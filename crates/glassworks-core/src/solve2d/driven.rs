//! Driven solves in an open domain: a line source, a plane wave or a port, and what monitors read.

use super::adjoint::density_gradient;
use super::dipole;
use super::open::{cell_load, line_response, open_dofs, open_matrix, open_mesh, open_output, point_load};
use super::{coefficients, free_positions, steps_in};
use crate::assemble::DofMap;
use crate::banded::Factors;
use crate::output::{FieldMeta, Output};
use crate::problem::{
    Analysis, Boundary, Layer, MeshSpec, Monitor, Polarization, Request, Request2d, Source2d,
};
use crate::space::{H1Triangle, Space};
use crate::sparse::{Ordering, SparseLu};
use crate::{C64, I, ZERO};
use std::f64::consts::PI;

/// Solves in a domain surrounded by a perfectly matched layer. At `k = k0 (1 + i / 2Q)`,
/// analyticity makes `F(k)` exactly the average of `F` at real frequencies over a Lorentzian
/// of half-width `Im k`, and a resonance's `F` stays finite and smooth at any `Q`.
pub(super) fn solve(
    req: &Request2d,
    wavelength: f64,
    source: Source2d,
    pml: f64,
    average_q: Option<f64>,
) -> Result<Output, String> {
    let mesh = open_mesh(req, pml)?;
    let space = H1Triangle::new(&mesh);
    let coef = coefficients(req);
    let k0 = 2.0 * PI / wavelength;
    let k = C64::new(k0, average_q.map_or(0.0, |q| k0 / (2.0 * q)));
    let dofs = open_dofs(req, &space, pml);
    let a = open_matrix(req, &space, &dofs, &coef, pml, k);
    let ordering = Ordering::reused(&a, || free_positions(&space, &dofs));
    let lu = SparseLu::factor_symmetric(&a, &ordering);
    let (mut ldos, mut measured) = (None, None);
    let u = match source {
        Source2d::Line { x, y } => {
            let (u, response) = line_response(&lu, a.n, &point_load(&space, &dofs, x + pml, y + pml)?, k);
            ldos = Some(2.0 * response.re);
            measured = Some([response.re, response.im]);
            dofs.expand(&u)
        }
        Source2d::Dipole { x, y, angle_deg } => {
            let load = dipole::load(req, &space, &dofs, &coef, pml, (x + pml, y + pml, angle_deg), k0)?;
            let (u, response) = line_response(&lu, a.n, &load, k);
            ldos = Some(2.0 * response.re);
            measured = Some([response.re, response.im]);
            dofs.expand(&u)
        }
        Source2d::PlaneWave { angle_deg } => plane_wave(&space, &dofs, &coef, &lu, k0, k, angle_deg, pml),
        Source2d::Port { x, y0, y1 } => port(req, &space, &dofs, &coef, &lu, k0, k, (x, y0, y1), pml)?,
    };
    let meta = FieldMeta { k0, ldos, response: measured, ..Default::default() };
    let mut out = open_output(req, &space, &dofs, pml, &u, meta);
    let mut bu = vec![ZERO; u.len()];
    out.channels.push((
        "monitors".into(),
        monitor_power(&space, &coef, &u, k0, req.cell, pml, &req.monitors, &mut bu)?,
    ));
    if let Some(density) = &req.density {
        let (mut gradient, mut imaginary) =
            density_gradient(req, density, &space, &dofs, &lu, pml, k, &u, &bu);
        if let Source2d::Dipole { x, y, angle_deg } = source {
            let at = (x + pml, y + pml, angle_deg);
            dipole::add_source_slope(req, density, &space, pml, at, k0, &u, &mut gradient, &mut imaginary)?;
        }
        out.channels.push(("gradient".into(), gradient));
        // For the LDOS, the two together are the slope of `-2 u_0`, and so of a resonance's `k`.
        out.channels.push(("gradient_im".into(), imaginary));
    }
    Ok(out)
}

/// Scattered-field form. The total field u = u_inc + u_s solves L u = 0, and the
/// incident wave solves the background's L_b u_inc = 0, so
/// L u_s = -(L - L_b) u_inc, which is nonzero only inside the shapes.
#[allow(clippy::too_many_arguments)]
fn plane_wave(
    space: &H1Triangle,
    dofs: &DofMap,
    coef: &[(C64, C64)],
    lu: &SparseLu,
    k0: f64,
    k: C64,
    angle_deg: f64,
    pml: f64,
) -> Vec<C64> {
    let (pb, qb) = coef[0];
    let kb = k0 * (qb.re / pb.re).sqrt();
    let (kx, ky) = (kb * angle_deg.to_radians().cos(), kb * angle_deg.to_radians().sin());
    let inc: Vec<C64> = space
        .dof_positions()
        .iter()
        .map(|p| C64::from_polar(1.0, kx * (p.0 - pml) + ky * (p.1 - pml)))
        .collect();
    let mut rhs = vec![ZERO; dofs.n_free];
    let (mut sx, mut sy, mut m, mut d) = (vec![0.0; 36], vec![0.0; 36], vec![0.0; 36], vec![0; 6]);
    for c in (0..space.n_cells()).filter(|&c| space.mesh.region[c] != 0) {
        let (p, q) = coef[space.mesh.region[c]];
        let (dp, dq) = (p - pb, q - qb);
        space.cell_dofs(c, &mut d);
        space.local_matrices(c, &mut sx, &mut sy, &mut m);
        for i in 0..6 {
            let Some((r, f)) = dofs.map[d[i]] else { continue };
            for j in 0..6 {
                let n = i * 6 + j;
                rhs[r] -= f.conj() * (dp * (sx[n] + sy[n]) - k * k * dq * m[n]) * inc[d[j]];
            }
        }
    }
    lu.solve(&mut rhs);
    let mut u = dofs.expand(&rhs);
    u.iter_mut().zip(inc).for_each(|(v, i)| *v += i);
    u
}

/// Modes are orthogonal with weight p, so the sheet `p phi` excites only this
/// one, with a field of `-(k / 2 beta) phi exp(i beta |x - x0|)`.
#[allow(clippy::too_many_arguments)]
fn port(
    req: &Request2d,
    space: &H1Triangle,
    dofs: &DofMap,
    coef: &[(C64, C64)],
    lu: &SparseLu,
    k0: f64,
    k: C64,
    (x, y0, y1): (f64, f64, f64),
    pml: f64,
) -> Result<Vec<C64>, String> {
    let profile = port_mode(req, k0, x, y0, y1)?;
    let mut rhs = vec![ZERO; dofs.n_free];
    // Only cells near the port need searching.
    let near: Vec<usize> = (0..space.n_cells())
        .filter(|&c| (space.mesh.centroid(c).0 - pml - x).abs() < 4.0 * req.cell)
        .collect();
    for pair in profile.windows(2).filter(|w| w[1].0 > w[0].0) {
        let y = (pair[0].0 + pair[1].0) / 2.0;
        let c = *near
            .iter()
            .find(|&&c| space.contains(c, x + pml, y + pml))
            .ok_or("the port lies outside the mesh")?;
        let weight = coef[space.mesh.region[c]].0 * (pair[0].1 + pair[1].1) / 2.0 * (pair[1].0 - pair[0].0);
        for (r, v) in cell_load(space, dofs, c, x + pml, y + pml) {
            rhs[r] += I * k * weight * v;
        }
    }
    lu.solve(&mut rhs);
    Ok(dofs.expand(&rhs))
}

/// The guided mode of the layers along `x` from `y0` to `y1`, from the 1D solver: its
/// propagation constant is unused here, because the sheet's own fields set it. Returns the
/// profile as `(y, value)` samples, peaking at 1.
fn port_mode(req: &Request2d, k0: f64, x: f64, y0: f64, y1: f64) -> Result<Vec<(f64, f64)>, String> {
    let steps = steps_in(y1 - y0, req.cell / 8.0);
    let dy = (y1 - y0) / steps as f64;
    let material = |y: f64| {
        req.shapes.iter().rev().find(|s| s.outline.level(x, y) <= 0.0).map_or(req.background, |s| s.material)
    };
    let mut layers: Vec<Layer> = vec![];
    let mut last = None;
    for i in 0..steps {
        let m = material(y0 + (i as f64 + 0.5) * dy);
        if last == Some((m.eps, m.mu)) {
            layers.last_mut().unwrap().thickness += dy;
        } else {
            layers.push(Layer { thickness: dy, eps: [m.eps, 0.0], mu: [m.mu, 0.0] });
            last = Some((m.eps, m.mu));
        }
    }
    let wavelength = 2.0 * PI / k0;
    let o = crate::solve1d::solve(&Request {
        layers,
        polarization: req.polarization,
        left: Boundary::Pec,
        right: Boundary::Pec,
        mesh: MeshSpec { order: 4, max_cell: wavelength / 8.0 },
        analysis: Analysis::Waveguide { wavelength, n_modes: 1 },
    })?;
    if o.meta.fields.is_empty() {
        return Err("no guided mode at the port: it needs a layer of higher index than both ends".into());
    }
    let channel = |name: &str| &o.channels.iter().find(|c| c.0 == name).unwrap().1;
    let name = if req.polarization == Polarization::Te { "field0/Ey/re" } else { "field0/Hy/re" };
    Ok(channel("field0/x").iter().zip(channel(name)).map(|(s, v)| (y0 + s, *v)).collect())
}

/// Time-averaged power crossing each monitor, `integral(p Im(conj(u) du/dn)) / (2 k0)`,
/// in units where a unit plane wave in vacuum carries 1/2 per unit length. A derivative is
/// least accurate at a cell's edge, where a line along the grid would sit. Power is conserved
/// across a thin strip, so this averages it over two cells, which cancels that error.
/// `P = u^H B u` for a Hermitian `B`, and the first monitor's `B u` goes into `bu`.
fn monitor_power(
    space: &H1Triangle,
    coef: &[(C64, C64)],
    u: &[C64],
    k0: f64,
    cell: f64,
    pml: f64,
    monitors: &[Monitor],
    bu: &mut [C64],
) -> Result<Vec<f64>, String> {
    let mut d = vec![0; 6];
    monitors
        .iter()
        .enumerate()
        .map(|(index, m)| {
            let (dx, dy) = (m.x1 - m.x0, m.y1 - m.y0);
            let length = dx.hypot(dy);
            let steps = steps_in(length, cell / 4.0);
            let (nx, ny) = (dy / length, -dx / length);
            // Only cells near the strip need searching.
            let near: Vec<usize> = (0..space.n_cells())
                .filter(|&c| {
                    let p = space.mesh.centroid(c);
                    let (rx, ry) = (p.0 - pml - m.x0, p.1 - pml - m.y0);
                    let along = (rx * dx + ry * dy) / length;
                    (rx * nx + ry * ny).abs() < 4.0 * cell
                        && along > -4.0 * cell
                        && along < length + 4.0 * cell
                })
                .collect();
            let mut total = 0.0;
            for i in 0..steps {
                let t = (i as f64 + 0.5) / steps as f64;
                for offset in [-0.75, -0.25, 0.25, 0.75] {
                    let (x, y) =
                        (m.x0 + t * dx + offset * cell * nx + pml, m.y0 + t * dy + offset * cell * ny + pml);
                    let c = *near
                        .iter()
                        .find(|&&c| space.contains(c, x, y))
                        .ok_or("a monitor lies outside the mesh")?;
                    space.cell_dofs(c, &mut d);
                    let (val, (gx, gy)) = (space.basis_at(c, x, y), space.gradients_at(c, x, y));
                    let (mut v, mut dn) = (ZERO, ZERO);
                    for j in 0..6 {
                        v += u[d[j]] * val[j];
                        dn += u[d[j]] * (gx[j] * nx + gy[j] * ny);
                    }
                    let weight = coef[space.mesh.region[c]].0.re * length / steps as f64 / 4.0;
                    total += weight * (v.conj() * dn).im;
                    // `(B u)_j` is the slope of `P` along `conj(u_j)`, from `Im z = (z - conj(z)) / 2i`.
                    if index == 0 {
                        for j in 0..6 {
                            bu[d[j]] += (v * -(gx[j] * nx + gy[j] * ny) + dn * val[j]) * weight
                                / (2.0 * k0)
                                / (I * 2.0);
                        }
                    }
                }
            }
            Ok(total / (2.0 * k0))
        })
        .collect()
}
