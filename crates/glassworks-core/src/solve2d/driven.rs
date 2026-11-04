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
            let placed = (x + pml, y + pml, angle_deg);
            dipole::add_source_slope(
                req,
                density,
                &space,
                pml,
                placed,
                k0,
                &u,
                &mut gradient,
                &mut imaginary,
            )?;
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
    let (p_background, q_background) = coef[0];
    let k_background = k0 * (q_background.re / p_background.re).sqrt();
    let (kx, ky) = (k_background * angle_deg.to_radians().cos(), k_background * angle_deg.to_radians().sin());
    let incident: Vec<C64> = space
        .dof_positions()
        .iter()
        .map(|p| C64::from_polar(1.0, kx * (p.0 - pml) + ky * (p.1 - pml)))
        .collect();
    let mut rhs = vec![ZERO; dofs.n_free];
    let (mut stiff_x, mut stiff_y, mut mass) = (vec![0.0; 36], vec![0.0; 36], vec![0.0; 36]);
    let mut cell_dofs = vec![0; 6];
    for cell in (0..space.n_cells()).filter(|&cell| space.mesh.region[cell] != 0) {
        let (p, q) = coef[space.mesh.region[cell]];
        let (dp, dq) = (p - p_background, q - q_background);
        space.cell_dofs(cell, &mut cell_dofs);
        space.local_matrices(cell, &mut stiff_x, &mut stiff_y, &mut mass);
        for i in 0..6 {
            let Some((row, factor)) = dofs.map[cell_dofs[i]] else { continue };
            for j in 0..6 {
                let local = i * 6 + j;
                let entry = dp * (stiff_x[local] + stiff_y[local]) - k * k * dq * mass[local];
                rhs[row] -= factor.conj() * entry * incident[cell_dofs[j]];
            }
        }
    }
    lu.solve(&mut rhs);
    let mut u = dofs.expand(&rhs);
    u.iter_mut().zip(incident).for_each(|(total, incident)| *total += incident);
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
        .filter(|&cell| (space.mesh.centroid(cell).0 - pml - x).abs() < 4.0 * req.cell)
        .collect();
    for pair in profile.windows(2).filter(|pair| pair[1].0 > pair[0].0) {
        let y = (pair[0].0 + pair[1].0) / 2.0;
        let cell = *near
            .iter()
            .find(|&&cell| space.contains(cell, x + pml, y + pml))
            .ok_or("the port lies outside the mesh")?;
        let weight =
            coef[space.mesh.region[cell]].0 * (pair[0].1 + pair[1].1) / 2.0 * (pair[1].0 - pair[0].0);
        for (row, value) in cell_load(space, dofs, cell, x + pml, y + pml) {
            rhs[row] += I * k * weight * value;
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
        let topmost = req.shapes.iter().rev().find(|shape| shape.outline.level(x, y) <= 0.0);
        topmost.map_or(req.background, |shape| shape.material)
    };
    let mut layers: Vec<Layer> = vec![];
    let mut previous = None;
    for i in 0..steps {
        let here = material(y0 + (i as f64 + 0.5) * dy);
        if previous == Some((here.eps, here.mu)) {
            layers.last_mut().unwrap().thickness += dy;
        } else {
            layers.push(Layer { thickness: dy, eps: [here.eps, 0.0], mu: [here.mu, 0.0] });
            previous = Some((here.eps, here.mu));
        }
    }
    let wavelength = 2.0 * PI / k0;
    let solved = crate::solve1d::solve(&Request {
        layers,
        polarization: req.polarization,
        left: Boundary::Pec,
        right: Boundary::Pec,
        mesh: MeshSpec { order: 4, max_cell: wavelength / 8.0 },
        analysis: Analysis::Waveguide { wavelength, n_modes: 1 },
    })?;
    if solved.meta.fields.is_empty() {
        return Err("no guided mode at the port: it needs a layer of higher index than both ends".into());
    }
    let channel = |name: &str| &solved.channels.iter().find(|channel| channel.0 == name).unwrap().1;
    let name = if req.polarization == Polarization::Te { "field0/Ey/re" } else { "field0/Hy/re" };
    Ok(channel("field0/x").iter().zip(channel(name)).map(|(offset, value)| (y0 + offset, *value)).collect())
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
    cell_size: f64,
    pml: f64,
    monitors: &[Monitor],
    bu: &mut [C64],
) -> Result<Vec<f64>, String> {
    let mut cell_dofs = vec![0; 6];
    monitors
        .iter()
        .enumerate()
        .map(|(index, monitor)| {
            let (dx, dy) = (monitor.x1 - monitor.x0, monitor.y1 - monitor.y0);
            let length = dx.hypot(dy);
            let steps = steps_in(length, cell_size / 4.0);
            // The unit normal, which power crossing the strip flows along.
            let (normal_x, normal_y) = (dy / length, -dx / length);
            // Only cells near the strip need searching.
            let near: Vec<usize> = (0..space.n_cells())
                .filter(|&cell| {
                    let p = space.mesh.centroid(cell);
                    let (rx, ry) = (p.0 - pml - monitor.x0, p.1 - pml - monitor.y0);
                    let along = (rx * dx + ry * dy) / length;
                    (rx * normal_x + ry * normal_y).abs() < 4.0 * cell_size
                        && along > -4.0 * cell_size
                        && along < length + 4.0 * cell_size
                })
                .collect();
            let mut total = 0.0;
            for i in 0..steps {
                let t = (i as f64 + 0.5) / steps as f64;
                for offset in [-0.75, -0.25, 0.25, 0.75] {
                    let x = monitor.x0 + t * dx + offset * cell_size * normal_x + pml;
                    let y = monitor.y0 + t * dy + offset * cell_size * normal_y + pml;
                    let cell = *near
                        .iter()
                        .find(|&&cell| space.contains(cell, x, y))
                        .ok_or("a monitor lies outside the mesh")?;
                    space.cell_dofs(cell, &mut cell_dofs);
                    let (values, (d_dx, d_dy)) = (space.basis_at(cell, x, y), space.gradients_at(cell, x, y));
                    let d_dn = |j: usize| d_dx[j] * normal_x + d_dy[j] * normal_y;
                    let (mut value, mut du_dn) = (ZERO, ZERO);
                    for j in 0..6 {
                        value += u[cell_dofs[j]] * values[j];
                        du_dn += u[cell_dofs[j]] * d_dn(j);
                    }
                    let weight = coef[space.mesh.region[cell]].0.re * length / steps as f64 / 4.0;
                    total += weight * (value.conj() * du_dn).im;
                    // `(B u)_j` is the slope of `P` along `conj(u_j)`, from `Im z = (z - conj(z)) / 2i`.
                    if index == 0 {
                        for j in 0..6 {
                            bu[cell_dofs[j]] +=
                                (value * -d_dn(j) + du_dn * values[j]) * weight / (2.0 * k0) / (I * 2.0);
                        }
                    }
                }
            }
            Ok(total / (2.0 * k0))
        })
        .collect()
}
