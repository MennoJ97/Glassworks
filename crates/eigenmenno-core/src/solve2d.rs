//! Scalar problems on a plane: structures uniform along z, where `u = E_z` has
//! `p = 1/mu`, `q = eps` and `u = eta0 H_z` has `p = 1/eps`, `q = mu` in
//!
//! ```text
//! -div(p grad u) - k0^2 q u = source
//! ```
//!
//! This is the 1D equation with `d/dx` replaced by the gradient, so the weak form
//! and `assemble` are unchanged. `Te` and `Tm` in the request keep their 1D meaning
//! for these coefficients: `Te` solves for `E_z` and `Tm` for `H_z`.

use crate::assemble::{assemble, assemble_split, Csr, DofMap};
use crate::banded::Factors;
use crate::eigen::eigs_with;
use crate::mesh::Mesh;
use crate::problem::{Analysis, Analysis2d, Boundary, Layer, MeshSpec, Monitor, Polarization, Request, Request2d, Source2d};
use crate::solve::{normalized, FieldMeta, Meta, Output};
use crate::space::{H1Triangle, Space};
use crate::sparse::{Ordering, SparseLu};
use crate::{C64, I, ONE, ZERO};
use std::collections::HashMap;
use std::f64::consts::PI;

const EIG_TOL: f64 = 1e-10;
/// Refuses meshes likely to exhaust WebAssembly's 4 GiB. An eigenmode solve, the heaviest,
/// took 3.2 GiB at 579,121 unknowns and 1.7 GiB at 314,721, and a driven one about half that.
const MAX_DOFS: usize = 400_000;
/// A pole search that needs more driven solves than this has lost its way.
const MAX_POLE_STEPS: usize = 30;

/// `(p, q)` of a region.
fn coefficients(req: &Request2d) -> Vec<(C64, C64)> {
    let materials = std::iter::once(req.background).chain(req.shapes.iter().map(|s| s.material));
    materials
        .map(|m| match req.polarization {
            Polarization::Te => (C64::new(1.0 / m.mu, 0.0), C64::new(m.eps, 0.0)),
            Polarization::Tm => (C64::new(1.0 / m.eps, 0.0), C64::new(m.mu, 0.0)),
        })
        .collect()
}

/// Dirichlet walls, or Bloch ties from the right and top edges to the left and
/// bottom. Corners tie straight to the origin, carrying both phases.
fn constraints(req: &Request2d, space: &H1Triangle, nx: usize, ny: usize, phase: [f64; 2]) -> DofMap {
    let pos = space.dof_positions();
    let (hx, hy) = (req.width / (2 * nx) as f64, req.height / (2 * ny) as f64);
    let on_right = |p: (f64, f64)| (p.0 / hx).round() as usize == 2 * nx;
    let on_top = |p: (f64, f64)| (p.1 / hy).round() as usize == 2 * ny;
    let on_left = |p: (f64, f64)| (p.0 / hx).round() as usize == 0;
    let on_bottom = |p: (f64, f64)| (p.1 / hy).round() as usize == 0;
    if req.walls != Boundary::Periodic {
        let wall = if req.polarization == Polarization::Te { Boundary::Pec } else { Boundary::Pmc };
        let zero: Vec<usize> = if req.walls == wall {
            (0..pos.len()).filter(|&d| on_left(pos[d]) || on_right(pos[d]) || on_bottom(pos[d]) || on_top(pos[d])).collect()
        } else {
            vec![]
        };
        return DofMap::new(pos.len(), &zero, &[]);
    }
    let key = |p: (f64, f64)| ((p.0 / hx).round() as usize % (2 * nx), (p.1 / hy).round() as usize % (2 * ny));
    let masters: HashMap<_, _> =
        (0..pos.len()).filter(|&d| !on_right(pos[d]) && !on_top(pos[d]) && (on_left(pos[d]) || on_bottom(pos[d]))).map(|d| (key(pos[d]), d)).collect();
    let ties: Vec<_> = (0..pos.len())
        .filter(|&d| on_right(pos[d]) || on_top(pos[d]))
        .map(|d| {
            let angle = if on_right(pos[d]) { phase[0] } else { 0.0 } + if on_top(pos[d]) { phase[1] } else { 0.0 };
            (d, masters[&key(pos[d])], C64::from_polar(1.0, angle))
        })
        .collect();
    DofMap::new(pos.len(), &[], &ties)
}

/// Where each free unknown is, on the low side of any periodic tie.
fn free_positions(space: &H1Triangle, dofs: &DofMap) -> Vec<(f64, f64)> {
    let all = space.dof_positions();
    let mut pos = vec![(f64::MAX, f64::MAX); dofs.n_free];
    for (d, m) in dofs.map.iter().enumerate() {
        if let Some((r, _)) = m {
            if all[d].0 + all[d].1 < pos[*r].0 + pos[*r].1 {
                pos[*r] = all[d];
            }
        }
    }
    pos
}

/// The `n` values of `(k0^2, x)` nearest `sigma`, for one set of constraints. Bloch phases
/// change only the matrix values, so a band diagram shares one `ordering` across its k-points.
fn modes(
    space: &H1Triangle,
    dofs: &DofMap,
    coef: &[(C64, C64)],
    n: usize,
    sigma: f64,
    ordering: &mut Option<Ordering>,
) -> Result<Vec<(f64, Vec<C64>)>, String> {
    let region = &space.mesh.region;
    let a = Csr::from_triplets(dofs.n_free, assemble(space, dofs, |c| (coef[region[c]].0, ZERO)));
    let b = Csr::from_triplets(dofs.n_free, assemble(space, dofs, |c| (ZERO, coef[region[c]].1)));
    let shifted = a.sub_scaled(sigma, &b);
    let ordering = ordering.get_or_insert_with(|| Ordering::new(&shifted, &free_positions(space, dofs)));
    eigs_with(&SparseLu::factor(&shifted, ordering) as &dyn Factors, &b, sigma, n, EIG_TOL)
}

/// `V = integral(q |u|^2) / max(q |u|^2)` over the cells `keep` selects, an area here since
/// nothing varies along z. For `E_z` this is the usual mode volume per unit length. For
/// `H_z` it measures magnetic energy rather than the electric energy usually quoted.
fn mode_volume(space: &H1Triangle, coef: &[(C64, C64)], u: &[C64], keep: impl Fn(usize) -> bool) -> f64 {
    let (mut sx, mut sy, mut m, mut d) = (vec![0.0; 36], vec![0.0; 36], vec![0.0; 36], vec![0; 6]);
    let (mut energy, mut peak) = (0.0, 0f64);
    for c in (0..space.n_cells()).filter(|&c| keep(c)) {
        let q = coef[space.mesh.region[c]].1.re;
        space.cell_dofs(c, &mut d);
        space.local_matrices(c, &mut sx, &mut sy, &mut m);
        for i in 0..6 {
            peak = peak.max(q * u[d[i]].norm_sqr());
            energy += q * (0..6).map(|j| (u[d[i]].conj() * u[d[j]]).re * m[i * 6 + j]).sum::<f64>();
        }
    }
    energy / peak
}

/// Steps of at most `step` that fill `size`. The tolerance keeps rounding noise from adding
/// one to a whole number of them, so a scaled structure gets the same mesh and sampling.
fn steps_in(size: f64, step: f64) -> usize {
    ((size / step - 1e-9).ceil() as usize).max(1)
}

fn grid_size(req: &Request2d, width: f64, height: f64) -> Result<(usize, usize), String> {
    let (nx, ny) = (steps_in(width, req.cell), steps_in(height, req.cell));
    if (2 * nx + 1) * (2 * ny + 1) > MAX_DOFS {
        return Err("this mesh is too fine: use larger cells or a smaller domain".into());
    }
    Ok((nx, ny))
}

fn validate(req: &Request2d) -> Result<(), String> {
    if !(req.width > 0.0 && req.height > 0.0 && req.cell > 0.0) {
        return Err("the domain and its cell size must be positive".into());
    }
    let positive = |m: crate::problem::Material| m.eps > 0.0 && m.mu > 0.0;
    if !positive(req.background) || !req.shapes.iter().all(|s| positive(s.material)) {
        return Err("eigenmode analyses support lossless materials with positive eps and mu".into());
    }
    if let Some(d) = &req.density {
        let plane_wave = matches!(req.analysis, Analysis2d::Driven { source: Source2d::PlaneWave { .. }, .. });
        if !matches!(req.analysis, Analysis2d::Driven { .. }) || plane_wave || req.polarization != Polarization::Te || req.monitors.is_empty() {
            return Err("a density map needs a driven E_z solve from a line source or port, and a monitor".into());
        }
        if d.rho.len() != d.nx * d.ny || !(d.width > 0.0 && d.height > 0.0 && d.eps_low > 0.0 && d.eps_high > 0.0) {
            return Err("the density map's size or materials are invalid".into());
        }
    }
    let driven = matches!(req.analysis, Analysis2d::Driven { .. } | Analysis2d::Resonance { .. });
    if !driven && matches!(req.walls, Boundary::Absorbing | Boundary::Port) {
        return Err("2D eigenmode analyses need PEC, PMC or periodic walls".into());
    }
    if let Analysis2d::Resonance { target_k0, x, y, pml } = req.analysis {
        if req.walls != Boundary::Absorbing || !(target_k0 > 0.0 && pml > 0.0) {
            return Err("2D leaky modes need absorbing walls, a positive target and a positive PML".into());
        }
        if !(0.0..=req.width).contains(&x) || !(0.0..=req.height).contains(&y) {
            return Err("the source lies outside the domain".into());
        }
    }
    if let Analysis2d::Driven { wavelength, pml, source, average_q } = req.analysis {
        if req.walls != Boundary::Absorbing || !(wavelength > 0.0 && pml > 0.0) {
            return Err("2D driven solves need absorbing walls, a positive wavelength and a positive PML".into());
        }
        if average_q.is_some_and(|q| !(q > 0.0)) || (average_q.is_some() && matches!(source, Source2d::PlaneWave { .. })) {
            return Err("averaging over frequency needs a line source and a positive Q".into());
        }
        match source {
            Source2d::Line { x, y } if !(0.0..=req.width).contains(&x) || !(0.0..=req.height).contains(&y) => {
                return Err("the source lies outside the domain".into());
            }
            // The incident wave is defined in the background, so scatterers can't reach the absorber.
            Source2d::PlaneWave { .. } if !req.shapes.iter().all(|s| s.outline.inside(req.width, req.height)) => {
                return Err("shapes must stay inside the domain for a plane wave".into());
            }
            Source2d::Port { x, y0, y1 } if !(0.0..=req.width).contains(&x) || !(0.0 <= y0 && y0 < y1 && y1 <= req.height) => {
                return Err("the port lies outside the domain".into());
            }
            _ => {}
        }
    }
    if matches!(req.analysis, Analysis2d::Bands { .. }) && req.walls != Boundary::Periodic {
        return Err("band diagrams need periodic walls".into());
    }
    Ok(())
}

pub fn solve(req: &Request2d) -> Result<Output, String> {
    validate(req)?;
    if let Analysis2d::Driven { wavelength, source, pml, average_q } = req.analysis {
        return driven(req, wavelength, source, pml, average_q);
    }
    if let Analysis2d::Resonance { target_k0, x, y, pml } = req.analysis {
        return resonance(req, target_k0, (x, y), pml);
    }
    let (nx, ny) = grid_size(req, req.width, req.height)?;
    let outlines: Vec<_> = req.shapes.iter().map(|s| s.outline).collect();
    let mesh = Mesh::fitted_grid(nx, ny, req.width, req.height, &outlines);
    let space = H1Triangle::new(&mesh);
    let coef = coefficients(req);
    let length = req.width.max(req.height);
    // A negative shift keeps the operator definite when searching up from k0 = 0.
    let sigma = -(PI / length).powi(2);
    let mut out = Output { meta: Meta { length: req.width, ..Default::default() }, channels: vec![] };
    match &req.analysis {
        Analysis2d::Cavity { n_modes, bloch, target_k0 } => {
            let dofs = constraints(req, &space, nx, ny, *bloch);
            out.meta.n_dofs = dofs.n_free;
            push_mesh(&space, &mut out);
            let field = if req.polarization == Polarization::Te { "Ez" } else { "Hz" };
            for (lambda, x) in modes(&space, &dofs, &coef, *n_modes, target_k0.map_or(sigma, |k| k * k), &mut None)? {
                let i = out.meta.fields.len();
                let u = normalized(dofs.expand(&x));
                out.channels.push((format!("field{i}/{field}/re"), u.iter().map(|v| v.re).collect()));
                out.channels.push((format!("field{i}/{field}/im"), u.iter().map(|v| v.im).collect()));
                let volume = Some(mode_volume(&space, &coef, &u, |_| true));
                out.meta.fields.push(FieldMeta { components: vec![field], k0: lambda.max(0.0).sqrt(), volume, ..Default::default() });
            }
        }
        Analysis2d::Driven { .. } | Analysis2d::Resonance { .. } => unreachable!("handled above"),
        Analysis2d::Bands { n_bands, n_k } => {
            // Gamma to X to M to Gamma, with `s` running from 0 to 3 along the path.
            let points = 3 * n_k + 1;
            let phases = |s: f64| match s {
                s if s <= 1.0 => [PI * s, 0.0],
                s if s <= 2.0 => [PI, PI * (s - 1.0)],
                s => [PI * (3.0 - s), PI * (3.0 - s)],
            };
            let path: Vec<f64> = (0..points).map(|i| i as f64 / *n_k as f64).collect();
            let mut table = vec![vec![]; *n_bands];
            let mut ordering = None;
            for &s in &path {
                let dofs = constraints(req, &space, nx, ny, phases(s));
                out.meta.n_dofs = dofs.n_free;
                for (band, (lambda, _)) in table.iter_mut().zip(modes(&space, &dofs, &coef, *n_bands, sigma, &mut ordering)?) {
                    band.push(lambda.max(0.0).sqrt());
                }
            }
            out.channels.push(("s".into(), path));
            out.channels.extend(table.into_iter().enumerate().map(|(i, b)| (format!("band{i}"), b)));
        }
    }
    Ok(out)
}

/// A mesh extended by `pml` on every side, whose outer edge is a wall.
fn open_mesh(req: &Request2d, pml: f64) -> Result<Mesh, String> {
    let (width, height) = (req.width + 2.0 * pml, req.height + 2.0 * pml);
    let (nx, ny) = grid_size(req, width, height)?;
    let outlines: Vec<_> = req.shapes.iter().map(|s| s.outline.shifted(pml, pml)).collect();
    Ok(Mesh::fitted_grid(nx, ny, width, height, &outlines))
}

fn open_dofs(req: &Request2d, space: &H1Triangle, pml: f64) -> DofMap {
    let (width, height) = (req.width + 2.0 * pml, req.height + 2.0 * pml);
    let pos = space.dof_positions();
    let edge = |p: (f64, f64)| p.0 < 1e-9 || p.1 < 1e-9 || p.0 > width - 1e-9 || p.1 > height - 1e-9;
    let zero: Vec<usize> = (0..pos.len()).filter(|&d| edge(pos[d])).collect();
    DofMap::new(pos.len(), &zero, &[])
}

fn cell_center(mesh: &Mesh, c: usize) -> (f64, f64) {
    let v = &mesh.cells[3 * c..3 * c + 3];
    let (x, y) = v.iter().fold((0.0, 0.0), |s, &i| (s.0 + mesh.vertex(i).0, s.1 + mesh.vertex(i).1));
    (x / 3.0, y / 3.0)
}

/// The pixel of `req.density` holding a cell on an extended mesh, if any.
fn pixel_of(req: &Request2d, mesh: &Mesh, pml: f64, c: usize) -> Option<usize> {
    let d = req.density.as_ref()?;
    let p = cell_center(mesh, c);
    let (u, v) = ((p.0 - pml - d.x0) / d.width, (p.1 - pml - d.y0) / d.height);
    ((0.0..1.0).contains(&u) && (0.0..1.0).contains(&v)).then(|| (v * d.ny as f64) as usize * d.nx + (u * d.nx as f64) as usize)
}

/// The operator at any `k0`, real or complex. Coordinates in the layer stretch by
/// `s = 1 + i sigma / k`, which makes outgoing waves decay. With `sigma` growing
/// quadratically into the layer, a wave crossing it and back loses a factor of about 1e-4.
/// Asking for more steepens the layer, which reflects more where rods run through it.
/// That stretch is equivalent to weighting the x and y stiffness differently, with the
/// mass scaled by both stretches. It's analytic in `k0` above the real axis, as causality
/// requires, and it still absorbs a mode decaying slowly in time, just below that axis.
fn open_matrix(req: &Request2d, space: &H1Triangle, dofs: &DofMap, coef: &[(C64, C64)], pml: f64, k0: C64) -> Csr {
    let (pb, qb) = coef[0];
    // The background's wavenumber, from eps mu = q / p whichever polarization this is.
    let k = k0 * (qb.re / pb.re).sqrt();
    let sigma_max = -3.0 * (1e-4f64).ln() / (2.0 * pml);
    let stretch = |d: f64, size: f64| {
        let depth = (pml - d).max(d - pml - size).max(0.0) / pml;
        C64::new(1.0, 0.0) + I * (sigma_max * depth * depth) / k
    };
    let t = assemble_split(space, dofs, |c| {
        let (mut p, mut q) = coef[space.mesh.region[c]];
        if let Some(pixel) = pixel_of(req, space.mesh, pml, c) {
            let d = req.density.as_ref().unwrap();
            (p, q) = (ONE, C64::new(d.eps_low + d.rho[pixel] * (d.eps_high - d.eps_low), 0.0));
        }
        let (x, y) = cell_center(space.mesh, c);
        let (sx, sy) = (stretch(x, req.width), stretch(y, req.height));
        (p * sy / sx, p * sx / sy, -k0 * k0 * q * sx * sy)
    });
    Csr::from_triplets(dofs.n_free, t)
}

fn locate(space: &H1Triangle, x: f64, y: f64) -> Option<usize> {
    (0..space.n_cells()).find(|&c| space.barycentric(c, x, y).iter().all(|&l| l >= -1e-9))
}

/// `e_0`, which reads the field at a point as `sum(e_0 u)`, and loads it there as a source.
fn point_load(space: &H1Triangle, dofs: &DofMap, x: f64, y: f64) -> Result<Vec<(usize, f64)>, String> {
    let cell = locate(space, x, y).ok_or("the source lies outside the mesh")?;
    Ok(cell_load(space, dofs, cell, x, y))
}

fn cell_load(space: &H1Triangle, dofs: &DofMap, cell: usize, x: f64, y: f64) -> Vec<(usize, f64)> {
    let mut d = vec![0; 6];
    space.cell_dofs(cell, &mut d);
    d.iter().zip(space.basis_at(cell, x, y)).filter_map(|(&g, v)| dofs.map[g].map(|(r, _)| (r, v))).collect()
}

/// `u` for `A(k) u = i k e_0`, and `-sum(e_0 u)`. Dividing by `i k` turns this into the
/// Green's function. So `2 Re(-u_0)` is `pi` times the LDOS per unit k, or `k / 2` in vacuum.
fn line_response(lu: &SparseLu, n: usize, load: &[(usize, f64)], k: C64) -> (Vec<C64>, C64) {
    let mut u = vec![ZERO; n];
    load.iter().for_each(|&(r, v)| u[r] += I * k * v);
    lu.solve(&mut u);
    let response = -load.iter().map(|&(r, v)| u[r] * v).sum::<C64>();
    (u, response)
}

/// The guided mode of the layers along `x` from `y0` to `y1`, from the 1D solver: its
/// propagation constant is unused here, because the sheet's own fields set it. Returns the
/// profile as `(y, value)` samples, peaking at 1.
fn port_mode(req: &Request2d, k0: f64, x: f64, y0: f64, y1: f64) -> Result<Vec<(f64, f64)>, String> {
    let steps = steps_in(y1 - y0, req.cell / 8.0);
    let dy = (y1 - y0) / steps as f64;
    let material = |y: f64| req.shapes.iter().rev().find(|s| s.outline.level(x, y) <= 0.0).map_or(req.background, |s| s.material);
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
    let o = crate::solve::solve(&Request {
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
/// `P = u^H B u` for a Hermitian `B`, and the first monitor's `B u` goes into `first`.
fn monitor_power(space: &H1Triangle, coef: &[(C64, C64)], u: &[C64], k0: f64, cell: f64, pml: f64, monitors: &[Monitor], first: &mut [C64]) -> Result<Vec<f64>, String> {
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
                    let p = cell_center(space.mesh, c);
                    let (rx, ry) = (p.0 - pml - m.x0, p.1 - pml - m.y0);
                    let along = (rx * dx + ry * dy) / length;
                    (rx * nx + ry * ny).abs() < 4.0 * cell && along > -4.0 * cell && along < length + 4.0 * cell
                })
                .collect();
            let mut total = 0.0;
            for i in 0..steps {
                let t = (i as f64 + 0.5) / steps as f64;
                for offset in [-0.75, -0.25, 0.25, 0.75] {
                    let (x, y) = (m.x0 + t * dx + offset * cell * nx + pml, m.y0 + t * dy + offset * cell * ny + pml);
                    let c = *near.iter().find(|&&c| space.barycentric(c, x, y).iter().all(|&l| l >= -1e-9)).ok_or("a monitor lies outside the mesh")?;
                    space.cell_dofs(c, &mut d);
                    let (val, (gx, gy)) = (space.basis_at(c, x, y), space.gradients_at(c, x, y));
                    let (mut v, mut dn) = (ZERO, ZERO);
                    for j in 0..6 {
                        v += u[d[j]] * val[j];
                        dn += u[d[j]] * (gx[j] * nx + gy[j] * ny);
                    }
                    let weight = coef[space.mesh.region[c]].0.re * length / steps as f64 / 4.0;
                    total += weight * (v.conj() * dn).im;
                    if index == 0 {
                        for j in 0..6 {
                            first[d[j]] += (v * -(gx[j] * nx + gy[j] * ny) + dn * val[j]) * weight / (2.0 * k0) / (I * 2.0);
                        }
                    }
                }
            }
            Ok(total / (2.0 * k0))
        })
        .collect()
}

/// Solves in a domain surrounded by a perfectly matched layer. At `k = k0 (1 + i / 2Q)`,
/// analyticity makes `F(k)` exactly the average of `F` at real frequencies over a Lorentzian
/// of half-width `Im k`, and a resonance's `F` stays finite and smooth at any `Q`.
fn driven(req: &Request2d, wavelength: f64, source: Source2d, pml: f64, average_q: Option<f64>) -> Result<Output, String> {
    let mesh = open_mesh(req, pml)?;
    let space = H1Triangle::new(&mesh);
    let coef = coefficients(req);
    let k0 = 2.0 * PI / wavelength;
    let k = C64::new(k0, average_q.map_or(0.0, |q| k0 / (2.0 * q)));
    let dofs = open_dofs(req, &space, pml);
    let a = open_matrix(req, &space, &dofs, &coef, pml, k);
    let ordering = Ordering::reused(&a, &free_positions(&space, &dofs));
    let lu = SparseLu::factor(&a, &ordering);
    let (pb, qb) = coef[0];
    let mut ldos = None;
    let u = match source {
        Source2d::Line { x, y } => {
            let (u, response) = line_response(&lu, a.n, &point_load(&space, &dofs, x + pml, y + pml)?, k);
            ldos = Some(2.0 * response.re);
            dofs.expand(&u)
        }
        Source2d::PlaneWave { angle_deg } => {
            // Scattered-field form. The total field u = u_inc + u_s solves L u = 0, and the
            // incident wave solves the background's L_b u_inc = 0, so
            // L u_s = -(L - L_b) u_inc, which is nonzero only inside the shapes.
            let kb = k0 * (qb.re / pb.re).sqrt();
            let (kx, ky) = (kb * angle_deg.to_radians().cos(), kb * angle_deg.to_radians().sin());
            let inc: Vec<C64> = space.dof_positions().iter().map(|p| C64::from_polar(1.0, kx * (p.0 - pml) + ky * (p.1 - pml))).collect();
            let mut rhs = vec![ZERO; dofs.n_free];
            let (mut sx, mut sy, mut m, mut d) = (vec![0.0; 36], vec![0.0; 36], vec![0.0; 36], vec![0; 6]);
            for c in (0..space.n_cells()).filter(|&c| mesh.region[c] != 0) {
                let (p, q) = coef[mesh.region[c]];
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
        Source2d::Port { x, y0, y1 } => {
            // Modes are orthogonal with weight p, so the sheet `p phi` excites only this
            // one, with a field of `-(k / 2 beta) phi exp(i beta |x - x0|)`.
            let profile = port_mode(req, k0, x, y0, y1)?;
            let mut rhs = vec![ZERO; dofs.n_free];
            // Only cells near the port need searching.
            let near: Vec<usize> = (0..space.n_cells()).filter(|&c| (cell_center(&mesh, c).0 - pml - x).abs() < 4.0 * req.cell).collect();
            for pair in profile.windows(2).filter(|w| w[1].0 > w[0].0) {
                let y = (pair[0].0 + pair[1].0) / 2.0;
                let c = *near.iter().find(|&&c| space.barycentric(c, x + pml, y + pml).iter().all(|&l| l >= -1e-9)).ok_or("the port lies outside the mesh")?;
                let weight = coef[mesh.region[c]].0 * (pair[0].1 + pair[1].1) / 2.0 * (pair[1].0 - pair[0].0);
                for (r, v) in cell_load(&space, &dofs, c, x + pml, y + pml) {
                    rhs[r] += I * k * weight * v;
                }
            }
            lu.solve(&mut rhs);
            dofs.expand(&rhs)
        }
    };
    let meta = FieldMeta { k0, ldos, ..Default::default() };
    let mut out = open_output(req, &space, &dofs, pml, &u, meta);
    let mut bu = vec![ZERO; u.len()];
    out.channels.push(("monitors".into(), monitor_power(&space, &coef, &u, k0, req.cell, pml, &req.monitors, &mut bu)?));
    if let Some(density) = &req.density {
        // With `dP = 2 Re(u^H B du)` and `A du = -dA u`, an adjoint `A^T w = conj(B u)` gives
        // `dP = -2 Re(w^T dA u)`. The operator is complex symmetric, even in the layer, so
        // `A^T = A` and the same factors solve it. Each pixel's `dA` is `-k^2 deps M`.
        let mut w = vec![ZERO; dofs.n_free];
        for (g, v) in dofs.map.iter().zip(&bu) {
            if let Some((r, _)) = g {
                w[*r] += v.conj();
            }
        }
        lu.solve(&mut w);
        let w = dofs.expand(&w);
        let mut gradient = vec![0.0; density.rho.len()];
        let (mut sx, mut sy, mut m, mut d) = (vec![0.0; 36], vec![0.0; 36], vec![0.0; 36], vec![0; 6]);
        for c in 0..space.n_cells() {
            let Some(pixel) = pixel_of(req, &mesh, pml, c) else { continue };
            space.cell_dofs(c, &mut d);
            space.local_matrices(c, &mut sx, &mut sy, &mut m);
            let product: C64 = (0..36).map(|n| w[d[n / 6]] * m[n] * u[d[n % 6]]).sum();
            gradient[pixel] += 2.0 * k0 * k0 * (density.eps_high - density.eps_low) * product.re;
        }
        out.channels.push(("gradient".into(), gradient));
    }
    Ok(out)
}

/// Finds a leaky mode with no eigensolver. The line source's response has a pole at the
/// mode's complex `k`, so its reciprocal has a zero there, which a secant iteration finds
/// in a few driven solves. Close to the pole, the field is the mode.
fn resonance(req: &Request2d, target_k0: f64, (x, y): (f64, f64), pml: f64) -> Result<Output, String> {
    let mesh = open_mesh(req, pml)?;
    let space = H1Triangle::new(&mesh);
    let coef = coefficients(req);
    let dofs = open_dofs(req, &space, pml);
    let load = point_load(&space, &dofs, x + pml, y + pml)?;
    let ordering = Ordering::new(&open_matrix(req, &space, &dofs, &coef, pml, C64::new(target_k0, 0.0)), &free_positions(&space, &dofs));
    let at = |k: C64| {
        let a = open_matrix(req, &space, &dofs, &coef, pml, k);
        line_response(&SparseLu::factor(&a, &ordering), a.n, &load, k)
    };
    let mut prev = C64::new(target_k0, 0.0);
    let mut k = prev * C64::new(1.0, -2e-3);
    let (mut g_prev, mut g) = (at(prev).1.inv(), at(k).1.inv());
    for _ in 0..MAX_POLE_STEPS {
        let next = k - g * (k - prev) / (g - g_prev);
        if !next.is_finite() {
            break;
        }
        (prev, g_prev) = (k, g);
        k = next;
        let (u, response) = at(k);
        g = response.inv();
        if (k - prev).norm() <= 1e-4 * k.im.abs() {
            // A search can wander off to a pole of the mesh, far above anything physical.
            if !(k.im < 0.0 && (k.re - target_k0).abs() < 0.25 * target_k0) {
                break;
            }
            let u = normalized(dofs.expand(&u));
            let inside = |c: usize| {
                let p = cell_center(&mesh, c);
                p.0 > pml && p.1 > pml && p.0 < pml + req.width && p.1 < pml + req.height
            };
            let meta = FieldMeta {
                k0: k.re,
                q: Some(k.re / (2.0 * -k.im)),
                volume: Some(mode_volume(&space, &coef, &u, inside)),
                ..Default::default()
            };
            return Ok(open_output(req, &space, &dofs, pml, &u, meta));
        }
    }
    Err("no leaky mode found near this wavelength: try a different target, or a source where the mode is strong".into())
}

fn open_output(req: &Request2d, space: &H1Triangle, dofs: &DofMap, pml: f64, u: &[C64], meta: FieldMeta) -> Output {
    let mut out = Output { meta: Meta { length: req.width, n_dofs: dofs.n_free, ..Default::default() }, channels: vec![] };
    push_mesh(space, &mut out);
    // Back to the caller's coordinates, so the layer sits outside `[0, width]`.
    for name in ["x", "y"] {
        out.channels.iter_mut().find(|c| c.0 == name).unwrap().1.iter_mut().for_each(|v| *v -= pml);
    }
    let field = if req.polarization == Polarization::Te { "Ez" } else { "Hz" };
    out.channels.push((format!("field0/{field}/re"), u.iter().map(|v| v.re).collect()));
    out.channels.push((format!("field0/{field}/im"), u.iter().map(|v| v.im).collect()));
    out.meta.fields.push(FieldMeta { components: vec![field], ..meta });
    out
}

/// Dof positions, and four straight triangles per cell for linear plotting.
fn push_mesh(space: &H1Triangle, out: &mut Output) {
    let pos = space.dof_positions();
    out.channels.push(("x".into(), pos.iter().map(|p| p.0).collect()));
    out.channels.push(("y".into(), pos.iter().map(|p| p.1).collect()));
    let mut triangles = Vec::with_capacity(space.n_cells() * 12);
    let mut d = vec![0; 6];
    for c in 0..space.n_cells() {
        space.cell_dofs(c, &mut d);
        for t in [[0, 3, 5], [3, 1, 4], [5, 4, 2], [3, 4, 5]] {
            triangles.extend(t.map(|k| d[k] as f64));
        }
    }
    out.channels.push(("triangles".into(), triangles));
}
