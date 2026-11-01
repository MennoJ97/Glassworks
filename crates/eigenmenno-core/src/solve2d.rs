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
use crate::problem::{Analysis2d, Boundary, Polarization, Request2d, Source2d};
use crate::solve::{normalized, FieldMeta, Meta, Output};
use crate::space::{H1Triangle, Space};
use crate::sparse::{SparseCholesky, SparseLu};
use crate::{C64, I, ZERO};
use faer::sparse::linalg::solvers::SymbolicLlt;
use std::collections::HashMap;
use std::f64::consts::PI;

const EIG_TOL: f64 = 1e-10;
/// Refuses meshes likely to exhaust a browser tab's memory: 638,401 unknowns did.
const MAX_DOFS: usize = 160_000;

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

/// The `n` lowest `(k0^2, x)` for one set of constraints. Bloch phases change only the
/// matrix values, so a band diagram shares one `analysis` across its k-points.
fn modes(
    space: &H1Triangle,
    dofs: &DofMap,
    coef: &[(C64, C64)],
    n: usize,
    sigma: f64,
    analysis: &mut Option<SymbolicLlt<usize>>,
) -> Result<Vec<(f64, Vec<C64>)>, String> {
    let region = &space.mesh.region;
    let a = Csr::from_triplets(dofs.n_free, assemble(space, dofs, |c| (coef[region[c]].0, ZERO)));
    let b = Csr::from_triplets(dofs.n_free, assemble(space, dofs, |c| (ZERO, coef[region[c]].1)));
    let factors = SparseCholesky::factor(&a.sub_scaled(sigma, &b), analysis)?;
    eigs_with(&factors as &dyn Factors, &b, sigma, n, EIG_TOL)
}

fn grid_size(req: &Request2d, width: f64, height: f64) -> Result<(usize, usize), String> {
    let (nx, ny) = (((width / req.cell).ceil() as usize).max(1), ((height / req.cell).ceil() as usize).max(1));
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
    let driven = matches!(req.analysis, Analysis2d::Driven { .. });
    if !driven && matches!(req.walls, Boundary::Absorbing | Boundary::Port) {
        return Err("2D eigenmode analyses need PEC, PMC or periodic walls".into());
    }
    if let Analysis2d::Driven { wavelength, pml, source } = req.analysis {
        if req.walls != Boundary::Absorbing || !(wavelength > 0.0 && pml > 0.0) {
            return Err("2D driven solves need absorbing walls, a positive wavelength and a positive PML".into());
        }
        match source {
            Source2d::Line { x, y } if !(0.0..=req.width).contains(&x) || !(0.0..=req.height).contains(&y) => {
                return Err("the source lies outside the domain".into());
            }
            // The incident wave is defined in the background, so scatterers can't reach the absorber.
            Source2d::PlaneWave { .. } if !req.shapes.iter().all(|s| s.outline.inside(req.width, req.height)) => {
                return Err("shapes must stay inside the domain for a plane wave".into());
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
    if let Analysis2d::Driven { wavelength, source, pml } = req.analysis {
        return driven(req, wavelength, source, pml);
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
        Analysis2d::Cavity { n_modes, bloch } => {
            let dofs = constraints(req, &space, nx, ny, *bloch);
            out.meta.n_dofs = dofs.n_free;
            push_mesh(&space, &mut out);
            let field = if req.polarization == Polarization::Te { "Ez" } else { "Hz" };
            for (lambda, x) in modes(&space, &dofs, &coef, *n_modes, sigma, &mut None)? {
                let i = out.meta.fields.len();
                let u = normalized(dofs.expand(&x));
                out.channels.push((format!("field{i}/{field}/re"), u.iter().map(|v| v.re).collect()));
                out.channels.push((format!("field{i}/{field}/im"), u.iter().map(|v| v.im).collect()));
                out.meta.fields.push(FieldMeta { components: vec![field], k0: lambda.max(0.0).sqrt(), n_eff: None, power: None });
            }
        }
        Analysis2d::Driven { .. } => unreachable!("handled above"),
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
            let mut analysis = None;
            for &s in &path {
                let dofs = constraints(req, &space, nx, ny, phases(s));
                out.meta.n_dofs = dofs.n_free;
                for (band, (lambda, _)) in table.iter_mut().zip(modes(&space, &dofs, &coef, *n_bands, sigma, &mut analysis)?) {
                    band.push(lambda.max(0.0).sqrt());
                }
            }
            out.channels.push(("s".into(), path));
            out.channels.extend(table.into_iter().enumerate().map(|(i, b)| (format!("band{i}"), b)));
        }
    }
    Ok(out)
}

/// Solves in a domain surrounded by a perfectly matched layer, on a mesh extended by
/// `pml` on every side. Coordinates in the layer stretch by `s = 1 + i sigma / k`, which
/// makes outgoing waves decay. With `sigma` growing quadratically into the layer, a wave
/// crossing it and back loses a factor of about 1e-6. That stretch is equivalent to
/// weighting the x and y stiffness differently, with the mass scaled by both stretches.
fn driven(req: &Request2d, wavelength: f64, source: Source2d, pml: f64) -> Result<Output, String> {
    let (width, height) = (req.width + 2.0 * pml, req.height + 2.0 * pml);
    let (nx, ny) = grid_size(req, width, height)?;
    let outlines: Vec<_> = req.shapes.iter().map(|s| s.outline.shifted(pml, pml)).collect();
    let mesh = Mesh::fitted_grid(nx, ny, width, height, &outlines);
    let space = H1Triangle::new(&mesh);
    let coef = coefficients(req);
    let k0 = 2.0 * PI / wavelength;
    let k2 = C64::new(k0 * k0, 0.0);
    let (pb, qb) = coef[0];
    // The background's wavenumber, from eps mu = q / p whichever polarization this is.
    let k = k0 * (qb.re / pb.re).sqrt();
    let sigma_max = -3.0 * (1e-6f64).ln() / (2.0 * pml);
    let stretch = |d: f64, size: f64| {
        let depth = (pml - d).max(d - pml - size).max(0.0) / pml;
        C64::new(1.0, sigma_max * depth * depth / k)
    };
    let pos = space.dof_positions();
    let edge = |p: (f64, f64)| p.0 < 1e-9 || p.1 < 1e-9 || p.0 > width - 1e-9 || p.1 > height - 1e-9;
    let zero: Vec<usize> = (0..pos.len()).filter(|&d| edge(pos[d])).collect();
    let dofs = DofMap::new(pos.len(), &zero, &[]);
    let cell_center = |c: usize| {
        let v = &mesh.cells[3 * c..3 * c + 3];
        let (x, y) = v.iter().fold((0.0, 0.0), |s, &i| (s.0 + mesh.vertex(i).0, s.1 + mesh.vertex(i).1));
        (x / 3.0, y / 3.0)
    };
    let t = assemble_split(&space, &dofs, |c| {
        let (p, q) = coef[mesh.region[c]];
        let (x, y) = cell_center(c);
        let (sx, sy) = (stretch(x, req.width), stretch(y, req.height));
        (p * sy / sx, p * sx / sy, -k2 * q * sx * sy)
    });
    let mut rhs = vec![ZERO; dofs.n_free];
    let mut d = vec![0; 6];
    let mut incident_wave: Option<Vec<C64>> = None;
    match source {
        Source2d::Line { x, y } => {
            let (x, y) = (x + pml, y + pml);
            let cell = (0..space.n_cells())
                .find(|&c| space.barycentric(c, x, y).iter().all(|&l| l >= -1e-9))
                .ok_or("the source lies outside the mesh")?;
            space.cell_dofs(cell, &mut d);
            for (&g, v) in d.iter().zip(space.basis_at(cell, x, y)) {
                if let Some((r, f)) = dofs.map[g] {
                    rhs[r] += f.conj() * I * k0 * v;
                }
            }
        }
        Source2d::PlaneWave { angle_deg } => {
            // Scattered-field form. The total field u = u_inc + u_s solves L u = 0, and the
            // incident wave solves the background's L_b u_inc = 0, so
            // L u_s = -(L - L_b) u_inc, which is nonzero only inside the shapes.
            let (kx, ky) = (k * angle_deg.to_radians().cos(), k * angle_deg.to_radians().sin());
            let inc: Vec<C64> = pos.iter().map(|p| C64::from_polar(1.0, kx * (p.0 - pml) + ky * (p.1 - pml))).collect();
            let (mut sx, mut sy, mut m) = (vec![0.0; 36], vec![0.0; 36], vec![0.0; 36]);
            for c in (0..space.n_cells()).filter(|&c| mesh.region[c] != 0) {
                let (p, q) = coef[mesh.region[c]];
                let (dp, dq) = (p - pb, q - qb);
                space.cell_dofs(c, &mut d);
                space.local_matrices(c, &mut sx, &mut sy, &mut m);
                for i in 0..6 {
                    let Some((r, f)) = dofs.map[d[i]] else { continue };
                    for j in 0..6 {
                        let k = i * 6 + j;
                        rhs[r] -= f.conj() * (dp * (sx[k] + sy[k]) - k2 * dq * m[k]) * inc[d[j]];
                    }
                }
            }
            incident_wave = Some(inc);
        }
    }
    let a = Csr::from_triplets(dofs.n_free, t);
    SparseLu::factor(&a)?.solve(&mut rhs);
    let mut u = dofs.expand(&rhs);
    if let Some(inc) = incident_wave {
        u.iter_mut().zip(inc).for_each(|(v, i)| *v += i);
    }
    let mut out = Output { meta: Meta { length: req.width, n_dofs: dofs.n_free, ..Default::default() }, channels: vec![] };
    push_mesh(&space, &mut out);
    // Back to the caller's coordinates, so the layer sits outside `[0, width]`.
    for name in ["x", "y"] {
        out.channels.iter_mut().find(|c| c.0 == name).unwrap().1.iter_mut().for_each(|v| *v -= pml);
    }
    let field = if req.polarization == Polarization::Te { "Ez" } else { "Hz" };
    out.channels.push((format!("field0/{field}/re"), u.iter().map(|v| v.re).collect()));
    out.channels.push((format!("field0/{field}/im"), u.iter().map(|v| v.im).collect()));
    out.meta.fields.push(FieldMeta { components: vec![field], k0, n_eff: None, power: None });
    Ok(out)
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
