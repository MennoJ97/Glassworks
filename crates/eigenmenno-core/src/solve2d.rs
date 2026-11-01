//! Scalar problems on a plane: structures uniform along z, where `u = E_z` has
//! `p = 1/mu`, `q = eps` and `u = eta0 H_z` has `p = 1/eps`, `q = mu` in
//!
//! ```text
//! -div(p grad u) = k0^2 q u
//! ```
//!
//! This is the 1D equation with `d/dx` replaced by the gradient, so the weak form
//! and `assemble` are unchanged. `Te` and `Tm` in the request keep their 1D meaning
//! for these coefficients: `Te` solves for `E_z` and `Tm` for `H_z`.

use crate::assemble::{assemble, Csr, DofMap};
use crate::banded::Factors;
use crate::eigen::eigs_with;
use crate::mesh::Mesh;
use crate::problem::{Analysis2d, Boundary, Polarization, Request2d};
use crate::solve::{normalized, FieldMeta, Meta, Output};
use crate::space::{H1Triangle, Space};
use crate::sparse::SparseLu;
use crate::{C64, ZERO};
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

/// The `n` lowest `(k0^2, x)` for one set of constraints.
fn modes(space: &H1Triangle, dofs: &DofMap, coef: &[(C64, C64)], n: usize, sigma: f64) -> Result<Vec<(f64, Vec<C64>)>, String> {
    let region = &space.mesh.region;
    let a = Csr::from_triplets(dofs.n_free, assemble(space, dofs, |c| (coef[region[c]].0, ZERO)));
    let b = Csr::from_triplets(dofs.n_free, assemble(space, dofs, |c| (ZERO, coef[region[c]].1)));
    let lu = SparseLu::factor(&a.sub_scaled(sigma, &b))?;
    eigs_with(&lu as &dyn Factors, &b, sigma, n, EIG_TOL)
}

fn validate(req: &Request2d, nx: usize, ny: usize) -> Result<(), String> {
    if !(req.width > 0.0 && req.height > 0.0 && req.cell > 0.0) {
        return Err("the domain and its cell size must be positive".into());
    }
    if (2 * nx + 1) * (2 * ny + 1) > MAX_DOFS {
        return Err("this mesh is too fine: use larger cells or a smaller domain".into());
    }
    let positive = |m: crate::problem::Material| m.eps > 0.0 && m.mu > 0.0;
    if !positive(req.background) || !req.shapes.iter().all(|s| positive(s.material)) {
        return Err("eigenmode analyses support lossless materials with positive eps and mu".into());
    }
    if matches!(req.walls, Boundary::Absorbing | Boundary::Port) {
        return Err("2D eigenmode analyses need PEC, PMC or periodic walls".into());
    }
    if matches!(req.analysis, Analysis2d::Bands { .. }) && req.walls != Boundary::Periodic {
        return Err("band diagrams need periodic walls".into());
    }
    Ok(())
}

pub fn solve(req: &Request2d) -> Result<Output, String> {
    let (nx, ny) = ((req.width / req.cell).ceil() as usize, (req.height / req.cell).ceil() as usize);
    let (nx, ny) = (nx.max(1), ny.max(1));
    validate(req, nx, ny)?;
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
            for (lambda, x) in modes(&space, &dofs, &coef, *n_modes, sigma)? {
                let i = out.meta.fields.len();
                let u = normalized(dofs.expand(&x));
                out.channels.push((format!("field{i}/{field}/re"), u.iter().map(|v| v.re).collect()));
                out.channels.push((format!("field{i}/{field}/im"), u.iter().map(|v| v.im).collect()));
                out.meta.fields.push(FieldMeta { components: vec![field], k0: lambda.max(0.0).sqrt(), n_eff: None, power: None });
            }
        }
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
            for &s in &path {
                let dofs = constraints(req, &space, nx, ny, phases(s));
                out.meta.n_dofs = dofs.n_free;
                for (band, (lambda, _)) in table.iter_mut().zip(modes(&space, &dofs, &coef, *n_bands, sigma)?) {
                    band.push(lambda.max(0.0).sqrt());
                }
            }
            out.channels.push(("s".into(), path));
            out.channels.extend(table.into_iter().enumerate().map(|(i, b)| (format!("band{i}"), b)));
        }
    }
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
