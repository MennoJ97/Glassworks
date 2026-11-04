//! Closed domains: PEC or PMC walls, or a periodic cell, where modes come from the eigensolver.

use super::{coefficients, free_positions, grid_size, mode_volume, push_mesh};
use crate::assemble::{assemble, Csr, DofMap};
use crate::banded::Factors;
use crate::eigen::eigs_factored;
use crate::mesh::Mesh;
use crate::output::{normalized, FieldMeta, Meta, Output};
use crate::problem::{Analysis2d, Boundary, Polarization, Request2d};
use crate::space::H1Triangle;
use crate::sparse::{Ordering, SparseLu};
use crate::{C64, ZERO};
use std::collections::HashMap;
use std::f64::consts::PI;

/// Relative residual at which the eigensolver accepts a mode.
const EIG_TOL: f64 = 1e-10;
/// A band diagram reads only eigenvalues, whose error goes as the square of the residual, so
/// this still gives `k0` to about 1e-10, far below the mesh's error.
const BAND_TOL: f64 = 1e-5;

/// Dirichlet walls, or Bloch ties from the right and top edges to the left and
/// bottom. Corners tie straight to the origin, carrying both phases.
fn constraints(req: &Request2d, space: &H1Triangle, nx: usize, ny: usize, phase: [f64; 2]) -> DofMap {
    let positions = space.dof_positions();
    // Dofs sit on vertices and edge midpoints, so they lie on a grid twice as fine as the cells.
    // On the outer edge, which never moves, each one's position on that grid is exact.
    let (hx, hy) = (req.width / (2 * nx) as f64, req.height / (2 * ny) as f64);
    let on_right = |p: (f64, f64)| (p.0 / hx).round() as usize == 2 * nx;
    let on_top = |p: (f64, f64)| (p.1 / hy).round() as usize == 2 * ny;
    let on_left = |p: (f64, f64)| (p.0 / hx).round() as usize == 0;
    let on_bottom = |p: (f64, f64)| (p.1 / hy).round() as usize == 0;
    if req.walls != Boundary::Periodic {
        let wall = if req.polarization == Polarization::Te { Boundary::Pec } else { Boundary::Pmc };
        let zero: Vec<usize> = if req.walls == wall {
            (0..positions.len())
                .filter(|&dof| {
                    on_left(positions[dof])
                        || on_right(positions[dof])
                        || on_bottom(positions[dof])
                        || on_top(positions[dof])
                })
                .collect()
        } else {
            vec![]
        };
        return DofMap::new(positions.len(), &zero, &[]);
    }
    // A boundary node's grid position, with the right and top edges folded onto the left and bottom.
    let wrapped =
        |p: (f64, f64)| ((p.0 / hx).round() as usize % (2 * nx), (p.1 / hy).round() as usize % (2 * ny));
    let masters: HashMap<_, _> = (0..positions.len())
        .filter(|&dof| {
            !on_right(positions[dof])
                && !on_top(positions[dof])
                && (on_left(positions[dof]) || on_bottom(positions[dof]))
        })
        .map(|dof| (wrapped(positions[dof]), dof))
        .collect();
    let ties: Vec<_> = (0..positions.len())
        .filter(|&dof| on_right(positions[dof]) || on_top(positions[dof]))
        .map(|dof| {
            let angle = if on_right(positions[dof]) { phase[0] } else { 0.0 }
                + if on_top(positions[dof]) { phase[1] } else { 0.0 };
            (dof, masters[&wrapped(positions[dof])], C64::from_polar(1.0, angle))
        })
        .collect();
    DofMap::new(positions.len(), &[], &ties)
}

/// The `n` values of `(k0^2, x)` nearest `sigma`, for one set of constraints.
fn modes(
    space: &H1Triangle,
    dofs: &DofMap,
    coef: &[(C64, C64)],
    n: usize,
    sigma: f64,
    tol: f64,
) -> Result<Vec<(f64, Vec<C64>)>, String> {
    let region = &space.mesh.region;
    let a = Csr::from_triplets(dofs.n_free, assemble(space, dofs, |c| (coef[region[c]].0, ZERO)));
    let b = Csr::from_triplets(dofs.n_free, assemble(space, dofs, |c| (ZERO, coef[region[c]].1)));
    let shifted = a.sub_scaled(sigma, &b);
    let ordering = Ordering::reused(&shifted, || free_positions(space, dofs));
    eigs_factored(&SparseLu::factor(&shifted, &ordering) as &dyn Factors, &b, sigma, n, tol)
}

/// Modes of a closed cell, or a band diagram: the frequencies at Bloch wavevectors along the
/// edge of the irreducible zone.
pub(super) fn solve(req: &Request2d) -> Result<Output, String> {
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
            for (lambda, x) in
                modes(&space, &dofs, &coef, *n_modes, target_k0.map_or(sigma, |k0| k0 * k0), EIG_TOL)?
            {
                let i = out.meta.fields.len();
                let u = normalized(dofs.expand(&x));
                out.channels.push((format!("field{i}/{field}/re"), u.iter().map(|value| value.re).collect()));
                out.channels.push((format!("field{i}/{field}/im"), u.iter().map(|value| value.im).collect()));
                let volume = Some(mode_volume(&space, &coef, &u, |_| true));
                out.meta.fields.push(FieldMeta {
                    components: vec![field],
                    k0: lambda.max(0.0).sqrt(),
                    volume,
                    ..Default::default()
                });
            }
        }
        Analysis2d::Driven { .. } | Analysis2d::Resonance { .. } | Analysis2d::Grating { .. } => {
            unreachable!("sent elsewhere by solve2d::solve")
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
                let found = modes(&space, &dofs, &coef, *n_bands, sigma, BAND_TOL)?;
                for (band, (lambda, _)) in table.iter_mut().zip(found) {
                    band.push(lambda.max(0.0).sqrt());
                }
            }
            out.channels.push(("s".into(), path));
            out.channels.extend(table.into_iter().enumerate().map(|(i, band)| (format!("band{i}"), band)));
        }
    }
    Ok(out)
}
