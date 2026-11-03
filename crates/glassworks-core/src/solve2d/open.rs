//! An open domain: the structure surrounded by a perfectly matched layer, whose outer edge is a
//! wall. Everything shared by driven solves and the search for leaky modes.

use super::{grid_size, push_mesh, side_layer};
use crate::assemble::{assemble_split, Csr, DofMap};
use crate::banded::Factors;
use crate::mesh::Mesh;
use crate::output::{FieldMeta, Meta, Output};
use crate::problem::{Polarization, Request2d};
use crate::space::{H1Triangle, Space};
use crate::sparse::SparseLu;
use crate::{C64, I, ONE, ZERO};

/// A mesh extended by `pml` on every side that has a layer, whose outer edge is a wall.
pub(super) fn open_mesh(req: &Request2d, pml: f64) -> Result<Mesh, String> {
    let side = side_layer(req, pml);
    let (width, height) = (req.width + 2.0 * side, req.height + 2.0 * pml);
    let (nx, ny) = grid_size(req, width, height)?;
    let outlines: Vec<_> = req.shapes.iter().map(|s| s.outline.shifted(side, pml)).collect();
    Ok(Mesh::fitted_grid(nx, ny, width, height, &outlines))
}

pub(super) fn open_dofs(req: &Request2d, space: &H1Triangle, pml: f64) -> DofMap {
    let (width, height) = (req.width + 2.0 * pml, req.height + 2.0 * pml);
    let pos = space.dof_positions();
    let edge = |p: (f64, f64)| p.0 < 1e-9 || p.1 < 1e-9 || p.0 > width - 1e-9 || p.1 > height - 1e-9;
    let zero: Vec<usize> = (0..pos.len()).filter(|&d| edge(pos[d])).collect();
    DofMap::new(pos.len(), &zero, &[])
}

/// The pixel of `req.density` holding a cell on an extended mesh, if any.
pub(super) fn pixel_of(req: &Request2d, mesh: &Mesh, pml: f64, c: usize) -> Option<usize> {
    let d = req.density.as_ref()?;
    let p = mesh.centroid(c);
    let (u, v) = ((p.0 - side_layer(req, pml) - d.x0) / d.width, (p.1 - pml - d.y0) / d.height);
    ((0.0..1.0).contains(&u) && (0.0..1.0).contains(&v))
        .then(|| (v * d.ny as f64) as usize * d.nx + (u * d.nx as f64) as usize)
}

/// `(p, q)` of a cell, where a density pixel mixes permittivities: `q = eps` for `E_z`, and
/// `p = 1 / eps` for `H_z`.
pub(super) fn cell_material(
    req: &Request2d,
    mesh: &Mesh,
    coef: &[(C64, C64)],
    pml: f64,
    c: usize,
) -> (C64, C64) {
    let Some(pixel) = pixel_of(req, mesh, pml, c) else { return coef[mesh.region[c]] };
    let d = req.density.as_ref().unwrap();
    let eps = C64::new(d.eps_low + d.rho[pixel] * (d.eps_high - d.eps_low), 0.0);
    match req.polarization {
        Polarization::Te => (ONE, eps),
        Polarization::Tm => (eps.inv(), ONE),
    }
}

/// The operator at any `k0`, real or complex. Coordinates in the layer stretch by
/// `s = 1 + i sigma / k`, which makes outgoing waves decay. With `sigma` growing
/// quadratically into the layer, a wave crossing it and back loses a factor of about 1e-4.
/// Asking for more steepens the layer, which reflects more where rods run through it.
/// That stretch is equivalent to weighting the x and y stiffness differently, with the
/// mass scaled by both stretches. It's analytic in `k0` above the real axis, as causality
/// requires, and it still absorbs a mode decaying slowly in time, just below that axis.
pub(super) fn open_matrix(
    req: &Request2d,
    space: &H1Triangle,
    dofs: &DofMap,
    coef: &[(C64, C64)],
    pml: f64,
    k0: C64,
) -> Csr {
    stretched_matrix(req, space, dofs, coef, pml, k0, ONE)
}

/// `open_matrix` with its stretch multiplied by `strength`, which a grating tunes to its orders.
pub(super) fn stretched_matrix(
    req: &Request2d,
    space: &H1Triangle,
    dofs: &DofMap,
    coef: &[(C64, C64)],
    pml: f64,
    k0: C64,
    strength: C64,
) -> Csr {
    let (pb, qb) = coef[0];
    // The background's wavenumber, from eps mu = q / p whichever polarization this is.
    let kb = k0 * (qb.re / pb.re).sqrt();
    // `sigma = sigma_max (depth / pml)^2`, so across the layer `integral(sigma) = sigma_max pml / 3`,
    // and a round trip at normal incidence keeps `exp(-2 sigma_max pml / 3) = 1e-4` of the amplitude.
    let sigma_max = -3.0 * (1e-4f64).ln() / (2.0 * pml);
    let stretch = |d: f64, size: f64, layer: f64| {
        let depth = (layer - d).max(d - layer - size).max(0.0) / pml;
        C64::new(1.0, 0.0) + I * strength * (sigma_max * depth * depth) / kb
    };
    let t = assemble_split(space, dofs, |c| {
        let (p, q) = cell_material(req, space.mesh, coef, pml, c);
        let (x, y) = space.mesh.centroid(c);
        let (sx, sy) = (stretch(x, req.width, side_layer(req, pml)), stretch(y, req.height, pml));
        (p * sy / sx, p * sx / sy, -k0 * k0 * q * sx * sy)
    });
    Csr::from_triplets(dofs.n_free, t)
}

pub(super) fn locate(space: &H1Triangle, x: f64, y: f64) -> Option<usize> {
    (0..space.n_cells()).find(|&c| space.contains(c, x, y))
}

/// `e_0`, which reads the field at a point as `sum(e_0 u)`, and loads it there as a source.
pub(super) fn point_load(
    space: &H1Triangle,
    dofs: &DofMap,
    x: f64,
    y: f64,
) -> Result<Vec<(usize, f64)>, String> {
    let cell = locate(space, x, y).ok_or("the source lies outside the mesh")?;
    Ok(cell_load(space, dofs, cell, x, y))
}

pub(super) fn cell_load(space: &H1Triangle, dofs: &DofMap, cell: usize, x: f64, y: f64) -> Vec<(usize, f64)> {
    let mut d = vec![0; 6];
    space.cell_dofs(cell, &mut d);
    d.iter().zip(space.basis_at(cell, x, y)).filter_map(|(&g, v)| dofs.map[g].map(|(r, _)| (r, v))).collect()
}

/// `u` for `A(k) u = i k e_0`, and `-sum(e_0 u)`. Dividing by `i k` turns this into the
/// Green's function. So `2 Re(-u_0)` is `pi` times the LDOS per unit k, or `k / 2` in vacuum.
pub(super) fn line_response(lu: &SparseLu, n: usize, load: &[(usize, f64)], k: C64) -> (Vec<C64>, C64) {
    let mut u = vec![ZERO; n];
    load.iter().for_each(|&(r, v)| u[r] += I * k * v);
    lu.solve(&mut u);
    let response = -load.iter().map(|&(r, v)| u[r] * v).sum::<C64>();
    (u, response)
}

pub(super) fn open_output(
    req: &Request2d,
    space: &H1Triangle,
    dofs: &DofMap,
    pml: f64,
    u: &[C64],
    meta: FieldMeta,
) -> Output {
    let mut out = Output {
        meta: Meta { length: req.width, n_dofs: dofs.n_free, ..Default::default() },
        channels: vec![],
    };
    push_mesh(space, &mut out);
    // Back to the caller's coordinates, so the layer sits outside `[0, width]`.
    for (name, layer) in [("x", side_layer(req, pml)), ("y", pml)] {
        out.channels.iter_mut().find(|c| c.0 == name).unwrap().1.iter_mut().for_each(|v| *v -= layer);
    }
    let field = if req.polarization == Polarization::Te { "Ez" } else { "Hz" };
    out.channels.push((format!("field0/{field}/re"), u.iter().map(|v| v.re).collect()));
    out.channels.push((format!("field0/{field}/im"), u.iter().map(|v| v.im).collect()));
    out.meta.fields.push(FieldMeta { components: vec![field], ..meta });
    out
}
