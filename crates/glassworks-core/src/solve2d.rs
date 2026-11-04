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

mod adjoint;
mod closed;
mod dipole;
mod driven;
mod grating;
mod leaky;
mod open;
mod validate;

use crate::assemble::DofMap;
use crate::output::Output;
use crate::problem::{Analysis2d, Boundary, Polarization, Request2d};
use crate::space::{H1Triangle, Space};
use crate::C64;

/// Refuses meshes likely to exhaust WebAssembly's 4 GiB. An eigenmode solve, the heaviest,
/// took 3.2 GiB at 579,121 unknowns and 1.7 GiB at 314,721, and a driven one about half that.
const MAX_DOFS: usize = 400_000;

pub fn solve(req: &Request2d) -> Result<Output, String> {
    validate::validate(req)?;
    match req.analysis {
        Analysis2d::Driven { wavelength, source, pml, average_q } => {
            driven::solve(req, wavelength, source, pml, average_q)
        }
        Analysis2d::Resonance { target_k0, x, y, pml, dipole_deg } => {
            leaky::solve(req, target_k0, (x, y), pml, dipole_deg)
        }
        Analysis2d::Cavity { .. } | Analysis2d::Bands { .. } => closed::solve(req),
        Analysis2d::Grating { .. } => grating::solve(req),
    }
}

/// How far an absorbing layer extends the domain sideways: not at all for a grating.
fn side_layer(req: &Request2d, pml: f64) -> f64 {
    if req.walls == Boundary::Periodic {
        0.0
    } else {
        pml
    }
}

/// `(p, q)` of a region.
fn coefficients(req: &Request2d) -> Vec<(C64, C64)> {
    let materials = std::iter::once(req.background).chain(req.shapes.iter().map(|shape| shape.material));
    materials
        .map(|material| match req.polarization {
            Polarization::Te => (C64::new(1.0 / material.mu, 0.0), C64::new(material.eps, 0.0)),
            Polarization::Tm => (C64::new(1.0 / material.eps, 0.0), C64::new(material.mu, 0.0)),
        })
        .collect()
}

/// Where each free unknown is, on the low side of any periodic tie.
fn free_positions(space: &H1Triangle, dofs: &DofMap) -> Vec<(f64, f64)> {
    let all = space.dof_positions();
    let mut positions = vec![(f64::MAX, f64::MAX); dofs.n_free];
    for (dof, entry) in dofs.map.iter().enumerate() {
        if let Some((index, _)) = entry {
            if all[dof].0 + all[dof].1 < positions[*index].0 + positions[*index].1 {
                positions[*index] = all[dof];
            }
        }
    }
    positions
}

/// `V = integral(q |u|^2) / max(q |u|^2)` over the cells `keep` selects, an area here since
/// nothing varies along z. For `E_z` this is the usual mode volume per unit length. For
/// `H_z` it measures magnetic energy rather than the electric energy usually quoted.
fn mode_volume(space: &H1Triangle, coef: &[(C64, C64)], u: &[C64], keep: impl Fn(usize) -> bool) -> f64 {
    let (mut stiff_x, mut stiff_y, mut mass) = (vec![0.0; 36], vec![0.0; 36], vec![0.0; 36]);
    let mut cell_dofs = [0; 6];
    let (mut energy, mut peak) = (0.0, 0f64);
    for cell in (0..space.n_cells()).filter(|&cell| keep(cell)) {
        let q = coef[space.mesh.region[cell]].1.re;
        space.cell_dofs(cell, &mut cell_dofs);
        space.local_matrices(cell, &mut stiff_x, &mut stiff_y, &mut mass);
        let u_cell = cell_dofs.map(|dof| u[dof]);
        for i in 0..6 {
            peak = peak.max(q * u_cell[i].norm_sqr());
            energy += q * (0..6).map(|j| (u_cell[i].conj() * u_cell[j]).re * mass[i * 6 + j]).sum::<f64>();
        }
    }
    energy / peak
}

/// Steps of at most `step` that fill `size`. The tolerance keeps rounding noise from adding
/// one to a whole number of them, so a scaled structure gets the same mesh and sampling.
fn steps_in(size: f64, step: f64) -> usize {
    ((size / step - 1e-9).ceil() as usize).max(1)
}

/// Cells along x and y. Second-order triangles put a dof on every vertex and edge midpoint,
/// which is a grid twice as fine: `(2 nx + 1) (2 ny + 1)` dofs.
fn grid_size(req: &Request2d, width: f64, height: f64) -> Result<(usize, usize), String> {
    let (nx, ny) = (steps_in(width, req.cell), steps_in(height, req.cell));
    if (2 * nx + 1) * (2 * ny + 1) > MAX_DOFS {
        return Err("this mesh is too fine: use larger cells or a smaller domain".into());
    }
    Ok((nx, ny))
}

/// Dof positions, and four straight triangles per cell for linear plotting.
fn push_mesh(space: &H1Triangle, out: &mut Output) {
    let positions = space.dof_positions();
    out.channels.push(("x".into(), positions.iter().map(|p| p.0).collect()));
    out.channels.push(("y".into(), positions.iter().map(|p| p.1).collect()));
    let mut triangles = Vec::with_capacity(space.n_cells() * 12);
    let mut cell_dofs = vec![0; 6];
    for cell in 0..space.n_cells() {
        space.cell_dofs(cell, &mut cell_dofs);
        // Local dofs 0-2 are the corners and 3-5 the edge midpoints.
        for corners in [[0, 3, 5], [3, 1, 4], [5, 4, 2], [3, 4, 5]] {
            triangles.extend(corners.map(|local| cell_dofs[local] as f64));
        }
    }
    out.channels.push(("triangles".into(), triangles));
}
