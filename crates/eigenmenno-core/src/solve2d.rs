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

mod closed;
mod driven;
mod gradient;
mod open;
mod resonance;
mod validate;

use crate::assemble::DofMap;
use crate::problem::{Analysis2d, Polarization, Request2d};
use crate::solve::Output;
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
        Analysis2d::Resonance { target_k0, x, y, pml } => resonance::solve(req, target_k0, (x, y), pml),
        Analysis2d::Cavity { .. } | Analysis2d::Bands { .. } => closed::solve(req),
    }
}

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
