//! An in-plane electric dipole, for `H_z`. A line source there is a magnetic current, while an
//! emitter in a planar device is an electric dipole, which reads `E.d = (i / k) p dH/dn`, with
//! `n` a quarter turn from `d`. Chapter 19 uses this to design cavities in `H_z`.

use super::open::{cell_load, cell_material, locate, pixel_of};
use crate::assemble::DofMap;
use crate::problem::{Density, Request2d};
use crate::space::{H1Triangle, Space};
use crate::C64;

/// A dipole's two points on the mesh, with each one's cell and weight. Its derivative is taken
/// across them, at mid-edge nodes on a grid, which needs no choice of a cell at the vertex between
/// them. The weight, without `p`, makes `F` in vacuum `k / 2`, as a line source's is.
fn points(
    req: &Request2d,
    space: &H1Triangle,
    (x, y, angle_deg): (f64, f64, f64),
    k0: f64,
) -> Result<[(f64, f64, usize, f64); 2], String> {
    let step = req.cell / 2.0;
    let (ny, nx) = angle_deg.to_radians().sin_cos();
    let weight = 2f64.sqrt() / k0 / (2.0 * step);
    let point = |sign: f64| {
        let (px, py) = (x - sign * step * ny, y + sign * step * nx);
        locate(space, px, py).map(|c| (px, py, c, sign * weight)).ok_or("the source lies outside the mesh")
    };
    Ok([point(1.0)?, point(-1.0)?])
}

/// The dipole as `e_0`, with `p` read at each of its points. `at` is on the mesh.
pub(super) fn load(
    req: &Request2d,
    space: &H1Triangle,
    dofs: &DofMap,
    coef: &[(C64, C64)],
    pml: f64,
    at: (f64, f64, f64),
    k0: f64,
) -> Result<Vec<(usize, f64)>, String> {
    let mut load = vec![];
    for (px, py, c, weight) in points(req, space, at, k0)? {
        let p = cell_material(req, space.mesh, coef, pml, c).0.re;
        load.extend(cell_load(space, dofs, c, px, py).into_iter().map(|(r, v)| (r, v * p * weight)));
    }
    Ok(load)
}

/// Since `e_0` holds `p`, a pixel at the dipole moves the source as well as the field. With
/// `A^-T e_0 = u / (i k)`, that adds `-2 u^T de_0` to the slope of the response, here in the
/// gradient's units of twice that, and with `dp = -deps / eps^2`.
#[allow(clippy::too_many_arguments)]
pub(super) fn add_source_slope(
    req: &Request2d,
    density: &Density,
    space: &H1Triangle,
    pml: f64,
    at: (f64, f64, f64),
    k0: f64,
    u: &[C64],
    gradient: &mut [f64],
    imaginary: &mut [f64],
) -> Result<(), String> {
    let mut d = vec![0; 6];
    for (px, py, c, weight) in points(req, space, at, k0)? {
        let Some(pixel) = pixel_of(req, space.mesh, pml, c) else { continue };
        space.cell_dofs(c, &mut d);
        let value: C64 = d.iter().zip(space.basis_at(c, px, py)).map(|(&g, v)| u[g] * v).sum();
        let deps = density.eps_high - density.eps_low;
        let eps = density.eps_low + density.rho[pixel] * deps;
        let change = value * 4.0 * weight * deps / (eps * eps);
        gradient[pixel] += change.re;
        imaginary[pixel] += change.im;
    }
    Ok(())
}
