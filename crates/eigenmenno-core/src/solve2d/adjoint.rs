//! The slope of a goal with respect to every density pixel, by the adjoint method.

use super::open::pixel_of;
use crate::assemble::DofMap;
use crate::banded::Factors;
use crate::problem::{Density, Request2d};
use crate::space::{H1Triangle, Space};
use crate::sparse::SparseLu;
use crate::{C64, I, ZERO};

/// Real and imaginary parts of each pixel's slope. `bu` is `B u` from the first monitor.
#[allow(clippy::too_many_arguments)]
pub(super) fn density_gradient(
    req: &Request2d,
    density: &Density,
    space: &H1Triangle,
    dofs: &DofMap,
    lu: &SparseLu,
    pml: f64,
    k: C64,
    u: &[C64],
    bu: &[C64],
) -> (Vec<f64>, Vec<f64>) {
    // With `dP = 2 Re(u^H B du)` and `A du = -dA u`, an adjoint `A^T w = conj(B u)` gives
    // `dP = -2 Re(w^T dA u)`. The operator is complex symmetric, even in the layer, so
    // `A^T = A` and the same factors solve it. Each pixel's `dA` is `-k^2 deps M`.
    // With no monitor, the goal is the LDOS `F = -2 Re(e_0^T u)`, for which `A^T w = -e_0`.
    // Since `A u = i k e_0`, that's `w = -u / (i k)`, with no further solve.
    let w = if req.monitors.is_empty() {
        u.iter().map(|v| -v / (I * k)).collect()
    } else {
        let mut w = vec![ZERO; dofs.n_free];
        for (g, v) in dofs.map.iter().zip(bu) {
            if let Some((r, _)) = g {
                w[*r] += v.conj();
            }
        }
        lu.solve(&mut w);
        dofs.expand(&w)
    };
    pixel_slopes(req, density, space, pml, k, &w, u)
}

/// `-2 w^T dA u` for each pixel, whose `dA` is `-k^2 deps M`. `w` and `u` are on all dofs.
pub(super) fn pixel_slopes(
    req: &Request2d,
    density: &Density,
    space: &H1Triangle,
    pml: f64,
    k: C64,
    w: &[C64],
    u: &[C64],
) -> (Vec<f64>, Vec<f64>) {
    let (mut gradient, mut imaginary) = (vec![0.0; density.rho.len()], vec![0.0; density.rho.len()]);
    let (mut sx, mut sy, mut m, mut d) = (vec![0.0; 36], vec![0.0; 36], vec![0.0; 36], vec![0; 6]);
    for c in 0..space.n_cells() {
        let Some(pixel) = pixel_of(req, space.mesh, pml, c) else { continue };
        space.cell_dofs(c, &mut d);
        space.local_matrices(c, &mut sx, &mut sy, &mut m);
        let product: C64 = (0..36).map(|n| w[d[n / 6]] * m[n] * u[d[n % 6]]).sum();
        let change = 2.0 * (density.eps_high - density.eps_low) * k * k * product;
        gradient[pixel] += change.re;
        imaginary[pixel] += change.im;
    }
    (gradient, imaginary)
}
