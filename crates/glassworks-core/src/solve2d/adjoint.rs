//! The slope of a goal with respect to every density pixel, by the adjoint method.

use super::open::pixel_of;
use crate::assemble::DofMap;
use crate::banded::Factors;
use crate::problem::{Density, Polarization, Request2d};
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
    // `A^T = A` and the same factors solve it.
    // With no monitor, the goal is the LDOS `F = -2 Re(e_0^T u)`, for which `A^T w = -e_0`.
    // Since `A u = i k e_0`, that's `w = -u / (i k)`, with no further solve.
    let w = if req.monitors.is_empty() {
        u.iter().map(|value| -value / (I * k)).collect()
    } else {
        let mut w = vec![ZERO; dofs.n_free];
        for (entry, bu_dof) in dofs.map.iter().zip(bu) {
            if let Some((row, _)) = entry {
                w[*row] += bu_dof.conj();
            }
        }
        lu.solve(&mut w);
        dofs.expand(&w)
    };
    pixel_slopes(req, density, space, pml, k, &w, u)
}

/// `-2 w^T dA u` for each pixel, with `w` and `u` on all dofs. A pixel's `dA` is `-k^2 deps M`
/// for `E_z`, where it mixes `q = eps`, and `-(deps / eps^2) S` for `H_z`, where `p = 1 / eps`.
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
    let (mut stiff_x, mut stiff_y, mut mass) = (vec![0.0; 36], vec![0.0; 36], vec![0.0; 36]);
    let mut cell_dofs = vec![0; 6];
    for cell in 0..space.n_cells() {
        let Some(pixel) = pixel_of(req, space.mesh, pml, cell) else { continue };
        space.cell_dofs(cell, &mut cell_dofs);
        space.local_matrices(cell, &mut stiff_x, &mut stiff_y, &mut mass);
        // `w^T local u`, for one of the cell's 6 x 6 matrices.
        let quadratic = |local: &[f64]| {
            (0..36).map(|n| w[cell_dofs[n / 6]] * local[n] * u[cell_dofs[n % 6]]).sum::<C64>()
        };
        let deps = density.eps_high - density.eps_low;
        let change = match req.polarization {
            Polarization::Te => 2.0 * deps * k * k * quadratic(&mass),
            Polarization::Tm => {
                let eps = density.eps_low + density.rho[pixel] * deps;
                (quadratic(&stiff_x) + quadratic(&stiff_y)) * 2.0 * deps / (eps * eps)
            }
        };
        gradient[pixel] += change.re;
        imaginary[pixel] += change.im;
    }
    (gradient, imaginary)
}
