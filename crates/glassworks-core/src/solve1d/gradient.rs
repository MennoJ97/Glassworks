//! The slope of R, T and A with respect to each layer's permittivity, for the mirror
//! designer on the website. Chapter 13 derives the same adjoint in 2D.

use super::{Driven, Setup};
use crate::banded::Factors;
use crate::output::Output;
use crate::problem::Polarization;
use crate::space::Space;
use crate::{C64, ONE, ZERO};

/// Adjoint gradients of R, T and A. For a real `J(u)` with `A u = b`,
/// `dJ/dθ = -2 Re(λᵀ (dA/dθ) u)`, where `Aᵀ λ = dJ/du`. Without periodic walls `A`
/// is complex symmetric, so the factors of `A` also solve for `λ`, and each output
/// costs one back-substitution however many layers there are.
///
/// Layers at an open end are left out: they also set the angle and the boundary
/// terms. The result is the exact gradient of the discrete problem on this mesh.
pub(super) fn power_gradients(setup: &Setup, k0: f64, solved: &Driven, out: &mut Output) {
    let last = setup.space.n_dofs() - 1;
    // λ for a `J` that depends only on the field at `dof`, where `dJ/du = dj_du` and is 0 elsewhere.
    let adjoint = |dof: usize, dj_du: C64| {
        let mut lambda = vec![ZERO; setup.dofs.n_free];
        if let Some((row, _)) = setup.dofs.map[dof] {
            lambda[row] = dj_du;
        }
        solved.lu.solve(&mut lambda);
        setup.dofs.expand(&lambda)
    };
    let dr = layer_sensitivity(setup, k0, solved, &adjoint(0, (solved.u[0] - ONE).conj()));
    let dt = layer_sensitivity(setup, k0, solved, &adjoint(last, solved.u[last].conj() * solved.t_scale));
    let da: Vec<C64> = dr.iter().zip(&dt).map(|(r, t)| -r - t).collect();
    let n_layers = setup.req.layers.len();
    let design = 1..if setup.req.right.is_open() { n_layers - 1 } else { n_layers };
    out.channels.push(("layer".into(), design.clone().map(|layer| layer as f64).collect()));
    for (name, per_layer) in [("R", &dr), ("T", &dt), ("A", &da)] {
        // `dε = i dε''` multiplies `λᵀ (dA/dε) u` by `i`, which turns Re into -Im.
        let re = design.clone().map(|layer| -2.0 * per_layer[layer].re).collect();
        let im = design.clone().map(|layer| 2.0 * per_layer[layer].im).collect();
        out.channels.push((format!("d{name}/deps_re"), re));
        out.channels.push((format!("d{name}/deps_im"), im));
    }
}

/// `λᵀ (dA/dε) u` for each layer. The derivative has the usual `a S + b M` form:
/// TE has `q = ε` in `b = kt^2 p - k0^2 q`, and TM has `p = 1/ε` in both `a` and `b`.
fn layer_sensitivity(setup: &Setup, k0: f64, solved: &Driven, lambda: &[C64]) -> Vec<C64> {
    let n = setup.space.local_size();
    let (mut stiff, mut unused, mut mass, mut cell_dofs) =
        (vec![0.0; n * n], vec![0.0; n * n], vec![0.0; n * n], vec![0; n]);
    let mut per_layer = vec![ZERO; setup.req.layers.len()];
    for cell in 0..setup.space.n_cells() {
        setup.space.cell_dofs(cell, &mut cell_dofs);
        setup.space.local_matrices(cell, &mut stiff, &mut unused, &mut mass);
        let (p, _) = setup.coef(cell);
        let (a, b) = match setup.req.polarization {
            Polarization::Te => (ZERO, C64::new(-k0 * k0, 0.0)),
            Polarization::Tm => (-p * p, -p * p * (solved.kt * solved.kt)),
        };
        let layer = setup.space.mesh.region[cell];
        for i in 0..n {
            for j in 0..n {
                let entry = a * stiff[i * n + j] + b * mass[i * n + j];
                per_layer[layer] += lambda[cell_dofs[i]] * entry * solved.u[cell_dofs[j]];
            }
        }
    }
    per_layer
}
