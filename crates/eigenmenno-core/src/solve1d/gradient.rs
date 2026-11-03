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
pub(super) fn power_gradients(s: &Setup, k0: f64, d: &Driven, out: &mut Output) {
    let last = s.space.n_dofs() - 1;
    let adjoint = |dof: usize, g: C64| {
        let mut x = vec![ZERO; s.dofs.n_free];
        if let Some((r, _)) = s.dofs.map[dof] {
            x[r] = g;
        }
        d.lu.solve(&mut x);
        s.dofs.expand(&x)
    };
    let dr = layer_sensitivity(s, k0, d, &adjoint(0, (d.u[0] - ONE).conj()));
    let dt = layer_sensitivity(s, k0, d, &adjoint(last, d.u[last].conj() * d.t_scale));
    let da: Vec<C64> = dr.iter().zip(&dt).map(|(r, t)| -r - t).collect();
    let n = s.req.layers.len();
    let design = 1..if s.req.right.is_open() { n - 1 } else { n };
    out.channels.push(("layer".into(), design.clone().map(|l| l as f64).collect()));
    for (name, w) in [("R", &dr), ("T", &dt), ("A", &da)] {
        // `dε = i dε''` multiplies `λᵀ (dA/dε) u` by `i`, which turns Re into -Im.
        out.channels.push((format!("d{name}/deps_re"), design.clone().map(|l| -2.0 * w[l].re).collect()));
        out.channels.push((format!("d{name}/deps_im"), design.clone().map(|l| 2.0 * w[l].im).collect()));
    }
}

/// `λᵀ (dA/dε) u` for each layer. The derivative has the usual `a S + b M` form:
/// TE has `q = ε` in `b = kt^2 p - k0^2 q`, and TM has `p = 1/ε` in both `a` and `b`.
fn layer_sensitivity(s: &Setup, k0: f64, d: &Driven, lambda: &[C64]) -> Vec<C64> {
    let n = s.space.local_size();
    let (mut stiff, mut unused, mut mass, mut dofs) =
        (vec![0.0; n * n], vec![0.0; n * n], vec![0.0; n * n], vec![0; n]);
    let mut w = vec![ZERO; s.req.layers.len()];
    for cell in 0..s.space.n_cells() {
        s.space.cell_dofs(cell, &mut dofs);
        s.space.local_matrices(cell, &mut stiff, &mut unused, &mut mass);
        let (p, _) = s.coef(cell);
        let (a, b) = match s.req.polarization {
            Polarization::Te => (ZERO, C64::new(-k0 * k0, 0.0)),
            Polarization::Tm => (-p * p, -p * p * (d.kt * d.kt)),
        };
        for i in 0..n {
            for j in 0..n {
                w[s.space.mesh.region[cell]] +=
                    lambda[dofs[i]] * (a * stiff[i * n + j] + b * mass[i * n + j]) * d.u[dofs[j]];
            }
        }
    }
    w
}
