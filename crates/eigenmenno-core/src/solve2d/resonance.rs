//! Leaky modes of an open domain, found from the poles of a line source's response.

use super::open::{cell_center, line_response, open_dofs, open_matrix, open_mesh, open_output, point_load};
use super::{coefficients, free_positions, mode_volume};
use crate::problem::Request2d;
use crate::solve::{normalized, FieldMeta, Output};
use crate::space::H1Triangle;
use crate::sparse::{Ordering, SparseLu};
use crate::C64;

/// A pole search that needs more driven solves than this has lost its way.
const MAX_POLE_STEPS: usize = 30;

/// Finds a leaky mode with no eigensolver. The line source's response has a pole at the
/// mode's complex `k`, so its reciprocal has a zero there, which a secant iteration finds
/// in a few driven solves. Close to the pole, the field is the mode.
pub(super) fn solve(req: &Request2d, target_k0: f64, (x, y): (f64, f64), pml: f64) -> Result<Output, String> {
    let mesh = open_mesh(req, pml)?;
    let space = H1Triangle::new(&mesh);
    let coef = coefficients(req);
    let dofs = open_dofs(req, &space, pml);
    let load = point_load(&space, &dofs, x + pml, y + pml)?;
    let at = |k: C64| {
        let a = open_matrix(req, &space, &dofs, &coef, pml, k);
        let ordering = Ordering::reused(&a, || free_positions(&space, &dofs));
        line_response(&SparseLu::factor(&a, &ordering), a.n, &load, k)
    };
    let mut prev = C64::new(target_k0, 0.0);
    let mut k = prev * C64::new(1.0, -2e-3);
    let (mut g_prev, mut g) = (at(prev).1.inv(), at(k).1.inv());
    for _ in 0..MAX_POLE_STEPS {
        let next = k - g * (k - prev) / (g - g_prev);
        if !next.is_finite() {
            break;
        }
        (prev, g_prev) = (k, g);
        k = next;
        let (u, response) = at(k);
        g = response.inv();
        if (k - prev).norm() <= 1e-4 * k.im.abs() {
            // A search can wander off to a pole of the mesh, far above anything physical.
            if !(k.im < 0.0 && (k.re - target_k0).abs() < 0.25 * target_k0) {
                break;
            }
            let u = normalized(dofs.expand(&u));
            let inside = |c: usize| {
                let p = cell_center(&mesh, c);
                p.0 > pml && p.1 > pml && p.0 < pml + req.width && p.1 < pml + req.height
            };
            let meta = FieldMeta {
                k0: k.re,
                q: Some(k.re / (2.0 * -k.im)),
                volume: Some(mode_volume(&space, &coef, &u, inside)),
                ..Default::default()
            };
            return Ok(open_output(req, &space, &dofs, pml, &u, meta));
        }
    }
    Err("no leaky mode found near this wavelength: try a different target, or a source where the mode is strong".into())
}
