//! A faster factorization when `A^T = A`, as an absorbing layer's operator is. Half of every
//! update to a front repeats the other half, so this skips it, and pivoting with it, as
//! `A = L D L^T`. Refinement covers its weak pivots, and a pivot too small to trust falls back
//! to the pivoted `factor`. Chapter 16 has the measurements.

use super::blocked::{subtract, BLOCK};
use super::{Csr, Ordering, SparseLu, C64, ZERO};

impl<'a> SparseLu<'a> {
    /// For `A^T = A`. About half the cost of `factor`, with the same accuracy after refinement.
    pub fn factor_symmetric(a: &'a Csr, ordering: &'a Ordering) -> Self {
        Self::factor_as(a, ordering, true).unwrap_or_else(|| Self::factor(a, ordering))
    }
}

/// As `partial_lu`, for a complex symmetric front, updating only its lower triangle and leaving
/// `U = D L^T` for the solves. Returns whether any pivot was weak, or `None` if one is unusable.
pub(super) fn partial_ldl(f: &mut [C64], m: usize, p: usize, scale: f64) -> Option<bool> {
    let mut weak = false;
    let mut k = 0;
    while k < p {
        let e = (k + BLOCK).min(p);
        for j in k..e {
            let d = f[j + j * m];
            let below = f[j * m + j + 1..(j + 1) * m].iter().fold(0.0f64, |b, v| b.max(v.norm_sqr()));
            if d.norm_sqr() < 1e-16 * below || d.norm_sqr() < 1e-28 * scale * scale {
                return None;
            }
            weak |= d.norm_sqr() < 1e-4 * below || d.norm_sqr() < 1e-16 * scale * scale;
            let inv = d.inv();
            f[j * m + j + 1..(j + 1) * m].iter_mut().for_each(|v| *v *= inv);
            let (left, right) = f.split_at_mut((j + 1) * m);
            for (c, col) in (j + 1..e).zip(right.chunks_exact_mut(m)) {
                let u = d * left[j * m + c];
                col[c..].iter_mut().zip(&left[j * m + c..]).for_each(|(y, x)| *y -= u * *x);
            }
        }
        // Update the later columns from this block.
        let (left, right) = f.split_at_mut(e * m);
        for (c, col) in (e..m).zip(right.chunks_exact_mut(m)) {
            let mut q = k;
            while q < e {
                let take = (e - q).min(4);
                let at = |o: usize| {
                    let from = if o < take { q + o } else { q };
                    &left[from * m + c..(from + 1) * m]
                };
                let coef = |o: usize| {
                    if o < take {
                        left[(q + o) * m + q + o] * left[(q + o) * m + c]
                    } else {
                        ZERO
                    }
                };
                subtract(&mut col[c..], [at(0), at(1), at(2), at(3)], [coef(0), coef(1), coef(2), coef(3)]);
                q += take;
            }
        }
        k = e;
    }
    for c in 1..m {
        for i in 0..c.min(p) {
            f[i + c * m] = f[i + i * m] * f[c + i * m];
        }
    }
    Some(weak)
}
