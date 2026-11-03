//! The speed in `partial_lu`: it eliminates a block of columns, then updates every later
//! column from the whole block at once, four columns at a time, while the block is in cache.
//! The arithmetic is the same as updating after each column, so a first reading can skip this.

use super::C64;
use crate::ZERO;

/// Columns eliminated together before updating the rest of a front, to stay in cache.
pub(super) const BLOCK: usize = 32;

/// `y -= sum(u[q] x[q])`, over four columns at once for speed.
#[inline]
pub(super) fn subtract(y: &mut [C64], x: [&[C64]; 4], u: [C64; 4]) {
    let n = y.len();
    let [x0, x1, x2, x3] = [&x[0][..n], &x[1][..n], &x[2][..n], &x[3][..n]];
    for i in 0..n {
        y[i] -= u[0] * x0[i] + u[1] * x1[i] + u[2] * x2[i] + u[3] * x3[i];
    }
}

/// Brings every column after `e` up to date with the eliminated columns `k..e`: first the
/// rows inside the block, then the rest of each column, four eliminated columns at a time.
#[inline]
pub(super) fn update_after(f: &mut [C64], m: usize, k: usize, e: usize) {
    let (left, right) = f.split_at_mut(e * m);
    for col in right.chunks_exact_mut(m) {
        for q in k..e {
            let u = col[q];
            for i in q + 1..e {
                col[i] -= u * left[q * m + i];
            }
        }
        let (head, tail) = col.split_at_mut(e);
        let mut q = k;
        while q < e {
            let take = (e - q).min(4);
            let at = |o: usize| {
                if o < take {
                    &left[(q + o) * m + e..(q + o + 1) * m]
                } else {
                    &left[q * m + e..(q + 1) * m]
                }
            };
            let coef = |o: usize| if o < take { head[q + o] } else { ZERO };
            subtract(tail, [at(0), at(1), at(2), at(3)], [coef(0), coef(1), coef(2), coef(3)]);
            q += take;
        }
    }
}
