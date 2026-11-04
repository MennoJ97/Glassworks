//! The speed in `partial_lu`: it eliminates a block of columns, then updates every later
//! column from the whole block at once, four columns at a time, while the block is in cache.
//! The arithmetic is the same as updating after each column, so a first reading can skip this.
//! Chapter 16 has what it saved.

use super::C64;
use crate::ZERO;

/// Columns eliminated together before updating the rest of a front, to stay in cache.
pub(super) const BLOCK: usize = 32;

/// `y -= sum(scales[q] x[q])`, over four columns at once for speed: `y` is read and
/// written once instead of four times.
#[inline]
pub(super) fn subtract(y: &mut [C64], x: [&[C64]; 4], scales: [C64; 4]) {
    let n = y.len();
    let [x0, x1, x2, x3] = [&x[0][..n], &x[1][..n], &x[2][..n], &x[3][..n]];
    for i in 0..n {
        y[i] -= scales[0] * x0[i] + scales[1] * x1[i] + scales[2] * x2[i] + scales[3] * x3[i];
    }
}

/// Brings every column after `block_end` up to date with the eliminated columns
/// `block_start..block_end` of a column-major `front` with `size` rows: first the rows inside
/// the block, then the rest of each column, four eliminated columns at a time.
#[inline]
pub(super) fn update_after(front: &mut [C64], size: usize, block_start: usize, block_end: usize) {
    let (left, right) = front.split_at_mut(block_end * size);
    for col in right.chunks_exact_mut(size) {
        // Within the block, this column's entries become rows of U, one after another.
        for q in block_start..block_end {
            let in_pivot_row = col[q];
            for i in q + 1..block_end {
                col[i] -= in_pivot_row * left[q * size + i];
            }
        }
        // Below the block, each entry loses `U[q] L[i][q]` for every eliminated column q.
        let (head, tail) = col.split_at_mut(block_end);
        let mut q = block_start;
        while q < block_end {
            let take = (block_end - q).min(4);
            // A short last group repeats column q with a zero scale, so `subtract` stays four wide.
            let multipliers = |lane: usize| {
                let from = if lane < take { q + lane } else { q };
                &left[from * size + block_end..(from + 1) * size]
            };
            let scale = |lane: usize| if lane < take { head[q + lane] } else { ZERO };
            subtract(
                tail,
                [multipliers(0), multipliers(1), multipliers(2), multipliers(3)],
                [scale(0), scale(1), scale(2), scale(3)],
            );
            q += take;
        }
    }
}
