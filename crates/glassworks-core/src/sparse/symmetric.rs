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
pub(super) fn partial_ldl(front: &mut [C64], size: usize, n_pivots: usize, scale: f64) -> Option<bool> {
    let mut weak = false;
    let mut block_start = 0;
    while block_start < n_pivots {
        let block_end = (block_start + BLOCK).min(n_pivots);
        for j in block_start..block_end {
            // Without row swaps, the pivot is whatever sits on the diagonal.
            let pivot = front[j + j * size];
            let below = front[j * size + j + 1..(j + 1) * size]
                .iter()
                .fold(0.0f64, |largest, entry| largest.max(entry.norm_sqr()));
            // A pivot far smaller than the column below it, or than the matrix, gives factors
            // too inaccurate for refinement to repair: `factor_symmetric` falls back to `factor`.
            if pivot.norm_sqr() < 1e-16 * below || pivot.norm_sqr() < 1e-28 * scale * scale {
                return None;
            }
            weak |= pivot.norm_sqr() < 1e-4 * below || pivot.norm_sqr() < 1e-16 * scale * scale;
            let inv_pivot = pivot.inv();
            front[j * size + j + 1..(j + 1) * size].iter_mut().for_each(|entry| *entry *= inv_pivot);
            let (left, right) = front.split_at_mut((j + 1) * size);
            // By symmetry, `U[j][column] = D[j] L[column][j]`, so U needn't be stored while factoring.
            for (column, col) in (j + 1..block_end).zip(right.chunks_exact_mut(size)) {
                let in_pivot_row = pivot * left[j * size + column];
                col[column..]
                    .iter_mut()
                    .zip(&left[j * size + column..])
                    .for_each(|(entry, multiplier)| *entry -= in_pivot_row * *multiplier);
            }
        }
        // Update the later columns from this block, on and below the diagonal only.
        let (left, right) = front.split_at_mut(block_end * size);
        for (column, col) in (block_end..size).zip(right.chunks_exact_mut(size)) {
            let mut q = block_start;
            while q < block_end {
                let take = (block_end - q).min(4);
                let multipliers = |lane: usize| {
                    let from = if lane < take { q + lane } else { q };
                    &left[from * size + column..(from + 1) * size]
                };
                let in_pivot_row = |lane: usize| {
                    if lane < take {
                        left[(q + lane) * size + q + lane] * left[(q + lane) * size + column]
                    } else {
                        ZERO
                    }
                };
                subtract(
                    &mut col[column..],
                    [multipliers(0), multipliers(1), multipliers(2), multipliers(3)],
                    [in_pivot_row(0), in_pivot_row(1), in_pivot_row(2), in_pivot_row(3)],
                );
                q += take;
            }
        }
        block_start = block_end;
    }
    // Fill in `U = D L^T` above the diagonal, where the solves shared with `factor` read it.
    for column in 1..size {
        for i in 0..column.min(n_pivots) {
            front[i + column * size] = front[i + i * size] * front[column + i * size];
        }
    }
    Some(weak)
}
