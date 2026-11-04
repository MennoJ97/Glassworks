//! Banded LU with partial pivoting, after reverse Cuthill-McKee reordering.
//! Cost is O(N b^2) for bandwidth b: ideal in 1D, adequate for small 2D meshes.

use crate::assemble::Csr;
use crate::{C64, ZERO};

/// Reverse Cuthill-McKee: breadth-first numbering keeps neighbours close, which
/// narrows the band. A periodic ring folds to bandwidth 2. `perm[new] = old`.
pub fn rcm(a: &Csr) -> Vec<usize> {
    let degree = |i: usize| a.ptr[i + 1] - a.ptr[i];
    let mut starts: Vec<usize> = (0..a.n).collect();
    starts.sort_by_key(|&i| degree(i));
    let (mut perm, mut seen) = (Vec::with_capacity(a.n), vec![false; a.n]);
    for start in starts {
        if seen[start] {
            continue;
        }
        seen[start] = true;
        // `perm` doubles as the breadth-first queue: `head` is the next node to visit.
        let mut head = perm.len();
        perm.push(start);
        while head < perm.len() {
            let node = perm[head];
            head += 1;
            let mut unseen: Vec<usize> =
                a.col[a.ptr[node]..a.ptr[node + 1]].iter().copied().filter(|&next| !seen[next]).collect();
            unseen.sort_by_key(|&next| degree(next));
            for next in unseen {
                seen[next] = true;
                perm.push(next);
            }
        }
    }
    perm.reverse();
    perm
}

/// Anything that can solve `A x = rhs` in place once `A` has been factored.
pub trait Factors {
    fn solve(&self, rhs: &mut [C64]);
}

/// LU factors in LAPACK band layout: entry (i, j) lives at offset `2b + i - j` of
/// column j, which leaves room for pivoting to grow the upper band from b to 2b.
pub struct BandedLu {
    n: usize,
    /// Bandwidth: how far from the diagonal the reordered matrix has entries.
    b: usize,
    /// Stored length of each column, `3b + 1`.
    col_len: usize,
    /// The columns of L and U, one after another.
    bands: Vec<C64>,
    /// Row swapped with row j when column j was eliminated.
    pivot_row: Vec<usize>,
    perm: Vec<usize>,
}

impl BandedLu {
    /// Reorders with `rcm`, copies the matrix into band storage, then eliminates column by
    /// column. Fails only on an exactly singular matrix.
    pub fn factor(a: &Csr) -> Result<Self, String> {
        let n = a.n;
        let perm = rcm(a);
        let mut new_index = vec![0; n];
        for (new, &old) in perm.iter().enumerate() {
            new_index[old] = new;
        }
        let mut b = 0;
        for i in 0..n {
            for &j in &a.col[a.ptr[i]..a.ptr[i + 1]] {
                b = b.max(new_index[i].abs_diff(new_index[j]));
            }
        }
        // `diag` is the row of each stored column that holds the diagonal.
        let (col_len, diag) = (3 * b + 1, 2 * b);
        let mut bands = vec![ZERO; col_len * n];
        for i in 0..n {
            for k in a.ptr[i]..a.ptr[i + 1] {
                let (row, col) = (new_index[i], new_index[a.col[k]]);
                bands[diag + row + col * col_len - col] = a.val[k];
            }
        }
        let mut pivot_row = vec![0; n];
        let mut last_col = 0;
        for j in 0..n {
            let below = b.min(n - 1 - j);
            let diagonal_at = j * col_len + diag;
            // Partial pivoting: the matrix is indefinite, so pick the largest entry.
            let pivot = (0..=below)
                .max_by(|&x, &y| {
                    bands[diagonal_at + x].norm_sqr().total_cmp(&bands[diagonal_at + y].norm_sqr())
                })
                .unwrap();
            if bands[diagonal_at + pivot] == ZERO {
                return Err(
                    "singular matrix: this frequency hits an exact resonance of a lossless structure".into(),
                );
            }
            pivot_row[j] = j + pivot;
            // The pivot row ends by last_col: the row swapped in from j + pivot reaches b
            // columns past its own diagonal, and earlier fill stops at the old last_col.
            last_col = last_col.max((j + pivot + b).min(n - 1));
            if pivot != 0 {
                for col in j..=last_col {
                    bands.swap(diag + j + col * col_len - col, diag + j + pivot + col * col_len - col);
                }
            }
            // Eliminate below the pivot, updating the columns to the right.
            let inv_pivot = bands[diagonal_at].inv();
            let (left, right) = bands.split_at_mut((j + 1) * col_len);
            let multipliers = &mut left[diagonal_at + 1..diagonal_at + 1 + below];
            multipliers.iter_mut().for_each(|x| *x *= inv_pivot);
            // a[r][c] -= multipliers[r] * a[j][c], for rows r below j and columns c up to last_col.
            for col in j + 1..=last_col {
                let start = (col - j - 1) * col_len + diag + j - col;
                let in_pivot_row = right[start];
                for (entry, m) in right[start + 1..start + 1 + below].iter_mut().zip(multipliers.iter()) {
                    *entry -= *m * in_pivot_row;
                }
            }
        }
        Ok(Self { n, b, col_len, bands, pivot_row, perm })
    }
}

impl Factors for BandedLu {
    fn solve(&self, rhs: &mut [C64]) {
        let (n, b, col_len, diag) = (self.n, self.b, self.col_len, 2 * self.b);
        let mut x: Vec<C64> = self.perm.iter().map(|&old| rhs[old]).collect();
        // Forward substitution with L, applying each row swap as factoring did.
        for j in 0..n {
            x.swap(j, self.pivot_row[j]);
            let (xj, below, diagonal_at) = (x[j], b.min(n - 1 - j), j * col_len + diag);
            for (xi, l) in x[j + 1..j + 1 + below].iter_mut().zip(&self.bands[diagonal_at + 1..]) {
                *xi -= *l * xj;
            }
        }
        // Back substitution with U, column by column.
        for j in (0..n).rev() {
            let (diagonal_at, above) = (j * col_len + diag, j.min(diag));
            x[j] /= self.bands[diagonal_at];
            let xj = x[j];
            for (xi, u) in x[j - above..j].iter_mut().zip(&self.bands[diagonal_at - above..diagonal_at]) {
                *xi -= *u * xj;
            }
        }
        for (value, &old) in x.iter().zip(&self.perm) {
            rhs[old] = *value;
        }
    }
}
