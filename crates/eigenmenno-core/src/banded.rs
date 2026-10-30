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
    for s in starts {
        if seen[s] {
            continue;
        }
        seen[s] = true;
        let mut head = perm.len();
        perm.push(s);
        while head < perm.len() {
            let v = perm[head];
            head += 1;
            let mut next: Vec<usize> = a.col[a.ptr[v]..a.ptr[v + 1]].iter().copied().filter(|&w| !seen[w]).collect();
            next.sort_by_key(|&w| degree(w));
            for w in next {
                seen[w] = true;
                perm.push(w);
            }
        }
    }
    perm.reverse();
    perm
}

/// LU factors in LAPACK band layout: entry (i, j) lives at offset `2b + i - j` of
/// column j, which leaves room for pivoting to grow the upper band from b to 2b.
pub struct BandedLu {
    n: usize,
    b: usize,
    ld: usize,
    ab: Vec<C64>,
    piv: Vec<usize>,
    perm: Vec<usize>,
}

impl BandedLu {
    pub fn factor(a: &Csr) -> Result<Self, String> {
        let n = a.n;
        let perm = rcm(a);
        let mut inv = vec![0; n];
        for (new, &old) in perm.iter().enumerate() {
            inv[old] = new;
        }
        let mut b = 0;
        for i in 0..n {
            for &j in &a.col[a.ptr[i]..a.ptr[i + 1]] {
                b = b.max(inv[i].abs_diff(inv[j]));
            }
        }
        let (ld, diag) = (3 * b + 1, 2 * b);
        let mut ab = vec![ZERO; ld * n];
        for i in 0..n {
            for k in a.ptr[i]..a.ptr[i + 1] {
                let (r, c) = (inv[i], inv[a.col[k]]);
                ab[diag + r + c * ld - c] = a.val[k];
            }
        }
        let mut piv = vec![0; n];
        let mut last_col = 0;
        for j in 0..n {
            let below = b.min(n - 1 - j);
            let col = j * ld + diag;
            // Partial pivoting: the matrix is indefinite, so pick the largest entry.
            let p = (0..=below).max_by(|&x, &y| ab[col + x].norm_sqr().total_cmp(&ab[col + y].norm_sqr())).unwrap();
            if ab[col + p] == ZERO {
                return Err("singular matrix: this frequency hits an exact resonance of a lossless structure".into());
            }
            piv[j] = j + p;
            last_col = last_col.max((j + p + b).min(n - 1));
            if p != 0 {
                for c in j..=last_col {
                    ab.swap(diag + j + c * ld - c, diag + j + p + c * ld - c);
                }
            }
            // Eliminate below the pivot, updating the columns to the right.
            let inv_pivot = ab[col].inv();
            let (left, right) = ab.split_at_mut((j + 1) * ld);
            let l = &mut left[col + 1..col + 1 + below];
            l.iter_mut().for_each(|x| *x *= inv_pivot);
            for c in j + 1..=last_col {
                let start = (c - j - 1) * ld + diag + j - c;
                let t = right[start];
                for (x, m) in right[start + 1..start + 1 + below].iter_mut().zip(l.iter()) {
                    *x -= *m * t;
                }
            }
        }
        Ok(Self { n, b, ld, ab, piv, perm })
    }

    /// Solves `A x = rhs` in place.
    pub fn solve(&self, rhs: &mut [C64]) {
        let (n, b, ld, diag) = (self.n, self.b, self.ld, 2 * self.b);
        let mut x: Vec<C64> = self.perm.iter().map(|&old| rhs[old]).collect();
        for j in 0..n {
            x.swap(j, self.piv[j]);
            let (t, below, col) = (x[j], b.min(n - 1 - j), j * ld + diag);
            for (xi, m) in x[j + 1..j + 1 + below].iter_mut().zip(&self.ab[col + 1..]) {
                *xi -= *m * t;
            }
        }
        for j in (0..n).rev() {
            let (col, above) = (j * ld + diag, j.min(diag));
            x[j] /= self.ab[col];
            let t = x[j];
            for (xi, u) in x[j - above..j].iter_mut().zip(&self.ab[col - above..col]) {
                *xi -= *u * t;
            }
        }
        for (w, &old) in x.iter().zip(&self.perm) {
            rhs[old] = *w;
        }
    }
}
