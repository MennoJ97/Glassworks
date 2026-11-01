//! Sparse LU for 2D, where the band in `banded` fills in too much. `faer` reorders
//! the matrix to limit that fill. See `docs/` and `plans/2d-solver.md` for timings.

use crate::assemble::Csr;
use crate::banded::Factors;
use crate::C64;
use faer::prelude::*;
use faer::sparse::linalg::solvers::Lu;
use faer::sparse::{SparseColMat, Triplet};
use faer::{c64, Mat};

pub struct SparseLu {
    lu: Lu<usize, c64>,
    n: usize,
}

impl SparseLu {
    pub fn factor(a: &Csr) -> Result<Self, String> {
        let mut triplets = Vec::with_capacity(a.val.len());
        for i in 0..a.n {
            for k in a.ptr[i]..a.ptr[i + 1] {
                triplets.push(Triplet::new(i, a.col[k], c64::new(a.val[k].re, a.val[k].im)));
            }
        }
        let m = SparseColMat::<usize, c64>::try_new_from_triplets(a.n, a.n, &triplets).map_err(|e| format!("{e:?}"))?;
        let lu = m.sp_lu().map_err(|_| "singular matrix: this frequency hits an exact resonance".to_string())?;
        Ok(Self { lu, n: a.n })
    }
}

impl Factors for SparseLu {
    fn solve(&self, rhs: &mut [C64]) {
        let b = Mat::<c64>::from_fn(self.n, 1, |i, _| c64::new(rhs[i].re, rhs[i].im));
        let x = self.lu.solve(&b);
        for (i, r) in rhs.iter_mut().enumerate() {
            *r = C64::new(x[(i, 0)].re, x[(i, 0)].im);
        }
    }
}
