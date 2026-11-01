//! Sparse Cholesky for 2D, where the band in `banded` fills in too much. `faer`
//! reorders the matrix to limit that fill. See `plans/2d-solver.md` for timings.
//!
//! Shift-invert below every eigenvalue makes `S - sigma Q` Hermitian positive definite,
//! which needs no pivoting and stores only one triangle. It took about a third of the
//! time of general sparse LU over a band diagram. Driven problems are indefinite, so
//! they'll need LU.

use crate::assemble::Csr;
use crate::banded::Factors;
use crate::C64;
use faer::prelude::*;
use faer::sparse::linalg::solvers::{Llt, SymbolicLlt};
use faer::sparse::{SparseColMat, Triplet};
use faer::{c64, Mat, Side};

pub struct SparseCholesky {
    llt: Llt<usize, c64>,
    n: usize,
}

impl SparseCholesky {
    /// `analysis` holds the elimination order and fill pattern, which depend only on
    /// where the nonzeros are. It's computed on the first call, then reused for later
    /// matrices, which must have the same pattern, as Bloch phases do.
    pub fn factor(a: &Csr, analysis: &mut Option<SymbolicLlt<usize>>) -> Result<Self, String> {
        let mut triplets = Vec::with_capacity(a.val.len());
        for i in 0..a.n {
            for k in a.ptr[i]..a.ptr[i + 1] {
                triplets.push(Triplet::new(i, a.col[k], c64::new(a.val[k].re, a.val[k].im)));
            }
        }
        let m = SparseColMat::<usize, c64>::try_new_from_triplets(a.n, a.n, &triplets).map_err(|e| format!("{e:?}"))?;
        if analysis.is_none() {
            *analysis = Some(SymbolicLlt::try_new(m.symbolic(), Side::Lower).map_err(|e| format!("{e:?}"))?);
        }
        let llt = Llt::try_new_with_symbolic(analysis.clone().unwrap(), m.as_ref(), Side::Lower).map_err(|e| format!("{e:?}"))?;
        Ok(Self { llt, n: a.n })
    }
}

impl Factors for SparseCholesky {
    fn solve(&self, rhs: &mut [C64]) {
        let b = Mat::<c64>::from_fn(self.n, 1, |i, _| c64::new(rhs[i].re, rhs[i].im));
        let x = self.llt.solve(&b);
        for (i, r) in rhs.iter_mut().enumerate() {
            *r = C64::new(x[(i, 0)].re, x[(i, 0)].im);
        }
    }
}
