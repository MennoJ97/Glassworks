//! Global sparse matrices and boundary constraints.

use crate::space::Space;
use crate::{C64, ONE, ZERO};

pub type Triplets = Vec<(usize, usize, C64)>;

/// Compressed sparse row matrix.
pub struct Csr {
    pub n: usize,
    pub ptr: Vec<usize>,
    pub col: Vec<usize>,
    pub val: Vec<C64>,
}

impl Csr {
    /// Builds from unsorted triplets, summing duplicates.
    pub fn from_triplets(n: usize, t: Triplets) -> Self {
        // Rows are short, so dealing triplets into rows beats sorting them all at once.
        let mut start = vec![0; n + 1];
        t.iter().for_each(|e| start[e.0 + 1] += 1);
        for i in 0..n {
            start[i + 1] += start[i];
        }
        let mut sorted = vec![(0, 0, ZERO); t.len()];
        let mut next = start.clone();
        for e in t {
            sorted[next[e.0]] = e;
            next[e.0] += 1;
        }
        for i in 0..n {
            sorted[start[i]..start[i + 1]].sort_unstable_by_key(|e| e.1);
        }
        let mut m = Csr { n, ptr: vec![0; n + 1], col: vec![], val: vec![] };
        let mut last = None;
        for (i, j, v) in sorted {
            if last == Some((i, j)) {
                *m.val.last_mut().unwrap() += v;
            } else {
                m.col.push(j);
                m.val.push(v);
                m.ptr[i + 1] += 1;
                last = Some((i, j));
            }
        }
        for i in 0..n {
            m.ptr[i + 1] += m.ptr[i];
        }
        m
    }

    fn triplets(&self) -> impl Iterator<Item = (usize, usize, C64)> + '_ {
        (0..self.n)
            .flat_map(move |i| (self.ptr[i]..self.ptr[i + 1]).map(move |k| (i, self.col[k], self.val[k])))
    }

    pub fn mul(&self, x: &[C64], y: &mut [C64]) {
        y.fill(ZERO);
        for (i, j, v) in self.triplets() {
            y[i] += v * x[j];
        }
    }

    pub fn mul_transposed(&self, x: &[C64], y: &mut [C64]) {
        y.fill(ZERO);
        for (i, j, v) in self.triplets() {
            y[j] += v * x[i];
        }
    }

    /// `self - sigma * other`.
    pub fn sub_scaled(&self, sigma: f64, other: &Csr) -> Csr {
        let scaled = other.triplets().map(|(i, j, v)| (i, j, v * -sigma));
        Csr::from_triplets(self.n, self.triplets().chain(scaled).collect())
    }
}

/// Maps each global dof to `(reduced index, factor)`, or `None` if it is fixed
/// to zero. Dirichlet walls and Bloch periodicity are both just entries here.
pub struct DofMap {
    pub map: Vec<Option<(usize, C64)>>,
    pub n_free: usize,
}

impl DofMap {
    /// `zero` lists Dirichlet dofs; each tie `(slave, master, phase)` sets
    /// `u[slave] = phase * u[master]`, and no master may itself be a slave.
    pub fn new(n: usize, zero: &[usize], ties: &[(usize, usize, C64)]) -> Self {
        let mut d = DofMap { map: vec![None; n], n_free: 0 };
        let mut fixed = vec![false; n];
        zero.iter().chain(ties.iter().map(|t| &t.0)).for_each(|&i| fixed[i] = true);
        for i in (0..n).filter(|&i| !fixed[i]) {
            d.map[i] = Some((d.n_free, ONE));
            d.n_free += 1;
        }
        for &(slave, master, phase) in ties {
            d.map[slave] = d.map[master].map(|(i, f)| (i, f * phase));
        }
        d
    }

    /// Expands a reduced vector back onto all dofs.
    pub fn expand(&self, reduced: &[C64]) -> Vec<C64> {
        self.map.iter().map(|m| m.map_or(ZERO, |(i, f)| reduced[i] * f)).collect()
    }
}

/// Assembles `sum_c (a_c S_c + b_c M_c)`, where `coef` gives `(a_c, b_c)` per
/// cell. Almost every operator in this solver has this form.
pub fn assemble(space: &dyn Space, dofs: &DofMap, coef: impl Fn(usize) -> (C64, C64)) -> Triplets {
    assemble_split(space, dofs, |c| {
        let (a, b) = coef(c);
        (a, a, b)
    })
}

/// `assemble` with separate x and y stiffness weights: `(a_x, a_y, b)`.
pub fn assemble_split(space: &dyn Space, dofs: &DofMap, coef: impl Fn(usize) -> (C64, C64, C64)) -> Triplets {
    let n = space.local_size();
    let (mut sx, mut sy, mut m, mut d) = (vec![0.0; n * n], vec![0.0; n * n], vec![0.0; n * n], vec![0; n]);
    let mut t = Vec::with_capacity(space.n_cells() * n * n);
    for c in 0..space.n_cells() {
        space.cell_dofs(c, &mut d);
        space.local_matrices(c, &mut sx, &mut sy, &mut m);
        let (ax, ay, b) = coef(c);
        for i in 0..n {
            let Some((ri, fi)) = dofs.map[d[i]] else { continue };
            for j in 0..n {
                let Some((rj, fj)) = dofs.map[d[j]] else { continue };
                let k = i * n + j;
                t.push((ri, rj, fi.conj() * (ax * sx[k] + ay * sy[k] + b * m[k]) * fj));
            }
        }
    }
    t
}
