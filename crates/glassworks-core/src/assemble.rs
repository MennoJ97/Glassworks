//! Global sparse matrices and boundary constraints.

use crate::space::Space;
use crate::{C64, ONE, ZERO};

/// Matrix entries as `(row, col, value)`, in any order, with repeats to be summed.
pub type Triplets = Vec<(usize, usize, C64)>;

/// Compressed sparse row matrix: row `i`'s entries are `val[ptr[i]..ptr[i + 1]]`, in the
/// columns `col[ptr[i]..ptr[i + 1]]`, sorted.
pub struct Csr {
    pub n: usize,
    pub ptr: Vec<usize>,
    pub col: Vec<usize>,
    pub val: Vec<C64>,
}

impl Csr {
    /// Builds from unsorted triplets, summing duplicates.
    pub fn from_triplets(n: usize, triplets: Triplets) -> Self {
        // Rows are short, so dealing triplets into rows beats sorting them all at once.
        // Count each row's entries, and add up the counts to find where each row starts.
        let mut row_start = vec![0; n + 1];
        triplets.iter().for_each(|entry| row_start[entry.0 + 1] += 1);
        for i in 0..n {
            row_start[i + 1] += row_start[i];
        }
        let mut by_row = vec![(0, 0, ZERO); triplets.len()];
        let mut next_slot = row_start.clone();
        for entry in triplets {
            by_row[next_slot[entry.0]] = entry;
            next_slot[entry.0] += 1;
        }
        for i in 0..n {
            by_row[row_start[i]..row_start[i + 1]].sort_unstable_by_key(|entry| entry.1);
        }
        // Repeats are now next to each other, so summing them is one pass.
        let mut matrix = Csr { n, ptr: vec![0; n + 1], col: vec![], val: vec![] };
        let mut previous = None;
        for (i, j, value) in by_row {
            if previous == Some((i, j)) {
                *matrix.val.last_mut().unwrap() += value;
            } else {
                matrix.col.push(j);
                matrix.val.push(value);
                matrix.ptr[i + 1] += 1;
                previous = Some((i, j));
            }
        }
        for i in 0..n {
            matrix.ptr[i + 1] += matrix.ptr[i];
        }
        matrix
    }

    fn triplets(&self) -> impl Iterator<Item = (usize, usize, C64)> + '_ {
        (0..self.n)
            .flat_map(move |i| (self.ptr[i]..self.ptr[i + 1]).map(move |k| (i, self.col[k], self.val[k])))
    }

    /// `y = self x`.
    pub fn mul(&self, x: &[C64], y: &mut [C64]) {
        y.fill(ZERO);
        for (i, j, value) in self.triplets() {
            y[i] += value * x[j];
        }
    }

    /// `self - sigma * other`.
    pub fn sub_scaled(&self, sigma: f64, other: &Csr) -> Csr {
        let scaled = other.triplets().map(|(i, j, value)| (i, j, value * -sigma));
        Csr::from_triplets(self.n, self.triplets().chain(scaled).collect())
    }
}

/// Maps each global dof to `(reduced index, factor)`, or `None` if it is fixed
/// to zero. Dirichlet walls and Bloch periodicity are both just entries here.
pub struct DofMap {
    pub map: Vec<Option<(usize, C64)>>,
    /// The number of unknowns left to solve for: the size of every assembled matrix.
    pub n_free: usize,
}

impl DofMap {
    /// `zero` lists Dirichlet dofs; each tie `(slave, master, phase)` sets
    /// `u[slave] = phase * u[master]`, and no master may itself be a slave.
    pub fn new(n_dofs: usize, zero: &[usize], ties: &[(usize, usize, C64)]) -> Self {
        let mut dofs = DofMap { map: vec![None; n_dofs], n_free: 0 };
        let mut fixed = vec![false; n_dofs];
        zero.iter().chain(ties.iter().map(|tie| &tie.0)).for_each(|&dof| fixed[dof] = true);
        for dof in (0..n_dofs).filter(|&dof| !fixed[dof]) {
            dofs.map[dof] = Some((dofs.n_free, ONE));
            dofs.n_free += 1;
        }
        for &(slave, master, phase) in ties {
            dofs.map[slave] = dofs.map[master].map(|(index, factor)| (index, factor * phase));
        }
        dofs
    }

    /// Expands a reduced vector back onto all dofs.
    pub fn expand(&self, reduced: &[C64]) -> Vec<C64> {
        self.map.iter().map(|entry| entry.map_or(ZERO, |(index, factor)| reduced[index] * factor)).collect()
    }
}

/// Assembles `sum_c (a_c S_c + b_c M_c)`, where `coef` gives `(a_c, b_c)` per
/// cell. Almost every operator in this solver has this form.
pub fn assemble(space: &dyn Space, dofs: &DofMap, coef: impl Fn(usize) -> (C64, C64)) -> Triplets {
    assemble_split(space, dofs, |cell| {
        let (a, b) = coef(cell);
        (a, a, b)
    })
}

/// `assemble` with separate x and y stiffness weights: `(a_x, a_y, b)`.
pub fn assemble_split(space: &dyn Space, dofs: &DofMap, coef: impl Fn(usize) -> (C64, C64, C64)) -> Triplets {
    let n = space.local_size();
    let (mut stiff_x, mut stiff_y, mut mass) = (vec![0.0; n * n], vec![0.0; n * n], vec![0.0; n * n]);
    let mut cell_dofs = vec![0; n];
    let mut triplets = Vec::with_capacity(space.n_cells() * n * n);
    for cell in 0..space.n_cells() {
        space.cell_dofs(cell, &mut cell_dofs);
        space.local_matrices(cell, &mut stiff_x, &mut stiff_y, &mut mass);
        let (ax, ay, b) = coef(cell);
        for i in 0..n {
            let Some((row, row_factor)) = dofs.map[cell_dofs[i]] else { continue };
            for j in 0..n {
                let Some((col, col_factor)) = dofs.map[cell_dofs[j]] else { continue };
                let local = i * n + j;
                // A tied dof is `factor` times a free one. The test function is conjugated,
                // so the row picks up `conj(row_factor)` and the column `col_factor`.
                let entry = ax * stiff_x[local] + ay * stiff_y[local] + b * mass[local];
                triplets.push((row, col, row_factor.conj() * entry * col_factor));
            }
        }
    }
    triplets
}
