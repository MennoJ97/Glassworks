//! Solving `A^T x = b` with the factors of `A`, which an adjoint needs when `A` isn't symmetric,
//! as Bloch ties make it. Chapter 15 uses this for gratings.

use super::{Csr, SparseLu, C64, ZERO};

impl Csr {
    /// `y = A^T x`, so refinement can check a transposed solve against the matrix.
    pub fn mul_transposed(&self, x: &[C64], y: &mut [C64]) {
        y.fill(ZERO);
        for i in 0..self.n {
            for k in self.ptr[i]..self.ptr[i + 1] {
                y[self.col[k]] += self.val[k] * x[i];
            }
        }
    }
}

impl SparseLu<'_> {
    /// `solve_once` for `A^T`: the transpose of each of its steps, in the opposite order.
    fn solve_transposed_once(&self, rhs: &mut [C64]) {
        let nodes = &self.ordering.nodes;
        let mut x = vec![ZERO; rhs.len()];
        for (old, &new) in self.ordering.new.iter().enumerate() {
            x[new] = rhs[old];
        }
        // Forward through the fronts with `U^T`, which is lower triangular.
        for (node, front) in nodes.iter().zip(&self.fronts) {
            let n_pivots = node.end - node.start;
            let size = n_pivots + node.boundary.len();
            let (head, tail) = x.split_at_mut(node.end);
            let own = &mut head[node.start..];
            for j in 0..n_pivots {
                let col = &front.lower[j * size..(j + 1) * size];
                let above: C64 = (0..j).map(|i| own[i] * col[i]).sum();
                own[j] = (own[j] - above) / col[j];
            }
            for (offset, &unknown) in node.boundary.iter().enumerate() {
                let col = &front.upper[offset * n_pivots..(offset + 1) * n_pivots];
                tail[unknown - node.end] -= own.iter().zip(col).map(|(xi, u)| xi * u).sum::<C64>();
            }
        }
        // Back through them with `L^T`, then undo each front's row swaps, last swap first.
        for (node, front) in nodes.iter().zip(&self.fronts).rev() {
            let n_pivots = node.end - node.start;
            let size = n_pivots + node.boundary.len();
            let (head, tail) = x.split_at_mut(node.end);
            let own = &mut head[node.start..];
            for j in (0..n_pivots).rev() {
                let col = &front.lower[j * size..(j + 1) * size];
                let below: C64 = (j + 1..n_pivots).map(|i| own[i] * col[i]).sum();
                let beyond: C64 = node
                    .boundary
                    .iter()
                    .enumerate()
                    .map(|(offset, &unknown)| tail[unknown - node.end] * col[n_pivots + offset])
                    .sum();
                own[j] -= below + beyond;
            }
            front.swaps.iter().enumerate().rev().for_each(|(j, &row)| own.swap(j, row as usize));
        }
        for (old, &new) in self.ordering.new.iter().enumerate() {
            rhs[old] = x[new];
        }
    }

    /// Solves `A^T x = rhs`, for adjoints when Bloch ties make `A` unsymmetric.
    pub fn solve_transposed(&self, rhs: &mut [C64]) {
        self.refined(rhs, Self::solve_transposed_once, Csr::mul_transposed);
    }
}
