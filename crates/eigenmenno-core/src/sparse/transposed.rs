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
        for (node, front) in nodes.iter().zip(&self.fronts) {
            let (p, m) = (node.end - node.start, node.end - node.start + node.boundary.len());
            let (head, tail) = x.split_at_mut(node.end);
            let piv = &mut head[node.start..];
            for j in 0..p {
                let col = &front.lower[j * m..(j + 1) * m];
                let above: C64 = (0..j).map(|i| piv[i] * col[i]).sum();
                piv[j] = (piv[j] - above) / col[j];
            }
            for (o, &g) in node.boundary.iter().enumerate() {
                let col = &front.upper[o * p..(o + 1) * p];
                tail[g - node.end] -= piv.iter().zip(col).map(|(y, u)| y * u).sum::<C64>();
            }
        }
        for (node, front) in nodes.iter().zip(&self.fronts).rev() {
            let (p, m) = (node.end - node.start, node.end - node.start + node.boundary.len());
            let (head, tail) = x.split_at_mut(node.end);
            let piv = &mut head[node.start..];
            for j in (0..p).rev() {
                let col = &front.lower[j * m..(j + 1) * m];
                let below: C64 = (j + 1..p).map(|i| piv[i] * col[i]).sum();
                let beyond: C64 =
                    node.boundary.iter().enumerate().map(|(o, &g)| tail[g - node.end] * col[p + o]).sum();
                piv[j] -= below + beyond;
            }
            front.swaps.iter().enumerate().rev().for_each(|(j, &r)| piv.swap(j, r as usize));
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
