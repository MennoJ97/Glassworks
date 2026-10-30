//! Finite element spaces. Assembly only sees the `Space` trait, so 2D Lagrange
//! and Nedelec elements can be added as new implementations.

use crate::mesh::Mesh;
use crate::{C64, ZERO};

pub trait Space {
    fn n_dofs(&self) -> usize;
    fn n_cells(&self) -> usize;
    /// Number of dofs on each cell.
    fn local_size(&self) -> usize;
    fn cell_dofs(&self, cell: usize, dofs: &mut [usize]);
    /// Unit-coefficient stiffness and mass matrices of a cell, row-major.
    fn local_matrices(&self, cell: usize, stiffness: &mut [f64], mass: &mut [f64]);
}

/// Gauss-Legendre points and weights on [-1, 1]: Newton's method on `P_n`.
fn gauss_legendre(n: usize) -> Vec<(f64, f64)> {
    let legendre = |x: f64| {
        let (mut p0, mut p1) = (1.0, x);
        for k in 2..=n {
            let k = k as f64;
            (p0, p1) = (p1, ((2.0 * k - 1.0) * x * p1 - (k - 1.0) * p0) / k);
        }
        (p1, n as f64 * (x * p1 - p0) / (x * x - 1.0))
    };
    (0..n)
        .map(|i| {
            let mut x = -(std::f64::consts::PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
            for _ in 0..50 {
                let (p, d) = legendre(x);
                x -= p / d;
            }
            let d = legendre(x).1;
            (x, 2.0 / ((1.0 - x * x) * d * d))
        })
        .collect()
}

/// Continuous Lagrange elements of order p on a line. Neighbouring cells share
/// their end dof, which is what makes the field continuous.
pub struct H1Line<'a> {
    pub mesh: &'a Mesh,
    pub order: usize,
    ref_stiffness: Vec<f64>,
    ref_mass: Vec<f64>,
}

impl<'a> H1Line<'a> {
    pub fn new(mesh: &'a Mesh, order: usize) -> Self {
        let n = order + 1;
        let mut s = Self { mesh, order, ref_stiffness: vec![0.0; n * n], ref_mass: vec![0.0; n * n] };
        for (x, w) in gauss_legendre(n) {
            let (val, der) = s.basis(x);
            for i in 0..n {
                for j in 0..n {
                    s.ref_stiffness[i * n + j] += w * der[i] * der[j];
                    s.ref_mass[i * n + j] += w * val[i] * val[j];
                }
            }
        }
        s
    }

    /// Lagrange basis values and derivatives at `xi` in [-1, 1], on evenly
    /// spaced nodes: basis `i` is 1 at node `i` and 0 at all the others.
    pub fn basis(&self, xi: f64) -> (Vec<f64>, Vec<f64>) {
        let n = self.order + 1;
        let node = |i: usize| -1.0 + 2.0 * i as f64 / self.order as f64;
        let (mut val, mut der) = (vec![1.0; n], vec![0.0; n]);
        for i in 0..n {
            for j in (0..n).filter(|&j| j != i) {
                let w = 1.0 / (node(i) - node(j));
                der[i] = der[i] * (xi - node(j)) * w + val[i] * w;
                val[i] *= (xi - node(j)) * w;
            }
        }
        (val, der)
    }

    pub fn cell_length(&self, cell: usize) -> f64 {
        self.mesh.coords[self.mesh.cells[2 * cell + 1]] - self.mesh.coords[self.mesh.cells[2 * cell]]
    }

    /// Solution value and x-derivative at `xi` in a cell.
    pub fn eval(&self, cell: usize, xi: f64, u: &[C64]) -> (C64, C64) {
        let (val, der) = self.basis(xi);
        let (mut v, mut d) = (ZERO, ZERO);
        for i in 0..=self.order {
            v += u[cell * self.order + i] * val[i];
            d += u[cell * self.order + i] * der[i] * 2.0 / self.cell_length(cell);
        }
        (v, d)
    }
}

impl Space for H1Line<'_> {
    fn n_dofs(&self) -> usize {
        self.n_cells() * self.order + 1
    }
    fn n_cells(&self) -> usize {
        self.mesh.n_cells()
    }
    fn local_size(&self) -> usize {
        self.order + 1
    }
    fn cell_dofs(&self, cell: usize, dofs: &mut [usize]) {
        for (i, d) in dofs.iter_mut().enumerate() {
            *d = cell * self.order + i;
        }
    }
    fn local_matrices(&self, cell: usize, stiffness: &mut [f64], mass: &mut [f64]) {
        // Mapping [-1, 1] onto a cell of length h scales d/dx by 2/h and dx by h/2.
        let h = self.cell_length(cell);
        for (s, r) in stiffness.iter_mut().zip(&self.ref_stiffness) {
            *s = r * 2.0 / h;
        }
        for (m, r) in mass.iter_mut().zip(&self.ref_mass) {
            *m = r * h / 2.0;
        }
    }
}
