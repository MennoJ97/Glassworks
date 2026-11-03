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
    /// Unit-coefficient stiffness and mass matrices of a cell, row-major. Stiffness comes
    /// in x and y parts, because an absorbing layer weights them differently.
    fn local_matrices(&self, cell: usize, stiff_x: &mut [f64], stiff_y: &mut [f64], mass: &mut [f64]);
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

/// A degree-4 rule on triangles ([Dunavant 1985]): barycentric coordinates and a weight, summing to 1.
const TRIANGLE_RULE: [([f64; 3], f64); 6] = [
    ([0.108103018168070, 0.445948490915965, 0.445948490915965], 0.223381589678011),
    ([0.445948490915965, 0.108103018168070, 0.445948490915965], 0.223381589678011),
    ([0.445948490915965, 0.445948490915965, 0.108103018168070], 0.223381589678011),
    ([0.816847572980459, 0.091576213509771, 0.091576213509771], 0.109951743655322),
    ([0.091576213509771, 0.816847572980459, 0.091576213509771], 0.109951743655322),
    ([0.091576213509771, 0.091576213509771, 0.816847572980459], 0.109951743655322),
];

/// The six second-order basis functions, at barycentric coordinates.
fn basis_values([l0, l1, l2]: [f64; 3]) -> [f64; 6] {
    [
        l0 * (2.0 * l0 - 1.0),
        l1 * (2.0 * l1 - 1.0),
        l2 * (2.0 * l2 - 1.0),
        4.0 * l0 * l1,
        4.0 * l1 * l2,
        4.0 * l2 * l0,
    ]
}

/// One derivative of the six second-order basis functions, given that derivative of each
/// barycentric coordinate.
fn basis_gradient([l0, l1, l2]: [f64; 3], g: &[f64; 3]) -> [f64; 6] {
    [
        (4.0 * l0 - 1.0) * g[0],
        (4.0 * l1 - 1.0) * g[1],
        (4.0 * l2 - 1.0) * g[2],
        4.0 * (l0 * g[1] + l1 * g[0]),
        4.0 * (l1 * g[2] + l2 * g[1]),
        4.0 * (l2 * g[0] + l0 * g[2]),
    ]
}

/// Continuous second-order triangles, with dofs on vertices and then edge midpoints.
pub struct H1Triangle<'a> {
    pub mesh: &'a Mesh,
    n_vertices: usize,
    edges: Vec<[usize; 2]>,
    /// The edge opposite each pair: `(v0, v1)`, `(v1, v2)`, `(v2, v0)`.
    cell_edges: Vec<[usize; 3]>,
}

impl<'a> H1Triangle<'a> {
    pub fn new(mesh: &'a Mesh) -> Self {
        let edges = mesh.edges();
        let index: std::collections::HashMap<[usize; 2], usize> =
            edges.iter().enumerate().map(|(i, e)| (*e, i)).collect();
        let id = |a: usize, b: usize| index[&[a.min(b), a.max(b)]];
        let cell_edges =
            mesh.cells.chunks(3).map(|c| [id(c[0], c[1]), id(c[1], c[2]), id(c[2], c[0])]).collect();
        Self { mesh, n_vertices: mesh.coords.len() / 2, edges, cell_edges }
    }

    /// The six basis values at a point in a cell.
    pub fn basis_at(&self, cell: usize, x: f64, y: f64) -> [f64; 6] {
        basis_values(self.barycentric(cell, x, y))
    }

    /// Whether a point is in a cell, or on its edge.
    pub fn contains(&self, cell: usize, x: f64, y: f64) -> bool {
        self.barycentric(cell, x, y).iter().all(|&l| l >= -1e-9)
    }

    /// A cell's signed doubled area, and the x and y gradients of its barycentric
    /// coordinates, which are constant on a straight triangle.
    fn shape(&self, cell: usize) -> (f64, [f64; 3], [f64; 3]) {
        let v = &self.mesh.cells[3 * cell..3 * cell + 3];
        let [p0, p1, p2] = [self.mesh.vertex(v[0]), self.mesh.vertex(v[1]), self.mesh.vertex(v[2])];
        let det = (p1.0 - p0.0) * (p2.1 - p0.1) - (p2.0 - p0.0) * (p1.1 - p0.1);
        (
            det,
            [(p1.1 - p2.1) / det, (p2.1 - p0.1) / det, (p0.1 - p1.1) / det],
            [(p2.0 - p1.0) / det, (p0.0 - p2.0) / det, (p1.0 - p0.0) / det],
        )
    }

    /// Barycentric coordinates of a point, all in `[0, 1]` if it's inside the cell.
    pub fn barycentric(&self, cell: usize, x: f64, y: f64) -> [f64; 3] {
        let v = &self.mesh.cells[3 * cell..3 * cell + 3];
        let [p0, p1, p2] = [self.mesh.vertex(v[0]), self.mesh.vertex(v[1]), self.mesh.vertex(v[2])];
        let det = (p1.0 - p0.0) * (p2.1 - p0.1) - (p2.0 - p0.0) * (p1.1 - p0.1);
        let l1 = ((x - p0.0) * (p2.1 - p0.1) - (p2.0 - p0.0) * (y - p0.1)) / det;
        let l2 = ((p1.0 - p0.0) * (y - p0.1) - (x - p0.0) * (p1.1 - p0.1)) / det;
        [1.0 - l1 - l2, l1, l2]
    }

    /// The six basis functions' x and y derivatives at a point in a cell.
    pub fn gradients_at(&self, cell: usize, x: f64, y: f64) -> ([f64; 6], [f64; 6]) {
        let (_, gx, gy) = self.shape(cell);
        let l = self.barycentric(cell, x, y);
        (basis_gradient(l, &gx), basis_gradient(l, &gy))
    }

    pub fn dof_positions(&self) -> Vec<(f64, f64)> {
        let vertices = (0..self.n_vertices).map(|v| self.mesh.vertex(v));
        let midpoints = self.edges.iter().map(|&[a, b]| {
            let (p, q) = (self.mesh.vertex(a), self.mesh.vertex(b));
            ((p.0 + q.0) / 2.0, (p.1 + q.1) / 2.0)
        });
        vertices.chain(midpoints).collect()
    }
}

impl Space for H1Triangle<'_> {
    fn n_dofs(&self) -> usize {
        self.n_vertices + self.edges.len()
    }
    fn n_cells(&self) -> usize {
        self.mesh.n_cells()
    }
    fn local_size(&self) -> usize {
        6
    }
    fn cell_dofs(&self, cell: usize, dofs: &mut [usize]) {
        dofs[..3].copy_from_slice(&self.mesh.cells[3 * cell..3 * cell + 3]);
        for (d, e) in dofs[3..].iter_mut().zip(self.cell_edges[cell]) {
            *d = self.n_vertices + e;
        }
    }
    fn local_matrices(&self, cell: usize, stiff_x: &mut [f64], stiff_y: &mut [f64], mass: &mut [f64]) {
        let (det, gx, gy) = self.shape(cell);
        let area = det.abs() / 2.0;
        stiff_x.fill(0.0);
        stiff_y.fill(0.0);
        mass.fill(0.0);
        for ([l0, l1, l2], w) in TRIANGLE_RULE {
            let val = basis_values([l0, l1, l2]);
            let (dx, dy) = (basis_gradient([l0, l1, l2], &gx), basis_gradient([l0, l1, l2], &gy));
            for i in 0..6 {
                for j in 0..6 {
                    stiff_x[i * 6 + j] += w * area * dx[i] * dx[j];
                    stiff_y[i * 6 + j] += w * area * dy[i] * dy[j];
                    mass[i * 6 + j] += w * area * val[i] * val[j];
                }
            }
        }
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
    fn local_matrices(&self, cell: usize, stiff_x: &mut [f64], stiff_y: &mut [f64], mass: &mut [f64]) {
        // Mapping [-1, 1] onto a cell of length h scales d/dx by 2/h and dx by h/2.
        let h = self.cell_length(cell);
        stiff_y.fill(0.0);
        for (s, r) in stiff_x.iter_mut().zip(&self.ref_stiffness) {
            *s = r * 2.0 / h;
        }
        for (m, r) in mass.iter_mut().zip(&self.ref_mass) {
            *m = r * h / 2.0;
        }
    }
}
