//! Throwaway 2D probe. P2 triangles on a structured grid, reusing eigenmenno-core.
use eigenmenno_core::assemble::{assemble, Csr, DofMap, Triplets};
use eigenmenno_core::banded::BandedLu;
use eigenmenno_core::mesh::Mesh;
use eigenmenno_core::space::Space;
use eigenmenno_core::{C64, ONE, ZERO};
use std::collections::HashMap;

/// Structured triangles; `coords` holds x, y per vertex, `cells` three ids per cell.
pub fn grid(nx: usize, ny: usize, lx: f64, ly: f64, x0: f64, y0: f64) -> Mesh {
    let mut m = Mesh { verts_per_cell: 3, coords: vec![], cells: vec![], region: vec![] };
    for j in 0..=ny { for i in 0..=nx {
        m.coords.extend([x0 + lx * i as f64 / nx as f64, y0 + ly * j as f64 / ny as f64]);
    } }
    let id = |i: usize, j: usize| i + j * (nx + 1);
    for j in 0..ny { for i in 0..nx {
        let (a, b, c, d) = (id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1));
        // Alternate the diagonal so the mesh has no preferred direction.
        if (i + j) % 2 == 0 { m.cells.extend([a, b, c, a, c, d]); } else { m.cells.extend([a, b, d, b, c, d]); }
        m.region.extend([0, 0]);
    } }
    m
}

fn edges_of(m: &Mesh) -> Vec<[usize; 2]> {
    let mut seen = std::collections::HashSet::new();
    let mut out = vec![];
    for c in m.cells.chunks(3) { for (a, b) in [(c[0], c[1]), (c[1], c[2]), (c[2], c[0])] {
        let e = [a.min(b), a.max(b)];
        if seen.insert(e) { out.push(e); }
    } }
    out
}

/// Moves the closer end of every edge crossing the level set onto it, so no edge
/// crosses any more and every cell lies on one side.
pub fn snap(m: &mut Mesh, phi: &dyn Fn(f64, f64) -> f64, project: &dyn Fn(f64, f64) -> (f64, f64)) -> usize {
    let mut moved = HashMap::new();
    for [a, b] in edges_of(m) {
        let f = |v: usize| phi(m.coords[2 * v], m.coords[2 * v + 1]);
        let (fa, fb) = (f(a), f(b));
        if fa * fb < 0.0 {
            moved.entry(if fa.abs() < fb.abs() { a } else { b }).or_insert(());
        }
    }
    for &v in moved.keys() {
        let (x, y) = project(m.coords[2 * v], m.coords[2 * v + 1]);
        m.coords[2 * v] = x;
        m.coords[2 * v + 1] = y;
    }
    moved.len()
}

pub fn classify(m: &mut Mesh, inside: &dyn Fn(f64, f64) -> usize) {
    for c in 0..m.region.len() {
        let v = &m.cells[3 * c..3 * c + 3];
        let cx = v.iter().map(|&i| m.coords[2 * i]).sum::<f64>() / 3.0;
        let cy = v.iter().map(|&i| m.coords[2 * i + 1]).sum::<f64>() / 3.0;
        m.region[c] = inside(cx, cy);
    }
}

pub fn min_angle(m: &Mesh) -> f64 {
    let mut worst = 180f64;
    for c in m.cells.chunks(3) {
        let p: Vec<(f64, f64)> = c.iter().map(|&i| (m.coords[2 * i], m.coords[2 * i + 1])).collect();
        for k in 0..3 {
            let (a, b, o) = (p[(k + 1) % 3], p[(k + 2) % 3], p[k]);
            let (u, v) = ((a.0 - o.0, a.1 - o.1), (b.0 - o.0, b.1 - o.1));
            let ang = ((u.0 * v.0 + u.1 * v.1) / (u.0.hypot(u.1) * v.0.hypot(v.1))).clamp(-1.0, 1.0).acos();
            worst = worst.min(ang.to_degrees());
        }
    }
    worst
}

/// Second-order Lagrange triangles with vertex and edge-midpoint dofs.
pub struct H1Tri<'a> { pub mesh: &'a Mesh, pub edge: HashMap<[usize; 2], usize>, pub nv: usize }

// Degree-4 rule (Dunavant, 6 points): (l0, l1, l2, weight summing to 1).
const Q: [(f64, f64, f64, f64); 6] = [
    (0.108103018168070, 0.445948490915965, 0.445948490915965, 0.223381589678011),
    (0.445948490915965, 0.108103018168070, 0.445948490915965, 0.223381589678011),
    (0.445948490915965, 0.445948490915965, 0.108103018168070, 0.223381589678011),
    (0.816847572980459, 0.091576213509771, 0.091576213509771, 0.109951743655322),
    (0.091576213509771, 0.816847572980459, 0.091576213509771, 0.109951743655322),
    (0.091576213509771, 0.091576213509771, 0.816847572980459, 0.109951743655322),
];

impl<'a> H1Tri<'a> {
    pub fn new(mesh: &'a Mesh) -> Self {
        let edge = edges_of(mesh).into_iter().enumerate().map(|(i, e)| (e, i)).collect();
        Self { mesh, edge, nv: mesh.coords.len() / 2 }
    }
    /// Unit-coefficient x and y stiffness and mass: anisotropy needs them apart.
    pub fn local(&self, c: usize, sxx: &mut [f64], syy: &mut [f64], mass: &mut [f64]) {
        let v = &self.mesh.cells[3 * c..3 * c + 3];
        let (x, y) = (|i: usize| self.mesh.coords[2 * v[i]], |i: usize| self.mesh.coords[2 * v[i] + 1]);
        let det = (x(1) - x(0)) * (y(2) - y(0)) - (x(2) - x(0)) * (y(1) - y(0));
        let gx = [(y(1) - y(2)) / det, (y(2) - y(0)) / det, (y(0) - y(1)) / det];
        let gy = [(x(2) - x(1)) / det, (x(0) - x(2)) / det, (x(1) - x(0)) / det];
        let area = det.abs() / 2.0;
        sxx.fill(0.0); syy.fill(0.0); mass.fill(0.0);
        for &(l0, l1, l2, w) in &Q {
            let l = [l0, l1, l2];
            let n = [l0 * (2.0 * l0 - 1.0), l1 * (2.0 * l1 - 1.0), l2 * (2.0 * l2 - 1.0), 4.0 * l0 * l1, 4.0 * l1 * l2, 4.0 * l2 * l0];
            let grad = |g: &[f64; 3]| [
                (4.0 * l[0] - 1.0) * g[0], (4.0 * l[1] - 1.0) * g[1], (4.0 * l[2] - 1.0) * g[2],
                4.0 * (l[0] * g[1] + l[1] * g[0]), 4.0 * (l[1] * g[2] + l[2] * g[1]), 4.0 * (l[2] * g[0] + l[0] * g[2]),
            ];
            let (dx, dy) = (grad(&gx), grad(&gy));
            for i in 0..6 { for j in 0..6 {
                sxx[i * 6 + j] += w * area * dx[i] * dx[j];
                syy[i * 6 + j] += w * area * dy[i] * dy[j];
                mass[i * 6 + j] += w * area * n[i] * n[j];
            } }
        }
    }
    pub fn dof_position(&self) -> Vec<(f64, f64)> {
        let c = &self.mesh.coords;
        let mut out = vec![(0.0, 0.0); self.n_dofs()];
        for v in 0..self.nv { out[v] = (c[2 * v], c[2 * v + 1]); }
        for (e, &i) in &self.edge { out[self.nv + i] = ((c[2 * e[0]] + c[2 * e[1]]) / 2.0, (c[2 * e[0] + 1] + c[2 * e[1] + 1]) / 2.0); }
        out
    }
}

impl Space for H1Tri<'_> {
    fn n_dofs(&self) -> usize { self.nv + self.edge.len() }
    fn n_cells(&self) -> usize { self.mesh.n_cells() }
    fn local_size(&self) -> usize { 6 }
    fn cell_dofs(&self, cell: usize, d: &mut [usize]) {
        let v = &self.mesh.cells[3 * cell..3 * cell + 3];
        let e = |a: usize, b: usize| self.nv + self.edge[&[a.min(b), a.max(b)]];
        d[0] = v[0]; d[1] = v[1]; d[2] = v[2];
        d[3] = e(v[0], v[1]); d[4] = e(v[1], v[2]); d[5] = e(v[2], v[0]);
    }
    fn local_matrices(&self, c: usize, s: &mut [f64], m: &mut [f64]) {
        let mut y = vec![0.0; 36];
        self.local(c, s, &mut y, m);
        s.iter_mut().zip(&y).for_each(|(a, b)| *a += b);
    }
}

/// Bloch ties for a unit cell, keyed on the unsnapped half-grid position of each dof.
pub fn bloch_map(space: &H1Tri, unsnapped: &Mesh, nx: usize, ny: usize, a: f64, phase: (f64, f64)) -> DofMap {
    let tmp = H1Tri { mesh: unsnapped, edge: space.edge.clone(), nv: space.nv };
    let pos = tmp.dof_position();
    let key = |(x, y): (f64, f64)| ((x / a * 2.0 * nx as f64).round() as usize, (y / a * 2.0 * ny as f64).round() as usize);
    let mut masters: HashMap<(usize, usize), usize> = HashMap::new();
    let mut map = vec![None; pos.len()];
    let mut n_free = 0;
    for (d, &p) in pos.iter().enumerate() {
        let (i, j) = key(p);
        let (wi, wj) = (i == 2 * nx, j == 2 * ny);
        if !wi && !wj { masters.insert((i, j), n_free); map[d] = Some((n_free, ONE)); n_free += 1; }
    }
    for (d, &p) in pos.iter().enumerate() {
        let (i, j) = key(p);
        let (wi, wj) = (i == 2 * nx, j == 2 * ny);
        if wi || wj {
            let m = masters[&(if wi { 0 } else { i }, if wj { 0 } else { j })];
            let ph = (if wi { phase.0 } else { 0.0 }) + (if wj { phase.1 } else { 0.0 });
            map[d] = Some((m, C64::from_polar(1.0, ph)));
        }
    }
    DofMap { map, n_free }
}

/// `-div(Lambda grad u) - k^2 s_x s_y u = f` with a quadratic-graded PML.
pub fn helmholtz_pml(space: &H1Tri, dofs: &DofMap, half: f64, t: f64, k: f64) -> Csr {
    let sigma_max = -3.0 * (1e-6f64).ln() / (2.0 * t);
    let stretch = |d: f64| {
        let depth = (d.abs() - half).max(0.0) / t;
        C64::new(1.0, sigma_max * depth * depth / k)
    };
    let (mut xx, mut yy, mut m, mut d) = (vec![0.0; 36], vec![0.0; 36], vec![0.0; 36], vec![0; 6]);
    let mut tr: Triplets = Vec::with_capacity(space.n_cells() * 36);
    for c in 0..space.n_cells() {
        space.cell_dofs(c, &mut d);
        space.local(c, &mut xx, &mut yy, &mut m);
        let v = &space.mesh.cells[3 * c..3 * c + 3];
        let cx = v.iter().map(|&i| space.mesh.coords[2 * i]).sum::<f64>() / 3.0;
        let cy = v.iter().map(|&i| space.mesh.coords[2 * i + 1]).sum::<f64>() / 3.0;
        let (sx, sy) = (stretch(cx), stretch(cy));
        let (ax, ay, b) = (sy / sx, sx / sy, -sx * sy * (k * k));
        for i in 0..6 { let Some((ri, _)) = dofs.map[d[i]] else { continue };
            for j in 0..6 { let Some((rj, _)) = dofs.map[d[j]] else { continue };
                tr.push((ri, rj, ax * xx[i * 6 + j] + ay * yy[i * 6 + j] + b * m[i * 6 + j]));
            } }
    }
    Csr::from_triplets(dofs.n_free, tr)
}

/// A Dirichlet-walled PML box of about `n` free dofs, for timing.
pub fn scaling_matrix(nx: usize) -> Csr {
    let mesh = grid(nx, nx, 4.0, 4.0, -2.0, -2.0);
    let space = H1Tri::new(&mesh);
    let pos = space.dof_position();
    let zero: Vec<usize> = pos.iter().enumerate().filter(|(_, p)| p.0.abs() >= 2.0 - 1e-9 || p.1.abs() >= 2.0 - 1e-9).map(|(i, _)| i).collect();
    let dofs = DofMap::new(space.n_dofs(), &zero, None);
    helmholtz_pml(&space, &dofs, 1.5, 0.5, 2.0 * std::f64::consts::PI)
}

pub fn rhs(n: usize) -> Vec<C64> { let mut b = vec![ZERO; n]; b[n / 2] = ONE; b }

pub fn solve_banded(a: &Csr) -> Vec<C64> {
    let mut x = rhs(a.n);
    BandedLu::factor(a).unwrap().solve(&mut x);
    x
}

#[cfg(feature = "sparse")]
pub fn solve_faer(a: &Csr) -> Vec<C64> {
    use faer::prelude::*;
    use faer::sparse::{SparseColMat, Triplet};
    let mut tr = Vec::with_capacity(a.val.len());
    for i in 0..a.n { for k in a.ptr[i]..a.ptr[i + 1] {
        tr.push(Triplet::new(i, a.col[k], faer::c64::new(a.val[k].re, a.val[k].im)));
    } }
    let m = SparseColMat::<usize, faer::c64>::try_new_from_triplets(a.n, a.n, &tr).unwrap();
    let b = faer::Mat::<faer::c64>::from_fn(a.n, 1, |i, _| if i == a.n / 2 { faer::c64::new(1.0, 0.0) } else { faer::c64::new(0.0, 0.0) });
    let y = m.sp_lu().unwrap().solve(&b);
    (0..a.n).map(|i| C64::new(y[(i, 0)].re, y[(i, 0)].im)).collect()
}

// Raw exports for timing in a WebAssembly engine.
#[no_mangle] pub extern "C" fn wasm_banded(nx: usize) -> f64 { let a = scaling_matrix(nx); solve_banded(&a)[a.n / 2].norm() }
#[cfg(feature = "sparse")]
#[no_mangle] pub extern "C" fn wasm_faer(nx: usize) -> f64 { let a = scaling_matrix(nx); solve_faer(&a)[a.n / 2].norm() }
#[no_mangle] pub extern "C" fn wasm_assemble(nx: usize) -> f64 { scaling_matrix(nx).n as f64 }
