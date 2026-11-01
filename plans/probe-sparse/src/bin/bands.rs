//! Where a band diagram's time goes: factoring, solving, or the rest of Lanczos. Each solver
//! runs one diagram of a lattice of rods, with the shift below every eigenvalue.
use eigenmenno_core::assemble::{assemble, Csr, DofMap};
use eigenmenno_core::banded::Factors;
use eigenmenno_core::eigen::eigs_with;
use eigenmenno_core::mesh::{Mesh, Outline};
use eigenmenno_core::space::H1Triangle;
use eigenmenno_core::sparse::{Ordering, SparseLu};
use eigenmenno_core::{C64, ONE, ZERO};
use probe_sparse::faer_reference::SparseCholesky;
use std::cell::{Cell, RefCell};
use std::collections::HashMap;
use std::f64::consts::PI;
use std::time::Instant;

struct Timed<'a> {
    inner: &'a dyn Factors,
    time: &'a Cell<f64>,
    count: &'a Cell<usize>,
}
impl Factors for Timed<'_> {
    fn solve(&self, rhs: &mut [C64]) {
        let t = Instant::now();
        self.inner.solve(rhs);
        self.time.set(self.time.get() + t.elapsed().as_secs_f64());
        self.count.set(self.count.get() + 1);
    }
}

struct Setup {
    s: Csr,
    q: Csr,
    pos: Vec<(f64, f64)>,
}

fn setup(n: usize, phase: [f64; 2]) -> Setup {
    let rod = [Outline::Circle { x: 0.5, y: 0.5, r: 0.2 }];
    let mesh = Mesh::fitted_grid(n, n, 1.0, 1.0, &rod);
    let space = H1Triangle::new(&mesh);
    let all = space.dof_positions();
    let h = 1.0 / (2 * n) as f64;
    let grid = |v: f64| (v / h).round() as usize;
    let key = |p: (f64, f64)| (grid(p.0) % (2 * n), grid(p.1) % (2 * n));
    let (hx, hy) = (|p: (f64, f64)| grid(p.0) == 2 * n, |p: (f64, f64)| grid(p.1) == 2 * n);
    let masters: HashMap<_, _> = (0..all.len()).filter(|&d| !hx(all[d]) && !hy(all[d]) && (grid(all[d].0) == 0 || grid(all[d].1) == 0)).map(|d| (key(all[d]), d)).collect();
    let ties: Vec<_> = (0..all.len())
        .filter(|&d| hx(all[d]) || hy(all[d]))
        .map(|d| {
            let angle = if hx(all[d]) { phase[0] } else { 0.0 } + if hy(all[d]) { phase[1] } else { 0.0 };
            (d, masters[&key(all[d])], C64::from_polar(1.0, angle))
        })
        .collect();
    let dofs = DofMap::new(all.len(), &[], &ties);
    let region = &mesh.region;
    let s = Csr::from_triplets(dofs.n_free, assemble(&space, &dofs, |_| (ONE, ZERO)));
    let q = Csr::from_triplets(dofs.n_free, assemble(&space, &dofs, |c| (ZERO, C64::new(if region[c] == 0 { 1.0 } else { 12.0 }, 0.0))));
    let mut pos = vec![(f64::MAX, f64::MAX); dofs.n_free];
    for (d, m) in dofs.map.iter().enumerate() {
        if let Some((r, _)) = m {
            if all[d].0 + all[d].1 < pos[*r].0 + pos[*r].1 {
                pos[*r] = all[d];
            }
        }
    }
    Setup { s, q, pos }
}

fn phases(s: f64) -> [f64; 2] {
    if s <= 1.0 { [PI * s, 0.0] } else if s <= 2.0 { [PI, PI * (s - 1.0)] } else { [PI * (3.0 - s), PI * (3.0 - s)] }
}

fn main() {
    let sigma = -PI * PI;
    for n in [25, 50, 100] {
        let setups: Vec<Setup> = (0..16).map(|i| setup(n, phases(i as f64 / 5.0))).collect();
        let shifted: Vec<Csr> = setups.iter().map(|e| e.s.sub_scaled(sigma, &e.q)).collect();
        let dofs = setups[0].s.n;
        for who in ["faer Cholesky", "ours"] {
            let (solve, count, factor) = (Cell::new(0.0), Cell::new(0), RefCell::new(0.0));
            let start = Instant::now();
            let mut analysis = None;
            let ordering = Ordering::new(&shifted[0], &setups[0].pos);
            let mut low = 0.0;
            for (e, a) in setups.iter().zip(&shifted) {
                let t = Instant::now();
                let (f1, f2);
                let inner: &dyn Factors = if who == "ours" {
                    f1 = SparseLu::factor(a, &ordering);
                    &f1
                } else {
                    f2 = SparseCholesky::factor(a, &mut analysis).unwrap();
                    &f2
                };
                *factor.borrow_mut() += t.elapsed().as_secs_f64();
                let timed = Timed { inner, time: &solve, count: &count };
                low += eigs_with(&timed as &dyn Factors, &e.q, sigma, 6, 1e-10).unwrap()[0].0;
            }
            let total = start.elapsed().as_secs_f64();
            let (fa, so) = (*factor.borrow(), solve.get());
            println!("{dofs:>6} {who:14}: total {:7.0} ms, factor {:6.0}, solves {:6.0} ({} at {:.2} ms), rest {:6.0}, check {low:.6}", total * 1e3, fa * 1e3, so * 1e3, count.get(), so * 1e3 / count.get() as f64, (total - fa - so) * 1e3);
        }
    }
}
