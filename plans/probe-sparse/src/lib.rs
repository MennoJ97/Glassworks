//! Throwaway benchmark of `multifrontal` against `faer`, on the solver's real 2D matrices.
//! The matrix builders repeat `solve2d.rs`, whose own are private.
pub mod faer_reference;
pub mod multifrontal;

use glassworks_core::assemble::{assemble, assemble_split, Csr, DofMap};
use glassworks_core::mesh::{Mesh, Outline};
use glassworks_core::space::H1Triangle;
use glassworks_core::{C64, I, ONE, ZERO};
use std::alloc::{GlobalAlloc, Layout, System};
use std::collections::HashMap;
use std::f64::consts::PI;
use std::sync::atomic::{AtomicUsize, Ordering::Relaxed};

struct Counting;
static NOW: AtomicUsize = AtomicUsize::new(0);
static PEAK: AtomicUsize = AtomicUsize::new(0);
unsafe impl GlobalAlloc for Counting {
    unsafe fn alloc(&self, l: Layout) -> *mut u8 {
        let now = NOW.fetch_add(l.size(), Relaxed) + l.size();
        PEAK.fetch_max(now, Relaxed);
        System.alloc(l)
    }
    unsafe fn dealloc(&self, p: *mut u8, l: Layout) {
        NOW.fetch_sub(l.size(), Relaxed);
        System.dealloc(p, l)
    }
}
#[global_allocator]
static COUNTER: Counting = Counting;
/// Restarts the peak, and returns how much memory is live.
pub fn reset_peak() -> usize {
    let now = NOW.load(Relaxed);
    PEAK.store(now, Relaxed);
    now
}
pub fn peak_mb(base: usize) -> f64 {
    (PEAK.load(Relaxed) - base) as f64 / 1e6
}

pub struct Case {
    pub name: &'static str,
    pub a: Csr,
    pub pos: Vec<(f64, f64)>,
    pub definite: bool,
}

fn circle(x: f64, y: f64, r: f64) -> Outline {
    Outline::Circle { x, y, r }
}

/// The matrix of `open_matrix` at complex `k`, for rods of `eps`, one region per material.
fn open(outlines: &[Outline], eps: f64, size: f64, cell: f64, pml: f64, k: C64) -> (Csr, Vec<(f64, f64)>) {
    let total = size + 2.0 * pml;
    let n = (total / cell).ceil() as usize;
    let shifted: Vec<_> = outlines.iter().map(|o| o.shifted(pml, pml)).collect();
    let mesh = Mesh::fitted_grid(n, n, total, total, &shifted);
    let space = H1Triangle::new(&mesh);
    let all = space.dof_positions();
    let edge = |p: (f64, f64)| p.0 < 1e-9 || p.1 < 1e-9 || p.0 > total - 1e-9 || p.1 > total - 1e-9;
    let zero: Vec<usize> = (0..all.len()).filter(|&d| edge(all[d])).collect();
    let dofs = DofMap::new(all.len(), &zero, &[]);
    let sigma_max = -3.0 * (1e-4f64).ln() / (2.0 * pml);
    let stretch = |d: f64| {
        let depth = (pml - d).max(d - pml - size).max(0.0) / pml;
        C64::new(1.0, 0.0) + I * (sigma_max * depth * depth) / k
    };
    let t = assemble_split(&space, &dofs, |c| {
        let v = &mesh.cells[3 * c..3 * c + 3];
        let (x, y) = v.iter().fold((0.0, 0.0), |s, &i| (s.0 + mesh.vertex(i).0 / 3.0, s.1 + mesh.vertex(i).1 / 3.0));
        let (sx, sy) = (stretch(x), stretch(y));
        let q = if mesh.region[c] == 0 { 1.0 } else { eps };
        (sy / sx, sx / sy, -k * k * q * sx * sy)
    });
    let mut pos = vec![(0.0, 0.0); dofs.n_free];
    for (d, m) in dofs.map.iter().enumerate() {
        if let Some((r, _)) = m {
            pos[*r] = all[d];
        }
    }
    (Csr::from_triplets(dofs.n_free, t), pos)
}

/// `S - shift Q` with Bloch ties on a periodic square of rods, as `modes` builds it.
fn periodic(outlines: &[Outline], eps: f64, size: f64, cell: f64, phase: [f64; 2], shift: f64) -> (Csr, Vec<(f64, f64)>) {
    let n = (size / cell).round() as usize;
    let mesh = Mesh::fitted_grid(n, n, size, size, outlines);
    let space = H1Triangle::new(&mesh);
    let all = space.dof_positions();
    let h = size / (2 * n) as f64;
    let grid = |v: f64| (v / h).round() as usize;
    let key = |p: (f64, f64)| (grid(p.0) % (2 * n), grid(p.1) % (2 * n));
    let (high_x, high_y) = (|p: (f64, f64)| grid(p.0) == 2 * n, |p: (f64, f64)| grid(p.1) == 2 * n);
    let masters: HashMap<_, _> = (0..all.len()).filter(|&d| !high_x(all[d]) && !high_y(all[d]) && (grid(all[d].0) == 0 || grid(all[d].1) == 0)).map(|d| (key(all[d]), d)).collect();
    let ties: Vec<_> = (0..all.len())
        .filter(|&d| high_x(all[d]) || high_y(all[d]))
        .map(|d| {
            let angle = if high_x(all[d]) { phase[0] } else { 0.0 } + if high_y(all[d]) { phase[1] } else { 0.0 };
            (d, masters[&key(all[d])], C64::from_polar(1.0, angle))
        })
        .collect();
    let dofs = DofMap::new(all.len(), &[], &ties);
    let region = &mesh.region;
    let t = assemble(&space, &dofs, |c| (ONE, C64::new(-shift * if region[c] == 0 { 1.0 } else { eps }, 0.0)));
    let mut pos = vec![(f64::MAX, f64::MAX); dofs.n_free];
    for (d, m) in dofs.map.iter().enumerate() {
        if let Some((r, _)) = m {
            if all[d].0 + all[d].1 < pos[*r].0 + pos[*r].1 {
                pos[*r] = all[d];
            }
        }
    }
    (Csr::from_triplets(dofs.n_free, t), pos)
}

/// Mekis's sharp bend, with rods carried on through the layer.
fn bend(cell: f64, pml: f64, freq: f64) -> (Csr, Vec<(f64, f64)>) {
    let (i0, j0, size) = (5i32, 3i32, 9i32);
    let reach = pml.ceil() as i32 + 1;
    let mut rods = vec![];
    for i in -reach..size + reach {
        for j in -reach..size + reach {
            if !((j == j0 && i <= i0) || (i == i0 && j >= j0)) {
                rods.push(circle(i as f64 + 0.5, j as f64 + 0.5, 0.18));
            }
        }
    }
    open(&rods, 3.4 * 3.4, size as f64, cell, pml, C64::new(2.0 * PI * freq, 0.0))
}

fn defect(n: usize) -> Vec<Outline> {
    let m = n / 2;
    (0..n * n).filter(|&i| i != m * n + m).map(|i| circle((i % n) as f64 + 0.5, (i / n) as f64 + 0.5, 0.2)).collect()
}

/// The leaky mode of a defect crystal, from the core's own search, to sit right on its pole.
fn pole(n: usize, cell: f64, pml: f64) -> C64 {
    use glassworks_core::problem::*;
    let c = n as f64 / 2.0;
    let shapes = defect(n).into_iter().map(|outline| Shape { outline, material: Material { eps: 12.0, mu: 1.0 } }).collect();
    let r = Request2d {
        width: n as f64, height: n as f64, cell, background: Material { eps: 1.0, mu: 1.0 }, shapes,
        polarization: Polarization::Te, walls: Boundary::Absorbing,
        analysis: Analysis2d::Resonance { target_k0: 2.0 * PI * 0.38, x: c, y: c, pml }, monitors: vec![],
    };
    let o = glassworks_core::solve2d::solve(&r).unwrap();
    let f = &o.meta.fields[0];
    C64::new(f.k0, -f.k0 / (2.0 * f.q.unwrap()))
}

/// The closed defect mode's `k0^2`, to shift right onto a lossless resonance.
fn closed_mode(cell: f64) -> f64 {
    use glassworks_core::problem::*;
    let shapes = defect(5).into_iter().map(|outline| Shape { outline, material: Material { eps: 12.0, mu: 1.0 } }).collect();
    let r = Request2d {
        width: 5.0, height: 5.0, cell, background: Material { eps: 1.0, mu: 1.0 }, shapes,
        polarization: Polarization::Te, walls: Boundary::Periodic,
        analysis: Analysis2d::Cavity { n_modes: 1, bloch: [0.0, 0.0], target_k0: Some(2.0 * PI * 0.35) }, monitors: vec![],
    };
    glassworks_core::solve2d::solve(&r).unwrap().meta.fields[0].k0.powi(2)
}

pub const CASES: usize = 9;

pub fn case(i: usize) -> Case {
    let rods = |n: usize| defect(n);
    let (name, (a, pos), definite) = match i {
        0 => ("bend a/6, 5a PML", bend(1.0 / 6.0, 5.0, 0.353), false),
        1 => ("bend a/8, 6a PML", bend(1.0 / 8.0, 6.0, 0.353), false),
        2 => ("7x7 cavity, real k at peak", {
            let k = pole(7, 0.1, 0.66);
            open(&rods(7), 12.0, 7.0, 0.1, 0.66, C64::new(k.re, 0.0))
        }, false),
        3 => ("13x13 cavity, on the pole", {
            let k = pole(13, 0.1, 0.66);
            open(&rods(13), 12.0, 13.0, 0.1, 0.66, k)
        }, false),
        4 => ("5x5 supercell, inside the gap", periodic(&rods(5), 12.0, 5.0, 1.0 / 16.0, [0.0, 0.0], (2.0 * PI * 0.35f64).powi(2)), false),
        5 => ("5x5 supercell, on the mode", periodic(&rods(5), 12.0, 5.0, 1.0 / 16.0, [0.0, 0.0], closed_mode(1.0 / 16.0)), false),
        6 => ("empty box, near k^2 = 2 pi^2", periodic(&[], 1.0, 1.0, 0.01, [0.7, 0.3], 2.0 * PI * PI * 1.0001), false),
        7 => ("rod lattice at M, 10,000", periodic(&[circle(0.5, 0.5, 0.2)], 12.0, 1.0, 0.02, [PI, PI], -PI * PI), true),
        _ => ("rod lattice at M, 40,000", periodic(&[circle(0.5, 0.5, 0.2)], 12.0, 1.0, 0.01, [PI, PI], -PI * PI), true),
    };
    Case { name, a, pos, definite }
}

/// `max|A x - b| / (max|A| max|x| + max|b|)`, the usual backward error.
pub fn residual(a: &Csr, x: &[C64], b: &[C64]) -> f64 {
    let mut r = vec![ZERO; a.n];
    a.mul(x, &mut r);
    let inf = |v: &[C64]| v.iter().fold(0.0f64, |m, c| m.max(c.norm()));
    let row = (0..a.n).map(|i| a.val[a.ptr[i]..a.ptr[i + 1]].iter().map(|v| v.norm()).sum::<f64>()).fold(0.0, f64::max);
    r.iter_mut().zip(b).for_each(|(v, w)| *v -= *w);
    inf(&r) / (row * inf(x) + inf(b))
}

/// A fixed pseudo-random right-hand side.
pub fn load(n: usize) -> Vec<C64> {
    let mut seed = 0x9E3779B97F4A7C15u64;
    let mut next = || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        (seed >> 11) as f64 / (1u64 << 53) as f64 - 0.5
    };
    (0..n).map(|_| C64::new(next(), next())).collect()
}

/// The WebAssembly entry points. `run.mjs` holds the clock, since these import nothing.
pub mod wasm {
    use super::*;
    use crate::multifrontal::{Analysis, Multifrontal};
    use glassworks_core::banded::Factors;
    use crate::faer_reference::SparseLu;
    static mut CASE: Option<Case> = None;
    static mut ANALYSIS: Option<Analysis> = None;
    #[allow(static_mut_refs)]
    fn current() -> &'static Case { unsafe { CASE.as_ref().unwrap() } }
    #[no_mangle] #[allow(static_mut_refs)]
    pub extern "C" fn build(i: usize) -> f64 { unsafe { CASE = Some(case(i)); ANALYSIS = None; } current().a.n as f64 }
    #[no_mangle] pub extern "C" fn cases() -> usize { CASES }
    #[no_mangle] #[allow(static_mut_refs)]
    pub extern "C" fn analyse() { unsafe { ANALYSIS = Some(Analysis::new(&current().a, &current().pos)); } }
    #[no_mangle] #[allow(static_mut_refs)]
    pub extern "C" fn run_mf() -> f64 {
        let (a, an) = (&current().a, unsafe { ANALYSIS.as_ref().unwrap() });
        let f = Multifrontal::factor(an, a);
        let b = load(a.n);
        let mut x = b.clone();
        f.solve(&mut x);
        residual(a, &x, &b)
    }
    #[no_mangle] pub extern "C" fn run_faer() -> f64 {
        let a = &current().a;
        let f = SparseLu::factor(a).unwrap();
        let b = load(a.n);
        let mut x = b.clone();
        f.solve(&mut x);
        residual(a, &x, &b)
    }
}
