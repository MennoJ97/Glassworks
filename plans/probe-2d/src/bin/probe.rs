use eigenmenno_core::assemble::{assemble, Csr, DofMap};
use eigenmenno_core::banded::rcm;
use eigenmenno_core::eigen::eigs_near;
use eigenmenno_core::space::Space;
use eigenmenno_core::{C64, ZERO};
use proto2d::*;
use std::f64::consts::PI;
use std::time::Instant;

fn analytic() {
    println!("== PEC square cavity, TM, against m^2 + n^2 (P2 on the structured grid)");
    for nx in [8, 16, 32] {
        let m = grid(nx, nx, 1.0, 1.0, 0.0, 0.0);
        let s = H1Tri::new(&m);
        let pos = s.dof_position();
        let zero: Vec<usize> = pos.iter().enumerate().filter(|(_, p)| p.0 < 1e-9 || p.1 < 1e-9 || p.0 > 1.0 - 1e-9 || p.1 > 1.0 - 1e-9).map(|(i, _)| i).collect();
        let d = DofMap::new(s.n_dofs(), &zero, None);
        let a = Csr::from_triplets(d.n_free, assemble(&s, &d, |_| (C64::new(1.0, 0.0), ZERO)));
        let b = Csr::from_triplets(d.n_free, assemble(&s, &d, |_| (ZERO, C64::new(1.0, 0.0))));
        let t = Instant::now();
        let e = eigs_near(&a, &b, -1.0, 6, 1e-11).unwrap();
        let want = [2.0, 5.0, 5.0, 8.0, 10.0, 10.0];
        let err = e.iter().zip(want).map(|((l, _), w)| (l / (PI * PI) - w).abs() / w).fold(0.0, f64::max);
        println!("  nx {nx:3}: dofs {:6}, worst relative k^2 error {err:.2e}, {:?}", d.n_free, t.elapsed());
    }
}

#[derive(Clone, Copy, PartialEq, Debug)]
enum Fit { Staircase, Snapped }

fn gap(nx: usize, fit: Fit, tm: bool) -> (f64, f64, f64, usize, f64) {
    // Square lattice of eps = 12 rods, r = 0.2 a: the MPB tutorial's structure.
    let (r, eps) = (0.2, 12.0);
    let base = grid(nx, nx, 1.0, 1.0, 0.0, 0.0);
    let mut mesh = grid(nx, nx, 1.0, 1.0, 0.0, 0.0);
    let dist = |x: f64, y: f64| (x - 0.5).hypot(y - 0.5) - r;
    if fit == Fit::Snapped {
        snap(&mut mesh, &dist, &|x, y| { let d = (x - 0.5).hypot(y - 0.5); (0.5 + (x - 0.5) * r / d, 0.5 + (y - 0.5) * r / d) });
    }
    classify(&mut mesh, &|x, y| usize::from(dist(x, y) < 0.0));
    let angle = min_angle(&mesh);
    let space = H1Tri::new(&mesh);
    let coef = |cell: usize| -> (C64, C64) {
        let e = if mesh.region[cell] == 1 { eps } else { 1.0 };
        if tm { (C64::new(1.0, 0.0), C64::new(e, 0.0)) } else { (C64::new(1.0 / e, 0.0), C64::new(1.0, 0.0)) }
    };
    let bands = |ph: (f64, f64), n: usize| -> Vec<f64> {
        let d = bloch_map(&space, &base, nx, nx, 1.0, ph);
        let a = Csr::from_triplets(d.n_free, assemble(&space, &d, |c| (coef(c).0, ZERO)));
        let b = Csr::from_triplets(d.n_free, assemble(&space, &d, |c| (ZERO, coef(c).1)));
        eigs_near(&a, &b, -1.0, n, 1e-10).unwrap().into_iter().map(|(l, _)| l.max(0.0).sqrt() / (2.0 * PI)).collect()
    };
    let t = Instant::now();
    let m = bands((PI, PI), 3);
    let per_k = t.elapsed().as_secs_f64() * 1000.0;
    let x = bands((PI, 0.0), 3);
    let n = bloch_map(&space, &base, nx, nx, 1.0, (PI, 0.0)).n_free;
    (m[0], x[1], per_k, n, angle)
}

fn bands() {
    println!("== Band edges, square lattice, eps 12, r 0.2a. MPB tutorial TM gap: 0.28262 to 0.41933");
    for tm in [true, false] {
        println!("  polarization {}", if tm { "TM (Ez)" } else { "TE (Hz)" });
        let (fine1, fine2, _, _, _) = gap(96, Fit::Snapped, tm);
        println!("    finest snapped, nx 96: band 1 at M {fine1:.5}, band 2 at X {fine2:.5}");
        for fit in [Fit::Staircase, Fit::Snapped] {
            for nx in [12, 24, 48] {
                let (a, b, ms, n, ang) = gap(nx, fit, tm);
                println!("    {fit:?} nx {nx:3}: dofs {n:6}, M1 {a:.5} (err {:.1e}), X2 {b:.5} (err {:.1e}), {ms:6.1} ms for 3 modes, min angle {ang:.1}", (a - fine1).abs(), (b - fine2).abs());
            }
        }
    }
}

fn bessel(x: f64) -> C64 {
    let n = 20000;
    let (mut j, mut y) = (0.0, 0.0);
    for i in 0..n {
        let t = (i as f64 + 0.5) / n as f64;
        j += (x * (PI * t).sin()).cos() / n as f64;
        let th = t * PI / 2.0;
        y += (x * th.cos()).cos() * (0.5772156649 + (2.0 * x * th.sin().powi(2)).ln()) * (PI / 2.0) / n as f64;
    }
    C64::new(j, y * 4.0 / (PI * PI))
}

fn pml() {
    println!("== Line source with PML against (i/4) H0(kr), wavelength 1, ring at r = 0.8 to 1.2");
    let k = 2.0 * PI;
    for (half, t, per_wave) in [(1.5, 0.25, 12), (1.5, 0.5, 12), (1.5, 1.0, 12), (1.5, 0.5, 6), (1.5, 0.5, 24)] {
        let width = 2.0 * (half + t);
        let nx = (width * per_wave as f64).round() as usize / 2 * 2;
        let mesh = grid(nx, nx, width, width, -width / 2.0, -width / 2.0);
        let s = H1Tri::new(&mesh);
        let pos = s.dof_position();
        let zero: Vec<usize> = pos.iter().enumerate().filter(|(_, p)| p.0.abs() > width / 2.0 - 1e-9 || p.1.abs() > width / 2.0 - 1e-9).map(|(i, _)| i).collect();
        let d = DofMap::new(s.n_dofs(), &zero, None);
        let a = helmholtz_pml(&s, &d, half, t, k);
        let src = pos.iter().position(|p| p.0.abs() < 1e-9 && p.1.abs() < 1e-9).unwrap();
        let mut b = vec![ZERO; d.n_free];
        b[d.map[src].unwrap().0] = C64::new(1.0, 0.0);
        let tt = Instant::now();
        eigenmenno_core::banded::BandedLu::factor(&a).unwrap().solve(&mut b);
        let ms = tt.elapsed().as_secs_f64() * 1000.0;
        let u = d.expand(&b);
        let (mut worst, mut count) = (0f64, 0);
        for (i, p) in pos.iter().enumerate().step_by(7) {
            let r = p.0.hypot(p.1);
            if (0.8..=1.2).contains(&r) {
                let exact = C64::new(0.0, 0.25) * bessel(k * r);
                worst = worst.max((u[i] - exact).norm() / exact.norm());
                count += 1;
            }
        }
        println!("  PML {t:4} wavelengths, {per_wave:2} P2 cells per wavelength: dofs {:6}, worst error {:.2}% over {count} points, {ms:.0} ms", d.n_free, 100.0 * worst);
    }
}

fn scaling() {
    println!("== Native solve of a 2D PML Helmholtz matrix");
    for nx in [25, 50, 100, 200, 400, 800] {
        let a = scaling_matrix(nx);
        let perm = rcm(&a);
        let mut inv = vec![0; a.n];
        perm.iter().enumerate().for_each(|(new, &old)| inv[old] = new);
        let mut bw = 0;
        for i in 0..a.n { for &j in &a.col[a.ptr[i]..a.ptr[i + 1]] { bw = bw.max(inv[i].abs_diff(inv[j])); } }
        let gb = (3 * bw + 1) as f64 * a.n as f64 * 16.0 / 1e9;
        let banded = if gb < 3.0 { let t = Instant::now(); let x = solve_banded(&a); format!("{:9.1} ms ({:.3})", t.elapsed().as_secs_f64() * 1000.0, x[a.n / 2].norm()) } else { "skipped".to_string() };
        #[cfg(feature = "sparse")]
        let sparse = { let t = Instant::now(); let x = solve_faer(&a); format!("{:9.1} ms ({:.3})", t.elapsed().as_secs_f64() * 1000.0, x[a.n / 2].norm()) };
        #[cfg(not(feature = "sparse"))]
        let sparse = "";
        println!("  nx {nx:4}: dofs {:7}, bandwidth {bw:5}, band storage {gb:6.2} GB, banded {banded}, faer {sparse}", a.n);
    }
}

fn main() {
    analytic();
    bands();
    pml();
    scaling();
}
