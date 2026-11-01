//! Checks the handwritten linear algebra against faer on the same matrices.

use faer::prelude::*;
use faer::sparse::{SparseColMat, Triplet};
use faer::Mat;
use eigenmenno_core::assemble::{assemble, Csr, DofMap};
use eigenmenno_core::banded::{rcm, BandedLu, Factors};
use eigenmenno_core::eigen::eigs_near;
use eigenmenno_core::mesh::{Mesh, Outline};
use eigenmenno_core::space::H1Triangle;
use eigenmenno_core::sparse::{Ordering, SparseLu};
use eigenmenno_core::C64;

fn ring(n: usize, shift: C64) -> Csr {
    let mut t = Vec::new();
    for i in 0..n {
        t.push((i, i, C64::new(2.0, 0.0) + shift));
        t.push((i, (i + 1) % n, C64::new(-1.0, 0.3)));
        t.push(((i + 1) % n, i, C64::new(-1.0, -0.3)));
    }
    Csr::from_triplets(n, t)
}

#[test]
fn rcm_folds_a_ring_to_a_narrow_band() {
    let a = ring(1000, C64::new(0.0, 0.0));
    let perm = rcm(&a);
    let mut inv = vec![0; a.n];
    for (new, &old) in perm.iter().enumerate() {
        inv[old] = new;
    }
    let mut b = 0;
    for i in 0..a.n {
        for &j in &a.col[a.ptr[i]..a.ptr[i + 1]] {
            b = b.max(inv[i].abs_diff(inv[j]));
        }
    }
    assert!(b <= 2, "bandwidth {b}");
}

#[test]
fn banded_lu_matches_faer_on_an_indefinite_ring() {
    // An indefinite diagonal forces real pivoting, and the wrap-around entries
    // test the RCM reordering.
    let n = 500;
    let a = ring(n, C64::new(-1.7, 0.01));
    let mut t = Vec::new();
    for i in 0..n {
        for k in a.ptr[i]..a.ptr[i + 1] {
            let v = a.val[k];
            t.push(Triplet::new(i, a.col[k], faer::c64::new(v.re, v.im)));
        }
    }
    let m = SparseColMat::<usize, faer::c64>::try_new_from_triplets(n, n, &t).unwrap();
    let rhs: Vec<C64> = (0..n).map(|i| C64::new((i as f64).sin(), (i as f64 * 0.7).cos())).collect();
    let b = Mat::<faer::c64>::from_fn(n, 1, |i, _| faer::c64::new(rhs[i].re, rhs[i].im));
    let reference = m.sp_lu().unwrap().solve(&b);

    let mut x = rhs.clone();
    BandedLu::factor(&a).unwrap().solve(&mut x);
    for i in 0..n {
        let r = reference[(i, 0)];
        assert!((x[i] - C64::new(r.re, r.im)).norm() < 1e-9, "dof {i}");
    }
}

#[test]
fn lanczos_matches_a_dense_eigensolver() {
    // Hermitian A with a clustered spectrum and a non-identity B.
    let n = 120;
    let a = ring(n, C64::new(0.0, 0.0));
    let b_t: Vec<_> = (0..n).map(|i| (i, i, C64::new(1.0 + 0.5 * (i as f64 * 0.1).sin(), 0.0))).collect();
    let b = Csr::from_triplets(n, b_t);

    // Reduce to a standard problem: B is diagonal, so C = B^-1/2 A B^-1/2.
    let d: Vec<f64> = (0..n).map(|i| b.val[i].re.sqrt().recip()).collect();
    let mut dense = Mat::<faer::c64>::zeros(n, n);
    for i in 0..n {
        for k in a.ptr[i]..a.ptr[i + 1] {
            let (j, v) = (a.col[k], a.val[k]);
            dense[(i, j)] = faer::c64::new(v.re * d[i] * d[j], v.im * d[i] * d[j]);
        }
    }
    let evd = dense.self_adjoint_eigen(faer::Side::Lower).unwrap();
    let mut reference: Vec<f64> = (0..n).map(|i| evd.S()[i].re).collect();
    reference.sort_by(f64::total_cmp);

    let sigma = 1.3;
    let k = 8;
    let found = eigs_near(&a, &b, sigma, k, 1e-11).unwrap();
    let mut nearest = reference.clone();
    nearest.sort_by(|p, q| (p - sigma).abs().total_cmp(&(q - sigma).abs()));
    nearest.truncate(k);
    nearest.sort_by(f64::total_cmp);
    for (f, r) in found.iter().zip(&nearest) {
        assert!((f.0 - r).abs() < 1e-9, "{} vs {r}", f.0);
    }
}

#[test]
fn lanczos_returns_an_entire_small_spectrum() {
    let n = 12;
    let a = ring(n, C64::new(0.0, 0.0));
    let b = Csr::from_triplets(n, (0..n).map(|i| (i, i, C64::new(1.0, 0.0))).collect());
    assert_eq!(eigs_near(&a, &b, -1.0, n, 1e-11).unwrap().len(), n);
}

/// A rod in a box of triangles at an indefinite, lossy shift, with positions for the analysis.
fn helmholtz(shift: C64) -> (Csr, Vec<(f64, f64)>) {
    let mesh = Mesh::fitted_grid(30, 30, 1.0, 1.0, &[Outline::Circle { x: 0.6, y: 0.45, r: 0.2 }]);
    let space = H1Triangle::new(&mesh);
    let all = space.dof_positions();
    let zero: Vec<usize> = (0..all.len()).filter(|&d| [all[d].0, all[d].1].iter().any(|v| *v < 1e-9 || *v > 1.0 - 1e-9)).collect();
    let dofs = DofMap::new(all.len(), &zero, &[]);
    let region = &mesh.region;
    let t = assemble(&space, &dofs, |c| (C64::new(1.0, 0.0), -shift * if region[c] == 0 { 1.0 } else { 12.0 }));
    let free: Vec<(f64, f64)> = (0..all.len()).filter(|&d| dofs.map[d].is_some()).map(|d| all[d]).collect();
    (Csr::from_triplets(dofs.n_free, t), free)
}

#[test]
fn sparse_lu_matches_faer_on_a_2d_matrix() {
    for shift in [C64::new(150.0, 0.0), C64::new(400.0, 3.0)] {
        let (a, pos) = helmholtz(shift);
        let rhs: Vec<C64> = (0..a.n).map(|i| C64::new(((i * 7) % 11) as f64 - 5.0, ((i * 3) % 7) as f64)).collect();
        let mut x = rhs.clone();
        SparseLu::factor(&a, &Ordering::new(&a, &pos)).solve(&mut x);
        let mut triplets = Vec::with_capacity(a.val.len());
        for i in 0..a.n {
            for k in a.ptr[i]..a.ptr[i + 1] {
                triplets.push(Triplet::new(i, a.col[k], faer::c64::new(a.val[k].re, a.val[k].im)));
            }
        }
        let m = SparseColMat::<usize, faer::c64>::try_new_from_triplets(a.n, a.n, &triplets).unwrap();
        let b = Mat::<faer::c64>::from_fn(a.n, 1, |i, _| faer::c64::new(rhs[i].re, rhs[i].im));
        let reference = m.sp_lu().unwrap().solve(&b);
        for i in 0..a.n {
            assert!((x[i] - C64::new(reference[(i, 0)].re, reference[(i, 0)].im)).norm() < 1e-8 * (1.0 + x[i].norm()), "unknown {i}");
        }
    }
}

#[test]
fn sparse_lu_recovers_when_a_part_sits_on_its_own_resonance() {
    // A chain of 129 cuts into a separator and parts, the left one holding unknowns 0 to 62.
    // Shifted to that part's exact eigenvalue, its front is singular, though the chain isn't.
    let n = 129;
    let lambda = 2.0 - 2.0 * (std::f64::consts::PI / 64.0).cos();
    let mut t = vec![];
    for i in 0..n {
        t.push((i, i, C64::new(2.0 - lambda, 0.0)));
        if i + 1 < n {
            t.push((i, i + 1, C64::new(-1.0, 0.0)));
            t.push((i + 1, i, C64::new(-1.0, 0.0)));
        }
    }
    let a = Csr::from_triplets(n, t);
    let pos: Vec<(f64, f64)> = (0..n).map(|i| (i as f64, 0.0)).collect();
    let rhs: Vec<C64> = (0..n).map(|i| C64::new(1.0 + (i % 5) as f64, 0.0)).collect();
    let mut x = rhs.clone();
    SparseLu::factor(&a, &Ordering::new(&a, &pos)).solve(&mut x);
    let mut residual = vec![C64::new(0.0, 0.0); n];
    a.mul(&x, &mut residual);
    let worst = residual.iter().zip(&rhs).fold(0.0f64, |m, (r, b)| m.max((r - b).norm()));
    let size = x.iter().fold(0.0f64, |m, v| m.max(v.norm()));
    eprintln!("chain residual {worst:.1e} on a solution of size {size:.1e}");
    assert!(worst < 1e-12 * (4.0 * size + 5.0), "residual {worst}");
    let mut direct = rhs.clone();
    BandedLu::factor(&a).unwrap().solve(&mut direct);
    for i in 0..n {
        assert!((x[i] - direct[i]).norm() < 1e-9 * size, "unknown {i}");
    }
}
