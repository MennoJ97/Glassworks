//! The hardest case for pivoting only inside each front: a matrix that's fine as a whole, but
//! whose subdomain, walled in by its separator, sits on a resonance of its own. Elimination
//! there meets a nearly singular block, and no later unknown can be used to pivot it.
use eigenmenno_core::assemble::{assemble, Csr, DofMap};
use eigenmenno_core::banded::Factors;
use eigenmenno_core::eigen::eigs_near;
use eigenmenno_core::mesh::{Mesh, Outline};
use eigenmenno_core::space::H1Triangle;
use probe_sparse::faer_reference::SparseLu;
use eigenmenno_core::{C64, ONE, ZERO};
use probe_sparse::multifrontal::{Analysis, Multifrontal};
use probe_sparse::{load, residual};

fn sub(a: &Csr, ids: &[usize]) -> Csr {
    let mut map = vec![usize::MAX; a.n];
    ids.iter().enumerate().for_each(|(i, &d)| map[d] = i);
    let mut t = vec![];
    for (i, &d) in ids.iter().enumerate() {
        for k in a.ptr[d]..a.ptr[d + 1] {
            if map[a.col[k]] != usize::MAX {
                t.push((i, map[a.col[k]], a.val[k]));
            }
        }
    }
    Csr::from_triplets(ids.len(), t)
}

fn main() {
    // Vacuum on the left, and a rod on the right so that the whole has different modes.
    let n = 50;
    let rod = [Outline::Circle { x: 0.75, y: 0.45, r: 0.12 }];
    let mesh = Mesh::fitted_grid(n, n, 1.0, 1.0, &rod);
    let space = H1Triangle::new(&mesh);
    let all = space.dof_positions();
    let edge = |p: (f64, f64)| p.0 < 1e-9 || p.1 < 1e-9 || p.0 > 1.0 - 1e-9 || p.1 > 1.0 - 1e-9;
    let zero: Vec<usize> = (0..all.len()).filter(|&d| edge(all[d])).collect();
    let dofs = DofMap::new(all.len(), &zero, &[]);
    let mut pos = vec![(0.0, 0.0); dofs.n_free];
    for (d, m) in dofs.map.iter().enumerate() {
        if let Some((r, _)) = m { pos[*r] = all[d]; }
    }
    let region = &mesh.region;
    let s = Csr::from_triplets(dofs.n_free, assemble(&space, &dofs, |_| (ONE, ZERO)));
    let m = Csr::from_triplets(dofs.n_free, assemble(&space, &dofs, |c| (ZERO, C64::new(if region[c] == 0 { 1.0 } else { 12.0 }, 0.0))));
    let an = Analysis::new(&s, &pos);
    println!("n {}", s.n);
    println!("depth | subdomain | delta | faer residual | ours residual | refined once | twice | replaced | weak | diff to faer | nearest global mode");
    for depth in [1, 2, 4] {
        let ids = an.subdomain(depth);
        let (ss, ms) = (sub(&s, &ids), sub(&m, &ids));
        let lambda = eigs_near(&ss, &ms, -1.0, 1, 1e-12).unwrap()[0].0;
        let global = eigs_near(&s, &m, lambda, 1, 1e-10).unwrap()[0].0;
        for delta in [1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-14, 0.0] {
            let a = s.sub_scaled(lambda * (1.0 + delta), &m);
            let b = load(a.n);
            let mut truth = b.clone();
            SparseLu::factor(&a).unwrap().solve(&mut truth);
            let f = Multifrontal::factor(&an, &a);
            let mut x = b.clone();
            f.solve(&mut x);
            let first = residual(&a, &x, &b);
            let mut refined = [0.0; 2];
            for r in &mut refined {
                let mut e = vec![ZERO; a.n];
                a.mul(&x, &mut e);
                e.iter_mut().zip(&b).for_each(|(v, w)| *v = *w - *v);
                f.solve(&mut e);
                x.iter_mut().zip(&e).for_each(|(v, d)| *v += *d);
                *r = residual(&a, &x, &b);
            }
            let norm = truth.iter().map(|v| v.norm_sqr()).sum::<f64>().sqrt();
            let diff = x.iter().zip(&truth).map(|(p, q)| (p - q).norm_sqr()).sum::<f64>().sqrt() / norm;
            println!("{depth} | {} | {delta:.0e} | {:.1e} | {first:.1e} | {:.1e} | {:.1e} | {} | {} | {diff:.1e} | {:.2e}", ids.len(), residual(&a, &truth, &b), refined[0], refined[1], f.replaced, f.weak, (global - lambda).abs() / lambda);
        }
    }
}
