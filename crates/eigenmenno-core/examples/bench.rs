//! Native timings: `cargo run --release -p eigenmenno-core --example bench`.

use eigenmenno_core::assemble::{assemble, Csr, DofMap};
use eigenmenno_core::banded::{BandedLu, Factors};
use eigenmenno_core::mesh::Mesh;
use eigenmenno_core::problem::{Analysis, Boundary, Layer, MeshSpec, Polarization, Request};
use eigenmenno_core::solve1d::solve;
use eigenmenno_core::space::{H1Line, Space};
use eigenmenno_core::C64;
use faer::prelude::*;
use faer::sparse::{SparseColMat, Triplet};
use faer::Mat;
use std::time::Instant;

fn main() {
    println!("linear solve: banded LU against faer sparse LU (complex Helmholtz matrix)");
    for (order, cells) in [(1, 100_000), (1, 1_000_000), (4, 250_000)] {
        let mesh = Mesh::layered_line(&[1.0], &[1.0 / cells as f64]);
        let space = H1Line::new(&mesh, order);
        let dofs = DofMap::new(space.n_dofs(), &[0, space.n_dofs() - 1], &[]);
        let k2 = C64::new(1.0e4, 10.0);
        let a = Csr::from_triplets(dofs.n_free, assemble(&space, &dofs, |_| (C64::new(1.0, 0.0), -k2)));
        let rhs: Vec<C64> = (0..a.n).map(|i| C64::new(1.0, (i % 7) as f64)).collect();

        let t = Instant::now();
        let mut x = rhs.clone();
        BandedLu::factor(&a).unwrap().solve(&mut x);
        let ours = t.elapsed();

        let mut trip = Vec::with_capacity(a.val.len());
        for i in 0..a.n {
            for k in a.ptr[i]..a.ptr[i + 1] {
                let v = a.val[k];
                trip.push(Triplet::new(i, a.col[k], faer::c64::new(v.re, v.im)));
            }
        }
        let m = SparseColMat::<usize, faer::c64>::try_new_from_triplets(a.n, a.n, &trip).unwrap();
        let b = Mat::<faer::c64>::from_fn(a.n, 1, |i, _| faer::c64::new(rhs[i].re, rhs[i].im));
        let t = Instant::now();
        let y = m.sp_lu().unwrap().solve(&b);
        let theirs = t.elapsed();

        let diff = (0..a.n).map(|i| (x[i] - C64::new(y[(i, 0)].re, y[(i, 0)].im)).norm()).fold(0.0, f64::max);
        println!(
            "  p={order} dofs={:>8}  banded {:>8.1?}  faer {:>8.1?}  max difference {diff:.1e}",
            a.n, ours, theirs
        );
    }

    println!("full problems");
    let layers = |n: usize| -> Vec<Layer> {
        (0..n)
            .map(|i| Layer { thickness: 0.2, eps: [if i % 2 == 0 { 1.0 } else { 6.0 }, 0.0], mu: [1.0, 0.0] })
            .collect()
    };
    let cases = [
        (
            "cavity, 10 modes, 20 layers",
            Analysis::Cavity { n_modes: 10, target_k0: None, bloch_phase: 0.0 },
            Boundary::Pec,
            Boundary::Pec,
            20,
        ),
        (
            "bands, 41 k x 6 bands",
            Analysis::Bands { n_bands: 6, n_k: 41 },
            Boundary::Periodic,
            Boundary::Periodic,
            2,
        ),
        (
            "sweep, 500 wavelengths",
            Analysis::Sweep {
                wavelengths: (0..500).map(|i| 0.5 + 0.002 * i as f64).collect(),
                angle_deg: 0.0,
            },
            Boundary::Port,
            Boundary::Absorbing,
            20,
        ),
    ];
    for (name, analysis, left, right, n) in cases {
        let r = Request {
            layers: layers(n),
            polarization: Polarization::Te,
            left,
            right,
            mesh: MeshSpec { order: 4, max_cell: 0.02 },
            analysis,
        };
        let t = Instant::now();
        let o = solve(&r).unwrap();
        println!("  {name:<30} dofs={:>6}  {:>8.1?}", o.meta.n_dofs, t.elapsed());
    }
}
