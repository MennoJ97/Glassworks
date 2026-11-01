use eigenmenno_core::banded::Factors;
use eigenmenno_core::sparse::{SparseCholesky, SparseLu};
use eigenmenno_core::C64;
use probe_sparse::multifrontal::{Analysis, Multifrontal};
use probe_sparse::*;
use std::time::Instant;

fn refine(c: &Case, f: &Multifrontal, x: &mut [C64], b: &[C64]) {
    let mut r = vec![eigenmenno_core::ZERO; c.a.n];
    c.a.mul(x, &mut r);
    r.iter_mut().zip(b).for_each(|(v, w)| *v = *w - *v);
    f.solve(&mut r);
    x.iter_mut().zip(&r).for_each(|(v, d)| *v += *d);
}

fn main() {
    println!("case | n | who | analyse ms | factor ms | solve ms | peak MB | residual | refined | replaced | weak | diff to faer");
    let only: Vec<usize> = std::env::args().skip(1).map(|s| s.parse().unwrap()).collect();
    for i in (0..CASES).filter(|i| only.is_empty() || only.contains(i)) {
        let c = case(i);
        let b = load(c.a.n);
        // faer's general LU.
        let base = reset_peak();
        let t = Instant::now();
        let lu = SparseLu::factor(&c.a).unwrap();
        let factor = t.elapsed().as_secs_f64() * 1e3;
        let peak = peak_mb(base);
        let t = Instant::now();
        let mut truth = b.clone();
        lu.solve(&mut truth);
        let solve = t.elapsed().as_secs_f64() * 1e3;
        println!("{} | {} | faer LU | - | {factor:.0} | {solve:.1} | {peak:.0} | {:.1e} | | | |", c.name, c.a.n, residual(&c.a, &truth, &b));
        drop(lu);
        if c.definite {
            let base = reset_peak();
            let t = Instant::now();
            let llt = SparseCholesky::factor(&c.a, &mut None).unwrap();
            let factor = t.elapsed().as_secs_f64() * 1e3;
            let peak = peak_mb(base);
            let t = Instant::now();
            let mut x = b.clone();
            llt.solve(&mut x);
            println!("{} | {} | faer LLT | - | {factor:.0} | {:.1} | {peak:.0} | {:.1e} | | | |", c.name, c.a.n, t.elapsed().as_secs_f64() * 1e3, residual(&c.a, &x, &b));
        }
        // Ours.
        let base = reset_peak();
        let t = Instant::now();
        let an = Analysis::new(&c.a, &c.pos);
        let analyse = t.elapsed().as_secs_f64() * 1e3;
        let t = Instant::now();
        let f = Multifrontal::factor(&an, &c.a);
        let factor = t.elapsed().as_secs_f64() * 1e3;
        let peak = peak_mb(base);
        let t = Instant::now();
        let mut x = b.clone();
        f.solve(&mut x);
        let solve = t.elapsed().as_secs_f64() * 1e3;
        let first = residual(&c.a, &x, &b);
        for _ in 0..2 {
            refine(&c, &f, &mut x, &b);
        }
        let norm = truth.iter().map(|v| v.norm_sqr()).sum::<f64>().sqrt();
        let diff = x.iter().zip(&truth).map(|(p, q)| (p - q).norm_sqr()).sum::<f64>().sqrt() / norm;
        println!("{} | {} | ours | {analyse:.0} | {factor:.0} | {solve:.1} | {peak:.0} | {first:.1e} | {:.1e} | {} | {} | {diff:.1e}", c.name, c.a.n, residual(&c.a, &x, &b), f.replaced, f.weak);
    }
}
