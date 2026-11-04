//! Shift-invert Lanczos for `A x = lambda B x`, with A Hermitian and B positive
//! definite. Inverting `A - sigma B` maps the eigenvalues nearest `sigma` to the
//! largest ones, which Lanczos finds first. With vector fields, this would also
//! skip the large null space of the curl-curl operator.

use crate::assemble::Csr;
use crate::banded::{BandedLu, Factors};
use crate::{C64, ZERO};

/// Largest Krylov basis, in vectors: the only limit on clustered spectra.
const MAX_BASIS: usize = 3000;

/// An eigenpair of the inverted operator, with `theta = 1 / (lambda - sigma)`.
pub struct Mode {
    pub theta: f64,
    /// B-orthonormal eigenvector, and B times it.
    pub x: Vec<C64>,
    bx: Vec<C64>,
}

/// `a^H b`, conjugating the first vector.
fn dot(a: &[C64], b: &[C64]) -> C64 {
    a.iter().zip(b).map(|(x, y)| x.conj() * y).sum()
}

/// `y += alpha x`, named after the BLAS routine.
fn axpy(alpha: C64, x: &[C64], y: &mut [C64]) {
    y.iter_mut().zip(x).for_each(|(yi, xi)| *yi += *xi * alpha);
}

/// Implicit QL for a symmetric tridiagonal matrix with diagonal `diag` and
/// off-diagonal `off`. Returns eigenvalues and vectors, with `vectors[k * n + i]` the
/// k-th component of vector i. Each sweep chases a bulge down with plane rotations
/// (`cos`, `sin`); `g`, `r`, `p`, `f` and `b` are the textbook's working values.
fn tridiagonal_eigen(diag: &[f64], off: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let n = diag.len();
    let (mut d, mut e) = (diag.to_vec(), off.to_vec());
    e.resize(n, 0.0);
    let mut vectors = vec![0.0; n * n];
    (0..n).for_each(|i| vectors[i * n + i] = 1.0);
    for l in 0..n {
        for _ in 0..200 {
            // The block d[l..=m] is unreduced: its last coupling e[m] is negligible.
            let m = (l..n - 1)
                .find(|&m| e[m].abs() <= f64::EPSILON * (d[m].abs() + d[m + 1].abs()))
                .unwrap_or(n - 1);
            if m == l {
                break;
            }
            let mut g = (d[l + 1] - d[l]) / (2.0 * e[l]);
            let mut r = g.hypot(1.0);
            g = d[m] - d[l] + e[l] / (g + r.copysign(g));
            let (mut sin, mut cos, mut p) = (1.0, 1.0, 0.0);
            for i in (l..m).rev() {
                let (f, b) = (sin * e[i], cos * e[i]);
                r = f.hypot(g).max(f64::MIN_POSITIVE);
                e[i + 1] = r;
                sin = f / r;
                cos = g / r;
                g = d[i + 1] - p;
                r = (d[i] - g) * sin + 2.0 * cos * b;
                p = sin * r;
                d[i + 1] = g + p;
                g = cos * r - b;
                for k in 0..n {
                    let f = vectors[k * n + i + 1];
                    vectors[k * n + i + 1] = sin * vectors[k * n + i] + cos * f;
                    vectors[k * n + i] = cos * vectors[k * n + i] - sin * f;
                }
            }
            d[l] -= p;
            e[l] = g;
            e[m] = 0.0;
        }
    }
    (d, vectors)
}

/// One Lanczos run, kept orthogonal to `found`. Returns the converged modes, aiming for
/// `want` of them. `lu` factors `A - sigma B`, and `seed` carries the random start between runs.
fn lanczos(lu: &dyn Factors, b: &Csr, found: &[Mode], want: usize, tol: f64, seed: &mut u64) -> Vec<Mode> {
    let n = b.n;
    let max_dim = (n - found.len()).min(MAX_BASIS);
    // `v` is the basis and `bv` is B times it. `alpha` and `beta` are the diagonal and
    // off-diagonal of the tridiagonal matrix `T`, the operator projected onto the basis.
    let (mut v, mut bv): (Vec<Vec<C64>>, Vec<Vec<C64>>) = (vec![], vec![]);
    let (mut alpha, mut beta): (Vec<f64>, Vec<f64>) = (vec![], vec![]);
    let mut next_check = want.max(8);
    // The next basis vector, before it's normalized. A deterministic random start means
    // repeated runs give identical results.
    let mut next: Vec<C64> = (0..n)
        .map(|_| {
            *seed ^= *seed << 13;
            *seed ^= *seed >> 7;
            *seed ^= *seed << 17;
            C64::new((*seed >> 11) as f64 / (1u64 << 53) as f64 - 0.5, 0.0)
        })
        .collect();
    found.iter().for_each(|mode| axpy(-dot(&mode.bx, &next), &mode.x, &mut next));
    let mut b_next = vec![ZERO; n];
    loop {
        b.mul(&next, &mut b_next);
        // The B-norm of what's left: beta_j, coupling the last vector to this one.
        let norm = dot(&next, &b_next).re.max(0.0).sqrt();
        let scale = alpha.iter().fold(0.0, |largest: f64, a| largest.max(a.abs()));
        // Out of room, or the basis already spans an invariant subspace.
        let stop = v.len() == max_dim || (!v.is_empty() && norm <= 1e-10 * scale);
        if !stop {
            if !v.is_empty() {
                beta.push(norm);
            }
            v.push(next.iter().map(|x| x / norm).collect());
            bv.push(b_next.iter().map(|x| x / norm).collect());
            let j = v.len() - 1;
            next.copy_from_slice(&bv[j]);
            lu.solve(&mut next);
            alpha.push(dot(&bv[j], &next).re);
            // Full reorthogonalization, twice: why the basis is stored at all.
            for _ in 0..2 {
                (0..=j).for_each(|i| axpy(-dot(&bv[i], &next), &v[i], &mut next));
                found.iter().for_each(|mode| axpy(-dot(&mode.bx, &next), &mode.x, &mut next));
            }
            if v.len() < next_check {
                continue;
            }
        }
        // Checking costs O(m^3), so do it at geometrically spaced sizes.
        let m = v.len();
        next_check = (m + 8).max(m * 5 / 4);
        let (theta, tri_vectors) = tridiagonal_eigen(&alpha, &beta);
        let mut order: Vec<usize> = (0..m).collect();
        order.sort_by(|&p, &q| theta[q].abs().total_cmp(&theta[p].abs()));
        // A Ritz pair's residual is the coupling to the next basis vector times
        // the last component of its eigenvector.
        b.mul(&next, &mut b_next);
        let coupling = dot(&next, &b_next).re.max(0.0).sqrt();
        let converged = |i: usize| coupling * tri_vectors[(m - 1) * m + i].abs() <= tol * theta[i].abs();
        if stop || order.iter().take(want).all(|&i| converged(i)) {
            return order
                .into_iter()
                .filter(|&i| converged(i))
                .map(|i| {
                    let mut x = vec![ZERO; n];
                    (0..m).for_each(|k| axpy(C64::new(tri_vectors[k * m + i], 0.0), &v[k], &mut x));
                    let mut bx = vec![ZERO; n];
                    b.mul(&x, &mut bx);
                    Mode { theta: theta[i], x, bx }
                })
                .collect();
        }
    }
}

/// The `k` eigenpairs `(lambda, x)` closest to `sigma`, by increasing lambda.
///
/// A single Krylov space holds only one vector per eigenspace, so degenerate
/// modes (Bloch bands at the zone edge) would be missed. Each run therefore
/// stays orthogonal to the modes already found, and we stop once a run finds
/// nothing closer than the current k-th mode.
pub fn eigs_near(a: &Csr, b: &Csr, sigma: f64, k: usize, tol: f64) -> Result<Vec<(f64, Vec<C64>)>, String> {
    eigs_factored(&BandedLu::factor(&a.sub_scaled(sigma, b))?, b, sigma, k, tol)
}

/// `eigs_near` for an `lu` that already factors `A - sigma B`.
pub fn eigs_factored(
    lu: &dyn Factors,
    b: &Csr,
    sigma: f64,
    k: usize,
    tol: f64,
) -> Result<Vec<(f64, Vec<C64>)>, String> {
    let k = k.min(b.n);
    let mut seed = 0x9E3779B97F4A7C15;
    let mut found: Vec<Mode> = vec![];
    while found.len() < b.n {
        let new = lanczos(lu, b, &found, k.saturating_sub(found.len()).max(1), tol, &mut seed);
        let kth = found.get(k.wrapping_sub(1)).map_or(0.0, |mode| mode.theta.abs());
        let finished = new.iter().all(|mode| mode.theta.abs() < kth * (1.0 - 1e-9));
        found.extend(new);
        found.sort_by(|p, q| q.theta.abs().total_cmp(&p.theta.abs()));
        if finished {
            break;
        }
    }
    if found.len() < k {
        return Err(format!("the eigensolver converged only {} of {k} modes", found.len()));
    }
    let mut pairs: Vec<_> =
        found.into_iter().take(k).map(|mode| (sigma + 1.0 / mode.theta, mode.x)).collect();
    pairs.sort_by(|p, q| p.0.total_cmp(&q.0));
    Ok(pairs)
}
