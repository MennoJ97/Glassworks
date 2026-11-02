//! Multifrontal sparse LU for 2D, where the band in `banded` fills in too much. Cutting the
//! domain in two along its longer side leaves two parts joined only through a separator, so
//! each can be eliminated without touching the other. Doing that recursively is nested
//! dissection ([George 1973]). Each separator's unknowns are then eliminated in a dense "front",
//! which also holds whatever its two parts left behind. See `plans/2d-solver.md` for timings.
//!
//! Pivoting stays inside each front's own unknowns, because a pivot from elsewhere would
//! change which fronts need what. A part walled in by its separator can sit on a resonance of
//! its own, though, and then its front has only a poor pivot, even if the whole matrix is
//! fine. So a front reports any weak pivot, and solves are then refined against the matrix.

use crate::assemble::Csr;
use crate::banded::Factors;
use crate::{C64, ZERO};
use std::cell::RefCell;
use std::rc::Rc;

/// Parts this small are eliminated whole.
const LEAF: usize = 16;
/// Columns eliminated together before updating the rest of a front, to stay in cache.
const BLOCK: usize = 32;
const REFINEMENT_STEPS: usize = 2;

struct Node {
    /// This node's pivots are `start..end` in the new numbering.
    start: usize,
    end: usize,
    children: Vec<usize>,
    /// New indices of later unknowns this front couples to, sorted.
    boundary: Vec<usize>,
    /// Matrix entries added here: `(index into a.val, row, column)`, within the front.
    entries: Vec<(u32, u32, u32)>,
}

/// The elimination order and fill, which depend only on where the nonzeros are. It can be
/// reused for matrices with the same pattern, as Bloch phases and wavelengths have.
pub struct Ordering {
    /// `new[old]` numbering.
    new: Vec<usize>,
    nodes: Vec<Node>,
}

/// Splits `ids` in two and recurses, pushing `(pivots, children)` with children first.
fn dissect(ids: Vec<usize>, a: &Csr, pos: &[(f64, f64)], side: &mut [u8], nodes: &mut Vec<(Vec<usize>, Vec<usize>)>) -> usize {
    if ids.len() <= LEAF {
        nodes.push((ids, vec![]));
        return nodes.len() - 1;
    }
    let extent = |f: &dyn Fn(usize) -> f64| {
        let (low, high) = ids.iter().fold((f64::MAX, f64::MIN), |b, &d| (b.0.min(f(d)), b.1.max(f(d))));
        high - low
    };
    let along_x = extent(&|d| pos[d].0) >= extent(&|d| pos[d].1);
    let key = |d: usize| if along_x { pos[d].0 } else { pos[d].1 };
    let mut left = ids;
    let mid = left.len() / 2;
    left.select_nth_unstable_by(mid, |&p, &q| key(p).total_cmp(&key(q)));
    let mut right = left.split_off(mid);
    right.iter().for_each(|&d| side[d] = 2);
    // The separator comes from the matrix, not the grid: whatever on the left touches the
    // right, whether across the cut, over a flipped edge, or through a periodic tie.
    let touches = |d: usize, s: &[u8], v: u8| a.col[a.ptr[d]..a.ptr[d + 1]].iter().any(|&c| s[c] == v);
    let (mut separator, left): (Vec<_>, Vec<_>) = left.into_iter().partition(|&d| touches(d, side, 2));
    left.iter().for_each(|&d| side[d] = 1);
    // Elements couple nodes two columns apart, which can make that two columns thick. A node
    // touching nothing further left can join the right, unless there's no left at all: that
    // would only rebuild the same set, and never finish.
    if !left.is_empty() {
        separator.retain(|&d| {
            let needed = touches(d, side, 1);
            if !needed {
                right.push(d);
            }
            needed
        });
    }
    left.iter().chain(&right).for_each(|&d| side[d] = 0);
    let children = [left, right].into_iter().filter(|p| !p.is_empty()).map(|p| dissect(p, a, pos, side, nodes)).collect();
    nodes.push((separator, children));
    nodes.len() - 1
}

impl Ordering {
    /// `pos` gives each unknown's position. `a` must be structurally symmetric.
    pub fn new(a: &Csr, pos: &[(f64, f64)]) -> Self {
        let n = a.n;
        let mut raw = vec![];
        dissect((0..n).collect(), a, pos, &mut vec![0; n], &mut raw);
        // Children come before parents, so numbering in this order eliminates them first.
        let (mut new, mut owner, mut next) = (vec![0; n], vec![0; n], 0);
        let mut nodes: Vec<Node> = vec![];
        for (i, (pivots, children)) in raw.iter().enumerate() {
            let start = next;
            for &d in pivots {
                new[d] = next;
                owner[next] = i;
                next += 1;
            }
            nodes.push(Node { start, end: next, children: children.clone(), boundary: vec![], entries: vec![] });
        }
        // A front couples to its own entries' later unknowns, and to what its children left.
        let mut seen = vec![usize::MAX; n];
        for i in 0..nodes.len() {
            let end = nodes[i].end;
            let own = raw[i].0.iter().flat_map(|&d| a.col[a.ptr[d]..a.ptr[d + 1]].iter().map(|&c| new[c]));
            let inherited = raw[i].1.iter().flat_map(|&c| nodes[c].boundary.iter().copied());
            let mut boundary: Vec<usize> = own.chain(inherited).filter(|&g| g >= end && std::mem::replace(&mut seen[g], i) != i).collect();
            boundary.sort_unstable();
            nodes[i].boundary = boundary;
        }
        // Each entry belongs to the first front that holds both its row and its column.
        for r in 0..n {
            for k in a.ptr[r]..a.ptr[r + 1] {
                let (i, j) = (new[r], new[a.col[k]]);
                let node = &mut nodes[owner[i.min(j)]];
                let local = |g: usize| if g < node.end { g - node.start } else { node.end - node.start + node.boundary.binary_search(&g).unwrap() };
                let entry = (k as u32, local(i) as u32, local(j) as u32);
                node.entries.push(entry);
            }
        }
        Self { new, nodes }
    }

    /// The last ordering made if this matrix has its nonzeros, which is cheap to check, and
    /// otherwise a new one. A design's iterations and a band diagram's k-points all reuse it.
    pub fn reused(a: &Csr, pos: impl FnOnce() -> Vec<(f64, f64)>) -> Rc<Self> {
        LAST.with_borrow_mut(|last| match last {
            Some((ptr, col, ordering)) if *ptr == a.ptr && *col == a.col => ordering.clone(),
            _ => {
                let ordering = Rc::new(Self::new(a, &pos()));
                *last = Some((a.ptr.clone(), a.col.clone(), ordering.clone()));
                ordering
            }
        })
    }
}

thread_local! {
    static LAST: RefCell<Option<(Vec<usize>, Vec<usize>, Rc<Ordering>)>> = const { RefCell::new(None) };
}

/// `y -= sum(u[q] x[q])`, over four columns at once for speed.
#[inline]
fn subtract(y: &mut [C64], x: [&[C64]; 4], u: [C64; 4]) {
    let n = y.len();
    let [x0, x1, x2, x3] = [&x[0][..n], &x[1][..n], &x[2][..n], &x[3][..n]];
    for i in 0..n {
        y[i] -= u[0] * x0[i] + u[1] * x1[i] + u[2] * x2[i] + u[3] * x3[i];
    }
}

/// Eliminates the first `p` of a front's `m` unknowns in place, swapping rows among those `p`.
/// `f` is column-major. Returns the swaps, and whether any pivot was weak.
fn partial_lu(f: &mut [C64], m: usize, p: usize, scale: f64) -> (Vec<u32>, bool) {
    let mut swaps = Vec::with_capacity(p);
    let mut weak = false;
    let mut k = 0;
    while k < p {
        let e = (k + BLOCK).min(p);
        for j in k..e {
            let col = &f[j * m..(j + 1) * m];
            let (mut r, mut best) = (j, 0.0);
            for i in j..p {
                if col[i].norm_sqr() > best {
                    (r, best) = (i, col[i].norm_sqr());
                }
            }
            let below = col[p..].iter().fold(0.0f64, |b, v| b.max(v.norm_sqr()));
            weak |= best < 1e-4 * below || best < 1e-16 * scale * scale;
            if r != j {
                for c in 0..m {
                    f.swap(j + c * m, r + c * m);
                }
            }
            swaps.push(r as u32);
            // Refinement corrects for this, unless the whole matrix is singular.
            if f[j + j * m].norm() < 1e-14 * scale {
                f[j + j * m] = C64::new(1e-14 * scale, 0.0);
            }
            let inv = f[j + j * m].inv();
            f[j * m + j + 1..(j + 1) * m].iter_mut().for_each(|v| *v *= inv);
            let (left, right) = f.split_at_mut((j + 1) * m);
            let l = &left[j * m + j + 1..];
            for col in right.chunks_exact_mut(m).take(e - j - 1) {
                let u = col[j];
                col[j + 1..].iter_mut().zip(l).for_each(|(y, x)| *y -= u * *x);
            }
        }
        // Update the later columns from this block.
        let (left, right) = f.split_at_mut(e * m);
        for col in right.chunks_exact_mut(m) {
            for q in k..e {
                let u = col[q];
                for i in q + 1..e {
                    col[i] -= u * left[q * m + i];
                }
            }
            let (head, tail) = col.split_at_mut(e);
            let mut q = k;
            while q < e {
                let take = (e - q).min(4);
                let at = |o: usize| if o < take { &left[(q + o) * m + e..(q + o + 1) * m] } else { &left[q * m + e..(q + 1) * m] };
                let coef = |o: usize| if o < take { head[q + o] } else { ZERO };
                subtract(tail, [at(0), at(1), at(2), at(3)], [coef(0), coef(1), coef(2), coef(3)]);
                q += take;
            }
        }
        k = e;
    }
    (swaps, weak)
}

struct Front {
    /// The front's first `p` columns: LU of the pivots, then `L` below them.
    lower: Vec<C64>,
    /// Rows `0..p` of the later columns: `U` to the boundary.
    upper: Vec<C64>,
    swaps: Vec<u32>,
}

pub struct SparseLu<'a> {
    matrix: &'a Csr,
    ordering: &'a Ordering,
    fronts: Vec<Front>,
    weak: bool,
}

impl<'a> SparseLu<'a> {
    pub fn factor(a: &'a Csr, ordering: &'a Ordering) -> Self {
        let scale = a.val.iter().fold(0.0f64, |m, v| m.max(v.norm_sqr())).sqrt();
        let mut local = vec![0; a.n];
        // What each finished front leaves for its parent.
        let mut updates: Vec<Vec<C64>> = vec![vec![]; ordering.nodes.len()];
        let mut fronts = Vec::with_capacity(ordering.nodes.len());
        let mut weak = false;
        for (i, node) in ordering.nodes.iter().enumerate() {
            let p = node.end - node.start;
            let m = p + node.boundary.len();
            let mut f = vec![ZERO; m * m];
            for &(k, r, c) in &node.entries {
                f[r as usize + c as usize * m] += a.val[k as usize];
            }
            (node.start..node.end).for_each(|g| local[g] = g - node.start);
            node.boundary.iter().enumerate().for_each(|(o, &g)| local[g] = p + o);
            for &child in &node.children {
                let update = std::mem::take(&mut updates[child]);
                let map: Vec<usize> = ordering.nodes[child].boundary.iter().map(|&g| local[g]).collect();
                for (c, column) in update.chunks_exact(map.len().max(1)).enumerate() {
                    let target = map[c] * m;
                    column.iter().zip(&map).for_each(|(v, &r)| f[target + r] += *v);
                }
            }
            let (swaps, front_weak) = partial_lu(&mut f, m, p, scale);
            weak |= front_weak;
            let (mut upper, mut update) = (Vec::with_capacity(p * (m - p)), Vec::with_capacity((m - p) * (m - p)));
            for column in f[p * m..].chunks_exact(m) {
                upper.extend_from_slice(&column[..p]);
                update.extend_from_slice(&column[p..]);
            }
            updates[i] = update;
            f.truncate(p * m);
            f.shrink_to_fit();
            fronts.push(Front { lower: f, upper, swaps });
        }
        Self { matrix: a, ordering, fronts, weak }
    }

    fn solve_once(&self, rhs: &mut [C64]) {
        let nodes = &self.ordering.nodes;
        let mut x = vec![ZERO; rhs.len()];
        for (old, &new) in self.ordering.new.iter().enumerate() {
            x[new] = rhs[old];
        }
        for (node, front) in nodes.iter().zip(&self.fronts) {
            let (p, m) = (node.end - node.start, node.end - node.start + node.boundary.len());
            let (head, tail) = x.split_at_mut(node.end);
            let piv = &mut head[node.start..];
            front.swaps.iter().enumerate().for_each(|(j, &r)| piv.swap(j, r as usize));
            for j in 0..p {
                let (v, col) = (piv[j], &front.lower[j * m..(j + 1) * m]);
                for i in j + 1..p {
                    piv[i] -= v * col[i];
                }
                for (o, &g) in node.boundary.iter().enumerate() {
                    tail[g - node.end] -= v * col[p + o];
                }
            }
        }
        for (node, front) in nodes.iter().zip(&self.fronts).rev() {
            let (p, m) = (node.end - node.start, node.end - node.start + node.boundary.len());
            let (head, tail) = x.split_at_mut(node.end);
            let piv = &mut head[node.start..];
            for (o, &g) in node.boundary.iter().enumerate() {
                let (v, col) = (tail[g - node.end], &front.upper[o * p..(o + 1) * p]);
                piv.iter_mut().zip(col).for_each(|(y, u)| *y -= v * *u);
            }
            for j in (0..p).rev() {
                let col = &front.lower[j * m..(j + 1) * m];
                piv[j] /= col[j];
                let v = piv[j];
                for i in 0..j {
                    piv[i] -= v * col[i];
                }
            }
        }
        for (old, &new) in self.ordering.new.iter().enumerate() {
            rhs[old] = x[new];
        }
    }
}

impl Factors for SparseLu<'_> {
    fn solve(&self, rhs: &mut [C64]) {
        let b = rhs.to_vec();
        self.solve_once(rhs);
        if self.weak {
            for _ in 0..REFINEMENT_STEPS {
                let mut r = vec![ZERO; rhs.len()];
                self.matrix.mul(rhs, &mut r);
                r.iter_mut().zip(&b).for_each(|(v, w)| *v = *w - *v);
                self.solve_once(&mut r);
                rhs.iter_mut().zip(&r).for_each(|(v, d)| *v += *d);
            }
        }
    }
}
