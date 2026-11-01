//! Multifrontal sparse LU with geometric nested dissection. Cutting the domain in two along
//! its longer side leaves two parts joined only through a separator, so each part can be
//! eliminated without touching the other. Each separator's unknowns are then eliminated in a
//! dense "front", holding everything its parts left behind.
//!
//! Pivoting stays inside each front's own unknowns, with no delayed pivots. A pivot too small
//! to use is replaced, and iterative refinement corrects for that.

use eigenmenno_core::assemble::Csr;
use eigenmenno_core::C64;

const ZERO: C64 = C64::new(0.0, 0.0);
const LEAF: usize = 16;
const BLOCK: usize = 32;
const TINY: f64 = 1e-14;

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
/// reused for any matrix with the same pattern.
pub struct Analysis {
    pub n: usize,
    /// `new[old]` numbering.
    new: Vec<usize>,
    nodes: Vec<Node>,
}

fn dissect(ids: Vec<usize>, a: &Csr, pos: &[(f64, f64)], side: &mut [u8], nodes: &mut Vec<(Vec<usize>, Vec<usize>)>) -> usize {
    if ids.len() <= LEAF {
        nodes.push((ids, vec![]));
        return nodes.len() - 1;
    }
    let bound = |f: &dyn Fn(usize) -> f64| ids.iter().fold((f64::MAX, f64::MIN), |b, &d| (b.0.min(f(d)), b.1.max(f(d))));
    let (xs, ys) = (bound(&|d| pos[d].0), bound(&|d| pos[d].1));
    let along_x = xs.1 - xs.0 >= ys.1 - ys.0;
    let key = |d: usize| if along_x { pos[d].0 } else { pos[d].1 };
    let mut sorted = ids;
    let mid = sorted.len() / 2;
    sorted.select_nth_unstable_by(mid, |&p, &q| key(p).total_cmp(&key(q)));
    let mut right = sorted.split_off(mid);
    let left = sorted;
    right.iter().for_each(|&d| side[d] = 2);
    // The separator: whatever on the left touches the right, across the cut or a periodic tie.
    let touches = |d: usize, s: &[u8], v: u8| a.col[a.ptr[d]..a.ptr[d + 1]].iter().any(|&c| s[c] == v);
    let (mut sep, left): (Vec<_>, Vec<_>) = left.into_iter().partition(|&d| touches(d, side, 2));
    left.iter().for_each(|&d| side[d] = 1);
    if !left.is_empty() {
        sep.retain(|&d| {
            let needed = touches(d, side, 1);
            if !needed {
                right.push(d);
            }
            needed
        });
    }
    left.iter().chain(&right).for_each(|&d| side[d] = 0);
    let children = [left, right].into_iter().filter(|p| !p.is_empty()).map(|p| dissect(p, a, pos, side, nodes)).collect();
    nodes.push((sep, children));
    nodes.len() - 1
}

impl Analysis {
    /// `pos` gives each unknown's position. `a` must be structurally symmetric.
    pub fn new(a: &Csr, pos: &[(f64, f64)]) -> Self {
        let n = a.n;
        let mut raw = vec![];
        dissect((0..n).collect(), a, pos, &mut vec![0; n], &mut raw);
        // Children come before parents, so numbering in this order eliminates them first.
        let mut new = vec![0; n];
        let mut next = 0;
        let mut nodes: Vec<Node> = vec![];
        let mut owner = vec![0; n];
        for (i, (pivots, children)) in raw.iter().enumerate() {
            let start = next;
            for &d in pivots {
                new[d] = next;
                owner[next] = i;
                next += 1;
            }
            nodes.push(Node { start, end: next, children: children.clone(), boundary: vec![], entries: vec![] });
        }
        let mut seen = vec![usize::MAX; n];
        for i in 0..nodes.len() {
            let end = nodes[i].end;
            let mut set = vec![];
            for &d in &raw[i].0 {
                for &c in &a.col[a.ptr[d]..a.ptr[d + 1]] {
                    if new[c] >= end && seen[new[c]] != i {
                        seen[new[c]] = i;
                        set.push(new[c]);
                    }
                }
            }
            for &child in &raw[i].1 {
                for &b in &nodes[child].boundary {
                    if b >= end && seen[b] != i {
                        seen[b] = i;
                        set.push(b);
                    }
                }
            }
            set.sort_unstable();
            nodes[i].boundary = set;
        }
        // Each entry belongs to the first front that holds both its row and its column.
        for r in 0..n {
            for k in a.ptr[r]..a.ptr[r + 1] {
                let (i, j) = (new[r], new[a.col[k]]);
                let node = &mut nodes[owner[i.min(j)]];
                let local = |g: usize| if g < node.end { g - node.start } else { node.end - node.start + node.boundary.binary_search(&g).unwrap() };
                let e = (k as u32, local(i) as u32, local(j) as u32);
                node.entries.push(e);
            }
        }
        Self { n, new, nodes }
    }
}

impl Analysis {
    /// The original unknowns eliminated under the `depth`-th child chain from the root,
    /// always taking the first child. Only the stress test needs this.
    pub fn subdomain(&self, depth: usize) -> Vec<usize> {
        let mut at = self.nodes.len() - 1;
        for _ in 0..depth {
            at = self.nodes[at].children[0];
        }
        let mut first = at;
        while let Some(&c) = self.nodes[first].children.first() {
            first = c;
        }
        let range = self.nodes[first].start..self.nodes[at].end;
        (0..self.n).filter(|&d| range.contains(&self.new[d])).collect()
    }
}

struct Factored {
    /// The front's first `p` columns: LU of the pivots, then `L` below them.
    lower: Vec<C64>,
    /// Rows `0..p` of the later columns: `U` to the boundary.
    upper: Vec<C64>,
    swaps: Vec<u32>,
}

pub struct Multifrontal<'a> {
    analysis: &'a Analysis,
    fronts: Vec<Factored>,
    /// Pivots replaced for being too small, and those small against their column.
    pub replaced: usize,
    pub weak: usize,
    pub stored: usize,
}

/// `y -= sum(u[q] x[q])`, over up to four columns at once for speed.
#[inline]
fn subtract(y: &mut [C64], x: [&[C64]; 4], u: [C64; 4]) {
    let n = y.len();
    let [x0, x1, x2, x3] = [&x[0][..n], &x[1][..n], &x[2][..n], &x[3][..n]];
    for i in 0..n {
        y[i] -= u[0] * x0[i] + u[1] * x1[i] + u[2] * x2[i] + u[3] * x3[i];
    }
}

/// Eliminates the first `p` of `m` unknowns in place, with row swaps among those `p`.
fn partial_lu(f: &mut [C64], m: usize, p: usize, swaps: &mut Vec<u32>, tiny: f64, replaced: &mut usize, weak: &mut usize) {
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
            if best < 1e-4 * below {
                *weak += 1;
            }
            if r != j {
                for c in 0..m {
                    f.swap(j + c * m, r + c * m);
                }
            }
            swaps.push(r as u32);
            if f[j + j * m].norm() < tiny {
                f[j + j * m] = C64::new(tiny, 0.0);
                *replaced += 1;
            }
            let inv = f[j + j * m].inv();
            f[j * m + j + 1..(j + 1) * m].iter_mut().for_each(|v| *v *= inv);
            let (left, right) = f.split_at_mut((j + 1) * m);
            let l = &left[j * m + j + 1..];
            for c in 0..e - j - 1 {
                let col = &mut right[c * m..(c + 1) * m];
                let u = col[j];
                col[j + 1..].iter_mut().zip(l).for_each(|(y, x)| *y -= u * *x);
            }
        }
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
}

impl<'a> Multifrontal<'a> {
    pub fn factor(analysis: &'a Analysis, a: &Csr) -> Self {
        let scale = a.val.iter().fold(0.0f64, |m, v| m.max(v.norm()));
        let mut local = vec![0; analysis.n];
        let mut updates: Vec<Option<Vec<C64>>> = (0..analysis.nodes.len()).map(|_| None).collect();
        let mut out = Self { analysis, fronts: Vec::with_capacity(analysis.nodes.len()), replaced: 0, weak: 0, stored: 0 };
        for (i, node) in analysis.nodes.iter().enumerate() {
            let p = node.end - node.start;
            let b = node.boundary.len();
            let m = p + b;
            let mut f = vec![ZERO; m * m];
            for &(k, r, c) in &node.entries {
                f[r as usize + c as usize * m] += a.val[k as usize];
            }
            for g in node.start..node.end {
                local[g] = g - node.start;
            }
            for (o, &g) in node.boundary.iter().enumerate() {
                local[g] = p + o;
            }
            // Add what each child left behind.
            for &child in &node.children {
                let update = updates[child].take().unwrap();
                let cb = &analysis.nodes[child].boundary;
                let map: Vec<usize> = cb.iter().map(|&g| local[g]).collect();
                for (c, column) in update.chunks_exact(cb.len().max(1)).enumerate() {
                    let target = map[c] * m;
                    for (v, &r) in column.iter().zip(&map) {
                        f[target + r] += *v;
                    }
                }
            }
            let mut swaps = Vec::with_capacity(p);
            partial_lu(&mut f, m, p, &mut swaps, TINY * scale, &mut out.replaced, &mut out.weak);
            let upper: Vec<C64> = (p..m).flat_map(|c| f[c * m..c * m + p].iter().copied()).collect();
            let update: Vec<C64> = (p..m).flat_map(|c| f[c * m + p..(c + 1) * m].iter().copied()).collect();
            f.truncate(p * m);
            f.shrink_to_fit();
            out.stored += f.len() + upper.len();
            updates[i] = Some(update);
            out.fronts.push(Factored { lower: f, upper, swaps });
        }
        out
    }

    pub fn solve(&self, rhs: &mut [C64]) {
        let mut x = vec![ZERO; self.analysis.n];
        for (old, &new) in self.analysis.new.iter().enumerate() {
            x[new] = rhs[old];
        }
        for (node, front) in self.analysis.nodes.iter().zip(&self.fronts) {
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
        for (node, front) in self.analysis.nodes.iter().zip(&self.fronts).rev() {
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
        for (old, &new) in self.analysis.new.iter().enumerate() {
            rhs[old] = x[new];
        }
    }
}
