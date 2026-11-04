//! Multifrontal sparse LU for 2D, where the band in `banded` fills in too much. Cutting the
//! domain in two along its longer side leaves two parts joined only through a separator, so
//! each can be eliminated without touching the other. Doing that recursively is nested
//! dissection ([George 1973]). Each separator's unknowns are then eliminated in a dense "front",
//! which also holds whatever its two parts left behind. Chapter 10 of the tutorial has its timings.
//!
//! Pivoting stays inside each front's own unknowns, because a pivot from elsewhere would
//! change which fronts need what. A part walled in by its separator can sit on a resonance of
//! its own, though, and then its front has only a poor pivot, even if the whole matrix is
//! fine. So a front reports any weak pivot, and solves are then refined against the matrix.
//!
//! This file has the method itself. Three refinements that a first reading can skip are in
//! `sparse/`: `blocked.rs` updates a front in cache-sized blocks, `symmetric.rs` is faster
//! when `A^T = A`, and `transposed.rs` solves `A^T x = b`.

use crate::assemble::Csr;
use crate::banded::Factors;
use crate::{C64, ZERO};
use std::cell::RefCell;
use std::rc::Rc;

mod blocked;
mod symmetric;
mod transposed;

use blocked::BLOCK;
use symmetric::partial_ldl;

/// Parts this small are eliminated whole.
const LEAF: usize = 16;
/// Corrections against the matrix that a solve gets when some pivot was weak.
const REFINEMENT_STEPS: usize = 2;

/// One separator, or one leaf part, in the dissection tree: a front to eliminate.
struct Node {
    /// This node's pivots are `start..end` in the new numbering.
    start: usize,
    end: usize,
    /// Nodes whose leftover updates this front adds in. They come earlier in `nodes`.
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
fn dissect(
    ids: Vec<usize>,
    a: &Csr,
    positions: &[(f64, f64)],
    side: &mut [u8],
    nodes: &mut Vec<(Vec<usize>, Vec<usize>)>,
) -> usize {
    if ids.len() <= LEAF {
        nodes.push((ids, vec![]));
        return nodes.len() - 1;
    }
    let extent = |coordinate: &dyn Fn(usize) -> f64| {
        let (low, high) = ids.iter().fold((f64::MAX, f64::MIN), |(low, high), &id| {
            (low.min(coordinate(id)), high.max(coordinate(id)))
        });
        high - low
    };
    let along_x = extent(&|id| positions[id].0) >= extent(&|id| positions[id].1);
    let key = |id: usize| if along_x { positions[id].0 } else { positions[id].1 };
    let mut left = ids;
    let mid = left.len() / 2;
    left.select_nth_unstable_by(mid, |&p, &q| key(p).total_cmp(&key(q)));
    let mut right = left.split_off(mid);
    right.iter().for_each(|&id| side[id] = 2);
    // The separator comes from the matrix, not the grid: whatever on the left touches the
    // right, whether across the cut, over a flipped edge, or through a periodic tie.
    let touches = |id: usize, sides: &[u8], which: u8| {
        a.col[a.ptr[id]..a.ptr[id + 1]].iter().any(|&neighbour| sides[neighbour] == which)
    };
    let (mut separator, left): (Vec<_>, Vec<_>) = left.into_iter().partition(|&id| touches(id, side, 2));
    left.iter().for_each(|&id| side[id] = 1);
    // Elements couple nodes two columns apart, which can make that two columns thick. A node
    // touching nothing further left can join the right, unless there's no left at all: that
    // would only rebuild the same set, and never finish.
    if !left.is_empty() {
        separator.retain(|&id| {
            let needed = touches(id, side, 1);
            if !needed {
                right.push(id);
            }
            needed
        });
    }
    left.iter().chain(&right).for_each(|&id| side[id] = 0);
    let children = [left, right]
        .into_iter()
        .filter(|part| !part.is_empty())
        .map(|part| dissect(part, a, positions, side, nodes))
        .collect();
    nodes.push((separator, children));
    nodes.len() - 1
}

impl Ordering {
    /// `positions` gives each unknown's position. `a` must be structurally symmetric.
    pub fn new(a: &Csr, positions: &[(f64, f64)]) -> Self {
        let n = a.n;
        let mut tree = vec![];
        dissect((0..n).collect(), a, positions, &mut vec![0; n], &mut tree);
        // Children come before parents, so numbering in this order eliminates them first.
        let (mut new, mut owner, mut next) = (vec![0; n], vec![0; n], 0);
        let mut nodes: Vec<Node> = vec![];
        for (i, (pivots, children)) in tree.iter().enumerate() {
            let start = next;
            for &old in pivots {
                new[old] = next;
                owner[next] = i;
                next += 1;
            }
            nodes.push(Node {
                start,
                end: next,
                children: children.clone(),
                boundary: vec![],
                entries: vec![],
            });
        }
        // A front couples to its own entries' later unknowns, and to what its children left.
        let mut seen = vec![usize::MAX; n];
        for i in 0..nodes.len() {
            let end = nodes[i].end;
            let own = tree[i]
                .0
                .iter()
                .flat_map(|&old| a.col[a.ptr[old]..a.ptr[old + 1]].iter().map(|&col| new[col]));
            let inherited = tree[i].1.iter().flat_map(|&child| nodes[child].boundary.iter().copied());
            let mut boundary: Vec<usize> = own
                .chain(inherited)
                .filter(|&unknown| unknown >= end && std::mem::replace(&mut seen[unknown], i) != i)
                .collect();
            boundary.sort_unstable();
            nodes[i].boundary = boundary;
        }
        // Each entry belongs to the first front that holds both its row and its column.
        for row in 0..n {
            for k in a.ptr[row]..a.ptr[row + 1] {
                let (i, j) = (new[row], new[a.col[k]]);
                let node = &mut nodes[owner[i.min(j)]];
                let local = |unknown: usize| {
                    if unknown < node.end {
                        unknown - node.start
                    } else {
                        node.end - node.start + node.boundary.binary_search(&unknown).unwrap()
                    }
                };
                let entry = (k as u32, local(i) as u32, local(j) as u32);
                node.entries.push(entry);
            }
        }
        Self { new, nodes }
    }

    /// The last ordering made if this matrix has its nonzeros, which is cheap to check, and
    /// otherwise a new one. A design's iterations and a band diagram's k-points all reuse it.
    pub fn reused(a: &Csr, positions: impl FnOnce() -> Vec<(f64, f64)>) -> Rc<Self> {
        LAST.with_borrow_mut(|last| match last {
            Some((ptr, col, ordering)) if *ptr == a.ptr && *col == a.col => ordering.clone(),
            _ => {
                let ordering = Rc::new(Self::new(a, &positions()));
                *last = Some((a.ptr.clone(), a.col.clone(), ordering.clone()));
                ordering
            }
        })
    }
}

thread_local! {
    static LAST: RefCell<Option<(Vec<usize>, Vec<usize>, Rc<Ordering>)>> = const { RefCell::new(None) };
}

/// Eliminates the first `n_pivots` of a front's `size` unknowns in place, swapping rows among
/// those `n_pivots`. `front` is column-major. Returns the swaps, and whether any pivot was weak.
fn partial_lu(front: &mut [C64], size: usize, n_pivots: usize, scale: f64) -> (Vec<u32>, bool) {
    let mut swaps = Vec::with_capacity(n_pivots);
    let mut weak = false;
    let mut block_start = 0;
    while block_start < n_pivots {
        let block_end = (block_start + BLOCK).min(n_pivots);
        for j in block_start..block_end {
            let col = &front[j * size..(j + 1) * size];
            let (mut pivot_row, mut best) = (j, 0.0);
            for i in j..n_pivots {
                if col[i].norm_sqr() > best {
                    (pivot_row, best) = (i, col[i].norm_sqr());
                }
            }
            let below = col[n_pivots..].iter().fold(0.0f64, |largest, v| largest.max(v.norm_sqr()));
            weak |= best < 1e-4 * below || best < 1e-16 * scale * scale;
            if pivot_row != j {
                for column in 0..size {
                    front.swap(j + column * size, pivot_row + column * size);
                }
            }
            swaps.push(pivot_row as u32);
            // Refinement corrects for this, unless the whole matrix is singular.
            if front[j + j * size].norm() < 1e-14 * scale {
                front[j + j * size] = C64::new(1e-14 * scale, 0.0);
            }
            let inv_pivot = front[j + j * size].inv();
            front[j * size + j + 1..(j + 1) * size].iter_mut().for_each(|entry| *entry *= inv_pivot);
            let (left, right) = front.split_at_mut((j + 1) * size);
            let multipliers = &left[j * size + j + 1..];
            for col in right.chunks_exact_mut(size).take(block_end - j - 1) {
                let in_pivot_row = col[j];
                col[j + 1..]
                    .iter_mut()
                    .zip(multipliers)
                    .for_each(|(entry, multiplier)| *entry -= in_pivot_row * *multiplier);
            }
        }
        // Update the later columns from this block, all at once to stay in cache.
        blocked::update_after(front, size, block_start, block_end);
        block_start = block_end;
    }
    (swaps, weak)
}

struct Front {
    /// The front's first `n_pivots` columns: LU of the pivots, then `L` below them.
    lower: Vec<C64>,
    /// Rows `0..n_pivots` of the later columns: `U` to the boundary.
    upper: Vec<C64>,
    /// Row swaps from pivoting, as `partial_lu` returns them.
    swaps: Vec<u32>,
}

/// The factors, one front per node of the ordering. The matrix is kept for refinement.
pub struct SparseLu<'a> {
    matrix: &'a Csr,
    ordering: &'a Ordering,
    fronts: Vec<Front>,
    /// Whether any front had a weak pivot, in which case solves are refined.
    weak: bool,
}

impl<'a> SparseLu<'a> {
    pub fn factor(a: &'a Csr, ordering: &'a Ordering) -> Self {
        Self::factor_as(a, ordering, false).expect("pivoting always finds a usable pivot")
    }

    /// Visits the fronts children first. Each gathers its matrix entries and its children's
    /// updates, eliminates its own pivots, and leaves an update for its parent. `None` means
    /// the symmetric factorization met an unusable pivot, and the caller falls back to LU.
    fn factor_as(a: &'a Csr, ordering: &'a Ordering, symmetric: bool) -> Option<Self> {
        let scale = a.val.iter().fold(0.0f64, |largest, v| largest.max(v.norm_sqr())).sqrt();
        // Where each unknown sits in the front being built.
        let mut local_index = vec![0; a.n];
        // What each finished front leaves for its parent.
        let mut updates: Vec<Vec<C64>> = vec![vec![]; ordering.nodes.len()];
        let mut fronts = Vec::with_capacity(ordering.nodes.len());
        let mut weak = false;
        for (i, node) in ordering.nodes.iter().enumerate() {
            let n_pivots = node.end - node.start;
            let size = n_pivots + node.boundary.len();
            let mut front = vec![ZERO; size * size];
            for &(k, row, col) in &node.entries {
                front[row as usize + col as usize * size] += a.val[k as usize];
            }
            (node.start..node.end).for_each(|unknown| local_index[unknown] = unknown - node.start);
            node.boundary
                .iter()
                .enumerate()
                .for_each(|(offset, &unknown)| local_index[unknown] = n_pivots + offset);
            for &child in &node.children {
                let update = std::mem::take(&mut updates[child]);
                let in_parent: Vec<usize> =
                    ordering.nodes[child].boundary.iter().map(|&unknown| local_index[unknown]).collect();
                for (col, column) in update.chunks_exact(in_parent.len().max(1)).enumerate() {
                    let target = in_parent[col] * size;
                    // A symmetric update holds only its lower triangle.
                    let from = if symmetric { col } else { 0 };
                    column[from..]
                        .iter()
                        .zip(&in_parent[from..])
                        .for_each(|(value, &row)| front[target + row] += *value);
                }
            }
            let (swaps, front_weak) = if symmetric {
                (vec![], partial_ldl(&mut front, size, n_pivots, scale)?)
            } else {
                partial_lu(&mut front, size, n_pivots, scale)
            };
            weak |= front_weak;
            // After elimination, the trailing block of `front` is the Schur complement, which the
            // parent adds in.
            let n_boundary = size - n_pivots;
            let (mut upper, mut update) =
                (Vec::with_capacity(n_pivots * n_boundary), Vec::with_capacity(n_boundary * n_boundary));
            for column in front[n_pivots * size..].chunks_exact(size) {
                upper.extend_from_slice(&column[..n_pivots]);
                update.extend_from_slice(&column[n_pivots..]);
            }
            updates[i] = update;
            front.truncate(n_pivots * size);
            front.shrink_to_fit();
            fronts.push(Front { lower: front, upper, swaps });
        }
        Some(Self { matrix: a, ordering, fronts, weak })
    }

    /// Forward substitution through the fronts in elimination order, then back substitution
    /// in reverse, all in the new numbering.
    fn solve_once(&self, rhs: &mut [C64]) {
        let nodes = &self.ordering.nodes;
        let mut x = vec![ZERO; rhs.len()];
        for (old, &new) in self.ordering.new.iter().enumerate() {
            x[new] = rhs[old];
        }
        for (node, front) in nodes.iter().zip(&self.fronts) {
            let n_pivots = node.end - node.start;
            let size = n_pivots + node.boundary.len();
            let (head, tail) = x.split_at_mut(node.end);
            let own = &mut head[node.start..];
            front.swaps.iter().enumerate().for_each(|(j, &row)| own.swap(j, row as usize));
            for j in 0..n_pivots {
                let (xj, col) = (own[j], &front.lower[j * size..(j + 1) * size]);
                for i in j + 1..n_pivots {
                    own[i] -= xj * col[i];
                }
                for (offset, &unknown) in node.boundary.iter().enumerate() {
                    tail[unknown - node.end] -= xj * col[n_pivots + offset];
                }
            }
        }
        for (node, front) in nodes.iter().zip(&self.fronts).rev() {
            let n_pivots = node.end - node.start;
            let size = n_pivots + node.boundary.len();
            let (head, tail) = x.split_at_mut(node.end);
            let own = &mut head[node.start..];
            for (offset, &unknown) in node.boundary.iter().enumerate() {
                let (later, col) =
                    (tail[unknown - node.end], &front.upper[offset * n_pivots..(offset + 1) * n_pivots]);
                own.iter_mut().zip(col).for_each(|(xi, u)| *xi -= later * *u);
            }
            for j in (0..n_pivots).rev() {
                let col = &front.lower[j * size..(j + 1) * size];
                own[j] /= col[j];
                let xj = own[j];
                for i in 0..j {
                    own[i] -= xj * col[i];
                }
            }
        }
        for (old, &new) in self.ordering.new.iter().enumerate() {
            rhs[old] = x[new];
        }
    }

    /// One solve, then refinement against the matrix if any pivot was weak.
    fn refined(
        &self,
        rhs: &mut [C64],
        once: impl Fn(&Self, &mut [C64]),
        mul: impl Fn(&Csr, &[C64], &mut [C64]),
    ) {
        let b = rhs.to_vec();
        once(self, rhs);
        if self.weak {
            for _ in 0..REFINEMENT_STEPS {
                let mut residual = vec![ZERO; rhs.len()];
                mul(self.matrix, rhs, &mut residual);
                residual.iter_mut().zip(&b).for_each(|(ri, bi)| *ri = *bi - *ri);
                once(self, &mut residual);
                rhs.iter_mut().zip(&residual).for_each(|(xi, correction)| *xi += *correction);
            }
        }
    }
}

impl Factors for SparseLu<'_> {
    fn solve(&self, rhs: &mut [C64]) {
        self.refined(rhs, Self::solve_once, Csr::mul);
    }
}
