# 10. Solving on a 2D mesh

Chapter 5's banded LU was close to ideal in 1D. In 2D, it's the first thing that stops working.

## Fill-in

Factoring a sparse matrix usually creates new nonzeros. Eliminating an unknown connects all its
remaining neighbours to each other, and those new entries are **fill-in**. In 1D, a node's
neighbours are next to it, so little happens.

A 2D mesh of $N$ unknowns is about $\sqrt{N}$ nodes across. Even after chapter 5's reordering, the
band can't be narrower than that, and it fills in completely. So a banded factorization stores about
$N^{3/2}$ entries and costs about $N^2$ operations. At 40,000 unknowns it took 10 s and 1.5 GB, natively
on an Apple M4 Pro, and 160,000 needed 12 GB, far more than this page can give it.

## Nested dissection

The order of elimination decides the fill, so let's choose it from the geometry
([George 1973](20-references.md#george-1973)).

Cut the mesh in half with a line of nodes. Those nodes are a **separator**: nothing on the left
touches anything on the right except through them. If we number the separator last, then eliminating
the left creates no fill on the right, and the two sides can be factored independently. Do the same to
each half, and each quarter, and so on, which makes a tree of ever smaller separators.
<!-- code: dissect "let children" .. "nodes.push((separator" -->

Now the expensive part is only the separators, and most of them are small. For a 2D mesh, the
factorization costs about $N^{3/2}$ operations rather than $N^2$, and storage about $N\log N$. No
ordering of such a mesh can do better than $N^{3/2}$ ([Hoffman et al. 1973](20-references.md#hoffman-1973)).

[`sparse.rs`](../crates/glassworks-core/src/sparse.rs) does this with a few decisions that keep it short:

1. **Split by position, along the longer side.** Whatever on the left touches the right becomes the
   separator. This reads the separator from the matrix rather than the grid, so flipped edges and
   periodic ties need no special care. A periodic tie just makes the separator two lines long.
2. **Thin the separator.** Second-order elements couple nodes two columns apart, so a straight cut
   can leave a separator two columns thick. A node touching nothing further left can join the right.
3. **Reuse the ordering.** It depends only on where the nonzeros are, not their values. So a band
   diagram finds it once for all its wavevectors, and a design for all its iterations.
<!-- code: dissect "let along_x" .. "partition(" -->
<!-- code: dissect "Elements couple nodes" .. "});" -->
<!-- code: Ordering::reused -->

## The multifrontal method

Eliminating a separator's unknowns is a dense problem, because by then they're all coupled to each
other. The **multifrontal** method makes this explicit ([Duff & Reid 1983](20-references.md#duff-reid-1983)). Each separator gets a dense matrix, its
**front**, for its own unknowns and whatever they couple to later. It starts from the matrix entries
that first appear there, and adds the block that each of its two children left behind. Each entry of a
child's block goes to the row and column of the same unknown in the parent's front, which is called
**extend-add**. Eliminating the separator's unknowns from its front leaves a smaller dense block, which
it passes up to its parent.
<!-- code: SparseLu::factor_as "for &child in &node.children" .. "front[target + row] += *value" -->

Take three unknowns in a row, 1 – 3 – 2, with 3 as the separator. The front of 1 holds unknowns 1 and 3:

$$\begin{pmatrix}a_{11} & a_{13}\\ a_{31} & 0\end{pmatrix}.$$

Its (3, 3) entry starts at 0, because $a_{33}$ first appears in the front of 3. Eliminating 1 leaves
$-a_{31}a_{13}/a_{11}$ there, and the front of 2 likewise leaves $-a_{32}a_{23}/a_{22}$. The front of 3
adds both to $a_{33}$, which gives exactly the pivot that plain elimination would. In general, a front
$F$ with its own unknowns first, split into blocks $F_{11}$ for those and $F_{22}$ for the rest, leaves
the **Schur complement** $F_{22} - F_{21}F_{11}^{-1}F_{12}$.
<!-- code: SparseLu::factor_as "let (mut upper" .. "updates[i] = update" -->

So the whole sparse factorization becomes a tree of small dense LUs, and dense loops are about the
fastest code a computer can run. On driven problems in WebAssembly, this was 2 to 6 times faster
than [`faer`](https://docs.rs/faer)'s general sparse LU, which has to find its ordering with no
knowledge of the geometry. `faer` remains in the test suite as an independent reference.

Fronts stay small until the last few, and a front's cost grows as the cube of its size, so most of the
work is in those last fronts. The very last is the top separator, about $\sqrt N$ unknowns, and its
dense front costs about $(\sqrt N)^3 = N^{3/2}$, which is where the estimate above comes from. That's why
thinning separators mattered as much as it did.

A solve with these factors walks the same tree twice. Forward substitution goes up from the leaves:
each front solves with its $L$ for its own unknowns, then subtracts their effect from the boundary
unknowns that later fronts hold. Back substitution comes down from the root: each front subtracts $U$
times the boundary values found above it, then solves for its own unknowns. Each factor entry is used
once, so a solve is far cheaper than factoring.
<!-- code: SparseLu::solve_once -->

## Pivoting, and its price

A pivot can only come from a front's own unknowns, because taking one from elsewhere would change
which fronts need what. That's a real restriction: a front can have only poor pivots even when the
whole matrix is fine. Multifrontal solvers can pass such a pivot up to the parent front, which
changes the fronts' sizes as the factorization runs ([Duff & Reid 1983](20-references.md#duff-reid-1983)).
<!-- code: partial_lu "for i in j..n_pivots" .. "weak |=" -->

`sparse.rs` takes the simpler route of **static pivoting** instead
([Li & Demmel 2003](20-references.md#li-demmel-2003)). It keeps the fronts as they are, replaces a
pivot too small to divide by, and notes any weak one. Each solve is then corrected with two steps of
**iterative refinement**: compute what the answer gets wrong, solve for that, and add it
([Golub & Van Loan 2013](20-references.md#golub-van-loan-2013); [Arioli et al. 1989](20-references.md#arioli-1989)).
The test `sparse_lu_recovers_when_a_part_sits_on_its_own_resonance` builds such a front on purpose,
from a region that sits on a resonance of its own. The error is 0.008 without refinement and
$10^{-12}$ with it.
<!-- code: SparseLu::refined "if self.weak" .. "rhs.iter_mut().zip(&residual)" -->

When $A^T = A$, as it will be with chapter 11's absorbing layer, `SparseLu::factor_symmetric` factors
$A = LDL^T$ instead, with no pivoting and about half the work. If a pivot is unusable, it falls back to
this pivoted LU. Chapter 16 has the details.

## How big can a mesh be?

The limit in a browser is memory, since WebAssembly can address at most 4 GiB per instance. So the solver
refuses meshes above 400,000 unknowns. Chapter 16 shows where that number comes from, and how the
factorization is made faster.
<!-- code: solve2d.rs "Refuses meshes" .. "MAX_DOFS: usize" -->

Next: [open boundaries and leaky modes](11-open-boundaries.md).
