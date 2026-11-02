# 10. Solving on a 2D mesh

Chapter 5's banded LU was close to ideal in 1D. In 2D, it's the first thing that stops working.

## Fill-in

Factoring a sparse matrix usually creates new nonzeros. Eliminating an unknown connects all its
remaining neighbours to each other, and those new entries are **fill-in**. In 1D, a node's
neighbours are next to it, so little happens.

A 2D mesh of $N$ unknowns is about $\sqrt{N}$ nodes across. Even after chapter 5's reordering, the
band can't be narrower than that, and it fills in completely. So a banded factorization stores about
$N^{3/2}$ entries and costs about $N^2$ operations. At 40,000 unknowns it took 9.7 s and 1.5 GB, and
160,000 needed 12 GB, far more than this page can give it.

## Nested dissection

The order of elimination decides the fill, so let's choose it from the geometry
([George 1973](15-references.md#george-1973)).

Cut the mesh in half with a line of nodes. Those nodes are a **separator**: nothing on the left
touches anything on the right except through them. If we number the separator last, then eliminating
the left creates no fill on the right, and the two sides can be factored independently. Do the same to
each half, and each quarter, and so on, which makes a tree of ever smaller separators.

Now the expensive part is only the separators, and most of them are small. For a 2D mesh, the
factorization costs about $N^{3/2}$ operations rather than $N^2$, and storage about $N\log N$
([George 1973](15-references.md#george-1973); [Davis 2006](15-references.md#davis-2006)).

[`sparse.rs`](../crates/eigenmenno-core/src/sparse.rs) does this with a few decisions that keep it short:

1. **Split by position, along the longer side.** Whatever on the left touches the right becomes the
   separator. This reads the separator from the matrix rather than the grid, so flipped edges and
   periodic ties need no special care. A periodic tie just makes the separator two lines long.
2. **Thin the separator.** Second-order elements couple nodes two columns apart, so a straight cut
   can leave a separator two columns thick. A node touching nothing further left can join the right.
3. **Reuse the ordering.** It depends only on where the nonzeros are, not their values. So a band
   diagram finds it once for all its wavevectors, and a design for all its iterations.

## The multifrontal method

Eliminating a separator's unknowns is a dense problem, because by then they're all coupled to each
other. The **multifrontal** method makes this explicit. Each separator gets a dense matrix, its
**front**, for its own unknowns and whatever they couple to later. It starts from the matrix entries
that first appear there, and adds the block that each of its two children left behind. Eliminating the
separator's unknowns from its front leaves a smaller dense block, which it passes up to its parent.

So the whole sparse factorization becomes a tree of small dense LUs, and dense loops are about the
fastest code a computer can run. On driven problems in WebAssembly, this was 2 to 6 times faster
than [`faer`](https://docs.rs/faer)'s general sparse LU, which has to find its ordering with no
knowledge of the geometry. `faer` remains in the test suite as an independent reference.

Fronts stay small until the last few, and a front's cost grows as the cube of its size, so most of the
work is in those last fronts. That's why thinning separators mattered as much as it did.

## Pivoting, and its price

A pivot can only come from a front's own unknowns, because taking one from elsewhere would change
which fronts need what. That's a real restriction. A part of the mesh walled in by its separator is a
small cavity, and it can happen to sit on a resonance of its own. Its front then has only poor pivots,
even when the whole matrix is fine.

So `sparse.rs` notes any weak pivot, and then corrects each solve with two steps of **iterative
refinement**: compute what the answer gets wrong, solve for that, and add it. The test
`sparse_lu_recovers_when_a_part_sits_on_its_own_resonance` builds this case on purpose. The error is
0.008 without refinement and $10^{-12}$ with it.

## An optimization that didn't help

A band diagram's matrix is Hermitian and positive definite after the shift, which allows Cholesky,
with no pivoting and half the storage. It's the textbook move, so it was tried: it halved factoring
but saved only 4 to 7% overall. A band diagram spends most of its time on back-substitutions and
Lanczos, and $L$ has as many entries as $U$. So it was left out. Thinning the separators did far more,
for every problem. Measure before you add.

## How big can a mesh be?

The limit in a browser is memory. WebAssembly can address at most 4 GiB per instance. An eigenmode
solve is the heaviest, because Lanczos keeps its vectors alongside the factors, and it took 3.2 GiB at
579,121 unknowns. So the solver refuses meshes above 400,000 unknowns, which leaves room to spare.

Next: [open boundaries and leaky modes](11-open-boundaries.md).
