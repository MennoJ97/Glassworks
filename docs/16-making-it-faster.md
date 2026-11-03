# 16. Making it faster

*Going further. Nothing here changes an answer. It builds on chapters 6 and 10, and is about how to choose what to speed up. The lesson that keeps returning: measure before you add.*

## An optimization that didn't help: Cholesky

A band diagram's matrix is Hermitian and positive definite after the shift, which allows Cholesky,
with no pivoting and half the storage. It's the textbook move, so it was tried: it halved factoring
but saved only 4 to 7% overall. A band diagram spends most of its time on back-substitutions and
Lanczos, and $L$ has as many entries as $U$. So it was left out. Thinning the separators did far more,
for every problem. Measure before you add.


## Blocks that stay in cache

Factoring took 77% of a native cavity run, so that's where speed matters. `partial_lu` eliminates a
front's columns 32 at a time, then updates every later column from all 32 at once, four at a time,
while they're still in cache. The arithmetic is the same as updating after each column. Other block
and leaf sizes were slower, and a hand-written SIMD kernel gained nothing, since the compiler already
vectorizes this loop. The code is in [`sparse/blocked.rs`](../crates/glassworks-core/src/sparse/blocked.rs).
<!-- code: update_after -->

## When the matrix is symmetric

An absorbing layer's operator is complex symmetric, $A^T = A$, so half of every update to a front repeats the other half. Factoring it as $A = L D L^T$ updates only the lower triangle, which took 54 ms down to 30 ms on a cavity-sized matrix, and a whole cavity design from 16 s to 11 s on the page.
<!-- code: partial_ldl "let inv = d.inv()" .. "*y -= u" -->

The price is pivoting. Swapping rows alone would break the symmetry, so these factors don't pivot at all, and weak pivots become about three times as common. Solves are then refined against the matrix, as for pivoting inside fronts. Over 8,612 solves scanned through resonances, the worst residual was $1.1 \times 10^{-11}$, against $6.4 \times 10^{-12}$ with pivoting. If a pivot is too small to trust, the whole matrix is factored the pivoted way instead.
<!-- code: partial_ldl "let below" .. "weak |=" -->
<!-- code: SparseLu::factor_symmetric -->

The code is in [`sparse/symmetric.rs`](../crates/glassworks-core/src/sparse/symmetric.rs), apart from the method it refines.

## A band diagram that asks for too much

Warm-starting Lanczos from the previous wavevector looked just as obvious as Cholesky, and it did nothing: 62 solves per
wavevector before and after. Lanczos gains digits geometrically, so a good start saves only the first one or two of the ten it
was asked for. The better question was why ten. A band diagram reads only eigenvalues, whose error goes as the square of the
residual ([Parlett 1998](20-references.md#parlett-1998)), so a tolerance of $10^{-5}$ still gives $k_0$ to about $10^{-10}$. That takes 47 to 51 solves per wavevector and 15 to 23%
off a diagram, for a few lines. Adding the warm start on top saved only 2 to 9% more of the solves, so it was left out too.
<!-- code: BAND_TOL -->

## Where the mesh limit comes from

An eigenmode solve is the heaviest, because Lanczos keeps its vectors alongside the factors, and it took 3.2 GiB at
579,121 unknowns. Counting allocations natively for ten modes, the mesh and matrices hold 0.6 GiB, the
factors 1.2 GiB, and the basis 1.0 GiB: 56 vectors, each stored with $B$ times it. So the solver refuses
meshes above 400,000 unknowns, which leaves room to spare.
<!-- code: solve2d.rs "Refuses meshes" .. "const MAX_DOFS" -->

Next: [designing a cavity](17-designing-a-cavity.md).
