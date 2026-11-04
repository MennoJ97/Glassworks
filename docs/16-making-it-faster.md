# 16. Making it faster

*Going further. Nothing here changes an answer. It builds on chapters 6 and 10, and covers three standard ways to speed up a sparse direct solver and an eigensolver: blocking for the cache, using symmetry, and not asking for more digits than the answer needs.*

Measure first. On a native cavity run, factoring took 77% of the time, so that's where speed matters most.

## Blocks that stay in cache

A front is a dense matrix, and dense LU spends nearly all its time updating the columns still to come. Done one column at a time, each update reads the whole remaining matrix from memory for very little arithmetic. A **blocked** factorization eliminates a panel of columns first, then updates everything after it from the whole panel at once, so each entry is reused many times while it's in cache. The arithmetic is the same, just reordered. This is how LAPACK factors dense matrices ([Golub & Van Loan 2013](20-references.md#golub-van-loan-2013); [Anderson et al. 1999](20-references.md#anderson-1999)).

`partial_lu` uses panels of 32 columns, and updates later columns four at a time. Other block and leaf sizes were slower here, and a hand-written SIMD kernel gained nothing, since the compiler already vectorizes this loop. The code is in [`sparse/blocked.rs`](../crates/glassworks-core/src/sparse/blocked.rs).
<!-- code: update_after -->

## When the matrix is symmetric

An absorbing layer's operator is complex symmetric, $A^T = A$ (chapter 11), so half of every update to a front repeats the other half. A symmetric factorization, $A = L D L^T$, updates only the lower triangle, for about half the work ([Golub & Van Loan 2013](20-references.md#golub-van-loan-2013); [Duff & Reid 1983](20-references.md#duff-reid-1983)). Here that took a cavity-sized matrix from 54 ms to 30 ms, and a whole cavity design from 16 s to 11 s on the page.
<!-- code: partial_ldl "let inv_pivot = pivot.inv()" .. "*entry -= in_pivot_row" -->

The price is pivoting. Swapping rows alone would break the symmetry. A stable symmetric factorization swaps rows and columns together, and sometimes takes a $2 \times 2$ block as its pivot ([Bunch & Kaufman 1977](20-references.md#bunch-kaufman-1977)). These factors take chapter 10's simpler route instead, static pivoting: they don't pivot at all, and refine each solve against the matrix ([Li & Demmel 2003](20-references.md#li-demmel-2003)). Weak pivots become about three times as common as with the pivoted LU. Over 8,612 solves scanned through resonances, the worst residual was $1.1 \times 10^{-11}$, against $6.4 \times 10^{-12}$ with pivoting. If a pivot is too small to trust, the whole matrix is factored the pivoted way instead.
<!-- code: partial_ldl "let below" .. "weak |=" -->
<!-- code: SparseLu::factor_symmetric -->

The code is in [`sparse/symmetric.rs`](../crates/glassworks-core/src/sparse/symmetric.rs), apart from the refinement, which it shares with chapter 10.

## Asking for the digits you need

A band diagram reads only eigenvalues, and an eigenvalue is far more accurate than its eigenvector. For a Hermitian problem, a Ritz value's error is bounded by the square of its residual divided by the gap to the next eigenvalue ([Parlett 1998](20-references.md#parlett-1998)). So a residual tolerance of $10^{-5}$ still gives $k_0$ to about $10^{-10}$. Lanczos gains accuracy at a roughly geometric rate ([Saad 2011](20-references.md#saad-2011)), so each digit not asked for saves steps. Here it took 62 solves per wavevector down to 47 to 51, and 15 to 23% off a diagram, for a few lines.
<!-- code: BAND_TOL -->

## Where the mesh limit comes from

An eigenmode solve is the heaviest, because Lanczos keeps its vectors alongside the factors, and it took 3.2 GiB at
579,121 unknowns. Counting allocations natively for ten modes, the mesh and matrices hold 0.6 GiB, the
factors 1.2 GiB, and the basis 1.0 GiB: 56 vectors, each stored with $B$ times it. WebAssembly can address at most 4 GiB, so the solver refuses
meshes above 400,000 unknowns, which leaves room to spare.
<!-- code: solve2d.rs "Refuses meshes" .. "const MAX_DOFS" -->

Next: [designing a cavity](17-designing-a-cavity.md).
