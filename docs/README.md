# Learn the finite element method by reading a small solver

This tutorial explains the finite element method (FEM) using Eigenmenno's own
source code. Each chapter takes one idea, derives it, and points at the few lines
that implement it. The essentials are about 3,000 lines of Rust, so you can read
all of them. The 1D chapters cover about 1,200 of those.

You need calculus (integration by parts) and basic linear algebra (matrices,
eigenvalues). Knowing Maxwell's equations helps but isn't required, because
chapter 1 reduces them to one ordinary differential equation.

## The essentials, in reading order

These build on each other, and together they take you from Maxwell's equations to a solver that designs its own structures.

1. [From Maxwell's equations to one equation](01-from-maxwell.md): the physics
2. [The weak form](02-weak-form.md): the idea FEM is built on
3. [Elements](03-elements.md): turning a function into a list of numbers
4. [Assembly and boundary conditions](04-assembly.md): building the matrices
5. [Solving linear systems](05-linear-systems.md): banded LU and reordering
6. [Eigenproblems](06-eigenproblems.md): finding modes with Lanczos
7. [Driven problems](07-driven.md): waves, ports, reflection and transmission
8. [Is the answer right?](08-verification.md): testing against exact results
9. [Fields on triangles](09-triangles.md): second-order elements, fitted meshes and band diagrams
10. [Solving on a 2D mesh](10-sparse-lu.md): fill-in, nested dissection and multifrontal LU
11. [Open boundaries and leaky modes](11-open-boundaries.md): absorbing layers, LDOS and Q
12. [Waveguides and bends](12-waveguides.md): ports, monitors and scale invariance
13. [Inverse design](13-inverse-design.md): adjoints, and letting the solver draw the structure

<a id="going-further"></a>

## Going further

Each of these stands on its own, so read any that interest you, or none.

14. [Beyond scalar fields](14-vector-fields.md): why vector fields are harder
15. [Diffraction gratings](15-gratings.md): Bloch ties and absorbing layers at once, and transposed solves
16. [Making it faster](16-making-it-faster.md): what was tried, and which of it was worth keeping
17. [Designing a cavity](17-designing-a-cavity.md): holding a resonance still while the optimizer works
18. [When numbers lie](18-when-numbers-lie.md): two results that looked wrong, or looked fine and weren't

The [references](19-references.md) list every source cited in the chapters.

## Where the code lives

The files follow the chapters. Whatever only a Going further chapter needs has a file of its own, so the essentials read straight through.

| Chapter | Source file |
|---|---|
| 1, 7 | [`solve1d.rs`](../crates/eigenmenno-core/src/solve1d.rs), which fills in [`output.rs`](../crates/eigenmenno-core/src/output.rs) |
| 3, 9 | [`mesh.rs`](../crates/eigenmenno-core/src/mesh.rs), [`space.rs`](../crates/eigenmenno-core/src/space.rs) |
| 4 | [`assemble.rs`](../crates/eigenmenno-core/src/assemble.rs) |
| 5 | [`banded.rs`](../crates/eigenmenno-core/src/banded.rs) |
| 6 | [`eigen.rs`](../crates/eigenmenno-core/src/eigen.rs) |
| 8 | [`tests/analytic.rs`](../crates/eigenmenno-core/tests/analytic.rs), [`tests/cross_check.rs`](../crates/eigenmenno-core/tests/cross_check.rs), and for 2D [`tests/analytic2d.rs`](../crates/eigenmenno-core/tests/analytic2d.rs) |
| 9 to 13 | [`solve2d.rs`](../crates/eigenmenno-core/src/solve2d.rs), which sends each request to one of the files below |
| 9 | [`solve2d/closed.rs`](../crates/eigenmenno-core/src/solve2d/closed.rs) |
| 10 | [`sparse.rs`](../crates/eigenmenno-core/src/sparse.rs) |
| 11 | [`solve2d/open.rs`](../crates/eigenmenno-core/src/solve2d/open.rs), [`solve2d/leaky.rs`](../crates/eigenmenno-core/src/solve2d/leaky.rs) |
| 12 | [`solve2d/driven.rs`](../crates/eigenmenno-core/src/solve2d/driven.rs) |
| 13 | [`solve2d/adjoint.rs`](../crates/eigenmenno-core/src/solve2d/adjoint.rs), and on the website, `web/js/optimizer.js` and `web/js/design.js` |

| Going further | Source file |
|---|---|
| 15 | [`solve2d/grating.rs`](../crates/eigenmenno-core/src/solve2d/grating.rs), [`sparse/transposed.rs`](../crates/eigenmenno-core/src/sparse/transposed.rs), [`tests/gratings.rs`](../crates/eigenmenno-core/tests/gratings.rs) |
| 16 | [`sparse/symmetric.rs`](../crates/eigenmenno-core/src/sparse/symmetric.rs) |
| 17 | `web/js/cavity-design.js`, on the website rather than in the solver |
