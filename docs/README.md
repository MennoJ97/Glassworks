# Learn the finite element method by reading a small solver

This tutorial explains the finite element method (FEM) using the source code
of Glassworks itself. Each chapter takes one idea, derives it, and points at the few lines
that implement it. The essentials are about 2,500 lines of Rust, so you can read
all of them. The 1D chapters cover about 900 of those. Another 500 lines turn requests
into structs, check them and package the results, and you can look those up when needed.

You need calculus (integration by parts) and basic linear algebra (matrices,
eigenvalues). Knowing Maxwell's equations helps but isn't required, because
chapter 1 reduces them to one ordinary differential equation.

On the website, each chapter's code sits beside it. Point at a name in `code` type, or at a
passage or equation marked in the margin, and the code shows the lines that do it. As you
scroll, the code follows the passage you're reading. If you're reading the Markdown files
instead, on GitHub for example, follow the links to the source files.

This teaches the method, and it doesn't design devices. Its 2D results use bulk-like indices rather than a real chip's effective index, and they can't see light leaking out of the plane.

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
16. [Making it faster](16-making-it-faster.md): blocking, symmetry, and asking only for the digits you need
17. [Designing a cavity](17-designing-a-cavity.md): optimizing the LDOS, and holding a resonance still while the optimizer works
18. [When numbers lie](18-when-numbers-lie.md): checking a result against an invariance and against energy conservation
19. [Designing in the other polarization](19-designing-in-hz.md): bends and cavities in $H_z$, where a chip's light usually is

The [references](20-references.md) list every source cited in the chapters.

## Where the code lives

The files follow the chapters. Whatever only a Going further chapter needs has a file of its own, so the essentials read straight through.

Some files are for looking up rather than reading: [`problem.rs`](../crates/glassworks-core/src/problem.rs) defines the requests, [`output.rs`](../crates/glassworks-core/src/output.rs) what a solve returns, and `solve1d/validate.rs` and [`solve2d/validate.rs`](../crates/glassworks-core/src/solve2d/validate.rs) check a request before any work starts.

| Chapter | Source file |
|---|---|
| 1, 7 | [`solve1d.rs`](../crates/glassworks-core/src/solve1d.rs) |
| 3, 9 | [`mesh.rs`](../crates/glassworks-core/src/mesh.rs), [`space.rs`](../crates/glassworks-core/src/space.rs) |
| 4 | [`assemble.rs`](../crates/glassworks-core/src/assemble.rs) |
| 5 | [`banded.rs`](../crates/glassworks-core/src/banded.rs) |
| 6 | [`eigen.rs`](../crates/glassworks-core/src/eigen.rs) |
| 8 | [`tests/analytic.rs`](../crates/glassworks-core/tests/analytic.rs), [`tests/cross_check.rs`](../crates/glassworks-core/tests/cross_check.rs), and for 2D [`tests/analytic2d.rs`](../crates/glassworks-core/tests/analytic2d.rs) |
| 9 to 13 | [`solve2d.rs`](../crates/glassworks-core/src/solve2d.rs), which sends each request to one of the files below |
| 9 | [`solve2d/closed.rs`](../crates/glassworks-core/src/solve2d/closed.rs) |
| 10 | [`sparse.rs`](../crates/glassworks-core/src/sparse.rs) |
| 11 | [`solve2d/open.rs`](../crates/glassworks-core/src/solve2d/open.rs), [`solve2d/leaky.rs`](../crates/glassworks-core/src/solve2d/leaky.rs) |
| 12 | [`solve2d/driven.rs`](../crates/glassworks-core/src/solve2d/driven.rs) |
| 13 | [`solve2d/adjoint.rs`](../crates/glassworks-core/src/solve2d/adjoint.rs), and on the website, `web/js/optimizer.js` and `web/js/design.js`. The website's 1D mirror designer uses [`solve1d/gradient.rs`](../crates/glassworks-core/src/solve1d/gradient.rs), the same adjoint for layers |

| Going further | Source file |
|---|---|
| 15 | [`solve2d/grating.rs`](../crates/glassworks-core/src/solve2d/grating.rs), [`sparse/transposed.rs`](../crates/glassworks-core/src/sparse/transposed.rs), [`tests/gratings.rs`](../crates/glassworks-core/tests/gratings.rs) |
| 16 | [`sparse/symmetric.rs`](../crates/glassworks-core/src/sparse/symmetric.rs), [`sparse/blocked.rs`](../crates/glassworks-core/src/sparse/blocked.rs) |
| 17 | `web/js/cavity-design.js`, on the website rather than in the solver |
| 19 | [`solve2d/dipole.rs`](../crates/glassworks-core/src/solve2d/dipole.rs), [`tests/hz_design.rs`](../crates/glassworks-core/tests/hz_design.rs), and on the website, `web/js/hz-design.js` |
