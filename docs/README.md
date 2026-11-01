# Learn the finite element method by reading a small solver

This tutorial explains the finite element method (FEM) using Eigenmenno's own
source code. Each chapter takes one idea, derives it, and points at the few lines
that implement it. The whole solver is about 1,900 lines of Rust, so you can read
all of it. The 1D chapters cover about 1,100 of those.

You need calculus (integration by parts) and basic linear algebra (matrices,
eigenvalues). Knowing Maxwell's equations helps but isn't required, because
chapter 1 reduces them to one ordinary differential equation.

## Chapters, in reading order

1. [From Maxwell's equations to one equation](01-from-maxwell.md): the physics
2. [The weak form](02-weak-form.md): the idea FEM is built on
3. [Elements](03-elements.md): turning a function into a list of numbers
4. [Assembly and boundary conditions](04-assembly.md): building the matrices
5. [Solving linear systems](05-linear-systems.md): banded LU and reordering
6. [Eigenproblems](06-eigenproblems.md): finding modes with Lanczos
7. [Driven problems](07-driven.md): waves, ports, reflection and transmission
8. [Is the answer right?](08-verification.md): testing against exact results
9. [Toward 2D](09-toward-2d.md): what changes, and why vector fields are harder

The [references](10-references.md) list every source cited in the chapters.

## Where the code lives

| Chapter | Source file |
|---|---|
| 1, 7 | [`solve.rs`](../crates/eigenmenno-core/src/solve.rs) |
| 3 | [`mesh.rs`](../crates/eigenmenno-core/src/mesh.rs), [`space.rs`](../crates/eigenmenno-core/src/space.rs) |
| 4 | [`assemble.rs`](../crates/eigenmenno-core/src/assemble.rs) |
| 5 | [`banded.rs`](../crates/eigenmenno-core/src/banded.rs) |
| 6 | [`eigen.rs`](../crates/eigenmenno-core/src/eigen.rs) |
| 8 | [`tests/analytic.rs`](../crates/eigenmenno-core/tests/analytic.rs), [`tests/cross_check.rs`](../crates/eigenmenno-core/tests/cross_check.rs) |
