# 9. Toward 2D

The solver was written with 2D in mind. Scalar eigenmodes on triangles now exist, in
[`solve2d.rs`](../crates/eigenmenno-core/src/solve2d.rs). Open boundaries and vector
fields don't yet, so those parts of this chapter describe what would change, and why
vector fields need new ideas.

## Step one: scalar fields on triangles

The easiest 2D problem still has one scalar unknown. A structure that's uniform along
$z$, with $u = E_z$ or $u = \eta_0 H_z$, gives

$$-\nabla\cdot(p\nabla u) - k_0^2 q\,u = f.$$

This is chapter 1's equation with $d/dx$ replaced by $\nabla$, so the weak form is
identical. Here's what changed:

| Module | Change |
|---|---|
| `mesh.rs` | Triangles, with `verts_per_cell = 3`. A structured grid has its nodes snapped onto shape boundaries, which needs no mesher and keeps opposite edges matched for periodic ties. A general mesher such as [Triangle](https://www.cs.cmu.edu/~quake/triangle.html) ([Shewchuk 1996](10-references.md#shewchuk-1996)) would handle harder geometry. |
| `space.rs` | A new `H1Triangle` implementing `Space`. Its local matrices use gradients on the triangle rather than the scale factor $h/2$. |
| `assemble.rs` | Assembly is unchanged, because it only sees the `Space` trait. `DofMap` takes many periodic ties rather than one. |
| `banded.rs` | Kept for 1D. 2D uses `SparseLu`, behind the same small `Factors` trait (see below). |
| `eigen.rs` | Nothing, apart from taking any `Factors`. |
| `solve2d.rs` | New: the 2D problems, with boundaries as edges rather than two end points. |
| Website | A table of shapes, and a triangle renderer alongside `FieldPlot`. |

The curved boundary limits accuracy. Straight cells cut the corners of a circle, so a
rod's band edges converge at second order, even with second-order elements. On a staircase
mesh with no snapping, they stop converging altogether.

## The solver has to change

In 2D, a mesh of $N$ unknowns has a bandwidth of about $\sqrt{N}$ even after RCM, so
banded LU costs $O(N^2)$. For 40,000 unknowns it took 9.7 s and 1.5 GB, and by 160,000 it
needed 12 GB, more than a browser tab can address.

Sparse solvers do better using **nested dissection**, cutting the mesh recursively so
elimination creates far less fill
([George 1973](10-references.md#george-1973)). `faer`'s sparse LU limits fill in the same spirit, and it took 0.31 s
on those 40,000 unknowns, about 30 times faster. So
[`sparse.rs`](../crates/eigenmenno-core/src/sparse.rs) wraps `faer` for 2D, at the cost of
more to download.

## Why vector fields are different

Scalar fields are the easy case. Waveguide modes with full vector fields, or in-plane
fields around sharp corners, need the curl-curl equation:

$$\nabla\times\left(\frac{1}{\mu}\nabla\times\mathbf{E}\right) = k_0^2\,\varepsilon\,\mathbf{E}.$$

The obvious approach is to use the Lagrange elements above for each component.
**This fails, badly.** It produces *spurious modes*: non-physical solutions scattered
among the real ones, at eigenvalues that look perfectly plausible. It's a classic
failure in the history of computational electromagnetics
([Sun et al. 1995](10-references.md#sun-1995)).

The cause is the curl's null space. Any gradient field satisfies
$\nabla\times\nabla\phi = 0$, so every gradient solves the curl-curl equation with
$k_0 = 0$. These static fields aren't resonances, but a discretization still has to
represent them consistently ([Sun et al. 1995](10-references.md#sun-1995)).

Lagrange elements can't. Take a continuous, piecewise polynomial $\phi$. Across a cell
boundary, the tangential part of its gradient is continuous but the normal part jumps.
A field with every component continuous can't do that, so most of these gradients
aren't in the Lagrange space. Their poor approximations don't stay at zero: they turn
up at nonzero eigenvalues, mixed in with the real modes. This happens even in a uniform
material, so it isn't only a problem at interfaces.

Interfaces do add a second problem. Physics only requires the **tangential** component
of $\mathbf{E}$ to be continuous across a material boundary, while the normal component
must jump, and Lagrange elements force both to be continuous.

## Nédélec edge elements

**Nédélec elements** fix both problems with one choice of continuity
([Nédélec 1980](10-references.md#nedelec-1980)). Their unknowns live on mesh *edges*,
and each one represents the tangential field along that edge. Tangential continuity is
then automatic, and the normal component is free to jump
([Jin 2014](10-references.md#jin-2014)).

That's the same continuity a gradient has. The gradient of every continuous, piecewise
polynomial of matching order lies exactly in the Nédélec space, so the static fields
are represented exactly and stay at zero ([Monk 2003](10-references.md#monk-2003)).

Each edge has a direction, so neighbouring triangles have to agree on a sign. That's
the one change `Space` would need: `cell_dofs` also has to give a sign per dof.

## The null space, and why shift-invert already helps

Representing those gradients exactly has a cost: the matrix has a huge null space, with
thousands of eigenvalues at exactly $\lambda = 0$, none of them a resonance.

An eigensolver hunting for the smallest eigenvalues would drown in them. Shift-invert
from chapter 6 avoids this neatly. With $\sigma$ near the modes of interest, the null
space sits far from the shift, at $\theta = -1/\sigma$, a small value, so Lanczos
ignores it. That's why the eigensolver's interface takes a shift, even though 1D
doesn't need this.

## Open boundaries: perfectly matched layers

Chapter 7's open end was exact because in 1D a wave has only one direction to go. In
2D, a local boundary condition can't be reflectionless for all angles at once
([Engquist & Majda 1977](10-references.md#engquist-majda-1977); [Johnson 2021](10-references.md#johnson-2021)).

This isn't implemented yet. The standard solution is a **perfectly matched layer** (PML): an artificial absorbing
region around the domain, equivalent to a complex stretch of the coordinates
([Bérenger 1994](10-references.md#berenger-1994); [Chew & Weedon 1994](10-references.md#chew-weedon-1994)). PML looks like a complex, anisotropic *material*
([Sacks et al. 1995](10-references.md#sacks-1995); [Johnson 2021](10-references.md#johnson-2021)). That nearly fits the existing assembly, but the anisotropy means
the stiffness has to split into separate $x$ and $y$ parts. It also makes the eigenproblem
non-Hermitian. Leaky modes could use a general eigensolver such as Krylov-Schur
([Stewart 2002](10-references.md#stewart-2002)), but a driven solve at complex frequency
gets the local density of states directly, with no eigensolver at all.

## When the browser isn't the right tool

For large 3D problems, compute in the browser isn't the right tool. Use
[Palace](https://github.com/awslabs/palace) (built on [MFEM](https://mfem.org), [Anderson et al. 2021](10-references.md#anderson-2021)), which
runs natively on clusters and supports eigenmode and driven simulations in 3D.

Next: [references](10-references.md).
