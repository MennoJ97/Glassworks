# 9. Toward 2D

The solver was written with 2D in mind, though none of this chapter is implemented
yet. It describes what would change, and why vector fields need new ideas.

## Step one: scalar fields on triangles

The easiest 2D problem still has one scalar unknown. A structure that's uniform along
$z$, with $u = E_z$ or $u = \eta_0 H_z$, gives

$$-\nabla\cdot(p\nabla u) - k_0^2 q\,u = f.$$

This is chapter 1's equation with $d/dx$ replaced by $\nabla$, so the weak form is
identical. What changes:

| Module | Change |
|---|---|
| `mesh.rs` | Triangles: `verts_per_cell = 3`. The storage already allows this. A mesher is needed, such as a compiled [Triangle](https://www.cs.cmu.edu/~quake/triangle.html) or a JavaScript Delaunay library. |
| `space.rs` | A new `H1Triangle` implementing `Space`. Its local matrices use a $2\times2$ Jacobian rather than the scale factor $h/2$. |
| `assemble.rs` | Nothing. It already only sees the `Space` trait. |
| `banded.rs` | Works for small meshes, then gets replaced (see below). |
| `eigen.rs` | Nothing, apart from swapping the linear solver. |
| `solve.rs` | Boundaries become edges rather than two end points. |
| Website | A geometry editor for shapes, and a WebGL renderer in place of `FieldPlot`. |

## The solver has to change

In 2D, a mesh of $N$ unknowns has a bandwidth of about $\sqrt{N}$ even after RCM, so
banded LU costs $O(N^2)$. Rough scale: around 10,000 unknowns is still fine, but
100,000 won't be interactive.

Sparse solvers do better using **nested dissection**, cutting the mesh recursively so
elimination creates far less fill. This is where `faer` comes back in. It's already
been shown to compile to WebAssembly here, and `BandedLu` is used in just a handful
of places, so the swap is small. For 2D problems `faer` should be used unless
measurements say otherwise.

## Why vector fields are different

Scalar fields are the easy case. Waveguide modes with full vector fields, or in-plane
fields around sharp corners, need the curl-curl equation:

$$\nabla\times\left(\frac{1}{\mu}\nabla\times\mathbf{E}\right) = k_0^2\,\varepsilon\,\mathbf{E}.$$

The obvious approach is to use the Lagrange elements above for each component.
**This fails, badly.** It produces *spurious modes*: non-physical solutions scattered
among the real ones, at eigenvalues that look perfectly plausible. It's a classic
failure in the history of computational electromagnetics.

The reason is continuity. Lagrange elements force *all* components to be continuous.
But physics only requires the **tangential** component of $\mathbf{E}$ to be continuous
across a material boundary, while the normal component must jump.

## Nédélec edge elements

**Nédélec elements** fix this by matching the physics. Their unknowns live on mesh
*edges*, and each one represents the tangential field along that edge. Tangential
continuity is then automatic, and the normal component is free to jump.

Each edge has a direction, so neighbouring triangles have to agree on a sign. That's
the one change `Space` would need: `cell_dofs` also has to give a sign per dof.

## The null space, and why shift-invert already helps

Any gradient field satisfies $\nabla\times\nabla\phi = 0$. So the curl-curl operator
has a huge null space: thousands of eigenvalues at exactly $\lambda = 0$, none of them
physical.

An eigensolver hunting for the smallest eigenvalues would drown in them. Shift-invert
from chapter 6 avoids this neatly. With $\sigma$ near the modes of interest, the null
space sits far from the shift, at $\theta = -1/\sigma$, a small value, so Lanczos
ignores it. That's why the eigensolver's interface takes a shift, even though 1D
doesn't need this.

## Open boundaries: perfectly matched layers

Chapter 7's open end was exact because in 1D a wave has only one direction to go. In
2D, a boundary condition can't be reflectionless for all angles at once.

The standard solution is a **perfectly matched layer** (PML): an artificial absorbing
region around the domain, equivalent to a complex stretch of the coordinates. The
useful part for us is that PML looks like a complex, anisotropic *material*, so it
slots into the existing assembly. It also makes the eigenproblem non-Hermitian, so
lossy and leaky modes would need `faer`'s general Krylov-Schur eigensolver.

## When the browser isn't the right tool

For large 3D problems, compute in the browser isn't the right tool. Use
[Palace](https://github.com/awslabs/palace) (built on [MFEM](https://mfem.org)), which
runs natively on clusters and supports eigenmode and driven simulations in 3D.
