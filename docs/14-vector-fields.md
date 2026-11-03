# 14. Beyond scalar fields

*Going further. This builds on chapter 9.*

Everything in the last five chapters had one scalar unknown, which is what makes the mesh, the solver
and the design tools so simple. This solver doesn't have vector fields yet. This chapter describes what
would change, and why they need new ideas.

## Why vector fields are different

Waveguide modes with full vector fields, or in-plane fields around sharp corners, need the curl-curl equation:

$$\nabla\times\left(\frac{1}{\mu}\nabla\times\mathbf{E}\right) = k_0^2\,\varepsilon\,\mathbf{E}.$$

The obvious approach is to use the Lagrange elements from chapter 9 for each component.
**This fails, badly.** It produces *spurious modes*: non-physical solutions scattered among the real ones,
at eigenvalues that look perfectly plausible. It's a classic failure in the history of computational
electromagnetics ([Sun et al. 1995](20-references.md#sun-1995)).

The cause is the curl's null space. Any gradient field satisfies $\nabla\times\nabla\phi = 0$, so every gradient
solves the curl-curl equation with $k_0 = 0$. These static fields aren't resonances, but a discretization still
has to represent them consistently ([Sun et al. 1995](20-references.md#sun-1995)).

Lagrange elements can't. Take a continuous, piecewise polynomial $\phi$. Across a cell boundary, the tangential
part of its gradient is continuous but the normal part jumps. A field with every component continuous can't do
that, so most of these gradients aren't in the Lagrange space. Their poor approximations don't stay at zero: they
turn up at nonzero eigenvalues, mixed in with the real modes. This happens even in a uniform material, so it isn't
only a problem at interfaces.

Interfaces do add a second problem. Physics only requires the **tangential** component of $\mathbf{E}$ to be
continuous across a material boundary, while the normal component must jump, and Lagrange elements force both to be
continuous.

## Nédélec edge elements

**Nédélec elements** fix both problems with one choice of continuity
([Nédélec 1980](20-references.md#nedelec-1980)). Their unknowns live on mesh *edges*, and each one represents the
tangential field along that edge. Tangential continuity is then automatic, and the normal component is free to
jump ([Jin 2014](20-references.md#jin-2014)).

That's the same continuity a gradient has. The gradient of every continuous, piecewise polynomial of matching order
lies exactly in the Nédélec space, so the static fields are represented exactly and stay at zero
([Monk 2003](20-references.md#monk-2003)).

Each edge has a direction, so neighbouring triangles have to agree on a sign. That's the one change `Space` would need:
`cell_dofs` also has to give a sign per dof.

## The null space, and why shift-invert already helps

Representing those gradients exactly has a cost: the matrix has a huge null space, with thousands of eigenvalues at
exactly $\lambda = 0$, none of them a resonance.

An eigensolver hunting for the smallest eigenvalues would drown in them. Shift-invert from chapter 6 avoids this
neatly. With $\sigma$ near the modes of interest, the null space sits far from the shift, at $\theta = -1/\sigma$, a
small value, so Lanczos ignores it. That's why the eigensolver's interface, `eigs_near`, takes a shift, even though 1D
doesn't need this.
<!-- code: eigs_near "pub fn eigs_near" .. "sub_scaled" -->

## When the browser isn't the right tool

Large 3D problems don't fit in a browser. For those, use
[Palace](https://github.com/awslabs/palace) (built on [MFEM](https://mfem.org), [Anderson et al. 2021](20-references.md#anderson-2021)), which
runs natively on clusters and supports eigenmode and driven simulations in 3D.

Next: [diffraction gratings](15-gratings.md).
