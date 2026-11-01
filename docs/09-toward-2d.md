# 9. Toward 2D

The solver was written with 2D in mind. Scalar modes, open scattering, leaky cavities and
waveguides on triangles now exist, in [`solve2d.rs`](../crates/eigenmenno-core/src/solve2d.rs). Vector fields don't
yet, so that part of this chapter describes what would change, and why they need new
ideas.

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
elimination creates far less fill ([George 1973](10-references.md#george-1973)).
[`sparse.rs`](../crates/eigenmenno-core/src/sparse.rs) does this in about 290 lines:

1. Split the unknowns in two by position, along the longer side. Whatever on the left touches
   the right becomes a **separator**. The two parts then don't interact, so each can be
   eliminated without touching the other. Reading separators from the matrix, not the grid,
   means flipped edges and periodic ties need no special care. A tie just makes the separator
   two lines long.
2. Number each part first and its separator last, recursively. This **ordering** depends only
   on where the nonzeros are, so a band diagram computes it once for all its wavevectors.
3. Eliminate each separator in a dense **front**, which holds its own matrix entries, plus
   whatever its two parts left behind. Eliminating it leaves a smaller dense block for
   its parent. This is the **multifrontal** method.

Fronts stay small until near the end, and dense loops run fast. Around 40,000 unknowns, this
factored in about 0.5 s in WebAssembly, against about 3.5 to 5.6 s for
[`faer`](https://docs.rs/faer)'s general sparse LU, which has to find its ordering with no
knowledge of the geometry. `faer`'s Cholesky is still the faster on band diagrams, where the
matrix is positive definite, and needs no pivoting and half the storage. This solver doesn't
use that yet, which makes band diagrams about 1.4 times slower than with `faer`.

The cost is in the pivoting. A pivot can only come from a front's own unknowns, because using
one from elsewhere would change the fronts that follow. A part walled in by its separator can
sit on a resonance of its own, and then its front has only a poor pivot, even though the whole
matrix is fine. `sparse.rs` notes any weak pivot, and then corrects each solve with two steps of
**iterative refinement**: find what the answer gets wrong, solve for that, and add it.
`sparse_lu_recovers_when_a_part_sits_on_its_own_resonance` builds this case on purpose: the
error is 0.12 without refinement and 9e-13 with it.

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

The standard solution is a **perfectly matched layer** (PML): an artificial absorbing
region around the domain, equivalent to a complex stretch of the coordinates
([Bérenger 1994](10-references.md#berenger-1994); [Chew & Weedon 1994](10-references.md#chew-weedon-1994)). PML looks like a complex, anisotropic *material*
([Sacks et al. 1995](10-references.md#sacks-1995); [Johnson 2021](10-references.md#johnson-2021)). With a stretch $s_x$, $s_y$ along each axis, the weights become

$$a_x = p\,\frac{s_y}{s_x}, \qquad a_y = p\,\frac{s_x}{s_y}, \qquad b = -k_0^2\,q\,s_x s_y.$$

So the stiffness had to split into separate $x$ and $y$ parts, in `assemble_split`. A plane
wave then enters through the **scattered-field** form: the incident wave already solves the
background problem, so only the difference between each shape and the background remains
as a source. The test `dielectric_cylinder_matches_the_mie_series` checks this against the
analytic solution.

## Leaky modes, with no eigensolver

PML also makes the eigenproblem non-Hermitian, and its coefficients depend on $k$. Leaky
modes could use a general eigensolver such as Krylov-Schur
([Stewart 2002](10-references.md#stewart-2002)). Driven solves turn out to be enough.

Put a line source at $\mathbf{r}_0$, so $A(k)\,u = ik\,e_0$, where $e_0$ reads the field
there. Dividing by $ik$ gives the Green's function, so

$$F(k) = 2\,\mathrm{Re}(-u_0)$$

is $\pi$ times the **local density of states** per unit $k$, how many modes the source can
radiate into. In vacuum it's $k/2$, which `vacuum_ldos_is_half_of_k` checks.

The response $u_0(k)$ is analytic above the real axis, because a causal system can't
respond before it's driven. The stretch keeps that property, since $s = 1 + i\sigma/k$ has
no singularity up there. For such a function, the value at $k + i\Gamma$ is the
Lorentzian average of its values on the real axis. So one solve, at
$k_0(1 + i/2Q)$, gives the LDOS averaged over a bandwidth, however sharp the resonances
inside it ([Liang & Johnson 2013](10-references.md#liang-johnson-2013)). On a fixed mesh this
is an exact identity, which `complex_frequency_is_a_lorentzian_average_of_real_ones` checks.

Two warnings about that average. In 2D, a Lorentzian's slow tail meets an LDOS that grows
with $k$, so the average carries a small term that grows logarithmically as cells shrink.
At $Q = 100$, vacuum reads about 1% high, and each halving of the cell adds about 0.2%. At
$Q = 20$ it's about 5%, adding 1%. And the identity's tail reaches the sharp modes far above the
mesh's cutoff, which is why the test combines two widths to get a tail that falls faster.

A leaky mode is a pole of that same response. Near it, $u_0 \approx R/(k - k_p)$, so $1/u_0$
crosses zero at the complex $k_p$, and a secant iteration finds it in 6 or 7 driven solves
(`resonance` in `solve2d.rs`). Then

$$Q = \frac{\mathrm{Re}\,k_p}{2\,|\mathrm{Im}\,k_p|}, \qquad
V = \frac{\int q\,|u|^2\,dA}{\max\,q\,|u|^2}.$$

The solve closest to the pole is dominated by the mode, so it gives the field for $V$, with
the integral stopping at the layer, whose artificial field isn't part of the mode's energy.
Because nothing varies along $z$, $V$ is an area. It only means the usual electric mode
volume for $E_z$: for $H_z$ it measures magnetic energy. Between degenerate modes it depends
on which mix the eigensolver returned.

The search needs a source where the mode is strong and a target within about a quarter of
the answer. It refuses a pole further away than that, because in an empty domain it
otherwise wanders up to a pole of the mesh and reports a plausible-looking $Q$.

### Example: one rod missing

Remove one rod from a square lattice of rods, with $\varepsilon = 12$ and $r = 0.2a$. Its
$E_z$ gap runs from $a/\lambda = 0.281$ to $0.417$, and one mode moves into it. Plane waves put
it at $0.3746$ in a periodic 5 × 5 supercell, and this solver converges to that from above.
With open walls, more crystal blocks more leakage, while the mode's shape stays put:

| Rods | $a/\lambda$ | $Q$ | $V / a^2$ |
|---:|---:|---:|---:|
| 3 × 3 | 0.3786 | 20 | 0.64 |
| 5 × 5 | 0.3784 | 198 | 0.64 |
| 7 × 7 | 0.3781 | 1,650 | 0.64 |
| 9 × 9 | 0.3781 | 12,700 | 0.65 |
| 13 × 13 | 0.3781 | 683,000 | 0.65 |

These use 12 cells per period, and finer meshes move $Q$ by about 4%. So in 2D, $Q/V$ still
changes mostly through $Q$ here, but $V$ now has a real meaning: this mode's area is smaller
than one lattice cell.

## Waveguides: power around a corner

To ask how much light gets around a bend, you need to launch one guided mode, and to measure
power. Both are short.

A **port** is a current sheet with the mode's shape. The 1D solver from chapter 7 already
finds a layered guide's modes, so `port_mode` reads the layers along the sheet and asks it.
Modes are orthogonal with weight $p$, so a sheet of $p\,\phi(y)$ excites only that mode, with

$$u = -\frac{k_0}{2\beta}\,\phi(y)\,e^{i\beta|x - x_0|},$$

going both ways. The half going backwards leaves through the absorbing layer.

A **monitor** integrates the Poynting flux across a line,

$$P = \frac{1}{2k_0}\int p\,\mathrm{Im}\!\left(u^*\frac{\partial u}{\partial n}\right)ds.$$

There's a trap here. A second-order element's derivative is least accurate at the edge of a
cell, which is exactly where a line along the grid sits: this read 2% high on a slab. Power
is conserved across a thin strip, so `monitor_power` averages over two cells instead. Around
a line source, the power crossing a closed loop then matches $F/4$, with the LDOS $F$ from
earlier, within 0.3%.

### Example: a sharp bend

[Mekis et al. 1996](10-references.md#mekis-1996) showed that a photonic crystal can steer light
around a corner of zero radius. Their crystal has rods of index 3.4 and radius $0.18a$ on a
square lattice, with a row removed to make the guide. No extended mode exists in the band gap
for light to radiate into, so what doesn't get through must come back. A line source in the
guide excites the one guided mode, and transmission is the power through a monitor, divided by
the same source's power through a straight guide:

| Check | Mekis et al. | This solver |
|---|---:|---:|
| Guided band at $k = 0$ | 0.312 | 0.312 |
| Top of the gap, reached at $k = 0.38$ | 0.443 | 0.445 |
| Sharp bend at $\omega = 0.353$ | 98% at most | 97.9% |
| Sharp bend at $\omega = 0.392$ | over 95% | 96.6% |
| Sharp bend at $\omega = 0.403$ | about 95% | 95.8% |
| Sharp bend at $\omega = 0.42$, read from their figure | about 90% | 90.8% |
| Stepped bend at $\omega = 0.353$ | 100% | 99.8% |
| Slab of index 3.5, sharp bend | 30% at most | 23% |

Frequencies are in units of $a/\lambda$. They computed these with a time-domain method, so
the agreement is between two quite different solvers. The slab row is the weakest: they don't
give that guide's width, so this uses a width of $0.143\lambda$, from their other bend. Thinner
slabs do exceed 30%, reaching 57% at $0.06\lambda$, because the mode is then mostly in the air.
A wider monitor also collects some of the radiated light, which reads about 18% to 23%.

Getting this right took a thicker layer than free space needs. Rods running through a thin
layer reflect the guided mode: 2% of its power at $2a$ of layer, which formed a resonance with
the bend and gave a transmission over 100%. At $6a$ that's gone at these frequencies, and asking
for a weaker absorption also helps. Close to the guide's cutoff, below about $\omega = 0.35$,
light travels slowly enough that reflection returns, so these numbers shouldn't be trusted there.

## When the browser isn't the right tool

For large 3D problems, compute in the browser isn't the right tool. Use
[Palace](https://github.com/awslabs/palace) (built on [MFEM](https://mfem.org), [Anderson et al. 2021](10-references.md#anderson-2021)), which
runs natively on clusters and supports eigenmode and driven simulations in 3D.

Next: [references](10-references.md).
