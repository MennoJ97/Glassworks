# 9. Fields on triangles

Everything so far has been 1D. Most of it carries over to 2D, and this chapter and the next
four follow what changes, from the mesh up to designing structures. They assume chapters 1
to 8. To try any of this, each 2D problem has a view on the page, such as **2D band diagram**,
under its problem at the top.

## One scalar field

The easiest 2D problem still has a single unknown. A structure that's uniform along $z$, with
$u = E_z$ or $u = \eta_0 H_z$, gives

$$-\nabla\cdot(p\nabla u) - k_0^2 q\,u = f,$$

with $p = 1/\mu$, $q = \varepsilon$ for $E_z$, and $p = 1/\varepsilon$, $q = \mu$ for $H_z$. The two
are chapter 1's two polarizations, and they behave quite differently, as bands below show. The
page lists them by which field points out of the plane, because the names TE and TM swap
between communities. In the code, `Te` solves for $E_z$ and `Tm` for $H_z$.

This is chapter 1's equation with $d/dx$ replaced by $\nabla$. Integrating by parts in the same
way gives the same weak form:

$$\int p\,\nabla u\cdot\nabla v\,dA - k_0^2\int q\,u\,v\,dA = \int f\,v\,dA.$$

So [`assemble.rs`](../crates/glassworks-core/src/assemble.rs) doesn't change. It only ever saw a
`Space` that gives it local matrices and dofs. Two things are new: a `Space` on triangles, and a
mesh of them. Chapter 14 explains why this doesn't extend to in-plane vector fields, where the
scalar picture breaks down.

## Second-order triangles

Chapter 3 built a basis on a line from polynomials through evenly spaced nodes. The natural
coordinates on a triangle are **barycentric**: $\lambda_0, \lambda_1, \lambda_2$ are each 1 at one
vertex and 0 on the opposite edge, and they add up to 1. A second-order element has six nodes,
at the three vertices and three edge midpoints, with basis functions

$$\lambda_i(2\lambda_i - 1) \quad\text{at a vertex}, \qquad 4\lambda_i\lambda_j \quad\text{at the midpoint between } i \text{ and } j.$$
<!-- code: basis_values -->

Each is 1 at its own node and 0 at the other five. Along an edge, each of these is a quadratic
fixed by that edge's three nodes, and two triangles sharing an edge share those nodes, so $u$ is
continuous across it with nothing else to enforce.

On a straight triangle, the gradients of the $\lambda_i$ are constant, so the stiffness integrand
is a polynomial of degree 2, and the mass integrand has degree 4. Chapter 3's quadrature idea
applies directly, with a small rule on the triangle ([Dunavant 1985](20-references.md#dunavant-1985)).
[`space.rs`](../crates/glassworks-core/src/space.rs) has all of this in `H1Triangle`.
<!-- code: TRIANGLE_RULE -->

Barycentric coordinates do the job of chapter 3's reference element. A point of the triangle with
corners $\mathbf{p}_0, \mathbf{p}_1, \mathbf{p}_2$ is
$\mathbf{x} = \mathbf{p}_0 + \lambda_1(\mathbf{p}_1 - \mathbf{p}_0) + \lambda_2(\mathbf{p}_2 - \mathbf{p}_0)$,
and `det` is the Jacobian determinant of that map, twice the area. Its sign only records whether the
corners run anticlockwise, and dividing by it gives the right gradients either way. Each $\lambda_i$
is linear, 0 on the opposite edge and 1 at vertex $i$, so $\nabla\lambda_i$ points from that edge
straight towards vertex $i$, with length one over the triangle's height there. On the triangle
(0, 0), (1, 0), (0, 1), `det` is 1 and the gradients are (−1, −1), (1, 0) and (0, 1). The rule's
weights add up to 1, so multiplying by the area turns the rule into an integral over the cell.
<!-- code: H1Triangle::shape -->
<!-- code: H1Triangle::local_matrices "let (det" .. "mass[i * 6 + j]" -->

Two checks need no integrals. The six basis functions add up to 1 everywhere, and a constant has no
gradient, so each row of the stiffness matrix sums to 0. For the same reason, all 36 entries of the
mass matrix add up to the triangle's area.

On a mesh that fits the geometry, eigenvalue error should fall as $h^{2r}$ for order $r$, as chapter 8
showed. Second-order elements give $h^4$, a factor of 16 per halving, and a square cavity's error falls
about 15 times.

### Example: an empty box

Take an empty 1 × 1.5 box with perfectly conducting walls, so $E_z = 0$ on them. Its modes are products
of 1D cavity modes like chapter 3's, $\sin(m\pi x)\sin(n\pi y/1.5)$, with $k_0^2 = \pi^2(m^2 + n^2/2.25)$ for
$m, n \ge 1$. The lowest is $k_0 = \pi\sqrt{1 + 1/2.25} = 3.776$. For $H_z$ the same walls make its
normal derivative 0 instead, so its modes are cosines, $m$ or $n$ can be 0, and a constant mode comes
first, at $k_0 = 0$. The test `rectangular_cavity_matches_mode_formula` checks the first six of each.
<!-- code: rectangular_cavity_matches_mode_formula -->

A box is also the bridge from 1D. A short one layered along $x$, with magnetic walls, has first modes
with no $y$ variation at all. Those are exactly the 1D solver's modes, and
`layered_box_matches_the_1d_solver` checks that the 2D solver finds them.
<!-- code: layered_box_matches_the_1d_solver -->

## A mesh that fits the shapes

A general mesher, such as [Triangle](https://www.cs.cmu.edu/~quake/triangle.html)
([Shewchuk 1996](20-references.md#shewchuk-1996)), could handle any geometry. Here a simpler idea
goes a long way, and it has a useful property that a general mesher lacks.

Start from a structured grid. On a **staircase** mesh, each cell takes the material at its centre,
so a circle becomes a jagged outline. Errors there stop falling with the mesh, because the shape's
outline is only ever approximated to first order, with a jump in the wrong place:

| Grid (unknowns) | $E_z$ staircase | $E_z$ fitted | $H_z$ staircase | $H_z$ fitted |
|---|---:|---:|---:|---:|
| 12 × 12 (576) | 1.3e-2 | 3.6e-3 | 3.9e-2 | 7.0e-3 |
| 24 × 24 (2,304) | 6.7e-4 | 9.0e-4 | 9.3e-3 | 1.6e-3 |
| 48 × 48 (9,216) | 5.9e-4 | 9.6e-5 | 9.7e-3 | 1.2e-4 |

These are errors in the band edges of the rod lattice at the end of this chapter. A staircase can land
close by luck at one mesh, as $E_z$ does at 24 × 24, but it stops improving. `Mesh::fitted_grid` in
[`mesh.rs`](../crates/glassworks-core/src/mesh.rs) fits the outline instead, in three steps:

1. **Snap nodes onto outlines** (`Mesh::snap`). Corners of a shape take the nearest grid node. Then,
   wherever a grid line crosses an outline, one of its two ends moves onto the outline. Nodes very
   close to an outline join it too, or they'd form slivers against the ones that moved.
2. **Choose diagonals afterwards.** Each grid square is cut into two triangles. With no grid line
   crossing an outline any more, at least one diagonal doesn't cross it either. Doing this after
   snapping matters: an early version snapped along diagonals as well, which pulled two neighbours
   across a curve onto almost the same point.
3. **Flip to Delaunay** (`Mesh::flip_to_delaunay`). Snapping still leaves some thin triangles along
   curves. Two triangles sharing an edge form a quadrilateral with two possible diagonals, and the
   better one is where the angles opposite the edge add up to less than 180° ([Shewchuk 1996](20-references.md#shewchuk-1996)). Repeating this removes
   slivers, and a flip that would cross an outline is refused, so every triangle stays in one material.
<!-- code: Mesh::fitted_grid "m.snap" .. "flip_to_delaunay" -->

A triangle's material is then read at its centroid. Over 200 random placements, the worst angle was
22° for circles and 16° for rectangles.
<!-- code: Mesh::fitted_grid "m.region =" .. ".collect()" -->

Nodes on the domain's outer edge never move. That's the useful property: opposite edges of the
domain have matching nodes, which is exactly what periodic boundaries need.
<!-- code: Mesh::snap "let free" .. "};" -->

Straight cells still cut the corners of a circle, so a curved outline's error falls at second order,
even with second-order elements. Fitting doesn't undo that: it just gets the outline's position right.

## Periodic ties and band diagrams

Chapter 4 handled a Bloch boundary by mapping the last dof onto the first, with a phase. In 2D this is
the same map with many entries. A periodic cell has two directions, with a phase across each set by
the wavevector $\mathbf{k}$:

$$u(x + L_x, y) = e^{ik_x L_x}\,u(x, y), \qquad u(x, y + L_y) = e^{ik_y L_y}\,u(x, y).$$
<!-- code: constraints "let angle" .. "from_polar" -->

Every dof on the right edge ties to its partner on the left, and every one on the top to its partner
on the bottom. Corners tie to the origin, carrying both phases. Because outer nodes never moved, each
partner exists, and `DofMap` needs nothing else: assembly already multiplies each entry by its dofs'
phases, as chapter 4 showed. The matrix is Hermitian, so Lanczos and shift-invert from chapter 6 work
unchanged.
<!-- code: constraints "let key" .. "&ties)" -->

A **band diagram** repeats that eigenproblem for many wavevectors. A phase of $2\pi$ changes nothing,
so each phase $k_xL_x$, $k_yL_y$ only needs to run from $-\pi$ to $\pi$. A square cell's symmetry, such as
a centred rod's, cuts that down to a triangle: $\Gamma$ = (0, 0), X = ($\pi$, 0), M = ($\pi$, $\pi$), as
phases per period. This triangle is called the irreducible Brillouin zone. The code walks its edge from
$\Gamma$ to X to M and back to $\Gamma$, with `s` from 0 to 3. For a square lattice, the extremes of the
bands are usually found on this path, though this isn't guaranteed
([Harrison et al. 2007](20-references.md#harrison-2007)). Only the phases change between
wavevectors, so where the nonzeros are stays the same, which the next chapter exploits.
<!-- code: closed.rs "let points" .. "band.push" -->

### Example: a lattice of rods

Take a square lattice of rods with $\varepsilon = 12$ and radius $r = 0.2a$, on the **2D band diagram**
view. This is checked against a plane-wave expansion ([Ho et al. 1990](20-references.md#ho-1990)), written
for this check, which shares nothing with this solver. Frequencies are in units of $a/\lambda$,
and M1 and X2 mean the first band at M and the second at X:

| Grid (unknowns) | $E_z$ M1 | $E_z$ X2 | $H_z$ M1 | $H_z$ X2 |
|---|---:|---:|---:|---:|
| 12 × 12 (576) | 0.28456 | 0.42346 | 0.50783 | 0.44746 |
| 24 × 24 (2,304) | 0.28138 | 0.41835 | 0.49786 | 0.44235 |
| 48 × 48 (9,216) | 0.28084 | 0.41747 | 0.49637 | 0.44155 |
| 96 × 96 (36,864) | 0.28071 | 0.41724 | 0.49600 | 0.44135 |
| Plane waves | 0.28066 | 0.41716 | 0.49588 | 0.44109, still rising |

The [MPB tutorial](https://mpb.readthedocs.io/en/latest/Python_Tutorial/) gives 0.28262 and 0.41933 for this
lattice's $E_z$ gap, at a resolution of 32 ([Johnson & Joannopoulos 2001](20-references.md#johnson-joannopoulos-2001)). That's 0.5 to 0.7% above
both of these. Its coarse grid may explain it, but we haven't run MPB to check.

The solver comes down from above, as chapter 2 predicts. For $E_z$, this lattice has a **band gap**
from 0.281 to 0.417: no $E_z$ wave at those frequencies can travel through it, in any direction in the plane. That gap
is what makes the guides and cavities of the following chapters work.

Next: [solving on a 2D mesh](10-sparse-lu.md).
