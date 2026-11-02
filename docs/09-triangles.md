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
between communities.

This is chapter 1's equation with $d/dx$ replaced by $\nabla$. Integrating by parts in the same
way gives the same weak form:

$$\int p\,\nabla u\cdot\nabla v\,dA - k_0^2\int q\,u\,v\,dA = \int f\,v\,dA.$$

So [`assemble.rs`](../crates/eigenmenno-core/src/assemble.rs) doesn't change. It only ever saw a
`Space` that gives it local matrices and dofs. Two things are new: a `Space` on triangles, and a
mesh of them. Chapter 14 explains why this doesn't extend to in-plane vector fields, where the
scalar picture breaks down.

## Second-order triangles

Chapter 3 built a basis on a line from polynomials through evenly spaced nodes. The natural
coordinates on a triangle are **barycentric**: $\lambda_0, \lambda_1, \lambda_2$ are each 1 at one
vertex and 0 on the opposite edge, and they add up to 1. A second-order element has six nodes,
at the three vertices and three edge midpoints, with basis functions

$$\lambda_i(2\lambda_i - 1) \quad\text{at a vertex}, \qquad 4\lambda_i\lambda_j \quad\text{at the midpoint between } i \text{ and } j.$$

Each is 1 at its own node and 0 at the other five. Along an edge, each of these is a quadratic
fixed by that edge's three nodes, and two triangles sharing an edge share those nodes, so $u$ is
continuous across it with nothing else to enforce.

On a straight triangle, the gradients of the $\lambda_i$ are constant, so the stiffness integrand
is a polynomial of degree 2, and the mass integrand has degree 4. Chapter 3's quadrature idea
applies directly, with a small rule on the triangle.
[`space.rs`](../crates/eigenmenno-core/src/space.rs) has all of this in `H1Triangle`.

On a mesh that fits the geometry, eigenvalue error should fall as $h^{2p}$, as chapter 8 showed.
With $p = 2$ that's a factor of 16 per halving, and a square cavity's error falls about 15 times.

## A mesh that fits the shapes

A general mesher, such as [Triangle](https://www.cs.cmu.edu/~quake/triangle.html)
([Shewchuk 1996](15-references.md#shewchuk-1996)), could handle any geometry. Here a simpler idea
goes a long way, and it has a useful property that a general mesher lacks.

Start from a structured grid. On a **staircase** mesh, each cell takes the material at its centre,
so a circle becomes a jagged outline. Errors there stop falling with the mesh, because the shape's
outline is only ever approximated to first order, with a jump in the wrong place:

| Grid (unknowns) | $E_z$ staircase | $E_z$ fitted | $H_z$ staircase | $H_z$ fitted |
|---|---:|---:|---:|---:|
| 12 × 12 (576) | 1.3e-2 | 3.6e-3 | 3.9e-2 | 7.0e-3 |
| 24 × 24 (2,304) | 6.7e-4 | 9.0e-4 | 9.3e-3 | 1.6e-3 |
| 48 × 48 (9,216) | 5.9e-4 | 9.6e-5 | 9.7e-3 | 1.2e-4 |

These are errors in a rod lattice's band edges, explained below. `Mesh::fitted_grid` in
[`mesh.rs`](../crates/eigenmenno-core/src/mesh.rs) fits the outline instead, in three steps:

1. **Snap nodes onto outlines.** Corners of a shape take the nearest grid node. Then, wherever a
   grid line crosses an outline, one of its two ends moves onto the outline. Nodes very close to an
   outline join it too, or they'd form slivers against the ones that moved.
2. **Choose diagonals afterwards.** Each grid square is cut into two triangles. With no grid line
   crossing an outline any more, at least one diagonal doesn't cross it either. Doing this after
   snapping matters: an early version snapped along diagonals as well, which pulled two neighbours
   across a curve onto almost the same point.
3. **Flip to Delaunay.** Snapping still leaves some thin triangles along curves. Two triangles
   sharing an edge form a quadrilateral with two possible diagonals, and the better one is where the
   angles opposite the edge add up to less than 180°. Repeating this removes slivers, and a flip that
   would cross an outline is refused, so every triangle stays in one material.

A triangle's material is then read at its centroid. Over 200 random placements, the worst angle was
22° for circles and 16° for rectangles.

Nodes on the domain's outer edge never move. That's the useful property: opposite edges of the
domain have matching nodes, which is exactly what periodic boundaries need.

Straight cells still cut the corners of a circle, so a curved outline's error falls at second order,
even with second-order elements. Fitting doesn't undo that: it just gets the outline's position right.

## Ties, and what rounding noise can do

A fitted mesh makes many exact ties, and a fixed grid under a regular structure makes more. A
rectangle's edge can lie exactly halfway between two rows, many nodes can be equally far from a
shape, and a domain can be exactly 84 cells wide. What decides each tie is then floating-point
noise. This looks harmless until the same structure gets a different mesh for no physical reason.

This solver has a test that made that visible: multiplying every length and the wavelength by the same
number shouldn't change anything (chapter 12 explains why). At first, it changed one structure's gradients by
up to 17%, even though no solve was wrong. It took four separate fixes to find all of these:

| Tie | What noise did |
|---|---|
| A count of steps along a grid, monitor or port is a whole number | `ceil(84.00000000000001)` adds a row |
| A shape's corner is exactly halfway between two nodes | It decides which node moves |
| Several nodes are equally far from a shape | It decides the order they snap in |
| A triangle's centroid lies exactly on an outline, where it meets the outer edge whose nodes can't move | It picks that cell's material |

The fix is the same each time: decide with a small tolerance, relative to the cell, so that a tie
always falls the same way. Results at different scales now agree to $10^{-11}$. The last one took
longest to find, because the field was only wrong in the last column of cells, next to the wall, and
its effect was $10^{-6}$. Looking at *where* two solutions differ located it, after several good
theories about the eigensolver and the linear solver had been ruled out by measuring them.

## Periodic ties and band diagrams

Chapter 4 handled a Bloch boundary by mapping the last dof onto the first, with a phase. In 2D this is
the same map with many entries. A periodic cell has two directions, with a phase across each set by
the wavevector $\mathbf{k}$:

$$u(x + L_x, y) = e^{ik_x L_x}\,u(x, y), \qquad u(x, y + L_y) = e^{ik_y L_y}\,u(x, y).$$

Every dof on the right edge ties to its partner on the left, and every one on the top to its partner
on the bottom. Corners tie to the origin, carrying both phases. Because outer nodes never moved, each
partner exists, and `DofMap` needs nothing else: assembly already multiplies each entry by its dofs'
phases, as chapter 4 showed. The matrix is Hermitian, so Lanczos and shift-invert from chapter 6 work
unchanged.

A **band diagram** repeats that eigenproblem for wavevectors around the edge of the irreducible
Brillouin zone: from $\Gamma$ ($\mathbf{k} = 0$) to X to M and back to $\Gamma$. For a square lattice, the
extremes of the bands are usually found on this path, though this isn't guaranteed
([Harrison et al. 2007](15-references.md#harrison-2007)). Only the phases change between
wavevectors, so where the nonzeros are stays the same, which the next chapter exploits.

### Example: a lattice of rods

Take a square lattice of rods with $\varepsilon = 12$ and radius $r = 0.2a$, on the **2D band diagram**
view. This is checked against an independent plane-wave expansion. Frequencies are in units of $a/\lambda$,
and M1 and X2 mean the first band at M and the second at X:

| Grid (unknowns) | $E_z$ M1 | $E_z$ X2 | $H_z$ M1 | $H_z$ X2 |
|---|---:|---:|---:|---:|
| 12 × 12 (576) | 0.28456 | 0.42346 | 0.50783 | 0.44746 |
| 24 × 24 (2,304) | 0.28138 | 0.41835 | 0.49786 | 0.44235 |
| 48 × 48 (9,216) | 0.28084 | 0.41747 | 0.49637 | 0.44155 |
| 96 × 96 (36,864) | 0.28071 | 0.41724 | 0.49600 | 0.44135 |
| Plane waves | 0.28066 | 0.41716 | 0.49588 | 0.44109, still rising |

The solver comes down from above, as chapter 2 predicts. For $E_z$, this lattice has a **band gap**
from 0.281 to 0.417: no wave at those frequencies can travel through it in any direction. That gap
is what makes the guides and cavities of the following chapters work.

Next: [solving on a 2D mesh](10-sparse-lu.md).
