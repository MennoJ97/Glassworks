# 18. When numbers lie

*Going further. Two stories of a result that looked wrong, or looked fine and wasn't, and how each was tracked down. They build on chapter 8's habit of checking, and on chapters 9 and 12.*

## Ties, and what rounding noise can do

From chapter 9's mesh.

A fitted mesh makes many exact ties, and a fixed grid under a regular structure makes more. A
rectangle's edge can lie exactly halfway between two rows, many nodes can be equally far from a
shape, and a domain can be exactly 84 cells wide. What decides each tie is then floating-point
noise. This looks harmless until the same structure gets a different mesh for no physical reason.

This solver has a test that made that visible,
`scaling_a_structure_with_its_wavelength_changes_nothing`. Multiplying every length and the wavelength
by the same number shouldn't change anything (chapter 12 explains why). At first, it changed one structure's gradients by
up to 17%, even though no solve was wrong. It took four separate fixes to find all of these:
<!-- code: scaling_a_structure_with_its_wavelength_changes_nothing "let rho" .. "1e-11" -->

| Tie | What noise did |
|---|---|
| A count of steps along a grid, monitor or port is a whole number | `ceil(84.00000000000001)` adds a row |
| A shape's corner is exactly halfway between two nodes | It decides which node moves |
| Several nodes are equally far from a shape | It decides the order they snap in |
| A triangle's centroid lies exactly on an outline, where it meets the outer edge whose nodes can't move | It picks that cell's material |
<!-- code: mesh.rs "Halfway between two nodes" .. "ny as f64 - 1e-9" -->
<!-- code: mesh.rs "A shape on a grid has many nodes" .. "order.sort_by" -->
<!-- code: mesh.rs "Where an outline meets" .. "rposition" -->

The fix is the same each time: decide with a small tolerance, relative to the cell, so that a tie
always falls the same way. Results at different scales now agree to $10^{-11}$. The last one took
longest to find, because the field was only wrong in the last column of cells, next to the wall, and
its effect was $10^{-6}$. Looking at *where* two solutions differ located it, after several good
theories about the eigensolver and the linear solver had been ruled out by measuring them.
<!-- code: steps_in -->

## When a transmission reads 465%

From chapter 12's bend. Open that view at a wavelength of 1 µm with the crystal's rods 1 µm apart, far outside the band gap, and an
early version of the page reported a transmission of 465%. No energy is being created. The number is a ratio of
two separate solves, and it can't promise to stay below 100%, because its reference can fail.

The two monitors read 0.156 and 0.0336. Dividing by the straight guide assumes that it delivers a steady, known
power, and here it doesn't. There's a better reference: the power the source itself emits, which the LDOS gives,
as $F/4$. Against that, energy is conserved in every case:

| Wavelength | Cells per µm | Straight, share of source | Bend, share of source | The ratio |
|---:|---:|---:|---:|---:|
| 1 µm | 6 | 3% | 16% | 465% |
| 1 µm | 9 | 16% | 12% | 91% |
| 2 µm | fine enough | 5% | 15% | 336% |
| 2.833 µm | 6 | 50% | 43% | 98% |
| 2.833 µm | 9 | 50% | 47% | 98% |
| 3.5 µm | 6 | 5% | 9% | 270% |

Two different things go wrong outside the band. The first is physics: light isn't guided there. It leaks
through the crystal, and each monitor catches whatever happens to cross its line, so the ratio means nothing.
Its value at 2 µm goes from 336% to 418% on a finer mesh, while the bend's share of the source's power stays near
14%. The share is meaningful, and the ratio isn't.

The second is the mesh. At 1 µm, a wavelength inside a rod is 0.29 µm, and the first row's cells are 0.167 µm, so
under 2 cells per wavelength. Nothing has converged, and the ratio there jumps from 465% to 55% to 91%.

The 50% in the good rows is itself a check. A line source in a guide sends half its power each way, so a straight
guide that really carries the light delivers about half of what's emitted, and it does: 41% to 50% from
$\lambda = 2.3a$ to $3.2a$, collapsing just outside, to 16% at $2.26a$ and 1.4% at $3.3a$. That needs no knowledge of the band, and works for any
structure. So the page tests it, and only shows a ratio when it's believable:

| Situation | The page shows |
|---|---|
| Too few cells per wavelength in the densest material, even at the mesh limit | "Not resolved" |
| The straight guide's output isn't near half of what's emitted | The bend's share of what's emitted |
| The ratio is well above 100% | "Unreliable" |
| Otherwise | The ratio |

The mesh also adjusts itself to a short wavelength, as far as the limit on unknowns allows, and never gets
coarser than was asked for. With up to 2.3 million unknowns, solved outside the browser, the crystal at 1 µm
converges to about 19% through the straight guide and 11% through the bend. The page's own mesh reads 11.6%.

A band diagram explains why: this crystal's $E_z$ gaps are at $a/\lambda$ of 0.304 to 0.446 and 0.741 to 0.766,
and a wavelength of 1 µm, with $a = 1$ µm, lands in overlapping bands. The second gap does guide light, reaching 36%
to 48% at 1.31 to 1.35 µm, but the sharp bend passes only about 4 to 9% of it at 1.33 and 1.35 µm.

Next: [designing in the other polarization](19-designing-in-hz.md).
