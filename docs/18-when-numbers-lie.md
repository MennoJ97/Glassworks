# 18. When numbers lie

*Going further. When there's no exact answer to compare with, a result can still be checked against something it must obey: an invariance, or the conservation of energy. This chapter applies one of each, and builds on chapter 8's checks and on chapters 9 and 12.*

## An invariance, and rounding at ties

From chapter 9's mesh.

Multiplying every length and the wavelength by the same number can't change anything (chapter 12), and that makes a test that needs no reference solution: `scaling_a_structure_with_its_wavelength_changes_nothing`. At first, it changed one structure's gradients by up to 17%, even though no solve was wrong.
<!-- code: scaling_a_structure_with_its_wavelength_changes_nothing "let rho" .. "1e-11" -->

The cause was rounding at exact ties. Binary floating point can't represent most decimals exactly, not even 0.1, and every operation rounds ([Goldberg 1991](20-references.md#goldberg-1991)). So a value that should be exact is often off in its last digit, and a comparison at an exact tie can go either way. Geometric code is especially exposed. A decision such as which side of a line a point lies on is made by computing a number and testing its sign, and near a tie, rounding can make that decision wrong or inconsistent ([Shewchuk 1997](20-references.md#shewchuk-1997)).

A fitted mesh makes many exact ties, and a regular structure on a fixed grid makes more. A rectangle's edge can lie exactly halfway between two rows, many nodes can be equally far from a shape, and a domain can be exactly 84 cells wide. Rounding then decides each tie, and differently at different scales, so the same structure gets a different mesh for no physical reason. The scaling test found four:

| Tie | What rounding did |
|---|---|
| A count of steps along a grid, monitor or port is a whole number | `ceil(84.00000000000001)` adds a row |
| A shape's corner is exactly halfway between two nodes | It decides which node moves |
| Several nodes are equally far from a shape | It decides the order they snap in |
| A triangle's centroid lies exactly on an outline, where it meets the outer edge whose nodes can't move | It picks that cell's material |
<!-- code: mesh.rs "Halfway between two nodes" .. "ny as f64 - 1e-9" -->
<!-- code: mesh.rs "A shape on a grid has many nodes" .. "order.sort_by" -->
<!-- code: mesh.rs "Where an outline meets" .. "rposition" -->

Geometric algorithms deal with ties in two standard ways: compute the decision exactly ([Shewchuk 1997](20-references.md#shewchuk-1997)), or break every tie by a fixed rule, so that it always falls the same way ([Edelsbrunner & Mücke 1990](20-references.md#edelsbrunner-mucke-1990)). This mesher takes the second, simpler route. It decides each of these with a small tolerance relative to the cell, so a near-tie always falls in the same direction, whatever the scale. Results at different scales now agree to $10^{-11}$.
<!-- code: steps_in -->

## Conservation, and a transmission of 465%

From chapter 12's bend.

Open that view at a wavelength of 1 µm, with the crystal's rods 1 µm apart, far outside the band gap. An early version of the page reported a transmission of 465%. No energy is being created. The number is a ratio of two separate solves, the bend's output over that of a straight guide driven the same way, a common normalization in waveguide design ([Christiansen & Sigmund 2021](20-references.md#christiansen-sigmund-2021)). That ratio is only a transmission if the straight guide delivers a steady, known power, which it does only when the light is guided.

The conserved reference is the power the source itself emits, which the LDOS gives as $F/4$ (chapter 12). No monitor can collect more than that. Against it, energy is conserved in every case:

| Wavelength | Cells per µm | Straight, share of source | Bend, share of source | The ratio |
|---:|---:|---:|---:|---:|
| 1 µm | 6 | 3% | 16% | 465% |
| 1 µm | 9 | 16% | 12% | 91% |
| 2 µm | fine enough | 5% | 15% | 336% |
| 2.833 µm | 6 | 50% | 43% | 98% |
| 2.833 µm | 9 | 50% | 47% | 98% |
| 3.5 µm | 6 | 5% | 9% | 270% |

Two different things go wrong outside the band. The first is physics. A crystal waveguide guides only inside a band gap ([Joannopoulos et al. 2008](20-references.md#joannopoulos-2008)), and this crystal's $E_z$ gaps are at $a/\lambda$ of 0.304 to 0.446 and 0.741 to 0.766. A wavelength of 1 µm, with $a = 1$ µm, lands in overlapping bands. Light leaks through the crystal, and each monitor catches whatever happens to cross its line, so the ratio means nothing. At 2 µm it goes from 336% to 418% on a finer mesh, while the bend's share of the source's power stays near 14%.

The second is the mesh, which has to resolve the shortest wavelength in the structure ([Jin 2014](20-references.md#jin-2014)). At 1 µm, a wavelength inside a rod is 0.29 µm, and the first row's cells are 0.167 µm, under 2 cells per wavelength. Nothing has converged, and the ratio jumps from 465% to 55% to 91% as the mesh is refined.

The 50% in the good rows is itself a check. By symmetry, a line source in a guide sends half its power each way, so a straight guide that really carries the light delivers about half of what's emitted. It does, from 41% to 50% for $\lambda$ from $2.3a$ to $3.2a$, collapsing just outside, to 16% at $2.26a$ and 1.4% at $3.3a$. That needs no knowledge of the band, and works for any structure. So the page tests it, and only shows a ratio when it's believable:

| Situation | The page shows |
|---|---|
| Too few cells per wavelength in the densest material, even at the mesh limit | "Not resolved" |
| The straight guide's output isn't near half of what's emitted | The bend's share of what's emitted |
| The ratio is well above 100% | "Unreliable" |
| Otherwise | The ratio |

The mesh also refines itself for a short wavelength, as far as the limit on unknowns allows, and it never gets coarser than was asked for.

Next: [designing in the other polarization](19-designing-in-hz.md).
