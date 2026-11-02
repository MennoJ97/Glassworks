# 12. Waveguides and bends

To ask how much light gets around a corner, we need to launch one guided mode, and to measure power.
Both are short, and the answer turns out to need more care than the calculation.

## Ports

A **port** is a current sheet with a mode's shape. A straight guide is layered across its width, so
chapter 7's 1D solver already finds its modes, and `port_mode` in
[`solve2d/driven.rs`](../crates/eigenmenno-core/src/solve2d/driven.rs) simply reads the layers along the sheet and asks it.

Modes are orthogonal with weight $p$, so a sheet of $p\,\phi(y)$ excites only the mode $\phi$, with

$$u = -\frac{k_0}{2\beta}\,\phi(y)\,e^{i\beta|x - x_0|},$$

going both ways. This is chapter 7's current sheet again, with a mode in place of the plane wave. The half
going backwards leaves through the absorbing layer. On an exact slab mode, the field is within 1.2% and the
power within 0.5%.

A port needs its line to cross a layered guide, which a photonic crystal's isn't. There, a line source does
the same job, since nothing else exists in the gap for it to excite.

## Monitors

A **monitor** integrates the Poynting flux across a line:

$$P = \frac{1}{2k_0}\int p\,\mathrm{Im}\!\left(u^*\frac{\partial u}{\partial n}\right)ds.$$

There's a trap here. A second-order element's derivative is least accurate at the edge of a cell, which is
exactly where a line along the grid sits: this read 2% high on a slab. Power is conserved across a thin strip,
so the monitor averages over two cells instead. Around a line source, the power crossing a closed loop then
matches $F/4$, with the LDOS $F$ from chapter 11, within 0.3%.

### Example: a sharp bend

[Mekis et al. 1996](15-references.md#mekis-1996) showed that a photonic crystal can steer light around a
corner of zero radius. Their crystal has rods of index 3.4 and radius $0.18a$ on a square lattice, with a row
removed to make the guide. No extended mode exists in the band gap for light to radiate into, so what doesn't
get through must come back.

The page's **2D waveguide bend** view, under **Driven response**, does this. A line source in the guide excites
the one guided mode, and transmission is the power through a monitor, divided by the same source's power
through a straight guide:

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

Frequencies are in units of $a/\lambda$. They computed these with a time-domain method, so the agreement is
between two quite different solvers. The slab row is the weakest: they don't give that guide's width, so this
uses a width of $0.143\lambda$, from their other bend. Thinner slabs do exceed 30%, reaching 57% at $0.06\lambda$,
because the mode is then mostly in the air. A wider monitor also collects some of the radiated light, which
reads about 18% to 23%.

Getting this right took a thicker layer than free space needs. The rods run through the layer, as they must,
and they reflect the guided mode a little: 1.6% of its power at $2a$ of layer, 0.06% at $4a$, and 0.001% at $6a$.
That small reflection formed a resonance with the bend, and gave a transmission over 100%. Asking for weaker
absorption also helps, because a steep layer reflects more. Close to the guide's cutoff, below about
$\omega = 0.35$, light travels slowly enough that reflection returns, so these numbers shouldn't be trusted there.

## When a transmission reads 465%

Open that view at a wavelength of 1 µm with the crystal's rods 1 µm apart, far outside the band gap, and an
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

A band diagram explains why: this crystal's complete gaps are at $a/\lambda$ of 0.304 to 0.446 and 0.741 to 0.766,
and a wavelength of 1 µm, with $a = 1$ µm, lands in overlapping bands. The second gap does guide light, reaching 36%
to 48% at 1.31 to 1.35 µm, but the sharp bend passes only about 4 to 9% of it at 1.33 and 1.35 µm.

## Scale invariance

How do we get a working bend at 1 µm, then? Look at the equation again, away from any source. Multiply every
length by a number $s$, so $\nabla \to \nabla/s$, and every wavelength too, so $k_0 \to k_0/s$:

$$-\frac{1}{s^2}\nabla\cdot(p\nabla u) - \frac{k_0^2}{s^2}\,q\,u = 0.$$

Both terms pick up the same factor. Maxwell's equations have no length of their own, so a structure scaled along
with its wavelength behaves exactly as it did. A crystal with rods 0.353 µm apart reads the same 97.8% at 1 µm that
the original reads at 2.833 µm, and at 0.55 µm and 10 µm too, to all the digits we printed. The mesh scales with it.
The test `scaling_a_structure_with_its_wavelength_changes_nothing` holds this to $10^{-7}$. Chapter 9 told the story of
what it took to get there. Resizing is what the page offers when a wavelength doesn't suit a structure.

This is a real design principle: a design for one wavelength is a design for any, so long as the materials' indices
stay the same. Real materials' indices don't, so that last condition is where the principle stops.

Next: [inverse design](13-inverse-design.md).
