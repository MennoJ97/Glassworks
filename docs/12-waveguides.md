# 12. Waveguides and bends

To ask how much light gets around a corner, we need to launch one guided mode, and to measure power.
Both are short, and the answer turns out to need more care than the calculation.

## Ports

A **port** is a current sheet with a mode's shape. A straight guide is layered across its width, so
chapter 7's 1D solver already finds its modes, and `port_mode` in
[`solve2d/driven.rs`](../crates/eigenmenno-core/src/solve2d/driven.rs) simply reads the layers along the sheet and asks it.

Modes are orthogonal with weight $p$, so a sheet of $p\,\phi(y)$ excites only the mode $\phi$, with

$$u = -\frac{k_0}{2\beta}\,\phi(y)\,e^{i\beta|x - x_0|},$$
<!-- code: port "for pair in profile.windows" .. "rhs[r] += I * k" -->

going both ways. This is chapter 7's current sheet again, with a mode in place of the plane wave. The half
going backwards leaves through the absorbing layer. On an exact slab mode, the field is within 1.2% and the
power within 0.5%.
<!-- code: port_launches_the_exact_slab_mode "for x in [0.0, 1.5, 3.0]" .. "close(monitors(&o)[0]" -->

A port needs its line to cross a layered guide, which a photonic crystal's isn't. There, a line source does
the same job, since nothing else exists in the gap for it to excite.

## Monitors

A **monitor** integrates the Poynting flux across a line:

$$P = \frac{1}{2k_0}\int p\,\mathrm{Im}\!\left(u^*\frac{\partial u}{\partial n}\right)ds.$$
<!-- code: monitor_power "let (val, (gx, gy))" .. "total += weight" -->

There's a trap here. A second-order element's derivative is least accurate at the edge of a cell, which is
exactly where a line along the grid sits: this read 2% high on a slab. Power is conserved across a thin strip,
so the monitor averages over two cells instead. Around a line source, the power crossing a closed loop then
matches $F/4$, with the LDOS $F$ from chapter 11, within 0.3%. `power_around_a_line_source_matches_its_ldos` checks this.
<!-- code: monitor_power "for i in 0..steps" .. "offset * cell * ny + pml" -->

### Example: a sharp bend

[Mekis et al. 1996](20-references.md#mekis-1996) showed that a photonic crystal can steer light around a
corner of zero radius. Their crystal has rods of index 3.4 and radius $0.18a$ on a square lattice, with a row
removed to make the guide. No extended mode exists in the band gap for light to radiate into, so what doesn't
get through must come back.

The page's **2D waveguide bend** view, under **Driven response**, does this. Outside the band gap a ratio like this stops meaning anything, and chapter 18 shows how the page knows. A line source in the guide excites
the one guided mode, and transmission is the power through a monitor, divided by the same source's power
through a straight guide:
<!-- code: sharp_bends_transmit_as_in_mekis "let (straight, bend)" .. "let t =" -->

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
<!-- code: guided_band_matches_mekis "band(0.0, 0.312)" .. "band(0.38, 0.443)" -->
<!-- code: sharp_bends_transmit_as_in_mekis "close(transmission(Guide::Sharp, 0.353)" .. "close(transmission(Guide::Stepped" -->
<!-- code: a_sharp_bend_in_a_slab_loses_most_of_its_power -->

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
<!-- code: lattice_guide "let (arm, clad, pml)" .. "let rods" -->

## Scale invariance

This crystal guides light only near a wavelength of 2.8 times its rod spacing. How do we get a working bend at 1 µm, then? Look at the equation again, away from any source. Multiply every
length by a number $s$, so $\nabla \to \nabla/s$, and every wavelength too, so $k_0 \to k_0/s$:

$$-\frac{1}{s^2}\nabla\cdot(p\nabla u) - \frac{k_0^2}{s^2}\,q\,u = 0.$$

Both terms pick up the same factor. Maxwell's equations have no length of their own, so a structure scaled along
with its wavelength behaves exactly as it did. A crystal with rods 0.353 µm apart reads the same 97.8% at 1 µm that
the original reads at 2.833 µm, and at 0.55 µm and 10 µm too, to all the digits we printed. The mesh scales with it, since `steps_in` divides each length into the same number of cells at any size.
The test `scaling_a_structure_with_its_wavelength_changes_nothing` holds this to $10^{-7}$. Chapter 18 tells the story of
what it took to get there. Resizing is what the page offers when a wavelength doesn't suit a structure.

This is a real design principle: a design for one wavelength is a design for any, so long as the materials' indices
stay the same. Real materials' indices don't, so that last condition is where the principle stops.

Next: [inverse design](13-inverse-design.md).
