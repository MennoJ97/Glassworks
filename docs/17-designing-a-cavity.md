# 17. Designing a cavity

*Going further. This builds on the LDOS from chapter 11 and the optimizer from chapter 13. It follows [Liang & Johnson 2013](20-references.md#liang-johnson-2013) for the goal, and [Shaker et al. 2026](20-references.md#shaker-2026) for why a sharp resonance is hard to optimize.*

Chapter 11 found the modes of a crystal someone drew. What if the solver draws the cavity? The code is in `web/js/cavity-design.js`, on the **2D cavity** view under **Inverse design**.

## The goal: LDOS at a complex frequency

A cavity's job is to make an emitter at its centre radiate more, so the goal is the LDOS at a source there. [Liang & Johnson 2013](20-references.md#liang-johnson-2013) made this a practical goal for topology optimization, by averaging the LDOS over a bandwidth, which chapter 11's complex frequency $k_0(1 + i/2Q)$ does in one solve. That rewards any resonance near the target, however sharp, rather than one that has to land exactly on it. A broad bandwidth gives a smooth goal that the optimizer can climb, and narrowing it as the optimization goes on leads towards a sharp resonance ([Shaker et al. 2026](20-references.md#shaker-2026) review these choices). Here $Q$ rises in stages. In 2D, the average reads a little high in vacuum, as chapter 11 noted: about 5% at $Q = 20$ and 1% at $Q = 100$.

## The LDOS adjoint needs no solve at all

With $A u = ik\,e_0$ and $F = -2\,\mathrm{Re}(e_0^T u)$, chapter 13's steps give

$$\delta F = 2\,\mathrm{Re}\big(e_0^T A^{-1}\delta A\,u\big) = 2\,\mathrm{Re}\big((A^{-1}e_0)^T \delta A\,u\big).$$

But $A^{-1}e_0$ is just $u/(ik)$, which we already have. Because $A$ is complex symmetric, the Green's function is reciprocal, and a source's adjoint is its own field, so the LDOS gradient needs no extra solve ([Shaker et al. 2026](20-references.md#shaker-2026), reviewing Liang & Johnson). So

$$\delta F = 2\,\mathrm{Re}\!\left(\frac{u^T \delta A\,u}{ik}\right),$$
<!-- code: density_gradient "With no monitor" .. "u.iter().map" -->

and the gradient costs one more sum over cells. `ldos_gradient_matches_finite_differences` checks this.
<!-- code: pixel_slopes "for cell in 0..space.n_cells()" .. "imaginary[pixel] +=" -->

## The setup

The benchmark is chapter 11's crystal of 7 × 7 rods with $\varepsilon = 12$, in a region of the same size, 2.67 wavelengths wide, in 64 × 64 pixels. A crystal is a good cavity, so this is a fair comparison. The design is kept symmetric under the square's eight rotations and mirrors, around the source.

The filter's radius sets the smallest feature a design can have (chapter 13), so it has to be small enough for the benchmark. A filter of 3 pixels, 0.125 wavelengths, is nearly as wide as these rods, and applied to the crystal itself it leaves a fill of 5.4% rather than the true 12.3%. So the page uses 2. A run has 8 stages of 20 iterations, ending at $Q = 500$ and $\beta = 128$.
<!-- code: CAVITY -->

## Why a sharp resonance is hard

With only the goal above, finished designs were strong, but their peaks sat 2 to 5% away from the target. At $Q = 1000$ the bandwidth is 0.1%, so they miss it completely. [Shaker et al. 2026](20-references.md#shaker-2026) explain why. As $Q$ grows, the LDOS becomes far more sensitive to a shift of the resonance than to anything else, so the goal becomes stiff along the direction in which the frequency changes. A step that's good for the cavity also moves its peak further than its own width.

Their fix is to read the LDOS at the real part of the mode's own frequency rather than at the target, finding that frequency with a few steps of a shift-invert eigensolver at each iteration. That doubles the number of factorizations. A loose bound keeps the mode from drifting far from the target. This page uses the same idea more cheaply, and adds a standard way of holding the frequency on target.

### Finding the mode and its slope

Near a pole, $1/u_0$ is nearly linear in $k$ (chapter 11), so two solves place the pole. The same complex quantity $u^T \delta A\,u / (ik)$ is the slope of $-u_0$, not only of its real part. Writing $1/u_0$ as a straight line, with the slope $d(1/u_0)/dk$ from those two solves, a change of $u_0$ moves the pole by
<!-- code: findResonance "const offset" .. "const points" -->

$$\delta k_p = \frac{\delta(1/u_0)}{-\,d(1/u_0)/dk}.$$
<!-- code: findResonance "const r = dipole" .. "moves: Float64Array" -->

So for every density $x_i$, we get a number $c_i = \partial\,\mathrm{Re}\,k_p/\partial x_i$. Measured on a true resonance, that predicts the real shift within 2%.

### Holding the frequency

Keeping $\mathrm{Re}\,k_p$ at the target is an equality constraint, and $c$ is its gradient. The classic way to improve a goal without breaking such a constraint is the **gradient projection** method ([Rosen 1960](20-references.md#rosen-1960); [Nocedal & Wright 2006](20-references.md#nocedal-wright-2006)). Every direction perpendicular to $c$ leaves the frequency alone, to first order, so remove that component from the goal's gradient $g$:

$$g \;\to\; g - \frac{g\cdot c}{c\cdot c}\,c.$$
<!-- code: optimizeCavity "const square" .. "goal.map((v, k) => v - along" -->

A step along that direction can't remove an error the frequency already has, and MMA's bounds bend the step away from it anyway. So a second, **range-space** step moves along $c$ itself, to remove a set fraction of the remaining error. Here that fraction is half. [Feppon et al. 2020](20-references.md#feppon-2020) combine the same two parts, a projected step and a correction along the constraint's gradient that shrinks the violation geometrically.
<!-- code: optimizeCavity "const moved" .. "x = next.map" -->

Following Shaker et al., each step also reads the LDOS at the mode's current wavelength rather than at the target. A Lorentzian shows why this matters. With detuning $D$, bandwidth $\Gamma$ and the mode's own decay rate $\gamma$,

$$F \propto \frac{\Gamma + \gamma}{D^2 + (\Gamma + \gamma)^2},$$

and its slope along $\gamma$ changes sign at $|D| = \Gamma + \gamma$. Closer than that, raising $F$ raises $Q$. Further out, on the peak's wing, raising $F$ lowers $Q$.
<!-- code: optimizeCavity "solve, own)" .. "own = (2 * Math.PI) / poleK" -->

The last stages of sharpening still shift the mode a little, because a pixel flip is a big step. Scale invariance (chapter 12) takes out what's left exactly: resize the finished design by $\lambda_\text{target}/\lambda_\text{mode}$.
<!-- code: optimizeCavity "tuned.set(rho" -->

One run proves little, because different starts differ by far more than most of these choices. Here are nine starts with and without the hold. Both versions read the LDOS at the target, so that only the hold differs, and both are compared at a bandwidth narrow enough for a strong mode to stand out:

| Version | Median LDOS at $Q = 10{,}000$ | Median mode $Q$ | Largest resizing |
|---|---:|---:|---:|
| Holding the frequency | 8,411 | 56,000 | 1.2% |
| Not holding it | 3,304 | 5,900 | 4.4% |

## The result

Everything is compared at $Q = 1000$, relative to vacuum, each at its own resonance. The crystal is tuned to its own, and every pattern goes through the same search for its strongest peak within 4%, then resizing:
<!-- code: optimizeCavity "const scan" .. "let peak" -->

| Structure | LDOS | $Q$ | $V$ (wavelengths²) |
|---|---:|---:|---:|
| 16 random symmetric patterns | 0 to 29 | | |
| Crystal of 7 × 7 rods | 299 | 1,680 | 0.093 |
| 18 optimized starts | 886 to 1,075, median 1,014 | 9,120 to 5,859,698 | 0.008 to 0.010 |

The LDOS is what's reported, rather than $Q/V$, because it measures what actually reaches the source. $V$ is measured at the field's strongest point, which needn't be at the source. Once a mode's own $Q$ is far above the bandwidth, the LDOS stops growing with it, and here the mode's volume caps it near 1,070.

Drawn as ordinary rectangles, and found with chapter 11's separate search for leaky modes, the default start has an LDOS of 1,066 at $Q = 1000$, a $Q$ of about 1,512,000 and $V = 0.0079$ wavelengths². Those change by under 1% from cells of $\lambda/24$ to $\lambda/72$, and $Q$ stays within 0.3% for absorbing layers of 0.5 to 1.5 wavelengths.
<!-- code: cavityShapes -->
<!-- code: optimizeCavity "Checks it as a drawn structure" .. "const found" -->

A $Q$ in the millions from so small a footprint deserves a check that shares nothing with the pole search. Outside the structure, the field is a sum of outgoing cylindrical waves, $u = \sum_m a_m H_m(kr) e^{im\theta}$ ([Bohren & Huffman 1983](20-references.md#bohren-huffman-1983)). Their power gives $Q = k^2 \int \epsilon |u|^2 \,/\, (4 \sum_m |a_m|^2)$, and that agrees within 0.5%, even at 5.9 million.

<a id="fabrication"></a>

## Fabrication

A real structure's edges come out slightly wrong, and the more sensitive a design, the more that matters. For a resonance, the sensitivity has a closed form. Moving a boundary shifts the frequency in proportion to the field on that boundary, by perturbation theory for shifting material boundaries ([Johnson et al. 2002](20-references.md#johnson-2002)). For $E_z$, which lies along every edge, that's

$$\frac{\delta k}{k} = -\frac{1}{2}\,\frac{\Delta\epsilon \oint u^2 \,dl}{\int \epsilon u^2\,dA}$$

per unit of outward move. To test a design against it, the edges have to move by less than a cell, so each pixel is drawn as its own rectangle and the mesh's nodes move onto the true edge (chapter 9). On a rod crystal, the test `moving_an_edge_shifts_a_resonance_as_perturbation_theory_says` agrees within 0.03%.

At a wavelength of 1 µm, moving every edge outwards shifts the optimized cavity's wavelength by 1.3% per nanometre, against a bandwidth of 0.1%. So a design that's 0.5 nm off reads an LDOS of 6 at its intended wavelength, rather than 1,000. The crystal shifts by 0.32% per nanometre, and misses as well, reading 18 at 1 nm. A cavity this sharp would need tuning after it was made. Topology optimization can also be made to anticipate such errors, by optimizing eroded and dilated copies of the design together ([Wang et al. 2011](20-references.md#wang-2011); [Christiansen & Sigmund 2021](20-references.md#christiansen-sigmund-2021)). Chapter 19 does something similar with random errors.

A design with fine lines is also harder to make. Removing the parts thinner than a given width from finished designs, medians over nine starts:

| Removed | LDOS | $Q$ |
|---|---:|---:|
| Nothing | 1,006 | 533,000 |
| Below 83 nm | 943 | 22,000 |
| Below 125 nm | 16 | 136 |

The crystal's rods are 151 nm across. Chapter 13's bend is far more forgiving. With every pixel edge moved at random, three starts keep 97.8% to 98.3% at a spread of 1 nm and 94.7% to 98.7% at 5 nm. Removing what's thinner than 83 nm keeps 95% to 98%, while removing what's thinner than 125 nm leaves 5% to 87%. A shift that moves a resonance 0.1% wide off its target barely matters to a bend that serves a band of 10%.

## What this leaves out

- **This is 2D.** A real slab also radiates out of the plane, and materials absorb. Either could cap $Q$ well below these values. [Albrechtsen et al. 2022](20-references.md#albrechtsen-2022) built measured fabrication limits into topology optimization of a cavity in a silicon membrane, made it, and report a $Q$ of about 1,100.
- **No minimum feature size is enforced**, beyond the filter's. Methods that guarantee one exist ([Wang et al. 2011](20-references.md#wang-2011)), and would cost some of the LDOS.

Next: [when numbers lie](18-when-numbers-lie.md).
