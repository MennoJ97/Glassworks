# 19. Designing in the other polarization

*Going further. This builds on chapters 11, 12 and 13, and on chapter 17 for the cavity. Try it by choosing $H_z$ on the **2D waveguide bend** and **2D cavity** views.*

Chapters 13 and 17 design for $E_z$, the electric field out of the plane. But chips usually carry the other one. A silicon-on-insulator waveguide normally runs on its TE-like mode, with the electric field in the plane of the chip, and in a top-down 2D model that's $H_z$. So how well does a design for one serve the other?

Not well. Here's a bend optimized for $E_z$, read in both, on cells of $\lambda/72$:

| Transmission at $\lambda$ of 0.95, 1, 1.05 | $E_z$ | $H_z$ |
|---|---:|---:|
| Optimized for $E_z$ | 98.3, 98.1, 97.5% | 48.0, 43.8, 38.8% |
| Plain sharp corner | 24.0, 29.2, 34.6% | 34.5, 27.0, 20.2% |

$H_z$ needs its own designs. It's still one scalar field, so this needs none of chapter 14's vector elements. It does change what a pixel means, and three things that $E_z$ let us ignore.

## The gradient

A density mixes permittivities, and in $E_z$ that's $q = \varepsilon$, in the mass term. In $H_z$ it's $p = 1/\varepsilon$, in the stiffness, so a pixel's change to the matrix is

$$\delta A = -\frac{\delta\varepsilon}{\varepsilon^2}\,S.$$
<!-- code: pixel_slopes "Polarization::Tm" .. "deps / (eps * eps)" -->

The adjoint from chapter 13 is otherwise unchanged. There's no $k_0^2$ here, and each slope is now a sum of $\nabla w\cdot\nabla u$ over the pixel's cells. `a_monitors_gradient_matches_finite_differences` in [`tests/hz_design.rs`](../crates/eigenmenno-core/tests/hz_design.rs) checks this.

## A bend needs a wider guide

Copying the $E_z$ setup gave only 79 to 85%, and more iterations reached only 85 to 94%. The optimizer wasn't the problem. In $H_z$ the electric field crosses the slab's walls, and it's weak inside a high index, so a thin slab barely holds the mode. At a width of $0.125\lambda$, the effective index is 2.73 for $E_z$, and only 1.32 for $H_z$.

Twice as wide gives $H_z$ an effective index of 2.92, which asks a comparable thing of the design. It also carries a weak second mode, with an index of 1.10, and a monitor's total power wouldn't show light converted into it. But two modes would beat along the guide, and over six starts, transmission varies by under 2% from 2.1 to 3.3 wavelengths along it, against about 0.5% for $E_z$. So very little light is in the second mode.
<!-- code: bendWidth -->

| | Plain sharp corner | Six starts |
|---|---:|---:|
| $H_z$, width $0.25\lambda$ | 33% | 94.1 to 97.8% |

These are means over the band on the page's mesh, with its closer monitor. $H_z$ also needs more iterations: 40 per stage, where $E_z$ takes 10.
<!-- code: bendIterations -->

## An emitter that's electric

A cavity's goal is the LDOS at its emitter. In $E_z$ a line source is an electric current, as an emitter should be. In $H_z$ a line source is a *magnetic* current, and its LDOS says little about an atom or a quantum dot. So we need an in-plane electric dipole. It reads the electric field along its direction $\mathbf{d}$, which $H_z$ gives through a derivative:

$$\mathbf{E}\cdot\mathbf{d} = \frac{i}{k}\,p\,\frac{\partial H_z}{\partial n},$$
<!-- code: dipole.rs "let mut load = vec![]" .. "Ok(load)" -->

with $\mathbf{n}$ a quarter turn from $\mathbf{d}$. So it's a new $e_0$ with the same algebra: $A u = ik\,e_0$, and $F = -2\,\mathrm{Re}(e_0^T u)$. [`solve2d/dipole.rs`](../crates/eigenmenno-core/src/solve2d/dipole.rs) takes the derivative across two points, and scales it so that $F = k/2$ in vacuum, as for a line source. `the_power_a_dipole_supplies_crosses_any_loop` checks it against energy balance, which knows nothing of that scale.
<!-- code: dipole.rs "let step" .. "Ok([point" -->

Chapter 17's shortcut, that the source is its own adjoint, still holds. But $e_0$ now holds $p$ at the emitter, so a pixel there moves the source as well as the field. That adds a term,

$$\delta F \supset -4\,\mathrm{Re}\big(u^T \delta e_0\big),$$
<!-- code: add_source_slope "let value" .. "imaginary[pixel]" -->

which the finite-difference test catches if it's removed: 1.26 against the right 4.99 in one pixel.

## Two things that broke the cavity

**The frequency hold lost its grip.** Chapter 17 places a pole with two solves, assuming that $1/u_0$ is a straight line in $k$. A dipole's response also holds a large, nearly constant background from its near field, since the derivative of a Green's function on itself diverges as the points come together. Then $u_0 \approx B + R/(k - k_p)$, which isn't a straight line until a design is very sharp. The mode swung by $\pm 3\%$ between iterations, and the LDOS fell to under 1. Three solves fit $B$ as well, and that fixes it.
<!-- code: fitPole -->

**Grey pixels matter more.** Sharpening went well, then removing the last grey pixels moved the mode, so that finished designs needed resizing by up to 10.4%, where $E_z$ needs 1.2%. That was outside the search for the finished peak. Here the field crosses each pixel's edge, so its value there changes the mode far more. The search now looks further. We also tried holding the frequency through two more sharpness stages. Three starts looked better and nine worse, with a median of 345 against 437, so it was dropped.
<!-- code: optimizeCavity "Grey pixels matter more" .. "0.96 + 0.004" -->

## Holding the frequency inside the step

Even with three solves, the chart swung by up to 100 times late in a run. The step was the cause. Chapter 17 removes the frequency's slope from the goal, lets MMA step, and then corrects along that slope. But MMA with only bounds acts nearly as a sign step: each pixel's model is monotone, so it moves to its limit however small its slope. Here 81 to 90% of pixels took the full move, and the step had a cosine of only 0.24 with the direction asked for. So removing the slope from the direction didn't remove it from the step. Each step moved the mode by 5 to 6 bandwidths, more than the correction could undo.

So in $H_z$ the frequency is held inside the step. MMA's step for a goal $g$ has a closed form for each pixel, so trying a goal of $g + \lambda f$, where $f$ is the frequency's slope, costs no solve. A bisection on $\lambda$ finds the step whose shift, carried through the filter and projection, is half the way to the target. The next iteration then reads the LDOS where that predicts the peak. Topology optimization has long held a volume this way, with the same bisected multiplier ([Andreassen et al. 2011](20-references.md#andreassen-2011)), and MMA treats a constraint the same way through its dual. Here it holds a frequency instead. Each pixel flips between its limits as its part of $g + \lambda f$ changes sign, so the shift comes in small jumps, and the hold is good to within one pixel's jump.
<!-- code: HeldMma "const stepFor" .. "const next = choose" -->
<!-- code: holdingStep -->

Others have found the same difficulty. The LDOS's sharpest curvature lies along the frequency's slope, so a high $Q$ makes the goal stiff in just that direction ([Shaker et al. 2025](20-references.md#shaker-2025)). They follow the mode with an eigensolve at each step instead, which costs about twice the solves.

A prediction through the projection is only as good as the projection's slope, and a $\tanh$ projection's slope grows without limit as it sharpens. A subpixel-smoothed projection ([Hammond et al. 2025](20-references.md#hammond-2025)) instead works out from the filtered field's value and slope how far an edge lies through each pixel. So an edge stays about a pixel wide at any sharpness, and its slope stays finite. Its second-order form ([Romano et al. 2026](20-references.md#romano-2026)) adds the curvature, which keeps it smooth where two edges meet. `smoothProject` and `holdingStep` in `web/js/hz-design.js` do both.
<!-- code: smoothProject "const gx" .. "rho[k] = smoothed" -->
<!-- code: smoothed -->

Each change helps alone. Over nine starts with disorder sampled as below, the median LDOS goes from 665 to 1,248 with the held step, and to 1,028 with the smoothed projection. Together they reach 1,348. The chart still drops on the first iteration of each stage, where sharpening moves the mode before any step could predict it, and it recovers on the next. $E_z$ keeps chapter 17's step, since none of this was measured there.

## The result

Rods are no benchmark here, since they read about 1. The usual crystal for this polarization is a triangular lattice of air holes, with the middle one filled. We scanned its period and hole size, in the same square, and the best is $a = 0.276\lambda$, $r = 0.4a$. All are at $Q = 1000$, relative to vacuum:
<!-- code: holeCrystal "const a" .. "if ((i || j) && inside)" -->

| Structure | LDOS | Mode $Q$ |
|---|---:|---:|
| Crystal of rods | 1 | |
| Crystal of air holes | 371 | |
| Nine optimized starts | 908 to 1,660, median 1,348 | 1,001 to 30,071 |

So the median beats the crystal by 3.6 times, and the weakest start by 2.4. Chapter 17's $E_z$ designs beat their own crystal by at least 2.9 times, though against a weaker kind of benchmark. On cells half the size, 18 other starts have a median of 1,370, and the crystal reads 382. An $H_z$ LDOS and an $E_z$ one measure different emitters, so comparing them directly means little.

<a id="what-this-doesnt-settle"></a>

## What this doesn't settle

- **The index is still bulk.** These designs use 3.5 and $\varepsilon = 12$. A real 2D model of a silicon slab uses its effective index, about 2.86 for 220 nm ([Nikkhah et al. 2024](20-references.md#nikkhah-2024)), and the third dimension leaks light that none of this sees.
- **These designs are better, not sturdier.** Over 18 starts, on half cells, with edges disturbed randomly and correlated over 100 nm, the median LDOS falls from 1,370 to 1,137 at 2 nm, and to 639 at 5 nm, against the crystal's 382, 335 and 269. At 5 nm, 17 of the 18 still beat it. But they keep less of their LDOS than the crystal does, 0.47 at 5 nm against 0.70, and they shift further, 0.96% per nanometre against 0.59, so they would still need tuning. They also tend to form slots: air channels about 87 nm wide crossing at the source, where $H_z$'s electric field is strongest. That's narrow for lithography, and minimum features weren't measured. The crystal is one structure with 18 draws, and its disorder is in each hole's radius, so this comparison is rough.
- **Sampling disorder while optimizing helped the earlier step, mostly as a better search.** The page takes three samples per iteration, with edges moved by 2 nm and correlated over 100 nm. With chapter 17's step, over 18 starts, it raised the median LDOS from 417 to 664 with no disorder, and from 197 to 363 at 5 nm. But the share kept at 5 nm only rose from 0.43 to 0.51. Whether it still helps with the held step wasn't measured. It makes a run take about twice as long.
- **The dipole's direction.** Without disorder, a design's eight-fold symmetry treats both in-plane directions alike: over 18 starts they agree within 1 part in 100,000. Disorder breaks that. With independent disorder on each edge, the two directions differ by a median of 6% at 2 nm, and 18% at 5 nm, and in the worst draws by 38% and 79%.
- **Few starts,** with spreads wide enough that small differences here shouldn't be trusted. A run is also chaotic, since one pixel flipping changes the path. Changing a gradient in its last digit, by adding its terms in another order, moved one start from 1,395 to 719. The same start reads 1,323 in Node and 1,333 in a browser, because their maths functions differ in the last digit. Each repeats exactly on its own engine. The numbers above are from Node, so the page's own will differ within that spread.

Next: [references](20-references.md).
