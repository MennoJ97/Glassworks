# 19. Designing in the other polarization

*Going further. This builds on chapters 11, 12 and 13, and on chapter 17 for the cavity. Try it by choosing $H_z$ on the **2D waveguide bend** and **2D cavity** views.*

Chapters 13 and 17 design for $E_z$, the electric field out of the plane. But chips usually carry the other one. A silicon-on-insulator waveguide normally runs on its TE-like mode, with the electric field in the plane of the chip, and in a top-down 2D model that's $H_z$ ([Nikkhah et al. 2024](20-references.md#nikkhah-2024)). Topology optimization treats both the same way, with $\varepsilon$ moving from the mass term to the stiffness ([Christiansen & Sigmund 2021](20-references.md#christiansen-sigmund-2021)), but a design for one doesn't serve the other. Here's a bend optimized for $E_z$, read in both, on cells of $\lambda/72$:

| Transmission at $\lambda$ of 0.95, 1, 1.05 | $E_z$ | $H_z$ |
|---|---:|---:|
| Optimized for $E_z$ | 98.3, 98.1, 97.5% | 48.0, 43.8, 38.8% |
| Plain sharp corner | 24.0, 29.2, 34.6% | 34.5, 27.0, 20.2% |

$H_z$ is still one scalar field, so this needs none of chapter 14's vector elements. It changes what a pixel means, what a guide can hold, and what an emitter is.

## The gradient

A density mixes permittivities, and in $E_z$ that's $q = \varepsilon$, in the mass term. In $H_z$ it's $p = 1/\varepsilon$, in the stiffness, so a pixel's change to the matrix is

$$\delta A = -\frac{\delta\varepsilon}{\varepsilon^2}\,S.$$
<!-- code: pixel_slopes "Polarization::Tm" .. "deps / (eps * eps)" -->

The adjoint from chapter 13 is otherwise unchanged. There's no $k_0^2$ here, and each slope is now a sum of $\nabla w\cdot\nabla u$ over the pixel's cells. `a_monitors_gradient_matches_finite_differences` in [`tests/hz_design.rs`](../crates/glassworks-core/tests/hz_design.rs) checks this.

## A bend needs a wider guide

In $H_z$ the electric field crosses the slab's walls, and since $\varepsilon E$ is continuous across them, the field is $\varepsilon$ times weaker inside the slab than just outside ([Jackson 1999](20-references.md#jackson-1999)). So a thin slab holds this polarization far more loosely. The slab dispersion equations of chapter 8 show it ([Yariv & Yeh 2007](20-references.md#yariv-yeh-2007)). At a width of $0.125\lambda$, the effective index is 2.73 for $E_z$ and only 1.32 for $H_z$. Twice as wide gives $H_z$ an effective index of 2.92, close to what $E_z$ had, so the page's $H_z$ bend uses that width.
<!-- code: bendWidth -->

The wider guide also carries a weak second mode, with an index of 1.10. A monitor's total power wouldn't show light converted into it, but two modes would beat along the guide, and transmission varies by under 2% from 2.1 to 3.3 wavelengths along it. So very little light is in the second mode. In chapter 13's tighter turn, with a curved bend of radius $0.25\lambda$ as the baseline, three starts reach 98.3% to 98.5%, against 93.2% for the curve and 33.8% for a sharp corner. $H_z$ needs more iterations, 40 per stage rather than 10.
<!-- code: bendIterations -->

## An emitter that's electric

A cavity's goal is the LDOS at its emitter. In $E_z$ a line source is an electric current, as an emitter should be. In $H_z$ a line source is a *magnetic* current, and its LDOS says little about an atom or a quantum dot. So we need an in-plane electric dipole. It reads the electric field along its direction $\mathbf{d}$, which $H_z$ gives through a derivative, from the curl equation of chapter 1:

$$\mathbf{E}\cdot\mathbf{d} = \frac{i}{k}\,p\,\frac{\partial H_z}{\partial n},$$
<!-- code: dipole.rs "let mut load = vec![]" .. "Ok(load)" -->

with $\mathbf{n}$ a quarter turn from $\mathbf{d}$. So it's a new $e_0$ with the same algebra: $A u = ik\,e_0$, and $F = -2\,\mathrm{Re}(e_0^T u)$. [`solve2d/dipole.rs`](../crates/glassworks-core/src/solve2d/dipole.rs) takes the derivative across two points, and scales it so that $F = k/2$ in vacuum, as for a line source. `the_power_a_dipole_supplies_crosses_any_loop` checks it against energy balance, which knows nothing of that scale.
<!-- code: dipole.rs "let step" .. "Ok([point" -->

Chapter 17's shortcut, that the source is its own adjoint, still holds. But $e_0$ now holds $p$ at the emitter, so a pixel there moves the source as well as the field. That adds a term,

$$\delta F \supset -4\,\mathrm{Re}\big(u^T \delta e_0\big),$$
<!-- code: add_source_slope "let value" .. "imaginary[pixel]" -->

which the finite-difference test catches if it's removed: 1.26 against the right 4.99 in one pixel.

## Finding the mode

Chapter 17 places a pole with two solves, assuming that $1/u_0$ is a straight line in $k$. A dipole's response also holds a large, nearly constant background from its near field. Near a resonance, a response is the pole's term plus a smooth, non-resonant background ([Bykov & Doskolovich 2013](20-references.md#bykov-doskolovich-2013)),

$$u_0 \approx B + \frac{R}{k - k_p},$$

which isn't a straight line in $1/u_0$ until a design is very sharp. So three solves fit $B$ as well.
<!-- code: fitPole -->

Removing the last grey pixels moves this mode further than it moves $E_z$'s. Here the field crosses each pixel's edge, and a boundary's effect on a field normal to it is much harder to capture than on a field along it ([Johnson et al. 2002](20-references.md#johnson-2002)). Finished designs needed resizing by up to 10.4%, where $E_z$ needs 1.2%, so the search for the finished peak looks further.
<!-- code: optimizeCavity "Grey pixels matter more" .. "0.96 + 0.004" -->

Rounding costs strength as well. Read at the last stage's bandwidth, $Q = 500$, just before and just after it, three starts keep 81 to 91% of their LDOS at a wavelength of 1 µm, and 37 to 86% at 1.5 µm, where $E_z$ keeps 98%. The page's chart ends at that bandwidth, before the rounding, and then reports the result at $Q = 1000$, which reads higher for a sharp mode, so on the chart the loss mostly hides. Its cause is the same as the shift's. Neither a grey pixel nor a rounded one is the edge between pixels that the optimizer was steering towards, and a mode this sharp needs its edges right to about a nanometre, where a pixel is $\lambda/24$, 42 nm at 1 µm. More iterations after rounding, edges drawn more finely, and grey pixels averaged differently along and across an edge, as [Farjadpour et al. 2006](20-references.md#farjadpour-2006) do, were all tried, and none closed the gap reliably.

## Holding the frequency inside the step

Chapter 17 holds the frequency by projecting the goal's gradient, then letting MMA step. That relies on the step following the direction it was given. MMA with only bounds doesn't, quite: each pixel's model is monotone, so a pixel with any slope at all moves to its limit, and the step behaves nearly as a sign step. In $H_z$, 81 to 90% of pixels took the full move, and the projection's work was mostly undone.

The standard remedy is to put the constraint inside the subproblem that produces the step. MMA does that through the dual of its subproblem ([Svanberg 1987](20-references.md#svanberg-1987)), and topology optimization's optimality-criteria update does it for a volume constraint, by bisecting on the constraint's Lagrange multiplier ([Andreassen et al. 2011](20-references.md#andreassen-2011)). A projection doesn't preserve volume, so the direct approach evaluates the constraint through the filter and projection at each bisection step ([Ferrari & Sigmund 2020](20-references.md#ferrari-sigmund-2020)).
<!-- code: HeldMma "const stepFor" .. "const next = choose" -->

The $H_z$ cavity does the same for its frequency. MMA's step for a goal $g$ has a closed form for each pixel, so trying a goal of $g + \lambda c$, with $c$ the frequency's slope, costs no solve. A bisection on $\lambda$ finds the step whose predicted shift, carried through the filter and projection, is half the way to the target, as in chapter 17's range-space step. The next iteration then reads the LDOS where that predicts the peak, which is [Shaker et al. 2026](20-references.md#shaker-2026)'s idea of following the mode, with a first-order prediction in place of their eigensolver.
<!-- code: holdingStep -->

Each part of this is published, but we haven't found them put together: a frequency, rather than a volume, held by bisection inside MMA's step, half the way at a time, through the projection below, at no extra solve. So it's this project's own small variant, not a method from those papers. It's also a little cruder than MMA's own dual, which mixes the two goals' coefficients where this mixes their gradients. Tried here, that held the frequency exactly but halved the median LDOS, so the cruder step stays.

A prediction through the projection is only as good as the projection's slope, and a $\tanh$ projection's slope grows without limit as it sharpens. A **subpixel-smoothed projection** ([Hammond et al. 2025](20-references.md#hammond-2025)) instead works out, from the filtered field's value and slope, how far an edge lies through each pixel. So an edge stays about a pixel wide at any sharpness, and its slope stays finite. Its second-order form ([Romano et al. 2026](20-references.md#romano-2026)) adds the curvature, which keeps it smooth where two edges meet. `smoothProject` in `web/js/hz-design.js` does both.
<!-- code: smoothProject "const gx" .. "rho[k] = smoothed" -->
<!-- code: smoothed -->

## Random fabrication errors

Each iteration also averages the gradient over three copies of the design, with every edge moved at random, by 2 nm with a correlation length of 100 nm. That's the stochastic form of robust topology optimization, which optimizes the expected performance over a random field of manufacturing errors ([Schevenels et al. 2011](20-references.md#schevenels-2011)). Each copy is read at the wavelength its own predicted shift gives, for the reason in chapter 17: read at the target, a shifted copy sits on its peak's wing. Three samples are far too few to estimate an expectation well, so this mostly steers the search, and it makes a run take about twice as long.

## The result

Rods are no benchmark here, since they read about 1. The usual crystal for this polarization is a triangular lattice of air holes ([Joannopoulos et al. 2008](20-references.md#joannopoulos-2008)), with the middle one filled. Its period and hole size were scanned in the same square, and the best is $a = 0.276\lambda$, $r = 0.4a$. All are at $Q = 1000$, relative to vacuum:
<!-- code: holeCrystal "const a" .. "if ((i || j) && inside)" -->

| Structure | LDOS | Mode $Q$ |
|---|---:|---:|
| Crystal of rods | 1 | |
| Crystal of air holes | 371 | |
| Nine optimized starts | 908 to 1,660, median 1,348 | 1,001 to 30,071 |

The median beats the crystal by 3.6 times, and the weakest start by 2.4. An $H_z$ LDOS and an $E_z$ one measure different emitters, so comparing them with chapter 17's directly means little.

<a id="what-this-doesnt-settle"></a>

## What this doesn't settle

- **The index is still bulk.** These designs use 3.5 and $\varepsilon = 12$. A real 2D model of a silicon slab uses its effective index, about 2.86 for 220 nm of silicon in oxide at a wavelength of 1.525 µm ([Nikkhah et al. 2024](20-references.md#nikkhah-2024)), and the third dimension leaks light that none of this sees.
- **These designs are better, not sturdier.** With edges disturbed randomly by 5 nm, the median LDOS over 18 starts falls to 639, against the crystal's 269, but the designs keep less of their LDOS than the crystal does, and shift further, 0.96% per nanometre against 0.59%. They would still need tuning. They also tend to form air slots about 87 nm wide at the source, where $H_z$'s electric field is strongest, which is narrow for lithography.
- **Few starts.** A run is chaotic, since one pixel flipping changes the path. Changing a gradient in its last digit, by adding its terms in another order, moved one start from 1,395 to 719, so small differences between versions here shouldn't be trusted. The numbers above are from Node, and a browser's maths functions differ in the last digit, so the page's own will differ within that spread.

Next: [references](20-references.md).
