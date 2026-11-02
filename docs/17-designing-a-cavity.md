# 17. Designing a cavity

*Going further. This builds on the LDOS from chapter 11 and the optimizer from chapter 13. It's closer to a research log than the rest of the tutorial, because the first attempts failed in instructive ways.*

Chapter 11 found the modes of a crystal someone drew. What if the solver draws the cavity? The code is in `web/js/cavity-design.js`, on the **2D cavity** view under **Inverse design**.

## The LDOS adjoint needs no solve at all

Chapter 11's LDOS is even kinder. With $A u = ik\,e_0$ and $F = -2\,\mathrm{Re}(e_0^T u)$, the same steps give

$$\delta F = 2\,\mathrm{Re}\big(e_0^T A^{-1}\delta A\,u\big) = 2\,\mathrm{Re}\big((A^{-1}e_0)^T \delta A\,u\big).$$

But $A^{-1}e_0$ is just $u/(ik)$, which we already have. A source's adjoint is its own field, because the Green's function is
reciprocal. So

$$\delta F = 2\,\mathrm{Re}\!\left(\frac{u^T \delta A\,u}{ik}\right),$$

and the gradient costs one more sum over cells. `ldos_gradient_matches_finite_differences` checks this.

## The first attempts

Maximize the LDOS at a line source in the middle of a square, with Liang & Johnson's idea
([Liang & Johnson 2013](20-references.md#liang-johnson-2013)): use chapter 11's complex frequency, and start at a low $Q$, whose
smooth average gives the optimizer something to climb. Raising $Q$ in stages narrows the bandwidth, and leads it toward a sharp
resonance. The gradient is the one that needs no extra solve.

Our benchmark is chapter 11's crystal of 7 × 7 rods with $\varepsilon = 12$, in a region of the same size, 2.67 wavelengths wide, in 64 × 64
pixels. The run has 8 stages of 20 iterations, ending at $Q = 500$ and $\beta = 128$. This is the honest comparison: a crystal is a good cavity. And it starts out badly. It's worth seeing how, since each failure
teaches something.

**The resonance moved.** Finished designs were strong, but with peaks 2 to 5% away from the target. At $Q = 1000$ the bandwidth is 0.1%, so
they miss it completely, and the optimizer spends its effort chasing a peak that moves whenever it takes a step.

**The filter was too wide.** A filter of 3 pixels, 0.125 wavelengths, is nearly as wide as these rods. Applied to the crystal itself,
it left a fill of 5.4% rather than the true 12.3%. The design space couldn't even hold our benchmark. 2 pixels can.

**Q/V wasn't the right number to read.** One early design had $Q = 1,889$ and $V = 0.0137\lambda^2$, more than 7 times the crystal's $Q/V$.
But its LDOS at its own resonance was only 15, against the crystal's 299, at $Q = 1000$. $V$ is measured at the field's strongest point, which
may not be at the source. The LDOS says what actually reaches the source, so that's what we report.

### Holding the frequency still

The first of those has a satisfying fix. The trouble is that a step good for the cavity also moves its resonance. So let's take steps
that don't.

Near a pole, $1/u_0$ is nearly linear in $k$, as chapter 11 showed, which gives us two things. First, two solves place the pole. Second, the
same complex quantity $u^T \delta A\,u / (ik)$ is the slope of $-u_0$, not only of its real part. Writing $1/u_0$ as a straight line, with the
slope $d(1/u_0)/dk$ from those two solves, a change of $u_0$ moves the pole by

$$\delta k_p = \frac{\delta(1/u_0)}{-\,d(1/u_0)/dk}.$$

So for every density $x_i$, we get a number $c_i = \partial\,\mathrm{Re}\,k_p/\partial x_i$. Measured on a true resonance, that predicts the real
shift within 2%. Together these numbers are a direction in the space of all densities: the one direction along which the frequency changes
fastest. Every direction perpendicular to it leaves the frequency alone, to first order. So remove that component from the gradient $g$:

$$g \;\to\; g - \frac{g\cdot c}{c\cdot c}\,c.$$

The step now improves the cavity without moving it. MMA has bounds, so its step isn't quite perpendicular, and any error builds up. So afterwards, a
correction moves along $c$ by just enough to aim half of the remaining frequency error back at the target. That makes "stay at this wavelength" a real
constraint on the design, rather than something to hope for.

Each step reads LDOS at the mode's current wavelength rather than at the target, for a reason the results below will show.

The last stages of sharpening still shift it a little, because a pixel flip is a big step. Scale invariance takes out what's left exactly: resize the
finished design by $\lambda_\text{target}/\lambda_\text{mode}$.

It's worth measuring, because one good run proves little here: starts differ by far more than most of these choices. This table covers nine
starts, at a bandwidth narrow enough for a strong mode to stand out. These runs read LDOS at the target, as the first version did.

| Version | Median LDOS at $Q = 10{,}000$ | Median mode $Q$ | Largest resizing |
|---|---:|---:|---:|
| Holding the frequency | 8,411 | 56,000 | 1.2% |
| Not holding it | 3,304 | 5,900 | 4.4% |

So the hold is worth about ten times in $Q$. Holding alone isn't enough, though: with no goal at all, LDOS at $Q = 1000$ reaches only 27 to 372. Removing only the correction looked worst of all on three starts, at
61 to 485 at $Q = 1000$, where the full method reads 110 to 1,061.

### The result

We compare everything at $Q = 1000$, relative to vacuum, and each at its own resonance. The crystal is tuned to its own, and every pattern goes through the same search for its strongest peak within 4%, then resizing:

| Structure | LDOS | $Q$ | $V$ (wavelengths²) |
|---|---:|---:|---:|
| 16 random symmetric patterns | 0 to 29 | | |
| Crystal of 7 × 7 rods | 299 | 1,680 | 0.093 |
| 18 optimized starts | 886 to 1,075, median 1,014 | 9,120 to 5,859,698 | 0.008 to 0.010 |

All 18 beat the crystal by at least 2.9 times. The crystal's $Q$ here is a little different from chapter 11's, because of its mesh and surroundings.
At $Q = 3000$, the page's default start reaches an LDOS of 3,193, and the crystal 508. On the earlier setup, raising $Q$ only as far as 500 while sharpening was the steadiest schedule: taking it
to 2,000 reached only 269 on one start.

Drawn as ordinary rectangles, and using the separate search for leaky modes from chapter 11, the default start has an LDOS of 1,066 at $Q = 1000$, $Q$ of about 1,512,000 and $V = 0.0079$ wavelengths². Those change by under 1% from cells of $\lambda/24$
to $\lambda/72$. Its $Q$ stays within 0.3% for layers of 0.5 to 1.5 wavelengths, and within 0.7% with 0.25 to 1.5 wavelengths of air around it.

### A goal that can work against itself

Over nine starts, flipping the goal's sign and minimizing the LDOS appeared to beat maximizing it. Nine more starts reversed the median, 951 against 1,001, so most of that was chance. The lesson is the sample size. But the flipped goal did fail less often, and looking into why turned up something real.

Near a pole, with detuning $D$, bandwidth $\Gamma$ and the mode's own decay $\gamma$,

$$F \propto \frac{\Gamma + \gamma}{D^2 + (\Gamma + \gamma)^2}.$$

Its slope along $\gamma$ changes sign at $|D| = \Gamma + \gamma$. Closer than that, raising $F$ raises $Q$. Further out, on the wing of the peak, raising $F$ lowers $Q$. Holding the frequency removes nearly all of the gradient, and the small remainder sets the step, so this matters. In 18 runs reading $F$ at the target, the remainder pointed towards more decay on 0.4% of iterations below half a linewidth, and 97% from one to three linewidths. Runs that spent more time past one linewidth also tended to finish weaker.

So we read $F$ at the mode's own wavelength, which stays close to the peak. Over those 18 starts, the median mode $Q$ rose from 113,000 to 374,000, and the worst LDOS at $Q = 1000$ from 110 to 886. These figures include a later change to the factorization, which moved three starts through round-off alone. This isn't the whole story: other runs lost their $Q$ while sharpening, even when close to the peak.

<a id="would-it-work-if-you-built-it"></a>

### Would it work if you built it?

A $Q$ in the millions from so small a footprint deserves suspicion, so we tested it. Everything here is at a wavelength of 1 µm, in 2D, and with no loss. Lengths scale with the wavelength.

**The $Q$ is real, in this model.** Outside the structure, the field splits into outgoing waves, $u = \sum_m a_m H_m(kr) e^{im\theta}$, and their power gives $Q = k^2 \int \epsilon |u|^2 \,/\, (4 \sum_m |a_m|^2)$. That shares nothing with the search for the pole, and it agrees within 0.5%, even at 5.9 million. Three things make the $a_m$ small. Eight-fold symmetry allows only $m = 0, 4, 8, 12, \dots$, and so small a structure barely couples to higher $m$, which leaves about four channels to quieten. Interference then quietens the lowest ones: for $m = 0$, the net radiation is $1.7 \times 10^{-4}$ of the sum of its parts, against $5.7 \times 10^{-2}$ for the crystal. And less field reaches the outside. A larger region raises $Q$ quickly, while $V$ and the LDOS stay almost the same.

**Moving an edge by less than a cell.** Fabrication errors are often quoted in single nanometres, against cells of 42 nm. Giving the edge's cell a fractional density looks natural, and is wrong by 15 to 39% in the shift it predicts. What works is drawing each pixel as its own rectangle, so that mesh nodes move onto the true edge. As a check, the shift should follow from perturbation theory, $\delta k / k = -\tfrac{1}{2} \Delta\epsilon \oint u^2 \,dl \,/\, \int \epsilon u^2$ per unit of outward move. On a rod crystal, the test `moving_an_edge_shifts_a_resonance_as_perturbation_theory_says` agrees within 0.03%, and on finished designs, within about 1 to 6%.

**An error moves the resonance out of its bandwidth.** Moving every edge shifts the wavelength by 1.3% per nanometre, against a bandwidth of 0.1%, so a design that misses by 0.5 nm has an LDOS of 6 at its intended wavelength rather than 1,000. The crystal shifts 0.32% per nanometre. It misses as well, with 18 at 1 nm.

The reason is where the source sits. Left alone, the optimizer does what you'd expect, and puts the source in high index: a rod about 83 nm across, with 87% of the energy in material. $F$ follows $|u_0|^2 / \int \epsilon |u|^2$, and that tiny rod makes a sharp peak cheaply. But its surface is 47 nm from the peak, and the shift goes with the field at the edges. So we also tried keeping a disc of 250 nm radius clear of material. The source then sits in air, where a broader peak is still cheap, since energy in air counts 12 times less. The LDOS falls 39%, and the shift falls 5 times, to 0.26% per nanometre. This is for $E_z$, which is continuous across edges. We haven't tried the other polarization, whose field jumps at an edge.

**Random disorder.** Each exposed pixel side now moves independently, with a Gaussian spread we chose, and a correlation length of one pixel. Over 18 starts, medians:

| Spread | Design | Mode $Q$ | LDOS, retuned | LDOS at the original wavelength |
|---:|---|---:|---:|---:|
| 0 | Rod | 393,000 | 1,005 | 1,005 |
| 0 | Clear disc | 306,000 | 599 | 599 |
| 1 nm | Rod | 216,000 | 1,001 | 38 |
| 1 nm | Clear disc | 87,000 | 595 | 525 |
| 5 nm | Rod | 22,000 | 928 | 2 |
| 5 nm | Clear disc | 12,000 | 530 | 153 |

The crystal, with each rod's radius varied instead, which is a different model, reads $Q$ of 1,373, 270 retuned and 23 at the original wavelength, at 5 nm. So $Q$ in the millions doesn't survive. The remaining $Q$ is still several times the crystal's. If the wavelength can be tuned, the rod is the better design. If it can't, the clear disc is far better, and random errors hurt it much less than a shared one, which still needs tuning.

**Optimizing for disorder.** Adding random samples to each step, each read at its own wavelength, raises the clear disc's $Q$ by about 2, 3 and 4 times at 1, 2 and 5 nm, on 12, 14 and 14 of 18 unseen starts. It doesn't reduce the shift. Reading them at the target wavelength instead destroys $Q$, through the wing effect above, and so does the textbook approach of optimizing eroded and dilated copies. It takes about twice the time.

**Feature sizes.** All of these designs are built from lines about 83 nm wide. Removing what's thinner than a given width from a finished design, medians over nine starts:

| Removed | Rod: LDOS | Rod: $Q$ | Clear disc: LDOS | Clear disc: $Q$ |
|---|---:|---:|---:|---:|
| Nothing | 1,006 | 533,000 | 614 | 296,000 |
| Below 83 nm | 943 | 22,000 | 601 | 56,000 |
| Below 125 nm | 16 | 136 | 7 | 79 |

Cleaning up a design afterwards isn't a fair test of designing with a limit, so we also widened the filter. It doesn't set a strict limit, because isolated thin parts remain, but it does remove most of them. With about 6% of the material thinner than 125 nm, a clear disc reaches a median LDOS of only 279, with half the starts below $Q = 10{,}000$. A larger disc of 417 nm does better, at 492, with 16% thinner than that, $Q$ of 18,000 and a shift of 0.44% per nanometre. The crystal's rods are 151 nm across. So how much of the advantage survives depends on how fine a process is, and on the optimizer as well as the physics.

### What we don't know

These results deserve some caution.

- **This is 2D.** A real slab also radiates out of the plane, and materials absorb. Either could cap $Q$ well below all of these.
- **The disorder is assumed.** Its spread and correlation length are our choices, not a measured process, and six random draws per start understate the worst cases.
- **Comparing these at $Q = 1000$ hides differences.** Once a mode's own $Q$ is far above the bandwidth, the LDOS stops growing with it, and is capped near 1,070 here by the mode's volume. That's why we also quote higher bandwidths.

The page designs the unconstrained rod. The clear disc and the disordered goal are experiments, in `plans/probe-tolerance`.

Next: [when numbers lie](18-when-numbers-lie.md).
