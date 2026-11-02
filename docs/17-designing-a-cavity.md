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
([Liang & Johnson 2013](19-references.md#liang-johnson-2013)): use chapter 11's complex frequency, and start at a low $Q$, whose
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

### What we don't know

These results deserve some caution.

- **A $Q$ in the millions from a footprint this small has no explanation here.** Checks against the mesh, layer and surroundings all hold, so it isn't an obvious
  numerical artifact, but we haven't found the physical cause.
- **No fabrication tolerance has been checked.** A mode this sharp would be very sensitive. A real device would need the optimization to be robust against small
  changes, which this one isn't asked to be.
- **Comparing these at $Q = 1000$ hides differences.** Once a mode's own $Q$ is far above the bandwidth, the LDOS stops growing with it, and is capped near 1,070 here by the mode's volume. That's why we also quote higher bandwidths.

Next: [beyond scalar fields](14-vector-fields.md).

Next: [when numbers lie](18-when-numbers-lie.md).
