# 13. Inverse design

So far we've analysed structures that someone drew. Instead, we can ask for the structure: state what the light
should do, and let the solver find the shape. This is **topology optimization**, and it works because a solver can
report not only an answer, but how that answer changes with every part of the structure, at almost no extra cost.

## Densities

Give each pixel of a region a density $\rho$ between 0 and 1, mixing two materials:

$$\varepsilon = \varepsilon_\text{low} + \rho\,(\varepsilon_\text{high} - \varepsilon_\text{low}).$$

That turns a question about shapes, which can't be differentiated, into a smooth function of thousands of numbers,
which can. We then need a final design that's only 0s and 1s, since a real structure has no grey material. The
rest of this chapter is about the gradient, then about getting from grey to sharp.

The page does this for layers on the **Mirror** and **Broadband coating** views, and in 2D on the **2D waveguide bend**
and **2D cavity** views under **Inverse design**.

## The adjoint method

Hundreds of pixels would take hundreds of solves by finite differences. The **adjoint method** needs one more,
however many pixels there are.

A monitor's power from chapter 12 is a quadratic form, $P = u^H B\,u$, with $B$ Hermitian. A change to the matrix
changes the field by $A\,\delta u = -\delta A\,u$, so

$$\delta P = 2\,\mathrm{Re}(u^H B\,\delta u) = -2\,\mathrm{Re}\big(u^H B\,A^{-1}\delta A\,u\big).$$

Computing $A^{-1}\delta A\,u$ for every pixel is exactly the expensive part. But read from the left instead: only the
one row vector $u^H B\,A^{-1}$ is needed, and that's a single solve. Taking the transpose gives

$$\delta P = -2\,\mathrm{Re}(w^T \delta A\,u), \qquad A^T w = \overline{B\,u}.$$

One solve for $w$ prices every pixel. It costs even less than that sounds. Chapter 11 noted that $A$ is complex symmetric,
even inside the absorbing layer, so $A^T = A$, and $w$ only needs a back-substitution with the factors we already have.
For a pixel, $\delta A$ is just $-k_0^2\,\delta\varepsilon$ times its mass matrix, so each pixel's slope is a small
sum over its cells. `density_gradient_matches_finite_differences` checks this.

### The LDOS adjoint needs no solve at all

Chapter 11's LDOS is even kinder. With $A u = ik\,e_0$ and $F = -2\,\mathrm{Re}(e_0^T u)$, the same steps give

$$\delta F = 2\,\mathrm{Re}\big(e_0^T A^{-1}\delta A\,u\big) = 2\,\mathrm{Re}\big((A^{-1}e_0)^T \delta A\,u\big).$$

But $A^{-1}e_0$ is just $u/(ik)$, which we already have. A source's adjoint is its own field, because the Green's function is
reciprocal. So

$$\delta F = 2\,\mathrm{Re}\!\left(\frac{u^T \delta A\,u}{ik}\right),$$

and the gradient costs one more sum over cells. `ldos_gradient_matches_finite_differences` checks this.

## From grey to sharp

All of this lives in `web/js/design.js`, outside the solver. Only the physics needs to be in Rust.

**Filter.** If every pixel were free, the optimizer would make patterns at the scale of a single pixel, relying on
the mesh's errors. A **cone filter** replaces each density by a weighted average of its neighbours, with a weight
falling linearly to zero at a chosen radius. That sets a smallest feature size. It's linear, so its gradient
is its transpose: the same weights, carrying slopes the other way.

**Projection.** Then push the filtered densities toward 0 or 1:

$$\hat\rho = \frac{\tanh(\beta/2) + \tanh\big(\beta(\tilde\rho - \tfrac12)\big)}{2\tanh(\beta/2)}.$$

At small $\beta$ this is nearly a straight line. As $\beta$ grows it becomes a step. Starting with a sharp step would make
the problem nearly impossible, since every pixel's slope would be zero or huge, so $\beta$ doubles in stages, and the
structure sharpens gradually. The slope of $\hat\rho$ joins the chain rule.

**Optimizer.** With thousands of unknowns and an expensive function, a good step matters. The **method of moving asymptotes**
([Svanberg 1987](15-references.md#svanberg-1987)) replaces the objective with a simple convex model for each variable,
built from its slope and two asymptotes on either side, and minimizing that has a closed form. A variable that keeps
oscillating gets tighter asymptotes, and a steady one gets looser ones. Since projection amplifies a change in density
by up to about $\beta/2$, the largest step allowed shrinks as $\beta$ grows.

**Pixels on cells.** Each pixel should be a whole number of mesh cells, and the solver refuses anything else. If a pixel's
edge runs through a cell, that cell's material depends on where its centre happens to fall. Optimizers are very good at
exploiting such things, and no real structure reproduces them. In an earlier version, a design's resonance moved by 2% when
drawn as ordinary rectangles, and it looked far better on the coarse mesh than it was. With pixels on cells, a finished design
is exactly the same matrix as the same shapes drawn by hand.

### Example: a corner in a slab

Take a slab of index 3.5 and width $0.125\lambda$, with a monitor far along the output. A plain sharp corner passes 23 to 34% of the
light. Now cut out a square one wavelength wide at that corner, which leaves two guides ending at an empty region that passes about
5%, and give that square to the optimizer, with three wavelengths to serve at once
([Jensen & Sigmund 2004](15-references.md#jensen-sigmund-2004) did this kind of thing for photonic crystal bends). A
pattern from a longer run gives:

| Wavelength | Empty region | Sharp corner, 72 cells per $\lambda$ | Optimized, 24 cells per $\lambda$ | Optimized, 48 cells per $\lambda$ |
|---:|---:|---:|---:|---:|
| 0.95 | 5.4% | 23.0% | 97.8% | 97.6% |
| 1.00 | 5.3% | 28.6% | 99.3% | 99.3% |
| 1.05 | 5.5% | 34.2% | 99.4% | 99.6% |

The page's monitor is closer to the corner, where it still catches some light radiated from the gap: moving only the monitor from there to the far position takes the empty region from 7.1% to 5.4%, and the page's sharp corner reads 35%. The page also uses fewer iterations, 10 for each of 8 sharpness stages, and reaches a mean of 98.1% to 98.4% over three random starts.
Giving it 25 per stage reached 96.9% to 99.0%, so more iterations can reach higher but aren't more dependable. Skipping every
other sharpness stage did worse than cutting iterations, reaching only 96.5% to 97.6%. How the design sharpens matters more
than how long it runs.

For comparison, chapter 12's crystal bend reaches 97.9% to 95.8% over a much wider band, about 13%, in a crystal several periods
across. So neither wins outright: this pattern is compact and needs no band gap, while the crystal gets its bandwidth from having
nowhere else for light to go. We only optimized over 10% here, and over 0.85 to 1.15 this reached about 80 to 85%.

Optimizers are known to exploit a mesh's errors, so `optimized_bend_survives_a_finer_mesh` checks the finished pattern at twice the
resolution, and it holds. Around it, the power entering and leaving a closed loop balances to 0.02%. Optimizing for a single wavelength
reached only 43% to 82% at the band's edges, which is why this asks for three.

The structure and its mesh scale with the middle of the band, so 1.9 to 2.1 µm gets the same design as 0.95 to 1.05 µm, at twice the
size. Chapter 12 explains why that's exact.

## Designing a cavity

Now the harder problem. Maximize the LDOS at a line source in the middle of a square, with Liang & Johnson's idea
([Liang & Johnson 2013](15-references.md#liang-johnson-2013)): use chapter 11's complex frequency, and start at a low $Q$, whose
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
