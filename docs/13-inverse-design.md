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
view under **Inverse design**. Chapter 17 designs a cavity the same way.

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

## From grey to sharp

All of this lives in `web/js/optimizer.js`, outside the solver. Only the physics needs to be in Rust.

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
([Svanberg 1987](19-references.md#svanberg-1987)) replaces the objective with a simple convex model for each variable,
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
([Jensen & Sigmund 2004](19-references.md#jensen-sigmund-2004) did this kind of thing for photonic crystal bends). A
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

That completes the essentials: a field, a mesh, a solver, open boundaries, and a design. The next chapters go further, and each stands on its own: [start with whichever interests you](README.md#going-further).

Next: [beyond scalar fields](14-vector-fields.md).
