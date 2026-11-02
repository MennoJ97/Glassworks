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

The adjoint from chapter 13 is otherwise unchanged. There's no $k_0^2$ here, and each slope is now a sum of $\nabla w\cdot\nabla u$ over the pixel's cells. `a_monitors_gradient_matches_finite_differences` in [`tests/hz_design.rs`](../crates/eigenmenno-core/tests/hz_design.rs) checks this.

## A bend needs a wider guide

Copying the $E_z$ setup gave only 79 to 85%, and more iterations reached only 85 to 94%. The optimizer wasn't the problem. In $H_z$ the electric field crosses the slab's walls, and it's weak inside a high index, so a thin slab barely holds the mode. At a width of $0.125\lambda$, the effective index is 2.73 for $E_z$, and only 1.32 for $H_z$.

Twice as wide gives $H_z$ an effective index of 2.92, which asks a comparable thing of the design. It also carries a weak second mode, with an index of 1.10, and a monitor's total power wouldn't show light converted into it. But two modes would beat along the guide, and over six starts, transmission varies by under 2% from 2.1 to 3.3 wavelengths along it, against about 0.5% for $E_z$. So very little light is in the second mode.

| | Plain sharp corner | Six starts |
|---|---:|---:|
| $H_z$, width $0.25\lambda$ | 33% | 94.1 to 97.8% |

These are means over the band on the page's mesh, with its closer monitor. $H_z$ also needs more iterations: 40 per stage, where $E_z$ takes 10.

## An emitter that's electric

A cavity's goal is the LDOS at its emitter. In $E_z$ a line source is an electric current, as an emitter should be. In $H_z$ a line source is a *magnetic* current, and its LDOS says little about an atom or a quantum dot. So we need an in-plane electric dipole. It reads the electric field along its direction $\mathbf{d}$, which $H_z$ gives through a derivative:

$$\mathbf{E}\cdot\mathbf{d} = \frac{i}{k}\,p\,\frac{\partial H_z}{\partial n},$$

with $\mathbf{n}$ a quarter turn from $\mathbf{d}$. So it's a new $e_0$ with the same algebra: $A u = ik\,e_0$, and $F = -2\,\mathrm{Re}(e_0^T u)$. [`solve2d/dipole.rs`](../crates/eigenmenno-core/src/solve2d/dipole.rs) takes the derivative across two points, and scales it so that $F = k/2$ in vacuum, as for a line source. `the_power_a_dipole_supplies_crosses_any_loop` checks it against energy balance, which knows nothing of that scale.

Chapter 17's shortcut, that the source is its own adjoint, still holds. But $e_0$ now holds $p$ at the emitter, so a pixel there moves the source as well as the field. That adds a term,

$$\delta F \supset -4\,\mathrm{Re}\big(u^T \delta e_0\big),$$

which the finite-difference test catches if it's removed: 1.26 against the right 4.99 in one pixel.

## Two things that broke the cavity

**The frequency hold lost its grip.** Chapter 17 places a pole with two solves, assuming that $1/u_0$ is a straight line in $k$. A dipole's response also holds a large, nearly constant background from its near field, since the derivative of a Green's function on itself diverges as the points come together. Then $u_0 \approx B + R/(k - k_p)$, which isn't a straight line until a design is very sharp. The mode swung by $\pm 3\%$ between iterations, and the LDOS fell to under 1. Three solves fit $B$ as well, and that fixes it.

**Grey pixels matter more.** Sharpening went well, then removing the last grey pixels moved the mode, so that finished designs needed resizing by up to 10.4%, where $E_z$ needs 1.2%. That was outside the search for the finished peak. Here the field crosses each pixel's edge, so its value there changes the mode far more. The search now looks further. We also tried holding the frequency through two more sharpness stages. Three starts looked better and nine worse, with a median of 345 against 437, so it was dropped.

## The result

Rods are no benchmark here, since they read about 1. The usual crystal for this polarization is a triangular lattice of air holes, with the middle one filled. We scanned its period and hole size, in the same square, and the best is $a = 0.276\lambda$, $r = 0.4a$. All are at $Q = 1000$, relative to vacuum:

| Structure | LDOS | Mode $Q$ |
|---|---:|---:|
| Crystal of rods | 1 | |
| Crystal of air holes | 371 | |
| Nine optimized starts | 200 to 681, median 437 | 760 to 18,345 |

So the median beats the crystal by only 1.2 times, and four of the nine fall short of it. That's a far smaller margin than chapter 17's in $E_z$, which beats its own crystal by at least 2.9 times, though against a weaker kind of benchmark. On cells half the size, the median holds at 432, and the crystal reads 378. An $H_z$ LDOS and an $E_z$ one measure different emitters, so comparing them directly means little.

## What this doesn't settle

- **The index is still bulk.** These designs use 3.5 and $\varepsilon = 12$. A real 2D model of a silicon slab uses its effective index, about 2.86 for 220 nm ([Nikkhah et al. 2024](20-references.md#nikkhah-2024)), and the third dimension leaks light that none of this sees.
- **Tolerance is unmeasured here.** The results in chapter 17 are for $E_z$, whose field is continuous across an edge. $H_z$'s isn't, so they may not carry over.
- **The dipole's direction.** A design here is optimized and read along one in-plane direction. Its eight-fold symmetry should treat the other direction alike, but we haven't checked that.
- **Few starts,** with spreads wide enough that small differences here shouldn't be trusted. A run is also chaotic: the same start reads 681 in Node and 507 in a browser, with mode $Q$ of 18,345 and 3,420, because their `Math.cosh` differs in the last digit. Each repeats exactly on its own engine. The numbers above are from Node, so the page's own may differ within that spread.

Next: [references](20-references.md).
