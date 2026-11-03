# 15. Diffraction gratings

*Going further. This builds on chapter 9's Bloch ties, chapter 11's absorbing layers and chapter 13's adjoint.*

A grating repeats forever along $x$, and it's open above and below. So it needs the two boundary conditions we already have at once: Bloch ties on the sides, and an absorbing layer on the top and bottom. Try it on the **2D grating** view, under **Driven response**. The code is in [`solve2d/grating.rs`](../crates/eigenmenno-core/src/solve2d/grating.rs).

## Diffraction orders

A plane wave at angle $\theta$ has $k_x = k\sin\theta$ along the grating. Since the structure repeats with period $L$, what comes back can only add multiples of $2\pi/L$:

$$k_{x,m} = k_x + \frac{2\pi m}{L}, \qquad k_{y,m} = \sqrt{k^2 - k_{x,m}^2}.$$
<!-- code: grating.rs "let orders" .. ".collect()" -->

These are the **diffraction orders**. Where $|k_{x,m}| < k$ an order propagates, and otherwise it's evanescent and carries no power. What we want from the solver is the share of power in each propagating order.

## Launching and reading

The plane wave comes from a sheet of current with the wave's own phase along $x$, as in chapter 12's port. This has one advantage over chapter 11's scattered-field form: a substrate can run out through the layer, because nothing assumes a uniform background.
<!-- code: grating.rs "The sheet" .. "lu.solve(&mut u)" -->

Each order's amplitude comes from a line across the period: the average of $u\,e^{-ik_{x,m}x}$ picks out that order, and all the others average to zero. Below the structure, that's the answer. Above it, a line holds both the reflected wave and the incident one, so four lines fit both at once, by least squares:
<!-- code: add_order -->

$$u_m(y) = a_m^{\uparrow} e^{ik_{y,m}y} + a_m^{\downarrow} e^{-ik_{y,m}y}.$$
<!-- code: split_order "let n" .. "(up, down)" -->

These use field values rather than derivatives, which chapter 12 found are least accurate at a cell's edge. The share of power in an order is then

$$\eta_m = \frac{p\,k_{y,m}\,|a_m|^2}{p_0\,k_{y,0}\,|a_0^{\downarrow}|^2}.$$
<!-- code: grating.rs "let (_, incident)" .. "let power_in" -->
<!-- code: grating.rs "let ratio" .. "let share" -->

Measuring the incident amplitude this way means nothing depends on how exactly the mesh launches the wave.

## What an absorbing layer can't do

The first version disagreed with its reference by far more than the mesh could explain, and refining the mesh changed nothing. There were two separate causes, and neither is a mesh error.

A layer absorbs at a rate of $k_y/k$ ([Johnson 2021](20-references.md#johnson-2021)). So an order that's nearly grazing barely notices it: at $k_y/k = 0.13$, a layer tuned to lose $10^{-4}$ on a round trip at normal incidence loses only a factor of 0.3. A thicker layer doesn't help, since the absorption is set for its whole thickness. Here we know every order in advance, so the layer can be set for the slowest one.

The second is an order that doesn't propagate at all. Chapter 11's stretch $s = 1 + i\sigma/k$ only rotates an evanescent wave's phase, without adding any decay. So an order just past its cutoff decays slowly, passes through the layer untouched, and returns from the wall. That needs a real part in the stretch ([Johnson 2021](20-references.md#johnson-2021)), which does speed up its decay, and which leaves a propagating wave alone apart from shortening its wavelength:

$$s = 1 + \left(\frac{k}{\kappa_{\min}} + i\,\frac{k}{k_{y,\min}}\right)\frac{\sigma}{k}.$$
<!-- code: stretched_matrix "let stretch" .. "I * strength" -->

Here $k_{y,\min}$ belongs to the most grazing propagating order and $\kappa_{\min}$ to the slowest evanescent one. Both factors are capped at 10 (`MAX_STRENGTH`), since a steeper layer stops being resolved. The first change alone took the most grazing reflected order from 0.00122 to the correct 0.00179, and the second took the zeroth transmitted order from 0.45894 to the correct 0.45646, on the same mesh.
<!-- code: grating.rs "let background" .. "MAX_STRENGTH));" -->

## Checking it

The reference is the Fourier modal method ([Moharam et al. 1995](20-references.md#moharam-1995)), with Li's rule for $H_z$ ([Li 1996](20-references.md#li-1996)), in `plans/probe-grating/rcwa.py`. It expands fields in plane waves, with no mesh and no absorbing layer, so it shares nothing with this solver. For a ridge of index 2.5 on glass, with $L = 1.5\lambda$ at $20°$, on cells of $\lambda/80$ and one wavelength of layer:

| Order | $E_z$ here | $E_z$ reference | $H_z$ here | $H_z$ reference |
|---:|---:|---:|---:|---:|
| Reflected $-2$ | 0.00254 | 0.00254 | 0.01012 | 0.01013 |
| Reflected $-1$ | 0.05361 | 0.05361 | 0.04710 | 0.04708 |
| Reflected $0$ | 0.13877 | 0.13877 | 0.12396 | 0.12398 |
| Transmitted $-2$ | 0.28440 | 0.28439 | 0.05527 | 0.05528 |
| Transmitted $-1$ | 0.21809 | 0.21810 | 0.16843 | 0.16844 |
| Transmitted $0$ | 0.11827 | 0.11828 | 0.26553 | 0.26549 |
| Transmitted $1$ | 0.18431 | 0.18431 | 0.32958 | 0.32960 |

[`tests/gratings.rs`](../crates/eigenmenno-core/tests/gratings.rs) also checks a flat interface against Fresnel's formula, and that a symmetric grating mirrors its orders.

A real limit remains. Where an order sits exactly on its cutoff, a **Rayleigh anomaly**, the shares change abruptly with wavelength, and the slowest rate goes to zero, past what the cap can cover. For a period of $1.5\lambda$, the zeroth transmitted order is about 1% off at that wavelength, under 0.2% within 0.1% of it, and about $10^{-4}$ off once 7% away.

## Why the adjoint needs a transposed solve

Chapter 13 relied on $A^T = A$. Bloch ties break this. With $Z$ carrying the phases, $A = Z^H A_\text{full} Z$, where $A_\text{full}$ is still complex symmetric, so

$$A^T = Z^T A_\text{full} \overline{Z},$$
<!-- code: grating.rs "lu.solve_transposed" .. "let full" -->

which is the same operator at the opposite phase. So the adjoint really needs $A^T w = g$. Happily, an LU factorization already holds the answer. A solve is a sequence of simple steps, so $A^{-T}$ is each step's transpose, in the reverse order, with no new factorization. [`sparse/transposed.rs`](../crates/eigenmenno-core/src/sparse/transposed.rs) mirrors the ordinary solve line by line.
<!-- code: SparseLu::solve_transposed_once -->

`a_grating_gradient_matches_finite_differences` checks the result at an angle, where the two would differ, and for an order on each side. The page doesn't optimize gratings yet, but the solver returns everything that would need.

Next: [making it faster](16-making-it-faster.md).
