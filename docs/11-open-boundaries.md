# 11. Open boundaries and leaky modes

Chapter 7's open end was exact, because in 1D a wave has only one direction to go. In 2D, waves
hit the edge of the domain at every angle at once, and no local boundary condition is reflectionless
for all of them ([Engquist & Majda 1977](20-references.md#engquist-majda-1977); [Johnson 2021](20-references.md#johnson-2021)).
The way out also changes how we can find resonances.

## Perfectly matched layers

A **perfectly matched layer** (PML) surrounds the domain with an absorbing region
([Bérenger 1994](20-references.md#berenger-1994)). The neatest way to understand it is as a change of
coordinates into the complex plane ([Chew & Weedon 1994](20-references.md#chew-weedon-1994)).

An outgoing wave $e^{ikx}$ is an analytic function of $x$. So it can be followed along a path that
leaves the real axis. Replace $x$ by a complex coordinate, with $dx \to s\,dx$ and

$$s = 1 + \frac{i\sigma(x)}{k}.$$
<!-- code: stretched_matrix "let kb" .. "/ kb" -->

Here $k = k_0 n_b$ is the background's wavenumber, for its index $n_b$. Then $e^{ikx}$ becomes $e^{ikx}\,e^{-\int\sigma\,dx}$: the same wave, decaying exponentially. Where
$\sigma = 0$ nothing has changed, and the solution inside the domain is still the exact solution of
the original equation. So the join is free of reflection at every angle and frequency, which no
ordinary absorbing material could manage. Dividing by $k$ makes the decay rate the same at every frequency.
A wall at the far side of the layer reflects what's left, which has been reduced by the round trip.

To get the equations, each derivative along $x$ becomes $\frac{1}{s_x}\partial_x$, and the area
element becomes $s_x s_y\,dA$. In the weak form's $x$ term, that's

$$\int p\,\frac{\partial_x u}{s_x}\,\frac{\partial_x v}{s_x}\,s_x s_y\,dA = \int p\,\frac{s_y}{s_x}\,\partial_x u\,\partial_x v\,dA,$$

and likewise for $y$, while the mass term only gains $s_x s_y$. So it's only a change of the weights:

$$a_x = p\,\frac{s_y}{s_x}, \qquad a_y = p\,\frac{s_x}{s_y}, \qquad b = -k_0^2\,q\,s_x s_y.$$
<!-- code: stretched_matrix "let t = assemble_split" .. "sx * sy)" -->

So a PML looks like a complex, anisotropic *material*
([Sacks et al. 1995](20-references.md#sacks-1995); [Johnson 2021](20-references.md#johnson-2021)).
[`open.rs`](../crates/eigenmenno-core/src/solve2d/open.rs) builds it cell by cell in `open_matrix`.
The stiffness had to split into separate $x$ and $y$ parts, which is `assemble_split`, and nothing else in
assembly changes. The matrix isn't Hermitian any more, but it's still complex **symmetric**, since each
weight only multiplies a symmetric local matrix. Chapter 13 will rely on that.

$\sigma$ grows quadratically into the layer, so that it starts gently, with a strength chosen for a set
loss over a round trip. With $\sigma = \sigma_\text{max}(d/L)^2$ at depth $d$ in a layer $L$ thick, $\int\sigma\,dx = \sigma_\text{max}L/3$,
so a wave crossing it and back at normal incidence keeps $e^{-2\sigma_\text{max}L/3}$ of its amplitude, which is set to $10^{-4}$.
<!-- code: stretched_matrix "`sigma = sigma_max" .. "let sigma_max" -->

On a discrete mesh, the join isn't perfect any more, so the layer can't be
arbitrarily steep or thin. In free space, half a wavelength of layer is a good start. A line source
then matches the exact $\frac{i}{4}H_0^{(1)}(kr)$ to 0.05% with one wavelength of layer at 24 cells
per wavelength, including at the corners.

### Structures run through the layer

The derivation assumes the medium continues unchanged into the layer. If a waveguide stopped at the
layer's edge, a wave would meet the end of the guide, reflect from it, and scatter, and the layer
couldn't help. So anything that reaches the layer must keep going to the outer wall, with the stretch
applied on top of its material.

That's easy to check. A straight slab's power through a monitor reads 0.005735 whether the layer is
0.5, 1 or 1.5 wavelengths thick. Stop the slab where the layer starts, and it reads 0.004354, 24% lower.

Rods running through a layer are subtler: a periodic structure isn't uniform along the path, so
stretching it isn't quite the same structure any more, and a thin layer reflects. Chapter 12 measures how
thick it has to be.

### Plane waves

A plane wave never decays, so it can't simply be a source inside the layer. It enters through the
**scattered-field** form ([Jin 2014](20-references.md#jin-2014)). The incident wave $u_i$ already solves the background problem, so writing the
total as $u = u_i + u_s$ leaves

$$L\,u_s = -(L - L_b)\,u_i,$$
<!-- code: plane_wave "region[c] != 0" .. "*v += i" -->

where $L_b$ is the background's operator. The right side is zero everywhere except inside the shapes,
and $u_s$ is purely outgoing, so the layer absorbs it. The **2D scattering** view does this, and
`dielectric_cylinder_matches_the_mie_series` checks it against the exact series solution for a
cylinder ([Bohren & Huffman 1983](20-references.md#bohren-huffman-1983)), within 1%.

## LDOS, at complex frequency

A line source at $\mathbf{r}_0$ gives $A(k)\,u = ik\,e_0$, where $e_0$ reads the field at that point.
Here $k$ is $k_0$, possibly complex. Dividing by $ik$ gives the Green's function, so
<!-- code: line_response -->

$$F(k) = 2\,\mathrm{Re}(-u_0)$$
<!-- code: driven.rs "Source2d::Line" .. "ldos = Some" -->

is $\pi$ times the **local density of states** (LDOS) per unit $k$. It says how many modes the source
can radiate into, and so how much power it emits: $\mathrm{Re}(-u_0)$ is the work the source does against
its own field, which is chapter 7's minus sign. In vacuum, $u = ik\cdot\frac{i}{4}H_0^{(1)}(kr) = -\frac{k}{4}\big(J_0(kr) + i\,Y_0(kr)\big)$.
At the source, $Y_0$ diverges, but only in the imaginary part, and $J_0(0) = 1$ gives $F = k/2$.
<!-- code: vacuum_ldos_is_half_of_k -->

The ratio to that is what's known as the Purcell enhancement in cavity physics. It's the number a cavity
is really trying to raise.

Now a useful fact. The response $u_0(k)$ is analytic above the real axis, because a causal system can't
respond before it's driven. The stretch keeps that property, since $s = 1 + i\sigma/k$ has no singularity
up there. For such a function, the value at $k + i\Gamma$ is exactly the Lorentzian average of its values
along the real axis, with half-width $\Gamma$.

So one solve at $k_0(1 + i/2Q)$ gives the LDOS averaged over a bandwidth of $1/Q$, however sharp the
resonances inside it ([Liang & Johnson 2013](20-references.md#liang-johnson-2013)). A real frequency would
need many solves to find a sharp peak, or miss it completely. On a fixed mesh this is an exact identity,
which `complex_frequency_is_a_lorentzian_average_of_real_ones` checks.
<!-- code: driven.rs "let k0 = 2.0 * PI" .. "let k = C64::new" -->

Two warnings about that average. In 2D, a Lorentzian's slow tail meets an LDOS that grows with $k$, so
the average carries a small term that grows logarithmically as cells shrink. At $Q = 100$, vacuum reads
about 1% high. [Chapter 17](17-designing-a-cavity.md) starts at lower $Q$, where it's larger. And
the tail reaches the sharp modes far above the mesh's cutoff, which a sum over real frequencies can't
resolve, so the test combines two widths to get a tail that falls off faster.

## Leaky modes, with no eigensolver

The layer makes the eigenproblem non-Hermitian, and its coefficients depend on $k$. A general eigensolver
such as Krylov-Schur could handle the first ([Stewart 2002](20-references.md#stewart-2002)). Driven solves
turn out to be enough, and [`leaky.rs`](../crates/eigenmenno-core/src/solve2d/leaky.rs) needs no new machinery at all.

A leaky mode is a pole of that same response, at a complex $k_p$ below the real axis. Near it,

$$u_0 \approx \frac{R}{k - k_p},$$

so $1/u_0$ is nearly a straight line in $k$, crossing zero at $k_p$. A **secant iteration** finds that zero:
evaluate $1/u_0$ at two points, follow the straight line through them, and repeat. It takes 6 or 7 driven
solves. The mode's sharpness and size then follow ([Joannopoulos et al. 2008](20-references.md#joannopoulos-2008)):
<!-- code: leaky.rs "let mut prev" .. "g = response.inv()" -->

$$Q = \frac{\mathrm{Re}\,k_p}{2\,|\mathrm{Im}\,k_p|}, \qquad
V = \frac{\int q\,|u|^2\,dA}{\max\,q\,|u|^2}.$$
<!-- code: mode_volume "let (mut energy" .. "energy / peak" -->

The mode varies in time as $e^{-ikct}$, so its energy decays as $e^{2\,\mathrm{Im}\,k\,ct}$, and $Q$ counts the
radians of oscillation per e-fold of energy. $V$ is the area the mode would fill at its peak strength. In an
empty box with walls, $\sin^2\sin^2$ averages to 1/4 and peaks at 1, so $V$ is a quarter of the area.
<!-- code: mode_area_of_a_box_is_a_quarter_of_it -->

The solve closest to the pole is dominated by the mode, so it provides the field for $V$. The integral
stops at the layer, whose artificial field isn't part of the mode's energy. Because nothing varies along
$z$, $V$ is an area. It only means the usual electric mode volume for $E_z$: for $H_z$ it measures magnetic
energy.
<!-- code: leaky.rs "let u = normalized" .. "volume: Some" -->

The search needs a source where the mode is strong, and a target within about a quarter of the answer. It
refuses a pole further away, because in an empty domain it otherwise wanders up to a pole of the mesh
and reports a plausible-looking $Q$.
<!-- code: leaky.rs "A search can wander" .. "break;" -->

### Example: one rod missing

Remove one rod from the lattice in chapter 9, and try it on the **2D resonances** view. One mode moves
into the gap. Chapter 9's plane-wave expansion puts it at $a/\lambda = 0.3746$ in a periodic 5 × 5 supercell, and this solver
converges to that from above. With open walls, more crystal blocks more leakage, while the mode's shape
stays put:

| Rods | $a/\lambda$ | $Q$ | $V / a^2$ |
|---:|---:|---:|---:|
| 3 × 3 | 0.3786 | 20 | 0.64 |
| 5 × 5 | 0.3784 | 198 | 0.64 |
| 7 × 7 | 0.3781 | 1,650 | 0.64 |
| 9 × 9 | 0.3781 | 12,700 | 0.65 |
| 13 × 13 | 0.3781 | 683,000 | 0.65 |

These use 12 cells per period, and finer meshes move $Q$ by about 4%. $Q$ is checked against the
half-height of the real-frequency LDOS, not a published value, since we couldn't find one for this
structure. So adding rods changes $Q/V$ only through $Q$. Chapter 17 asks what a free choice of shape
can do in the same footprint.

Next: [waveguides and bends](12-waveguides.md).
