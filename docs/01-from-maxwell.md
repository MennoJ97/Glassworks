# 1. From Maxwell's equations to one equation

FEM solves differential equations, so first we need the equation. Starting from
Maxwell's equations, we'll use the symmetry of a layered structure to reduce six
field components to one unknown function of one variable.

## Maxwell's curl equations

We assume every field oscillates at one angular frequency, as $e^{-i\omega t}$.
This choice of sign matters later: with it, a wave moving toward $+x$ looks like
$e^{+ikx}$, and a lossy material has a *positive* imaginary permittivity
([Jackson 1999](20-references.md#jackson-1999)).

With no sources, the curl equations become

$$\nabla \times \mathbf{E} = i\omega\mu_0\mu\,\mathbf{H}, \qquad \nabla \times \mathbf{H} = -i\omega\varepsilon_0\varepsilon\,\mathbf{E}.$$

Here $\varepsilon(x)$ and $\mu(x)$ are the relative permittivity and permeability of
the material. A material's refractive index $n$ and absorption $k$ give
$\varepsilon = (n + ik)^2$ when $\mu = 1$.

We use two shorthands. The vacuum wavenumber is $k_0 = \omega/c = 2\pi/\lambda$, and
$\eta_0 = \sqrt{\mu_0/\varepsilon_0} \approx 377\,\Omega$ converts $\mathbf{H}$ into
the same units as $\mathbf{E}$. Plots in the app show $\eta_0 \mathbf{H}$, so both
fields share one axis.

## Using the symmetry

In a stack of flat layers, the material only changes along $x$. Nothing changes
along $y$, so $\partial_y = 0$. Along $z$, the structure is uniform, so the field
can only pick up a phase: it varies as $e^{ik_t z}$ for some constant $k_t$.

Writing out the curl equations with those two facts, the six components split into
two independent groups ([Yeh 2005](20-references.md#yeh-2005)):

| Polarization | Components | Solved for |
|---|---|---|
| TE | $E_y, H_x, H_z$ | $u = E_y$ |
| TM | $H_y, E_x, E_z$ | $u = \eta_0 H_y$ |

Each group needs only one unknown, because the other two components follow from it.
For TE:

$$\eta_0 H_x = -\frac{k_t}{k_0\mu}E_y, \qquad \eta_0 H_z = -\frac{i}{k_0\mu}\frac{dE_y}{dx}.$$
<!-- code: Setup::push_field "let (second, third)" .. "-I * d * scale" -->

For TM:

$$E_x = \frac{k_t}{k_0\varepsilon}\,\eta_0 H_y, \qquad E_z = \frac{i}{k_0\varepsilon}\frac{d(\eta_0 H_y)}{dx}.$$
<!-- code: Setup::push_field "} else {" .. "I * d * scale" -->

These are exactly the lines in `Setup::push_field` in
[`solve1d.rs`](../crates/glassworks-core/src/solve1d.rs).

## The one equation

Eliminating the other components gives one equation. For TE, use
$\omega\mu_0 = k_0\eta_0$ and $\omega\varepsilon_0 = k_0/\eta_0$. Since $\partial_z$ becomes
$ik_t$, the $y$ part of the second curl equation, times $\eta_0$, is

$$ik_t\,\eta_0 H_x - \frac{d}{dx}(\eta_0 H_z) = -ik_0\varepsilon E_y.$$

Substitute $\eta_0 H_x$ and $\eta_0 H_z$ from above and multiply by $ik_0$:

$$-\frac{d}{dx}\left(\frac{1}{\mu}\frac{dE_y}{dx}\right) + \frac{k_t^2}{\mu}E_y - k_0^2\varepsilon E_y = 0.$$

TM works the same way, from the $y$ part of the first curl equation, with
$\varepsilon$ and $\mu$ swapped. So both polarizations share one form
([Jin 2014](20-references.md#jin-2014)):

$$-\frac{d}{dx}\left(p\,\frac{du}{dx}\right) + k_t^2\,p\,u - k_0^2\,q\,u = f$$

| | $u$ | $p$ | $q$ |
|---|---|---|---|
| TE | $E_y$ | $1/\mu$ | $\varepsilon$ |
| TM | $\eta_0 H_y$ | $1/\varepsilon$ | $\mu$ |
<!-- code: layer_coef "match pol" .. "Polarization::Tm" -->

The right side $f$ is a source, which is zero unless we add a current. A current
density $\mathbf{J}$ adds $+\mathbf{J}$ to the right of the second curl equation. In TE,
a current sheet $J_y = J\,\delta(x - x_0)$ goes through the same steps and gives
$f = ik_0\eta_0 J\,\delta(x - x_0)$. Chapter 7 uses a unit sheet, $\eta_0 J = 1$.
<!-- code: driven "Source::Line" -->

Everything Glassworks calculates is this one equation with different questions asked
of it:

| Analysis | Known | Unknown |
|---|---|---|
| Resonances | $k_t = 0$ | $k_0$ and $u$ (an eigenproblem) |
| Band diagram | Bloch phase across the cell | $k_0$ and $u$ (an eigenproblem) |
| Guided modes | $k_0$ | $\beta = k_t$ and $u$ (an eigenproblem) |
| Spectrum, single wavelength | $k_0$, $k_t$, $f$ | $u$ (a linear system) |

In `solve1d.rs`, `cavity_modes`, `bands`, `waveguide` and `driven` each answer one
row of this table.

One more thing to notice is the $p$ inside the derivative. It can jump at an
interface between materials. Where $p$ jumps, $u'$ must jump too, because $p u'$ is
continuous. From the equations above, $p u' = ik_0\,\eta_0 H_z$ in TE and
$p u' = -ik_0 E_z$ in TM. Both are fields tangential to the interface, and tangential
fields are continuous across an interface that carries no current sheet. So $u'$ isn't even defined at an interface, and
$u''$ certainly isn't. That's a real problem for any method that needs second
derivatives, and it's the reason for the next chapter.

Next: [the weak form](02-weak-form.md).
