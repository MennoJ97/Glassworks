# 7. Driven problems

In a driven problem, we know the frequency and the source and solve for the field:
one linear system rather than an eigenproblem. The interesting part is how to model
a domain that's open to infinity, where waves leave and never come back.

## Open ends

Beyond the structure there's a uniform medium with only an outgoing wave. To the
right of the last layer, $u = A e^{ik_x x}$ with
$k_x = \sqrt{k_0^2 \varepsilon\mu - k_t^2}$, so at the end

$$u'(b) = ik_x\,u(b).$$

Now substitute that into the boundary term from the weak form:

$$-\big[p u' v\big]_{x=b} = -i\,p\,k_x\,u(b)\,v(b).$$

The term is proportional to $u(b)$, so it adds one entry to the matrix, at the
position of the last dof. Abbreviating $Y = p k_x$, that entry is $-iY$.

In 1D this is **exact**: the condition is perfectly reflectionless at every frequency
and angle ([Johnson 2021](20-references.md#johnson-2021)). Higher dimensions aren't so lucky, which is why chapter 11 needs perfectly
matched layers.

The square root has two branches, and physics picks one. The wave must carry energy
away or decay, so $\operatorname{Im} k_x \ge 0$. If the sideways wavevector $k_t$ is
large enough, $k_x$ becomes purely imaginary and the wave is evanescent, which is
total internal reflection ([Born & Wolf 1999](20-references.md#born-wolf-1999)).

## The incident wave

On the left there's both an incoming wave of amplitude 1 and a reflected wave:

$$u = e^{ik_x x} + r\,e^{-ik_x x}.$$

Eliminating the unknown $r$ between $u(a)$ and $u'(a)$ gives

$$u'(a) = ik_x\big(2 - u(a)\big).$$

Substituting into the boundary term, the $u(a)$ part goes in the matrix like before,
and the known part moves to the right-hand side:

$$\text{matrix: } -iY_a \text{ at the first dof}, \qquad \text{right side: } -2iY_a.$$

The full system is

$$\left(S + k_t^2 P - k_0^2 Q - iY_a\mathbf{e}_a\mathbf{e}_a^T - iY_b\mathbf{e}_b\mathbf{e}_b^T\right)\mathbf{u} = -2iY_a\mathbf{e}_a.$$

This matrix is complex and not Hermitian, so only the general LU of chapter 5 will
do. It's also why lossy materials work here but not in the eigensolver of chapter 6:
a complex $\varepsilon$ just makes the entries complex.

## Reflection and transmission

Once we have $\mathbf{u}$, the answers read straight off the two end values:

$$r = u(a) - 1, \qquad t = u(b).$$

Power involves the wave's admittance, since a wave in glass carries more power than
one of equal amplitude in air:

$$R = |r|^2, \qquad T = |t|^2\,\frac{\operatorname{Re} Y_b}{\operatorname{Re} Y_a}, \qquad A = 1 - R - T.$$

For a lossless structure, $A$ should come out as zero, which makes it a good built-in
check. In the app it shows as 0.0%.

## Oblique incidence

A wave arriving at angle $\theta$ in a medium of index $n_1$ has $k_t = k_0 n_1\sin\theta$,
and $k_t$ is the same in every layer. That's just **Snell's law**, arising from the
translation symmetry along $z$ ([Born & Wolf 1999](20-references.md#born-wolf-1999)). No separate physics is needed: the solver adds
$k_t^2 P$ and the correct $k_x$ at the ports.

## A current source

A sheet of current at $x_0$ adds a delta function to the equation. With the
$e^{-i\omega t}$ convention, the source term for a unit sheet is

$$-(u')' - k_0^2 u = ik_0\,\delta(x - x_0).$$

The weak form handles the delta function gracefully, since
$\int \delta(x - x_0)\,v\,dx = v(x_0)$. So the right-hand side is just the basis
functions evaluated at $x_0$. See `add_point_load` in
[`solve1d.rs`](../crates/eigenmenno-core/src/solve1d.rs).

The exact solution in free space is $u = -\tfrac12 e^{ik_0|x - x_0|}$. Note the minus
sign: the field opposes the current at the sheet, which is exactly what lets the
source do work on the field. An earlier version of this project's test had the sign
wrong, and the solver was right. Chapter 8 tells that story.

Next: [is the answer right?](08-verification.md)
