# 7. Driven problems

In a driven problem, we know the frequency and the source and solve for the field:
one linear system rather than an eigenproblem. The interesting part is how to model
a domain that's open to infinity, where waves leave and never come back.

## Open ends

Beyond the structure there's a uniform medium with only an outgoing wave. To the
right of the last layer, $u = c\,e^{ik_x x}$ with
$k_x = \sqrt{k_0^2 \varepsilon\mu - k_t^2}$, so at the end

$$u'(b) = ik_x\,u(b).$$
<!-- code: driven "`p kx` at an end" .. "this entry." -->

Now substitute that into the boundary term from the weak form:

$$-\big[p u' v\big]_{x=b} = -i\,p\,k_x\,u(b)\,v(b).$$

The term is proportional to $u(b)$, so it adds one entry to the matrix, at the
position of the last dof. Abbreviating $Y = p k_x$, that entry is $-iY$.
<!-- code: driven "let last = s.space" .. "t.push" -->

In 1D this is **exact**: the condition is perfectly reflectionless at every frequency
and angle ([Johnson 2021](20-references.md#johnson-2021)). Higher dimensions aren't so lucky, which is why chapter 11 needs perfectly
matched layers.

The square root has two branches, and physics picks one. The wave must carry energy
away or decay, so $\operatorname{Im} k_x \ge 0$. If the sideways wavevector $k_t$ is
large enough, $k_x$ becomes purely imaginary and the wave is evanescent, which is
total internal reflection ([Born & Wolf 1999](20-references.md#born-wolf-1999)).
<!-- code: driven "let admittance" .. "if w.im < 0.0" -->

## The incident wave

On the left there's both an incoming wave of amplitude 1 and a reflected wave:

$$u = e^{ik_x x} + r\,e^{-ik_x x}.$$

Put the left end at $x = 0$, as the mesh does. Then $u(0) = 1 + r$ and
$u'(0) = ik_x(1 - r)$. Eliminating the unknown $r$ gives

$$u'(0) = ik_x\big(2 - u(0)\big).$$

The lower limit of the boundary term enters with a plus sign:

$$+p\,u'(0)\,v(0) = -iY_a\,u(0)\,v(0) + 2iY_a\,v(0).$$

The first part is the matrix entry, like before. The second is known, so it moves to
the right-hand side as $-2iY_a$:

$$\text{matrix: } -iY_a \text{ at the first dof}, \qquad \text{right side: } -2iY_a.$$
<!-- code: driven "Source::PlaneWave =>" -->

The full system is

$$\left(S + k_t^2 P - k_0^2 Q - iY_a\mathbf{e}_a\mathbf{e}_a^T - iY_b\mathbf{e}_b\mathbf{e}_b^T\right)\mathbf{u} = -2iY_a\mathbf{e}_a.$$
<!-- code: driven "let mut t = assemble" .. "});" -->

This matrix is complex and not Hermitian, so only the general LU of chapter 5 will
do. It's also why lossy materials work here but not in the eigensolver of chapter 6:
a complex $\varepsilon$ just makes the entries complex.
<!-- code: driven "BandedLu::factor" .. "lu.solve" -->

## Reflection and transmission

Once we have $\mathbf{u}$, the answers read straight off the two end values:

$$r = u(0) - 1, \qquad t = u(b).$$
<!-- code: driven "let r = (u[0]" .. "let tr" -->

Power involves the wave's admittance, since a wave in glass carries more power than
one of equal amplitude in air:

$$R = |r|^2, \qquad T = |t|^2\,\frac{\operatorname{Re} Y_b}{\operatorname{Re} Y_a}, \qquad A = 1 - R - T.$$
<!-- code: driven "let t_scale" .. "1.0 - r - tr" -->

For a lossless structure, $A$ should come out as zero, which makes it a good built-in
check. In the app it shows as 0.0%.

## Oblique incidence

A wave arriving at angle $\theta$ in a medium of index $n_1$ has $k_t = k_0 n_1\sin\theta$,
and $k_t$ is the same in every layer. That's just **Snell's law**, arising from the
translation symmetry along $z$ ([Born & Wolf 1999](20-references.md#born-wolf-1999)). No separate physics is needed: the solver adds
$k_t^2 P$ and the correct $k_x$ at the ports.
<!-- code: driven "let kt = if plane" .. "let w = (k2" -->

## A current source

A sheet of current at $x_0$ adds a delta function to the equation. With the
$e^{-i\omega t}$ convention, the source term for a unit sheet is

$$-(u')' - k_0^2 u = ik_0\,\delta(x - x_0).$$
<!-- code: driven "Source::Line" -->

The weak form handles the delta function gracefully, since
$\int \delta(x - x_0)\,v\,dx = v(x_0)$. So the right-hand side is just the basis
functions evaluated at $x_0$. See `add_point_load` in
[`solve1d.rs`](../crates/eigenmenno-core/src/solve1d.rs).
<!-- code: add_point_load "let xi" .. "rhs[r] +=" -->

The exact solution in free space is $u = -\tfrac12 e^{ik_0|x - x_0|}$. The minus sign
matters: the field opposes the current at the sheet, which is exactly what lets the
source do work on the field ([Jackson 1999](20-references.md#jackson-1999)).

You can check the sign. Integrate the equation across the sheet. The $k_0^2u$ term
drops out, since $u$ is continuous, so $u'$ jumps by $-ik_0$. The solution
$u = -\tfrac12 e^{ik_0|x - x_0|}$ has slope $+\tfrac12 ik_0$ just left of the sheet and
$-\tfrac12 ik_0$ just right of it, a jump of $-ik_0$. With a plus sign in front, the
jump would be $+ik_0$. An earlier version of this project's test had that plus sign,
and the solver was right.
<!-- code: line_source_radiates_the_free_space_green_function "-u'' - k^2 u" .. "let exact" -->

Next: [is the answer right?](08-verification.md)
