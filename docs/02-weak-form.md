# 2. The weak form

The equation from chapter 1 has a second derivative, but at a material interface
that derivative doesn't exist. The weak form fixes this, and it's the central idea
of FEM ([Strang & Fix 2008](20-references.md#strang-fix-2008); [Brenner & Scott 2008](20-references.md#brenner-scott-2008)).

## Testing instead of matching

Asking the equation to hold at every point is too demanding. Instead, we ask for
something that looks weaker: multiply the equation by an arbitrary smooth *test
function* $v(x)$ and require the integral to vanish, for every possible $v$.

$$\int \left[-(p u')' + k_t^2 p u - k_0^2 q u\right] v \, dx = \int f v \, dx$$

If this holds for every $v$, the residual can't be nonzero anywhere it matters,
since we could pick a $v$ concentrated right there. So we've given up nothing for a
smooth solution.

## Integrating by parts

Now the key step. Move one derivative off $u$ and onto $v$:

$$-\int (p u')' v \, dx = \int p u' v' \, dx - \big[p u' v\big]_a^b$$

This gives the **weak form**:

$$\int p\,u'v'\,dx + k_t^2\int p\,u v\,dx - k_0^2\int q\,u v\,dx - \big[p u' v\big]_a^b = \int f v\,dx$$

Three things just happened, and each one is useful.

**Only first derivatives remain.** $u$ now only needs a first derivative that can be
integrated. A function with a kink is fine, so a field with a jump in slope at an
interface fits naturally. The jump is handled without any special code.

**Interface conditions came for free.** We never wrote "$p u'$ is continuous" in the
code. The weak form enforces it automatically, as long as cell boundaries line up
with material boundaries. That's why
[`Mesh::layered_line`](../crates/eigenmenno-core/src/mesh.rs) puts a vertex on every
interface.
<!-- code: Mesh::layered_line "let n = (t / h)" .. "start += t" -->

**The boundary term appeared.** $[p u' v]_a^b$ is the only place the domain's ends
enter the equation, and it's exactly where boundary conditions go. Chapter 4 uses
it for walls and chapter 7 uses it for open ends.

## Galerkin's method

So far $u$ and $v$ are functions, which have infinitely many degrees of freedom. A
computer needs finitely many, so we choose a finite set of basis functions
$\phi_1, \dots, \phi_N$ and write

$$u(x) \approx \sum_{j=1}^N u_j \phi_j(x).$$
<!-- code: H1Line::eval "let (mut v, mut d)" .. "v += u[" -->

The unknowns are now just the numbers $u_j$. For test functions, Galerkin's method
uses the same basis, $v = \phi_i$ ([Gander & Wanner 2012](20-references.md#gander-wanner-2012)). That gives one equation per basis function:

$$\sum_j \left[\underbrace{\int p\,\phi_i'\phi_j'\,dx}_{S_{ij}} + k_t^2\underbrace{\int p\,\phi_i\phi_j\,dx}_{P_{ij}} - k_0^2\underbrace{\int q\,\phi_i\phi_j\,dx}_{Q_{ij}}\right]u_j = \int f\phi_i\,dx$$

In matrix form, this is

$$\left(S + k_t^2 P - k_0^2 Q\right)\mathbf{u} = \mathbf{f}.$$
<!-- code: driven "let mut t = assemble" .. "});" -->

$S$ is called the **stiffness** matrix and $P$ and $Q$ are **mass** matrices, names
inherited from FEM's origins in structural mechanics
([Turner et al. 1956](20-references.md#turner-1956)).

A differential equation has become a matrix equation. All three of $S$, $P$ and $Q$
are of the same two kinds, so one function can build them all. See `assemble` in
[`assemble.rs`](../crates/eigenmenno-core/src/assemble.rs), which builds
$\sum_c (a_c S_c + b_c M_c)$, where the coefficients $a_c$ and $b_c$ are per cell.
For example, `cavity_modes` in [`solve1d.rs`](../crates/eigenmenno-core/src/solve1d.rs)
gets $S$ with $(a_c, b_c) = (p, 0)$ and $Q$ with $(0, q)$.
<!-- code: cavity_modes "let a = s.matrix" .. "let b = s.matrix" -->

## Why Galerkin is a good choice

When the problem is real, symmetric and positive definite, Galerkin's method has a strong property. The
error $u - u_h$ is orthogonal to the whole basis, so the solution is the *best
possible* approximation available in that basis, measured in the problem's natural
energy norm ([Céa 1964](20-references.md#cea-1964)).

There's a matching property for eigenvalues: the approximate eigenvalues are never
below the true ones for this kind of problem, and their error converges twice as
fast as the error in the field
([Strang & Fix 2008](20-references.md#strang-fix-2008); [Boffi 2010](20-references.md#boffi-2010)). We measure both effects in chapter 8.

Next: [elements](03-elements.md), which choose the $\phi_i$.
