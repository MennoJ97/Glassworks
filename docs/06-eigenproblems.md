# 6. Eigenproblems

A cavity's resonances are the values of $k_0$ where the equation has a nonzero
solution with no source. In matrix form:

$$S\mathbf{u} = k_0^2\,Q\,\mathbf{u}$$
<!-- code: cavity_modes "let a" .. "let b" -->

This is a **generalized eigenproblem**, $A\mathbf{x} = \lambda B\mathbf{x}$. Here $A$
is Hermitian and $B$ is positive definite, which guarantees real eigenvalues
([Parlett 1998](20-references.md#parlett-1998)).

Guided modes have the same form. Rearranging the operator from chapter 2 with
$k_t = \beta$:

$$(k_0^2 Q - S)\,\mathbf{u} = \beta^2\,P\,\mathbf{u}.$$
<!-- code: waveguide "let a" .. "let b" -->

The matrix may have thousands of eigenvalues and we want only a few. So we need an
iterative method.

## The power method, and its problem

Repeatedly multiplying a vector by a matrix makes it line up with the eigenvector of
largest $|\lambda|$ ([Golub & Van Loan 2013](20-references.md#golub-van-loan-2013)). It's simple, but it finds the *largest* eigenvalue, and we want
the smallest ones, or the ones near a particular value. The largest eigenvalues of
an FEM matrix are also the mesh's least accurate modes.

## Shift and invert

The trick is to change the operator. Pick a target $\sigma$ near the modes you want,
and consider

$$(A - \sigma B)^{-1} B\,\mathbf{x} = \theta\,\mathbf{x}, \qquad \theta = \frac{1}{\lambda - \sigma}.$$
<!-- code: eigs_factored "let mut pairs" -->

It has the same eigenvectors, but eigenvalues close to $\sigma$ now have the *largest*
$|\theta|$. The modes you want have become the easiest to find
([Ericsson & Ruhe 1980](20-references.md#ericsson-ruhe-1980)).

Applying the operator means solving a linear system, but we only factor
$A - \sigma B$ once and reuse the factors each iteration. Chapter 5's fast solver is
what makes this cheap.
<!-- code: eigs_near "sub_scaled" -->
<!-- code: lanczos "u.copy_from_slice" .. "lu.solve" -->

Each analysis chooses its own shift:
- **Resonances:** a small negative $\sigma = -(\pi/L)^2$, for a structure of length
  $L$, unless you ask for modes near a given $k_0$. A zero shift would make $S$
  singular whenever a static mode exists, as in the TM cavity of chapter 4. A negative
  shift keeps $S - \sigma Q$ positive definite, and the lowest modes are still the
  closest to it.
- **Guided modes:** just above $k_0^2 n_\text{max}^2$, the largest $\beta^2$ physically
  possible, since the most tightly confined modes are the ones we want.
<!-- code: cavity_modes "A zero shift" .. "let sigma" -->
<!-- code: waveguide "let n2: Vec" .. "n2_max * 1.001" -->

## Lanczos

The power method throws away every vector but the last. **Lanczos** keeps them all,
building an orthonormal basis $V_m$ of the *Krylov space*
$\{\mathbf{v}, C\mathbf{v}, C^2\mathbf{v}, \dots\}$ of the shifted operator
$C = (A - \sigma B)^{-1}B$, and solving the eigenproblem restricted to that space ([Lanczos 1950](20-references.md#lanczos-1950)).

For a Hermitian operator, a remarkable simplification happens. The restricted matrix
is **tridiagonal**, so each new vector needs only a three-term recurrence:

$$C\mathbf{v}_j = \beta_{j-1}\mathbf{v}_{j-1} + \alpha_j\mathbf{v}_j + \beta_j\mathbf{v}_{j+1}$$

These $\alpha_j$ and $\beta_j$ are Lanczos coefficients, not the propagation constant.
Call the restricted matrix $T_m$. It has $\alpha_j$ on its diagonal and $\beta_j$
beside it. The code runs one step like this:

```text
u = C v_j
alpha_j = <v_j, u>_B
subtract from u its projection on every v_i, twice
beta_j = |u|_B, then v_{j+1} = u / beta_j
```

The full projection includes the two recurrence terms, $\alpha_j\mathbf{v}_j$ and
$\beta_{j-1}\mathbf{v}_{j-1}$, so what's left is $\beta_j\mathbf{v}_{j+1}$. In the
code, $\beta_j$ is the norm computed at the top of the next pass.
<!-- code: lanczos "b.mul(&u, &mut bu)" .. "found.iter().for_each" -->

The eigenvalues of $T_m$, called *Ritz values*, are very good approximations to the outermost
eigenvalues of $C$ after just a few steps
([Parlett 1998](20-references.md#parlett-1998); [Saad 2011](20-references.md#saad-2011)). `tridiagonal_eigen` in
[`eigen.rs`](../crates/eigenmenno-core/src/eigen.rs) finds them with the QL algorithm
([Bowdler et al. 1968](20-references.md#bowdler-1968)). It's the standard method, like LAPACK's
`dsteqr`, so you can treat it as a black box.
<!-- code: lanczos "let (theta, z)" .. "order.sort_by" -->

Three details in the code are worth pointing out.

**The inner product uses $B$.** The operator is only Hermitian in the $B$ inner
product, $\langle \mathbf{x}, \mathbf{y}\rangle_B = \mathbf{x}^H B\,\mathbf{y}$, so the
basis is $B$-orthonormal ([Ericsson & Ruhe 1980](20-references.md#ericsson-ruhe-1980)). That's why the code stores both `v` and `bv`.
<!-- code: lanczos "b.mul(&u, &mut bu)" .. "let norm" -->
<!-- code: lanczos "v.push(u.iter()" .. "bv.push" -->

To see why $C$ is Hermitian in this product, write it out:

$$\langle \mathbf{x}, C\mathbf{y}\rangle_B = \mathbf{x}^H B(A - \sigma B)^{-1}B\,\mathbf{y}.$$

The matrix $B(A - \sigma B)^{-1}B$ is Hermitian, because $B$ and $(A - \sigma B)^{-1}$ are,
so $C$ is self-adjoint in the $B$ inner product. That's also why $\alpha_j$ and $\beta_j$
are real, and the code keeps only the real part, `.re`.
<!-- code: lanczos "alpha.push" -->

**Full reorthogonalization.** In exact arithmetic, the three-term recurrence keeps the
basis orthogonal. In floating point, it doesn't, and copies of converged eigenvectors
creep back in ([Paige 1972](20-references.md#paige-1972)). So each new vector is orthogonalized again against the whole basis.
This costs memory and $O(m^2 N)$ time in total, but it makes the code short and
reliable.
<!-- code: lanczos "Full reorthogonalization" .. "found.iter().for_each" -->

**The stopping test is nearly free.** The residual of a Ritz pair equals the coupling
$\beta_m$ to the next basis vector times the last component of its eigenvector in the
tridiagonal problem. So we can tell what has converged without computing a single
actual eigenvector ([Parlett 1998](20-references.md#parlett-1998)).
<!-- code: lanczos "A Ritz pair's residual" .. "let converged" -->

## Degenerate modes

A Krylov space started from one vector contains only one direction per eigenspace.
So when two modes share an eigenvalue, Lanczos can only find one of them. This
actually happens: a uniform periodic structure has pairs of modes at the same frequency, such as
waves going left and right at $k = 0$.

`eigs_near` handles this by **locking** ([Bai et al. 2000](20-references.md#bai-2000)). Its loop, in `eigs_factored`, runs Lanczos, keeps the converged modes,
then runs it again from a fresh random vector that's kept orthogonal to them. The
second copy of each degenerate mode is then visible. It stops when a run finds
nothing closer than the current $k$-th best mode. The test
`bloch_modes_of_uniform_medium_include_degenerate_pairs` checks this.
<!-- code: eigs_factored "while found.len()" .. "break" -->
<!-- code: lanczos "found.iter().for_each" -->

## The known weak spot

Lanczos converges quickly when the wanted eigenvalues are well separated. A guided
mode close to cutoff in a wide cladding is the bad case: it sits among a dense cluster
of *box modes*, which are artefacts of closing the domain with walls. Those need a
large basis. In the browser on an Apple M4 Pro, a waveguide with 20 µm of cladding took 187 ms for
only 1,431 unknowns, against 8 ms for the resonances of the same layers on the same mesh.
<!-- code: MAX_BASIS -->

The standard fix is **thick restarting**, as in Krylov-Schur, which bounds the basis
size ([Wu & Simon 2000](20-references.md#wu-simon-2000); [Stewart 2002](20-references.md#stewart-2002)). We left it out to keep the code readable, and `faer` provides it if we ever
need it.

Next: [driven problems](07-driven.md).
