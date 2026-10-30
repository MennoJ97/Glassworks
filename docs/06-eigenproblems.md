# 6. Eigenproblems

A cavity's resonances are the values of $k_0$ where the equation has a nonzero
solution with no source. In matrix form:

$$S\mathbf{u} = k_0^2\,Q\,\mathbf{u}$$

This is a **generalized eigenproblem**, $A\mathbf{x} = \lambda B\mathbf{x}$. Here $A$
is Hermitian and $B$ is positive definite, which guarantees real eigenvalues
([Parlett 1998](10-references.md#parlett-1998)).

Guided modes have the same form. Rearranging the operator from chapter 2 with
$k_t = \beta$:

$$(k_0^2 Q - S)\,\mathbf{u} = \beta^2\,P\,\mathbf{u}.$$

The matrix may have thousands of eigenvalues and we want only a few. So we need an
iterative method.

## The power method, and its problem

Repeatedly multiplying a vector by a matrix makes it line up with the eigenvector of
largest $|\lambda|$ ([Golub & Van Loan 2013](10-references.md#golub-van-loan-2013)). It's simple, but it finds the *largest* eigenvalue, and we want
the smallest ones, or the ones near a particular value. The largest eigenvalues of
an FEM matrix are also the mesh's least accurate modes.

## Shift and invert

The trick is to change the operator. Pick a target $\sigma$ near the modes you want,
and consider

$$(A - \sigma B)^{-1} B\,\mathbf{x} = \theta\,\mathbf{x}, \qquad \theta = \frac{1}{\lambda - \sigma}.$$

It has the same eigenvectors, but eigenvalues close to $\sigma$ now have the *largest*
$|\theta|$. The modes you want have become the easiest to find
([Ericsson & Ruhe 1980](10-references.md#ericsson-ruhe-1980)).

Applying the operator means solving a linear system, but we only factor
$A - \sigma B$ once and reuse the factors each iteration. Chapter 5's fast solver is
what makes this cheap.

Each analysis chooses its own shift:
- **Resonances:** a small negative $\sigma$, which finds the lowest modes and keeps
  the operator definite.
- **Guided modes:** just above $k_0^2 n_\text{max}^2$, the largest $\beta^2$ physically
  possible, since the most tightly confined modes are the ones we want.

## Lanczos

The power method throws away every vector but the last. **Lanczos** keeps them all,
building an orthonormal basis $V_m$ of the *Krylov space*
$\{\mathbf{v}, C\mathbf{v}, C^2\mathbf{v}, \dots\}$ and solving the eigenproblem
restricted to that space ([Lanczos 1950](10-references.md#lanczos-1950)).

For a Hermitian operator, a remarkable simplification happens. The restricted matrix
is **tridiagonal**, so each new vector needs only a three-term recurrence:

$$C\mathbf{v}_j = \beta_{j-1}\mathbf{v}_{j-1} + \alpha_j\mathbf{v}_j + \beta_j\mathbf{v}_{j+1}$$

Its eigenvalues, called *Ritz values*, are very good approximations to the outermost
eigenvalues of $C$ after just a few steps
([Parlett 1998](10-references.md#parlett-1998); [Saad 2011](10-references.md#saad-2011)). `tridiagonal_eigen` in
[`eigen.rs`](../crates/eigenmenno-core/src/eigen.rs) finds them with the QL algorithm
([Bowdler et al. 1968](10-references.md#bowdler-1968)).

Three details in the code are worth pointing out.

**The inner product uses $B$.** The operator is only Hermitian in the $B$ inner
product, $\langle \mathbf{x}, \mathbf{y}\rangle_B = \mathbf{x}^H B\,\mathbf{y}$, so the
basis is $B$-orthonormal ([Ericsson & Ruhe 1980](10-references.md#ericsson-ruhe-1980)). That's why the code stores both `v` and `bv`.

**Full reorthogonalization.** In exact arithmetic, the three-term recurrence keeps the
basis orthogonal. In floating point, it doesn't, and copies of converged eigenvectors
creep back in ([Paige 1972](10-references.md#paige-1972)). So each new vector is orthogonalized again against the whole basis.
This costs memory and $O(m^2 N)$ time in total, but it makes the code short and
reliable.

**The stopping test is nearly free.** The residual of a Ritz pair equals the coupling
to the next basis vector times the last component of its eigenvector in the
tridiagonal problem. So we can tell what has converged without computing a single
actual eigenvector ([Parlett 1998](10-references.md#parlett-1998)).

## Degenerate modes

A Krylov space started from one vector contains only one direction per eigenspace.
So when two modes share an eigenvalue, Lanczos can only find one of them. This
actually happens: a uniform periodic structure has two modes at each frequency, for
waves going left and right.

`eigs_near` handles this by **locking** ([Bai et al. 2000](10-references.md#bai-2000)). It runs Lanczos, keeps the converged modes,
then runs it again from a fresh random vector that's kept orthogonal to them. The
second copy of each degenerate mode is then visible. It stops when a run finds
nothing closer than the current $k$-th best mode. The test
`bloch_modes_of_uniform_medium_include_degenerate_pairs` checks this.

## The known weak spot

Lanczos converges quickly when the wanted eigenvalues are well separated. A guided
mode close to cutoff in a wide cladding is the bad case: it sits among a dense cluster
of *box modes*, which are artefacts of closing the domain with walls. Those need a
large basis. A waveguide with 20 µm of cladding took 454 ms for only 2,119 unknowns,
compared with 23 ms for a 2,799-unknown resonance problem.

The standard fix is **thick restarting**, as in Krylov-Schur, which bounds the basis
size ([Wu & Simon 2000](10-references.md#wu-simon-2000); [Stewart 2002](10-references.md#stewart-2002)). We left it out to keep the code readable, and `faer` provides it if we ever
need it.

Next: [driven problems](07-driven.md).
