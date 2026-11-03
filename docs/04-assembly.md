# 4. Assembly and boundary conditions

We can now compute a small matrix for each cell. This chapter combines them into one
global matrix, then applies the boundary conditions.

## Assembly

Each global entry $S_{ij}$ is an integral over the whole domain, which is a sum over
cells. So we loop over cells and add each local matrix into the global positions
given by its dof numbers (dof means *degree of freedom*, one unknown):

```text
for each cell c:
    compute local matrices S_c and M_c
    for each pair of local dofs (i, j):
        A[dof(i), dof(j)] += a_c * S_c[i, j] + b_c * M_c[i, j]
```
<!-- code: assemble_split "for c in 0..space.n_cells()" .. "t.push" -->

Where two cells share a dof, their contributions add. You could see this in the
worked example in chapter 3, where the middle node received $1/h$ from each side.

This loop is the function `assemble` in
[`assemble.rs`](../crates/eigenmenno-core/src/assemble.rs), and it doesn't know about
1D. It only sees the `Space` trait, which supplies dof numbers and local matrices.
That's what makes 2D an addition rather than a rewrite. The loop itself is in
`assemble_split`, which also lets 2D weight the $x$ and $y$ derivatives separately.

## Storing a sparse matrix

The loop produces a list of `(row, column, value)` triplets, including duplicates.
`Csr::from_triplets` deals them into their rows, sorts each short row by column, and
adds the duplicates together. That produces
**compressed sparse row** storage: all the column numbers and values in one long
array, and a pointer to where each row starts ([Davis 2006](20-references.md#davis-2006)). Multiplying by a vector then takes
one pass through the nonzeros, as in `Csr::mul`.
<!-- code: Csr -->
<!-- code: Csr::from_triplets "for (i, j, v) in sorted" .. "last = Some" -->

Take chapter 3's two cells of order 1, but with natural ends, so all 3 dofs stay. The
two cells give 8 triplets, two of them at (1, 1), and merging leaves 7 entries. For
$S$ with $h = 1/2$, row $i$ runs from `ptr[i]` up to `ptr[i + 1]`:

```text
ptr = [0, 2, 5, 7]
col = [0, 1, 0, 1, 2, 1, 2]
val = [2, -2, -2, 4, -2, -2, 2]
```

## Two kinds of boundary condition

Recall the boundary term from the weak form: $-[p u' v]_a^b$.

**Essential (Dirichlet) conditions** fix the value of $u$, for example $u = 0$. You
have to impose these explicitly, by removing that unknown from the system. Removing
the dof also removes its test function. Every remaining $\phi_i$ is zero at the wall,
so the boundary term vanishes there too, whatever $pu'$ is.

**Natural (Neumann) conditions** specify $p u'$. These need no special handling: if
$p u' = 0$ at a boundary, the boundary term vanishes, and leaving it out of the
matrix *is* the boundary condition. Doing nothing imposes it ([Strang & Fix 2008](20-references.md#strang-fix-2008)).

Physically:

| Wall | TE ($u = E_y$) | TM ($u = \eta_0 H_y$) |
|---|---|---|
| Perfect conductor (PEC) | $E_y = 0$: essential | $dH_y/dx = 0$: natural |
| Magnetic wall (PMC) | $dE_y/dx = 0$: natural | $H_y = 0$: essential |
<!-- code: Setup::new "let wall" .. "zero.push(last)" -->

A perfect conductor forces the tangential electric field to zero
([Jackson 1999](20-references.md#jackson-1999)). In TE that's the
unknown itself. In TM the same wall forces $E_z \propto du/dx$ to zero, which is a
condition on the derivative, so it's natural. The test
`tm_polarization_swaps_pec_and_pmc` checks this: a TM cavity with PEC walls has a
static mode at $k_0 = 0$, as the theory says.

## One mechanism for all constraints

`DofMap` in [`assemble.rs`](../crates/eigenmenno-core/src/assemble.rs) gives each
global dof a reduced index and a factor, or marks it as fixed to zero:

$$u_\text{global} = \text{factor} \times u_\text{reduced}$$
<!-- code: DofMap::expand "self.map.iter()" -->

A Dirichlet wall sets a dof to `None`. That's the whole mechanism.
<!-- code: DofMap::new "let mut fixed" .. "d.n_free += 1" -->

**Periodic (Bloch) boundaries** use the same map. In a periodic structure of period
$L$, Bloch's theorem says the field at the end equals the field at the start times a
phase: $u(L) = e^{ikL}u(0)$
([Bloch 1929](20-references.md#bloch-1929); [Floquet 1883](20-references.md#floquet-1883); [Joannopoulos et al. 2008](20-references.md#joannopoulos-2008)). So the last dof simply maps onto the first, with factor
$e^{ikL}$.
<!-- code: Setup::new "let tie" .. "DofMap::new" -->

With a Bloch tie, the basis function for dof 0 is $\psi = \phi_0 + e^{ikL}\phi_N$,
where $N$ is the last dof, and $\psi$ is complex. Complex functions pair as
$\int u\,\overline{v}\,dx$, so we test with $\overline\psi$. That puts $\overline{f}$
on the row side. During assembly, each cell $c$ adds

$$A_{r_ir_j} \mathrel{+}= \overline{f_i}\,\left(a_c S_c + b_c M_c\right)_{ij}\,f_j,$$
<!-- code: assemble_split "let Some((ri, fi))" .. "t.push" -->

where $i$ and $j$ are the cell's local dofs, $r_i$ is the reduced index of $i$, and
$f_i$ its factor. Testing with $\overline\psi$ also cancels the boundary terms at the
two ends: $\overline\psi$ is 1 at $x = 0$ and $e^{-ikL}$ at $x = L$, while
$u'(L) = e^{ikL}u'(0)$.

The conjugate keeps the matrix Hermitian, so its eigenvalues stay real, as a lossless
structure's frequencies must. Without it, the eigensolver of chapter 6 wouldn't apply.

Take order 1 and three cells of length $h$, with dofs 0 to 3 and $u_3 = e^{i\theta}u_0$,
where $\theta = kL$. The third cell's entries for dof 3 land in row and column 0, so

$$S = \frac{1}{h}\begin{pmatrix}2 & -1 & -e^{-i\theta}\\-1 & 2 & -1\\-e^{i\theta} & -1 & 2\end{pmatrix}.$$
<!-- code: DofMap::new "for &(slave, master, phase)" .. "d.map[slave]" -->

The two corners are conjugates, so $S$ is Hermitian. Entry $(0, 0)$ still sums to
$2/h$, because the third cell adds $\overline{e^{i\theta}}\,e^{i\theta}/h = 1/h$. But
the corners break the band: tying the last dof to the first connects the two ends of
the line. Chapter 5 deals with that.

Next: [solving linear systems](05-linear-systems.md).
