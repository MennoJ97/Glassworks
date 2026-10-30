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

Where two cells share a dof, their contributions add. You could see this in the
worked example in chapter 3, where the middle node received $1/h$ from each side.

This loop is the function `assemble` in
[`assemble.rs`](../crates/eigenmenno-core/src/assemble.rs), and it doesn't know about
1D. It only sees the `Space` trait, which supplies dof numbers and local matrices.
That's what makes 2D an addition rather than a rewrite.

## Storing a sparse matrix

The loop produces a list of `(row, column, value)` triplets, including duplicates.
`Csr::from_triplets` sorts them and adds the duplicates together, producing
**compressed sparse row** storage: all the column numbers and values in one long
array, and a pointer to where each row starts ([Davis 2006](10-references.md#davis-2006)). Multiplying by a vector then takes
one pass through the nonzeros.

## Two kinds of boundary condition

Recall the boundary term from the weak form: $-[p u' v]_a^b$.

**Essential (Dirichlet) conditions** fix the value of $u$, for example $u = 0$. You
have to impose these explicitly, by removing that unknown from the system.

**Natural (Neumann) conditions** specify $p u'$. These need no special handling: if
$p u' = 0$ at a boundary, the boundary term vanishes, and leaving it out of the
matrix *is* the boundary condition. Doing nothing imposes it ([Strang & Fix 2008](10-references.md#strang-fix-2008)).

Physically:

| Wall | TE ($u = E_y$) | TM ($u = \eta_0 H_y$) |
|---|---|---|
| Perfect conductor (PEC) | $E_y = 0$: essential | $dH_y/dx = 0$: natural |
| Magnetic wall (PMC) | $dE_y/dx = 0$: natural | $H_y = 0$: essential |

A perfect conductor forces the tangential electric field to zero
([Jackson 1999](10-references.md#jackson-1999)). In TE that's the
unknown itself. In TM the same wall forces $E_z \propto du/dx$ to zero, which is a
condition on the derivative, so it's natural. The test
`tm_polarization_swaps_pec_and_pmc` checks this: a TM cavity with PEC walls has a
static mode at $k_0 = 0$, as the theory says.

## One mechanism for all constraints

`DofMap` in [`assemble.rs`](../crates/eigenmenno-core/src/assemble.rs) gives each
global dof a reduced index and a factor, or marks it as fixed to zero:

$$u_\text{global} = \text{factor} \times u_\text{reduced}$$

A Dirichlet wall sets a dof to `None`. That's the whole mechanism.

**Periodic (Bloch) boundaries** use the same map. In a periodic structure of period
$L$, Bloch's theorem says the field at the end equals the field at the start times a
phase: $u(L) = e^{ikL}u(0)$
([Bloch 1929](10-references.md#bloch-1929); [Floquet 1883](10-references.md#floquet-1883); [Joannopoulos et al. 2008](10-references.md#joannopoulos-2008)). So the last dof simply maps onto the first, with factor
$e^{ikL}$.

During assembly, each entry becomes

$$A_{ij} = \overline{f_i}\,\left(a S + b M\right)_{ij}\,f_j.$$

The complex conjugate on the test side comes from using a conjugated test function.
It keeps the matrix Hermitian, so its eigenvalues stay real, as a lossless structure's
frequencies must. Without it, the eigensolver of chapter 6 wouldn't apply.

One wrinkle: tying the last dof to the first connects the two ends of the line, so
the matrix is no longer banded. Chapter 5 deals with that.

Next: [solving linear systems](05-linear-systems.md).
