# 5. Solving linear systems

FEM's output is a sparse matrix equation $A\mathbf{x} = \mathbf{b}$. A general dense
solver costs $O(N^3)$, which for a million unknowns is hopeless. Structure is what
makes it fast.

## Banded matrices

With the numbering from chapter 3, $A_{ij}$ is zero unless $i$ and $j$ share a cell.
So all nonzeros sit within a narrow band around the diagonal. For order $r$ in 1D,
the half-width is $b = r$.
<!-- code: BandedLu::factor "let mut b = 0" .. "b = b.max" -->

Without row swaps, Gaussian elimination never creates nonzeros outside the band, since
it only combines rows that already overlap. So LU factorization on a band costs

$$O(N b^2)$$
<!-- code: BandedLu::factor "// Eliminate below the pivot" .. "*x -= *m * t" -->

which is **linear** in $N$ ([Golub & Van Loan 2013](20-references.md#golub-van-loan-2013)). Doubling the mesh doubles the time.

## Pivoting

Our matrices are often *indefinite*: $S - k_0^2 Q$ has both positive and negative
eigenvalues. Plain elimination can then hit a tiny or zero pivot and lose all its
accuracy, so each step picks the largest available entry in the column and swaps it
into place. That's **partial pivoting** ([Golub & Van Loan 2013](20-references.md#golub-van-loan-2013)).
<!-- code: BandedLu::factor "// Partial pivoting" .. "ab.swap" -->

If every candidate is zero, the matrix is singular. In a lossless structure that
happens exactly when $k_0$ is a resonance, where a driven field has no finite answer,
so `BandedLu::factor` reports it as an error.
<!-- code: BandedLu::factor "if ab[col + p] == ZERO" .. "singular matrix" -->

Swapping rows can push entries above the original band, up to $2b$: a row swapped up
from as far as $b$ rows below reaches $b$ columns past its own diagonal. The storage in
[`banded.rs`](../crates/glassworks-core/src/banded.rs) keeps a column of $3b + 1$
entries to leave room, the same layout LAPACK uses
([Anderson et al. 1999](20-references.md#anderson-1999)). Entry $(i, j)$ is stored at
row $2b + i - j$ of column $j$. The diagonal sits at row $2b$, $L$ in the $b$ rows
below it, and $U$ with room for fill in the $2b$ rows above. For a $5 \times 5$
matrix with $b = 1$, that's 4 rows of 5 columns:

```text
          col 0   col 1   col 2   col 3   col 4
row 0       .       .     (0,2)   (1,3)   (2,4)    zero at first, room for fill
row 1       .     (0,1)   (1,2)   (2,3)   (3,4)    upper diagonal
row 2     (0,0)   (1,1)   (2,2)   (3,3)   (4,4)    diagonal
row 3     (1,0)   (2,1)   (3,2)   (4,3)     .      lower diagonal, then L
```
<!-- code: BandedLu::factor "let (ld, diag)" .. "ab[diag + r + c * ld - c]" -->

## Reordering

The band only stays narrow if neighbours have nearby numbers. Two cases break that:
periodic boundaries in 1D, where the last dof connects to the first, and every 2D
mesh, whose numbering is arbitrary.

**Reverse Cuthill-McKee** (RCM) renumbers the unknowns to fix it. Start at a node
with few neighbours, then number nodes breadth-first, then reverse the order. Numbers
now spread out in rings, so neighbours stay close
([Cuthill & McKee 1969](20-references.md#cuthill-mckee-1969); [Liu & Sherman 1976](20-references.md#liu-sherman-1976); [George & Liu 1981](20-references.md#george-liu-1981)).
<!-- code: rcm -->

`rcm` is about 30 lines, and it elegantly fixes the periodic case. A ring of $N$ nodes
folds flat, with both sides of the ring numbered alternately, giving a bandwidth of
2. On a ring of 6 nodes, breadth-first from node 0 visits 0, 1, 5, 2, 4, 3, so ring
neighbours end up at most 2 apart, and reversing keeps that. The test
`rcm_folds_a_ring_to_a_narrow_band` confirms this.

`BandedLu::factor` then factors the renumbered matrix $PAP^T$, where the permutation
$P$ moves old dof `perm[k]` to position $k$. Since $A\mathbf{x} = \mathbf{b}$ is the same
as $PAP^T(P\mathbf{x}) = P\mathbf{b}$, `solve` permutes the right-hand side on the way
in and undoes it on the way out.
<!-- code: BandedLu::factor "let perm = rcm(a)" .. "inv[old] = new" -->

## How fast is this?

Here's our banded LU against [`faer`](https://docs.rs/faer), an excellent
general-purpose sparse library, on the same 1D complex matrix. These are native,
single-threaded timings on an Apple M4 Pro, the median of three runs of `cargo run --release -p glassworks-core --example bench`.
They vary between runs, so treat them as rough:

| Unknowns | Order | Banded LU | `faer` sparse LU |
|---|---|---|---|
| 100,000 | 1 | 8 ms | 15 ms |
| 1,000,000 | 1 | 65 ms | 142 ms |
| 1,000,000 | 4 | 111 ms | 175 ms |

A specialised solver beats a general one here because it needs no symbolic analysis
of the sparsity. That's specific to 1D, though. In 2D, the bandwidth of a mesh with
$N$ unknowns grows like $\sqrt{N}$, so the cost becomes $O(N^2)$, and fill-reducing
sparse solvers win by a wide margin ([George 1973](20-references.md#george-1973)). Chapter 10 covers that switch.

In the browser, a complete driven solve with a million unknowns took 158 ms on the same
machine, including assembly. Expect a slower result on a phone.

Next: [eigenproblems](06-eigenproblems.md).
