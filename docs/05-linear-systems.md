# 5. Solving linear systems

FEM's output is a sparse matrix equation $A\mathbf{x} = \mathbf{b}$. A general dense
solver costs $O(N^3)$, which for a million unknowns is hopeless. Structure is what
makes it fast.

## Banded matrices

With the numbering from chapter 3, $A_{ij}$ is zero unless $i$ and $j$ share a cell.
So all nonzeros sit within a narrow band around the diagonal. For order $p$ in 1D,
the half-width is $b = p$.

Gaussian elimination never creates nonzeros outside the band, since it only combines
rows that already overlap. So LU factorization on a band costs

$$O(N b^2)$$

which is **linear** in $N$ ([Golub & Van Loan 2013](19-references.md#golub-van-loan-2013)). Doubling the mesh doubles the time.

## Pivoting

Our matrices are often *indefinite*: $S - k_0^2 Q$ has both positive and negative
eigenvalues. Plain elimination can then hit a tiny or zero pivot and lose all its
accuracy, so each step picks the largest available entry in the column and swaps it
into place. That's **partial pivoting** ([Golub & Van Loan 2013](19-references.md#golub-van-loan-2013)).

Swapping rows can push entries above the original band, up to $2b$. The storage in
[`banded.rs`](../crates/eigenmenno-core/src/banded.rs) keeps a column of $3b + 1$
entries to leave room, the same layout LAPACK uses
([Anderson et al. 1999](19-references.md#anderson-1999)).

## Reordering

The band only stays narrow if neighbours have nearby numbers. Two cases break that:
periodic boundaries in 1D, where the last dof connects to the first, and every 2D
mesh, whose numbering is arbitrary.

**Reverse Cuthill-McKee** (RCM) renumbers the unknowns to fix it. Start at a node
with few neighbours, then number nodes breadth-first, then reverse the order. Numbers
now spread out in rings, so neighbours stay close
([Cuthill & McKee 1969](19-references.md#cuthill-mckee-1969); [Liu & Sherman 1976](19-references.md#liu-sherman-1976); [George & Liu 1981](19-references.md#george-liu-1981)).

It's about 30 lines, and it elegantly fixes the periodic case. A ring of $N$ nodes
folds flat, with both sides of the ring numbered alternately, giving a bandwidth of
2. The test `rcm_folds_a_ring_to_a_narrow_band` confirms this.

## How fast is this?

Here's our banded LU against [`faer`](https://docs.rs/faer), an excellent
general-purpose sparse library, on the same 1D complex matrix. These are native,
single-threaded timings from `cargo run --release -p eigenmenno-core --example bench`.
They vary between runs, so treat them as rough:

| Unknowns | Order | Banded LU | `faer` sparse LU |
|---|---|---|---|
| 100,000 | 1 | 8 ms | 16 ms |
| 1,000,000 | 1 | 97 ms | 190 ms |
| 1,000,000 | 4 | 164 ms | 265 ms |

A specialised solver beats a general one here because it needs no symbolic analysis
of the sparsity. That's specific to 1D, though. In 2D, the bandwidth of a mesh with
$N$ unknowns grows like $\sqrt{N}$, so the cost becomes $O(N^2)$, and fill-reducing
sparse solvers win by a wide margin ([George 1973](19-references.md#george-1973)). Chapter 10 covers that switch.

In the browser, a complete driven solve with a million unknowns took 217 ms on an
Apple M4 Pro, including assembly. Expect a slower result on a phone.

Next: [eigenproblems](06-eigenproblems.md).
