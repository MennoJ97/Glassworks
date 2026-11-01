# 2D solver plan

The agreed plan for adding 2D. Pick up at the first milestone that isn't marked done,
and update this file as things land. Delete it, and `plans/probe-2d/`, once the last
milestone is finished and the tutorial covers what's here.

**Status (2026-11-01):** nothing is implemented. The numbers below come from the
throwaway probe in [`probe-2d/`](probe-2d/). **Next up: milestone 1.**

## Decisions

These are settled. Reopen one only if a milestone's measurements contradict it.

| Topic | Decision | Why |
|---|---|---|
| Physics | One scalar field, `u = Ez` (`p = 1/mu`, `q = eps`) or `u = Hz` (`p = 1/eps`, `q = mu`) | Chapter 1's weak form carries over, and spurious modes can't arise |
| Vector fields | No Nédélec elements yet | Only full-vector waveguide cross-sections need them, and that's after this plan |
| Elements | `H1Tri`: second-order triangles implementing `Space` | Cavity error falls about 15× per mesh halving. Other orders can wait for a measured reason |
| Mesh | Structured grid, with the nearest nodes snapped onto shape boundaries | About 40 lines with no dependency, periodic edges match, and it can serve 2D topology optimization. See open problem 1 |
| Mesh fallback | An unstructured mesher, such as [`spade`](https://github.com/Stoeoef/spade) (MIT/Apache) | Only if snapping can't be fixed. It makes periodic edges and own-solver separators harder |
| Linear solver | Add a small factorization trait. 1D keeps `banded.rs`, 2D uses `faer` | See the measurements. Banded LU can't scale in 2D |
| Open boundaries | PML | Exact local conditions don't exist in 2D (chapter 9) |
| Leaky modes | No non-Hermitian eigensolver. Use one driven solve at complex frequency | Avoids Krylov-Schur. Gives LDOS directly, which is what Q/V is a proxy for |
| Upper size | The page refuses above about 160,000 unknowns | 638,401 failed in WebAssembly, probably out of memory |
| Target size | 300–400 more lines of Rust, about 1,500 in total (estimate) | State the cost before adding more, and update the tutorial's "about 1,000 lines" |
| Own sparse solver | Not now. Decide after milestone 1 (see below) | `faer` first, behind the trait |

## Measurements

From `plans/probe-2d`, using second-order triangles. Native on an Apple Silicon Mac,
single-threaded. The WebAssembly runs were in Node while the native probe was also
running, so they may read slow.

Solving a 2D PML Helmholtz matrix once. Native times are the solve alone, and
WebAssembly times include assembly:

| Grid | Unknowns | Banded, native | `faer`, native | Banded, Wasm | `faer`, Wasm |
|---|---:|---:|---:|---:|---:|
| 25 × 25 | 2,401 | 47 ms | 8.6 ms | 238 ms | 34 ms |
| 50 × 50 | 9,801 | 591 ms | 36 ms | 911 ms | 100 ms |
| 100 × 100 | 39,601 | 9.7 s | 0.31 s | 16.6 s | 0.83 s |
| 200 × 200 | 159,201 | Needs 12 GB | 3.0 s | Failed | 8.7 s |
| 400 × 400 | 638,401 | Needs 98 GB | 27 s | Failed | Failed |
| 800 × 800 | 2,556,801 | Needs 784 GB | 345 s | Failed | Failed |

Download size, gzipped: the current page's solver is about 110 KB, and `faer` adds about
88 KB. Measured on the probe's own build (33 KB without `faer`, 121 KB with it).

Eigensolves: three modes at one k-point took 14 ms at 576 unknowns, about 200 ms at
2,304 and about 2.7 s at 9,216, all on banded LU. At 9,216 that's several times one
factorization, because Lanczos adds back-substitutions and reorthogonalization.

Band edges of the rod lattice (eps 12, r = 0.2a, square lattice), as error against the
probe's own 96 × 96 snapped mesh. M1 and X2 are the two edges of the TM gap:

| Grid (unknowns) | TM staircase | TM snapped | TE staircase | TE snapped |
|---|---:|---:|---:|---:|
| 12 × 12 (576) | 1.3e-2 | 3.6e-3 | 3.9e-2 | 7.0e-3 |
| 24 × 24 (2,304) | 6.7e-4 | 9.0e-4 | 9.3e-3 | 1.6e-3 |
| 48 × 48 (9,216) | 5.9e-4 | 9.6e-5 | 9.7e-3 | 1.2e-4 |

The reference values were TM 0.28068 (M1) and 0.41721 (X2), and TE 0.49594 (band 1 at M)
and 0.44131 (band 2 at X). These are worst-of-two errors per cell, and they only show
self-consistency until open problem 2 is settled.

A line source against the exact `(i/4) H0(kr)`, on rings of r = 0.8 to 1.2 wavelengths:

| PML thickness | Cells per vacuum wavelength | Worst error |
|---:|---:|---:|
| 0.25 λ | 12 | 2.24% |
| 0.5 λ | 6 | 2.26% |
| 0.5 λ | 12 | 0.37% |
| 0.5 λ | 24 | 0.09% |
| 1 λ | 12 | 0.14% |

So 0.5 λ of PML at 12 cells per vacuum wavelength is the starting default. Cells inside
high-index material will need refining, which the probe didn't test.

## Open problems

1. **Slivers from snapping.** The smallest triangle angle falls from 23° to 4° to 0.3°
   on the 12, 24 and 48 grids. Accuracy held up, but this isn't safe to build on. Try
   limiting how far a node may move, then flipping bad diagonals. Milestone 1 can't
   finish without this.
2. **Unexplained gap to MPB.** The probe's TM gap is 0.28068 to 0.41721. The
   [MPB tutorial](https://mpb.readthedocs.io/en/latest/Python_Tutorial/) gives 0.28262
   to 0.41933, about 0.5–0.7% higher, from a run at resolution 32 that may simply be
   under-resolved. Settle it by getting a converged reference, such as a
   high-resolution MPB run or an independent plane-wave calculation.
3. **The `faer` setup is untuned.** Its default ordering, symmetric factorization for
   `Aᵀ = A`, and reusing symbolic analysis across solves are all untested.
4. **Transposed solves through `faer`.** Bloch-periodic driven problems lose `Aᵀ = A`,
   so their adjoints need `Aᵀ` solves. Not needed before periodic gratings.

## Code changes

| Module | Change |
|---|---|
| `mesh.rs` | A structured grid and snapping. `coords` holds x, y per vertex. Finally reads `verts_per_cell` |
| `space.rs` | `H1Tri`. `local_matrices` returns `S_xx`, `S_yy` and `M` separately, with `S_yy = 0` in 1D |
| `assemble.rs` | `coef` returns `(a_x, a_y, b)` for `a_x S_xx + a_y S_yy + b M`, which PML needs. `DofMap` takes many ties, including corners with two Bloch phases |
| A new solver trait | One factor-and-solve interface, used by `eigs_near` and `driven`. `banded.rs` and `faer` implement it |
| `problem.rs` | A tagged geometry: the current layers, or a 2D domain with background, shapes, lattice and PML |
| `solve.rs` | Boundaries become edges. Field output adds node positions and triangles |
| Website | A shape editor and a triangle field renderer alongside `FieldPlot`. The `design.js` optimizer keeps projection and MMA, with a 2D filter |
| Tutorial | Correct chapter 9 and add 2D chapters as each milestone lands (see "Tutorial changes") |

The probe's `H1Tri`, `snap`, `bloch_map` and `helmholtz_pml` are good starting points,
but they aren't production code. `bloch_map` keys on the unsnapped positions, because
snapped nodes never sit on the periodic edges.

PML, as used in the probe, with `s = 1 + i sigma(d) / k` and a quadratic `sigma`:

```text
a_x = s_y / s_x,  a_y = s_x / s_y,  b = -k^2 s_x s_y
sigma_max = -3 ln(1e-6) / (2 T),  sigma(d) = sigma_max (d / T)^2
```

## Milestones

Each one adds a capability and ends with a check against a known answer. Add the check
as a test in `crates/eigenmenno-core/tests/analytic.rs` before calling it done.

### Milestone 1: bands on a good mesh

- [ ] Fix the slivers (open problem 1), with a test on the minimum angle.
- [ ] Settle the MPB gap (open problem 2) and write the result here.
- [ ] Add `H1Tri`, the grid and snapping, the many-tie `DofMap` and the solver trait to
      the core, with `faer` as a real dependency.
- [ ] **Example 1: square cavity.** Check against the exact `k^2 = pi^2 (m^2 + n^2)`.
- [ ] **Example 2: rod-lattice band diagram**, TM and TE along Γ-X-M-Γ. Check against
      the converged reference from the MPB work.
- [ ] Add a triangle field renderer and a minimal shape input on the page.
- [ ] Decide on the own sparse solver (see below).

Done when a visitor can see a correct band diagram and click a point to see its mode.
Aim for about 2,300 unknowns per k-point at first, which is accurate to about 0.1–0.3%.

### Milestone 2: open structures

- [ ] Add the split-stiffness PML, and a line source.
- [ ] Add a plane wave through a scattered-field formulation.
- [ ] **Example 3: dielectric cylinder scattering.** Check against the analytic Mie
      series, and the line source against `(i/4) H0(kr)`.

### Milestone 3: cavities

- [ ] **Example 4: point-defect cavity.** One rod removed from a supercell. First check
      that the mode falls inside the gap.
- [ ] Add a driven solve at complex frequency, for frequency-averaged LDOS. With
      `A u = i k e_0` and `F = 2 Re(-u_0)` in 1D units, the adjoint needs no extra solve,
      because `lambda = u / (i k)`. Check that LDOS against a Lorentzian average of
      real-frequency solves.
- [ ] Add Q and V for a finite crystal. This is where Q/V becomes meaningful, because
      in 1D V barely changes and Q/V is nearly just Q.

### Milestone 4: waveguides

- [ ] Add a waveguide port, taking the mode profile from the existing 1D solver.
- [ ] Add power monitors along lines.
- [ ] **Example 5: sharp bend in a photonic-crystal waveguide.** Check against
      [Mekis et al. 1996](https://doi.org/10.1103/PhysRevLett.77.3787) (98% transmission,
      against 30% in a conventional guide). Read the paper for its exact parameters
      first: only the abstract-level numbers have been checked so far.

### Milestone 5: 2D inverse design

- [ ] Reuse `design.js` with a 2D filter, and density per cell. Add 2D adjoint gradients
      in the solver, and keep `A` symmetric for PML problems.
- [ ] **Example 6: topology-optimized bend.** Compare against the hand-designed bend
      from example 5. [Jensen & Sigmund 2004](https://doi.org/10.1063/1.1688450) is the
      precedent.
- [ ] Optionally, optimize the cavity from example 4 by raising the LDOS bandwidth's Q
      in stages, as in [Liang & Johnson 2013](https://doi.org/10.1364/OE.21.030812).

### Later, not part of this plan

Periodic gratings (a port plus Bloch ties, which needs open problem 4), and full-vector
waveguide cross-sections, which are the first real need for Nédélec elements.

## Own sparse solver

Decide this at the end of milestone 1, by benchmarking a prototype against `faer` on
real 2D matrices. Keep `faer` as a test reference either way.

**Recommended candidate:** a multifrontal LU with geometric nested dissection, estimated
at 250–400 lines. On our structured grid the separators are just rows and columns of
nodes, and the elimination tree can be reused across wavelengths, k-points and optimizer
iterations.

- **Adopt it if** it's within a few times of `faer` in WebAssembly and stays reliable on
  hard matrices, such as lossless structures near resonance.
- **Main risk:** pivoting. Pivoting only inside each dense block can fail on those
  matrices, and deferring bad pivots is where the complexity lives.
- **Blocker:** an unstructured mesh removes the easy separators, so this depends on
  open problem 1.
- **Why bother:** about 88 KB less to download, no major runtime dependency, readers
  can still read the whole solver, and it would make a strong tutorial chapter following
  chapter 5's RCM.

A simple entry-by-entry sparse LU (about 200 lines) is likely several times slower than
`faer`. Matching `faer` outright isn't worth attempting.

## Tutorial changes

Chapter 9 says a few things these measurements correct. Fix them in milestone 1:

- Banded LU is already slow at 10,000 unknowns: about 0.6 s native, and more once an
  eigensolve multiplies it.
- PML doesn't quite slot into the existing assembly, because it needs the split stiffness.
- Leaky modes don't need `faer`'s Krylov-Schur for what's planned here.
- The README and `docs/README.md` give the solver's size and describe `faer` as tests-only.

## Running the probe

The probe isn't part of the workspace, so it can't affect normal builds. The full run
takes about 10 minutes, most of it in the largest solves at the end.

```bash
cd plans/probe-2d && CARGO_TARGET_DIR=/tmp/eigenmenno-probe cargo run --release --features sparse --bin probe
```

For the WebAssembly timings, build with the rustup toolchain and run `run.mjs`. It expects
the output copied next to it as `sparse.wasm`:

```bash
cd plans/probe-2d && CARGO_TARGET_DIR=/tmp/eigenmenno-probe RUSTFLAGS="-C target-feature=+simd128" cargo build --release --target wasm32-unknown-unknown --lib --features sparse
```
