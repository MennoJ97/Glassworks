# 2D solver plan

The agreed plan for adding 2D. Pick up at the first milestone that isn't marked done,
and update this file as things land. Delete it, with `plans/probe-2d/` and `plans/probe-sparse/`, once the last
milestone is finished and the tutorial covers what's here.

**Status (2026-11-01):** milestones 1 to 4 are built. Cavities, band diagrams, open
scattering, leaky modes and waveguide bends run on the page, tested against exact and
independent results. **Next up: milestone 5.** The solver numbers below come from the throwaway probe in
[`probe-2d/`](probe-2d/).

## Decisions

These are settled. Reopen one only if a milestone's measurements contradict it.

| Topic | Decision | Why |
|---|---|---|
| Physics | One scalar field, `u = Ez` (`p = 1/mu`, `q = eps`) or `u = Hz` (`p = 1/eps`, `q = mu`) | Chapter 1's weak form carries over, and spurious modes can't arise |
| Vector fields | No Nédélec elements yet | Only full-vector waveguide cross-sections need them, and that's after this plan |
| Elements | `H1Tri`: second-order triangles implementing `Space` | Cavity error falls about 15× per mesh halving. Other orders can wait for a measured reason |
| Mesh | Structured grid, with nodes snapped onto shape boundaries, diagonals chosen afterwards, then constrained Delaunay flips | No dependency, periodic edges match, and it can serve 2D topology optimization. Worst angle over 200 random placements: 22° for circles, 16° for rectangles |
| Mesh fallback | An unstructured mesher, such as [`spade`](https://github.com/Stoeoef/spade) (MIT/Apache) | Only if snapping can't be fixed. It makes periodic edges and own-solver separators harder |
| Linear solver | Add a small factorization trait. 1D keeps `banded.rs`, 2D uses `sparse.rs` | See the measurements. Banded LU can't scale in 2D |
| Open boundaries | PML | Exact local conditions don't exist in 2D (chapter 9) |
| Leaky modes | No non-Hermitian eigensolver. One driven solve at complex frequency for LDOS, and a secant search on `1/u_0` for Q | Avoids Krylov-Schur. Gives LDOS directly, which is what Q/V is a proxy for. It took 6 or 7 solves, and found Q up to 680,000 |
| Upper size | The page refuses above about 160,000 unknowns | 638,401 failed in WebAssembly, probably out of memory |
| Size | The 300–400 line estimate was low: milestone 1 alone added about 600 lines of Rust, for about 1,700 at that point, about 1,900 after milestone 2, about 2,000 after milestone 3, and about 2,150 after milestone 4 | State the cost before adding more, and keep the tutorial's line count current |
| Own sparse solver | Adopted in `sparse.rs`. `faer` is now only for tests (see below) | It's 2–6× faster than `faer` on driven problems in WebAssembly and within 8% on band diagrams, with 99 KB less to download. Refinement covers its one real risk |

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

Download size, gzipped: the page's solver went from about 110 KB to 222 KB with milestone
1, and to 270 KB with sparse LU in milestone 2, more than the 88 KB the probe's own build suggested (33 KB without `faer`, 121 KB with).

Eigensolves: three modes at one k-point took about 200 ms at 2,304 unknowns and 2.7 s at
9,216 on banded LU. A full band diagram (16 k-points, 6 bands) took, natively:

| Unknowns | `faer` LU | LU, reused analysis | Cholesky, reused analysis |
|---:|---:|---:|---:|
| 2,500 | 1.28 s | 1.18 s | 0.48 s |
| 10,000 | 9.3 s | 6.8 s | 2.9 s |

Cholesky factors about 8× faster than LU at 2,500 unknowns, and its back-substitutions
are about 3× faster. Reusing the analysis saves little once that's in place, about 10%.
Lanczos now dominates: about 62 back-substitutions per k-point, across two runs, with
reorthogonalization about as costly as the solves. On the page, a band diagram at 2,500
unknowns went from about 1.8 s with LU to about 0.6 s.

Band edges of the rod lattice (eps 12, r = 0.2a, square lattice) from the probe, as error
against its own 96 × 96 snapped mesh. M1 and X2 are the two edges of the TM gap. The
staircase columns show why fitting matters, and the final mesher's values follow:

| Grid (unknowns) | TM staircase | TM snapped | TE staircase | TE snapped |
|---|---:|---:|---:|---:|
| 12 × 12 (576) | 1.3e-2 | 3.6e-3 | 3.9e-2 | 7.0e-3 |
| 24 × 24 (2,304) | 6.7e-4 | 9.0e-4 | 9.3e-3 | 1.6e-3 |
| 48 × 48 (9,216) | 5.9e-4 | 9.6e-5 | 9.7e-3 | 1.2e-4 |

These are worst-of-two errors per cell. With the final mesher, the band edges converge
from above as straight cells cut the corners of the curved rod:

| Grid (unknowns) | Ez M1 | Ez X2 | Hz M1 | Hz X2 |
|---|---:|---:|---:|---:|
| 12 × 12 (576) | 0.28456 | 0.42346 | 0.50783 | 0.44746 |
| 24 × 24 (2,304) | 0.28138 | 0.41835 | 0.49786 | 0.44235 |
| 48 × 48 (9,216) | 0.28084 | 0.41747 | 0.49637 | 0.44155 |
| 96 × 96 (36,864) | 0.28071 | 0.41724 | 0.49600 | 0.44135 |
| Plane waves (2,401), see below | 0.28066 | 0.41716 | 0.49588 | 0.44109, still rising |

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

1. **Slivers from snapping: solved.** Slivers came from diagonal edges forcing two
   neighbours across a curve onto almost the same point. Now only grid lines force a
   move, and each cell's diagonal is chosen afterwards. With no grid line crossing an
   outline, one diagonal is always clear. Flips then clean up, but may not cross an
   outline.
2. **Gap to MPB: settled.** An independent plane-wave expansion gives a TM gap of 0.28066
   to 0.41716, which this solver converges to. The
   [MPB tutorial](https://mpb.readthedocs.io/en/latest/Python_Tutorial/)'s 0.28262 to
   0.41933 is 0.5–0.7% higher. The likely cause is its coarse resolution of 32, but no
   MPB run has confirmed that, so it isn't used as a reference.
3. **Band diagram speed.** Eigenproblems use the LU with a reused ordering, within 8% of `faer`'s
   Cholesky. Warm-starting Lanczos from the previous k-point is the next speed-up to try.
4. **Transposed solves.** Bloch-periodic driven problems lose `Aᵀ = A`,
   so their adjoints need `Aᵀ` solves. Not needed before periodic gratings.

## Code changes

| Module | Change |
|---|---|
| `solve2d.rs`, later | Gradients |
| Website, later | Dragging shapes, and the `design.js` optimizer with a 2D filter |
| Tutorial | Add 2D chapters as each milestone lands. Chapter 9 is current |

Milestone 1 added `Mesh::fitted_grid`, `H1Triangle`, many-tie `DofMap`, the `Factors` trait
with `SparseLu`, and `solve2d.rs`. 2D requests are a separate `Request2d`, sent to
`solve_2d`, rather than a tagged geometry. On the page, `FieldMap` sits alongside
`FieldPlot`. Milestone 2 added `assemble_split` for `a_x S_xx + a_y S_yy + b M`, `SparseLu`
and `Analysis2d::Driven`. The absorbing layer sits outside the caller's domain, and fields
there are returned at negative or beyond-the-domain coordinates, which the plot clips.

PML, as used in the probe, with `s = 1 + i sigma(d) / k` and a quadratic `sigma`:

```text
a_x = s_y / s_x,  a_y = s_x / s_y,  b = -k^2 s_x s_y
sigma_max = -3 ln(1e-4) / (2 T),  sigma(d) = sigma_max (d / T)^2
```

## Milestones

Each one adds a capability and ends with a check against a known answer. Add the check
as a test in `crates/eigenmenno-core/tests/analytic.rs` before calling it done.

### Milestone 1: bands on a good mesh

- [x] Fix the slivers (open problem 1), with a test on the minimum angle.
- [x] Settle the MPB gap (open problem 2) and write the result here.
- [x] Add `H1Triangle`, the grid and snapping, the many-tie `DofMap` and the solver trait
      to the core, with `faer` as a real dependency.
- [x] **Example 1: square cavity.** Checked against the exact mode formula, and a layered
      box against the 1D solver.
- [x] **Example 2: rod-lattice band diagram**, Ez and Hz along Γ-X-M-Γ. Checked against
      the plane-wave reference.
- [x] Add a triangle field renderer and a minimal shape input on the page.
- [x] Decide on the own sparse solver: keep `faer` for now, and revisit after milestone 4.

The page defaults to about 2,500 unknowns per k-point, accurate to about 0.3–0.5%.
The plane-wave script wasn't kept: it's about 30 lines
of numpy, using `eps(G) = delta eps_b + f (eps_a - eps_b) 2 J1(Gr) / (Gr)`, and the inverse
rule for Hz.

### Milestone 2: open structures

- [x] Add the split-stiffness PML, and a line source.
- [x] Add a plane wave through a scattered-field formulation.
- [x] **Example 3: dielectric cylinder scattering.** Checked against the analytic Mie
      series, within 1% at 32 cells per wavelength for both polarizations and at an oblique
      angle. The line source matches `(i/4) H0(kr)` to 0.02% with 1 λ of PML at 24 cells
      per wavelength, including at the corners, and to 0.05% since milestone 4 weakened the
      layer's target.
- [x] Add a scattering view to the page, with the source and travel direction marked.

Materials are still real and positive, so metals and loss aren't available in 2D yet. A
plane wave needs its shapes inside the domain, because the incident wave is only defined
in the background. Waveguides that run into the layer will need the ports in milestone 4.

### Milestone 3: cavities

- [x] **Example 4: point-defect cavity.** One rod removed from a supercell. The mode falls
      inside the gap, at a/lambda 0.3746 by an independent plane-wave expansion (4489 waves,
      using the same 30 lines of numpy as milestone 1). The solver gives 0.37593, 0.37503 and
      0.37485 at 16, 24 and 32 cells per period, for the 5 x 5 supercell.
- [x] Add a driven solve at complex frequency, for frequency-averaged LDOS, with
      `A u = i k e_0` and `F = 2 Re(-u_0)`. For milestone 5, the adjoint needs no extra solve,
      because `lambda = u / (i k)`. Checked against real-frequency solves to 2.4e-4, using
      `2 F(gamma) - F(2 gamma)`, whose kernel falls off as 1 / k^4. One Lorentzian alone
      missed by 2%, which went away with that faster tail, so the cause was the sharp modes
      above the mesh's cutoff, which quadrature can't resolve.
- [x] Add Q and V for a finite crystal. Q is checked against the half-height of the
      real-frequency LDOS, for 3 x 3 and 5 x 5 crystals.

Q was within 5% across meshes from 10 to 24 cells per period, and changed by under 0.2%
between 0.25 and 1 wavelength of layer. Mesh alignment moves a/lambda by about 0.0005. V
stays at about 0.65 a^2 from 3 x 3 rods on, so adding rods changes Q/V only through Q. V
should matter when comparing different cavities, as milestone 5 would. No published Q has
been compared yet, so Q rests on the internal LDOS check.

On the page, 5 x 5 rods at 16,129 unknowns take 1 to 4 s, for 6 or 7 factorizations, slowest just after loading. Reusing
the analysis across them is the obvious speed-up. Download: 277 KB gzipped.

### Milestone 4: waveguides

- [x] Add a waveguide port, taking the mode profile from the existing 1D solver. Checked
      against the exact slab mode: field within 1.2%, from error in beta that builds up with
      distance, and power within 0.5%.
- [x] Add power monitors along lines, averaged over two cells, since a grid line is where
      a derivative is worst. Checked against `F / 4` around a line source.
- [x] **Example 5: sharp bend in a photonic-crystal waveguide.** Checked against
      [Mekis et al. 1996](https://doi.org/10.1103/PhysRevLett.77.3787), read in full. Chapter
      9 has the table. Sharp bend: 97.9% at 0.353 against their 98% maximum. Stepped bend:
      99.8% against 100%. Band: 0.312 and 0.445, against 0.312 and 0.443.

The slab comparison is weak, because the paper gives no width for its sharp bend. It reads 23%
at their other bend's width over wavelength, against "at most 30%", but 57% at 0.06. Near the
guide's cutoff, below about 0.35, reflection from the layer returns: with 4a of layer, T reads
1.02 at 0.33. The paper reports trouble there too. T was stable to about 0.3% from 8 to 16
cells per period at 0.353, and with cladding from 3 to 5 rows.

The layer needed most of the work. Rods running through it reflect the guided mode: 1.6% of
its power at 2a, 0.06% at 4a, and 0.001% at 6a, at 0.353. A weaker target of 1e-3 helps
further at 4a, but fails the line-source test. So the target went from 1e-6 to 1e-4 for
everything, which leaves all earlier tests passing. The test uses 6a. Monitors and
ports both assume a straight line, and a port needs its line to cross a layered guide, so a
port can't launch a photonic-crystal guide's mode: use a line source there, which excites only
that mode because nothing else exists in the gap.

On the page: 5a of layer at a/6 cells, 51,529 unknowns and two solves in about 5 s with `faer`,
within 0.3 percentage points of the test's values. The solver adopted afterwards should take
about 2.5 s by the ratio in Node, but that hasn't been timed on the page on wall power. The first solve after loading took about 20 s, in a
hidden browser pane, and milestone 3's first solve was also about 4 times its later ones. The
cause hasn't been measured: WebAssembly compiling its hot loops is the likely one. The slab takes about 3 s. Spectra would need many
solves, so the page shows one wavelength at a time, and T can't yet be plotted against
frequency there. Download: 283 KB gzipped.

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

**Adopted 2026-11-01, in `sparse.rs`, and `faer` is now only a test reference.** It's a
multifrontal LU with geometric nested dissection, 303 lines against the 250–400 estimated, and the
benchmarks are in [`probe-sparse/`](probe-sparse/). It takes separators from the matrix's graph,
not the grid, so Delaunay flips and periodic ties need no special care.

This was to be adopted "if it's within a few times of `faer` in WebAssembly and stays reliable
on hard matrices". It does better than that:

- **Speed:** factoring is 5–29× faster than `faer`'s LU in WebAssembly, and whole driven and leaky
  requests are 2–6× faster. Band diagrams, which had `faer`'s Cholesky, are 2–8% slower.
- **Reliability:** all 38 physics tests pass unchanged, including the Mekis bends and the Mie
  series. Q and T agree with `faer` to about 1e-13, including Q = 740,000.
- **Download:** 283 KB to 184 KB gzipped, 99 KB less. `faer` itself is 114 KB, not the 88 KB
  estimated, measured by building with a stub in its place.
- **Memory:** peak use during factoring is 15–30% of `faer` LU's, which could raise the 160,000
  limit. That hasn't been tried.
- **Cost:** about 240 lines of Rust net. Readers can now read everything.

Iterative refinement runs only when a front reports a weak pivot, and it's needed: see the risk
below. On real matrices it's rare, with the a/6 bend the only benchmark case that needed it. Open
problem 4 should also be simple now, since transposed solves are in our own code.

### No Hermitian path

Band diagrams factor a Hermitian positive definite matrix, so a Cholesky path, with no pivoting and
half the storage, looked like the way to win them back. It was built and measured, and **left out**:

| Band diagram, WebAssembly | `faer` | LU, first version | LU, thin separators | Thin, with Cholesky |
|---:|---:|---:|---:|---:|
| 2,500 unknowns | 761 ms | 1,049 ms | 840 ms | 802 ms |
| 10,000 unknowns | 3,809 ms | 5,467 ms | 3,891 ms | 3,600 ms |

Cholesky halved factoring, but solves don't get faster, since `L` and `U` have the same number of
entries. It saved 4–7% for 68 lines and 1.6 KB, and changes nothing else. `bands.rs` showed the
real gap: the first version's solves took 1.6× as long as `faer`'s, which meant too many entries.

A second-order element couples nodes two columns apart, so a cut through an element's middle made a
separator two columns thick. Fronts cost the cube of their size. Moving any separator node with
nothing further left into the right part, about 10 lines, and shrinking leaves from 64 to 16,
did far more than Cholesky could, for every problem: it halved peak memory, and sped up the leaky
13 x 13 by another 2×. Leaf sizes from 8 to 24 made little difference.

If band diagrams ever need more, Lanczos itself costs about as much as the solves, so
warm-starting it from the previous k-point is the next thing to try.

### Measurements

These are all from one run on battery, which reads about 1.5× slower than wall power. Each case is
a real matrix from the solver, on an Apple Silicon Mac, single-threaded. "Ours" includes its
ordering. WebAssembly is in Node, and includes one solve and a residual check. Peak memory counts
allocations during factoring.

| Matrix | Unknowns | Factor, native: `faer` LU / ours | Peak MB: `faer` LU / ours | Wasm: `faer` LU / ours |
|---|---:|---:|---:|---:|
| Sharp bend, a/6, 5a layer | 51,529 | 921 / 262 ms | 620 / 116 | 2,239 / 311 ms |
| Sharp bend, a/8, 6a layer | 112,225 | 2,693 / 746 ms | 1,339 / 262 | 6,866 / 851 ms |
| 7 x 7 cavity, real k at its peak | 27,889 | 271 / 103 ms | 181 / 55 | 621 / 114 ms |
| 13 x 13 cavity, on its pole | 82,369 | 1,693 / 455 ms | 785 / 186 | 4,364 / 529 ms |
| 5 x 5 supercell, shift inside the gap | 25,600 | 1,178 / 148 ms | 376 / 65 | 3,308 / 156 ms |
| 5 x 5 supercell, shift on the mode | 25,600 | 1,173 / 147 ms | 376 / 65 | 3,280 / 157 ms |
| Empty periodic box near a mode | 40,000 | 3,001 / 277 ms | 678 / 103 | 8,577 / 297 ms |
| Rod lattice at M, definite | 10,000 | 226 / 43 ms | 103 / 23 | 605 / 45 ms |
| Rod lattice at M, definite | 40,000 | 1,918 / 277 ms | 636 / 103 | 5,406 / 303 ms |

On the two definite cases, `faer`'s Cholesky took 16 and 96 ms, and 18 and 73 MB, still ahead.
Backward errors were 7e-17 to 2e-15 for both, apart from 4e-14 on the a/6 bend before refinement,
which brought it to 1e-16. No pivot was replaced.

Whole requests in WebAssembly, settled after warming up:

| Request | Unknowns | `faer` | Ours |
|---|---:|---:|---:|
| Band diagram, the page's default | 2,500 | 778 ms | 830 ms |
| Band diagram | 10,000 | 3,794 ms | 3,860 ms |
| Box, 6 lowest modes | 2,401 | 55 ms | 60 ms |
| Supercell near a target | 10,000 | 1,438 ms | 520 ms |
| Leaky 5 x 5, the page's default | 16,129 | 1,531 ms | 476 ms |
| Leaky 13 x 13 | 82,369 | 31,646 ms | 5,122 ms |
| Plane-wave scattering, the page's default | 19,321 | 400 ms | 101 ms |
| Sharp bend a/6, one of its two solves | 51,529 | 3,747 ms | 1,909 ms |

About 1.6 s of that last 1.9 s is outside the factor and solve. Meshing around 420 rods, assembly
and monitors weren't separated, and they're the next thing to look at for bends. Warm-up is
real too: the first factor in a fresh process took 818 ms against 311 ms afterwards, which
probably explains the slow first solve on the page.

### The main risk: pivoting only inside each front

A matrix can be fine as a whole while a subdomain, walled in by its separator, sits on a
resonance of its own. Elimination there meets a nearly singular block, and no later unknown can
help. `src/bin/stress.rs` builds this on a 9,801-unknown closed box, shifted to within `delta` of
a subdomain's eigenvalue. The nearest global mode is 0.1–6% away:

| `delta` | `faer` residual | Ours | Refined once | Refined twice |
|---:|---:|---:|---:|---:|
| 1e-2 | 3e-16 | 2e-16 to 6e-16 | 2e-16 | 1e-16 to 2e-16 |
| 1e-6 | 3e-16 | 5e-13 to 6e-13 | 1e-16 to 2e-16 | 1e-16 |
| 1e-10 | 3e-16 | 4e-9 to 1e-8 | 2e-16 to 6e-16 | 1e-16 to 2e-16 |
| 1e-12 | 3e-16 | 2e-7 to 9e-7 | 5e-13 to 1e-11 | 2e-16 to 5e-16 |
| 1e-14 | 3e-16 | 1e-5 to 3e-5 | 7e-10 to 8e-8 | 1e-14 to 4e-11 |
| 0 | 3e-16 | 8e-6 to 9e-5 | 6e-10 to 6e-8 | 2e-14 to 9e-11 |

So error grows as about `1e-18 / delta` without refinement, and two steps recover it, to about
1e-10 at worst when tuned within 1e-14. Reaching that takes deliberately setting a shift from a
subdomain's own eigenvalue, and a real structure would need that coincidence on one of its many
blocks. Delayed pivots would remove the risk altogether, at a real cost in complexity, and
this doesn't justify them. The stress test covers `E_z` on one mesh, so `H_z` and meshes
tuned this way remain untested, though the physics tests cover `H_z` in ordinary use.

## Tutorial changes

Chapter 9, the README and `docs/README.md` were corrected in milestone 1. A proper 2D
chapter, covering snapping and Bloch ties, is still to write.

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

For the sparse solver comparison, from `plans/probe-sparse`. Give an output directory outside
the repository. `bench` covers the native table, `bands` splits a band diagram's time into its
parts, and `stress` covers the risk above. The probe has its own copy of the solver, and reaches
`faer` through its own `faer_reference.rs`:

```bash
CARGO_TARGET_DIR=/tmp/eigenmenno-sparse cargo run --release --bin bench
```

```bash
CARGO_TARGET_DIR=/tmp/eigenmenno-sparse cargo run --release --bin bands
```

```bash
CARGO_TARGET_DIR=/tmp/eigenmenno-sparse cargo run --release --bin stress
```

```bash
CARGO_TARGET_DIR=/tmp/eigenmenno-sparse RUSTFLAGS="-C target-feature=+simd128" cargo build --release --target wasm32-unknown-unknown --lib && node run.mjs /tmp/eigenmenno-sparse/wasm32-unknown-unknown/release/probe_sparse.wasm
```

The whole-request table compares a build from before the swap, commit `7a0706a` plus milestone 4,
against the core as it is.
