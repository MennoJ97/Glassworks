# 2D solver plan

The agreed plan for adding 2D. Pick up at the first milestone that isn't marked done,
and update this file as things land. Delete it, with `plans/probe-2d/` and `plans/probe-sparse/`, once the last
milestone is finished and the tutorial covers what's here.

**Status (2026-11-01):** milestones 1 to 5 are built, apart from one optional item. Cavities,
band diagrams, open scattering, leaky modes, waveguide bends and inverse design run on the page,
tested against exact and independent results. **What's left is a proper 2D tutorial chapter,
and optionally optimizing the cavity.** The solver numbers below come from the throwaway probe in
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
| Open boundaries | PML | Exact local conditions don't exist in 2D (chapter 11) |
| Leaky modes | No non-Hermitian eigensolver. One driven solve at complex frequency for LDOS, and a secant search on `1/u_0` for Q | Avoids Krylov-Schur. Gives LDOS directly, which is what Q/V is a proxy for. It took 6 or 7 solves, and found Q up to 680,000 |
| Upper size | The page refuses above 400,000 unknowns | Measured below. 638,401 failed with `faer` LU, whose memory use was several times ours |
| Size | The 300–400 line estimate was low: milestone 1 alone added about 600 lines of Rust, for about 1,700 at that point, about 1,900 after milestone 2, about 2,000 after milestone 3, about 2,150 after milestone 4, about 2,400 with our own sparse solver, about 2,460 after milestone 5, and about 2,530 after speeding up design runs and scaling, about 2,560 after the optimized cavity | State the cost before adding more, and keep the tutorial's line count current |
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
3. **Band diagram speed: settled for now.** Warm-starting Lanczos was measured and left out. Loosening
   the tolerance took 15 to 23% off instead: see "Faster band diagrams" under "Done after milestone 6".
4. **Transposed solves: done.** Bloch-periodic driven problems lose `Aᵀ = A`, and
   `sparse/transposed.rs` solves `Aᵀ` with the same factors. See "Periodic gratings" below.

## Code changes

| Module | Change |
|---|---|
| Tutorial | Chapters 9 to 14 cover 2D |

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
within 0.3 percentage points of the test's values. With the solver adopted afterwards and the
mesh fix, two solves take about 0.85 s in Node, settled, but that hasn't been timed on the page. The first solve after loading took about 20 s, in a
hidden browser pane, and milestone 3's first solve was also about 4 times its later ones. The
cause hasn't been measured: WebAssembly compiling its hot loops is the likely one. The slab takes about 3 s. Spectra would need many
solves, so the page shows one wavelength at a time, and T can't yet be plotted against
frequency there. Download: 283 KB gzipped.

### Milestone 5: 2D inverse design

- [x] Reuse `design.js` with a 2D filter, and density per cell. Add 2D adjoint gradients
      in the solver, and keep `A` symmetric for PML problems. Checked against finite differences,
      to 7 digits. It took about 75 lines of Rust.
- [x] **Example 6: topology-optimized bend.** Compared against the hand-designed bend from
      example 5, in chapter 13. [Jensen & Sigmund 2004](https://doi.org/10.1063/1.1688450) is the
      precedent, but this is a slab bend rather than a crystal one, so it isn't a repeat of theirs.
- [x] Optimize the cavity from example 4 by raising the LDOS bandwidth's Q in stages, as in
      [Liang & Johnson 2013](https://doi.org/10.1364/OE.21.030812). See "A cavity, held on its frequency".

The example is a slab of index 3.5 and width 0.125, through a square one wavelength wide of 24 x 24
pixels. Over wavelengths 0.95, 1 and 1.05, far along the guide with a narrow monitor:

| Wavelength | Empty region | Optimized, 1/24 | Optimized, 1/48 |
|---:|---:|---:|---:|
| 0.95 | 5.4% | 97.8% | 97.6% |
| 1 | 5.3% | 99.3% | 99.3% |
| 1.05 | 5.5% | 99.4% | 99.6% |

It's a saved pattern in `tests/analytic.rs`, since the optimizer is JavaScript. So that test guards
the pattern's physics, not the optimizer, which `tests/design.test.mjs` covers (see below). The page's own run reached 98.6% on average, in 105 s at 28,000 unknowns, on battery.

That run now takes about 11 s, also on battery, and reaches 98.4%. A run is timed best on a quiet
machine: one with another process holding a full core took 51.7 s. Most of the gain isn't in the solver:

| Change | Effect | Cost |
|---|---|---|
| 10 iterations per sharpness instead of 25 | 2.5× fewer solves | one number |
| One worker per wavelength, up to three | 339 to 112 ms per iteration | about 20 lines of JS |
| Sorting triplets by row first | assembly 26 to 12 ms per solve | about 15 lines |
| Reusing the ordering when the nonzeros match, copying finished fronts by column, and searching only near the port | with the row above, 471 to 339 ms per iteration on one worker | about 30 lines |

Over three starts, T for 0.95 to 1.05 is 98.1 to 98.4% at 10 iterations per stage, and 96.9 to 99.0%
at 25, so more iterations can reach higher but aren't more dependable. Before the mesh tie-break below,
which changed these slightly, 7 per stage fell to 95.8% on one start, and skipping every other sharpness
reached only 96.5 to 97.6% even at 25 per stage. Over 0.85 to 1.15, 10 did as well as 25 (81 and 85%
against 80 and 84%), and over 0.99 to 1.01 slightly worse (99.0 and 97.2% against 99.5 and 98.3%).

Measured and left out: a handwritten simd128 kernel for `subtract` gained 6% for 22 lines of unsafe
code, because the kernel already reaches about 0.5 ns per complex update. Not tried: reusing the
fronts outside the design region between iterations, which needs a solver that keeps state per wavelength.

### A cavity, held on its frequency

The goal is LDOS relative to vacuum, from a line source at the centre of a square of 64 x 64 pixels,
2.67 wavelengths wide and with eps 12, as big as the 7 x 7 rods it's compared with. It cost about 25 lines of
Rust: LDOS gradients with no extra solve, the complex response, and the gradient's imaginary part.
The optimizer is about 180 lines of JavaScript.

It failed at first, in three separate ways, in the order they were found:

1. **The resonance drifted.** Peaks ended 2 to 5% off target, where a bandwidth of 0.1% can't see them.
2. **The filter was too wide.** A radius of 3 pixels at 24 per wavelength left 5.4% of a crystal's
   12.3% fill, so the design space couldn't hold a crystal at all. Now it's 2 pixels.
3. **Pixels didn't line up with cells.** The optimizer exploited that, and a binary design's mode moved
   2% when drawn as rectangles. The solver now refuses a misaligned density map.

Following the resonance around helped less than holding it still, which was the idea that worked:

- Near a pole, `1 / u_0` is nearly linear in `k`, so two solves place the pole. The gradient's real and
  imaginary parts together give each density's effect on it, checked within 2% on a true resonance.
- Each step has its slope along the frequency removed, and a correction then aims half of what's left
  of the error back at the target.
- Sharpening still leaves up to about 2%, because a pixel flip is a big step. Resizing the finished
  design takes it out exactly, by scale invariance.

The hold's effect, over nine starts on the page's setup. Three starts leave more to chance, so only the
first two rows settle anything. LDOS at Q 1,000 hides differences, because it saturates near 1,060 once a
mode's own Q is far above that bandwidth, which is why this also reads it at Q 10,000:

| Version | Starts | LDOS at Q 1,000 | Median at Q 10,000 | Median mode Q | Largest resizing |
|---|---:|---|---:|---:|---:|
| With the hold | 9 | 110 to 1,061, median 917 | 8,411 | 56,000 | 1.2% |
| Without projection or correction | 9 | 114 to 888, median 765 | 3,304 | 5,900 | 4.4% |
| Without the projection | 3 | 355 to 954 | | | |
| Without the correction | 3 | 61 to 485 | | | |
| With no goal, only holding | 9 | 27 to 372, median 64 | | | |
| With the goal's sign flipped | 9 | 904 to 1,033, median 990 | 9,483 | 416,000 | 3.1% |

An earlier claim here, of 55 against 1,059, came from one start on an older setup, with the projection
removed. The commit that added this view repeats it, and it overstated the effect.

LDOS at Q 1,000, after the same scan and resizing for each:

| Structure | LDOS | Q | V (wavelengths²) |
|---|---:|---:|---:|
| 16 random symmetric patterns | 0 to 29 | | |
| Crystal of 7 x 7 rods | 299 | 1,680 | 0.093 |
| 9 optimized starts | 110 to 1,061, median 917 | 253 to 449,000 | 0.008 to 0.017 |

Eight of nine beat the crystal by at least 1.5 times, and the page's default start is the best. As
rectangles, its LDOS and Q hold to under 1% from cells of 1/24 to 1/72 of a wavelength. On a mesh of
1/48, its Q is 328,000 at layers of 0.5, 1 and 1.5 wavelengths, and on 1/24, 329,000 with 0.25 to 1.5
wavelengths of air around it. Its LDOS at Q 3,000 is 3,162, against the crystal's 508. Stopping Q at
500 while sharpening was the steadiest schedule: raising it to 2,000 reached 269 on one start.

Scatter between starts is larger than these choices, so each was compared over nine:

| Stages, layer | Solves | Cost of one | LDOS at Q 1,000 | Median |
|---|---:|---:|---|---:|
| 11, 0.5 | 471 | 117 ms | 300 to 1,059 | 884 |
| 8, 0.5 | 351 | 121 ms | 303 to 1,059 | 990 |
| 8, 0.25 (the page's) | 351 | 87 ms | 110 to 1,061 | 917 |

So the page's choice costs about 45% less for the same median, but has the lowest worst start.
Fifteen iterations per stage instead of 20 reached only 329 to 844 on three starts.

Three cautions. Q in the hundreds of thousands from a footprint this small needs a cause, which we haven't
found, and no fabrication tolerance has been checked: a mode this sharp would be very sensitive. The
table above also made minimizing look better than maximizing, which the next section corrects.

### Why the flipped goal looked better

Mostly it didn't: nine starts misled us. Over 18 starts, with 10 to 18 unseen by any earlier choice:

| Version | LDOS at Q 1,000: worst, median, best | Below 900 | Median mode Q | Below Q 10,000 | Beats the page on |
|---|---|---:|---:|---:|---:|
| The page: maximize, correct half the error | 110, 985, 1,061 | 6 | 113,000 | 6 | |
| Flipped goal | 708, 956, 1,033 | 3 | 172,000 | 3 | 7 of 18 |
| Correct 0.75 of the error | 103, 1,002, 1,061 | 5 | 97,000 | 6 | 11 of 18 |
| Correct all of the error | 342, 1,016, 1,064 | 3 | 268,000 | 3 | 13 of 18 |
| LDOS read at the mode's own wavelength | 814, 1,014, 1,075 | 1 | 374,000 | 0 | 14 of 18 |

The first nine alone gave flipped 990 against 917; the last nine give 951 against 1,001. So its lead in
the median was chance. It does fail less often, and the wing of the resonance is a likely reason.

Near a pole, with detuning D, bandwidth G and the mode's own decay g, LDOS goes as
A (G + g) / (D² + (G + g)²). Its slope along g changes sign at |D| = G + g. Inside that, raising LDOS
raises Q. Outside it, on the wing, raising LDOS lowers Q, and lowering LDOS raises it.

The hold removes the frequency part of the goal, which is nearly all of it: what's left has a median share
of 0.5% in the page's runs, and that remainder sets the step. It behaves as predicted. It points towards
more decay on 0.4% of iterations below half a linewidth, 15% from a half to one, 97% from one to three and
99% beyond. The flipped goal mirrors that, and spends 99% of its iterations outside one linewidth. Across the
18 page runs, more iterations past one linewidth goes with a weaker result (Spearman -0.76, from Q 100 on).

That's an association plus a confirmed sign, and these limits apply:

- Q doesn't always follow on the very next iteration, because the correction and the sharpening also move
  it. With full correction, Q rises on 61% of its wing iterations, opposite to the sign above.
- Full correction's failures (starts 12, 16 and 18) lose their Q at the last sharpening stages. Start 16
  starts falling while still at 0.05 linewidths, so the wing can't be the whole cause.
- Flipped with full correction tests nothing: the correction is clamped far from the target, and eight of nine
  starts reproduce the plain flipped result exactly.
- Maximizing, then flipping (median 621), or flipping, then maximizing (803), over nine starts, both did
  worse than either sign throughout. So building a blocking surround first isn't what helps. Without the hold, the
  flipped goal drives LDOS to zero and no leaky mode remains.
- Lowering decay as the only goal wasn't a usable comparison. As coded it raised decay, and Q collapsed to about 60.
- Earlier tables in this file, and chapter 13 until it's corrected, describe the half-correction method read at the target.

The page now reads LDOS at the mode's own wavelength, which needs the same solves and keeps the true goal.
Its default start has LDOS 1,066 at Q 1,000, 3,193 at Q 3,000, Q of 1,512,000 and V of 0.0079. Drawn as
rectangles, starts 1, 3, 5 and 16 change by about 1% or less in LDOS and under 2.7% in Q from 1/24 to 1/72. For the default, Q stays
within 0.3% for layers of 0.5 to 1.5 wavelengths, and within 0.7% with 0.25 to 1.5 of air. No fabrication tolerance has been checked,
so these Q values are numerical results, not achievable ones, and they still lack a physical explanation.

### Symmetric factors for open domains

Factorization took 77% of a native cavity run. Other block and leaf sizes were slower, and an explicit
SIMD kernel gained nothing, so the compiler was already vectorizing well. An absorbing layer's operator
is complex symmetric, which `factor_symmetric` uses: `L D L^T` with no pivoting, for 72 lines of Rust.

| Run, default start | Pivoted LU | Symmetric | Saved |
|---|---:|---:|---:|
| Cavity on the page, visible | 16.0 s | 10.95 s | 32% |
| Cavity, one thread, WebAssembly | 28.35 s | 18.92 s | 33% |
| Bend, one thread, WebAssembly | 18.4 s | 13.2 s | 28% |

Across 8,612 line-source solves at real `k`, scanned through resonances, on coarse uniform meshes where a
diagonal can land near zero, and with permittivities up to 1,000, the worst residual was 1.1e-11, against
6.4e-12 pivoted, and solutions differed by at most 1.2e-11. Those ran before the fallback was added.
Without pivoting, weak pivots are about 3 times as common overall, so more solves are refined. A pivot under
1e-8 of its column, or 1e-14 of the matrix, now falls back to `factor`. Ports and plane waves are covered
only by the test suite, not by this scan.

Of 18 starts, 15 give the same design to the digits printed. Three change through round-off, because the
optimizer amplifies differences near 1e-12: 5 from 814 to 886 at Q 1,000, with Q 13,484 to 9,120, 2 from
1,011 to 996, and 12 from 1,017 to 1,011. The medians hold at 1,014 and 374,000, but one start now
finishes below Q 10,000. The table above is from before this change.

One thread runs all of this in 28 s in Node, on battery. The page uses two workers, and paints only a few
times a second. A hidden pane ran the older, longer schedule far slower, at 229 s, so timing it needs a visible page.

### A finer mesh, and the crystal at 1 µm

Peak memory in WebAssembly, where one instance can have at most 4 GiB:

| Unknowns | Driven | Time, driven | Eigenmodes | Time, eigenmodes |
|---:|---:|---:|---:|---:|
| 51,529 | 0.12 GiB | 1.0 s | | |
| 143,641 | 0.34 GiB | 3.7 s | | |
| 160,801 | | | 0.7 GiB | 14 s |
| 314,721 | | | 1.7 GiB | 31 s |
| 337,561 | 0.82 GiB | 10 s | | |
| 576,081 | 1.51 GiB | 24 s | | |
| 579,121 | | | 3.2 GiB | 62 s |
| 900,601 | 2.46 GiB | 39 s | | |

Those are on battery. A limit of 400,000 leaves eigenmodes near 2.1 GiB, in case of fragmentation
and because a page can run up to three workers. A phone may still run out; the page then reports
that the solver stopped. A finer mesh is only used when it's asked for or a short wavelength needs it.

The original crystal at 1 µm, solved natively with the limit removed:

| Cell | Unknowns | Native peak | Straight, share of source | Bend, share of source | Ratio |
|---:|---:|---:|---:|---:|---:|
| 0.0654 | 337,561 | 1.5 GiB | 19.6% | 11.6% | 58.1% |
| 0.05 | 576,081 | 2.7 GiB | 20.1% | 10.2% | 47.1% |
| 0.04 | 900,601 | 3.9 GiB | 19.8% | 11.9% | 57.8% |
| 0.032 | 1,408,969 | 5.9 GiB | 19.1% | 11.0% | 56.6% |
| 0.025 | 2,307,361 | 6.7 GiB | 18.9% | 11.0% | 57.4% |

So about 11% of the source's power reaches the bend's output, and the page's 11.6% is close. It isn't
guided: the crystal's complete gaps are at $a/\lambda$ of 0.304 to 0.446 and 0.741 to 0.766 on a 1/24
cell, and bands overlap at 1. The second gap, at 1.31 to 1.35 µm here, does guide light, with 36 to 48%
reaching the straight guide's output on two meshes, but the sharp bend passes 76% at 1.31, then 8.5%
and 4% at 1.33 and 1.35.

### Tests of the optimizer

`npm test` runs `tests/design.test.mjs`, in 34 s, against the page's WebAssembly. It checks each filter's
adjoint, the projection's slope, that MMA finds a bowl's minimum within its bounds and move limit, the
whole chain's gradient against finite differences for all three goals, that whole runs end sharp and
beat their references, and that a band at twice the size gets the same gradient. Breaking `design.js`
seven ways was caught each time: either gradient's sign, the projection's slope, the filter's
adjoint, MMA's move limit, the final sharpening, and two iterations per stage. A bend's pass mark is
97%, against 98.1 to 98.4% over three starts.

### Bends at any wavelength

The straight-guide reference was only ever valid where the guide carries light, and the page showed
465% at the crystal's default of 1 µm. Its straight guide's output gets about half of what the source
emits when it works (41 to 50% from 2.3 to 3.2 periods), and 16% or less outside that, so:

| Situation | The page shows |
|---|---|
| Fewer than 4 cells per wavelength in the densest material, even at the solver's limit | no percentage, and an offer to resize |
| Straight guide's output outside 40 to 60% of what's emitted | the bend's share of what's emitted, and that offer |
| Ratio above 102% | "Unreliable", and that offer. From 100 to 102% it shows 100%, which is the mesh's error |
| Otherwise | the ratio |

The mesh aims for 4.5 cells per wavelength in the densest material, as far as the limit on unknowns
allows, and never coarsens what was asked for. It doesn't change the Mekis numbers at 2.833 µm, but
does refine below about 2.55 µm: at 2.4 µm the sharp bend's T goes from 95.6% to 91.9%.

Resizing is exact: the sharp bend reads 97.7544% at 0.55, 1, 2.833 and 10 µm. A test holds this to 1e-7,
and a density map's gradient to 1e-11. Getting there took
four fixes, because whole numbers of cells and points meant rounding noise was choosing: the count of
steps along a grid, monitor or port, a corner halfway between two nodes, the order of nodes at equal
distances from a shape, and the material of a cell whose centroid lies exactly on a shape's edge. That
last one happens where a slab meets the outer wall, whose nodes can't move. The slab otherwise fills
the layer all the way there, 6 of 6 cells in each of 12 columns, and the wall's own column has 5 of 7.
A straight guide's power reads 0.005735 at layers of 0.5, 1 and 1.5, and 0.004354 if the slab is
deliberately stopped where the layer starts. The tie-break moved the slab bend's T only from 22.770% to 22.767%.
The optimized bend scales with its band's centre, so 0.475 to 0.525, 0.95 to 1.05 and 1.9 to 2.1 µm
all finish at 98.37%. Its mesh can double per pixel, to serve a band's short end. A band reaching below
0.84 of its centre gets a warning, because the slab carries a second mode there. The original crystal at 1 µm is answered below.

A first run monitored a wider line near the corner, which read 98.8% against 97.7% far along the
guide, because it counted radiated light. A narrower monitor closed that gap, to under 1%. Other
findings: a single-wavelength optimization reached only 43% and 82% at the band's edges, and a
square half a wavelength wide reached 96% at the centre, but 81% and 92% at the edges.

Against the crystal bend, this one is smaller and needs no band gap, while the crystal's T stays
above 95% over about 13% of bandwidth against the 10% tried here. How wide this design could go
wasn't measured: a run over 0.9 to 1.1 with a 1.5-wavelength square reached 99.6% while optimizing,
but crashed on its finer-mesh check, for exceeding the mesh limit, so it isn't verified.

Only `E_z` has gradients, since `q = eps` is simple there. `H_z` would need `p = 1/eps` in the
stiffness. The adjoint here needs `A^T = A`, so a Bloch-periodic structure would need open problem 4.

### Milestone 6: would it work if you built it?

The page reports Q values in the millions from a very small footprint, with no physical explanation
and no check against fabrication errors. This milestone answers both, in that order, because
step 1 decides how much of the rest is worth doing. Nothing goes on the website until its results have
been seen and agreed.

- [x] **Measure tolerance, with no solver changes.** Grow and shrink the finished rectangles of nine
      starts, the crystal and the bend, by 1, 2, 5 and 10 nm at a wavelength of 1 µm, and read Q, LDOS and the
      mode's wavelength. Pixels are 42 nm there, so these moves are sub-pixel and must be drawn: grow by
      enlarging each material rectangle, and shrink by drawing the square in material and enlarging
      each air rectangle over it. The mesh has to resolve the move, so run this natively, beyond the page's limit if needed,
      and check convergence with the mesh, as the rest of this plan does.
- [x] **Explain the high Q.** Look at the mode's field, where it radiates, and how Q changes with the
      surrounding air and with the size of the design region. Compare the best starts with those that stay
      below Q 10,000. Finding no cause is an acceptable result, as long as it's stated.
- [x] **Decide what robustness should mean.** A shifted edge moves the resonance, and at a bandwidth
      of Q 1,000, likely by more than its width. So one wavelength may be unreachable for all versions of a
      design. Read each version at its own wavelength, and report how far they spread as a separate
      result, which is what a real device would have to tune away. Check this against step 1's numbers.
- [x] **Optimize for robustness, if step 1 shows a need.** Prototype only, not on the page. Project the filtered densities at three
      thresholds, for eroded, nominal and dilated designs ([Wang, Lazarov & Sigmund 2011], to be
      checked before citing), and raise the weakest. This belongs in `design.js`, with no new Rust expected;
      state the line cost if that changes. A threshold's shift in pixels sets the edge move, so measure that
      against step 1. It takes six solves per iteration instead of two, so about three times as
      much computing, for a run that now takes 11 s on the page. Nine starts aren't enough, as this
      plan learned: use 18 at least, with unseen ones for any final comparison.
- [x] **Teach it.** In chapter 17, where the cavity moved, with a note on the bend in chapter 13. First planned as a section in chapter 13, with a robust design beside the sharpest one, if there's a
      difference to show. A negative result is worth teaching too.

Done means the page can say what happens to its designs under a stated edge error, and either offers a
robust one or explains why not.

The check against a known answer: first-order perturbation theory, from the pole's slope with each density,
should predict the resonance shift of a small drawn edge move. Add that as a test.

#### Step 1 results: tolerance, measured 2026-11-02

Everything here is at a wavelength of 1 µm, with `eps` 12 and lossless, in 2D. It's scratch work, and none is
on the website. The harness isn't in the repo yet.

**How to move an edge by less than a cell.** Cells are 42 nm at the page's resolution, so 1 nm is 1/40 of
one. What was tried:

| Method | Result |
|---|---|
| Fractional density in the edge's cell | Rejected: its shift is off by 15 to 39%, and differs between growing and shrinking |
| Whole runs of pixels, grown or shrunk | Good for uniform moves, but can't move edges independently |
| Thin patches on each edge | Silently does nothing: `snap` places each node once, and the first shape takes it |
| Each pixel its own rectangle, with exposed sides moved and others joined | Used: nodes move onto the true edge, and random offsets work |

The last was checked on the default start. Unmoved, it gives Q 1,518,579 against 1,518,599. From 1/48 to 1/72,
the wavelength at ±1 nm agrees within 0.03% and Q within 2%. Boundary perturbation theory, from the field alone,
predicts a shift of 1.35% per nm, and this measures 1.33%. Random offsets, for the same draws on both meshes, mostly
agree within 10% in Q, and at worst about 30%. Q at ±5 nm isn't converged (376,000 and 634,000 at +5 nm), so
it's not quoted. At whole cells, moved runs reproduce an exact grid of 1/96 to every printed digit, on one mesh.

**A uniform edge error moves the resonance far outside its bandwidth.**

| Structure | Shift per nm | LDOS at the original wavelength, Q 1,000 |
|---|---:|---|
| 9 optimized starts | 1.30 to 1.34% | About 1,000 unmoved, 6 at ±0.5 nm, 1.5 at ±1 nm |
| Crystal | 0.32% | 289 at its own wavelength, 18 at ±1 nm, 5 at ±2 nm |

The bandwidth is 0.1%. So neither would work at a fixed wavelength without tuning, and the optimized
designs are about 4 times as sensitive. That confirms the concern below about one wavelength for all versions.

**Random disorder.** Each exposed pixel side moves independently, from a Gaussian, which makes the correlation
length one pixel, 42 nm. Nine starts, six draws each, on a mesh of 1/48, with medians across starts of each
start's median, and its range of those medians:

| Sigma | Mode Q | Spread of wavelength | LDOS at own wavelength | At the original wavelength |
|---:|---:|---:|---:|---:|
| 0 | 533,000 (9,300 to 5.9 million) | | 1,006 | 1,006 |
| 0.5 nm | 370,000 (10,700 to 1.4 million) | 0.14% | 1,006 | 171 |
| 1 nm | 254,000 (8,100 to 816,000) | 0.32% | 1,005 | 30 |
| 2 nm | 118,000 (12,000 to 259,000) | 0.62% | 992 | 34 |
| 5 nm | 26,000 (2,800 to 56,000) | 1.06% | 949 (lowest 615) | 3 |

So:

- **Q in the millions doesn't survive.** At 5 nm, the median falls to about 26,000. It's still about 16 times the crystal's 1,552,
  whose own Q hardly moves. The start with the highest Q loses the most, so ranking starts by unperturbed Q misleads.
- **LDOS holds, if the wavelength can be retuned.** It depends mainly on the small mode volume, and stays at about 3 times the crystal's 270 at
  5 nm. Without retuning, it's lost at under 1 nm.
- **Random disorder costs more Q than the same uniform error.** Up to ±2 nm, a uniform move changes Q by 0.16 to 2.1 times, and
  some designs improve. A uniform move keeps the eight-fold symmetry, and random ones don't. That suggests the
  high Q relies on symmetry, but this doesn't test it.

**The bend is far more tolerant.** Three starts, mean T over its band, with arms fixed, and about 1% mesh error:

| Error | T, starting from 98% |
|---|---:|
| Random, 1 nm | 97.8 to 98.3% |
| Random, 5 nm | 94.7 to 98.7% |
| Random, 10 nm | 84.9 to 96.8% |
| Uniform, ±5 nm | 96.7 to 97.6% |
| Uniform, -20 nm and +20 nm | 58 to 79%, and 85 to 91% |

Moving the arms as well gave a kink of about 2% at ±1 nm, likely because the port's 1D mode samples the guide in steps of 5 nm, far
coarser than the move. So those runs aren't used, and a test that moves a port's guide by a sub-cell amount shouldn't be trusted.

Limits on all of this:

- **It's 2D.** A real slab also radiates out of the plane, and materials absorb. Either could cap Q well below these values.
- **The disorder model is assumed.** Sigma and correlation length here are chosen, and not from a measured process. The crystal's random
  case changes each rod's radius, a different model from the designs', so its comparison is rough.
- **Few samples.** Six draws per start, and three on the finer mesh. Medians are indicative, and worst cases are understated.
- **Feature size isn't checked.** Designs may hold features of about 42 nm that couldn't be made at all.

**What this changes for the next steps.** Three thresholds, as first planned, test uniform erosion and dilation, which
is the milder of the two for Q, and a uniform error's main effect, the shift, can only be tuned away. So a robust goal
should sample random edges too, or aim at what they do to Q, which step 2 needs to explain first.

#### Step 2 results: why Q gets so high

**The Q is real, at least within this 2D model.** On a circle in the air, the field splits into
`u = sum a_m H_m(k r) exp(i m t)`, whose outgoing power gives `Q = k^2 integral(eps |u|^2) / (4 sum |a_m|^2)`. That
uses only the field, and shares nothing with the search for the pole. It agrees within 0.5% unmoved, including at 5.9
million, and within about 6% under disorder. The two circles tried give the same answer.

**Only a few channels radiate, which is what makes this possible.** Eight-fold symmetry allows only `m = 0, 4, 8, 12, ...`.
A footprint this small barely couples to higher `m`. Unmoved, the other channels carry under 0.1%, and `m >= 16` under 0.5%.
So there are about four channels to quieten, with a few hundred free pixels.

**The low channels are quiet by interference.** Each `a_m` is also an integral of the material's currents, so comparing its
size with the sum of the sizes of its parts measures cancellation. Using the circle's `a_m`, for the default start and the crystal:

| Channel | Default start | Crystal |
|---:|---:|---:|
| `m = 0` | 1.7e-4 | 5.7e-2 |
| `m = 4` | 8.4e-4 | 3.1e-2 |
| `m = 8` | 2.5e-3 | 2.5e-2 |
| `m = 12` | 9.2e-2 | 4.2e-2 |

So cancellation in the lowest channels is about 10 to 340 times stronger than the crystal's, while `m = 12` is just weak, with no
special cancellation. The integral's own value at `m = 0` is 2 to 5 times off the circle's for high-Q designs, because so much
cancels. The others agree within a few percent, and the crystal's `m = 0` within 5%.

**It's also confinement.** Relative to the peak, the field at a radius of 1.2 wavelengths is 1.7e-3, against the crystal's 1.1e-2, and
the outgoing field at 2.0 is 2.2e-4, against 8.3e-3. Those are about 38 times apart, or 1,400 in power, near the ratio of Q.
Size helps too, over nine starts each:

| Side, in wavelengths | Worst Q | Median Q | Best Q | Median LDOS | Median V |
|---:|---:|---:|---:|---:|---:|
| 1.67 | 7,800 | 38,000 | 93,000 | 986 | 0.0083 |
| 2.00 | 461 | 47,000 | 310,000 | 994 | 0.0082 |
| 2.67 (the page's) | 9,100 | 514,000 | 5.9 million | 1,015 | 0.0083 |
| 3.33 | 688,000 | 4.6 million | 56 million | 1,009 | 0.0084 |
| 4.00 | 377 | 1.9 million | 305 million | 1,021 | 0.0083 |

Q grows quickly with size, up to 3.33 at least, but V and LDOS don't change at all. So the mode's size and the Q are nearly separate questions:
the small region near the source sets V, and the rest controls leakage. Q above about 10 million hasn't been checked for numerical limits.

**Why disorder hurts.** Random edges break the symmetry, and the forbidden channels then carry 29 to 91% of the power. The test of
equal offsets on all eight images is less clear. It protected start 3's Q about 40 times at 2 nm, but barely helped starts 6 and 8, because those offsets also
disturb the cancellation in the allowed channels, and shift the wavelength further. Those runs can show forbidden power too, since the mesh isn't quite
symmetric once nodes move, which matters when Q is in the millions. So symmetry is part of the answer, not all of it.

Which leaves, for step 3 and after:

- **Q beyond the bandwidth buys little.** LDOS at Q 1,000 is close to its cap once a mode's Q passes about 10,000, and that survives 5 nm for most starts.
- **What disorder really costs is the wavelength,** and no robust goal removes a uniform shift. Across these designs, shift per nm barely varies. It's 1.3%
  at V 0.008 and 0.32% at the crystal's V 0.09, so a small mode may simply be a sensitive one. That's untested beyond these two groups.

#### Step 3 results: a trade-off, and a cheap way to take it

The optimizer always settles at the same mode, so these force others: an air disc around the source, a solid disc, or wider
features, with six starts each at first. `V` isn't comparable between them, because it's measured against the peak of `q |u|^2`, which
depends on what the source sits in. LDOS is the fair measure of how tightly a mode is held at the source. Medians:

| Forced | LDOS, Q 1,000 | Shift per nm | Mode Q |
|---|---:|---:|---:|
| Nothing, 15 starts | 1,005 | 1.30% | 403,000 |
| Air disc, radius 3 to 5 pixels | 154 to 325 | 0.38 to 1.50%, erratic | 1,300 to 4,200 |
| **Air disc, radius 6 pixels (250 nm), 15 starts** | **611** | **0.26%** | **215,000** |
| Air disc, radius 7 and 8 | 598, 589 | 0.27%, 0.26% | 110,000, 142,000 |
| Air disc, radius 10 | 543 | 0.48% | 108,000 |
| Air disc, radius 20 | 300 | 0.14% | 38,000 |
| Solid disc, radius 3 | 352 | 0.56% | 18,000 |
| Solid disc, radius 6, 10, 15, 20 | 245, 174, 134, 85 | 0.32, 0.24, 0.18, 0.14% | 814,000, 790,000, 254,000, 19,000 |
| Features, filter radius 3 | 1,036 | 1.33% | 259,000 |
| Features, filter radius 4, 6, 8 | 92, 190, 27 | 0.65, 0.23, 0.19% | 369, 2,600, 520 |
| Crystal, for comparison | 289 | 0.32% | 1,552 |

So a trade-off exists, but it isn't one law: at the same LDOS, an air disc is far less sensitive than a solid one, and from 1,005 to 611 costs 39% of
LDOS for 5 times less shift. A clear disc of 6 pixels beats the crystal on LDOS, shift and Q at once. Smaller discs are unreliable, for a reason we haven't
looked into. It held on nine unseen starts: LDOS 570 to 621 over all 15, shift 0.26 to 0.32%, Q 22,000 to 2.4 million, and none below Q 10,000.

Random disorder, nine starts, six draws each, against the unconstrained designs from step 1:

| Sigma | Spread of wavelength | LDOS at the original wavelength | LDOS, retuned | Mode Q |
|---:|---|---|---|---|
| 1 nm | 0.022% against 0.32% | **536** against 30, crystal 132 | 612 against 1,005 | 166,000 against 254,000 |
| 2 nm | 0.045% against 0.62% | **420** against 34, crystal 62 | 599 against 992 | 55,000 against 118,000 |
| 5 nm | 0.107% against 1.06% | **154** against 3, crystal 23 | 492 against 949 | 7,500 against 26,000 |

- **At a fixed wavelength, it survives random disorder.** Its spread stays inside the 0.1% bandwidth up to about 5 nm, and the unconstrained designs miss it by 0.5 nm. The spread improves about 10 to 15 times, more than the 5 for
  a uniform shift. Field spread around a ring averaging over many independent edges would explain that, but it hasn't been tested.
- **If the wavelength can be retuned, the unconstrained design is still better,** with nearly twice the LDOS and a higher Q.
- **A uniform bias still needs tuning.** At ±1 nm, LDOS at the original wavelength is about 20, like the crystal's 18, against 1.5.

Checks. Theory predicts 0.287% per nm from the field, on all three meshes, and 1/72 measures 0.270%, within 6%. Between 1/48 and 1/72, wavelengths at ±2 nm agree within 0.03%, Q within 16%
at worst and usually 7%, and random draws agree in shift within 0.01% absolute. At ±5 nm and beyond, wavelengths still agree within 0.4%, but Q doesn't always converge.

**Where the source sits.** A fair question is why air would work, since the field normally concentrates in high index. Left alone, it does: every unconstrained
design checked puts the source in a rod of 2 x 2 pixels, about 83 nm across, with 87% of the energy in material. LDOS follows `|E_0|^2 / integral(eps |E|^2)`:

| Design | Source in | Energy in material | Peak per energy | Nearest edge | Edge field per energy | Shift per nm |
|---|---|---:|---:|---:|---:|---:|
| Unconstrained | Rod | 87% | 10.5 | 47 nm | 2.50 | 1.30% |
| Solid disc, 3 pixels | Material | 95% | 3.9 | 178 nm | 1.14 | 0.56% |
| Solid disc, 6 pixels | Material | 97% | 2.4 | 326 nm | 0.65 | 0.32% |
| Air disc, 6 pixels | Air | 22% | 6.1 | 369 nm | 0.52 | 0.26% |

So the rod wins on enhancement, and the air disc is a compromise, not a better design. Its broad peak is cheap because energy in air counts 12 times less. The rod's edge is
47 nm from the peak, and the ratio of edge field, 4.8, matches the ratio of shift, 5. This is `E_z` only, which is continuous across edges. The other polarization has the slot effect and may behave differently. In 2D, a
line source in uniform material reads LDOS 1.00 at any `eps`, so comparing sources in air and in material is fair here. A 3D emitter's baseline depends on its surroundings.

#### Step 4 results: optimizing for robustness

Each iteration adds three samples, with edges moved by fractions of a pixel, a Gaussian of 2 nm. The step follows their mean gradient, and the frequency hold stays on the unmoved design.
All of it is judged afterwards by drawn edges, not by this model of them. In the scratchpad only, as `plans/probe-tolerance/robust-goal.patch`. On the first nine starts:

| Goal | LDOS | Shift per nm | Mode Q, unmoved | Mode Q at 5 nm |
|---|---:|---:|---:|---:|
| None | 1,006 | 1.30% | 533,000 | 26,000 |
| Samples at the target wavelength, 2 nm | 694 | 1.23% | 2,400 | 1,000 |
| Samples at the target wavelength, 5 nm | 286 | 1.26% | 455 | 477 |
| Eroded and dilated by 2 nm, at the target wavelength | 250 | 1.10% | 1,600 | 705 |
| **Samples at their own wavelength** | 1,029 | 1.29% | 1.1 million | 39,000 |
| The same, with a bandwidth 10 times sharper | 875 | 1.23% | 37,000 | 11,000 |

Samples at the target sit on their resonance's wing, where raising LDOS lowers Q, as found above, so Q collapses. The three-threshold method first planned does the same. Reading
each at its own wavelength, predicted from `moves`, avoids that. Eight variants, with and without the disc, were tried on those nine starts, so the two that looked best went to 18 unseen starts, each paired with a plain run from the same start:

| | No disc | No disc, robust | Air disc | Air disc, robust |
|---|---:|---:|---:|---:|
| LDOS | 1,005 | 1,035 | 599 | 600 |
| Lowest LDOS | 989 | 997 | 452 | 564 |
| Shift per nm | 1.31% | 1.30% | 0.27% | 0.26% |
| Mode Q, unmoved | 393,000 | 504,000 | 306,000 | 417,000 |
| Mode Q, 1 nm | 216,000 | 437,000 | 87,000 | 276,000 |
| Mode Q, 2 nm | 116,000 | 269,000 | 43,000 | 137,000 |
| Mode Q, 5 nm | 22,000 | 43,000 | 12,000 | 40,000 |
| LDOS retuned, 5 nm | 928 | 943 | 530 | 566 |
| LDOS at the original wavelength, 1 nm | 38 | 42 | 525 | 551 |
| LDOS at the original wavelength, 5 nm | 2 | 2 | 153 | 197 |

Paired, the air disc's Q improves 1.9, 2.9 and 4.4 times at 1, 2 and 5 nm, on 12, 14 and 14 of 18 starts, with unmoved Q about level. Without the disc, Q doubles or so, but it also does
when unmoved, so that looks more like a steadier optimization than robustness. Only the 14 of 18 are beyond what chance would readily give, at about 3% by a sign test, before allowing for how many things are compared. It doesn't remove the shift, as expected. At a fixed wavelength it
changes little (at 2 nm, 0.96 times with the disc). So robustness helps Q and worst cases, and doesn't replace tuning. It takes five solves per iteration instead of two, about twice the time on one thread (77 s against 41 s).


#### Feature sizes

For each pixel, the largest square of its own phase that holds it gives a local feature size. Pixels are 41.7 nm at a wavelength of 1 µm. Medians over nine starts, for what's thinner than each width:

| Design | Below 83 nm, material and air | Below 125 nm | Below 167 nm |
|---|---:|---:|---:|
| Rod | 5%, 3% | 59%, 40% | 94%, 59% |
| Air disc | 7%, 2% | 66%, 17% | 87%, 50% |
| Rod, robust | 3%, 4% | 61%, 55% | 84%, 83% |
| Air disc, robust | 3%, 3% | 52%, 21% | 81%, 62% |

Removing those features from finished designs, by opening and closing with a square. LDOS and Q after a search within 0.8 to 1.2 wavelengths, since the mode moves. Later columns aren't all the same mode:

| Design | Nothing | Below 83 nm | Below 125 nm |
|---|---:|---:|---:|
| Rod | 1,006, Q 533,000 | 943, Q 22,000 | 16, Q 136 |
| Air disc | 614, Q 296,000 | 601, Q 56,000 | 7, Q 79 |
| Rod, robust | 1,034, Q 536,000 | 1,002, Q 287,000 | 6, Q 44 |
| Air disc, robust | 599, Q 215,000 | 563, Q 21,000 | 9, Q 80 |

Designing with a wider filter from the start, six starts each, which removes most thin parts but sets no strict limit, since isolated thin parts remain:

| Design | LDOS: worst, median, best | Median Q | Below Q 10,000 | Shift per nm | Material below 125 nm |
|---|---:|---:|---:|---:|---:|
| No disc, filter 3 | 819, 1,036, 1,069 | 259,000 | 1 | 1.33% | 56% |
| No disc, filter 4 | 6, 92, 370 | 369 | 5 | 0.65% | 4% |
| Air disc, filter 3 | 579, 607, 611 | 643,000 | 0 | 0.26% | 37% |
| Air disc, filter 4 | 79, 279, 517 | 7,700 | 3 | 0.27% | 6% |
| Air disc, filter 6 | 8, 102, 345 | 1,600 | 5 | 0.24% | 1% |
| Air disc 417 nm, filter 4 | 374, 492, 540 | 18,000 | 3 | 0.44% | 16% |
| Air disc 417 nm, filter 6 | 128, 198, 316 | 2,100 | 5 | 0.37% | 6% |

For comparison, the crystal's rods are 151 nm across, with LDOS 289. Bends over three starts, from 98%: removing what's below 83 nm leaves 95 to 98%, and below 125 nm, 5 to 87%. Filter sizes are in pixels. Six starts per row is few, and a poor result can be the optimizer's as much as the physics'. A strict minimum feature size hasn't been tried.

The merge of the other session's work into this plan moved the cavity to chapter 17 and renamed `gradient.rs` and `resonance.rs`, and those names above are historical.

### After milestone 6, not yet ordered

None of these is planned in detail.

- **Designing gratings on the page.** The solver returns each order's gradient, so this is JavaScript only.
- **Vector fields.** Full-vector waveguide cross-sections are the first real need for Nédélec elements. See
  "Vector fields: a line estimate" below, which recommends not starting yet.

### Done after milestone 6

**Tutorial clean-up.** Chapters 9 to 14 were read line by line, and the items flagged when they were written were tested.

- **Joannopoulos, for band extremes on the zone's edge:** the book's text couldn't be reached, so chapter 9 now cites
  Harrison et al. 2007, whose abstract says this directly.
- **Davis 2006, for storage cost:** also unreachable, and not needed. George 1973's abstract states $O(n^3)$ operations and
  $O(n^2 \log n)$ storage, and that no ordering beats $O(n^3)$. Davis 2006 stays in chapter 4, for compressed rows.
- **Why eigenmode solves use the most memory: holds.** Counting allocations natively at 579,121 unknowns, for ten modes:
  0.63 GiB for mesh and matrices, 1.18 GiB held by the factors, 1.85 GiB at the peak of factoring, then 0.97 GiB more for
  Lanczos, whose 56 vectors each store $B$ times themselves. That ends at 2.78 GiB, or 2.50 GiB for four modes. At 160,801
  unknowns it reads 0.72 GiB, matching the 0.7 GiB in WebAssembly.
- **Why the plain corner reads 7.1% on the page and 5.4% in the table: holds.** On the same mesh, moving only the monitor from
  2.1 to 2.85 takes it from 7.10% to 5.40%, and to the table's value at each wavelength.
- **A real error, found while testing that.** The "plain corner" wasn't a corner. The guides stop at the design region, so it's
  two guides ending at an empty square a wavelength wide. With guides that meet, at 72 cells per wavelength:

| Width | Sharp corner, 0.95, 1, 1.05 | Empty square, 0.95, 1, 1.05 |
|---:|---:|---:|
| 0.06 | 57.4, 57.0, 56.6% | 6.3, 6.9, 7.6% |
| 0.1 | 45.3, 48.9, 51.5% | 5.0, 5.4, 5.8% |
| 0.125 | 23.0, 28.6, 34.2% | 5.4, 5.3, 5.5% |
| 0.143 | 18.7, 17.2, 19.4% | 8.0, 6.7, 6.1% |
| 0.17 | 46.1, 39.4, 29.0% | 9.6, 10.9, 10.7% |
| 0.2 | 39.8, 43.6, 46.0% | 7.2, 8.1, 9.4% |

  So the corner that chapter 13's bend replaces passes 23 to 34%, not 5%, and it now agrees with chapter 12's 18 to 23% at
  0.143. The optimizer really starts from the empty square, and the page now shows both references. The sharp corner's dip
  near 0.143 isn't explained.

Smaller corrections: projection amplifies by up to about $\beta/2$, not $\beta$, in chapter 13 and a comment in `design.js`,
though the step limit is unchanged. "At least 2.7 times" predated the factorization change, and is now 2.9.

**Faster band diagrams: warm start measured, and left out.** The Lanczos starting from the sum of the previous k-point's modes,
about 15 lines, over the path of 16 k-points with 6 bands, and solves per k-point:

| Unknowns | Tolerance 1e-10, cold | 1e-10, warm | 1e-5, cold | 1e-5, warm |
|---:|---:|---:|---:|---:|
| 2,500 | 62.5 | 62.5 | 47.5 | 46.5 |
| 10,000 | 62.0 | 62.5 | 46.5 | 45.5 |
| 40,000 | 62.0 | 61.0 | 51.0 | 46.5 |

With 12 bands at 10,000 unknowns, 88.8, 81.3, 70.8 and 66.2. Warm starts did no better on a path twice as fine. So band diagrams
now use a tolerance of 1e-5, for 6 net lines, and solves that return fields keep 1e-10. A whole diagram in WebAssembly goes from 741 to 576 ms at 2,500 unknowns, 3,444 to 2,692 ms at 10,000 and
16,827 to 14,383 ms at 40,000. The largest change in $k_0$ is 6e-11 away from $\Gamma$, and 9e-9 at $\Gamma$, where $k_0 = 0$ turns
round-off into a square root. Tolerances down to 1e-3 didn't lose or repeat a band, so locking doesn't need the tighter one. Factoring is
now about a third of what's left at 40,000 unknowns, and solves and Lanczos's own work are the rest.

**Dragging shapes.** On views with shapes the visitor drew: drag inside to move, an edge or corner of a rectangle or a circle's
rim to resize, and the source too. `grabShape` and `dragShape` are pure functions in `model.js`, tested in `tests/shapes.test.mjs`.
Mouse drags were checked in a browser. Touch was checked only as far as which touches claim the gesture, not on a phone.

**Periodic gratings.** `solve2d/grating.rs` ties one period along x, with absorbing layers above and below, launches a plane wave
from a sheet, and reads each diffraction order's share of the power from lines of field values. There's a view on the page under
**Driven response**, and chapter 15 teaches it. The estimate was about 190 lines, and it took 390 (2,965 to 3,355), so the estimate was
low by a factor of 2: `grating.rs` 227, the transposed solve 50, and most of the rest generalizing `open.rs` and `validate.rs`.

- **Reference:** `plans/probe-grating/rcwa.py`, an independent Fourier modal method with no mesh or layer, steady to 1e-6 at 320
  orders. Both polarizations agree within 4e-5 at cells of λ/80, and power balances within 2e-5.
- **Two absorbing-layer errors found by it, neither visible to mesh refinement.** A grazing order absorbs at a rate of `ky / k`, so
  the layer's strength is set for the slowest one. An evanescent order just past cutoff isn't absorbed at all and returns from the wall,
  so the stretch gets a real part, set the same way. A first guess, that the steep layer reflected into the transmission reading, was
  wrong: measured directly, the up-going part was about 1e-5.
- **Rayleigh anomalies:** about 1% off exactly on a cutoff, under 0.2% within 0.1% of it, and about 1e-4 at 7% away.
- **A boundary off the grid is harmless.** Edge nodes can't move, so a layer crossing the periodic edge off a grid row is stepped in the
  last column. It put under 4e-6 of the power in other orders.
- **Moharam et al. 1995:** Crossref confirms the journal, volume, issue and first page. Its author order differs from the usual
  citation's, which the references follow, and it gives no last page.

**The tutorial and code now separate the essentials.** Chapters 1 to 13 read straight through, and 14 to 18 are "Going further", each
standing alone: gratings are 15, and 16 to 18 hold what used to interrupt chapters 9, 10, 12 and 13. References are 19. So where this plan
says chapter 13 about a cavity, read 17, and chapter 10 about speed, read 16. Essentials are about 3,000 lines of Rust, with 343 apart from them
in `solve2d/grating.rs` and `sparse/`. Renames: `solve.rs` to `solve1d.rs` with `output.rs` split out, `gradient.rs` to `adjoint.rs`,
`resonance.rs` to `leaky.rs`, `eigs_with` to `eigs_factored`, and in JavaScript `optimizer.js` and `cavity-design.js` out of `design.js`.
Tests split the same way: `analytic.rs` for 1D, `analytic2d.rs` and `gratings.rs`.

### Vector fields: a line estimate

The target is the first real need: full-vector modes of a waveguide's cross-section, at a fixed wavelength, in closed walls.
This is an estimate from reading the code, with no prototype behind it. Gratings ran at twice their estimate, so the right column allows for that.

| Part | What it needs | First guess | Allowing for overrun |
|---|---|---:|---:|
| Nédélec triangles | Basis, curl, local matrices with signs. `local_matrices` assumes a scalar's x and y stiffness, so this needs its own form | 110 | 200 |
| Signs in assembly | An edge's direction decides a dof's sign. Each caller of `cell_dofs` must honour it, or a new path must be added | 30 | 60 |
| Mesh edges | Global edge directions, and ties for periodic walls if wanted | 15 | 30 |
| Mode problem | The tangential and axial blocks, with the Lee transformation so that `B` stays positive semi-definite | 110 | 220 |
| Eigensolver | `lanczos` assumes `B` is positive definite, and `B` here is singular on the axial block unless transformed. May need care | 20 | 60 |
| Request, validation, output | A new analysis, and three field components each on two kinds of dof | 50 | 90 |
| **Solver total** | | **335** | **660** |
| Tests | A hollow metal guide's exact modes, no spurious modes in a uniform guide, and a dielectric guide against a reference | 150 | 300 |
| Page | A new view, and drawing a vector field | 100 | 200 |

So **about 350 to 650 lines of Rust**, or 10 to 20% on the solver, all of it in Going further files apart from the signs in assembly. That's
modes only. Driven vector problems need ports and a tensor absorbing layer as well, for perhaps 200 to 400 more, and 3D is out of reach
in a browser anyway, as chapter 14 says.

Recommendation: **don't start yet.** It's real new capability, and chapter 14 would teach from working code, but it's the only item here that reaches into the
essentials, through `Space` and `assemble`, which every reader goes through. Nothing else planned depends on it. If it's wanted,
start with a throwaway probe of the element alone against a hollow guide, as `plans/probe-sparse/` did for the solver. It's about 150
lines, and it would replace this table's guesses with a measurement.

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
- **Memory:** peak use during factoring is 15–30% of `faer` LU's, which raised the limit from
  160,000 to 400,000. See "A finer mesh, and the crystal at 1 µm".
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

Lanczos itself costs about as much as the solves. A warm start from the previous k-point didn't
help, and a looser tolerance did: see "Faster band diagrams" under "Done after milestone 6".

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

That last row is from before the mesh fix below, which takes it to about 420 ms. Warm-up is
real too: the first factor in a fresh process took 818 ms against 311 ms afterwards, which
probably explains the slow first solve on the page.

### Where a bend's time went

About 1.6 s of that 1.9 s was outside the factor and solve, and nearly all of it was meshing.
The sharp bend has 418 rods and 51,529 unknowns (a/6, 5a layer). These are from one run on
battery, in ms, as medians over three rounds alternating old and new builds:

| Step | Native, before / after | Wasm, before / after |
|---|---:|---:|
| `fitted_grid`: snapping | 158 / 0.8 | 425 / 0.9 |
| `fitted_grid`: choosing diagonals (crossing checks) | 26 / 17 | 476 / 15 |
| `fitted_grid`: Delaunay flips | 16 / 12 | 159 / 12 |
| `fitted_grid`: region lookup | 22 / 10 | 511 / 15 |
| `H1Triangle::new` | 5 / 5 | 5 / 5 |
| `assemble_split` | 10 / 10 | 10 / 10 |
| `Csr::from_triplets` | 38 / 38 | 42 / 42 |
| `Ordering::new` | 25 / 25 | 31 / 31 |
| Factor | 237 / 237 | 268 / 246 |
| Solve, with refinement | 27 / 27 | 31 / 31 |
| `monitor_power` | 1.6 / 1.6 | 4.6 / 4.7 |
| **Whole request** | **567 / 385** | **1,966 / 416** |

The first request in a fresh instance went from 2,666 to 1,020 ms, 733 ms of it the factor warming
up. The a/8 bend, with 112,225 unknowns and 504 rods, went from 4,735 to 1,042 ms, and 1,446 to
984 ms natively. Snapping and `hypot` had the same effect there.

Two things caused it, and neither was the outer loop over rods as such:

- **`f64::hypot`.** WebAssembly has no such instruction, so Rust links a careful software version,
  which avoids overflow that coordinates can't reach. `Outline::level` took about 47 ns a call
  with it. Using `sqrt`, a single instruction, made the loops over every rod about 30× faster.
  That was one line, and it's the difference between this table's columns: native only gained 1.5–2×.
- **Snapping levelled and sorted the whole grid for every rod.** Nodes not yet placed are still
  on the grid, so only those within a step of a shape's bounding box can cross it. Looking only
  there took 15 net lines, and turned 184 ms into 1 ms even after the `sqrt` change.

The same `hypot` was in `SparseLu::factor`, through `Complex::norm` on every entry to find the
matrix's scale. `norm_sqr` costs no lines and saves 22 ms, or 8%, in WebAssembly.

In all that's 16 lines net. The meshes are bit-for-bit identical to the old ones on 11 grids,
including rods, overlaps, and shapes cut by the domain or outside it, and monitors match in every
digit. Factoring is now 59% of a bend, which is about where it should be.

Left alone, since none of them is worth its lines: crossing checks and region lookup still test
every outline, 42 ms in all, and a grid of buckets would take about 15 lines, as an estimate, for at most 10%.
`from_triplets` sorts about 900,000 triplets in 42 ms (10%), and is the next thing to try if bends
need more, since sorting by row and then within rows might save about half, though that's untried. Assembly and
`monitor_power` are too small to matter, so the search in `monitor_power` over nearby cells stays.

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

Chapters 9 to 14 now teach 2D: triangles and fitted meshes, sparse LU, open boundaries, waveguides,
inverse design, and why vector fields are harder. They replaced one chapter that had grown into a list of results.

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
