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
3. **Band diagram speed.** Eigenproblems use the LU with a reused ordering, within 8% of `faer`'s
   Cholesky. Warm-starting Lanczos from the previous k-point is the next speed-up to try.
4. **Transposed solves.** Bloch-periodic driven problems lose `Aᵀ = A`,
   so their adjoints need `Aᵀ` solves. Not needed before periodic gratings.

## Code changes

| Module | Change |
|---|---|
| Website, later | Dragging shapes |
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
- [ ] **Explain the high Q.** Look at the mode's field, where it radiates, and how Q changes with the
      surrounding air and with the size of the design region. Compare the best starts with those that stay
      below Q 10,000. Finding no cause is an acceptable result, as long as it's stated.
- [ ] **Decide what robustness should mean.** A shifted edge moves the resonance, and at a bandwidth
      of Q 1,000, likely by more than its width. So one wavelength may be unreachable for all versions of a
      design. Read each version at its own wavelength, and report how far they spread as a separate
      result, which is what a real device would have to tune away. Check this against step 1's numbers.
- [ ] **Optimize for robustness, if step 1 shows a need.** Project the filtered densities at three
      thresholds, for eroded, nominal and dilated designs ([Wang, Lazarov & Sigmund 2011], to be
      checked before citing), and raise the weakest. This belongs in `design.js`, with no new Rust expected;
      state the line cost if that changes. A threshold's shift in pixels sets the edge move, so measure that
      against step 1. It takes six solves per iteration instead of two, so about three times as
      much computing, for a run that now takes 11 s on the page. Nine starts aren't enough, as this
      plan learned: use 18 at least, with unseen ones for any final comparison.
- [ ] **Teach it.** A section in chapter 13, with a robust design beside the sharpest one, if there's a
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

### After milestone 6, not yet ordered

None of these is planned in detail.

- **Tutorial clean-up.** Small, and best done first, whatever comes next: check the two citations written from
  memory, test the two inferred explanations in chapters 10 and 13, and read chapters 9 to 14 line by line.
- **Faster band diagrams.** Warm-start Lanczos from the previous k-point (open problem 3). Symmetric factors
  don't help these, because Bloch phases make the matrix Hermitian.
- **Dragging shapes.** The remaining website item, in JavaScript only.
- **Periodic gratings.** A port plus Bloch ties, which needs open problem 4.
- **Vector fields.** Full-vector waveguide cross-sections are the first real need for Nédélec elements, and
  the way to full polarization or 3D. This is the largest of these, and works against keeping the solver lean,
  so it needs a line estimate first.

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
