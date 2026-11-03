# 8. Is the answer right?

A wrong field looks just as convincing as a right one, so a solver is only useful
if its answers are checked. Here's how this one is checked, including what those
checks turned up while it was being built.

## Test against exact results, not stored numbers

Each test in [`tests/analytic.rs`](../crates/eigenmenno-core/tests/analytic.rs)
compares against an independently known answer:

| Test | What it's compared to |
|---|---|
| [Cavity modes](../crates/eigenmenno-core/tests/analytic.rs "pec_cavity_matches_mode_formula") | $k_0 = m\pi/L$ |
| [Layered cavity](../crates/eigenmenno-core/tests/analytic.rs "layered_cavity_satisfies_transcendental_equation") | The transcendental equation from matching fields at the interface |
| [Bragg stack bands](../crates/eigenmenno-core/tests/analytic.rs "bloch_bands_satisfy_bragg_stack_dispersion") | The analytic Bloch dispersion relation ([Yeh et al. 1977](20-references.md#yeh-1977)) |
| [Degenerate Bloch modes](../crates/eigenmenno-core/tests/analytic.rs "bloch_modes_of_uniform_medium_include_degenerate_pairs") | The double modes of a uniform medium |
| [Slab waveguide, TE and TM](../crates/eigenmenno-core/tests/analytic.rs "slab_waveguide_te_and_tm_satisfy_dispersion") | The slab dispersion equations ([Yariv & Yeh 2007](20-references.md#yariv-yeh-2007)) |
| [Single interface at 40°](../crates/eigenmenno-core/tests/analytic.rs "single_interface_matches_fresnel_at_oblique_incidence") | The Fresnel equations ([Born & Wolf 1999](20-references.md#born-wolf-1999)) |
| [Multilayer spectrum](../crates/eigenmenno-core/tests/analytic.rs "multilayer_spectrum_matches_transfer_matrix"), [lossy layers](../crates/eigenmenno-core/tests/analytic.rs "lossy_layer_absorbs_and_matches_transfer_matrix") | A [separately written transfer-matrix calculation](../crates/eigenmenno-core/tests/analytic.rs "transfer_matrix") ([Yeh 2005](20-references.md#yeh-2005)) |
| [Current sheet](../crates/eigenmenno-core/tests/analytic.rs "line_source_radiates_the_free_space_green_function") | The free-space solution |
| [$E$ and $H$ fields](../crates/eigenmenno-core/tests/analytic.rs "derived_fields_satisfy_faraday_in_a_cavity") | $E_y = \sin\pi x$ and $\eta_0 H_z = -i\cos\pi x$ in a cavity |

[`tests/cross_check.rs`](../crates/eigenmenno-core/tests/cross_check.rs) also compares
the handwritten LUs and Lanczos against `faer`, used only in tests as an independent
reference.

## Test the convergence rate

Testing one mesh can pass by luck. A much stronger check is that the error shrinks at
the theoretically predicted **rate** ([Oberkampf & Roy 2010](20-references.md#oberkampf-roy-2010)).
According to [Strang & Fix 2008](20-references.md#strang-fix-2008), eigenvalue error should scale as

$$|\lambda_h - \lambda| \propto h^{2p}.$$
<!-- code: eigenvalue_error_converges_at_rate_2p "let (e1, e2)" .. "assert!" -->

So halving $h$ should divide the error by $2^{2p}$. This table is real output, for the
lowest mode of a unit cavity, showing the error in $k_0^2$:

| Order | 2 cells | 4 cells | 8 cells | 16 cells | Measured rate | Predicted |
|---|---|---|---|---|---|---|
| 1 | $2.13$ | $5.17\times10^{-1}$ | $1.27\times10^{-1}$ | $3.17\times10^{-2}$ | 2.0 | 2 |
| 2 | $7.42\times10^{-2}$ | $5.05\times10^{-3}$ | $3.23\times10^{-4}$ | $2.03\times10^{-5}$ | 4.0 | 4 |
| 3 | $1.35\times10^{-3}$ | $2.25\times10^{-5}$ | $3.57\times10^{-7}$ | $5.60\times10^{-9}$ | 6.0 | 6 |
| 4 | $1.35\times10^{-5}$ | $5.53\times10^{-8}$ | $2.19\times10^{-10}$ | $1.16\times10^{-12}$ | 8.0 | 8 |
| 5 | $8.50\times10^{-8}$ | $8.66\times10^{-11}$ | $1.18\times10^{-12}$ | $6.52\times10^{-13}$ | 9.9 | 10 |

The rates for orders 4 and 5 are taken before the error hits **round-off**, around
$10^{-12}$. Past that point, a finer mesh can't help, because the limit is now
floating-point arithmetic rather than the discretization.
<!-- code: eigenvalue_error_converges_at_rate_2p "Coarse enough" .. "let (e1, e2)" -->

## A real debugging story

When the tests first ran, 4 of 14 failed, with all four errors being small. The
tempting move is to loosen the tolerances. Instead, each one was studied by varying
one parameter and watching the error respond.

**Cavity test.** The error was $2\times10^{-8}$, against a tolerance of $10^{-8}$.
Halving $h$ divided it by about 250, close to the predicted $2^8 = 256$. So it was
ordinary discretization error on a too-coarse mesh. *Fix: a finer mesh in the test,
with the tolerance unchanged.*
<!-- code: pec_cavity_matches_mode_formula -->

**TM cavity test.** The test expected $k_0 = 0$ and got $3\times10^{-7}$, which looks
inaccurate. But the solver computes $k_0^2$, which was $10^{-13}$, right at round-off.
The square root amplified it. *Fix: compare $k_0^2$, the quantity actually solved for.*
<!-- code: tm_polarization_swaps_pec_and_pmc "Compare k0^2" .. "close(" -->

**Slab waveguide test.** Halving $h$ changed nothing, so the mesh wasn't the cause.
Doubling the cladding width shrank the error about 100 times. That pointed at the
walls: the near-cutoff mode decays so slowly that its tail reaches the conductor at
the domain's edge, with an error scaling as $e^{-2\gamma d}$ for decay rate $\gamma$
and cladding width $d$ ([Yariv & Yeh 2007](20-references.md#yariv-yeh-2007)). This is a real modelling error, and **no mesh refinement will
ever fix it**. *Fix: wider cladding, and a bound in the test based on that exponential.*
<!-- code: slab_waveguide_te_and_tm_satisfy_dispersion "Wide cladding" .. "vec![layer" -->
<!-- code: slab_waveguide_te_and_tm_satisfy_dispersion "Modes near cutoff" .. "close(res" -->

**Current sheet test.** The result had the opposite sign to the test's formula. A
derivation by hand showed the *test* was wrong, as chapter 7 describes. The solver was
right. Tests can have bugs too.
<!-- code: line_source_radiates_the_free_space_green_function "-u'' - k^2 u" .. "let exact" -->

Investigating those failures also exposed **a real bug in the solver**. Trying wider
cladding made the eigensolver report "converged only 2 of 3 modes". The Krylov basis
had a fixed small limit, too small for the cluster of box modes around the third mode.
Loosening tolerances would never have revealed this. *Fix: let the basis grow until
the wanted modes converge.*
<!-- code: MAX_BASIS -->
<!-- code: lanczos "let max_dim" -->

## What these tests don't cover

- **The website's behaviour**, like dragging and sharing, was checked by hand in a
  browser rather than by automated tests.
- **Material data**: the materials in the app use one constant index at all
  wavelengths, so the "gold-like metal" is a rough stand-in rather than real gold.
- **Very large or extreme inputs**: there are size limits, but the solver hasn't been
  stress-tested for every combination, such as very high contrast at very fine meshes.

Next: [fields on triangles](09-triangles.md).
