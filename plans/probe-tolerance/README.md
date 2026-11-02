# Tolerance probe

Moves a finished design's edges by less than a mesh cell, for milestone 6 in `plans/2d-solver.md`. Not part of the
solver or the website, and not tidied up.

Run these from a folder holding `pkg/` and `des/`, with `pkg/` a WebAssembly build of the solver with `MAX_DOFS` raised.

- `gen.mjs seed` optimizes a cavity and saves it to `des/`.
- `pix.mjs` draws each material pixel as its own rectangle, which is the method that works. See the plan for those that didn't.
- `m1.mjs seed uniform 3` and `m1.mjs seed random 2 sigma draws` measure a design.
- `m2.mjs` measures the crystal, and `b2.mjs seed` the bend.
- `sym.mjs` compares free and symmetric disorder, `ch2.mjs` splits a mode's radiation into channels, and `where.mjs` shows where its field sits.
- `sz.mjs` scales the design region, `tr.mjs` forces a disc around the source, and `m3.mjs` measures named designs.
- `robust-goal.patch` is the prototype robust goal, against `web/js/design.js`. `rob.mjs` runs and judges it, and `sumr.py` summarizes.

`hz/` measures `H_z` designs against the air-hole crystal, for "Tolerance in `H_z`" in the plan. Run its scripts from a folder holding them with `lib.mjs`, `pix.mjs`, `des/`, `pkg/` and `src/`, a checkout of the repository. `genhz.mjs seed` optimizes, `tol.mjs seed draws` and `corr.mjs seed L draws` measure, and `summary.py` and `summary_corr.py` read the logs.
