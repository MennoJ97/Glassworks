# Tolerance probe

Moves a finished design's edges by less than a mesh cell, for milestone 6 in `plans/2d-solver.md`. Not part of the
solver or the website, and not tidied up.

Run these from a folder holding `pkg/` and `des/`, with `pkg/` a WebAssembly build of the solver with `MAX_DOFS` raised.

- `gen.mjs seed` optimizes a cavity and saves it to `des/`.
- `pix.mjs` draws each material pixel as its own rectangle, which is the method that works. See the plan for those that didn't.
- `m1.mjs seed uniform 3` and `m1.mjs seed random 2 sigma draws` measure a design.
- `m2.mjs` measures the crystal, and `b2.mjs seed` the bend.
