# Glassworks

A small finite element solver for Maxwell eigenmode and driven problems. It runs
entirely in the visitor's browser, so a website hosting it only serves static files
and does no computing.

It handles 1D layered structures, and 2D eigenmodes, scattering, leaky cavities and waveguide bends, with [more of 2D planned](plans/2d-solver.md):

- **Resonances** of closed structures
- **Photonic band diagrams** of periodic structures
- **Guided modes** of slab waveguides
- **2D resonances and band diagrams** of shapes in a box or lattice, for either polarization
- **2D scattering** of a plane wave or line source by shapes in open space, with the local
  density of states at a source, optionally averaged over a bandwidth
- **2D leaky modes**, with their quality factor Q and mode area V
- **2D waveguide bends**, from a port or line source to a power monitor
- **Reflection and transmission spectra**, including lossy layers and oblique incidence
- **Fields at a single wavelength**, from a plane wave or a current sheet
- **Inverse design** of mirrors, broadband coatings, and 2D waveguide bends and cavities in
  either polarization, using adjoint gradients of reflection, transmission and absorption, of
  power through a line, or of the local density of states

**This teaches the method, and it doesn't design devices.** The 2D designs use bulk-like indices (3.5, or a permittivity of 12) rather than a real chip's effective index, and a 2D model can't see light leaking out of the plane. So they aren't silicon-on-insulator designs, whichever polarization they're in.

Visitors build their own structure on the page: drag boundaries and shapes, split layers, pick
materials, repeat blocks, and share the result as a link.

**New to the finite element method?** The [tutorial](docs/README.md) explains it using
this solver's source code. The website includes it too, generated from `docs/` at build
time.

## Run it

You need Node.js and a Rust toolchain with the `wasm32-unknown-unknown` target. The script
installs the matching `wasm-bindgen` into `.tools/` and the tutorial's build tools into
`node_modules/`. Neither ships to visitors.

```bash
./scripts/build-web.sh
```

```bash
python3 -m http.server 8123 --directory web
```

Then open <http://localhost:8123>.

Homebrew's Rust has no WebAssembly target, so on a Mac with Homebrew Rust, install
`rustup` as well. The build script finds its toolchain automatically:

```bash
brew install rustup && $(brew --prefix rustup)/bin/rustup toolchain install stable --profile minimal --target wasm32-unknown-unknown
```

## Test it

```bash
cargo test --release -p glassworks-core
```

The optimizer has its own tests, which need the build above:

```bash
npm test
```

```bash
cargo run --release -p glassworks-core --example bench
```

## Put it on a website

After building, upload the contents of `web/` to any static host. There is no server
code, and no server load per visitor beyond serving the files. The host needs to serve
`.wasm` files as `application/wasm`, which most do by default. No special headers are
needed, because the solver is single-threaded.

Visitors' structures stay on their own computers. They only leave if a visitor shares a
link, since the link itself contains the structure. The page loads nothing from third
parties.

`plans/` holds development notes, measurements and probes. It isn't part of the website,
and nothing on the site links to it.

## How it fits together

```text
crates/glassworks-core   The solver, with no browser dependencies
  mesh.rs                Geometry
  space.rs               Elements (the Space trait is the 2D extension point)
  assemble.rs            Sparse matrices and boundary constraints
  banded.rs              Banded LU and RCM reordering, for 1D
  sparse.rs              Multifrontal sparse LU with nested dissection, for 2D
  sparse/blocked.rs      A faster update of each front, with the same result
  sparse/symmetric.rs    A faster factorization when the matrix is symmetric
  sparse/transposed.rs   Transposed solves, which Bloch adjoints need
  eigen.rs               Shift-invert Lanczos
  problem.rs             The JSON problem description
  output.rs              What a solve returns
  solve1d.rs             The physics: from a request to fields, in 1D
  solve1d/gradient.rs    The adjoint for layer designs
  solve1d/validate.rs    Checks on a 1D request
  solve2d.rs             The same for 2D problems, sending each request on
  solve2d/closed.rs      Walls and periodic cells: modes and band diagrams
  solve2d/open.rs        The absorbing layer, shared by the next two
  solve2d/driven.rs      Line sources, plane waves, ports and monitors
  solve2d/leaky.rs       Leaky modes
  solve2d/adjoint.rs     The adjoint: each pixel's slope
  solve2d/grating.rs     Periodic gratings and their diffraction orders
  solve2d/dipole.rs      An in-plane electric dipole, for H_z cavities
  solve2d/validate.rs    Checks on a request
crates/glassworks-wasm   The thin browser interface
web/                     The page
  worker.js              Runs the solver off the main thread
  js/model.js            Geometry, examples, and building solver requests
  js/plot.js             Canvas drawing (the only code that draws)
  js/design.js           Inverse design of layers and bends
  js/optimizer.js        Filter, projection and the method of moving asymptotes
  js/cavity-design.js    Cavity design, with the frequency held on target
  js/hz-design.js        What H_z cavities add: smoothed projection and a held step
  js/solver.js           Worker client; drops stale solves
  js/app.js              Connects the controls
  learn/                 Tutorial pages, generated from docs/
  learn.js               The code panel beside each tutorial chapter
  THIRD_PARTY.md         Licences of the libraries compiled into pkg/
docs/                    The tutorial, in Markdown
scripts/build-docs.mjs   Turns docs/ into web/learn, with equations as MathML
scripts/code-index.mjs   Finds the lines each name in the tutorial points at
plans/                   Development notes and probes, not published with the site
```

## Why it's built this way

All the linear algebra is handwritten, so you can read the whole solver. In 1D a banded LU
is faster than a general sparse library, as [chapter 5](docs/05-linear-systems.md) measures.
In 2D, where bandwidth grows, `sparse.rs` uses the mesh's geometry to cut the fill. In
WebAssembly it beat `faer`'s general sparse LU on these matrices, and came within 8% of its
Cholesky on band diagrams. `faer` is a dependency of the tests only, as an independent reference.

| If you need | Use |
|---|---|
| A readable solver, or modes and spectra in a browser | This project |
| Fast sparse linear algebra or eigensolvers in Rust | [`faer`](https://docs.rs/faer), which this project tests against |
| A serious FEM library with Nédélec elements | [MFEM](https://mfem.org), [NGSolve](https://ngsolve.org), or [FEniCSx](https://fenicsproject.org) |
| Large 3D electromagnetic eigenmode or driven simulations | [Palace](https://github.com/awslabs/palace) |

## Known limits

- 2D materials are real and positive, so there are no metals or loss in 2D yet.
- Eigenmode analyses support lossless materials only. Use a driven solve for loss or
  metals.
- Materials have one constant index at all wavelengths, so "gold-like metal" is a rough
  stand-in.
- Guided modes near cutoff feel the walls at the domain's edge unless the cladding is
  wide, and wide cladding slows the eigensolver. See
  [chapter 6](docs/06-eigenproblems.md) and [chapter 8](docs/08-verification.md).
- Very large problems are limited by the visitor's memory: 1D requests over 2 million
  unknowns, and 2D requests over 400,000, are refused.

## License

The code is released under the [MIT License](LICENSE). The tutorial's text in `docs/` is released
under [CC BY 4.0](docs/LICENSE): reuse it freely, with credit. Code it quotes stays under MIT.
The libraries compiled into the solver keep their own licences, listed in
[web/THIRD_PARTY.md](web/THIRD_PARTY.md), which ships with the site. Update it when `Cargo.lock` changes.
