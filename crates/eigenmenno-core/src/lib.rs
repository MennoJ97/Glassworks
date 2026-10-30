//! A small finite element solver for Maxwell eigenmode and driven problems.
//! See `docs/` for a tutorial that walks through these modules in order:
//! `mesh`, `space`, `assemble`, `banded`, `eigen`, `problem`, `solve`.

pub mod assemble;
pub mod banded;
pub mod eigen;
pub mod mesh;
pub mod problem;
pub mod solve;
pub mod space;

pub use num_complex::Complex64 as C64;
pub const ZERO: C64 = C64::new(0.0, 0.0);
pub const ONE: C64 = C64::new(1.0, 0.0);
pub const I: C64 = C64::new(0.0, 1.0);
