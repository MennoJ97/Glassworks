//! The JSON problem description shared with the website.

use crate::mesh::Outline;
use serde::Deserialize;

#[derive(Deserialize, Clone)]
pub struct Layer {
    pub thickness: f64,
    /// Relative permittivity as `[re, im]`. Positive `im` is loss.
    #[serde(default = "unit")]
    pub eps: [f64; 2],
    #[serde(default = "unit")]
    pub mu: [f64; 2],
}

fn unit() -> [f64; 2] {
    [1.0, 0.0]
}

/// `Te` solves for `E_y`, `Tm` for `H_y`.
#[derive(Deserialize, Clone, Copy, PartialEq)]
#[serde(rename_all = "lowercase")]
pub enum Polarization {
    Te,
    Tm,
}

/// `Absorbing` is a reflectionless end; `Port` also launches the incident wave.
#[derive(Deserialize, Clone, Copy, PartialEq)]
#[serde(rename_all = "lowercase")]
pub enum Boundary {
    Pec,
    Pmc,
    Periodic,
    Absorbing,
    Port,
}

impl Boundary {
    pub fn is_open(self) -> bool {
        matches!(self, Boundary::Absorbing | Boundary::Port)
    }
}

/// `order` is the element order (1 to 5). `max_cell` is the largest cell in
/// vacuum; each layer scales it by its refractive index.
#[derive(Deserialize, Clone)]
pub struct MeshSpec {
    pub order: usize,
    pub max_cell: f64,
}

/// `Line` is a current sheet at position `x`.
#[derive(Deserialize, Clone, Copy)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum Source {
    PlaneWave,
    Line { x: f64 },
}

#[derive(Deserialize, Clone)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum Analysis {
    /// Resonances: find `k0`. Lowest modes unless `target_k0` is given.
    Cavity {
        n_modes: usize,
        #[serde(default)]
        target_k0: Option<f64>,
        #[serde(default)]
        bloch_phase: f64,
    },
    /// Photonic band diagram over `k` in `[0, pi / L]`.
    Bands { n_bands: usize, n_k: usize },
    /// Guided modes: find `beta` at a fixed wavelength.
    Waveguide { wavelength: f64, n_modes: usize },
    /// The field at one wavelength.
    Driven {
        wavelength: f64,
        #[serde(default)]
        angle_deg: f64,
        source: Source,
    },
    /// Reflection and transmission spectrum of a plane wave.
    Sweep {
        wavelengths: Vec<f64>,
        #[serde(default)]
        angle_deg: f64,
    },
    /// A plane wave's field and power at one wavelength, plus the gradients of
    /// that power with respect to the permittivity of each layer not at an open end.
    DrivenGradient {
        wavelength: f64,
        #[serde(default)]
        angle_deg: f64,
    },
}

/// A real, lossless material.
#[derive(Deserialize, Clone, Copy)]
pub struct Material {
    pub eps: f64,
    #[serde(default = "one")]
    pub mu: f64,
}

fn one() -> f64 {
    1.0
}

/// An outline filled with a material. Later shapes sit on top of earlier ones.
#[derive(Deserialize, Clone)]
pub struct Shape {
    #[serde(flatten)]
    pub outline: Outline,
    pub material: Material,
}

#[derive(Deserialize, Clone)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum Analysis2d {
    /// The lowest modes, or those nearest `target_k0`, at a Bloch phase per period if
    /// the walls are periodic.
    Cavity {
        n_modes: usize,
        #[serde(default)]
        bloch: [f64; 2],
        #[serde(default)]
        target_k0: Option<f64>,
    },
    /// Band diagram of a rectangular lattice, along Gamma, X, M and back.
    Bands { n_bands: usize, n_k: usize },
    /// The field at one wavelength in open space. An absorbing layer of thickness `pml`
    /// surrounds the domain, so walls must be `absorbing`. `average_q` solves at
    /// `k0 (1 + i / 2Q)`, which averages a line source's LDOS over a Lorentzian of that Q.
    Driven {
        wavelength: f64,
        source: Source2d,
        pml: f64,
        #[serde(default)]
        average_q: Option<f64>,
    },
    /// The leaky mode nearest `target_k0` in open space, excited by a line source that
    /// must overlap it. Walls must be `absorbing`.
    Resonance { target_k0: f64, x: f64, y: f64, pml: f64 },
}

/// `Line` is a current filament at a point. `PlaneWave` travels at `angle_deg` from the
/// x-axis through the background, with amplitude 1, and the field returned is the total.
/// `Port` launches a layered guide's strongest mode both ways along x, from a current
/// sheet on the line at `x` between `y0` and `y1`, across which the guide must be layered.
#[derive(Deserialize, Clone, Copy)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum Source2d {
    Line { x: f64, y: f64 },
    PlaneWave { angle_deg: f64 },
    Port { x: f64, y0: f64, y1: f64 },
}

/// A straight line that power crossing it is measured along. Power counts as positive
/// when it crosses from the left of the direction `(x0, y0)` to `(x1, y1)` to the right.
#[derive(Deserialize, Clone, Copy)]
pub struct Monitor {
    pub x0: f64,
    pub y0: f64,
    pub x1: f64,
    pub y1: f64,
}

/// A 2D problem on `[0, width] x [0, height]`. `cell` is the largest grid spacing.
#[derive(Deserialize, Clone)]
pub struct Request2d {
    pub width: f64,
    pub height: f64,
    pub cell: f64,
    pub background: Material,
    pub shapes: Vec<Shape>,
    pub polarization: Polarization,
    pub walls: Boundary,
    pub analysis: Analysis2d,
    /// Lines to measure power along, for driven solves.
    #[serde(default)]
    pub monitors: Vec<Monitor>,
}

#[derive(Deserialize, Clone)]
pub struct Request {
    pub layers: Vec<Layer>,
    pub polarization: Polarization,
    pub left: Boundary,
    pub right: Boundary,
    pub mesh: MeshSpec,
    pub analysis: Analysis,
}
