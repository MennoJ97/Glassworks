//! The JSON problem description shared with the website.

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
