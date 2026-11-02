//! What every solve hands back to the website, in 1D and 2D alike.

use crate::{C64, ZERO};
use serde::Serialize;

/// Numbers travel as `channels`, so large fields stay cheap to hand to
/// JavaScript. Fields are named `field{i}/x` and `field{i}/{component}/re|im`.
pub struct Output {
    pub meta: Meta,
    pub channels: Vec<(String, Vec<f64>)>,
}

#[derive(Serialize, Default)]
pub struct Meta {
    pub n_dofs: usize,
    pub length: f64,
    pub fields: Vec<FieldMeta>,
}

#[derive(Serialize, Default)]
pub struct FieldMeta {
    pub components: Vec<&'static str>,
    pub k0: f64,
    /// Effective index of a waveguide mode.
    pub n_eff: Option<f64>,
    /// Reflected, transmitted and absorbed power of a plane wave.
    pub power: Option<[f64; 3]>,
    /// Quality factor of a leaky 2D mode.
    pub q: Option<f64>,
    /// Effective area of a 2D mode: energy over its peak density.
    pub volume: Option<f64>,
    /// `2 Re(-u)` at a 2D line source, which is pi times the LDOS per unit k.
    pub ldos: Option<f64>,
    /// `-u_0`, whose reciprocal is nearly linear in `k` close to a resonance, as `[re, im]`.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub response: Option<[f64; 2]>,
}

/// Scales a mode to a unit real peak, so plots are comparable.
pub(crate) fn normalized(mut u: Vec<C64>) -> Vec<C64> {
    let peak = u.iter().fold(ZERO, |m, x| if x.norm() > m.norm() { *x } else { m });
    if peak != ZERO {
        u.iter_mut().for_each(|x| *x /= peak);
    }
    u
}
