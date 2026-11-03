//! What the test files share. Not every test file uses all of it.
#![allow(dead_code)]

use glassworks_core::output::Output;
use glassworks_core::problem::{Analysis, Boundary, Layer, MeshSpec, Polarization, Request};

pub fn channel<'a>(o: &'a Output, name: &str) -> &'a [f64] {
    &o.channels.iter().find(|(n, _)| n == name).unwrap_or_else(|| panic!("no channel {name}")).1
}

pub fn close(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a} differs from {b} by more than {tol}");
}

pub fn layer(thickness: f64, n: f64) -> Layer {
    Layer { thickness, eps: [n * n, 0.0], mu: [1.0, 0.0] }
}

pub fn request(
    layers: Vec<Layer>,
    pol: Polarization,
    b: (Boundary, Boundary),
    order: usize,
    h: f64,
    analysis: Analysis,
) -> Request {
    Request {
        layers,
        polarization: pol,
        left: b.0,
        right: b.1,
        mesh: MeshSpec { order, max_cell: h },
        analysis,
    }
}

pub fn cavity(n_modes: usize) -> Analysis {
    Analysis::Cavity { n_modes, target_k0: None, bloch_phase: 0.0 }
}
