//! Checks a 1D request before any work starts, so that errors say what to change.

use super::complex;
use crate::problem::{Analysis, Boundary, Layer, Request, Source};

/// Refuses meshes that would exhaust a browser tab's memory.
const MAX_DOFS: f64 = 2.0e6;

/// The first problem with a request, as the message the app shows, or `Ok` if there's none.
pub(super) fn validate(req: &Request) -> Result<(), String> {
    let err = |message: &str| Err(message.to_string());
    let index = |layer: &Layer| (complex(layer.eps) * complex(layer.mu)).norm().sqrt().max(1.0);
    if req.layers.is_empty() || req.layers.iter().any(|layer| !(layer.thickness > 0.0)) {
        return err("every layer needs a positive thickness");
    }
    if !(1..=5).contains(&req.mesh.order) || !(req.mesh.max_cell > 0.0) {
        return err("element order must be 1 to 5, and cell size positive");
    }
    // The mesh `solve` would build, which shrinks cells by each layer's index, counted
    // before it's built.
    let cells: f64 = req.layers.iter().map(|layer| layer.thickness * index(layer) / req.mesh.max_cell).sum();
    if cells * req.mesh.order as f64 > MAX_DOFS {
        return err("this mesh is too fine: lower the resolution or shrink the structure");
    }
    let periodic = req.left == Boundary::Periodic;
    if periodic != (req.right == Boundary::Periodic) {
        return err("periodic boundaries must be set on both sides");
    }
    let open = req.left.is_open() || req.right.is_open();
    let valid_wavelength = |wavelength: f64| wavelength > 0.0 && wavelength.is_finite();
    match &req.analysis {
        Analysis::Driven { wavelength, source, .. } => {
            if !valid_wavelength(*wavelength)
                || periodic
                || (matches!(source, Source::PlaneWave) && req.left != Boundary::Port)
            {
                return err("driven solves need a positive wavelength, no periodic walls, and a left port for a plane wave");
            }
        }
        Analysis::DrivenGradient { wavelength, .. } => {
            if !valid_wavelength(*wavelength) || periodic || req.left != Boundary::Port {
                return err("gradients need a positive wavelength, no periodic walls, and a left port");
            }
        }
        Analysis::Sweep { wavelengths, .. } => {
            if !wavelengths.iter().copied().all(valid_wavelength) || periodic || req.left != Boundary::Port {
                return err("spectra need positive wavelengths and a port on the left");
            }
        }
        // Lanczos (chapter 6) needs a Hermitian `A` and a positive definite `B`, which only
        // closed walls and real, positive materials give.
        eigen => {
            if open {
                return err("eigenmode analyses need closed boundaries: use PEC, PMC or periodic");
            }
            if req.layers.iter().any(|layer| {
                layer.eps[1] != 0.0 || layer.mu[1] != 0.0 || layer.eps[0] <= 0.0 || layer.mu[0] <= 0.0
            }) {
                return err("eigenmode analyses support lossless materials only: use a driven solve for loss or metals");
            }
            let bad = match eigen {
                Analysis::Bands { n_k, .. } => !periodic || *n_k < 2,
                Analysis::Waveguide { wavelength, .. } => periodic || !valid_wavelength(*wavelength),
                _ => false,
            };
            if bad {
                return err("band diagrams need periodic walls; waveguides need PEC or PMC walls and a positive wavelength");
            }
        }
    }
    Ok(())
}
