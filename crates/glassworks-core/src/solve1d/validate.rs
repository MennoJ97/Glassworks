//! Checks a 1D request before any work starts, so that errors say what to change.

use super::complex;
use crate::problem::{Analysis, Boundary, Layer, Request, Source};

/// Refuses meshes that would exhaust a browser tab's memory.
const MAX_DOFS: f64 = 2.0e6;

pub(super) fn validate(req: &Request) -> Result<(), String> {
    let err = |m: &str| Err(m.to_string());
    let index = |layer: &Layer| (complex(layer.eps) * complex(layer.mu)).norm().sqrt().max(1.0);
    if req.layers.is_empty() || req.layers.iter().any(|l| !(l.thickness > 0.0)) {
        return err("every layer needs a positive thickness");
    }
    if !(1..=5).contains(&req.mesh.order) || !(req.mesh.max_cell > 0.0) {
        return err("element order must be 1 to 5, and cell size positive");
    }
    let cells: f64 = req.layers.iter().map(|l| l.thickness * index(l) / req.mesh.max_cell).sum();
    if cells * req.mesh.order as f64 > MAX_DOFS {
        return err("this mesh is too fine: lower the resolution or shrink the structure");
    }
    let periodic = req.left == Boundary::Periodic;
    if periodic != (req.right == Boundary::Periodic) {
        return err("periodic boundaries must be set on both sides");
    }
    let open = req.left.is_open() || req.right.is_open();
    let ok_wl = |w: f64| w > 0.0 && w.is_finite();
    match &req.analysis {
        Analysis::Driven { wavelength, source, .. } => {
            if !ok_wl(*wavelength)
                || periodic
                || (matches!(source, Source::PlaneWave) && req.left != Boundary::Port)
            {
                return err("driven solves need a positive wavelength, no periodic walls, and a left port for a plane wave");
            }
        }
        Analysis::DrivenGradient { wavelength, .. } => {
            if !ok_wl(*wavelength) || periodic || req.left != Boundary::Port {
                return err("gradients need a positive wavelength, no periodic walls, and a left port");
            }
        }
        Analysis::Sweep { wavelengths, .. } => {
            if !wavelengths.iter().copied().all(ok_wl) || periodic || req.left != Boundary::Port {
                return err("spectra need positive wavelengths and a port on the left");
            }
        }
        eigen => {
            if open {
                return err("eigenmode analyses need closed boundaries: use PEC, PMC or periodic");
            }
            if req
                .layers
                .iter()
                .any(|l| l.eps[1] != 0.0 || l.mu[1] != 0.0 || l.eps[0] <= 0.0 || l.mu[0] <= 0.0)
            {
                return err("eigenmode analyses support lossless materials only: use a driven solve for loss or metals");
            }
            let bad = match eigen {
                Analysis::Bands { n_k, .. } => !periodic || *n_k < 2,
                Analysis::Waveguide { wavelength, .. } => periodic || !ok_wl(*wavelength),
                _ => false,
            };
            if bad {
                return err("band diagrams need periodic walls; waveguides need PEC or PMC walls and a positive wavelength");
            }
        }
    }
    Ok(())
}
