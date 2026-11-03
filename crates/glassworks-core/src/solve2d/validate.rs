//! Checks a 2D request before any work starts, so that errors say what to change.

use super::{side_layer, steps_in};
use crate::problem::{Analysis2d, Boundary, Polarization, Request2d, Source2d};

pub(super) fn validate(req: &Request2d) -> Result<(), String> {
    if !(req.width > 0.0 && req.height > 0.0 && req.cell > 0.0) {
        return Err("the domain and its cell size must be positive".into());
    }
    let positive = |m: crate::problem::Material| m.eps > 0.0 && m.mu > 0.0;
    if !positive(req.background) || !req.shapes.iter().all(|s| positive(s.material)) {
        return Err("eigenmode analyses support lossless materials with positive eps and mu".into());
    }
    if let Some(d) = &req.density {
        let plane_wave =
            matches!(req.analysis, Analysis2d::Driven { source: Source2d::PlaneWave { .. }, .. });
        let line = matches!(
            req.analysis,
            Analysis2d::Driven { source: Source2d::Line { .. } | Source2d::Dipole { .. }, .. }
        );
        let grating = matches!(req.analysis, Analysis2d::Grating { goal: Some(_), .. });
        if !(matches!(req.analysis, Analysis2d::Driven { .. }) || grating)
            || plane_wave
            || (req.monitors.is_empty() && !line && !grating)
        {
            return Err("a density map needs a driven solve from a line source, dipole or port, with a monitor for a port, or a grating with a goal".into());
        }
        if d.rho.len() != d.nx * d.ny
            || !(d.width > 0.0 && d.height > 0.0 && d.eps_low > 0.0 && d.eps_high > 0.0)
        {
            return Err("the density map's size or materials are invalid".into());
        }
        // A monitor reads `p` by shape, which isn't `H_z`'s `p` inside the pixels.
        let inside =
            |x: f64, y: f64| (d.x0..=d.x0 + d.width).contains(&x) && (d.y0..=d.y0 + d.height).contains(&y);
        if req.polarization == Polarization::Tm
            && req.monitors.iter().any(|m| inside(m.x0, m.y0) || inside(m.x1, m.y1))
        {
            return Err("an H_z monitor must stay outside a density map".into());
        }
        // A pixel edge through a cell would make its material depend on where the cell's centre
        // falls, which an optimizer can exploit, and which no drawn structure would reproduce.
        let (Analysis2d::Driven { pml, .. } | Analysis2d::Grating { pml, .. }) = req.analysis else {
            unreachable!("checked above")
        };
        let side = side_layer(req, pml);
        let (hx, hy) = (
            (req.width + 2.0 * side) / steps_in(req.width + 2.0 * side, req.cell) as f64,
            (req.height + 2.0 * pml) / steps_in(req.height + 2.0 * pml, req.cell) as f64,
        );
        let on_grid = |at: f64, h: f64| (at / h - (at / h).round()).abs() < 1e-6;
        let aligned = (0..=d.nx).all(|i| on_grid(side + d.x0 + i as f64 * d.width / d.nx as f64, hx))
            && (0..=d.ny).all(|j| on_grid(pml + d.y0 + j as f64 * d.height / d.ny as f64, hy));
        if !aligned {
            return Err("the density map's pixels must line up with the mesh: choose a cell that divides each pixel exactly".into());
        }
    }
    if let Analysis2d::Grating { wavelength, angle_deg, pml, source_y, reflect_y, transmit_y, .. } =
        req.analysis
    {
        if req.walls != Boundary::Periodic || !(wavelength > 0.0 && pml > 0.0 && angle_deg.abs() < 90.0) {
            return Err(
                "gratings need periodic walls, a positive wavelength and PML, and an angle under 90°".into(),
            );
        }
        if !(0.0 <= transmit_y && transmit_y < reflect_y && reflect_y < source_y && source_y <= req.height) {
            return Err(
                "a grating needs its lines in order inside the domain, from the bottom: transmit, reflect, source. A taller domain or a shorter wavelength makes room".into()
            );
        }
    }
    let driven = matches!(req.analysis, Analysis2d::Driven { .. } | Analysis2d::Resonance { .. });
    if !driven && matches!(req.walls, Boundary::Absorbing | Boundary::Port) {
        return Err("2D eigenmode analyses need PEC, PMC or periodic walls".into());
    }
    if let Analysis2d::Resonance { target_k0, x, y, pml, dipole_deg } = req.analysis {
        if dipole_deg.is_some() && req.polarization != Polarization::Tm {
            return Err("a dipole in the plane needs H_z, since its field has no E_z".into());
        }
        if req.walls != Boundary::Absorbing || !(target_k0 > 0.0 && pml > 0.0) {
            return Err("2D leaky modes need absorbing walls, a positive target and a positive PML".into());
        }
        if !(0.0..=req.width).contains(&x) || !(0.0..=req.height).contains(&y) {
            return Err("the source lies outside the domain".into());
        }
    }
    if let Analysis2d::Driven { wavelength, pml, source, average_q } = req.analysis {
        if req.walls != Boundary::Absorbing || !(wavelength > 0.0 && pml > 0.0) {
            return Err(
                "2D driven solves need absorbing walls, a positive wavelength and a positive PML".into()
            );
        }
        if average_q.is_some_and(|q| !(q > 0.0))
            || (average_q.is_some() && matches!(source, Source2d::PlaneWave { .. }))
        {
            return Err("averaging over frequency needs a line source and a positive Q".into());
        }
        match source {
            Source2d::Line { x, y } | Source2d::Dipole { x, y, .. }
                if !(0.0..=req.width).contains(&x) || !(0.0..=req.height).contains(&y) =>
            {
                return Err("the source lies outside the domain".into());
            }
            Source2d::Dipole { .. } if req.polarization != Polarization::Tm => {
                return Err("a dipole in the plane needs H_z, since its field has no E_z".into());
            }
            // The incident wave is defined in the background, so scatterers can't reach the absorber.
            Source2d::PlaneWave { .. }
                if !req.shapes.iter().all(|s| s.outline.inside(req.width, req.height)) =>
            {
                return Err("shapes must stay inside the domain for a plane wave".into());
            }
            Source2d::Port { x, y0, y1 }
                if !(0.0..=req.width).contains(&x) || !(0.0 <= y0 && y0 < y1 && y1 <= req.height) =>
            {
                return Err("the port lies outside the domain".into());
            }
            _ => {}
        }
    }
    if matches!(req.analysis, Analysis2d::Bands { .. }) && req.walls != Boundary::Periodic {
        return Err("band diagrams need periodic walls".into());
    }
    Ok(())
}
