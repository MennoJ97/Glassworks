//! A grating: one period, tied along x with a Bloch phase, between absorbing layers above and
//! below. The ties cost it `A^T = A`, so its adjoint needs a transposed solve.

use super::adjoint::pixel_slopes;
use super::open::{open_mesh, open_output, stretched_matrix};
use super::{coefficients, free_positions, grid_size, steps_in};
use crate::assemble::DofMap;
use crate::banded::Factors;
use crate::output::{FieldMeta, Output};
use crate::problem::{Analysis2d, Request2d};
use crate::space::{H1Triangle, Space};
use crate::sparse::{Ordering, SparseLu};
use crate::{C64, I, ZERO};
use std::collections::HashMap;
use std::f64::consts::PI;

/// A linear reading of the field, as `(dof, coefficient)` pairs on all dofs: `read` sums
/// `coefficient * u[dof]`.
type Reading = Vec<(usize, C64)>;

/// Lines per reflection reading, spread between `reflect_y` and the sheet.
const LINES: usize = 4;

/// The most that the layers are strengthened for a grazing order. Past this, the layer is too
/// steep for the mesh, and an order this close to its cutoff is reflected anyway.
const MAX_STRENGTH: f64 = 10.0;

/// Walls above and below, and each dof on the right tied to its partner on the left.
fn bloch_dofs(req: &Request2d, space: &H1Triangle, pml: f64, phase: f64) -> Result<DofMap, String> {
    let (width, height) = (req.width, req.height + 2.0 * pml);
    let (_, ny) = grid_size(req, width, height)?;
    let positions = space.dof_positions();
    // A dof's row on the grid of dofs, twice as fine as the cells, which pairs left with right.
    let grid_row = |y: f64| (y / (height / (2 * ny) as f64)).round() as usize;
    let is_at = |coord: f64, end: f64| (coord - end).abs() < 1e-9 * width;
    let on_wall = |p: (f64, f64)| grid_row(p.1) == 0 || grid_row(p.1) == 2 * ny;
    let zero: Vec<usize> = (0..positions.len()).filter(|&dof| on_wall(positions[dof])).collect();
    let left: HashMap<_, _> = (0..positions.len())
        .filter(|&dof| is_at(positions[dof].0, 0.0) && !on_wall(positions[dof]))
        .map(|dof| (grid_row(positions[dof].1), dof))
        .collect();
    // `u(width, y) = e^{i phase} u(0, y)`: one period further along, the Bloch wave's phase.
    let phase_factor = C64::from_polar(1.0, phase);
    let ties: Vec<_> = (0..positions.len())
        .filter(|&dof| is_at(positions[dof].0, width) && !on_wall(positions[dof]))
        .map(|dof| (dof, left[&grid_row(positions[dof].1)], phase_factor))
        .collect();
    Ok(DofMap::new(positions.len(), &zero, &ties))
}

/// Points evenly along a line across the period at height `y` on the mesh: each one's `x`,
/// its cell, and the basis there.
fn line_points(req: &Request2d, space: &H1Triangle, y: f64) -> Result<Vec<(f64, usize, [f64; 6])>, String> {
    let near: Vec<usize> = (0..space.n_cells())
        .filter(|&cell| (space.mesh.centroid(cell).1 - y).abs() < 4.0 * req.cell)
        .collect();
    let steps = steps_in(req.width, req.cell / 4.0);
    (0..steps)
        .map(|i| {
            let x = (i as f64 + 0.5) * req.width / steps as f64;
            let cell = *near
                .iter()
                .find(|&&cell| space.contains(cell, x, y))
                .ok_or("a line lies outside the mesh")?;
            Ok((x, cell, space.basis_at(cell, x, y)))
        })
        .collect()
}

/// The material along a line, as `(p, q)`. It has to be uniform for orders to mean anything.
fn line_material(
    space: &H1Triangle,
    coef: &[(C64, C64)],
    points: &[(f64, usize, [f64; 6])],
) -> Result<(f64, f64), String> {
    let region = space.mesh.region[points[0].1];
    if points.iter().any(|point| space.mesh.region[point.1] != region) {
        return Err("a grating's reflect and transmit lines must lie in uniform material".into());
    }
    Ok((coef[region].0.re, coef[region].1.re))
}

/// Adds `weight` times the mean over a line of `u exp(-i kx x)`, which picks out one order.
fn add_order(space: &H1Triangle, points: &[(f64, usize, [f64; 6])], kx: f64, weight: C64, out: &mut Reading) {
    let mut cell_dofs = vec![0; 6];
    for (x, cell, basis) in points {
        space.cell_dofs(*cell, &mut cell_dofs);
        let factor = C64::from_polar(1.0, -kx * x) * weight / points.len() as f64;
        cell_dofs.iter().zip(basis).for_each(|(&dof, &value)| out.push((dof, factor * value)));
    }
}

fn read(reading: &Reading, u: &[C64]) -> C64 {
    reading.iter().map(|&(dof, coefficient)| coefficient * u[dof]).sum()
}

/// Readings of one order's upward and downward amplitudes, from lines spaced `gap` apart.
/// Each line holds `up exp(i ky y) + down exp(-i ky y)`, fitted by least squares.
fn split_order(
    space: &H1Triangle,
    lines: &[Vec<(f64, usize, [f64; 6])>],
    kx: f64,
    ky: f64,
    gap: f64,
) -> (Reading, Reading) {
    let n = lines.len() as f64;
    let wave = |j: usize, sign: f64| C64::from_polar(1.0, sign * ky * gap * j as f64);
    // The 2 x 2 normal equations have `n` on the diagonal and `cross` off it, so their
    // inverse is written out, and each line's reading gets its own weight.
    let cross: C64 = (0..lines.len()).map(|j| wave(j, -2.0)).sum();
    let det = n * n - cross.norm_sqr();
    let (mut up, mut down) = (vec![], vec![]);
    for (j, points) in lines.iter().enumerate() {
        add_order(space, points, kx, (wave(j, -1.0) * n - cross * wave(j, 1.0)) / det, &mut up);
        add_order(space, points, kx, (wave(j, 1.0) * n - cross.conj() * wave(j, -1.0)) / det, &mut down);
    }
    (up, down)
}

/// The share of the incident power in each reflected and transmitted order (chapter 15). With
/// a density and a goal order, also the slope of that order's share for each pixel.
pub(super) fn solve(req: &Request2d) -> Result<Output, String> {
    let Analysis2d::Grating { wavelength, angle_deg, pml, source_y, reflect_y, transmit_y, goal } =
        req.analysis
    else {
        unreachable!("sent here by solve2d::solve")
    };
    let mesh = open_mesh(req, pml)?;
    let space = H1Triangle::new(&mesh);
    let coef = coefficients(req);
    let k0 = 2.0 * PI / wavelength;
    let gap = (source_y - reflect_y) / LINES as f64;
    let lines = (0..LINES)
        .map(|j| line_points(req, &space, reflect_y + pml + j as f64 * gap))
        .collect::<Result<Vec<_>, _>>()?;
    let below = line_points(req, &space, transmit_y + pml)?;
    let above = line_material(&space, &coef, &lines[0])?;
    for line in &lines[1..] {
        if line_material(&space, &coef, line)? != above {
            return Err("a grating's reflect lines must all lie in one material".into());
        }
    }
    let under = line_material(&space, &coef, &below)?;
    let wavenumber = |(p, q): (f64, f64)| k0 * (q / p).sqrt();
    let span = source_y - reflect_y;
    if !(0.1..=1.0).contains(&(span * wavenumber(above) / (2.0 * PI))) {
        return Err("a grating's source must be 0.1 to 1 wavelengths above its reflect line".into());
    }
    let kx = wavenumber(above) * angle_deg.to_radians().sin();
    // The orders that carry power, in a medium with wavenumber `k`, and each one's `ky`.
    let orders = |k: f64| -> Vec<(i32, f64, f64)> {
        let reach = ((k + kx.abs()) * req.width / (2.0 * PI)).ceil() as i32;
        (-reach..=reach)
            .map(|m| (m, kx + 2.0 * PI * m as f64 / req.width))
            .filter(|&(_, km)| km.abs() < k)
            .map(|(m, km)| (m, km, (k * k - km * km).sqrt()))
            .collect()
    };
    // A layer absorbs a propagating order at a rate of `ky / k`, and an evanescent one only
    // through a real stretch, at a rate of `kappa / k`. Both are set for the slowest order: the
    // real part of `strength` for the first, and a negative imaginary part for the second.
    let background = wavenumber((coef[0].0.re, coef[0].1.re));
    let (mut grazing, mut decay) = (background, background);
    for k in [wavenumber(above), wavenumber(under)] {
        let reach = ((k + kx.abs()) * req.width / (2.0 * PI)).ceil() as i32 + 1;
        for km in (-reach..=reach).map(|m| kx + 2.0 * PI * m as f64 / req.width) {
            let rate = (k * k - km * km).abs().sqrt();
            if km.abs() < k {
                grazing = grazing.min(rate);
            } else {
                decay = decay.min(rate);
            }
        }
    }
    let strength =
        C64::new((background / grazing).min(MAX_STRENGTH), -(background / decay).min(MAX_STRENGTH));
    let dofs = bloch_dofs(req, &space, pml, kx * req.width)?;
    let a = stretched_matrix(req, &space, &dofs, &coef, pml, C64::new(k0, 0.0), strength);
    let ordering = Ordering::reused(&a, || free_positions(&space, &dofs));
    let lu = SparseLu::factor(&a, &ordering);
    // The sheet, with the wave's own phase along x.
    let mut u = vec![ZERO; dofs.n_free];
    let mut cell_dofs = vec![0; 6];
    for (x, cell, basis) in line_points(req, &space, source_y + pml)? {
        space.cell_dofs(cell, &mut cell_dofs);
        for (&dof, &value) in cell_dofs.iter().zip(&basis) {
            if let Some((row, factor)) = dofs.map[dof] {
                u[row] += factor.conj() * I * k0 * value * C64::from_polar(1.0, kx * x);
            }
        }
    }
    lu.solve(&mut u);
    let u = dofs.expand(&u);
    let (_, incident) = split_order(&space, &lines, kx, (wavenumber(above).powi(2) - kx * kx).sqrt(), gap);
    let amplitude = read(&incident, &u);
    let power_in = above.0 * (wavenumber(above).powi(2) - kx * kx).sqrt() * amplitude.norm_sqr();
    let mut channels: Vec<(String, Vec<f64>)> = vec![];
    let mut goal_reading = None;
    for (transmitted, name, material) in [(false, "reflected", above), (true, "transmitted", under)] {
        let (mut labels, mut shares) = (vec![], vec![]);
        for (m, km, ky) in orders(wavenumber(material)) {
            // Below the grating only downward waves remain, so one line reads an order.
            let reading = if transmitted {
                let mut reading = vec![];
                add_order(&space, &below, km, C64::new(1.0, 0.0), &mut reading);
                reading
            } else {
                split_order(&space, &lines, km, ky, gap).0
            };
            let ratio = material.0 * ky / power_in;
            let share = ratio * read(&reading, &u).norm_sqr();
            labels.push(m as f64);
            shares.push(share);
            if goal.is_some_and(|wanted| wanted.transmitted == transmitted && wanted.order == m) {
                goal_reading = Some((reading, ratio));
            }
        }
        channels.push((format!("{name}_orders"), labels));
        channels.push((name.into(), shares));
    }
    let mut out = open_output(req, &space, &dofs, pml, &u, FieldMeta { k0, ..Default::default() });
    out.channels.extend(channels);
    if let (Some(density), Some(wanted)) = (&req.density, goal) {
        let (reading, ratio) = goal_reading.ok_or(format!("order {} carries no power here", wanted.order))?;
        // With `s = r |a|^2 / |i|^2`, `ds = 2 Re(g^T du)`, and `A^T w = g` gives `ds = -2 Re(w^T dA u)`.
        // The incident amplitude `i` is in `g` because the mesh's own version of it moves slightly.
        let value = read(&reading, &u);
        let mut w = vec![ZERO; dofs.n_free];
        let mut load = |terms: &Reading, scale: C64| {
            for &(dof, coefficient) in terms {
                if let Some((row, factor)) = dofs.map[dof] {
                    w[row] += factor * coefficient * scale;
                }
            }
        };
        load(&reading, value.conj() * ratio);
        let norm = amplitude.norm_sqr();
        load(&incident, -amplitude.conj() * ratio * value.norm_sqr() / norm);
        lu.solve_transposed(&mut w);
        // Back onto all dofs: the ties enter `dA` through `Z^H`, so `w` picks up conjugated phases.
        let full: Vec<C64> =
            dofs.map.iter().map(|entry| entry.map_or(ZERO, |(row, factor)| w[row] * factor.conj())).collect();
        let (gradient, _) = pixel_slopes(req, density, &space, pml, C64::new(k0, 0.0), &full, &u);
        out.channels.push(("gradient".into(), gradient));
    }
    Ok(out)
}
