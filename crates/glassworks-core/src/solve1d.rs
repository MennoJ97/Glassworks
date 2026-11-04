//! The physics: turns a `Request` into matrices, solves them, and packages fields.
//!
//! With time dependence `e^{-iwt}` and fields varying as `u(x) e^{i kt z}`, both
//! polarizations reduce to one scalar equation:
//!
//! ```text
//! -(p u')' + kt^2 p u - k0^2 q u = source
//! ```
//!
//! TE has `u = E_y`, `p = 1/mu`, `q = eps`. TM has `u = eta0 H_y`, `p = 1/eps`,
//! `q = mu`. Multiplying by a test function and integrating by parts gives the
//! matrices `S` (from `p u' v'`) and `M` (from `u v`) that `assemble` builds.

mod gradient;
mod validate;

use crate::assemble::{assemble, Csr, DofMap};
use crate::banded::{BandedLu, Factors};
use crate::eigen::eigs_near;
use crate::mesh::Mesh;
use crate::output::{normalized, FieldMeta, Meta, Output};
use crate::problem::{Analysis, Boundary, Layer, Polarization, Request, Source};
use crate::space::{H1Line, Space};
use crate::{C64, I, ONE, ZERO};
use std::f64::consts::PI;

const EIG_TOL: f64 = 1e-11;
/// Field output cap, so huge meshes cannot flood the plot.
const MAX_FIELD_POINTS: usize = 20_000;

/// A complex number from the `[re, im]` pair a request carries.
fn complex(pair: [f64; 2]) -> C64 {
    C64::new(pair[0], pair[1])
}

/// `(p, q)` of a layer.
fn layer_coef(polarization: Polarization, layer: &Layer) -> (C64, C64) {
    match polarization {
        Polarization::Te => (complex(layer.mu).inv(), complex(layer.eps)),
        Polarization::Tm => (complex(layer.eps).inv(), complex(layer.mu)),
    }
}

/// Everything shared by the analyses for one request.
struct Setup<'a> {
    req: &'a Request,
    space: H1Line<'a>,
    dofs: DofMap,
}

impl<'a> Setup<'a> {
    fn new(req: &'a Request, mesh: &'a Mesh, bloch_phase: f64) -> Self {
        let space = H1Line::new(mesh, req.mesh.order);
        let last = space.n_dofs() - 1;
        // A wall that fixes the solved field is a Dirichlet condition. The other
        // wall type holds naturally: its boundary term in the weak form vanishes.
        let wall = if req.polarization == Polarization::Te { Boundary::Pec } else { Boundary::Pmc };
        let mut zero = vec![];
        if req.left == wall {
            zero.push(0);
        }
        if req.right == wall {
            zero.push(last);
        }
        let tie = (req.left == Boundary::Periodic).then(|| (last, 0, C64::from_polar(1.0, bloch_phase)));
        let dofs = DofMap::new(space.n_dofs(), &zero, tie.as_slice());
        Self { req, space, dofs }
    }

    fn coef(&self, cell: usize) -> (C64, C64) {
        layer_coef(self.req.polarization, &self.req.layers[self.space.mesh.region[cell]])
    }

    /// The matrix of `a S + b M`, with `weights(p, q) = (a, b)` per cell.
    fn matrix(&self, weights: impl Fn(C64, C64) -> (C64, C64)) -> Csr {
        let triplets = assemble(&self.space, &self.dofs, |cell| {
            let (p, q) = self.coef(cell);
            weights(p, q)
        });
        Csr::from_triplets(self.dofs.n_free, triplets)
    }

    /// Samples the field and the components that follow from Maxwell's curl
    /// equations, appending them as channels and metadata.
    fn push_field(
        &self,
        out: &mut Output,
        u: &[C64],
        k0: f64,
        kt: f64,
        n_eff: Option<f64>,
        power: Option<[f64; 3]>,
    ) {
        let mesh = self.space.mesh;
        let te = self.req.polarization == Polarization::Te;
        let per_cell = (2 * self.req.mesh.order).max(2);
        let stride = (mesh.n_cells() * per_cell / MAX_FIELD_POINTS).max(1);
        let mut x = vec![];
        // Real and imaginary parts of the three components, in turn.
        let mut samples = vec![vec![]; 6];
        let scale = if k0 > 0.0 { 1.0 / k0 } else { 0.0 };
        for cell in (0..mesh.n_cells()).step_by(stride) {
            let layer = &self.req.layers[mesh.region[cell]];
            let (x0, h) = (mesh.coords[mesh.cells[2 * cell]], self.space.cell_length(cell));
            // Each cell repeats its end points, so jumps at interfaces plot as jumps.
            for sample in 0..per_cell {
                let xi = -1.0 + 2.0 * sample as f64 / (per_cell - 1) as f64;
                let (value, du) = self.space.eval(cell, xi, u);
                let (second, third) = if te {
                    // eta0 Hx and eta0 Hz.
                    (-value * kt * scale / complex(layer.mu), -I * du * scale / complex(layer.mu))
                } else {
                    // Ex and Ez.
                    (value * kt * scale / complex(layer.eps), I * du * scale / complex(layer.eps))
                };
                x.push(x0 + h * (xi + 1.0) / 2.0);
                for (k, component) in [value, second, third].into_iter().enumerate() {
                    samples[2 * k].push(component.re);
                    samples[2 * k + 1].push(component.im);
                }
            }
        }
        let components = if te { ["Ey", "Hx", "Hz"] } else { ["Hy", "Ex", "Ez"] };
        let i = out.meta.fields.len();
        out.channels.push((format!("field{i}/x"), x));
        for (k, channel) in samples.into_iter().enumerate() {
            out.channels.push((format!("field{i}/{}/{}", components[k / 2], ["re", "im"][k % 2]), channel));
        }
        out.meta.fields.push(FieldMeta {
            components: components.to_vec(),
            k0,
            n_eff,
            power,
            ..Default::default()
        });
    }
}

pub fn solve(req: &Request) -> Result<Output, String> {
    validate::validate(req)?;
    let thickness: Vec<f64> = req.layers.iter().map(|layer| layer.thickness).collect();
    // Scale cells by local index so each layer has the same cells per wavelength.
    let index = |layer: &Layer| (complex(layer.eps) * complex(layer.mu)).norm().sqrt().max(1.0);
    let max_cell: Vec<f64> = req.layers.iter().map(|layer| req.mesh.max_cell / index(layer)).collect();
    let mesh = Mesh::layered_line(&thickness, &max_cell);
    let length: f64 = thickness.iter().sum();
    let bloch_phase = if let Analysis::Cavity { bloch_phase, .. } = req.analysis { bloch_phase } else { 0.0 };
    let setup = Setup::new(req, &mesh, bloch_phase);
    let meta = Meta { length, n_dofs: setup.dofs.n_free, ..Default::default() };
    let mut out = Output { meta, channels: vec![] };
    let k0_of = |wavelength: f64| 2.0 * PI / wavelength;
    match &req.analysis {
        Analysis::Cavity { n_modes, target_k0, .. } => {
            for (lambda, x) in cavity_modes(&setup, *n_modes, *target_k0, length)? {
                setup.push_field(
                    &mut out,
                    &normalized(setup.dofs.expand(&x)),
                    lambda.max(0.0).sqrt(),
                    0.0,
                    None,
                    None,
                );
            }
        }
        Analysis::Bands { n_bands, n_k } => bands(&setup, *n_bands, *n_k, &mut out)?,
        Analysis::Waveguide { wavelength, n_modes } => {
            waveguide(&setup, k0_of(*wavelength), *n_modes, &mut out)?
        }
        Analysis::Driven { wavelength, angle_deg, source } => {
            let k0 = k0_of(*wavelength);
            let solved = driven(&setup, k0, *angle_deg, *source)?;
            setup.push_field(&mut out, &solved.u, k0, solved.kt, None, solved.power);
        }
        Analysis::DrivenGradient { wavelength, angle_deg } => {
            let k0 = k0_of(*wavelength);
            let solved = driven(&setup, k0, *angle_deg, Source::PlaneWave)?;
            setup.push_field(&mut out, &solved.u, k0, solved.kt, None, solved.power);
            gradient::power_gradients(&setup, k0, &solved, &mut out);
        }
        Analysis::Sweep { wavelengths, angle_deg } => {
            let mut series = [vec![], vec![], vec![]];
            for &wavelength in wavelengths {
                let power = driven(&setup, k0_of(wavelength), *angle_deg, Source::PlaneWave)?.power.unwrap();
                (0..3).for_each(|i| series[i].push(power[i]));
            }
            out.channels.push(("wavelength".into(), wavelengths.clone()));
            for (name, data) in ["R", "T", "A"].into_iter().zip(series) {
                out.channels.push((name.into(), data));
            }
        }
    }
    Ok(out)
}

/// `S u = k0^2 Q u`, with `Q` the mass matrix weighted by `q`.
fn cavity_modes(
    setup: &Setup,
    n_modes: usize,
    target_k0: Option<f64>,
    length: f64,
) -> Result<Vec<(f64, Vec<C64>)>, String> {
    let a = setup.matrix(|p, _| (p, ZERO));
    let b = setup.matrix(|_, q| (ZERO, q));
    // A zero shift would make S singular whenever a static mode exists, as in a TM
    // cavity with PEC walls. A negative shift keeps S - sigma Q positive definite,
    // and the lowest modes are still the ones closest to it.
    let sigma = target_k0.map_or(-(PI / length).powi(2), |k0| k0 * k0);
    eigs_near(&a, &b, sigma, n_modes, EIG_TOL)
}

fn bands(setup: &Setup, n_bands: usize, n_k: usize, out: &mut Output) -> Result<(), String> {
    let length = out.meta.length;
    let ks: Vec<f64> = (0..n_k).map(|i| PI / length * i as f64 / (n_k - 1) as f64).collect();
    let mut table = vec![vec![]; n_bands];
    for &k in &ks {
        let at_k = Setup::new(setup.req, setup.space.mesh, k * length);
        for (band, (lambda, _)) in table.iter_mut().zip(cavity_modes(&at_k, n_bands, None, length)?) {
            band.push(lambda.max(0.0).sqrt());
        }
    }
    out.channels.push(("k".into(), ks));
    out.channels.extend(table.into_iter().enumerate().map(|(i, band)| (format!("band{i}"), band)));
    Ok(())
}

/// `(k0^2 Q - S) u = beta^2 P u`, with `P` the mass matrix weighted by `p`.
fn waveguide(setup: &Setup, k0: f64, n_modes: usize, out: &mut Output) -> Result<(), String> {
    let k02 = C64::new(k0 * k0, 0.0);
    let a = setup.matrix(|p, q| (-p, k02 * q));
    let b = setup.matrix(|p, _| (ZERO, p));
    let n2: Vec<f64> = setup.req.layers.iter().map(|layer| layer.eps[0] * layer.mu[0]).collect();
    let n2_max = n2.iter().copied().fold(0.0, f64::max);
    let n2_clad = n2[0].max(n2[n2.len() - 1]);
    // Highest effective index first. Modes at or below the cladding index are box
    // modes of the truncated domain, not guided modes, so they are dropped.
    for (lambda, x) in eigs_near(&a, &b, k0 * k0 * n2_max * 1.001, n_modes, EIG_TOL)?.into_iter().rev() {
        if lambda > k0 * k0 * n2_clad {
            let beta = lambda.sqrt();
            setup.push_field(out, &normalized(setup.dofs.expand(&x)), k0, beta, Some(beta / k0), None);
        }
    }
    Ok(())
}

/// A solved driven problem, with `power` set for a plane wave. The factors and
/// `t_scale` (T per `|u|^2` at the right end) stay around for adjoint solves.
struct Driven {
    u: Vec<C64>,
    kt: f64,
    power: Option<[f64; 3]>,
    lu: BandedLu,
    t_scale: f64,
}

fn driven(setup: &Setup, k0: f64, angle_deg: f64, source: Source) -> Result<Driven, String> {
    let req = setup.req;
    let n2 = |layer: &Layer| complex(layer.eps) * complex(layer.mu);
    let plane_wave = matches!(source, Source::PlaneWave);
    let (first, last_layer) = (&req.layers[0], &req.layers[req.layers.len() - 1]);
    let kt = if plane_wave { k0 * n2(first).sqrt().re * angle_deg.to_radians().sin() } else { 0.0 };
    let (k02, kt2) = (C64::new(k0 * k0, 0.0), C64::new(kt * kt, 0.0));
    let mut triplets = assemble(&setup.space, &setup.dofs, |cell| {
        let (p, q) = setup.coef(cell);
        (p, kt2 * p - k02 * q)
    });
    // `p kx` at an end, with `kx` on the outgoing or decaying branch. An open end
    // has `u' = i kx u`, which turns the weak form's boundary term into this entry.
    let admittance = |layer: &Layer| {
        let kx = (k02 * n2(layer) - kt2).sqrt();
        let outgoing = if kx.im < 0.0 || (kx.im == 0.0 && kx.re < 0.0) { -kx } else { kx };
        layer_coef(req.polarization, layer).0 * outgoing
    };
    let last = setup.space.n_dofs() - 1;
    let (y_left, y_right) = (admittance(first), admittance(last_layer));
    for (is_open, dof, y) in [(req.left.is_open(), 0, y_left), (req.right.is_open(), last, y_right)] {
        if is_open {
            let (row, _) = setup.dofs.map[dof].unwrap();
            triplets.push((row, row, -I * y));
        }
    }
    let mut rhs = vec![ZERO; setup.dofs.n_free];
    match source {
        Source::PlaneWave => rhs[setup.dofs.map[0].unwrap().0] = -I * y_left * 2.0,
        Source::Line { x } => add_point_load(setup, x, I * k0, &mut rhs)?,
    }
    let lu = BandedLu::factor(&Csr::from_triplets(setup.dofs.n_free, triplets))?;
    lu.solve(&mut rhs);
    let u = setup.dofs.expand(&rhs);
    let t_scale = if req.right.is_open() { y_right.re / y_left.re } else { 0.0 };
    let power = plane_wave.then(|| {
        let reflected = (u[0] - ONE).norm_sqr();
        let transmitted = u[last].norm_sqr() * t_scale;
        [reflected, transmitted, 1.0 - reflected - transmitted]
    });
    Ok(Driven { u, kt, power, lu, t_scale })
}

/// Adds `amp * v(x)` to the right-hand side: the weak form of a current sheet.
fn add_point_load(setup: &Setup, x: f64, amp: C64, rhs: &mut [C64]) -> Result<(), String> {
    let mesh = setup.space.mesh;
    let cell = (0..mesh.n_cells())
        .find(|&cell| mesh.coords[mesh.cells[2 * cell]] <= x && x <= mesh.coords[mesh.cells[2 * cell + 1]])
        .ok_or("the source lies outside the structure")?;
    let xi = 2.0 * (x - mesh.coords[mesh.cells[2 * cell]]) / setup.space.cell_length(cell) - 1.0;
    let mut cell_dofs = vec![0; setup.space.local_size()];
    setup.space.cell_dofs(cell, &mut cell_dofs);
    for (&dof, value) in cell_dofs.iter().zip(setup.space.basis(xi).0) {
        if let Some((row, factor)) = setup.dofs.map[dof] {
            rhs[row] += factor.conj() * amp * value;
        }
    }
    Ok(())
}
