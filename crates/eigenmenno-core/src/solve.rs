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

use crate::assemble::{assemble, Csr, DofMap};
use crate::banded::BandedLu;
use crate::eigen::eigs_near;
use crate::mesh::Mesh;
use crate::problem::{Analysis, Boundary, Layer, Polarization, Request, Source};
use crate::space::{H1Line, Space};
use crate::{C64, I, ONE, ZERO};
use serde::Serialize;
use std::f64::consts::PI;

const EIG_TOL: f64 = 1e-11;
/// Field output cap, so huge meshes cannot flood the plot.
const MAX_FIELD_POINTS: usize = 20_000;
/// Refuses meshes that would exhaust a browser tab's memory.
const MAX_DOFS: f64 = 2.0e6;

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

fn c(v: [f64; 2]) -> C64 {
    C64::new(v[0], v[1])
}

/// `(p, q)` of a layer.
fn layer_coef(pol: Polarization, l: &Layer) -> (C64, C64) {
    match pol {
        Polarization::Te => (c(l.mu).inv(), c(l.eps)),
        Polarization::Tm => (c(l.eps).inv(), c(l.mu)),
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

    /// The matrix of `a S + b M`, with `f(p, q) = (a, b)` per cell.
    fn matrix(&self, f: impl Fn(C64, C64) -> (C64, C64)) -> Csr {
        let t = assemble(&self.space, &self.dofs, |cell| {
            let (p, q) = self.coef(cell);
            f(p, q)
        });
        Csr::from_triplets(self.dofs.n_free, t)
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
        let mut data = vec![vec![]; 6];
        let scale = if k0 > 0.0 { 1.0 / k0 } else { 0.0 };
        for cell in (0..mesh.n_cells()).step_by(stride) {
            let layer = &self.req.layers[mesh.region[cell]];
            let (x0, h) = (mesh.coords[mesh.cells[2 * cell]], self.space.cell_length(cell));
            // Each cell repeats its end points, so jumps at interfaces plot as jumps.
            for s in 0..per_cell {
                let xi = -1.0 + 2.0 * s as f64 / (per_cell - 1) as f64;
                let (v, d) = self.space.eval(cell, xi, u);
                let (second, third) = if te {
                    (-v * kt * scale / c(layer.mu), -I * d * scale / c(layer.mu))
                } else {
                    (v * kt * scale / c(layer.eps), I * d * scale / c(layer.eps))
                };
                x.push(x0 + h * (xi + 1.0) / 2.0);
                for (k, w) in [v, second, third].into_iter().enumerate() {
                    data[2 * k].push(w.re);
                    data[2 * k + 1].push(w.im);
                }
            }
        }
        let components = if te { ["Ey", "Hx", "Hz"] } else { ["Hy", "Ex", "Ez"] };
        let i = out.meta.fields.len();
        out.channels.push((format!("field{i}/x"), x));
        for (k, d) in data.into_iter().enumerate() {
            out.channels.push((format!("field{i}/{}/{}", components[k / 2], ["re", "im"][k % 2]), d));
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

/// Scales a mode to a unit real peak, so plots are comparable.
pub(crate) fn normalized(mut u: Vec<C64>) -> Vec<C64> {
    let peak = u.iter().fold(ZERO, |m, x| if x.norm() > m.norm() { *x } else { m });
    if peak != ZERO {
        u.iter_mut().for_each(|x| *x /= peak);
    }
    u
}

fn validate(req: &Request) -> Result<(), String> {
    let err = |m: &str| Err(m.to_string());
    let index = |l: &Layer| (c(l.eps) * c(l.mu)).norm().sqrt().max(1.0);
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

pub fn solve(req: &Request) -> Result<Output, String> {
    validate(req)?;
    let thickness: Vec<f64> = req.layers.iter().map(|l| l.thickness).collect();
    // Scale cells by local index so each layer has the same cells per wavelength.
    let max_cell: Vec<f64> =
        req.layers.iter().map(|l| req.mesh.max_cell / (c(l.eps) * c(l.mu)).norm().sqrt().max(1.0)).collect();
    let mesh = Mesh::layered_line(&thickness, &max_cell);
    let length: f64 = thickness.iter().sum();
    let bloch_phase = if let Analysis::Cavity { bloch_phase, .. } = req.analysis { bloch_phase } else { 0.0 };
    let s = Setup::new(req, &mesh, bloch_phase);
    let meta = Meta { length, n_dofs: s.dofs.n_free, ..Default::default() };
    let mut out = Output { meta, channels: vec![] };
    let k0_of = |wavelength: f64| 2.0 * PI / wavelength;
    match &req.analysis {
        Analysis::Cavity { n_modes, target_k0, .. } => {
            for (lambda, x) in cavity_modes(&s, *n_modes, *target_k0, length)? {
                s.push_field(
                    &mut out,
                    &normalized(s.dofs.expand(&x)),
                    lambda.max(0.0).sqrt(),
                    0.0,
                    None,
                    None,
                );
            }
        }
        Analysis::Bands { n_bands, n_k } => bands(&s, *n_bands, *n_k, &mut out)?,
        Analysis::Waveguide { wavelength, n_modes } => waveguide(&s, k0_of(*wavelength), *n_modes, &mut out)?,
        Analysis::Driven { wavelength, angle_deg, source } => {
            let k0 = k0_of(*wavelength);
            let d = driven(&s, k0, *angle_deg, *source)?;
            s.push_field(&mut out, &d.u, k0, d.kt, None, d.power);
        }
        Analysis::DrivenGradient { wavelength, angle_deg } => {
            let k0 = k0_of(*wavelength);
            let d = driven(&s, k0, *angle_deg, Source::PlaneWave)?;
            s.push_field(&mut out, &d.u, k0, d.kt, None, d.power);
            power_gradients(&s, k0, &d, &mut out);
        }
        Analysis::Sweep { wavelengths, angle_deg } => {
            let mut series = [vec![], vec![], vec![]];
            for &w in wavelengths {
                let power = driven(&s, k0_of(w), *angle_deg, Source::PlaneWave)?.power.unwrap();
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
    s: &Setup,
    n: usize,
    target_k0: Option<f64>,
    length: f64,
) -> Result<Vec<(f64, Vec<C64>)>, String> {
    let a = s.matrix(|p, _| (p, ZERO));
    let b = s.matrix(|_, q| (ZERO, q));
    // A negative shift keeps the operator definite when searching up from k0 = 0.
    let sigma = target_k0.map_or(-(PI / length).powi(2), |k| k * k);
    eigs_near(&a, &b, sigma, n, EIG_TOL)
}

fn bands(s: &Setup, n_bands: usize, n_k: usize, out: &mut Output) -> Result<(), String> {
    let length = out.meta.length;
    let ks: Vec<f64> = (0..n_k).map(|i| PI / length * i as f64 / (n_k - 1) as f64).collect();
    let mut table = vec![vec![]; n_bands];
    for &k in &ks {
        let s = Setup::new(s.req, s.space.mesh, k * length);
        for (band, (lambda, _)) in table.iter_mut().zip(cavity_modes(&s, n_bands, None, length)?) {
            band.push(lambda.max(0.0).sqrt());
        }
    }
    out.channels.push(("k".into(), ks));
    out.channels.extend(table.into_iter().enumerate().map(|(i, b)| (format!("band{i}"), b)));
    Ok(())
}

/// `(k0^2 Q - S) u = beta^2 P u`, with `P` the mass matrix weighted by `p`.
fn waveguide(s: &Setup, k0: f64, n_modes: usize, out: &mut Output) -> Result<(), String> {
    let k02 = C64::new(k0 * k0, 0.0);
    let a = s.matrix(|p, q| (-p, k02 * q));
    let b = s.matrix(|p, _| (ZERO, p));
    let n2: Vec<f64> = s.req.layers.iter().map(|l| l.eps[0] * l.mu[0]).collect();
    let n2_max = n2.iter().copied().fold(0.0, f64::max);
    let n2_clad = n2[0].max(n2[n2.len() - 1]);
    // Highest effective index first. Modes at or below the cladding index are box
    // modes of the truncated domain, not guided modes, so they are dropped.
    for (lambda, x) in eigs_near(&a, &b, k0 * k0 * n2_max * 1.001, n_modes, EIG_TOL)?.into_iter().rev() {
        if lambda > k0 * k0 * n2_clad {
            let beta = lambda.sqrt();
            s.push_field(out, &normalized(s.dofs.expand(&x)), k0, beta, Some(beta / k0), None);
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

fn driven(s: &Setup, k0: f64, angle_deg: f64, source: Source) -> Result<Driven, String> {
    let req = s.req;
    let n2 = |l: &Layer| c(l.eps) * c(l.mu);
    let plane = matches!(source, Source::PlaneWave);
    let (first, last_layer) = (&req.layers[0], &req.layers[req.layers.len() - 1]);
    let kt = if plane { k0 * n2(first).sqrt().re * angle_deg.to_radians().sin() } else { 0.0 };
    let (k2, kt2) = (C64::new(k0 * k0, 0.0), C64::new(kt * kt, 0.0));
    let mut t = assemble(&s.space, &s.dofs, |cell| {
        let (p, q) = s.coef(cell);
        (p, kt2 * p - k2 * q)
    });
    // `p kx` at an end, with `kx` on the outgoing or decaying branch. An open end
    // has `u' = i kx u`, which turns the weak form's boundary term into this entry.
    let admittance = |l: &Layer| {
        let w = (k2 * n2(l) - kt2).sqrt();
        layer_coef(req.polarization, l).0 * if w.im < 0.0 || (w.im == 0.0 && w.re < 0.0) { -w } else { w }
    };
    let last = s.space.n_dofs() - 1;
    let (y_left, y_right) = (admittance(first), admittance(last_layer));
    for (is_open, dof, y) in [(req.left.is_open(), 0, y_left), (req.right.is_open(), last, y_right)] {
        if is_open {
            let (d, _) = s.dofs.map[dof].unwrap();
            t.push((d, d, -I * y));
        }
    }
    let mut rhs = vec![ZERO; s.dofs.n_free];
    match source {
        Source::PlaneWave => rhs[s.dofs.map[0].unwrap().0] = -I * y_left * 2.0,
        Source::Line { x } => add_point_load(s, x, I * k0, &mut rhs)?,
    }
    let lu = BandedLu::factor(&Csr::from_triplets(s.dofs.n_free, t))?;
    lu.solve(&mut rhs);
    let u = s.dofs.expand(&rhs);
    let t_scale = if req.right.is_open() { y_right.re / y_left.re } else { 0.0 };
    let power = plane.then(|| {
        let r = (u[0] - ONE).norm_sqr();
        let tr = u[last].norm_sqr() * t_scale;
        [r, tr, 1.0 - r - tr]
    });
    Ok(Driven { u, kt, power, lu, t_scale })
}

/// Adjoint gradients of R, T and A. For a real `J(u)` with `A u = b`,
/// `dJ/dθ = -2 Re(λᵀ (dA/dθ) u)`, where `Aᵀ λ = dJ/du`. Without periodic walls `A`
/// is complex symmetric, so the factors of `A` also solve for `λ`, and each output
/// costs one back-substitution however many layers there are.
///
/// Layers at an open end are left out: they also set the angle and the boundary
/// terms. The result is the exact gradient of the discrete problem on this mesh.
fn power_gradients(s: &Setup, k0: f64, d: &Driven, out: &mut Output) {
    let last = s.space.n_dofs() - 1;
    let adjoint = |dof: usize, g: C64| {
        let mut x = vec![ZERO; s.dofs.n_free];
        if let Some((r, _)) = s.dofs.map[dof] {
            x[r] = g;
        }
        d.lu.solve(&mut x);
        s.dofs.expand(&x)
    };
    let dr = layer_sensitivity(s, k0, d, &adjoint(0, (d.u[0] - ONE).conj()));
    let dt = layer_sensitivity(s, k0, d, &adjoint(last, d.u[last].conj() * d.t_scale));
    let da: Vec<C64> = dr.iter().zip(&dt).map(|(r, t)| -r - t).collect();
    let n = s.req.layers.len();
    let design = 1..if s.req.right.is_open() { n - 1 } else { n };
    out.channels.push(("layer".into(), design.clone().map(|l| l as f64).collect()));
    for (name, w) in [("R", &dr), ("T", &dt), ("A", &da)] {
        // `dε = i dε''` multiplies `λᵀ (dA/dε) u` by `i`, which turns Re into -Im.
        out.channels.push((format!("d{name}/deps_re"), design.clone().map(|l| -2.0 * w[l].re).collect()));
        out.channels.push((format!("d{name}/deps_im"), design.clone().map(|l| 2.0 * w[l].im).collect()));
    }
}

/// `λᵀ (dA/dε) u` for each layer. The derivative has the usual `a S + b M` form:
/// TE has `q = ε` in `b = kt^2 p - k0^2 q`, and TM has `p = 1/ε` in both `a` and `b`.
fn layer_sensitivity(s: &Setup, k0: f64, d: &Driven, lambda: &[C64]) -> Vec<C64> {
    let n = s.space.local_size();
    let (mut stiff, mut unused, mut mass, mut dofs) =
        (vec![0.0; n * n], vec![0.0; n * n], vec![0.0; n * n], vec![0; n]);
    let mut w = vec![ZERO; s.req.layers.len()];
    for cell in 0..s.space.n_cells() {
        s.space.cell_dofs(cell, &mut dofs);
        s.space.local_matrices(cell, &mut stiff, &mut unused, &mut mass);
        let (p, _) = s.coef(cell);
        let (a, b) = match s.req.polarization {
            Polarization::Te => (ZERO, C64::new(-k0 * k0, 0.0)),
            Polarization::Tm => (-p * p, -p * p * (d.kt * d.kt)),
        };
        for i in 0..n {
            for j in 0..n {
                w[s.space.mesh.region[cell]] +=
                    lambda[dofs[i]] * (a * stiff[i * n + j] + b * mass[i * n + j]) * d.u[dofs[j]];
            }
        }
    }
    w
}

/// Adds `amp * v(x)` to the right-hand side: the weak form of a current sheet.
fn add_point_load(s: &Setup, x: f64, amp: C64, rhs: &mut [C64]) -> Result<(), String> {
    let mesh = s.space.mesh;
    let cell = (0..mesh.n_cells())
        .find(|&c| mesh.coords[mesh.cells[2 * c]] <= x && x <= mesh.coords[mesh.cells[2 * c + 1]])
        .ok_or("the source lies outside the structure")?;
    let xi = 2.0 * (x - mesh.coords[mesh.cells[2 * cell]]) / s.space.cell_length(cell) - 1.0;
    let mut dofs = vec![0; s.space.local_size()];
    s.space.cell_dofs(cell, &mut dofs);
    for (&g, v) in dofs.iter().zip(s.space.basis(xi).0) {
        if let Some((r, f)) = s.dofs.map[g] {
            rhs[r] += f.conj() * amp * v;
        }
    }
    Ok(())
}
