//! Every test compares against an exact result, not against stored numbers.

use eigenmenno_core::mesh::{Mesh, Outline};
use eigenmenno_core::problem::{Analysis, Analysis2d, Boundary, Density, Layer, Material, MeshSpec, Polarization, Request, Request2d, Shape, Source, Source2d, Monitor};
use eigenmenno_core::solve::{solve, Output};
use eigenmenno_core::C64;
use std::f64::consts::PI;

fn layer(thickness: f64, n: f64) -> Layer {
    Layer { thickness, eps: [n * n, 0.0], mu: [1.0, 0.0] }
}

fn request(layers: Vec<Layer>, pol: Polarization, b: (Boundary, Boundary), order: usize, h: f64, analysis: Analysis) -> Request {
    Request { layers, polarization: pol, left: b.0, right: b.1, mesh: MeshSpec { order, max_cell: h }, analysis }
}

fn channel<'a>(o: &'a Output, name: &str) -> &'a [f64] {
    &o.channels.iter().find(|(n, _)| n == name).unwrap_or_else(|| panic!("no channel {name}")).1
}

fn close(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a} differs from {b} by more than {tol}");
}

fn cavity(n_modes: usize) -> Analysis {
    Analysis::Cavity { n_modes, target_k0: None, bloch_phase: 0.0 }
}

#[test]
fn pec_cavity_matches_mode_formula() {
    let r = request(vec![layer(2.0, 1.0)], Polarization::Te, (Boundary::Pec, Boundary::Pec), 4, 0.05, cavity(6));
    let o = solve(&r).unwrap();
    assert_eq!(o.meta.fields.len(), 6);
    for (m, f) in o.meta.fields.iter().enumerate() {
        close(f.k0, (m + 1) as f64 * PI / 2.0, 1e-8);
    }
}

#[test]
fn eigenvalue_error_converges_at_rate_2p() {
    for order in 1..=4 {
        let err = |h: f64| {
            let r = request(vec![layer(1.0, 1.0)], Polarization::Te, (Boundary::Pec, Boundary::Pec), order, h, cavity(1));
            let k = solve(&r).unwrap().meta.fields[0].k0;
            (k * k - PI * PI).abs()
        };
        // Coarse enough that high orders stay above round-off.
        let (e1, e2) = (err(0.5), err(0.25));
        let rate = (e1 / e2).log2();
        let expected = 2.0 * order as f64;
        assert!((rate - expected).abs() < 0.6, "order {order}: rate {rate}, expected {expected}");
    }
}

#[test]
fn tm_polarization_swaps_pec_and_pmc() {
    // A TM field with PEC walls has a free H_y at the walls, including k0 = 0.
    let r = request(vec![layer(1.0, 1.0)], Polarization::Tm, (Boundary::Pec, Boundary::Pec), 3, 0.05, cavity(3));
    let o = solve(&r).unwrap();
    // Compare k0^2: it is what gets solved for, and a square root would amplify
    // round-off in the zero eigenvalue.
    for (m, f) in o.meta.fields.iter().enumerate() {
        close(f.k0 * f.k0, (m as f64 * PI).powi(2), 1e-6);
    }
}

#[test]
fn layered_cavity_satisfies_transcendental_equation() {
    let (d1, n1, d2, n2) = (0.7, 1.0, 0.4, 2.5);
    let r = request(vec![layer(d1, n1), layer(d2, n2)], Polarization::Te, (Boundary::Pec, Boundary::Pec), 5, 0.05, cavity(5));
    for f in &solve(&r).unwrap().meta.fields {
        let (a, b) = (f.k0 * n1 * d1, f.k0 * n2 * d2);
        // Continuity of E and H at the interface, cleared of poles.
        let res = n2 * a.sin() * b.cos() + n1 * a.cos() * b.sin();
        close(res, 0.0, 1e-8);
    }
}

#[test]
fn bloch_modes_of_uniform_medium_include_degenerate_pairs() {
    // Zone centre of a uniform cell: k0 = 0 once, then 2 pi m / L twice each.
    let r = request(vec![layer(1.0, 1.0)], Polarization::Te, (Boundary::Periodic, Boundary::Periodic), 4, 0.05, cavity(5));
    let o = solve(&r).unwrap();
    let want = [0.0, 2.0 * PI, 2.0 * PI, 4.0 * PI, 4.0 * PI];
    for (f, w) in o.meta.fields.iter().zip(want) {
        close(f.k0, w, 1e-6);
    }
}

#[test]
fn bloch_bands_satisfy_bragg_stack_dispersion() {
    let (d1, n1, d2, n2) = (0.5, 1.0, 0.5, 3.0);
    let r = request(
        vec![layer(d1, n1), layer(d2, n2)],
        Polarization::Te,
        (Boundary::Periodic, Boundary::Periodic),
        4,
        0.05,
        Analysis::Bands { n_bands: 4, n_k: 9 },
    );
    let o = solve(&r).unwrap();
    let ks = channel(&o, "k");
    for b in 0..4 {
        for (i, &k0) in channel(&o, &format!("band{b}")).iter().enumerate() {
            let (a, c) = (k0 * n1 * d1, k0 * n2 * d2);
            let rhs = a.cos() * c.cos() - 0.5 * (n1 / n2 + n2 / n1) * a.sin() * c.sin();
            close((ks[i] * (d1 + d2)).cos(), rhs, 1e-6);
        }
    }
}

#[test]
fn slab_waveguide_te_and_tm_satisfy_dispersion() {
    let (n_core, n_clad, d, wl) = (2.0, 1.0, 0.6, 1.0);
    for pol in [Polarization::Te, Polarization::Tm] {
        let r = request(
            // Wide cladding: the third mode is near cutoff and decays slowly, so the
            // PEC truncation error scales as exp(-2 gamma * cladding).
            vec![layer(10.0, n_clad), layer(d, n_core), layer(10.0, n_clad)],
            pol,
            (Boundary::Pec, Boundary::Pec),
            5,
            0.05,
            Analysis::Waveguide { wavelength: wl, n_modes: 3 },
        );
        let o = solve(&r).unwrap();
        assert_eq!(o.meta.fields.len(), 3, "a slab this thick guides three modes");
        let k0 = 2.0 * PI / wl;
        let ratio = if pol == Polarization::Te { 1.0 } else { (n_core / n_clad).powi(2) };
        for (m, f) in o.meta.fields.iter().enumerate() {
            let beta = f.n_eff.unwrap() * k0;
            let kappa = (k0 * k0 * n_core * n_core - beta * beta).sqrt();
            let gamma = (beta * beta - k0 * k0 * n_clad * n_clad).sqrt();
            // kappa d = m pi + 2 atan(ratio gamma / kappa)
            let res = kappa * d - m as f64 * PI - 2.0 * (ratio * gamma / kappa).atan();
            // Modes near cutoff feel the PEC walls through their evanescent tails.
            close(res, 0.0, 1e-6 + 5.0 * (-2.0 * gamma * 10.0).exp());
        }
    }
}

/// Independent transfer-matrix reference for layered reflection.
fn transfer_matrix(layers: &[Layer], pol: Polarization, k0: f64, angle: f64) -> (f64, f64) {
    let n2: Vec<C64> = layers.iter().map(|l| C64::new(l.eps[0], l.eps[1]) * C64::new(l.mu[0], l.mu[1])).collect();
    let kt = k0 * n2[0].sqrt().re * angle.sin();
    let kx = |i: usize| {
        let w = (n2[i] * k0 * k0 - kt * kt).sqrt();
        if w.im < 0.0 { -w } else { w }
    };
    let p = |i: usize| match pol {
        Polarization::Te => C64::new(1.0, 0.0) / C64::new(layers[i].mu[0], layers[i].mu[1]),
        Polarization::Tm => C64::new(1.0, 0.0) / C64::new(layers[i].eps[0], layers[i].eps[1]),
    };
    let n = layers.len();
    let i = C64::new(0.0, 1.0);
    // March (u, p u') from the transmitted side, where u = 1.
    let (mut u, mut f) = (C64::new(1.0, 0.0), i * p(n - 1) * kx(n - 1));
    for l in (1..n - 1).rev() {
        let (k, d) = (kx(l), layers[l].thickness);
        let y = p(l) * k;
        let (c, s) = ((k * d).cos(), (k * d).sin());
        let (nu, nf) = (u * c - f * s / y, f * c + u * y * s);
        u = nu;
        f = nf;
    }
    let y0 = i * p(0) * kx(0);
    let inc = (u + f / y0) / 2.0;
    let refl = (u - f / y0) / 2.0;
    let r = (refl / inc).norm_sqr();
    let t = (C64::new(1.0, 0.0) / inc).norm_sqr() * (p(n - 1) * kx(n - 1)).re / (p(0) * kx(0)).re;
    (r, t)
}

#[test]
fn single_interface_matches_fresnel_at_oblique_incidence() {
    let (n1, n2, wl) = (1.0, 1.5, 1.0);
    let th = 40f64.to_radians();
    let th2 = (n1 * th.sin() / n2).asin();
    for (pol, rs) in [
        (Polarization::Te, (n1 * th.cos() - n2 * th2.cos()) / (n1 * th.cos() + n2 * th2.cos())),
        (Polarization::Tm, (n2 * th.cos() - n1 * th2.cos()) / (n2 * th.cos() + n1 * th2.cos())),
    ] {
        let r = request(
            vec![layer(0.5, n1), layer(0.5, n2)],
            pol,
            (Boundary::Port, Boundary::Absorbing),
            5,
            0.05,
            Analysis::Driven { wavelength: wl, angle_deg: 40.0, source: Source::PlaneWave },
        );
        let o = solve(&r).unwrap();
        let f = &o.meta.fields[0];
        close(f.power.unwrap()[0], rs * rs, 1e-7);
        close(f.power.unwrap()[1], 1.0 - rs * rs, 1e-7);
    }
}

#[test]
fn multilayer_spectrum_matches_transfer_matrix() {
    let layers = vec![layer(0.3, 1.0), layer(0.17, 2.3), layer(0.27, 1.38), layer(0.17, 2.3), layer(0.3, 1.5)];
    for (pol, angle) in [(Polarization::Te, 0.0), (Polarization::Te, 35.0), (Polarization::Tm, 35.0)] {
        let wavelengths: Vec<f64> = (0..21).map(|i| 0.5 + 0.05 * i as f64).collect();
        let r = request(
            layers.clone(),
            pol,
            (Boundary::Port, Boundary::Absorbing),
            5,
            0.05,
            Analysis::Sweep { wavelengths: wavelengths.clone(), angle_deg: angle },
        );
        let o = solve(&r).unwrap();
        for (i, &wl) in wavelengths.iter().enumerate() {
            let (rr, tt) = transfer_matrix(&layers, pol, 2.0 * PI / wl, f64::to_radians(angle));
            close(channel(&o, "R")[i], rr, 1e-7);
            close(channel(&o, "T")[i], tt, 1e-7);
        }
    }
}

#[test]
fn lossy_layer_absorbs_and_matches_transfer_matrix() {
    let mut lossy = layer(0.3, 1.0);
    lossy.eps = [4.0, 0.8];
    let layers = vec![layer(0.3, 1.0), lossy, layer(0.3, 1.0)];
    let r = request(
        layers.clone(),
        Polarization::Te,
        (Boundary::Port, Boundary::Absorbing),
        5,
        0.05,
        Analysis::Driven { wavelength: 1.0, angle_deg: 0.0, source: Source::PlaneWave },
    );
    let o = solve(&r).unwrap();
    let f = &o.meta.fields[0];
    let (rr, tt) = transfer_matrix(&layers, Polarization::Te, 2.0 * PI, 0.0);
    close(f.power.unwrap()[0], rr, 1e-7);
    close(f.power.unwrap()[1], tt, 1e-7);
    assert!(f.power.unwrap()[2] > 0.1, "the layer should absorb");
}

fn gradient(layers: &[Layer], pol: Polarization, right: Boundary, order: usize, h: f64, angle: f64) -> Output {
    let analysis = Analysis::DrivenGradient { wavelength: 1.0, angle_deg: angle };
    solve(&request(layers.to_vec(), pol, (Boundary::Port, right), order, h, analysis)).unwrap()
}

/// Central differences of `power(layers)` over each `eps` listed in `layer`,
/// compared with the adjoint channels. `power` returns `[R, T, A]`.
fn check_gradient(o: &Output, layers: &[Layer], tol: f64, power: impl Fn(&[Layer]) -> [f64; 3]) {
    const STEP: f64 = 1e-5;
    for (row, &l) in channel(o, "layer").iter().enumerate() {
        for (part, suffix) in ["re", "im"].into_iter().enumerate() {
            let at = |delta: f64| {
                let mut moved = layers.to_vec();
                moved[l as usize].eps[part] += delta;
                power(&moved)
            };
            let (up, down) = (at(STEP), at(-STEP));
            for (i, name) in ["R", "T", "A"].into_iter().enumerate() {
                let want = (up[i] - down[i]) / (2.0 * STEP);
                close(channel(o, &format!("d{name}/deps_{suffix}"))[row], want, tol);
            }
        }
    }
}

#[test]
fn adjoint_gradient_matches_transfer_matrix_differences() {
    let mut lossy = layer(0.2, 1.0);
    lossy.eps = [4.0, 0.8];
    let layers = vec![layer(0.3, 1.0), layer(0.17, 2.3), lossy, layer(0.27, 1.38), layer(0.3, 1.5)];
    for (pol, angle) in [(Polarization::Te, 0.0), (Polarization::Te, 35.0), (Polarization::Tm, 35.0)] {
        let o = gradient(&layers, pol, Boundary::Absorbing, 5, 0.05, angle);
        assert_eq!(channel(&o, "layer"), [1.0, 2.0, 3.0]);
        assert_eq!(o.meta.fields.len(), 1);
        check_gradient(&o, &layers, 1e-5, |moved| {
            let (r, t) = transfer_matrix(moved, pol, 2.0 * PI, f64::to_radians(angle));
            [r, t, 1.0 - r - t]
        });
    }
}

#[test]
fn adjoint_gradient_is_exact_on_a_coarse_mesh() {
    let parsed: Analysis = serde_json::from_str(r#"{"type": "driven_gradient", "wavelength": 1.0}"#).unwrap();
    assert!(matches!(parsed, Analysis::DrivenGradient { .. }));
    // The mesh is coarse enough to differ from the exact answer, so this only passes
    // if the adjoint differentiates the discrete problem itself. The small steps
    // never change a layer's cell count.
    let mut lossy = layer(0.21, 1.0);
    lossy.eps = [4.0, 0.8];
    let layers = vec![layer(0.33, 1.0), layer(0.17, 2.3), lossy, layer(0.27, 1.38)];
    let cases = [
        (Polarization::Te, Boundary::Absorbing, vec![1.0, 2.0]),
        (Polarization::Tm, Boundary::Absorbing, vec![1.0, 2.0]),
        // A closed right end leaves the last layer free to change, and T at zero.
        (Polarization::Te, Boundary::Pec, vec![1.0, 2.0, 3.0]),
        (Polarization::Tm, Boundary::Pec, vec![1.0, 2.0, 3.0]),
    ];
    for (pol, right, design) in cases {
        let o = gradient(&layers, pol, right, 2, 0.13, 25.0);
        assert_eq!(channel(&o, "layer"), design);
        check_gradient(&o, &layers, 1e-6, |moved| {
            let analysis = Analysis::Driven { wavelength: 1.0, angle_deg: 25.0, source: Source::PlaneWave };
            let r = request(moved.to_vec(), pol, (Boundary::Port, right), 2, 0.13, analysis);
            solve(&r).unwrap().meta.fields[0].power.unwrap()
        });
    }
}

#[test]
fn line_source_radiates_the_free_space_green_function() {
    let r = request(
        vec![layer(2.0, 1.0)],
        Polarization::Te,
        (Boundary::Absorbing, Boundary::Absorbing),
        5,
        0.05,
        Analysis::Driven { wavelength: 1.0, angle_deg: 0.0, source: Source::Line { x: 1.0 } },
    );
    let o = solve(&r).unwrap();
    let k = 2.0 * PI;
    let (x, re, im) = (channel(&o, "field0/x"), channel(&o, "field0/Ey/re"), channel(&o, "field0/Ey/im"));
    for i in 0..x.len() {
        // -u'' - k^2 u = i k delta(x - x0) has u = -(1/2) e^{ik|x - x0|}. The field
        // opposes the current at the sheet, so the source delivers power.
        let exact = C64::new(0.0, k * (x[i] - 1.0).abs()).exp() * -0.5;
        close(re[i], exact.re, 1e-5);
        close(im[i], exact.im, 1e-5);
    }
}

#[test]
fn derived_fields_satisfy_faraday_in_a_cavity() {
    // Fundamental PEC mode: E_y = sin(pi x), so eta0 H_z = -(i/k0) E_y' = -i cos(pi x).
    let r = request(vec![layer(1.0, 1.0)], Polarization::Te, (Boundary::Pec, Boundary::Pec), 5, 0.05, cavity(1));
    let o = solve(&r).unwrap();
    let (x, ey, hz) = (channel(&o, "field0/x"), channel(&o, "field0/Ey/re"), channel(&o, "field0/Hz/im"));
    let sign = ey[x.len() / 2].signum();
    for i in 0..x.len() {
        close(ey[i] * sign, (PI * x[i]).sin(), 1e-5);
        close(hz[i] * sign, -(PI * x[i]).cos(), 1e-5);
    }
}

#[test]
fn invalid_requests_return_errors() {
    let bad = request(vec![layer(1.0, 1.0)], Polarization::Te, (Boundary::Periodic, Boundary::Pec), 1, 0.1, cavity(1));
    assert!(solve(&bad).is_err());
    let open = request(vec![layer(1.0, 1.0)], Polarization::Te, (Boundary::Port, Boundary::Absorbing), 1, 0.1, cavity(1));
    assert!(solve(&open).is_err());
    let mut lossy = layer(1.0, 1.0);
    lossy.eps[1] = 0.1;
    let loss = request(vec![lossy], Polarization::Te, (Boundary::Pec, Boundary::Pec), 1, 0.1, cavity(1));
    assert!(solve(&loss).is_err());
}

#[test]
fn request_parses_from_json() {
    let json = r#"{"layers":[{"thickness":1.0,"eps":[4,0]}],"polarization":"te","left":"pec","right":"pec",
        "mesh":{"order":3,"max_cell":0.1},"analysis":{"type":"cavity","n_modes":2}}"#;
    let r: Request = serde_json::from_str(json).unwrap();
    close(solve(&r).unwrap().meta.fields[0].k0, PI / 2.0, 1e-6);
}

fn rods(cell: f64, pol: Polarization, analysis: Analysis2d) -> Request2d {
    Request2d {
        width: 1.0,
        height: 1.0,
        cell,
        background: Material { eps: 1.0, mu: 1.0 },
        shapes: vec![Shape { outline: Outline::Circle { x: 0.5, y: 0.5, r: 0.2 }, material: Material { eps: 12.0, mu: 1.0 } }],
        polarization: pol,
        walls: Boundary::Periodic,
        analysis,
        monitors: vec![],
        density: None,
    }
}

/// Frequency `a / lambda` of the lowest modes of the rod lattice at one Bloch phase.
fn rod_bands(cell: f64, pol: Polarization, bloch: [f64; 2]) -> Vec<f64> {
    let o = eigenmenno_core::solve2d::solve(&rods(cell, pol, Analysis2d::Cavity { n_modes: 3, bloch, target_k0: None })).unwrap();
    o.meta.fields.iter().map(|f| f.k0 / (2.0 * PI)).collect()
}

#[test]
fn fitted_mesh_has_no_slivers_and_follows_the_outline() {
    let circle = Outline::Circle { x: 0.5, y: 0.5, r: 0.2 };
    let rect = Outline::Rect { x0: 0.13, y0: 0.27, x1: 0.81, y1: 0.42 };
    let shifted = Outline::Circle { x: 0.47, y: 0.53, r: 0.31 };
    for (outlines, area) in [(vec![circle], PI * 0.04), (vec![shifted], PI * 0.0961), (vec![rect], 0.68 * 0.15)] {
        for n in [12, 24, 48, 96, 192] {
            let m = Mesh::fitted_grid(n, n, 1.0, 1.0, &outlines);
            let mut inside = 0.0;
            let mut worst = 180f64;
            for (c, v) in m.cells.chunks(3).enumerate() {
                let p: Vec<(f64, f64)> = v.iter().map(|&i| m.vertex(i)).collect();
                let det = (p[1].0 - p[0].0) * (p[2].1 - p[0].1) - (p[2].0 - p[0].0) * (p[1].1 - p[0].1);
                assert!(det > 0.0, "cell {c} is inverted on the {n} grid");
                if m.region[c] == 1 {
                    inside += det / 2.0;
                }
                for k in 0..3 {
                    let (o, a, b) = (p[k], p[(k + 1) % 3], p[(k + 2) % 3]);
                    let (ux, uy, vx, vy) = (a.0 - o.0, a.1 - o.1, b.0 - o.0, b.1 - o.1);
                    worst = worst.min((ux * vy - uy * vx).atan2(ux * vx + uy * vy).to_degrees());
                }
            }
            assert!(worst > 14.0, "smallest angle {worst} on the {n} grid");
            // A polygon inscribed in a circle loses area like h^2.
            close(inside, area, 3.0 / (n * n) as f64);
        }
    }
}

#[test]
fn rectangular_cavity_matches_mode_formula() {
    let box_of = |pol, walls| Request2d {
        width: 1.0,
        height: 1.5,
        cell: 0.05,
        background: Material { eps: 1.0, mu: 1.0 },
        shapes: vec![],
        polarization: pol,
        walls,
        analysis: Analysis2d::Cavity { n_modes: 6, bloch: [0.0, 0.0], target_k0: None },
        monitors: vec![],
        density: None,
    };
    let k2 = |m: f64, n: f64| PI * PI * (m * m + n * n / 2.25);
    // E_z vanishes on a perfect conductor, so the mode numbers start from 1.
    let ez = eigenmenno_core::solve2d::solve(&box_of(Polarization::Te, Boundary::Pec)).unwrap();
    let want = [k2(1.0, 1.0), k2(1.0, 2.0), k2(2.0, 1.0), k2(1.0, 3.0), k2(2.0, 2.0), k2(2.0, 3.0)];
    for (f, w) in ez.meta.fields.iter().zip(want) {
        close(f.k0 * f.k0, w, 2e-4 * w);
    }
    // H_z has zero normal derivative there, so a static mode comes first.
    let hz = eigenmenno_core::solve2d::solve(&box_of(Polarization::Tm, Boundary::Pec)).unwrap();
    let want = [0.0, k2(0.0, 1.0), k2(1.0, 0.0), k2(1.0, 1.0), k2(0.0, 2.0), k2(1.0, 2.0)];
    for (f, w) in hz.meta.fields.iter().zip(want) {
        close(f.k0 * f.k0, w, 2e-4 * w.max(1.0));
    }
}

#[test]
fn layered_box_matches_the_1d_solver() {
    // With magnetic walls top and bottom, E_z modes with no y variation are exactly
    // the 1D problem's. The first few modes of this short box have none.
    let dielectric = Material { eps: 6.25, mu: 1.0 };
    let r = Request2d {
        width: 1.1,
        height: 0.1,
        cell: 0.02,
        background: Material { eps: 1.0, mu: 1.0 },
        shapes: vec![Shape { outline: Outline::Rect { x0: 0.7, y0: -1.0, x1: 2.0, y1: 2.0 }, material: dielectric }],
        polarization: Polarization::Te,
        walls: Boundary::Pmc,
        analysis: Analysis2d::Cavity { n_modes: 4, bloch: [0.0, 0.0], target_k0: None },
        monitors: vec![],
        density: None,
    };
    let o = request(vec![layer(0.7, 1.0), layer(0.4, 2.5)], Polarization::Te, (Boundary::Pmc, Boundary::Pmc), 5, 0.02, cavity(4));
    let (two, one) = (eigenmenno_core::solve2d::solve(&r).unwrap(), solve(&o).unwrap());
    for (a, b) in two.meta.fields.iter().zip(&one.meta.fields) {
        close(a.k0, b.k0, 2e-3 * b.k0.max(1.0));
    }
}

#[test]
fn rod_lattice_matches_a_plane_wave_expansion() {
    // Square lattice, eps 12, r = 0.2a. Reference: an independent plane-wave expansion
    // with 2401 waves, using the inverse rule for H_z, whose X value was still rising by
    // about 4e-5 per step. See plans/2d-solver.md. Straight cells around the curved rod
    // limit this to second order, and 48 cells per period leave about 3e-4.
    let cell = 1.0 / 48.0;
    let (m, x) = ([PI, PI], [PI, 0.0]);
    close(rod_bands(cell, Polarization::Te, m)[0], 0.28066, 3e-4);
    close(rod_bands(cell, Polarization::Te, x)[1], 0.41716, 4e-4);
    close(rod_bands(cell, Polarization::Tm, m)[0], 0.49588, 6e-4);
    close(rod_bands(cell, Polarization::Tm, x)[1], 0.44109, 6e-4);
}

#[test]
fn bloch_phase_gives_a_wave_travelling_along_k() {
    // An empty lattice's lowest mode is the plane wave exp(i k.r), with phase k L per
    // period. So the sign and direction of the tie appear directly in the field.
    let (width, height) = (1.0, 2.0);
    let bloch = [PI / 2.0, PI / 4.0];
    let r = Request2d {
        width,
        height,
        cell: 0.05,
        background: Material { eps: 1.0, mu: 1.0 },
        shapes: vec![],
        polarization: Polarization::Te,
        walls: Boundary::Periodic,
        analysis: Analysis2d::Cavity { n_modes: 1, bloch, target_k0: None },
        monitors: vec![],
        density: None,
    };
    let o = eigenmenno_core::solve2d::solve(&r).unwrap();
    let (kx, ky) = (bloch[0] / width, bloch[1] / height);
    close(o.meta.fields[0].k0, kx.hypot(ky), 1e-5);
    let (x, y) = (channel(&o, "x"), channel(&o, "y"));
    let (re, im) = (channel(&o, "field0/Ez/re"), channel(&o, "field0/Ez/im"));
    let reference = C64::new(re[0], im[0]) / C64::from_polar(1.0, kx * x[0] + ky * y[0]);
    for i in 0..x.len() {
        let want = reference * C64::from_polar(1.0, kx * x[i] + ky * y[i]);
        assert!((C64::new(re[i], im[i]) - want).norm() < 1e-4, "dof {i}");
    }
}

#[test]
fn band_path_visits_the_lattice_symmetry_points() {
    let bands = eigenmenno_core::solve2d::solve(&rods(0.05, Polarization::Te, Analysis2d::Bands { n_bands: 2, n_k: 2 })).unwrap();
    let (s, lowest) = (channel(&bands, "s"), channel(&bands, "band0"));
    assert_eq!(s, [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]);
    // The lowest band starts and ends at zero at Gamma, and peaks at M.
    close(lowest[0], 0.0, 1e-6);
    close(lowest[6], 0.0, 1e-6);
    close(lowest[4] / (2.0 * PI), rod_bands(0.05, Polarization::Te, [PI, PI])[0], 1e-9);
}

#[test]
fn request_2d_parses_from_json() {
    let json = r#"{"width": 1, "height": 1, "cell": 0.1, "background": {"eps": 1}, "polarization": "te", "walls": "periodic",
        "shapes": [{"type": "circle", "x": 0.5, "y": 0.5, "r": 0.2, "material": {"eps": 12}}],
        "analysis": {"type": "bands", "n_bands": 4, "n_k": 3}}"#;
    let r: Request2d = serde_json::from_str(json).unwrap();
    assert!(eigenmenno_core::solve2d::solve(&r).is_ok());
}

/// `J_n(x)` and `Y_n(x)` for x > 0, from their integral representations.
fn bessel(n: i32, x: f64) -> (f64, f64) {
    let m = 2000;
    let (mut j, mut y) = (0.0, 0.0);
    for i in 0..m {
        let t = PI * (i as f64 + 0.5) / m as f64;
        j += (n as f64 * t - x * t.sin()).cos() / m as f64;
        y += (x * t.sin() - n as f64 * t).sin() / m as f64;
    }
    // The remaining part of Y_n, by Simpson's rule.
    let (top, steps) = ((80.0 / x).asinh(), 4000);
    let h = top / steps as f64;
    let sign = if n % 2 == 0 { 1.0 } else { -1.0 };
    let f = |t: f64| ((n as f64 * t).exp() + sign * (-n as f64 * t).exp()) * (-x * t.sinh()).exp();
    let tail: f64 = (0..=steps).map(|i| f(i as f64 * h) * if i == 0 || i == steps { 1.0 } else if i % 2 == 1 { 4.0 } else { 2.0 }).sum();
    (j, y - tail * h / (3.0 * PI))
}

fn hankel(n: i32, x: f64) -> C64 {
    let (j, y) = bessel(n, x);
    C64::new(j, y)
}

fn open_box(shapes: Vec<Shape>, pol: Polarization, cell: f64, pml: f64, source: Source2d) -> Request2d {
    Request2d {
        width: 3.0,
        height: 3.0,
        cell,
        background: Material { eps: 1.0, mu: 1.0 },
        shapes,
        polarization: pol,
        walls: Boundary::Absorbing,
        analysis: Analysis2d::Driven { wavelength: 1.0, source, pml, average_q: None },
        monitors: vec![],
        density: None,
    }
}

fn field_of(o: &Output, name: &str) -> Vec<(f64, f64, C64)> {
    let (x, y) = (channel(o, "x"), channel(o, "y"));
    let (re, im) = (channel(o, &format!("field0/{name}/re")), channel(o, &format!("field0/{name}/im")));
    (0..x.len()).map(|i| (x[i], y[i], C64::new(re[i], im[i]))).collect()
}

#[test]
fn line_source_radiates_the_2d_green_function() {
    let k = 2.0 * PI;
    for (cell, pml, tol) in [(1.0 / 12.0, 0.5, 0.01), (1.0 / 24.0, 1.0, 0.002)] {
        let r = open_box(vec![], Polarization::Te, cell, pml, Source2d::Line { x: 1.5, y: 1.5 });
        let o = eigenmenno_core::solve2d::solve(&r).unwrap();
        let mut worst = 0f64;
        for (x, y, u) in field_of(&o, "Ez") {
            let d = (x - 1.5).hypot(y - 1.5);
            // The corners are the hardest part of the layer, so go right out to them.
            if d >= 0.6 && (0.0..=3.0).contains(&x) && (0.0..=3.0).contains(&y) {
                // -(u'' ) - k^2 u = i k delta gives u = i k (i/4) H0(k r).
                let want = C64::new(0.0, k) * C64::new(0.0, 0.25) * hankel(0, k * d);
                worst = worst.max((u - want).norm() / want.norm());
            }
        }
        eprintln!("line source cell {cell:.3} pml {pml}: worst relative error {worst:.2e}");
        assert!(worst < tol, "error {worst}");
    }
}

/// The analytic field of a plane wave along +x on a cylinder at the origin.
fn mie(r: f64, theta: f64, k: f64, radius: f64, n1: f64, pol: Polarization) -> C64 {
    let k1 = k * n1;
    // Field and p du/dr are continuous, with p = 1 for E_z and 1/eps for H_z.
    let (pb, p1) = match pol {
        Polarization::Te => (1.0, 1.0),
        Polarization::Tm => (1.0, 1.0 / (n1 * n1)),
    };
    let diff = |n: i32, x: f64, z: &dyn Fn(i32, f64) -> C64| z(n - 1, x) - z(n, x) * (n as f64 / x);
    let (j, h) = (|n: i32, x: f64| C64::new(bessel(n, x).0, 0.0), |n: i32, x: f64| hankel(n, x));
    let (x, x1) = (k * radius, k1 * radius);
    // Outside, the incident wave is exact. Its series would need order past k r, while
    // the scattered part falls off quickly with order.
    let mut u = if r >= radius { C64::from_polar(1.0, k * r * theta.cos()) } else { C64::new(0.0, 0.0) };
    for n in -14..=14i32 {
        let m = n.abs();
        // Symmetry: order -n matches order n in these coefficients with this phase.
        let (jb, djb, hb, dhb) = (j(m, x), diff(m, x, &j), h(m, x), diff(m, x, &h));
        let (ja, dja) = (j(m, x1), diff(m, x1, &j));
        let num = dja * jb * (p1 * k1) - djb * ja * (pb * k);
        let den = dja * hb * (p1 * k1) - dhb * ja * (pb * k);
        let a = -num / den;
        let spin = C64::from_polar(1.0, m as f64 * PI / 2.0 + n as f64 * theta);
        if r >= radius {
            u += spin * a * h(m, k * r);
        } else {
            let c = (jb + a * hb) / ja;
            u += spin * c * j(m, k1 * r);
        }
    }
    u
}

#[test]
fn dielectric_cylinder_matches_the_mie_series() {
    let (k, radius, n1) = (2.0 * PI, 0.5, 2.0);
    let rod = || vec![Shape { outline: Outline::Circle { x: 1.5, y: 1.5, r: radius }, material: Material { eps: n1 * n1, mu: 1.0 } }];
    for (pol, name, angle) in [(Polarization::Te, "Ez", 0.0), (Polarization::Tm, "Hz", 0.0), (Polarization::Te, "Ez", 30.0f64)] {
        let r = open_box(rod(), pol, 1.0 / 32.0, 1.0, Source2d::PlaneWave { angle_deg: angle });
        let o = eigenmenno_core::solve2d::solve(&r).unwrap();
        let a = angle.to_radians();
        // The incident phase at the cylinder's centre, relative to the domain's origin.
        let centre = C64::from_polar(1.0, k * (1.5 * a.cos() + 1.5 * a.sin()));
        let mut worst = 0f64;
        for (x, y, u) in field_of(&o, name).into_iter().step_by(307) {
            let (dx, dy) = (x - 1.5, y - 1.5);
            let d = dx.hypot(dy);
            let inside = (0.0..=3.0).contains(&x) && (0.0..=3.0).contains(&y);
            // Skip the interface, where straight cells cut the curve.
            if !inside || d < 1e-3 || (d - radius).abs() < 0.04 {
                continue;
            }
            let want = centre * mie(d, dy.atan2(dx) - a, k, radius, n1, pol);
            worst = worst.max((u - want).norm());
        }
        eprintln!("mie {name} angle {angle}: worst absolute error {worst:.2e}");
        assert!(worst < 0.02, "error {worst}");
    }
}

#[test]
fn absorbing_requests_are_checked() {
    let rod = |x: f64| vec![Shape { outline: Outline::Circle { x, y: 1.5, r: 0.5 }, material: Material { eps: 4.0, mu: 1.0 } }];
    let wave = Source2d::PlaneWave { angle_deg: 0.0 };
    assert!(eigenmenno_core::solve2d::solve(&open_box(rod(0.2), Polarization::Te, 0.2, 0.5, wave)).is_err());
    assert!(eigenmenno_core::solve2d::solve(&open_box(vec![], Polarization::Te, 0.2, 0.5, Source2d::Line { x: 9.0, y: 1.0 })).is_err());
    assert!(eigenmenno_core::solve2d::solve(&open_box(vec![], Polarization::Te, 0.2, 0.0, wave)).is_err());
    let mut walled = open_box(vec![], Polarization::Te, 0.2, 0.5, wave);
    walled.walls = Boundary::Pec;
    assert!(eigenmenno_core::solve2d::solve(&walled).is_err());
}

/// A square lattice of rods (eps 12, r = 0.2a), `n` by `n`, with the middle one missing.
fn defect_crystal(n: usize, cell: f64, walls: Boundary, analysis: Analysis2d) -> Request2d {
    let middle = n / 2;
    let shapes = (0..n * n)
        .filter(|&i| i != middle * n + middle)
        .map(|i| {
            let (x, y) = ((i % n) as f64 + 0.5, (i / n) as f64 + 0.5);
            Shape { outline: Outline::Circle { x, y, r: 0.2 }, material: Material { eps: 12.0, mu: 1.0 } }
        })
        .collect();
    Request2d {
        width: n as f64,
        height: n as f64,
        cell,
        background: Material { eps: 1.0, mu: 1.0 },
        shapes,
        polarization: Polarization::Te,
        walls,
        analysis,
        monitors: vec![],
        density: None,
    }
}

#[test]
fn removing_a_rod_puts_one_mode_in_the_gap() {
    let target = 2.0 * PI * 0.35;
    let r = defect_crystal(5, 1.0 / 16.0, Boundary::Periodic, Analysis2d::Cavity { n_modes: 3, bloch: [0.0, 0.0], target_k0: Some(target) });
    let o = eigenmenno_core::solve2d::solve(&r).unwrap();
    let f: Vec<f64> = o.meta.fields.iter().map(|m| m.k0 / (2.0 * PI)).collect();
    eprintln!("defect supercell a/lambda {f:?}, areas {:?}", o.meta.fields.iter().map(|m| m.volume).collect::<Vec<_>>());
    assert_eq!(f.iter().filter(|&&v| v > 0.28066 && v < 0.41716).count(), 1, "{f:?}");
    // Reference: an independent plane-wave expansion of this supercell, with 4489 waves.
    // See plans/2d-solver.md. Straight cells around the rods keep this mesh above it.
    close(f[2], 0.3746, 1.5e-3);
}

#[test]
fn mode_area_of_a_box_is_a_quarter_of_it() {
    let r = Request2d {
        width: 1.0,
        height: 1.5,
        cell: 0.05,
        background: Material { eps: 2.0, mu: 1.0 },
        shapes: vec![],
        polarization: Polarization::Te,
        walls: Boundary::Pec,
        analysis: Analysis2d::Cavity { n_modes: 1, bloch: [0.0, 0.0], target_k0: None },
        monitors: vec![],
        density: None,
    };
    // sin^2 sin^2 averages to 1/4 and peaks at 1, whatever eps is.
    let o = eigenmenno_core::solve2d::solve(&r).unwrap();
    close(o.meta.fields[0].volume.unwrap(), 1.5 / 4.0, 1e-4);
}

fn ldos(r: &Request2d) -> f64 {
    eigenmenno_core::solve2d::solve(r).unwrap().meta.fields[0].ldos.unwrap()
}

fn ldos_request(shapes: Vec<Shape>, k: f64, average_q: Option<f64>) -> Request2d {
    Request2d {
        width: 2.0,
        height: 2.0,
        cell: 0.1,
        background: Material { eps: 1.0, mu: 1.0 },
        shapes,
        polarization: Polarization::Te,
        walls: Boundary::Absorbing,
        analysis: Analysis2d::Driven { wavelength: 2.0 * PI / k, source: Source2d::Line { x: 1.0, y: 1.0 }, pml: 0.5, average_q },
        monitors: vec![],
        density: None,
    }
}

#[test]
fn vacuum_ldos_is_half_of_k() {
    let r = Request2d { cell: 1.0 / 24.0, ..ldos_request(vec![], 2.0 * PI, None) };
    close(ldos(&r), PI, 2e-2 * PI);
}

#[test]
fn complex_frequency_is_a_lorentzian_average_of_real_ones() {
    let rod = vec![Shape { outline: Outline::Circle { x: 1.0, y: 1.25, r: 0.25 }, material: Material { eps: 12.0, mu: 1.0 } }];
    let (k, q) = (4.0, 20.0);
    let gamma = k / (2.0 * q);
    let averaged = 2.0 * ldos(&ldos_request(rod.clone(), k, Some(q))) - ldos(&ldos_request(rod.clone(), k, Some(q / 2.0)));
    // One Lorentzian falls off as 1 / k^2. That weights the sharp modes far above the mesh's
    // cutoff, which no quadrature can resolve. Two Lorentzians, of half-widths gamma and
    // 2 gamma, combine to fall off as 1 / k^4. With k' = k + gamma tan(t), their kernel has
    // the weight below in t, on the whole real axis. Negative k' has the same F.
    let n = 300;
    let dt = PI / n as f64;
    let direct: f64 = (0..n)
        .map(|i| {
            let t = (i as f64 + 0.5) * dt - PI / 2.0;
            let weight = 6.0 / PI / (t.tan().powi(2) + 4.0);
            weight * ldos(&ldos_request(rod.clone(), (k + gamma * t.tan()).abs(), None)) * dt
        })
        .sum();
    eprintln!("LDOS averaged {averaged:.6}, by quadrature {direct:.6}");
    close(averaged, direct, 1e-3 * direct.abs());
}

#[test]
fn leaky_mode_matches_the_ldos_peak_and_q_grows_with_the_crystal() {
    let mode = |n: usize| {
        let c = n as f64 / 2.0;
        let analysis = Analysis2d::Resonance { target_k0: 2.0 * PI * 0.38, x: c, y: c, pml: 1.3 };
        let o = eigenmenno_core::solve2d::solve(&defect_crystal(n, 1.0 / 12.0, Boundary::Absorbing, analysis)).unwrap();
        let f = &o.meta.fields[0];
        eprintln!("{n} x {n} crystal: a/lambda {:.5}, Q {:.1}, V {:.3}", f.k0 / (2.0 * PI), f.q.unwrap(), f.volume.unwrap());
        (f.k0, f.q.unwrap(), f.volume.unwrap())
    };
    let (k3, q3, _) = mode(3);
    let (k5, q5, v5) = mode(5);
    let (_, q7, v7) = mode(7);
    assert!(q5 > 5.0 * q3 && q7 > 5.0 * q5, "Q {q3} {q5} {q7}");
    // Past a few rods, more crystal only stops leakage, and the mode's shape settles.
    close(v7, v5, 0.03 * v5);
    assert!(v7 < 2.0, "area {v7}");
    // A Lorentzian at real frequencies falls to half height at k (1 +- 1 / 2Q).
    for (n, k, q) in [(3, k3, q3), (5, k5, q5)] {
        let driven = |k: f64| {
            let source = Source2d::Line { x: n as f64 / 2.0, y: n as f64 / 2.0 };
            ldos(&defect_crystal(n, 1.0 / 12.0, Boundary::Absorbing, Analysis2d::Driven { wavelength: 2.0 * PI / k, source, pml: 1.3, average_q: None }))
        };
        let half = k / (2.0 * q);
        let (peak, low, high) = (driven(k), driven(k - half), driven(k + half));
        eprintln!("{n} x {n} LDOS: {low:.4} {peak:.4} {high:.4}");
        assert!(peak > low && peak > high);
        close(0.5 * (low + high) / peak, 0.5, 0.05);
    }
}

#[test]
fn an_empty_domain_has_no_leaky_mode() {
    let r = Request2d {
        analysis: Analysis2d::Resonance { target_k0: 4.0, x: 1.0, y: 1.0, pml: 0.5 },
        ..ldos_request(vec![], 4.0, None)
    };
    assert!(eigenmenno_core::solve2d::solve(&r).is_err());
}

fn monitors(o: &Output) -> &[f64] {
    channel(o, "monitors")
}

#[test]
fn power_around_a_line_source_matches_its_ldos() {
    // What the source supplies, `-Re(u_0) / 2 = F / 4`, must cross any loop around it.
    let (w, m) = (3.0, 0.33);
    let r = Request2d {
        width: w,
        height: w,
        cell: 1.0 / 24.0,
        background: Material { eps: 1.0, mu: 1.0 },
        shapes: vec![Shape { outline: Outline::Circle { x: 1.3, y: 1.4, r: 0.3 }, material: Material { eps: 6.0, mu: 1.0 } }],
        polarization: Polarization::Te,
        walls: Boundary::Absorbing,
        analysis: Analysis2d::Driven { wavelength: 1.0, source: Source2d::Line { x: 1.8, y: 1.5 }, pml: 0.6, average_q: None },
        // Each side runs with the outside on its right.
        monitors: vec![
            Monitor { x0: w - m, y0: m, x1: w - m, y1: w - m },
            Monitor { x0: w - m, y0: w - m, x1: m, y1: w - m },
            Monitor { x0: m, y0: w - m, x1: m, y1: m },
            Monitor { x0: m, y0: m, x1: w - m, y1: m },
        ],
        density: None,
    };
    let o = eigenmenno_core::solve2d::solve(&r).unwrap();
    let total: f64 = monitors(&o).iter().sum();
    close(total, o.meta.fields[0].ldos.unwrap() / 4.0, 3e-3 * total);
}

/// A slab of index 3.5 and width `w` along x through y = 1.5, or bent sharply up at x = 2.5,
/// in a 4 x 4 domain, launched by a port at x = 0.5. Wavelength 1.
fn slab(w: f64, bent: bool, cell: f64, span: f64) -> Request2d {
    let (cx, cy, far) = (2.5, 1.5, 6.0);
    let core = Material { eps: 3.5 * 3.5, mu: 1.0 };
    let along = |x1: f64| Shape { outline: Outline::Rect { x0: -far, y0: cy - w / 2.0, x1, y1: cy + w / 2.0 }, material: core };
    let (shapes, monitor) = if bent {
        let up = Shape { outline: Outline::Rect { x0: cx - w / 2.0, y0: cy - w / 2.0, x1: cx + w / 2.0, y1: far }, material: core };
        (vec![along(cx + w / 2.0), up], Monitor { x0: cx + span, y0: 3.5, x1: cx - span, y1: 3.5 })
    } else {
        (vec![along(far)], Monitor { x0: 3.5, y0: cy - span, x1: 3.5, y1: cy + span })
    };
    Request2d {
        width: 4.0,
        height: 4.0,
        cell,
        background: Material { eps: 1.0, mu: 1.0 },
        shapes,
        polarization: Polarization::Te,
        walls: Boundary::Absorbing,
        analysis: Analysis2d::Driven { wavelength: 1.0, source: Source2d::Port { x: 0.5, y0: cy - 1.2, y1: cy + 1.2 }, pml: 1.0, average_q: None },
        monitors: vec![monitor],
        density: None,
    }
}

#[test]
fn port_launches_the_exact_slab_mode() {
    let (w, k) = (0.12, 2.0 * PI);
    // An even mode has kappa tan(kappa w / 2) = gamma, with kappa^2 + gamma^2 = k^2 (n^2 - 1).
    let v2 = k * k * (3.5f64 * 3.5 - 1.0);
    let (mut lo, mut hi) = (1e-6, (PI / w).min(v2.sqrt()) - 1e-9);
    for _ in 0..80 {
        let kappa = (lo + hi) / 2.0;
        if kappa * (kappa * w / 2.0).tan() > (v2 - kappa * kappa).sqrt() {
            hi = kappa;
        } else {
            lo = kappa;
        }
    }
    let kappa = lo;
    let gamma = (v2 - kappa * kappa).sqrt();
    let beta = (k * k * 3.5 * 3.5 - kappa * kappa).sqrt();
    let o = eigenmenno_core::solve2d::solve(&slab(w, false, 1.0 / 32.0, 1.2)).unwrap();
    let on_axis = |x: f64| field_of(&o, "Ez").into_iter().min_by(|a, b| (a.0 - x).hypot(a.1 - 1.5).total_cmp(&(b.0 - x).hypot(b.1 - 1.5))).unwrap();
    // The sheet's field is `-(k / 2 beta) phi exp(i beta |x - x0|)`, on both sides. Four
    // cells across the core leave a small error in beta, which builds up with distance.
    for x in [0.0, 1.5, 3.0] {
        let (px, _, u) = on_axis(x);
        let want = C64::from_polar(-k / (2.0 * beta), beta * (px - 0.5).abs());
        let error = (u - want).norm() / want.norm();
        eprintln!("port field at {px}: relative error {error:.2e}");
        assert!(error < 2e-2, "at {px}: {u} against {want}");
    }
    // So its power is `(k / 8 beta) integral(phi^2)`, with phi peaking at 1.
    let c = (kappa * w / 2.0).cos();
    let integral = w / 2.0 + (kappa * w).sin() / (2.0 * kappa) + c * c / gamma;
    close(monitors(&o)[0], k / (8.0 * beta) * integral, 5e-3 * k / (8.0 * beta) * integral);
}

#[test]
fn ports_are_checked() {
    let source = Source2d::Port { x: 0.5, y0: 0.0, y1: 9.0 };
    let mut r = slab(0.12, false, 0.2, 1.0);
    r.analysis = Analysis2d::Driven { wavelength: 1.0, source, pml: 1.0, average_q: None };
    assert!(eigenmenno_core::solve2d::solve(&r).is_err());
    let mut empty = slab(0.12, false, 0.2, 1.0);
    empty.shapes.clear();
    assert!(eigenmenno_core::solve2d::solve(&empty).is_err());
}

/// Mekis et al. 1996: rods of index 3.4 and radius 0.18a on a square lattice, with rows
/// removed. The bends are sharp (`L = 0.33 sqrt(2) a`) or take one step (`1.33 sqrt(2) a`).
#[derive(Clone, Copy, PartialEq)]
enum Guide {
    Straight,
    Sharp,
    Stepped,
}

fn lattice_guide(kind: Guide, freq: f64) -> Request2d {
    let (arm, clad, pml) = (5i32, 3i32, 6.0);
    let (size, i0, j0) = (arm + 1 + clad, arm, clad);
    let reach = pml as i32 + 1;
    let rods = (-reach..size + reach).flat_map(|i| (-reach..size + reach).map(move |j| (i, j)));
    let shapes = rods
        .filter(|&(i, j)| match kind {
            Guide::Straight => j != j0,
            Guide::Sharp => !((j == j0 && i <= i0) || (i == i0 && j >= j0)),
            Guide::Stepped => !((j == j0 && i < i0) || (i == i0 - 1 && j == j0 + 1) || (i == i0 && j > j0)),
        })
        .map(|(i, j)| Shape { outline: Outline::Circle { x: i as f64 + 0.5, y: j as f64 + 0.5, r: 0.18 }, material: Material { eps: 3.4 * 3.4, mu: 1.0 } })
        .collect();
    let (cx, cy, span, far) = (i0 as f64 + 0.5, j0 as f64 + 0.5, 3.0, size as f64 - 1.03);
    let output = if kind == Guide::Straight {
        Monitor { x0: far, y0: cy - span, x1: far, y1: cy + span }
    } else {
        Monitor { x0: cx + span, y0: far, x1: cx - span, y1: far }
    };
    Request2d {
        width: size as f64,
        height: size as f64,
        cell: 1.0 / 8.0,
        background: Material { eps: 1.0, mu: 1.0 },
        shapes,
        polarization: Polarization::Te,
        walls: Boundary::Absorbing,
        analysis: Analysis2d::Driven { wavelength: 1.0 / freq, source: Source2d::Line { x: 1.5, y: cy }, pml, average_q: None },
        monitors: vec![output, Monitor { x0: 3.03, y0: cy - span, x1: 3.03, y1: cy + span }],
        density: None,
    }
}

#[test]
fn guided_band_matches_mekis() {
    // A guide with one rod per period, between lattices, with k along it.
    let band = |phase: f64, want: f64| {
        let shapes = (0..7).filter(|&j| j != 3).map(|j| Shape { outline: Outline::Circle { x: 0.5, y: j as f64 + 0.5, r: 0.18 }, material: Material { eps: 3.4 * 3.4, mu: 1.0 } }).collect();
        let r = Request2d {
            width: 1.0,
            height: 7.0,
            cell: 1.0 / 24.0,
            background: Material { eps: 1.0, mu: 1.0 },
            shapes,
            polarization: Polarization::Te,
            walls: Boundary::Periodic,
            analysis: Analysis2d::Cavity { n_modes: 3, bloch: [2.0 * PI * phase, 0.0], target_k0: Some(2.0 * PI * want) },
            monitors: vec![],
            density: None,
        };
        let o = eigenmenno_core::solve2d::solve(&r).unwrap();
        let nearest = o.meta.fields.iter().map(|f| f.k0 / (2.0 * PI)).min_by(|a, b| (a - want).abs().total_cmp(&(b - want).abs())).unwrap();
        close(nearest, want, 2e-3);
    };
    // "The band appears at 0.312 when k = 0 and reaches the top of the gap when k = 0.38",
    // and that top is 0.443.
    band(0.0, 0.312);
    band(0.38, 0.443);
}

#[test]
fn sharp_bends_transmit_as_in_mekis() {
    let solve_guide = |kind, freq| eigenmenno_core::solve2d::solve(&lattice_guide(kind, freq)).unwrap();
    let transmission = |kind, freq| {
        let (straight, bend) = (solve_guide(Guide::Straight, freq), solve_guide(kind, freq));
        // Nothing radiates in a band gap, so all power passing the input reaches the output.
        close(monitors(&bend)[0], monitors(&bend)[1], 2e-3 * monitors(&bend)[0]);
        let t = monitors(&bend)[0] / monitors(&straight)[0];
        eprintln!("frequency {freq}: T {t:.4}");
        t
    };
    // The abstract's 98% maximum for the sharp bend, and over 95% below 0.403.
    close(transmission(Guide::Sharp, 0.353), 0.98, 0.01);
    assert!(transmission(Guide::Sharp, 0.392) > 0.95);
    assert!(transmission(Guide::Sharp, 0.403) > 0.95);
    // Figure 4 has its reflection rising to about 0.1 by 0.42.
    close(transmission(Guide::Sharp, 0.42), 0.9, 0.03);
    // The stepped bend was tuned to have no reflection at 0.353.
    close(transmission(Guide::Stepped, 0.353), 1.0, 0.005);
}

#[test]
fn a_sharp_bend_in_a_slab_loses_most_of_its_power() {
    // Mekis et al. only say that this reaches at most 30%, with index 3.5. Their width isn't
    // given for a sharp bend, so this uses the width over wavelength of their other bend.
    let (w, cell, span) = (0.143, 1.0 / 32.0, 1.0);
    let (straight, bend) = (slab(w, false, cell, span), slab(w, true, cell, span));
    let t = monitors(&eigenmenno_core::solve2d::solve(&bend).unwrap())[0] / monitors(&eigenmenno_core::solve2d::solve(&straight).unwrap())[0];
    eprintln!("slab bend T {t:.4}");
    assert!(t > 0.1 && t < 0.3, "T {t}");
}

/// A slab bent around a square of pixels. Everything sits on the grid of cells 1/24.
fn density_bend(rho: Vec<f64>) -> Request2d {
    let (w, n, core) = (0.125, 12, Material { eps: 3.5 * 3.5, mu: 1.0 });
    let (x0, y0, size) = (1.0, 0.5, 0.5);
    let (cx, cy) = (x0 + size / 2.0, y0 + size / 2.0);
    let rect = |x0, y0, x1, y1| Shape { outline: Outline::Rect { x0, y0, x1, y1 }, material: core };
    Request2d {
        width: 2.5,
        height: 2.5,
        cell: 1.0 / 24.0,
        background: Material { eps: 1.0, mu: 1.0 },
        shapes: vec![rect(-3.0, cy - w / 2.0, x0, cy + w / 2.0), rect(cx - w / 2.0, y0 + size, cx + w / 2.0, 5.0)],
        polarization: Polarization::Te,
        walls: Boundary::Absorbing,
        analysis: Analysis2d::Driven { wavelength: 1.0, source: Source2d::Port { x: 0.4, y0: 0.0, y1: 1.5 }, pml: 0.5, average_q: None },
        monitors: vec![Monitor { x0: cx + 0.6, y0: 2.1, x1: cx - 0.6, y1: 2.1 }],
        density: Some(Density { x0, y0, width: size, height: size, nx: n, ny: n, eps_low: 1.0, eps_high: 3.5 * 3.5, rho }),
    }
}

#[test]
fn density_gradient_matches_finite_differences() {
    let n = 12 * 12;
    let rho: Vec<f64> = (0..n).map(|i| 0.5 + 0.4 * ((i * 37 % 101) as f64 / 101.0 - 0.5)).collect();
    let power = |rho: &[f64]| {
        let o = eigenmenno_core::solve2d::solve(&density_bend(rho.to_vec())).unwrap();
        (monitors(&o)[0], channel(&o, "gradient").to_vec())
    };
    let (_, gradient) = power(&rho);
    for pixel in [0, 40, 77, 143] {
        let step = 1e-5;
        let (mut up, mut down) = (rho.clone(), rho.clone());
        up[pixel] += step;
        down[pixel] -= step;
        let difference = (power(&up).0 - power(&down).0) / (2.0 * step);
        eprintln!("pixel {pixel}: adjoint {:.6e}, difference {difference:.6e}", gradient[pixel]);
        close(gradient[pixel], difference, 1e-5 * difference.abs().max(1e-3));
    }
}

#[test]
fn density_requests_are_checked() {
    let mut wrong = density_bend(vec![0.5; 144]);
    wrong.polarization = Polarization::Tm;
    assert!(eigenmenno_core::solve2d::solve(&wrong).is_err());
    assert!(eigenmenno_core::solve2d::solve(&density_bend(vec![0.5; 3])).is_err());
    let mut unmeasured = density_bend(vec![0.5; 144]);
    unmeasured.monitors.clear();
    assert!(eigenmenno_core::solve2d::solve(&unmeasured).is_err());
}

/// The page's optimizer's result for a bend over wavelengths 0.95 to 1.05, top row first.
const OPTIMIZED_BEND: [&str; 24] = [
    "........######....##.###",
    "........#####......#####",
    "........#####......#####",
    ".......######......#####",
    "#......#####......######",
    "##.....#####......######",
    "......#####.....#..#####",
    "...#######....####..####",
    "..#######.....####...###",
    "########.....####....###",
    "#######.....#####.....##",
    "######......#####.....##",
    "####.......#####......##",
    "#........#######.....###",
    "........######.......###",
    "......#######.........##",
    ".....#######...........#",
    "#########..............#",
    "#######................#",
    "..###...............####",
    "............#......#####",
    "........#######...#####.",
    "...##############.####..",
    "...##################...",
];

/// A slab of width 0.125 turning a corner through that pattern, as the page builds it, or an
/// empty region, or a straight guide. `extra` lengthens the output.
fn optimized_bend(pattern: Option<&[&str; 24]>, straight: bool, wavelength: f64, cell: f64, extra: f64, monitors: Vec<Monitor>) -> Request2d {
    let core = Material { eps: 3.5 * 3.5, mu: 1.0 };
    let (x0, y0, size, w) = (1.0, 0.5, 1.0, 0.125);
    let (cx, cy) = (x0 + size / 2.0, y0 + size / 2.0);
    let rect = |x0, y0, x1, y1| Shape { outline: Outline::Rect { x0, y0, x1, y1 }, material: core };
    let shapes = if straight {
        vec![rect(-3.0, cy - w / 2.0, 9.0, cy + w / 2.0)]
    } else {
        vec![rect(-3.0, cy - w / 2.0, x0, cy + w / 2.0), rect(cx - w / 2.0, y0 + size, cx + w / 2.0, 9.0)]
    };
    let rho = |p: &[&str; 24]| p.iter().rev().flat_map(|row| row.chars().map(|c| if c == '#' { 1.0 } else { 0.0 })).collect::<Vec<f64>>();
    let density = (!straight).then(|| Density {
        x0, y0, width: size, height: size, nx: 24, ny: 24, eps_low: 1.0, eps_high: 3.5 * 3.5,
        rho: pattern.map_or(vec![0.0; 576], rho),
    });
    Request2d {
        width: 2.5,
        height: 2.5 + extra,
        cell,
        background: Material { eps: 1.0, mu: 1.0 },
        shapes,
        polarization: Polarization::Te,
        walls: Boundary::Absorbing,
        analysis: Analysis2d::Driven { wavelength, source: Source2d::Port { x: 0.4, y0: cy - 0.7, y1: cy + 0.7 }, pml: 0.5, average_q: None },
        monitors,
        density,
    }
}

#[test]
fn optimized_bend_survives_a_finer_mesh() {
    let (cx, cy) = (1.5, 1.0);
    let power = |r: Request2d| monitors(&eigenmenno_core::solve2d::solve(&r).unwrap())[0];
    for wavelength in [0.95, 1.0, 1.05] {
        for cell in [1.0 / 24.0, 1.0 / 48.0] {
            // Far along the guide, and narrow, so that no light radiated at the corner counts.
            let (far, span) = (2.85, 0.35);
            let straight = power(optimized_bend(None, true, wavelength, cell, 0.75, vec![Monitor { x0: 2.1, y0: cy - span, x1: 2.1, y1: cy + span }]));
            let output = vec![Monitor { x0: cx + span, y0: far, x1: cx - span, y1: far }];
            let t = power(optimized_bend(Some(&OPTIMIZED_BEND), false, wavelength, cell, 0.75, output.clone())) / straight;
            let plain = power(optimized_bend(None, false, wavelength, cell, 0.75, output)) / straight;
            eprintln!("wavelength {wavelength} cell {cell:.4}: T {t:.4}, empty region {plain:.4}");
            assert!(t > 0.97 && t < 1.005, "T {t}");
            assert!(plain < 0.15, "empty region {plain}");
        }
    }
}

#[test]
fn power_is_conserved_around_an_optimized_bend() {
    // The port is outside this loop, and everything inside it is lossless.
    let (a, b, c, d) = (0.7, 2.3, 0.2, 2.3);
    let around = vec![
        Monitor { x0: b, y0: c, x1: b, y1: d },
        Monitor { x0: b, y0: d, x1: a, y1: d },
        Monitor { x0: a, y0: d, x1: a, y1: c },
        Monitor { x0: a, y0: c, x1: b, y1: c },
    ];
    let o = eigenmenno_core::solve2d::solve(&optimized_bend(Some(&OPTIMIZED_BEND), false, 1.0, 1.0 / 48.0, 0.0, around)).unwrap();
    let m = monitors(&o);
    let net: f64 = m.iter().sum();
    let largest = m.iter().fold(0f64, |x, v| x.max(v.abs()));
    eprintln!("around the bend: {m:?}, net {net:.2e}");
    assert!(net.abs() < 2e-3 * largest, "net {net} of {largest}");
}
