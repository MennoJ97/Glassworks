//! The 1D solver against exact results: every test compares with a formula, not stored numbers.
//! Chapter 8 explains why. `analytic2d.rs` does the same for 2D, and `gratings.rs` for gratings.

mod common;

use common::{cavity, channel, close, layer, request};
use glassworks_core::output::Output;
use glassworks_core::problem::{Analysis, Boundary, Layer, Polarization, Request, Source};
use glassworks_core::solve1d::solve;
use glassworks_core::C64;
use std::f64::consts::PI;

#[test]
fn pec_cavity_matches_mode_formula() {
    let r =
        request(vec![layer(2.0, 1.0)], Polarization::Te, (Boundary::Pec, Boundary::Pec), 4, 0.05, cavity(6));
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
            let r = request(
                vec![layer(1.0, 1.0)],
                Polarization::Te,
                (Boundary::Pec, Boundary::Pec),
                order,
                h,
                cavity(1),
            );
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
    let r =
        request(vec![layer(1.0, 1.0)], Polarization::Tm, (Boundary::Pec, Boundary::Pec), 3, 0.05, cavity(3));
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
    let r = request(
        vec![layer(d1, n1), layer(d2, n2)],
        Polarization::Te,
        (Boundary::Pec, Boundary::Pec),
        5,
        0.05,
        cavity(5),
    );
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
    let r = request(
        vec![layer(1.0, 1.0)],
        Polarization::Te,
        (Boundary::Periodic, Boundary::Periodic),
        4,
        0.05,
        cavity(5),
    );
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
    let n2: Vec<C64> =
        layers.iter().map(|l| C64::new(l.eps[0], l.eps[1]) * C64::new(l.mu[0], l.mu[1])).collect();
    let kt = k0 * n2[0].sqrt().re * angle.sin();
    let kx = |i: usize| {
        let w = (n2[i] * k0 * k0 - kt * kt).sqrt();
        if w.im < 0.0 {
            -w
        } else {
            w
        }
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
    let layers =
        vec![layer(0.3, 1.0), layer(0.17, 2.3), layer(0.27, 1.38), layer(0.17, 2.3), layer(0.3, 1.5)];
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

fn gradient(
    layers: &[Layer],
    pol: Polarization,
    right: Boundary,
    order: usize,
    h: f64,
    angle: f64,
) -> Output {
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
    let r =
        request(vec![layer(1.0, 1.0)], Polarization::Te, (Boundary::Pec, Boundary::Pec), 5, 0.05, cavity(1));
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
    let bad = request(
        vec![layer(1.0, 1.0)],
        Polarization::Te,
        (Boundary::Periodic, Boundary::Pec),
        1,
        0.1,
        cavity(1),
    );
    assert!(solve(&bad).is_err());
    let open = request(
        vec![layer(1.0, 1.0)],
        Polarization::Te,
        (Boundary::Port, Boundary::Absorbing),
        1,
        0.1,
        cavity(1),
    );
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
