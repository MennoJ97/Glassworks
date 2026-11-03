//! Gratings against Fresnel, an independent Fourier modal method, symmetry and finite differences.
//! Chapter 15 explains them.

mod common;

use common::{channel, close};
use glassworks_core::mesh::Outline;
use glassworks_core::output::Output;
use glassworks_core::problem::{
    Analysis2d, Boundary, Density, Material, OrderGoal, Polarization, Request2d, Shape,
};

/// A ridge on a substrate, one period of 1.5, between layers of `pml`. Lines sit on grid rows.
fn grating_request(
    n_ridge: f64,
    x1: f64,
    angle_deg: f64,
    pol: Polarization,
    cell: f64,
    pml: f64,
) -> Request2d {
    let glass = Material { eps: 2.25, mu: 1.0 };
    let rect = |x0, y0, x1, y1, material| Shape { outline: Outline::Rect { x0, y0, x1, y1 }, material };
    Request2d {
        width: 1.5,
        height: 2.0,
        cell,
        background: Material { eps: 1.0, mu: 1.0 },
        shapes: vec![
            rect(-1.0, -3.0, 2.5, 0.5, glass),
            rect(0.375, 0.5, x1, 1.0, Material { eps: n_ridge * n_ridge, mu: 1.0 }),
        ],
        polarization: pol,
        walls: Boundary::Periodic,
        analysis: Analysis2d::Grating {
            wavelength: 1.0,
            angle_deg,
            pml,
            source_y: 1.75,
            reflect_y: 1.25,
            transmit_y: 0.25,
            goal: None,
        },
        monitors: vec![],
        density: None,
    }
}

fn grating_orders(o: &Output, name: &str) -> Vec<(i32, f64)> {
    let channel = |n: &str| &o.channels.iter().find(|c| c.0 == n).unwrap().1;
    channel(&format!("{name}_orders")).iter().map(|m| *m as i32).zip(channel(name).iter().copied()).collect()
}

fn share(orders: &[(i32, f64)], m: i32) -> f64 {
    orders.iter().find(|e| e.0 == m).unwrap_or_else(|| panic!("no order {m}")).1
}

#[test]
fn a_flat_interface_at_an_angle_matches_fresnel() {
    // A grating with no ridge has one order each way, whatever its period.
    for (pol, angle) in [(Polarization::Te, 0.0), (Polarization::Te, 35.0), (Polarization::Tm, 35.0)] {
        let mut r = grating_request(1.0, 0.375, angle, pol, 1.0 / 40.0, 1.0);
        r.shapes.truncate(1);
        r.width = 0.5;
        let o = glassworks_core::solve2d::solve(&r).unwrap();
        let (cos_i, sin_i) = ((angle as f64).to_radians().cos(), (angle as f64).to_radians().sin());
        let cos_t = (1.0 - (sin_i / 1.5f64).powi(2)).sqrt();
        let amplitude = if pol == Polarization::Te {
            (cos_i - 1.5 * cos_t) / (cos_i + 1.5 * cos_t)
        } else {
            (1.5 * cos_i - cos_t) / (1.5 * cos_i + cos_t)
        };
        close(share(&grating_orders(&o, "reflected"), 0), amplitude * amplitude, 1e-5);
        close(share(&grating_orders(&o, "transmitted"), 0), 1.0 - amplitude * amplitude, 1e-5);
    }
}

#[test]
fn a_grating_matches_the_fourier_modal_method() {
    // From a separate Fourier modal method, which has no mesh or layer, and is steady to 1e-6
    // at 320 orders.
    let expected: [(Polarization, [f64; 3], [f64; 4]); 2] = [
        (Polarization::Te, [0.002538, 0.053613, 0.138772], [0.284394, 0.218096, 0.118276, 0.184311]),
        (Polarization::Tm, [0.010134, 0.047082, 0.123978], [0.055278, 0.168438, 0.265488, 0.329601]),
    ];
    for (pol, reflected, transmitted) in expected {
        let o = glassworks_core::solve2d::solve(&grating_request(2.5, 0.975, 20.0, pol, 1.0 / 80.0, 1.0))
            .unwrap();
        let (r, t) = (grating_orders(&o, "reflected"), grating_orders(&o, "transmitted"));
        assert_eq!(r.iter().map(|e| e.0).collect::<Vec<_>>(), [-2, -1, 0]);
        assert_eq!(t.iter().map(|e| e.0).collect::<Vec<_>>(), [-2, -1, 0, 1]);
        reflected.iter().zip(&r).for_each(|(want, e)| close(e.1, *want, 1e-4));
        transmitted.iter().zip(&t).for_each(|(want, e)| close(e.1, *want, 1e-4));
        close(r.iter().chain(&t).map(|e| e.1).sum(), 1.0, 1e-4);
    }
}

#[test]
fn a_symmetric_grating_mirrors_its_orders() {
    let solve = |angle| {
        glassworks_core::solve2d::solve(&grating_request(
            1.5,
            1.125,
            angle,
            Polarization::Te,
            1.0 / 40.0,
            1.0,
        ))
        .unwrap()
    };
    let (right, left) = (solve(15.0), solve(-15.0));
    for name in ["reflected", "transmitted"] {
        let (a, b) = (grating_orders(&right, name), grating_orders(&left, name));
        assert_eq!(a.len(), b.len());
        a.iter().for_each(|e| close(e.1, share(&b, -e.0), 1e-4));
    }
}

#[test]
fn scaling_a_grating_with_its_wavelength_changes_nothing() {
    let shares = |scale: f64| {
        let mut r = grating_request(2.5, 0.975, 20.0, Polarization::Te, scale / 40.0, scale);
        r.width *= scale;
        r.height *= scale;
        for shape in &mut r.shapes {
            let Outline::Rect { x0, y0, x1, y1 } = shape.outline else { unreachable!("all rectangles") };
            shape.outline = Outline::Rect { x0: x0 * scale, y0: y0 * scale, x1: x1 * scale, y1: y1 * scale };
        }
        r.analysis = Analysis2d::Grating {
            wavelength: scale,
            angle_deg: 20.0,
            pml: scale,
            source_y: 1.75 * scale,
            reflect_y: 1.25 * scale,
            transmit_y: 0.25 * scale,
            goal: None,
        };
        channel(&glassworks_core::solve2d::solve(&r).unwrap(), "transmitted").to_vec()
    };
    shares(1.0).iter().zip(shares(2.833)).for_each(|(a, b)| close(*a, b, 1e-9));
}

#[test]
fn a_grating_gradient_matches_finite_differences() {
    // At an angle, so the ties really do break A^T = A, and for orders on both sides.
    for (transmitted, order, pol) in [
        (true, 1, Polarization::Te),
        (false, -1, Polarization::Te),
        (true, 0, Polarization::Te),
        (true, 1, Polarization::Tm),
        (false, -1, Polarization::Tm),
    ] {
        let make = |rho: &[f64]| {
            let mut r = grating_request(1.0, 0.375, 20.0, pol, 1.0 / 40.0, 1.0);
            r.shapes.truncate(1);
            r.density = Some(Density {
                x0: 0.0,
                y0: 0.5,
                width: 1.5,
                height: 0.5,
                nx: 6,
                ny: 2,
                eps_low: 1.0,
                eps_high: 6.25,
                rho: rho.to_vec(),
            });
            if let Analysis2d::Grating { goal, .. } = &mut r.analysis {
                *goal = Some(OrderGoal { transmitted, order });
            }
            r
        };
        let value = |rho: &[f64]| {
            let o = glassworks_core::solve2d::solve(&make(rho)).unwrap();
            share(&grating_orders(&o, if transmitted { "transmitted" } else { "reflected" }), order)
        };
        let rho = [0.2, 0.7, 0.9, 0.4, 0.1, 0.6, 0.8, 0.3, 0.5, 0.15, 0.65, 0.35];
        let o = glassworks_core::solve2d::solve(&make(&rho)).unwrap();
        let gradient = channel(&o, "gradient");
        for pixel in [0, 3, 7, 10] {
            let (mut up, mut down) = (rho, rho);
            up[pixel] += 1e-5;
            down[pixel] -= 1e-5;
            let slope = (value(&up) - value(&down)) / 2e-5;
            assert!(
                (gradient[pixel] - slope).abs() < 1e-5 + 1e-4 * slope.abs(),
                "pixel {pixel}: {} and {slope}",
                gradient[pixel]
            );
        }
    }
}

#[test]
fn a_grating_refuses_lines_it_cannot_read() {
    let lines = |source_y: f64, reflect_y: f64, transmit_y: f64| {
        let mut r = grating_request(1.5, 1.125, 0.0, Polarization::Te, 1.0 / 20.0, 1.0);
        r.analysis = Analysis2d::Grating {
            wavelength: 1.0,
            angle_deg: 0.0,
            pml: 1.0,
            source_y,
            reflect_y,
            transmit_y,
            goal: None,
        };
        glassworks_core::solve2d::solve(&r).map(|_| ()).unwrap_err()
    };
    assert!(lines(1.0, 1.25, 0.25).contains("in order"));
    assert!(lines(1.26, 1.25, 0.25).contains("wavelengths above"));
    // Across the ridge, where no single order is defined.
    assert!(lines(1.75, 0.75, 0.25).contains("uniform"));
    // A patch across only an upper reflect line, or a layer filling one, changes that line alone.
    let crossing = |x0: f64, x1: f64| {
        let mut r = grating_request(1.5, 1.125, 0.0, Polarization::Te, 1.0 / 20.0, 1.0);
        let outline = Outline::Rect { x0, y0: 1.5, x1, y1: 1.7 };
        r.shapes.push(Shape { outline, material: Material { eps: 4.0, mu: 1.0 } });
        glassworks_core::solve2d::solve(&r).map(|_| ()).unwrap_err()
    };
    assert!(crossing(0.1, 0.4).contains("uniform"));
    assert!(crossing(-1.0, 2.5).contains("one material"));
}
