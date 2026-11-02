//! Designing in `H_z`: its gradients against finite differences, and an in-plane electric dipole
//! against vacuum, energy balance and a drawn structure. Chapter 19 explains the dipole.

mod common;

use common::{channel, close};
use eigenmenno_core::mesh::Outline;
use eigenmenno_core::output::Output;
use eigenmenno_core::problem::{
    Analysis2d, Boundary, Density, Material, Monitor, Polarization, Request2d, Shape, Source2d,
};
use eigenmenno_core::solve2d::solve;

const CORE: Material = Material { eps: 12.25, mu: 1.0 };

fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Shape {
    Shape { outline: Outline::Rect { x0, y0, x1, y1 }, material: CORE }
}

fn pixels(rho: &[f64]) -> Density {
    Density {
        x0: 1.0,
        y0: 0.5,
        width: 0.5,
        height: 0.5,
        nx: 12,
        ny: 12,
        eps_low: 1.0,
        eps_high: 12.25,
        rho: rho.to_vec(),
    }
}

fn open(source: Source2d, q: Option<f64>) -> Request2d {
    Request2d {
        width: 2.5,
        height: 2.5,
        cell: 1.0 / 24.0,
        background: Material { eps: 1.0, mu: 1.0 },
        shapes: vec![],
        polarization: Polarization::Tm,
        walls: Boundary::Absorbing,
        analysis: Analysis2d::Driven { wavelength: 1.0, source, pml: 0.5, average_q: q },
        monitors: vec![],
        density: None,
    }
}

fn start() -> Vec<f64> {
    (0..144).map(|i| 0.5 + 0.4 * ((i * 37 % 101) as f64 / 101.0 - 0.5)).collect()
}

fn run(r: &Request2d) -> Output {
    solve(r).unwrap_or_else(|e| panic!("{e}"))
}

/// Checks a gradient channel against central differences of `value`.
fn check(make: impl Fn(&[f64]) -> Request2d, value: impl Fn(&Output) -> f64, name: &str, chosen: &[usize]) {
    let rho = start();
    let exact = channel(&run(&make(&rho)), name).to_vec();
    for &pixel in chosen {
        let (mut up, mut down, h) = (rho.clone(), rho.clone(), 1e-5);
        (up[pixel], down[pixel]) = (rho[pixel] + h, rho[pixel] - h);
        let numeric = (value(&run(&make(&up))) - value(&run(&make(&down)))) / (2.0 * h);
        close(exact[pixel], numeric, 2e-5 * numeric.abs().max(1e-3));
    }
}

#[test]
fn a_monitors_gradient_matches_finite_differences() {
    let make = |rho: &[f64]| {
        let mut r = open(Source2d::Port { x: 0.4, y0: 0.0, y1: 1.5 }, None);
        r.shapes = vec![rect(-3.0, 0.6875, 1.0, 0.8125), rect(1.1875, 1.0, 1.3125, 5.0)];
        r.monitors = vec![Monitor { x0: 1.85, y0: 2.1, x1: 0.65, y1: 2.1 }];
        r.density = Some(pixels(rho));
        r
    };
    check(make, |o| channel(o, "monitors")[0], "gradient", &[0, 40, 77, 143]);
}

#[test]
fn a_line_sources_gradient_matches_finite_differences() {
    let make = |rho: &[f64]| {
        let mut r = open(Source2d::Line { x: 1.25, y: 0.75 }, Some(50.0));
        r.density = Some(pixels(rho));
        r
    };
    check(make, |o| o.meta.fields[0].ldos.unwrap(), "gradient", &[0, 40, 77, 143]);
}

#[test]
fn a_dipoles_gradient_matches_finite_differences_even_at_its_own_pixels() {
    // At the pixels' shared corner, so pixels 65, 66, 77 and 78 hold or surround the dipole's points.
    for angle_deg in [0.0, 90.0, 30.0] {
        let make = |rho: &[f64]| {
            let mut r = open(Source2d::Dipole { x: 1.25, y: 0.75, angle_deg }, Some(50.0));
            r.density = Some(pixels(rho));
            r
        };
        let chosen = [0, 40, 65, 66, 77, 78, 143];
        check(make, |o| o.meta.fields[0].ldos.unwrap(), "gradient", &chosen);
        // The imaginary part too, which holds a resonance's frequency still.
        check(make, |o| 2.0 * o.meta.fields[0].response.unwrap()[1], "gradient_im", &chosen);
    }
}

#[test]
fn a_dipole_in_vacuum_reads_one() {
    // Its exact `F` is `k / 2` once scaled. The derivative across two points reads a little low.
    for angle_deg in [0.0, 45.0, 90.0] {
        let r = open(Source2d::Dipole { x: 1.25, y: 1.25, angle_deg }, None);
        let o = run(&r);
        let k = o.meta.fields[0].k0;
        close(o.meta.fields[0].ldos.unwrap() / (k / 2.0), 1.0, 0.02);
    }
}

#[test]
fn the_power_a_dipole_supplies_crosses_any_loop() {
    // Energy balance, which knows nothing of how the source was scaled, beside a rod.
    let (w, m) = (2.5, 0.33);
    let mut r = open(Source2d::Dipole { x: 1.7, y: 1.3, angle_deg: 30.0 }, None);
    r.shapes = vec![Shape { outline: Outline::Circle { x: 1.1, y: 1.2, r: 0.3 }, material: CORE }];
    r.monitors = vec![
        Monitor { x0: w - m, y0: m, x1: w - m, y1: w - m },
        Monitor { x0: w - m, y0: w - m, x1: m, y1: w - m },
        Monitor { x0: m, y0: w - m, x1: m, y1: m },
        Monitor { x0: m, y0: m, x1: w - m, y1: m },
    ];
    let o = run(&r);
    let total: f64 = channel(&o, "monitors").iter().sum();
    close(total, o.meta.fields[0].ldos.unwrap() / 4.0, 5e-3 * total);
}

#[test]
fn a_finished_design_is_the_same_as_its_drawn_shapes() {
    // Whole pixels, including the dipole's, so the two make the very same matrix.
    let on = |i: usize, j: usize| (3..6).contains(&i) && (4..9).contains(&j);
    let rho: Vec<f64> = (0..144).map(|n| if on(n % 12, n / 12) { 1.0 } else { 0.0 }).collect();
    let source = Source2d::Dipole { x: 1.25, y: 0.75, angle_deg: 0.0 };
    let mut with_pixels = open(source, None);
    with_pixels.density = Some(pixels(&rho));
    let h = 0.5 / 12.0;
    let mut drawn = open(source, None);
    drawn.shapes = vec![rect(1.0 + 3.0 * h, 0.5 + 4.0 * h, 1.0 + 6.0 * h, 0.5 + 9.0 * h)];
    let read = |r: &Request2d| run(r).meta.fields[0].ldos.unwrap();
    close(read(&with_pixels), read(&drawn), 1e-9);
}

#[test]
fn a_dipole_needs_hz_and_monitors_stay_out_of_hz_pixels() {
    let mut wrong = open(Source2d::Dipole { x: 1.25, y: 0.75, angle_deg: 0.0 }, None);
    wrong.polarization = Polarization::Te;
    assert!(solve(&wrong).err().unwrap().contains("needs H_z"));
    let mut through = open(Source2d::Port { x: 0.4, y0: 0.0, y1: 1.5 }, None);
    through.shapes = vec![rect(-3.0, 0.6875, 1.0, 0.8125)];
    through.monitors = vec![Monitor { x0: 1.25, y0: 0.75, x1: 2.0, y1: 0.75 }];
    through.density = Some(pixels(&start()));
    assert!(solve(&through).err().unwrap().contains("outside a density map"));
}
