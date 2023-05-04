use heatphase::design::analyse::{analyse, build_mesh};
use heatphase::design::parse::parse;
use heatphase::thermal::material::SILICON;

const STACK: &str = "grid dx=0.25 dy=0.25\nwindow width=40 ymin=-5 ymax=5\nlayer si -5 0\nrect si -0.25 0.25 2 2.25\n";
const HEATER: &str = "heater h1 xmin=-1 xmax=1 ymin=4 ymax=4.25 power=0.1\nprobe wg xmin=-0.25 xmax=0.25 ymin=2 ymax=2.25\nprobe far x=10.125 y=2.125\n";

fn run(text: &str) -> heatphase::design::analyse::Analysis {
    analyse(&parse(text).unwrap()).unwrap()
}

#[test]
fn mesh_follows_the_file() {
    let m = build_mesh(&parse(STACK).unwrap()).unwrap();
    assert_eq!((m.nx, m.ny), (160, 40));
    assert_eq!(m.cells.iter().filter(|c| **c == SILICON).count(), 160 * 20 + 2);
}

#[test]
fn response_temperature_and_pi_power_are_consistent() {
    let a = run(&format!("{STACK}{HEATER}"));
    assert!(a.response[0][0] > 100.0 && a.response[0][0] < 250.0);
    assert!((a.temperature[0] - 0.1 * a.response[0][0]).abs() < 1e-12);
    assert!((a.pi_power[0][0] - 1.55 / (2.0 * 1.86e-4 * a.response[0][0])).abs() < 1e-9);
    assert!(a.pi_power[0][0] > 15.0 && a.pi_power[0][0] < 40.0);
    assert!(a.response[1][0] < 0.2 * a.response[0][0]);
    // the hottest point is in the heater
    assert!(a.peak > a.temperature[0] && a.peak_position.0.abs() < 1.0 && (a.peak_position.1 - 4.125).abs() < 0.01);
    assert!((a.injected - 100.0).abs() < 1e-9 && ((a.to_sink + a.to_surface) / a.injected - 1.0).abs() < 1e-10);
    assert_eq!(a.iterations, 0);
}

#[test]
fn fields_superpose() {
    let two = "heater h2 xmin=4 xmax=6 ymin=4 ymax=4.25 power=0.05\n";
    let both = run(&format!("{STACK}{HEATER}{two}"));
    let first = run(&format!("{STACK}{HEATER}"));
    let second = run(&format!("{STACK}{}{two}", HEATER.replace("power=0.1", "power=0")));
    for p in 0..2 {
        assert!((both.temperature[p] - first.temperature[p] - second.temperature[p]).abs() < 1e-9);
    }
    assert!((both.injected - 150.0).abs() < 1e-9);
}

#[test]
fn doubling_the_power_doubles_the_field_but_not_the_response() {
    let a = run(&format!("{STACK}{HEATER}"));
    let b = run(&format!("{STACK}{}", HEATER.replace("power=0.1", "power=0.2")));
    assert!((b.peak / a.peak - 2.0).abs() < 1e-10);
    assert!((b.response[0][0] - a.response[0][0]).abs() < 1e-9);
}

#[test]
fn conjugate_gradients_agree_and_count_iterations() {
    let a = run(&format!("{STACK}{HEATER}"));
    let b = run(&format!("{STACK}{HEATER}solver cg\n"));
    assert!((b.response[0][0] / a.response[0][0] - 1.0).abs() < 1e-8);
    assert!(b.iterations > 20);
}

#[test]
fn transient_reports_delays() {
    let a = run(&format!("{STACK}{HEATER}transient heater=h1 dt=1 steps=300\n"));
    let t = a.transient.unwrap();
    assert_eq!((t.heater.as_str(), t.time.len(), t.traces.len()), ("h1", 300, 2));
    assert!(t.tau[0] > 2.0 && t.tau[0] < 30.0);
    assert!(t.tau[1] > t.tau[0] && t.rise[0] > t.tau[0]);
    // a run that is too short cannot give a time constant
    let short = run(&format!("{STACK}{HEATER}transient heater=h1 dt=0.1 steps=5\n"));
    assert!(short.transient.unwrap().tau[0].is_nan());
}

#[test]
fn impossible_files_are_refused() {
    for (text, needle) in [
        (STACK.to_string(), "no heater"),
        (format!("{STACK}heater h xmin=50 xmax=52 ymin=4 ymax=4.25\n"), "outside the window"),
        (format!("{STACK}{HEATER}probe lost x=99 y=0\n"), "probe 'lost' lies outside"),
        (format!("window width=40.1 ymin=-5 ymax=5\n{HEATER}"), "whole number of cells"),
    ] {
        let err = analyse(&parse(&text).unwrap()).err().unwrap();
        assert!(err.contains(needle), "{err}");
    }
}
