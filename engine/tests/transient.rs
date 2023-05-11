use heatphase::analytic::slab_time_constant;
use heatphase::thermal::material::SILICA;
use heatphase::thermal::mesh::{Mesh, Rect};
use heatphase::thermal::steady::{source, Conduction};
use heatphase::thermal::system::Top;
use heatphase::thermal::transient::{crossing, step_response};

#[test]
fn crossing_interpolates() {
    let (t, v) = ([1.0, 2.0, 3.0], [0.2, 0.6, 1.0]);
    assert!((crossing(&t, &v, 1.0, 0.4).unwrap() - 1.5).abs() < 1e-12);
    assert!((crossing(&t, &v, 1.0, 0.1).unwrap() - 0.5).abs() < 1e-12);
    assert_eq!(crossing(&t, &v, 2.0, 0.9), None);
}

#[test]
fn the_step_response_settles_on_the_steady_field() {
    let m = Mesh::new(4.0, -2.0, 0.0, 0.5, 0.125, SILICA).unwrap();
    let q = source(&m, &Rect { xmin: -2.0, xmax: 2.0, ymin: -0.125, ymax: 0.0 }, 0.1);
    let steady = Conduction::new(&m, Top::Adiabatic).unwrap().solve(&q);
    let cell = m.index(3, m.ny - 1);
    let run = step_response(&m, Top::Adiabatic, &q, &[vec![(cell, 1.0)]], 0.5, 200).unwrap();
    assert_eq!((run.time.len(), run.traces[0].len()), (200, 200));
    assert!((run.time[199] - 100.0).abs() < 1e-12);
    // monotone rise to the steady value
    assert!(run.traces[0].windows(2).all(|w| w[1] >= w[0]));
    assert!((run.traces[0][199] / steady[cell] - 1.0).abs() < 1e-6);
}

#[test]
fn the_slowest_decay_matches_the_slab_time_constant() {
    let d = 2.0;
    let m = Mesh::new(2.0, -d, 0.0, 0.5, 0.0625, SILICA).unwrap();
    let q = source(&m, &Rect { xmin: -1.0, xmax: 1.0, ymin: -0.0625, ymax: 0.0 }, 0.1);
    let steady = Conduction::new(&m, Top::Adiabatic).unwrap().solve(&q);
    let cell = m.index(1, m.ny - 1);
    let dt = 0.01;
    let run = step_response(&m, Top::Adiabatic, &q, &[vec![(cell, 1.0)]], dt, 1000).unwrap();
    // late in the response only the slowest mode is left: the gap to the steady value decays as exp(-t / tau)
    let gap = |i: usize| steady[cell] - run.traces[0][i];
    let tau = (run.time[900] - run.time[500]) / (gap(500) / gap(900)).ln();
    let expected = slab_time_constant(d, SILICA.diffusivity());
    assert!((tau / expected - 1.0).abs() < 0.01, "{tau} {expected}");
}

#[test]
fn bad_arguments() {
    let m = Mesh::new(2.0, -1.0, 0.0, 0.5, 0.25, SILICA).unwrap();
    let q = vec![0.0; m.len()];
    assert!(step_response(&m, Top::Adiabatic, &q, &[], 0.0, 10).is_err());
    assert!(step_response(&m, Top::Adiabatic, &q, &[], 1.0, 0).is_err());
}
