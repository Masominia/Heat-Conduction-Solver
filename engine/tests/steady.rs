use heatphase::analytic::{line_source_surface, stack_temperature};
use heatphase::thermal::material::{Material, SILICA, SILICON};
use heatphase::thermal::mesh::{Mesh, Rect};
use heatphase::thermal::steady::{probe, source, Conduction};
use heatphase::thermal::system::{outflow, Top};

#[test]
fn source_carries_the_stated_power() {
    let m = Mesh::new(10.0, 0.0, 2.0, 0.25, 0.25, SILICA).unwrap();
    let q = source(&m, &Rect { xmin: -0.9, xmax: 1.1, ymin: 1.3, ymax: 1.6 }, 0.2);
    assert!((q.iter().sum::<f64>() - 200.0).abs() < 1e-9);
    // half outside the window: half the power is lost
    let q = source(&m, &Rect { xmin: 4.0, xmax: 6.0, ymin: 1.0, ymax: 2.0 }, 0.2);
    assert!((q.iter().sum::<f64>() - 100.0).abs() < 1e-9);
}

#[test]
fn one_dimensional_stack_is_exact() {
    // silicon 0..1, oxide 1..3, heated in the top row over the whole width
    let mut m = Mesh::new(6.0, 0.0, 3.0, 0.5, 0.25, SILICA).unwrap();
    m.paint(&Rect { xmin: -9.0, xmax: 9.0, ymin: 0.0, ymax: 1.0 }, SILICON);
    let c = Conduction::new(&m, Top::Adiabatic).unwrap();
    let t = c.solve(&source(&m, &Rect { xmin: -3.0, xmax: 3.0, ymin: 2.75, ymax: 3.0 }, 0.3));
    // at the centre of the row just under the heater: all the silicon and 1.625 um of oxide below
    let expected = stack_temperature(&[(1.0, 148.0), (1.625, 1.38)], 0.3, 6.0);
    let got = t[m.index(4, 10)];
    assert!((got / expected - 1.0).abs() < 1e-12, "{got} {expected}");
    // the same in every column
    assert!((t[m.index(0, 10)] - got).abs() < 1e-10);
}

#[test]
fn heat_in_equals_heat_out() {
    let mut m = Mesh::new(20.0, -4.0, 4.0, 0.25, 0.25, SILICA).unwrap();
    m.paint(&Rect { xmin: -99.0, xmax: 99.0, ymin: -4.0, ymax: 0.0 }, SILICON);
    for top in [Top::Adiabatic, Top::Convective(2500.0)] {
        let c = Conduction::new(&m, top).unwrap();
        let t = c.solve(&source(&m, &Rect { xmin: -1.0, xmax: 1.0, ymin: 3.0, ymax: 3.25 }, 0.1));
        let (sink, surface) = outflow(&m, top, &t);
        assert!(((sink + surface) / 100.0 - 1.0).abs() < 1e-10);
        assert!(t.iter().all(|v| *v > 0.0));
    }
}

#[test]
fn conjugate_gradients_give_the_same_field() {
    let mut m = Mesh::new(12.0, -2.0, 3.0, 0.25, 0.25, SILICA).unwrap();
    m.paint(&Rect { xmin: -99.0, xmax: 99.0, ymin: -2.0, ymax: 0.0 }, SILICON);
    let c = Conduction::new(&m, Top::Adiabatic).unwrap();
    let q = source(&m, &Rect { xmin: -1.0, xmax: 1.0, ymin: 2.0, ymax: 2.25 }, 0.1);
    let direct = c.solve(&q);
    let (iterative, iterations) = c.solve_cg(&q, 1e-13).unwrap();
    let peak = direct.iter().cloned().fold(0.0, f64::max);
    let worst = direct.iter().zip(&iterative).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
    assert!(worst < 1e-9 * peak && iterations > 10);
    assert_eq!(c.matrix().size(), m.len());
}

#[test]
fn line_heater_follows_the_closed_form() {
    let d = 2.0;
    let m = Mesh::new(40.0, -d, 0.0, 0.125, 0.125, SILICA).unwrap();
    let c = Conduction::new(&m, Top::Adiabatic).unwrap();
    let t = c.solve(&source(&m, &Rect { xmin: -0.125, xmax: 0.125, ymin: -0.125, ymax: 0.0 }, 0.1));
    for x in [1.0625, 2.0625, 4.0625] {
        let got = probe(&m, &Rect { xmin: x, xmax: x, ymin: -0.0625, ymax: -0.0625 }, &t).unwrap();
        let expected = line_source_surface(0.1, 1.38, d, x);
        assert!((got / expected - 1.0).abs() < 0.03, "{x}: {got} {expected}");
    }
}

#[test]
fn probes_average_by_area_or_take_a_point() {
    let m = Mesh::new(4.0, 0.0, 2.0, 1.0, 1.0, Material { conductivity: 1.0, capacity: 1.0 }).unwrap();
    let field: Vec<f64> = (0..m.len()).map(|i| i as f64).collect();
    // half of cell 1 and half of cell 2
    let mean = probe(&m, &Rect { xmin: -0.5, xmax: 0.5, ymin: 0.0, ymax: 1.0 }, &field).unwrap();
    assert!((mean - 1.5).abs() < 1e-12);
    assert_eq!(probe(&m, &Rect { xmin: 1.2, xmax: 1.2, ymin: 1.5, ymax: 1.5 }, &field), Some(7.0));
    assert_eq!(probe(&m, &Rect { xmin: 9.0, xmax: 9.0, ymin: 1.5, ymax: 1.5 }, &field), None);
}
