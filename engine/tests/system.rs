use heatphase::thermal::material::{Material, SILICA, SILICON};
use heatphase::thermal::mesh::{Mesh, Rect};
use heatphase::thermal::system::{assemble, bottom_conductance, capacities, outflow, top_conductance, Ordering, Top};

#[test]
fn numbering_follows_the_shorter_side() {
    let wide = Mesh::new(10.0, 0.0, 2.0, 0.5, 0.5, SILICA).unwrap();
    let o = Ordering::new(&wide);
    assert!(!o.by_rows && o.half_width() == 4 && o.index(2, 1) == 9);
    let tall = Mesh::new(2.0, 0.0, 10.0, 0.5, 0.5, SILICA).unwrap();
    let o = Ordering::new(&tall);
    assert!(o.by_rows && o.half_width() == 4 && o.index(1, 2) == 9);
    let v: Vec<f64> = (0..wide.len()).map(|i| i as f64).collect();
    let round_trip = Ordering::new(&wide).to_mesh(&Ordering::new(&wide).from_mesh(&v));
    assert_eq!(round_trip, v);
}

#[test]
fn conductances_between_cells() {
    let mut m = Mesh::new(2.0, 0.0, 1.0, 0.5, 0.25, SILICA).unwrap();
    m.paint(&Rect { xmin: 0.0, xmax: 1.0, ymin: 0.0, ymax: 1.0 }, SILICON);
    let a = assemble(&m, Top::Adiabatic);
    let o = Ordering::new(&m);
    // oxide to oxide, sideways: k dy / dx
    assert!((a.get(o.index(0, 1), o.index(1, 1)) + 1.38 * 0.25 / 0.5).abs() < 1e-12);
    // oxide to silicon: harmonic mean
    let harmonic = 2.0 * 1.38 * 148.0 / (1.38 + 148.0);
    assert!((a.get(o.index(1, 1), o.index(2, 1)) + harmonic * 0.5).abs() < 1e-12);
    // upwards: k dx / dy
    assert!((a.get(o.index(0, 1), o.index(0, 2)) + 1.38 * 2.0).abs() < 1e-12);
    // rows of an interior cell sum to zero; the bottom row carries the sink
    let row_sum = |ix: usize, iy: usize| -> f64 { (0..m.len()).map(|j| a.get(o.index(ix, iy), j)).sum() };
    assert!(row_sum(1, 2).abs() < 1e-9);
    assert!((row_sum(1, 0) - bottom_conductance(&m, 1)).abs() < 1e-9);
    assert!((bottom_conductance(&m, 0) - 1.38 * 0.5 / 0.125).abs() < 1e-12);
}

#[test]
fn film_in_series_with_half_a_cell() {
    let m = Mesh::new(2.0, 0.0, 1.0, 0.5, 0.25, SILICA).unwrap();
    assert_eq!(top_conductance(&m, 0, Top::Adiabatic), 0.0);
    let h = 1000.0;
    let expected = 1.0 / (1.0 / (h * 0.5e-6) + 0.125 / (1.38 * 0.5));
    assert!((top_conductance(&m, 0, Top::Convective(h)) - expected).abs() < 1e-18);
    // a very large film coefficient leaves only the half cell
    assert!((top_conductance(&m, 0, Top::Convective(1e15)) - 1.38 * 0.5 / 0.125).abs() < 1e-6);
}

#[test]
fn capacities_and_outflow() {
    let m = Mesh::new(2.0, 0.0, 1.0, 0.5, 0.25, Material { conductivity: 2.0, capacity: 1e6 }).unwrap();
    assert!(capacities(&m).iter().all(|c| (c - 1e6 * 0.125e-12).abs() < 1e-20));
    let field = vec![3.0; m.len()];
    let (sink, surface) = outflow(&m, Top::Adiabatic, &field);
    assert!((sink - 4.0 * 3.0 * 2.0 * 0.5 / 0.125).abs() < 1e-9);
    assert_eq!(surface, 0.0);
    assert!(outflow(&m, Top::Convective(500.0), &field).1 > 0.0);
}
