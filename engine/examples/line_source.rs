//! The solver against the closed form for a line heater on a slab.
//!
//!     cargo run --release --example line_source

use heatphase::analytic::line_source_surface;
use heatphase::thermal::material::SILICA;
use heatphase::thermal::mesh::{Mesh, Rect};
use heatphase::thermal::steady::{probe, source, Conduction};
use heatphase::thermal::system::Top;

fn main() {
    let d = 4.0;
    let mesh = Mesh::new(64.0, -d, 0.0, 0.125, 0.125, SILICA).unwrap();
    let conduction = Conduction::new(&mesh, Top::Adiabatic).unwrap();
    let field = conduction.solve(&source(&mesh, &Rect { xmin: -0.125, xmax: 0.125, ymin: -0.125, ymax: 0.0 }, 0.1));
    println!("{:>8} {:>12} {:>12} {:>8}", "x, um", "solver, K", "formula, K", "ratio");
    for k in 1..=12 {
        let x = k as f64 + 0.0625;
        let got = probe(&mesh, &Rect { xmin: x, xmax: x, ymin: -0.0625, ymax: -0.0625 }, &field).unwrap();
        let want = line_source_surface(0.1, 1.38, d, x);
        println!("{x:>8.3} {got:>12.5} {want:>12.5} {:>8.4}", got / want);
    }
}
