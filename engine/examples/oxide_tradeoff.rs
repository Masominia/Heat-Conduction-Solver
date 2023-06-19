//! Efficiency against speed as the buried oxide gets thicker.
//!
//!     cargo run --release --example oxide_tradeoff

use heatphase::design::analyse::analyse;
use heatphase::design::parse::parse;

fn main() {
    println!("{:>8} {:>12} {:>10} {:>16}", "box, um", "P_pi, mW", "tau, us", "P_pi x tau, nJ");
    for b in [1.0, 2.0, 3.0, 4.0] {
        let text = format!(
            "grid dx=0.25 dy=0.25\nwindow width=60 ymin=-10 ymax={top}\nlayer si -10 0\nrect si -0.25 0.25 {b} {c}\nrect tin -1 1 {h} {h2}\n\
             heater h xmin=-1 xmax=1 ymin={h} ymax={h2} power=0.1\nprobe wg xmin=-0.25 xmax=0.25 ymin={b} ymax={c}\ntransient heater=h dt=0.5 steps=500\n",
            top = b + 3.0,
            c = b + 0.25,
            h = b + 2.0,
            h2 = b + 2.25
        );
        let a = analyse(&parse(&text).unwrap()).unwrap();
        let tau = a.transient.as_ref().unwrap().tau[0];
        println!("{b:>8.0} {:>12.2} {tau:>10.2} {:>16.1}", a.pi_power[0][0], a.pi_power[0][0] * tau);
    }
}
