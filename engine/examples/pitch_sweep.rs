//! Crosstalk between two phase shifters against their pitch.
//!
//!     cargo run --release --example pitch_sweep

use heatphase::design::analyse::analyse;
use heatphase::design::parse::parse;

fn main() {
    println!("{:>10} {:>14} {:>14} {:>10}", "pitch, um", "self, K/(mW/um)", "neighbour", "crosstalk");
    for pitch in [4.0, 6.0, 8.0, 10.0, 15.0, 20.0, 30.0] {
        let text = format!(
            "grid dx=0.25 dy=0.25\nwindow width=80 ymin=-10 ymax=5\nlayer si -10 0\nrect si -0.25 0.25 2 2.25\nrect tin -1 1 4 4.25\n\
             heater h xmin=-1 xmax=1 ymin=4 ymax=4.25 power=0.1\nprobe wg xmin=-0.25 xmax=0.25 ymin=2 ymax=2.25\n\
             probe next xmin={} xmax={} ymin=2 ymax=2.25\n",
            pitch - 0.25,
            pitch + 0.25
        );
        let a = analyse(&parse(&text).unwrap()).unwrap();
        println!("{pitch:>10.0} {:>14.2} {:>14.3} {:>9.2}%", a.response[0][0], a.response[1][0], 100.0 * a.response[1][0] / a.response[0][0]);
    }
}
