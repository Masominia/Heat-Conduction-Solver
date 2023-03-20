//! Command line front end.
//!
//!     heatphase run <cross-section>     results, checks, the temperature field and the step response
//!     heatphase check <cross-section>   results and checks; exit status 1 if a check fails
//!
//! Output lines: `key value`, `probe <name> <K>`, `response <probe> <heater> <K per mW/um> <pi power mW>`,
//! `step <probe> <tau us> <rise us>`, `check <PASS|FAIL> <name> :: <detail>`, and for `run`:
//! `conductivity <row...>` and `field <row...>` per mesh row from the bottom, `trace <t us> <K per probe...>`.

use std::env;
use std::fs;
use std::process;

use heatphase::design::analyse::analyse;
use heatphase::design::expect::evaluate;
use heatphase::design::parse::parse;

fn fail(msg: &str) -> ! {
    eprintln!("heatphase: {msg}");
    process::exit(2)
}

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    if args.len() != 2 || (args[0] != "run" && args[0] != "check") {
        fail("usage: heatphase <run|check> <cross-section file>");
    }
    let text = fs::read_to_string(&args[1]).unwrap_or_else(|e| fail(&format!("cannot read {}: {e}", args[1])));
    let spec = parse(&text).unwrap_or_else(|e| fail(&e));
    let a = analyse(&spec).unwrap_or_else(|e| fail(&e));

    println!("title {}", spec.title);
    println!("nx {}\nny {}\ndx {}\ndy {}\nx0 {}\ny0 {}", a.mesh.nx, a.mesh.ny, a.mesh.dx, a.mesh.dy, a.mesh.x0, a.mesh.y0);
    println!("peak {:.10}\npeak_x {:.4}\npeak_y {:.4}", a.peak, a.peak_position.0, a.peak_position.1);
    println!("injected {:.10}\nto_sink {:.10}\nto_surface {:.10}\niterations {}", a.injected, a.to_sink, a.to_surface, a.iterations);
    for (p, probe) in spec.probes.iter().enumerate() {
        println!("probe {} {:.12e}", probe.name, a.temperature[p]);
        for (h, heater) in spec.heaters.iter().enumerate() {
            println!("response {} {} {:.12e} {:.10e}", probe.name, heater.name, a.response[p][h], a.pi_power[p][h]);
        }
    }
    if let Some(t) = &a.transient {
        for (p, probe) in spec.probes.iter().enumerate() {
            println!("step {} {:.8} {:.8}", probe.name, t.tau[p], t.rise[p]);
        }
    }
    let mut failed = 0;
    for c in evaluate(&spec, &a) {
        println!("check {} {} :: {}", if c.pass { "PASS" } else { "FAIL" }, c.name, c.detail);
        failed += usize::from(!c.pass);
    }
    if args[0] == "check" {
        process::exit(if failed == 0 { 0 } else { 1 });
    }
    for iy in 0..a.mesh.ny {
        let row = &a.mesh.cells[iy * a.mesh.nx..(iy + 1) * a.mesh.nx];
        println!("conductivity {}", row.iter().map(|m| format!("{}", m.conductivity)).collect::<Vec<_>>().join(" "));
        let row = &a.field[iy * a.mesh.nx..(iy + 1) * a.mesh.nx];
        println!("field {}", row.iter().map(|t| format!("{t:.6e}")).collect::<Vec<_>>().join(" "));
    }
    if let Some(t) = &a.transient {
        let stride = (t.time.len() / 500).max(1);
        for i in (0..t.time.len()).step_by(stride) {
            println!("trace {:.6} {}", t.time[i], t.traces.iter().map(|tr| format!("{:.8e}", tr[i])).collect::<Vec<_>>().join(" "));
        }
    }
}
