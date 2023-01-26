//! Step response by implicit Euler time stepping.

use super::mesh::Mesh;
use super::system::{assemble, capacities, Ordering, Top};

pub struct StepResponse {
    /// Times in microseconds, starting at the first step.
    pub time: Vec<f64>,
    /// One row per probe cell set: mean temperature at each time.
    pub traces: Vec<Vec<f64>>,
}

/// Switches a source on at t = 0 from a cold start and records weighted means of the temperature.
/// `probes` are lists of (cell, weight) with weights summing to one. `dt` is in microseconds.
pub fn step_response(mesh: &Mesh, top: Top, source: &[f64], probes: &[Vec<(usize, f64)>], dt: f64, steps: usize) -> Result<StepResponse, String> {
    if !(dt > 0.0) || steps == 0 {
        return Err("the transient needs a positive step and at least one step".into());
    }
    let order = Ordering::new(mesh);
    let c_over_dt: Vec<f64> = capacities(mesh).iter().map(|c| c / (dt * 1e-6)).collect();
    // (C / dt + A) T_new = C / dt T_old + q
    let mut lhs = assemble(mesh, top);
    for iy in 0..mesh.ny {
        for ix in 0..mesh.nx {
            let k = order.index(ix, iy);
            lhs.add(k, k, c_over_dt[mesh.index(ix, iy)]);
        }
    }
    lhs.factorise()?;
    let mut temperature = vec![0.0; mesh.len()];
    let mut out = StepResponse { time: Vec::with_capacity(steps), traces: vec![Vec::with_capacity(steps); probes.len()] };
    for step in 1..=steps {
        let rhs: Vec<f64> = (0..mesh.len()).map(|i| c_over_dt[i] * temperature[i] + source[i]).collect();
        temperature = order.to_mesh(&lhs.solve(&order.from_mesh(&rhs)));
        out.time.push(step as f64 * dt);
        for (trace, cells) in out.traces.iter_mut().zip(probes) {
            trace.push(cells.iter().map(|(cell, w)| temperature[*cell] * w).sum());
        }
    }
    Ok(out)
}

/// First time at which a rising trace reaches `fraction` of `final_value`, by linear interpolation; the trace
/// starts from zero at t = 0. None if it never gets there.
pub fn crossing(time: &[f64], trace: &[f64], final_value: f64, fraction: f64) -> Option<f64> {
    let level = fraction * final_value;
    let (mut t0, mut v0) = (0.0, 0.0);
    for (t, v) in time.iter().zip(trace) {
        if *v >= level {
            return Some(if *v > v0 { t0 + (level - v0) / (v - v0) * (t - t0) } else { *t });
        }
        t0 = *t;
        v0 = *v;
    }
    None
}
