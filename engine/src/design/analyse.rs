//! Everything computed for one cross-section.

use super::spec::{Solver, Spec};
use crate::analytic::pi_power;
use crate::thermal::mesh::Mesh;
use crate::thermal::steady::{probe, source, Conduction};
use crate::thermal::system::outflow;
use crate::thermal::transient::{crossing, step_response};

pub struct TransientResult {
    pub heater: String,
    pub time: Vec<f64>,
    /// One trace per probe.
    pub traces: Vec<Vec<f64>>,
    /// Time to 1 - 1/e of the steady value, microseconds; NaN if the run is too short.
    pub tau: Vec<f64>,
    /// Time from 10% to 90% of the steady value, microseconds.
    pub rise: Vec<f64>,
}

pub struct Analysis {
    pub mesh: Mesh,
    /// response[p][h]: temperature rise of probe p in K per (mW per um) dissipated in heater h.
    pub response: Vec<Vec<f64>>,
    /// Power for a phase shift of pi, mW, same layout.
    pub pi_power: Vec<Vec<f64>>,
    /// Field with every heater at its stated power, K above the sink, mesh order.
    pub field: Vec<f64>,
    /// Probe temperatures in that field.
    pub temperature: Vec<f64>,
    pub peak: f64,
    pub peak_position: (f64, f64),
    /// Heat entering (W / m), leaving through the sink and through the top surface.
    pub injected: f64,
    pub to_sink: f64,
    pub to_surface: f64,
    /// Iterations of the conjugate-gradient solver, summed over heaters; 0 for the direct solver.
    pub iterations: usize,
    pub transient: Option<TransientResult>,
}

pub fn build_mesh(spec: &Spec) -> Result<Mesh, String> {
    let mut mesh = Mesh::new(spec.width, spec.ymin, spec.ymax, spec.dx, spec.dy, spec.background)?;
    for (rect, material) in &spec.rects {
        mesh.paint(rect, *material);
    }
    Ok(mesh)
}

pub fn analyse(spec: &Spec) -> Result<Analysis, String> {
    if spec.heaters.is_empty() {
        return Err("the cross-section has no heater".into());
    }
    let mesh = build_mesh(spec)?;
    let conduction = Conduction::new(&mesh, spec.top)?;
    let mut iterations = 0;
    let mut unit_fields = Vec::new();
    for heater in &spec.heaters {
        let q = source(&mesh, &heater.rect, 1.0);
        if q.iter().all(|v| *v == 0.0) {
            return Err(format!("heater '{}' lies outside the window", heater.name));
        }
        unit_fields.push(match spec.solver {
            Solver::Direct => conduction.solve(&q),
            Solver::ConjugateGradient => {
                let (field, n) = conduction.solve_cg(&q, 1e-13)?;
                iterations += n;
                field
            }
        });
    }
    let mut response = Vec::new();
    for p in &spec.probes {
        let row: Result<Vec<f64>, String> = unit_fields.iter().map(|f| probe(&mesh, &p.rect, f).ok_or_else(|| format!("probe '{}' lies outside the window", p.name))).collect();
        response.push(row?);
    }
    let mut field = vec![0.0; mesh.len()];
    let mut injected = 0.0;
    for (heater, unit) in spec.heaters.iter().zip(&unit_fields) {
        for (t, u) in field.iter_mut().zip(unit) {
            *t += heater.power * u;
        }
        injected += source(&mesh, &heater.rect, heater.power).iter().sum::<f64>();
    }
    let temperature: Vec<f64> = response.iter().map(|row| row.iter().zip(&spec.heaters).map(|(r, h)| r * h.power).sum()).collect();
    let (mut peak, mut at) = (f64::MIN, 0);
    for (i, t) in field.iter().enumerate() {
        if *t > peak {
            peak = *t;
            at = i;
        }
    }
    let (to_sink, to_surface) = outflow(&mesh, spec.top, &field);

    let transient = match &spec.transient {
        None => None,
        Some(t) => {
            let h = spec.heaters.iter().position(|x| x.name == t.heater).ok_or("the transient names a heater that does not exist")?;
            let heater = &spec.heaters[h];
            let mut cells = Vec::new();
            for p in &spec.probes {
                let covered = mesh.coverage(&p.rect);
                let total: f64 = covered.iter().map(|c| c.1).sum();
                cells.push(if total > 0.0 {
                    covered.iter().map(|(c, a)| (*c, a / total)).collect()
                } else {
                    vec![(mesh.locate(p.rect.xmin, p.rect.ymin).ok_or("a probe lies outside the window")?, 1.0)]
                });
            }
            let run = step_response(&mesh, spec.top, &source(&mesh, &heater.rect, heater.power), &cells, t.dt, t.steps)?;
            let (mut tau, mut rise) = (Vec::new(), Vec::new());
            for (p, trace) in run.traces.iter().enumerate() {
                let steady = heater.power * response[p][h];
                tau.push(crossing(&run.time, trace, steady, 1.0 - (-1.0f64).exp()).unwrap_or(f64::NAN));
                rise.push(match (crossing(&run.time, trace, steady, 0.1), crossing(&run.time, trace, steady, 0.9)) {
                    (Some(a), Some(b)) => b - a,
                    _ => f64::NAN,
                });
            }
            Some(TransientResult { heater: t.heater.clone(), time: run.time, traces: run.traces, tau, rise })
        }
    };

    let pi = response.iter().map(|row| row.iter().map(|r| pi_power(spec.wavelength, spec.dn_dt, *r)).collect()).collect();
    Ok(Analysis {
        peak_position: (mesh.x(at % mesh.nx), mesh.y(at / mesh.nx)),
        mesh,
        response,
        pi_power: pi,
        field,
        temperature,
        peak,
        injected,
        to_sink,
        to_surface,
        iterations,
        transient,
    })
}
