//! Expectations of a cross-section file, checked against its analysis.
//!
//! ```text
//! temperature <probe> <K> <reltol>              with every heater at its stated power
//! response <probe> <heater> <K per mW/um> <reltol>
//! pi_power <probe> <heater> <mW> <reltol>
//! crosstalk <probe> <reference probe> <heater> <ratio> <reltol>   response of one probe over another's
//! hotter <probe> <probe>                        the first is warmer than the second
//! peak <K> <reltol>                             highest temperature anywhere
//! balance <tol>                                 heat out equals heat in, relative
//! surface_share <fraction> <tol>                share of the heat leaving through the top
//! tau <probe> <us> <reltol>                     time to 1 - 1/e of the step response
//! rise <probe> <us> <reltol>                    time from 10% to 90%
//! cells <count>
//! ```

use super::analyse::Analysis;
use super::spec::{Expectation, Spec};

#[derive(Clone, Debug, PartialEq)]
pub struct Check {
    pub name: String,
    pub pass: bool,
    pub detail: String,
}

fn arity(e: &Expectation, n: usize) -> Result<(), String> {
    if e.args.len() == n {
        Ok(())
    } else {
        Err(format!("'{}' takes {n} argument{}", e.kind, if n == 1 { "" } else { "s" }))
    }
}

fn num(e: &Expectation, i: usize) -> Result<f64, String> {
    e.args[i].parse::<f64>().map_err(|_| format!("'{}' is not a number", e.args[i]))
}

fn probe(spec: &Spec, name: &str) -> Result<usize, String> {
    spec.probes.iter().position(|p| p.name == name).ok_or_else(|| format!("no probe named '{name}'"))
}

fn heater(spec: &Spec, name: &str) -> Result<usize, String> {
    spec.heaters.iter().position(|h| h.name == name).ok_or_else(|| format!("no heater named '{name}'"))
}

fn relative(got: f64, want: f64, tol: f64) -> (bool, String) {
    ((got / want - 1.0).abs() <= tol, format!("{got:.8}, expected {want:.8}, relative error {:.2e}", got / want - 1.0))
}

fn one(spec: &Spec, a: &Analysis, e: &Expectation) -> Result<(bool, String), String> {
    match e.kind.as_str() {
        "temperature" => {
            arity(e, 3)?;
            Ok(relative(a.temperature[probe(spec, &e.args[0])?], num(e, 1)?, num(e, 2)?))
        }
        "response" | "pi_power" => {
            arity(e, 4)?;
            let (p, h) = (probe(spec, &e.args[0])?, heater(spec, &e.args[1])?);
            let got = if e.kind == "response" { a.response[p][h] } else { a.pi_power[p][h] };
            Ok(relative(got, num(e, 2)?, num(e, 3)?))
        }
        "crosstalk" => {
            arity(e, 5)?;
            let (p, r, h) = (probe(spec, &e.args[0])?, probe(spec, &e.args[1])?, heater(spec, &e.args[2])?);
            Ok(relative(a.response[p][h] / a.response[r][h], num(e, 3)?, num(e, 4)?))
        }
        "hotter" => {
            arity(e, 2)?;
            let (p, q) = (a.temperature[probe(spec, &e.args[0])?], a.temperature[probe(spec, &e.args[1])?]);
            Ok((p > q, format!("{p:.6} K against {q:.6} K")))
        }
        "peak" => {
            arity(e, 2)?;
            Ok(relative(a.peak, num(e, 0)?, num(e, 1)?))
        }
        "balance" => {
            arity(e, 1)?;
            let error = ((a.to_sink + a.to_surface) / a.injected - 1.0).abs();
            Ok((error <= num(e, 0)?, format!("{:.6} W/m in, {:.6} W/m out, relative difference {error:.2e}", a.injected, a.to_sink + a.to_surface)))
        }
        "surface_share" => {
            arity(e, 2)?;
            let share = a.to_surface / (a.to_sink + a.to_surface);
            Ok(((share - num(e, 0)?).abs() <= num(e, 1)?, format!("{share:.6} of the heat leaves through the top")))
        }
        "tau" | "rise" => {
            arity(e, 3)?;
            let t = a.transient.as_ref().ok_or_else(|| format!("'{}' needs a transient statement", e.kind))?;
            let p = probe(spec, &e.args[0])?;
            let got = if e.kind == "tau" { t.tau[p] } else { t.rise[p] };
            if got.is_nan() {
                return Ok((false, "the transient is too short to reach that level".into()));
            }
            Ok(relative(got, num(e, 1)?, num(e, 2)?))
        }
        "cells" => {
            arity(e, 1)?;
            Ok((a.mesh.len() as f64 == num(e, 0)?, format!("{} cells ({} by {})", a.mesh.len(), a.mesh.nx, a.mesh.ny)))
        }
        other => Err(format!("unknown expectation '{other}'")),
    }
}

/// A malformed expectation becomes a failed check with its line number; it does not stop the others.
pub fn evaluate(spec: &Spec, analysis: &Analysis) -> Vec<Check> {
    spec.expectations
        .iter()
        .map(|e| {
            let name = std::iter::once(e.kind.clone()).chain(e.args.iter().cloned()).collect::<Vec<_>>().join(" ");
            match one(spec, analysis, e) {
                Ok((pass, detail)) => Check { name, pass, detail },
                Err(msg) => Check { name, pass: false, detail: format!("line {}: {msg}", e.line) },
            }
        })
        .collect()
}
