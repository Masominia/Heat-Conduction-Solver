//! Parser for cross-section files.
//!
//! One statement per line, `#` starts a comment:
//!
//! ```text
//! title <text>
//! grid dx=<um> dy=<um>
//! window width=<um> ymin=<um> ymax=<um>
//! material <name> k=<W/(m K)> c=<J/(m^3 K)>        define or redefine a material
//! background <material>
//! layer <material> <ymin> <ymax>
//! rect <material> <xmin> <xmax> <ymin> <ymax>
//! heater <name> xmin= xmax= ymin= ymax= power=<mW per um>
//! probe <name> x= y=   |   probe <name> xmin= xmax= ymin= ymax=
//! top adiabatic | top convective h=<W/(m^2 K)>
//! solver direct | cg
//! optics wavelength=<um> dndt=<1/K>
//! transient heater=<name> dt=<us> steps=<n>
//! expect <kind> <args...>
//! ```
//!
//! Built-in materials: si, sio2, sin, tin, al, air.

use std::collections::HashMap;

use super::spec::{Expectation, Heater, Probe, Solver, Spec, Transient};
use crate::thermal::material::{by_name, Material};
use crate::thermal::mesh::Rect;
use crate::thermal::system::Top;

fn number(text: &str) -> Result<f64, String> {
    text.parse::<f64>().map_err(|_| format!("'{text}' is not a number"))
}

fn pairs(words: &[&str], allowed: &[&str], owner: &str) -> Result<HashMap<String, String>, String> {
    let mut out = HashMap::new();
    for w in words {
        let (k, v) = w.split_once('=').filter(|(k, v)| !k.is_empty() && !v.is_empty()).ok_or_else(|| format!("expected key=value, got '{w}'"))?;
        if !allowed.contains(&k) {
            return Err(format!("{owner} does not have '{k}'"));
        }
        out.insert(k.to_string(), v.to_string());
    }
    Ok(out)
}

fn need(p: &HashMap<String, String>, key: &str, owner: &str) -> Result<f64, String> {
    number(p.get(key).ok_or_else(|| format!("{owner} needs {key}=..."))?)
}

fn rect_of(p: &HashMap<String, String>, owner: &str) -> Result<Rect, String> {
    let r = Rect { xmin: need(p, "xmin", owner)?, xmax: need(p, "xmax", owner)?, ymin: need(p, "ymin", owner)?, ymax: need(p, "ymax", owner)? };
    if !(r.xmax > r.xmin && r.ymax > r.ymin) {
        return Err(format!("{owner} needs xmin < xmax and ymin < ymax"));
    }
    Ok(r)
}

struct Parser {
    spec: Spec,
    custom: HashMap<String, Material>,
}

impl Parser {
    fn material(&self, name: &str) -> Result<Material, String> {
        self.custom.get(name).copied().or_else(|| by_name(name)).ok_or_else(|| format!("unknown material '{name}'"))
    }

    fn name_is_free(&self, name: &str) -> Result<(), String> {
        if name.contains('=') {
            return Err("a name is missing".into());
        }
        if self.spec.heaters.iter().any(|h| h.name == name) || self.spec.probes.iter().any(|p| p.name == name) {
            return Err(format!("'{name}' is defined twice"));
        }
        Ok(())
    }

    fn statement(&mut self, line: &str, line_number: usize) -> Result<(), String> {
        let w: Vec<&str> = line.split_whitespace().collect();
        let arg = |i: usize| w.get(i).copied().ok_or_else(|| format!("'{}' needs more arguments", w[0]));
        match w[0] {
            "title" => self.spec.title = line[5..].trim().to_string(),
            "grid" => {
                let p = pairs(&w[1..], &["dx", "dy"], "grid")?;
                self.spec.dx = p.get("dx").map(|v| number(v)).transpose()?.unwrap_or(self.spec.dx);
                self.spec.dy = p.get("dy").map(|v| number(v)).transpose()?.unwrap_or(self.spec.dy);
                if !(self.spec.dx > 0.0 && self.spec.dy > 0.0) {
                    return Err("cell sizes must be positive".into());
                }
            }
            "window" => {
                let p = pairs(&w[1..], &["width", "ymin", "ymax"], "window")?;
                self.spec.width = p.get("width").map(|v| number(v)).transpose()?.unwrap_or(self.spec.width);
                self.spec.ymin = p.get("ymin").map(|v| number(v)).transpose()?.unwrap_or(self.spec.ymin);
                self.spec.ymax = p.get("ymax").map(|v| number(v)).transpose()?.unwrap_or(self.spec.ymax);
                if !(self.spec.width > 0.0 && self.spec.ymax > self.spec.ymin) {
                    return Err("the window must have a positive size".into());
                }
            }
            "material" => {
                let name = arg(1)?;
                let p = pairs(&w[2..], &["k", "c"], "material")?;
                let m = Material { conductivity: need(&p, "k", "material")?, capacity: need(&p, "c", "material")? };
                if !(m.conductivity > 0.0 && m.capacity > 0.0) {
                    return Err("conductivity and capacity must be positive".into());
                }
                self.custom.insert(name.to_string(), m);
            }
            "background" => self.spec.background = self.material(arg(1)?)?,
            "layer" => {
                if w.len() != 4 {
                    return Err("usage: layer <material> <ymin> <ymax>".into());
                }
                let r = Rect { xmin: -1e9, xmax: 1e9, ymin: number(w[2])?, ymax: number(w[3])? };
                if !(r.ymax > r.ymin) {
                    return Err("a layer needs ymin < ymax".into());
                }
                self.spec.rects.push((r, self.material(w[1])?));
            }
            "rect" => {
                if w.len() != 6 {
                    return Err("usage: rect <material> <xmin> <xmax> <ymin> <ymax>".into());
                }
                let r = Rect { xmin: number(w[2])?, xmax: number(w[3])?, ymin: number(w[4])?, ymax: number(w[5])? };
                if !(r.xmax > r.xmin && r.ymax > r.ymin) {
                    return Err("a rectangle needs xmin < xmax and ymin < ymax".into());
                }
                self.spec.rects.push((r, self.material(w[1])?));
            }
            "heater" => {
                let name = arg(1)?;
                self.name_is_free(name)?;
                let p = pairs(&w[2..], &["xmin", "xmax", "ymin", "ymax", "power"], "heater")?;
                let power = p.get("power").map(|v| number(v)).transpose()?.unwrap_or(0.1);
                if power < 0.0 {
                    return Err("heater power cannot be negative".into());
                }
                self.spec.heaters.push(Heater { name: name.to_string(), rect: rect_of(&p, "heater")?, power });
            }
            "probe" => {
                let name = arg(1)?;
                self.name_is_free(name)?;
                let p = pairs(&w[2..], &["x", "y", "xmin", "xmax", "ymin", "ymax"], "probe")?;
                let rect = if p.contains_key("x") || p.contains_key("y") {
                    let (x, y) = (need(&p, "x", "probe")?, need(&p, "y", "probe")?);
                    Rect { xmin: x, xmax: x, ymin: y, ymax: y }
                } else {
                    rect_of(&p, "probe")?
                };
                self.spec.probes.push(Probe { name: name.to_string(), rect });
            }
            "top" => {
                self.spec.top = match arg(1)? {
                    "adiabatic" => Top::Adiabatic,
                    "convective" => {
                        let h = need(&pairs(&w[2..], &["h"], "top")?, "h", "top convective")?;
                        if !(h > 0.0) {
                            return Err("the film coefficient must be positive".into());
                        }
                        Top::Convective(h)
                    }
                    _ => return Err("top must be adiabatic or convective".into()),
                }
            }
            "solver" => {
                self.spec.solver = match arg(1)? {
                    "direct" => Solver::Direct,
                    "cg" => Solver::ConjugateGradient,
                    _ => return Err("solver must be direct or cg".into()),
                }
            }
            "optics" => {
                let p = pairs(&w[1..], &["wavelength", "dndt"], "optics")?;
                self.spec.wavelength = p.get("wavelength").map(|v| number(v)).transpose()?.unwrap_or(self.spec.wavelength);
                self.spec.dn_dt = p.get("dndt").map(|v| number(v)).transpose()?.unwrap_or(self.spec.dn_dt);
                if !(self.spec.wavelength > 0.0 && self.spec.dn_dt > 0.0) {
                    return Err("wavelength and dndt must be positive".into());
                }
            }
            "transient" => {
                let mut heater = None;
                let mut rest = Vec::new();
                for word in &w[1..] {
                    match word.strip_prefix("heater=") {
                        Some(name) if !name.is_empty() => heater = Some(name.to_string()),
                        _ => rest.push(*word),
                    }
                }
                let p = pairs(&rest, &["dt", "steps"], "transient")?;
                let (dt, steps) = (need(&p, "dt", "transient")?, need(&p, "steps", "transient")?);
                if !(dt > 0.0 && steps >= 1.0) {
                    return Err("transient needs a positive dt and at least one step".into());
                }
                self.spec.transient = Some(Transient { heater: heater.ok_or("transient needs heater=<name>")?, dt, steps: steps as usize });
            }
            "expect" => self.spec.expectations.push(Expectation { kind: arg(1)?.to_string(), args: w[2..].iter().map(|s| s.to_string()).collect(), line: line_number }),
            other => return Err(format!("unknown statement '{other}'")),
        }
        Ok(())
    }
}

/// Parses a cross-section. Errors carry the line number.
pub fn parse(text: &str) -> Result<Spec, String> {
    let mut parser = Parser { spec: Spec::default(), custom: HashMap::new() };
    for (i, raw) in text.lines().enumerate() {
        let line = raw.split('#').next().unwrap_or("").trim();
        if line.is_empty() {
            continue;
        }
        parser.statement(line, i + 1).map_err(|e| format!("line {}: {e}", i + 1))?;
    }
    let spec = parser.spec;
    if let Some(t) = &spec.transient {
        if !spec.heaters.iter().any(|h| h.name == t.heater) {
            return Err(format!("the transient names a heater that does not exist: '{}'", t.heater));
        }
    }
    Ok(spec)
}
