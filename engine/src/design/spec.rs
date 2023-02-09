//! What a cross-section file says.

use crate::thermal::material::{Material, SILICA};
use crate::thermal::mesh::Rect;
use crate::thermal::system::Top;

#[derive(Clone, Debug, PartialEq)]
pub struct Heater {
    pub name: String,
    pub rect: Rect,
    /// Dissipated power in mW per um of heater length.
    pub power: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Probe {
    pub name: String,
    /// A rectangle to average over; a point probe has zero size.
    pub rect: Rect,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Solver {
    Direct,
    ConjugateGradient,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Transient {
    pub heater: String,
    /// Time step in microseconds.
    pub dt: f64,
    pub steps: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Expectation {
    pub kind: String,
    pub args: Vec<String>,
    pub line: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Spec {
    pub title: String,
    pub dx: f64,
    pub dy: f64,
    /// Full width of the window, centred on x = 0.
    pub width: f64,
    pub ymin: f64,
    pub ymax: f64,
    pub background: Material,
    /// Painted in order, later ones on top.
    pub rects: Vec<(Rect, Material)>,
    pub heaters: Vec<Heater>,
    pub probes: Vec<Probe>,
    pub top: Top,
    pub solver: Solver,
    pub wavelength: f64,
    /// Thermo-optic coefficient of the waveguide core, 1 / K.
    pub dn_dt: f64,
    pub transient: Option<Transient>,
    pub expectations: Vec<Expectation>,
}

impl Default for Spec {
    fn default() -> Self {
        Spec {
            title: String::new(),
            dx: 0.25,
            dy: 0.25,
            width: 60.0,
            ymin: -20.0,
            ymax: 5.0,
            background: SILICA,
            rects: Vec::new(),
            heaters: Vec::new(),
            probes: Vec::new(),
            top: Top::Adiabatic,
            solver: Solver::Direct,
            wavelength: 1.55,
            dn_dt: 1.86e-4,
            transient: None,
            expectations: Vec::new(),
        }
    }
}
