//! Steady and transient heat conduction in the cross-section of a photonic chip.
//!
//! Lengths are in micrometres, temperatures in kelvin above the heat sink, heater powers in milliwatts per
//! micrometre of heater length, times in microseconds.

pub mod analytic;
pub mod design;
pub mod linalg;
pub mod thermal;

pub use design::spec::Spec;
