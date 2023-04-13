//! Closed forms the solver is checked against.

use std::f64::consts::PI;

/// Temperature at the top of a stack of layers (thickness in um, conductivity in W/(m K)), heated uniformly
/// over its whole width with `power` mW per um of length spread over `width` um, bottom at the sink: the
/// one-dimensional result q'' times the sum of t / k.
pub fn stack_temperature(layers: &[(f64, f64)], power: f64, width: f64) -> f64 {
    let flux = 1e3 * power / width; // W/m per um of width
    flux * layers.iter().map(|(t, k)| t / k).sum::<f64>()
}

/// Surface temperature at a distance x from a line heater on top of a uniform slab of thickness d whose bottom
/// is at the sink and whose top is adiabatic: P' / (pi k) ln coth(pi |x| / (4 d)).
/// `power` in mW per um, lengths in um, k in W/(m K).
pub fn line_source_surface(power: f64, conductivity: f64, thickness: f64, x: f64) -> f64 {
    let u = PI * x.abs() / (4.0 * thickness);
    1e3 * power / (PI * conductivity) * (1.0 / u.tanh()).ln()
}

/// Mean surface temperature under a strip heater of the given width on the same slab, by integrating the line
/// source over the strip and averaging over it (the logarithmic singularity is integrable).
pub fn strip_surface_mean(power: f64, conductivity: f64, thickness: f64, width: f64) -> f64 {
    // mean over source and observation points of the line-source kernel: weight (1 - s / w) on separation s
    let n = 20000;
    let mut sum = 0.0;
    for i in 0..n {
        let s = (i as f64 + 0.5) / n as f64 * width;
        sum += (1.0 - s / width) * line_source_surface(1.0, conductivity, thickness, s);
    }
    power * 2.0 * sum / n as f64
}

/// Decay length of the surface temperature far from a heater on that slab: T ~ exp(-x / L), L = 2 d / pi.
pub fn crosstalk_decay_length(thickness: f64) -> f64 {
    2.0 * thickness / PI
}

/// Slowest time constant of a uniform slab with an adiabatic top and the bottom at the sink:
/// 4 d^2 / (pi^2 alpha), in microseconds for d in um and the diffusivity in um^2 per microsecond.
pub fn slab_time_constant(thickness: f64, diffusivity: f64) -> f64 {
    4.0 * thickness * thickness / (PI * PI * diffusivity)
}

/// Heater power for a phase shift of pi, mW: lambda / (2 dn/dT R), with R the temperature rise of the
/// waveguide in K per (mW per um). Independent of the heater length.
pub fn pi_power(wavelength: f64, dn_dt: f64, response: f64) -> f64 {
    wavelength / (2.0 * dn_dt * response)
}
