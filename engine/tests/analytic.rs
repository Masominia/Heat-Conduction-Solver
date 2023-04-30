use heatphase::analytic::{crosstalk_decay_length, line_source_surface, pi_power, slab_time_constant, stack_temperature, strip_surface_mean};

#[test]
fn stack_is_flux_times_resistance() {
    // 100 W/m over 10 um is 10 W/m per um; 2 um of k = 2 is 1 um K / W ... 10 K
    assert!((stack_temperature(&[(2.0, 2.0)], 0.1, 10.0) - 10.0).abs() < 1e-12);
    assert!((stack_temperature(&[(2.0, 2.0), (3.0, 1.5)], 0.1, 10.0) - 30.0).abs() < 1e-12);
}

#[test]
fn line_source_limits() {
    let (p, k, d) = (0.1, 1.38, 2.0);
    // near the source: logarithmic, P / (pi k) ln(4 d / (pi x))
    let x = 0.01;
    let near = 1e3 * p / (std::f64::consts::PI * k) * (4.0 * d / (std::f64::consts::PI * x)).ln();
    assert!((line_source_surface(p, k, d, x) / near - 1.0).abs() < 1e-4);
    // far from it: 2 P / (pi k) exp(-pi x / (2 d))
    let x = 12.0;
    let far = 2e3 * p / (std::f64::consts::PI * k) * (-x / crosstalk_decay_length(d)).exp();
    assert!((line_source_surface(p, k, d, x) / far - 1.0).abs() < 1e-3);
    // even in x, proportional to power
    assert_eq!(line_source_surface(p, k, d, -3.0), line_source_surface(p, k, d, 3.0));
    assert!((line_source_surface(2.0 * p, k, d, 3.0) / line_source_surface(p, k, d, 3.0) - 2.0).abs() < 1e-12);
}

#[test]
fn a_wide_strip_tends_to_the_one_dimensional_result() {
    let (k, d, w) = (1.38, 0.5, 40.0);
    let strip = strip_surface_mean(0.1, k, d, w);
    let flat = stack_temperature(&[(d, k)], 0.1, w);
    assert!((strip / flat - 1.0).abs() < 0.03, "{strip} {flat}");
    // a narrow strip is hotter than the same power spread wide
    assert!(strip_surface_mean(0.1, k, 2.0, 1.0) > strip_surface_mean(0.1, k, 2.0, 4.0));
}

#[test]
fn time_constant_and_pi_power() {
    assert!((slab_time_constant(2.0, 1.0) - 16.0 / std::f64::consts::PI.powi(2)).abs() < 1e-12);
    assert!((slab_time_constant(4.0, 1.0) / slab_time_constant(2.0, 1.0) - 4.0).abs() < 1e-12);
    // 1.55 um, 1.86e-4 per K, 170 K per mW/um: about 24.5 mW
    assert!((pi_power(1.55, 1.86e-4, 170.0) - 24.51).abs() < 0.01);
}
