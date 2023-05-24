use heatphase::thermal::material::{by_name, AIR, SILICA, SILICON};

#[test]
fn names() {
    assert_eq!(by_name("si"), Some(SILICON));
    assert_eq!(by_name("oxide"), Some(SILICA));
    assert_eq!(by_name("air"), Some(AIR));
    assert!(by_name("tin").is_some() && by_name("al").is_some() && by_name("sin").is_some());
    assert_eq!(by_name("wood"), None);
}

#[test]
fn silicon_conducts_about_a_hundred_times_better_than_oxide() {
    let ratio = SILICON.conductivity / SILICA.conductivity;
    assert!(ratio > 90.0 && ratio < 120.0);
}

#[test]
fn diffusivities() {
    // silicon about 90 um^2 per microsecond, oxide just under 1
    assert!((SILICON.diffusivity() - 90.8).abs() < 0.5);
    assert!((SILICA.diffusivity() - 0.857).abs() < 0.01);
    assert!(AIR.diffusivity() > SILICA.diffusivity());
}
