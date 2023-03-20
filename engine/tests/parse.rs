use heatphase::design::parse::parse;
use heatphase::design::spec::Solver;
use heatphase::thermal::material::{SILICA, SILICON};
use heatphase::thermal::system::Top;

#[test]
fn defaults() {
    let s = parse("").unwrap();
    assert_eq!((s.dx, s.dy, s.width, s.ymin, s.ymax), (0.25, 0.25, 60.0, -20.0, 5.0));
    assert_eq!(s.background, SILICA);
    assert_eq!((s.top, s.solver), (Top::Adiabatic, Solver::Direct));
    assert!(s.heaters.is_empty() && s.transient.is_none());
}

#[test]
fn a_full_file() {
    let text = "# note\ntitle My shifter\ngrid dx=0.5 dy=0.125\nwindow width=40 ymin=-5 ymax=6\nmaterial polymer k=0.2 c=1.8e6\nbackground polymer\n\
                layer si -5 0\nlayer sio2 0 3\nrect si -0.25 0.25 2 2.25 # core\nheater h1 xmin=-1 xmax=1 ymin=4 ymax=4.25 power=0.2\nheater h2 xmin=4 xmax=6 ymin=4 ymax=4.25\n\
                probe wg xmin=-0.25 xmax=0.25 ymin=2 ymax=2.25\nprobe point x=3 y=1\ntop convective h=50\nsolver cg\noptics wavelength=1.31 dndt=1.9e-4\n\
                transient heater=h2 dt=0.5 steps=100\nexpect peak 30 0.5\n";
    let s = parse(text).unwrap();
    assert_eq!(s.title, "My shifter");
    assert_eq!((s.dx, s.dy, s.width, s.ymin, s.ymax), (0.5, 0.125, 40.0, -5.0, 6.0));
    assert_eq!(s.background.conductivity, 0.2);
    assert_eq!(s.rects.len(), 3);
    assert_eq!(s.rects[0].1, SILICON);
    assert!(s.rects[0].0.xmin < -1e8);
    assert_eq!((s.heaters[0].power, s.heaters[1].power), (0.2, 0.1));
    assert_eq!(s.heaters[1].rect.xmax, 6.0);
    assert_eq!(s.probes[1].rect.xmin, s.probes[1].rect.xmax);
    assert_eq!((s.top, s.solver), (Top::Convective(50.0), Solver::ConjugateGradient));
    assert_eq!((s.wavelength, s.dn_dt), (1.31, 1.9e-4));
    let t = s.transient.unwrap();
    assert_eq!((t.heater.as_str(), t.dt, t.steps), ("h2", 0.5, 100));
    assert_eq!((s.expectations[0].kind.as_str(), s.expectations[0].line), ("peak", 18));
}

#[test]
fn errors_name_the_line() {
    for (text, needle) in [
        ("grid dx=0.5\nfrobnicate", "line 2: unknown statement"),
        ("grid dx=0", "cell sizes"),
        ("grid dz=1", "does not have 'dz'"),
        ("window width=10 ymin=2 ymax=1", "positive size"),
        ("material glass k=1", "needs c="),
        ("material glass k=-1 c=1", "must be positive"),
        ("background wood", "unknown material 'wood'"),
        ("layer si 2 1", "ymin < ymax"),
        ("layer si 1", "usage: layer"),
        ("rect si 1 0 0 1", "xmin < xmax"),
        ("rect si a 1 0 1", "not a number"),
        ("heater h xmin=0 xmax=1 ymin=0", "needs ymax"),
        ("heater h xmin=0 xmax=1 ymin=0 ymax=1 power=-1", "cannot be negative"),
        ("heater xmin=0 xmax=1 ymin=0 ymax=1", "name is missing"),
        ("heater h xmin=0 xmax=1 ymin=0 ymax=1\nprobe h x=0 y=0", "defined twice"),
        ("probe p x=1", "needs y"),
        ("top open", "adiabatic or convective"),
        ("top convective", "needs h"),
        ("top convective h=0", "must be positive"),
        ("solver magic", "direct or cg"),
        ("optics wavelength=0", "must be positive"),
        ("transient dt=1 steps=10", "needs heater"),
        ("heater h xmin=0 xmax=1 ymin=0 ymax=1\ntransient heater=h dt=0 steps=10", "positive dt"),
        ("transient heater=ghost dt=1 steps=10", "does not exist"),
        ("expect", "needs more arguments"),
    ] {
        let err = parse(text).unwrap_err();
        assert!(err.contains(needle), "{text:?} gave {err:?}");
    }
}
