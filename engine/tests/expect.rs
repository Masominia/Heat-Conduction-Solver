use heatphase::design::analyse::analyse;
use heatphase::design::expect::{evaluate, Check};
use heatphase::design::parse::parse;

const BASE: &str = "grid dx=0.25 dy=0.25\nwindow width=40 ymin=-5 ymax=5\nlayer si -5 0\nheater h1 xmin=-1 xmax=1 ymin=4 ymax=4.25 power=0.1\n\
                    probe wg xmin=-0.25 xmax=0.25 ymin=2 ymax=2.25\nprobe far x=8.125 y=2.125\ntransient heater=h1 dt=1 steps=200\n";

fn checks(text: &str) -> Vec<Check> {
    let spec = parse(&format!("{BASE}{text}")).unwrap();
    evaluate(&spec, &analyse(&spec).unwrap())
}

#[test]
fn passing_and_failing() {
    let c = checks("expect hotter wg far\nexpect hotter far wg\nexpect balance 1e-9\nexpect cells 6400\nexpect cells 5\n");
    assert_eq!(c.iter().map(|x| x.pass).collect::<Vec<_>>(), vec![true, false, true, true, false]);
    assert_eq!(c[0].name, "hotter wg far");
}

#[test]
fn every_kind_is_understood() {
    let c = checks(
        "expect temperature wg 17 0.5\nexpect response wg h1 170 0.5\nexpect pi_power wg h1 24 0.5\nexpect crosstalk far wg h1 0.1 0.9\nexpect peak 40 0.5\n\
         expect surface_share 0 1e-12\nexpect tau wg 10 0.9\nexpect rise wg 20 0.9\n",
    );
    assert_eq!(c.len(), 8);
    for x in &c {
        assert!(x.pass, "{}: {}", x.name, x.detail);
    }
}

#[test]
fn malformed_expectations_fail_with_their_line() {
    let c = checks("expect peak 40\nexpect wobble 1\nexpect peak abc 0.1\nexpect temperature ghost 1 0.1\nexpect response wg ghost 1 0.1\n");
    assert!(c.iter().all(|x| !x.pass));
    assert!(c[0].detail.contains("line 8") && c[0].detail.contains("takes 2 arguments"));
    assert!(c[1].detail.contains("unknown expectation"));
    assert!(c[2].detail.contains("not a number"));
    assert!(c[3].detail.contains("no probe named 'ghost'"));
    assert!(c[4].detail.contains("no heater named 'ghost'"));
}

#[test]
fn time_expectations_need_a_transient() {
    let spec = parse("window width=40 ymin=-5 ymax=5\nheater h xmin=-1 xmax=1 ymin=4 ymax=4.25\nprobe wg x=0.125 y=2.125\nexpect tau wg 10 0.5\n").unwrap();
    let c = evaluate(&spec, &analyse(&spec).unwrap());
    assert!(!c[0].pass && c[0].detail.contains("needs a transient"));
}
