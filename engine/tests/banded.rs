use heatphase::linalg::banded::SymBanded;

fn sample(n: usize, m: usize) -> SymBanded {
    let mut a = SymBanded::new(n, m);
    for i in 0..n {
        a.add(i, i, 4.0 + 0.1 * i as f64);
        if i + 1 < n {
            a.add(i, i + 1, -1.0);
        }
        if i + m < n && m > 1 {
            a.add(i + m, i, -0.5);
        }
    }
    a
}

#[test]
fn storage_is_symmetric() {
    let a = sample(10, 3);
    assert_eq!(a.get(2, 3), -1.0);
    assert_eq!(a.get(3, 2), -1.0);
    assert_eq!(a.get(5, 2), -0.5);
    assert_eq!(a.get(2, 5), -0.5);
    assert_eq!(a.get(0, 9), 0.0);
    assert_eq!((a.size(), a.half_width()), (10, 3));
}

#[test]
fn multiply_matches_the_definition() {
    let a = sample(8, 3);
    let x: Vec<f64> = (0..8).map(|i| 1.0 + (i * i % 5) as f64).collect();
    let y = a.multiply(&x);
    for i in 0..8 {
        let direct: f64 = (0..8).map(|j| a.get(i, j) * x[j]).sum();
        assert!((y[i] - direct).abs() < 1e-12);
    }
}

#[test]
fn cholesky_solves_the_system() {
    for (n, m) in [(1, 0), (5, 1), (40, 3), (60, 7)] {
        let a = sample(n, m);
        let x: Vec<f64> = (0..n).map(|i| (0.3 * i as f64).sin() + 2.0).collect();
        let b = a.multiply(&x);
        let mut f = a.clone();
        f.factorise().unwrap();
        let got = f.solve(&b);
        let worst = got.iter().zip(&x).map(|(g, w)| (g - w).abs()).fold(0.0, f64::max);
        assert!(worst < 1e-12, "{n} {m}: {worst}");
    }
}

#[test]
fn indefinite_matrices_are_refused() {
    let mut a = SymBanded::new(3, 1);
    a.add(0, 0, 1.0);
    a.add(1, 1, 1.0);
    a.add(2, 2, 1.0);
    a.add(0, 1, 2.0);
    assert!(a.factorise().unwrap_err().contains("not positive definite"));
}

#[test]
fn a_band_wider_than_the_matrix_is_clipped() {
    assert_eq!(SymBanded::new(4, 50).half_width(), 3);
}

#[test]
#[should_panic(expected = "outside the band")]
fn writing_outside_the_band_panics() {
    SymBanded::new(6, 1).add(0, 3, 1.0);
}

#[test]
#[should_panic(expected = "before factorise")]
fn solving_before_factorising_panics() {
    sample(4, 1).solve(&[1.0; 4]);
}
