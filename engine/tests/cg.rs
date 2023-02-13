use heatphase::linalg::banded::SymBanded;
use heatphase::linalg::cg;

fn laplacian(n: usize) -> SymBanded {
    let mut a = SymBanded::new(n, 1);
    for i in 0..n {
        a.add(i, i, 2.0 + 1e-3 * i as f64);
        if i + 1 < n {
            a.add(i, i + 1, -1.0);
        }
    }
    a
}

#[test]
fn agrees_with_the_direct_solver() {
    let a = laplacian(200);
    let b: Vec<f64> = (0..200).map(|i| ((i * 7) % 11) as f64 - 3.0).collect();
    let mut f = a.clone();
    f.factorise().unwrap();
    let direct = f.solve(&b);
    let r = cg::solve(&a, &b, 1e-13, 5000);
    assert!(r.converged && r.residual <= 1e-13);
    let worst = r.x.iter().zip(&direct).map(|(p, q)| (p - q).abs()).fold(0.0, f64::max);
    let scale = direct.iter().fold(0.0f64, |m, v| m.max(v.abs()));
    assert!(worst < 1e-9 * scale, "{worst}");
}

#[test]
fn converges_in_at_most_n_iterations_on_a_small_system() {
    let a = laplacian(12);
    let r = cg::solve(&a, &[1.0; 12], 1e-12, 100);
    assert!(r.converged && r.iterations <= 13, "{}", r.iterations);
}

#[test]
fn a_zero_right_hand_side_needs_no_iteration() {
    let r = cg::solve(&laplacian(5), &[0.0; 5], 1e-12, 10);
    assert!(r.converged && r.iterations == 0 && r.x.iter().all(|v| *v == 0.0));
}

#[test]
fn too_few_iterations_are_reported() {
    let r = cg::solve(&laplacian(300), &[1.0; 300], 1e-14, 3);
    assert!(!r.converged && r.iterations == 3 && r.residual > 1e-14);
}
