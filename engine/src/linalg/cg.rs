//! Conjugate gradients with a diagonal preconditioner, as an independent way to solve the same system.

use super::banded::SymBanded;

pub struct CgResult {
    pub x: Vec<f64>,
    pub iterations: usize,
    /// Norm of the residual over the norm of the right-hand side.
    pub residual: f64,
    pub converged: bool,
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// Solves A x = b for a symmetric positive definite band matrix that has not been factorised.
pub fn solve(a: &SymBanded, b: &[f64], tolerance: f64, max_iterations: usize) -> CgResult {
    let n = a.size();
    let inverse_diagonal: Vec<f64> = (0..n).map(|i| 1.0 / a.get(i, i)).collect();
    let norm_b = dot(b, b).sqrt();
    let mut x = vec![0.0; n];
    if norm_b == 0.0 {
        return CgResult { x, iterations: 0, residual: 0.0, converged: true };
    }
    let mut r = b.to_vec();
    let mut z: Vec<f64> = r.iter().zip(&inverse_diagonal).map(|(r, d)| r * d).collect();
    let mut p = z.clone();
    let mut rz = dot(&r, &z);
    let mut residual = 1.0;
    for iteration in 1..=max_iterations {
        let ap = a.multiply(&p);
        let alpha = rz / dot(&p, &ap);
        for i in 0..n {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }
        residual = dot(&r, &r).sqrt() / norm_b;
        if residual <= tolerance {
            return CgResult { x, iterations: iteration, residual, converged: true };
        }
        for i in 0..n {
            z[i] = r[i] * inverse_diagonal[i];
        }
        let rz_new = dot(&r, &z);
        let beta = rz_new / rz;
        rz = rz_new;
        for i in 0..n {
            p[i] = z[i] + beta * p[i];
        }
    }
    CgResult { x, iterations: max_iterations, residual, converged: false }
}
