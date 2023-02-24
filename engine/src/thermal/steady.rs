//! Steady temperature fields.

use super::mesh::{Mesh, Rect};
use super::system::{assemble, Ordering, Top};
use crate::linalg::banded::SymBanded;
use crate::linalg::cg;

/// Source vector (W / m, mesh order) of a heater dissipating `power` mW per um uniformly over a rectangle.
/// The part of the rectangle outside the window is lost; the rest keeps its share.
pub fn source(mesh: &Mesh, rect: &Rect, power: f64) -> Vec<f64> {
    let mut q = vec![0.0; mesh.len()];
    let area = rect.area();
    if area > 0.0 {
        for (cell, covered) in mesh.coverage(rect) {
            q[cell] += 1e3 * power * covered / area;
        }
    }
    q
}

/// A factorised conduction problem that can be solved for many sources.
pub struct Conduction {
    pub order: Ordering,
    pub top: Top,
    matrix: SymBanded,
    factor: SymBanded,
}

impl Conduction {
    pub fn new(mesh: &Mesh, top: Top) -> Result<Conduction, String> {
        let matrix = assemble(mesh, top);
        let mut factor = matrix.clone();
        factor.factorise()?;
        Ok(Conduction { order: Ordering::new(mesh), top, matrix, factor })
    }

    /// Temperature (K above the sink, mesh order) for a source in mesh order, by the direct solver.
    pub fn solve(&self, source: &[f64]) -> Vec<f64> {
        self.order.to_mesh(&self.factor.solve(&self.order.from_mesh(source)))
    }

    /// The same by conjugate gradients; returns the field and the number of iterations.
    pub fn solve_cg(&self, source: &[f64], tolerance: f64) -> Result<(Vec<f64>, usize), String> {
        let result = cg::solve(&self.matrix, &self.order.from_mesh(source), tolerance, 40 * self.matrix.size());
        if !result.converged {
            return Err(format!("conjugate gradients stopped at a residual of {:.1e}", result.residual));
        }
        Ok((self.order.to_mesh(&result.x), result.iterations))
    }

    pub fn matrix(&self) -> &SymBanded {
        &self.matrix
    }
}

/// Area-weighted mean of a field over a rectangle, or the value of the cell under its centre if it covers no area.
pub fn probe(mesh: &Mesh, rect: &Rect, field: &[f64]) -> Option<f64> {
    let covered = mesh.coverage(rect);
    let total: f64 = covered.iter().map(|c| c.1).sum();
    if total > 0.0 {
        Some(covered.iter().map(|(cell, a)| field[*cell] * a).sum::<f64>() / total)
    } else {
        mesh.locate(0.5 * (rect.xmin + rect.xmax), 0.5 * (rect.ymin + rect.ymax)).map(|cell| field[cell])
    }
}
