//! The conduction matrix of a mesh: finite volumes with harmonic-mean conductances between cells.

use super::mesh::Mesh;
use crate::linalg::banded::SymBanded;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Top {
    /// No heat leaves through the top surface.
    Adiabatic,
    /// Film coefficient in W / (m^2 K) to surroundings at the sink temperature.
    Convective(f64),
}

/// Cells are numbered along the shorter side of the mesh so that the band is as narrow as possible.
#[derive(Clone, Copy, Debug)]
pub struct Ordering {
    pub nx: usize,
    pub ny: usize,
    pub by_rows: bool,
}

impl Ordering {
    pub fn new(mesh: &Mesh) -> Ordering {
        Ordering { nx: mesh.nx, ny: mesh.ny, by_rows: mesh.nx <= mesh.ny }
    }

    pub fn index(&self, ix: usize, iy: usize) -> usize {
        if self.by_rows {
            iy * self.nx + ix
        } else {
            ix * self.ny + iy
        }
    }

    pub fn half_width(&self) -> usize {
        if self.by_rows {
            self.nx
        } else {
            self.ny
        }
    }

    /// Rearranges a vector indexed in solver order into mesh order (row by row).
    pub fn to_mesh(&self, v: &[f64]) -> Vec<f64> {
        let mut out = vec![0.0; v.len()];
        for iy in 0..self.ny {
            for ix in 0..self.nx {
                out[iy * self.nx + ix] = v[self.index(ix, iy)];
            }
        }
        out
    }

    pub fn from_mesh(&self, v: &[f64]) -> Vec<f64> {
        let mut out = vec![0.0; v.len()];
        for iy in 0..self.ny {
            for ix in 0..self.nx {
                out[self.index(ix, iy)] = v[iy * self.nx + ix];
            }
        }
        out
    }
}

/// Conductance per unit length (W / (m K)) from a boundary cell to the sink below it, and to the surroundings
/// above the top surface.
pub fn bottom_conductance(mesh: &Mesh, ix: usize) -> f64 {
    mesh.cells[mesh.index(ix, 0)].conductivity * mesh.dx / (0.5 * mesh.dy)
}

pub fn top_conductance(mesh: &Mesh, ix: usize, top: Top) -> f64 {
    match top {
        Top::Adiabatic => 0.0,
        Top::Convective(h) => {
            let k = mesh.cells[mesh.index(ix, mesh.ny - 1)].conductivity;
            // half a cell of conduction in series with the film; dx in metres for the film
            1.0 / (1.0 / (h * mesh.dx * 1e-6) + 0.5 * mesh.dy / (k * mesh.dx))
        }
    }
}

/// The matrix A of A T = q, with T in kelvin above the sink and q in watts per metre of length.
/// The bottom is held at the sink temperature, the sides are adiabatic.
pub fn assemble(mesh: &Mesh, top: Top) -> SymBanded {
    let order = Ordering::new(mesh);
    let mut a = SymBanded::new(mesh.len(), order.half_width());
    let harmonic = |p: f64, q: f64| 2.0 * p * q / (p + q);
    for iy in 0..mesh.ny {
        for ix in 0..mesh.nx {
            let me = order.index(ix, iy);
            let k = mesh.cells[mesh.index(ix, iy)].conductivity;
            if ix + 1 < mesh.nx {
                let g = harmonic(k, mesh.cells[mesh.index(ix + 1, iy)].conductivity) * mesh.dy / mesh.dx;
                let other = order.index(ix + 1, iy);
                a.add(me, me, g);
                a.add(other, other, g);
                a.add(me, other, -g);
            }
            if iy + 1 < mesh.ny {
                let g = harmonic(k, mesh.cells[mesh.index(ix, iy + 1)].conductivity) * mesh.dx / mesh.dy;
                let other = order.index(ix, iy + 1);
                a.add(me, me, g);
                a.add(other, other, g);
                a.add(me, other, -g);
            }
            if iy == 0 {
                a.add(me, me, bottom_conductance(mesh, ix));
            }
            if iy + 1 == mesh.ny {
                a.add(me, me, top_conductance(mesh, ix, top));
            }
        }
    }
    a
}

/// Heat capacity of each cell per metre of length, J / (m K), in mesh order.
pub fn capacities(mesh: &Mesh) -> Vec<f64> {
    mesh.cells.iter().map(|m| m.capacity * mesh.dx * mesh.dy * 1e-12).collect()
}

/// Heat leaving through the sink and through the top surface, W / m, for a temperature field in mesh order.
pub fn outflow(mesh: &Mesh, top: Top, temperature: &[f64]) -> (f64, f64) {
    let (mut sink, mut surface) = (0.0, 0.0);
    for ix in 0..mesh.nx {
        sink += bottom_conductance(mesh, ix) * temperature[mesh.index(ix, 0)];
        surface += top_conductance(mesh, ix, top) * temperature[mesh.index(ix, mesh.ny - 1)];
    }
    (sink, surface)
}
