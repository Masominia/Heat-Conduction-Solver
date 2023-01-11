//! A uniform mesh of rectangular cells, each of one material.

use super::material::Material;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Rect {
    pub xmin: f64,
    pub xmax: f64,
    pub ymin: f64,
    pub ymax: f64,
}

impl Rect {
    pub fn area(&self) -> f64 {
        (self.xmax - self.xmin).max(0.0) * (self.ymax - self.ymin).max(0.0)
    }
}

#[derive(Clone, Debug)]
pub struct Mesh {
    pub nx: usize,
    pub ny: usize,
    pub dx: f64,
    pub dy: f64,
    pub x0: f64,
    pub y0: f64,
    /// Row by row from the bottom: `cells[iy * nx + ix]`.
    pub cells: Vec<Material>,
}

impl Mesh {
    /// A window filled with one material. The sizes must be whole numbers of cells.
    pub fn new(width: f64, ymin: f64, ymax: f64, dx: f64, dy: f64, fill: Material) -> Result<Mesh, String> {
        if !(dx > 0.0 && dy > 0.0) {
            return Err("cell sizes must be positive".into());
        }
        let whole = |length: f64, step: f64, what: &str| -> Result<usize, String> {
            let count = length / step;
            let n = count.round();
            if (count - n).abs() > 1e-6 || n < 2.0 {
                Err(format!("the window is not a whole number of cells in {what}"))
            } else {
                Ok(n as usize)
            }
        };
        let (nx, ny) = (whole(width, dx, "x")?, whole(ymax - ymin, dy, "y")?);
        Ok(Mesh { nx, ny, dx, dy, x0: -0.5 * width, y0: ymin, cells: vec![fill; nx * ny] })
    }

    pub fn len(&self) -> usize {
        self.nx * self.ny
    }

    pub fn is_empty(&self) -> bool {
        self.cells.is_empty()
    }

    pub fn x(&self, ix: usize) -> f64 {
        self.x0 + (ix as f64 + 0.5) * self.dx
    }

    pub fn y(&self, iy: usize) -> f64 {
        self.y0 + (iy as f64 + 0.5) * self.dy
    }

    pub fn index(&self, ix: usize, iy: usize) -> usize {
        iy * self.nx + ix
    }

    /// Gives a material to every cell whose centre lies inside the rectangle.
    pub fn paint(&mut self, rect: &Rect, material: Material) {
        for iy in 0..self.ny {
            let y = self.y(iy);
            if y < rect.ymin || y > rect.ymax {
                continue;
            }
            for ix in 0..self.nx {
                let x = self.x(ix);
                if x >= rect.xmin && x <= rect.xmax {
                    let k = self.index(ix, iy);
                    self.cells[k] = material;
                }
            }
        }
    }

    /// Fraction of each cell covered by the rectangle, as (cell index, covered area in um^2) pairs.
    pub fn coverage(&self, rect: &Rect) -> Vec<(usize, f64)> {
        let mut out = Vec::new();
        for iy in 0..self.ny {
            let (lo, hi) = (self.y0 + iy as f64 * self.dy, self.y0 + (iy + 1) as f64 * self.dy);
            let fy = (hi.min(rect.ymax) - lo.max(rect.ymin)).max(0.0);
            if fy <= 0.0 {
                continue;
            }
            for ix in 0..self.nx {
                let (lo, hi) = (self.x0 + ix as f64 * self.dx, self.x0 + (ix + 1) as f64 * self.dx);
                let fx = (hi.min(rect.xmax) - lo.max(rect.xmin)).max(0.0);
                if fx > 0.0 {
                    out.push((self.index(ix, iy), fx * fy));
                }
            }
        }
        out
    }

    /// The cell containing a point, if it lies inside the window.
    pub fn locate(&self, x: f64, y: f64) -> Option<usize> {
        let (fx, fy) = ((x - self.x0) / self.dx, (y - self.y0) / self.dy);
        if fx < 0.0 || fy < 0.0 || fx >= self.nx as f64 || fy >= self.ny as f64 {
            None
        } else {
            Some(self.index(fx as usize, fy as usize))
        }
    }
}
