//! Symmetric positive definite banded matrices: storage of the lower band, Cholesky factorisation, solves.

#[derive(Clone, Debug)]
pub struct SymBanded {
    n: usize,
    m: usize,
    /// Row i holds elements (i, i - m) ... (i, i): `data[i * (m + 1) + (j + m - i)]`.
    data: Vec<f64>,
    factorised: bool,
}

impl SymBanded {
    pub fn new(n: usize, half_width: usize) -> SymBanded {
        let m = half_width.min(n.saturating_sub(1));
        SymBanded { n, m, data: vec![0.0; n * (m + 1)], factorised: false }
    }

    pub fn size(&self) -> usize {
        self.n
    }

    pub fn half_width(&self) -> usize {
        self.m
    }

    fn at(&self, i: usize, j: usize) -> usize {
        i * (self.m + 1) + (j + self.m - i)
    }

    /// Element (i, j) for |i - j| within the band; either triangle may be asked for.
    pub fn get(&self, i: usize, j: usize) -> f64 {
        let (i, j) = if i >= j { (i, j) } else { (j, i) };
        if i - j > self.m {
            0.0
        } else {
            self.data[self.at(i, j)]
        }
    }

    /// Adds to element (i, j) and, by symmetry, to (j, i). Panics outside the band.
    pub fn add(&mut self, i: usize, j: usize, value: f64) {
        let (i, j) = if i >= j { (i, j) } else { (j, i) };
        assert!(i - j <= self.m, "element outside the band");
        let k = self.at(i, j);
        self.data[k] += value;
    }

    pub fn multiply(&self, x: &[f64]) -> Vec<f64> {
        let mut y = vec![0.0; self.n];
        for i in 0..self.n {
            let lo = i.saturating_sub(self.m);
            for j in lo..=i {
                let a = self.data[self.at(i, j)];
                y[i] += a * x[j];
                if j != i {
                    y[j] += a * x[i];
                }
            }
        }
        y
    }

    /// Replaces the matrix by its Cholesky factor L (A = L L^T). Fails if the matrix is not positive definite.
    pub fn factorise(&mut self) -> Result<(), String> {
        let (n, m, w) = (self.n, self.m, self.m + 1);
        for i in 0..n {
            let lo = i.saturating_sub(m);
            for j in lo..=i {
                // sum over k of L(i, k) L(j, k), k from max(lo, j - m) to j - 1
                let start = lo.max(j.saturating_sub(m));
                let mut sum = self.data[i * w + (j + m - i)];
                for k in start..j {
                    sum -= self.data[i * w + (k + m - i)] * self.data[j * w + (k + m - j)];
                }
                if i == j {
                    if sum <= 0.0 {
                        return Err(format!("the matrix is not positive definite (row {i})"));
                    }
                    self.data[i * w + m] = sum.sqrt();
                } else {
                    self.data[i * w + (j + m - i)] = sum / self.data[j * w + m];
                }
            }
        }
        self.factorised = true;
        Ok(())
    }

    /// Solves A x = b with the factor from `factorise`.
    pub fn solve(&self, b: &[f64]) -> Vec<f64> {
        assert!(self.factorised, "solve called before factorise");
        let (n, m, w) = (self.n, self.m, self.m + 1);
        let mut x = b.to_vec();
        for i in 0..n {
            let lo = i.saturating_sub(m);
            let mut s = x[i];
            for k in lo..i {
                s -= self.data[i * w + (k + m - i)] * x[k];
            }
            x[i] = s / self.data[i * w + m];
        }
        for i in (0..n).rev() {
            let hi = (i + m).min(n - 1);
            let mut s = x[i];
            for k in i + 1..=hi {
                s -= self.data[k * w + (i + m - k)] * x[k];
            }
            x[i] = s / self.data[i * w + m];
        }
        x
    }
}
