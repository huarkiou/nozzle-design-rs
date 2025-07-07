#[derive(Debug, Clone, Copy)]
pub struct RootBracket {
    pub a: f64,
    pub b: f64,
    pub fa: f64,
    pub fb: f64,
    pub iterations: usize,
    pub converged: bool,
}

impl RootBracket {
    pub fn bracket(&self) -> (f64, f64) {
        (self.a, self.b)
    }

    pub fn root(&self) -> f64 {
        self.a / 2. + self.b / 2.
    }

    pub fn has_root(&self) -> bool {
        self.fa * self.fb <= 0.0
    }

    pub fn is_converged(&self) -> bool {
        self.converged
    }
}
