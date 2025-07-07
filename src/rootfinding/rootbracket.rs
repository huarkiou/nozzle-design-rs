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
    pub fn is_converged(&self) -> bool {
        self.converged
    }
}
