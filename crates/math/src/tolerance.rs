/// 误差容限 tol_total = tol_abs + tol_rel * diff.abs()
#[derive(Clone, Copy, Debug)]
pub struct Tolerance {
    /// 绝对误差容限
    pub abs: f64,
    /// 相对误差容限
    pub rel: f64,
}

impl Tolerance {
    pub fn new(abs: f64, rel: f64) -> Self {
        Self { abs, rel }
    }

    pub fn to_f64(&self, raw_value: f64) -> f64 {
        self.abs + self.rel * raw_value
    }

    pub fn approx_eq(&self, a: f64, b: f64) -> bool {
        let diff = (a - b).abs();
        let max_val = a.abs().max(b.abs());
        diff <= self.abs + self.rel * max_val
    }
}
