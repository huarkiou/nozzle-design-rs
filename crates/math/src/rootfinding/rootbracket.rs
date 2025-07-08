/// 一维区间求根结果
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
    /// 获取根所在的区间(a,b)
    pub fn bracket(&self) -> (f64, f64) {
        (self.a, self.b)
    }

    /// 获取根的数值解(a+b)/2
    pub fn root(&self) -> f64 {
        self.a / 2. + self.b / 2.
    }

    /// 区间(a,b)内是否一定存在根
    pub fn has_root(&self) -> bool {
        self.fa * self.fb <= 0.0
    }

    /// 计算是否收敛
    pub fn is_converged(&self) -> bool {
        self.converged
    }
}
