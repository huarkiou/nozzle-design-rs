pub mod bisection;
mod error;
pub mod newton2d;
pub mod secant;
pub mod toms748;

pub use bisection::{max_iterations, solve_bracket};
pub use error::RootFindingError;
pub use newton2d::newton_2d;

/// 求根结果区间（所有 rootfinding 方法共用）。
///
/// `T = f64` 用于连续域，`T = i64` 用于整数域。
#[derive(Debug, Clone, Copy)]
pub struct RootBracket<T> {
    pub lo: T,
    pub hi: T,
    pub flo: f64,
    pub fhi: f64,
    pub iterations: usize,
    pub converged: bool,
}

impl RootBracket<f64> {
    /// 取区间中点作为根的近似解
    pub fn root(&self) -> f64 {
        self.lo / 2.0 + self.hi / 2.0
    }
}

impl<T> RootBracket<T> {
    /// 区间内是否一定存在根
    pub fn has_root(&self) -> bool {
        self.flo * self.fhi <= 0.0
    }

    /// 是否已收敛
    pub fn is_converged(&self) -> bool {
        self.converged
    }
}

impl RootBracket<i64> {
    /// 取两端点中函数值更接近零的索引
    pub fn best_index(&self) -> i64 {
        if self.flo.abs() < self.fhi.abs() {
            self.lo
        } else {
            self.hi
        }
    }

    /// 取两端点中函数值更接近零的值
    pub fn best_value(&self) -> f64 {
        if self.flo.abs() < self.fhi.abs() {
            self.flo
        } else {
            self.fhi
        }
    }
}
