use num_traits::{Num, ToPrimitive};

use crate::rootfinding::RootFindingError;

/// 二分搜索结果（泛型域）。
///
/// `T = f64` 用于连续域根查找，`T = i64` 用于整数域索引搜索。
#[derive(Debug, Clone, Copy)]
pub struct BisectBracket<T> {
    pub lo: T,
    pub hi: T,
    pub flo: f64,
    pub fhi: f64,
    pub iterations: usize,
    pub converged: bool,
}

impl BisectBracket<f64> {
    /// 取区间中点作为根的近似解
    pub fn root(&self) -> f64 {
        self.lo / 2.0 + self.hi / 2.0
    }
}

impl<T> BisectBracket<T> {
    /// 区间内是否一定存在根
    pub fn has_root(&self) -> bool {
        self.flo * self.fhi <= 0.0
    }

    /// 是否已收敛
    pub fn is_converged(&self) -> bool {
        self.converged
    }
}

impl BisectBracket<i64> {
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

/// 预测达到指定精度所需的最大二分迭代次数
pub fn max_iterations<T>(interval: T, tol: T, ndivide: usize) -> usize
where
    T: Num + ToPrimitive,
{
    assert!(ndivide >= 2);
    let ret = ((interval.to_f64().unwrap() / tol.to_f64().unwrap()).log2()
        / (ndivide as f64).log2())
    .round() as usize
        + 1;
    if ret < 1 {
        1
    } else {
        ret
    }
}

/// 通用二分法查找函数根所在区间。
///
/// 参数化中点计算和收敛判定，同时支持连续域（f64）和整数域（i64）。
///
/// # 参数
/// - `lo`, `hi`: 搜索区间，要求 `f(lo) * f(hi) <= 0`
/// - `f`: 目标函数 `T -> f64`
/// - `midpoint`: 中点计算，如 `|a, b| (a + b) / 2`
/// - `is_done`: 收敛判定，如 `|a, b| b - a <= 1`（整数）或 tolerance 比较（浮点）
/// - `max_iter`: 最大迭代次数
pub fn solve_bracket<T, F, M, C>(
    lo: T,
    hi: T,
    f: &F,
    midpoint: M,
    is_done: C,
    max_iter: usize,
) -> Result<BisectBracket<T>, RootFindingError>
where
    T: Copy,
    F: Fn(T) -> f64,
    M: Fn(T, T) -> T,
    C: Fn(T, T) -> bool,
{
    let mut left = lo;
    let mut right = hi;
    let mut f_left = f(left);
    let mut f_right = f(right);

    if f_left.is_nan() || f_right.is_nan() {
        return Err(RootFindingError::FunctionNotBracketed);
    }
    if f_left * f_right > 0.0 {
        return Err(RootFindingError::FunctionNotBracketed);
    }

    let mut iterations = 0;
    while !is_done(left, right) && iterations < max_iter {
        let mid = midpoint(left, right);
        let f_mid = f(mid);

        if f_mid.is_nan() {
            return Err(RootFindingError::FunctionNotBracketed);
        }

        if f_left * f_mid <= 0.0 {
            right = mid;
            f_right = f_mid;
        } else {
            left = mid;
            f_left = f_mid;
        }
        iterations += 1;
    }

    let converged = iterations < max_iter && is_done(left, right);
    Ok(BisectBracket {
        lo: left,
        hi: right,
        flo: f_left,
        fhi: f_right,
        iterations,
        converged,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Tolerance;

    const MAX_ITER: usize = 100;

    // ── f64 ──

    fn float_bracket(
        f: impl Fn(f64) -> f64,
        a: f64,
        b: f64,
        tol: Tolerance,
        max_iter: usize,
    ) -> Result<BisectBracket<f64>, RootFindingError> {
        // 端点为精确根时直接返回
        let fa = f(a);
        if fa == 0.0 {
            return Ok(BisectBracket {
                lo: a,
                hi: a,
                flo: fa,
                fhi: fa,
                iterations: 0,
                converged: true,
            });
        }
        let fb = f(b);
        if fb == 0.0 {
            return Ok(BisectBracket {
                lo: b,
                hi: b,
                flo: fb,
                fhi: fb,
                iterations: 0,
                converged: true,
            });
        }
        solve_bracket(
            a,
            b,
            &f,
            |x, y| (x + y) / 2.0,
            |x, y| tol.approx_eq(x, y),
            max_iter,
        )
    }

    #[test]
    fn test_max_iterations_predict() {
        assert_eq!(max_iterations(10 - 1, 1, 2), 4);
        assert_eq!(max_iterations(9.2, 1.01, 2), 4);
    }

    #[test]
    fn test_float_cos_x() {
        let f = |x: f64| x.cos() - x;
        let tol = Tolerance::new(1e-10, 1e-10);
        let r = float_bracket(f, 0.0, 1.0, tol, MAX_ITER).unwrap();
        assert!(r.converged);
        assert!((r.hi - r.lo) < 2.0 * tol.to_f64(1.0));
        assert!(r.has_root());
    }

    #[test]
    fn test_float_not_bracketed() {
        let f = |x: f64| x * x + 1.0;
        let r = float_bracket(f, 0.0, 2.0, Tolerance::new(1e-8, 1e-8), 100);
        assert!(matches!(r, Err(RootFindingError::FunctionNotBracketed)));
    }

    #[test]
    fn test_float_convergence_failure() {
        let f = |x: f64| x.powi(3) - 2.0 * x - 5.0;
        let tol = Tolerance::new(1e-12, 1e-12);
        let r = float_bracket(f, 2.0, 3.0, tol, 3).unwrap();
        assert!(!r.converged);
        assert!(r.has_root());
    }

    #[test]
    fn test_float_oscillatory() {
        let f = |x: f64| (1.0 / x).sin();
        let tol = Tolerance::new(1e-10, 1e-10);
        let r = float_bracket(f, 0.1, 1.0, tol, MAX_ITER).unwrap();
        assert!(r.converged);
        assert!(r.has_root());
    }

    #[test]
    fn test_float_steep() {
        let f = |x: f64| (x - 1.0).powi(2) * x.log10();
        let tol = Tolerance::new(1e-10, 1e-10);
        let r = float_bracket(f, 0.5, 2.0, tol, MAX_ITER).unwrap();
        assert!(r.converged);
        assert!(r.has_root());
    }

    #[test]
    fn test_float_high_precision() {
        let f = |x: f64| x.cos() - x;
        let tol = Tolerance::new(1e-15, 1e-15);
        let r = float_bracket(f, 0.0, 1.0, tol, MAX_ITER).unwrap();
        assert!(r.converged);
    }

    #[test]
    fn test_float_low_max_iter() {
        let f = |x: f64| x.cos() - x;
        let tol = Tolerance::new(1e-10, 1e-10);
        let r = float_bracket(f, 0.0, 1.0, tol, 3).unwrap();
        assert!(!r.converged);
        assert_eq!(r.iterations, 3);
        assert!(r.has_root());
    }

    #[test]
    fn test_float_multiple_roots() {
        let f = |x: f64| (x - 1.0) * (x - 2.0) * (x - 3.0);
        let tol = Tolerance::new(1e-10, 1e-10);
        let r = float_bracket(f, 0.5, 3.1, tol, MAX_ITER).unwrap();
        assert!(r.converged);
        assert!(r.has_root());
    }

    #[test]
    fn test_float_discontinuous() {
        let f = |x: f64| if x < 0.0 { -1.0 } else { 1.0 };
        let tol = Tolerance::new(1e-10, 1e-10);
        let r = float_bracket(f, -1.0, 1.0, tol, MAX_ITER).unwrap();
        assert!(r.converged);
        assert!(r.has_root());
    }

    // ── i64 ──

    fn int_bracket(
        f: impl Fn(i64) -> f64,
        lo: i64,
        hi: i64,
        max_iter: usize,
    ) -> Result<BisectBracket<i64>, RootFindingError> {
        if lo >= hi {
            return Err(RootFindingError::InvalidInterval);
        }
        solve_bracket(lo, hi, &f, |a, b| (a + b) / 2, |a, b| b - a <= 1, max_iter)
    }

    #[test]
    fn test_int_basic() {
        let f = |i: i64| (i - 5) as f64;
        let r = int_bracket(f, 0, 10, 20).unwrap();
        assert!(r.has_root());
        assert!(r.hi - r.lo <= 1);
        assert!(r.lo <= 5 && 5 <= r.hi);
    }

    #[test]
    fn test_int_negative() {
        let f = |i: i64| (i + 3) as f64;
        let r = int_bracket(f, -10, 0, 20).unwrap();
        assert!(r.lo <= -3 && -3 <= r.hi);
    }

    #[test]
    fn test_int_large() {
        let f = |i: i64| (i - 500) as f64;
        let r = int_bracket(f, 0, 1000, 20).unwrap();
        assert!(r.lo <= 500 && 500 <= r.hi);
    }

    #[test]
    fn test_int_best_index() {
        let f = |i: i64| (i - 7) as f64;
        let r = int_bracket(f, 0, 10, 20).unwrap();
        assert_eq!(r.best_index(), 7);
    }

    #[test]
    fn test_int_not_bracketed() {
        let f = |i: i64| (i * i) as f64;
        let r = int_bracket(f, -5, 5, 20);
        assert!(r.is_err());
    }

    #[test]
    fn test_int_invalid() {
        let f = |i: i64| i as f64;
        let r = int_bracket(f, 5, 3, 20);
        assert!(r.is_err());
    }
}
