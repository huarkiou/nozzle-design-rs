use num_traits::{Num, ToPrimitive};

use crate::{
    rootfinding::{RootBracket, RootFindingError},
    Tolerance,
};

// ── 整数域二分搜索 ────────────────────────────────────────────

/// 整数域二分搜索结果
#[derive(Debug, Clone, Copy)]
pub struct IntBracket {
    /// 区间下界
    pub lo: i64,
    /// 区间上界
    pub hi: i64,
    /// f(lo)
    pub flo: f64,
    /// f(hi)
    pub fhi: f64,
    /// 迭代次数
    pub iterations: usize,
}

impl IntBracket {
    /// 区间内是否一定存在根 (flo * fhi <= 0)
    pub fn has_root(&self) -> bool {
        self.flo * self.fhi <= 0.0
    }

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

/// 整数域二分法查找函数根所在区间。
///
/// 在闭区间 `[lo, hi]` 上搜索使 `f` 变号的最小子区间 `[left, right]`，
/// 满足 `right - left <= 1`。
///
/// # 参数
/// - `lo`, `hi`: 整数搜索区间，要求 `f(lo) * f(hi) <= 0`
/// - `f`: 目标函数 `i64 -> f64`
/// - `max_iter`: 最大迭代次数
///
/// # 返回
/// - `Ok(IntBracket)`: 包含根的相邻整数区间
/// - `Err(RootFindingError)`: 区间无效、未包围根、或迭代中出现 NaN
pub fn solve_bracket_int<F>(
    lo: i64,
    hi: i64,
    f: &F,
    max_iter: usize,
) -> Result<IntBracket, RootFindingError>
where
    F: Fn(i64) -> f64,
{
    if lo >= hi {
        return Err(RootFindingError::InvalidInterval);
    }

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
    while right - left > 1 && iterations < max_iter {
        let mid = (left + right) / 2;
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

    Ok(IntBracket {
        lo: left,
        hi: right,
        flo: f_left,
        fhi: f_right,
        iterations,
    })
}

// ── 连续域二分搜索 ────────────────────────────────────────────

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

/// 二分搜索法查找函数根
///
/// # 参数
/// - `a`, `b`: 初始区间，要求函数f在区间端点的函数值符号相反
/// - `f`: 连续函数
/// - `tol`: 误差容限
/// - `max_iter`: 最大迭代次数
pub fn solve_bracket<F>(
    a: f64,
    b: f64,
    f: &F,
    tol: Tolerance,
    max_iter: usize,
) -> Result<RootBracket, RootFindingError>
where
    F: Fn(f64) -> f64,
{
    if a >= b {
        return Err(RootFindingError::InvalidInterval);
    }

    let mut left = a;
    let mut right = b;
    let mut f_left = f(left);
    let mut f_right = f(right);

    // 检查是否包围了根
    if f_left * f_right > 0.0 {
        return Err(RootFindingError::FunctionNotBracketed);
    }

    // 检查端点是否本身就是根
    if f_left == 0.0 {
        return Ok(RootBracket {
            a: left,
            b: left,
            fa: f_left,
            fb: f_left,
            iterations: 0,
            converged: true,
        });
    }
    if f_right == 0.0 {
        return Ok(RootBracket {
            a: right,
            b: right,
            fa: f_right,
            fb: f_right,
            iterations: 0,
            converged: true,
        });
    }

    let mut iter = 0;
    while !tol.approx_eq(left, right) && iter < max_iter {
        let mid = (left + right) / 2.0;
        let f_mid = f(mid);

        if f_mid == 0.0 {
            return Ok(RootBracket {
                a: mid,
                b: mid,
                fa: f_mid,
                fb: f_mid,
                iterations: iter + 1,
                converged: true,
            });
        }

        if f_left * f_mid < 0.0 {
            right = mid;
            f_right = f_mid;
        } else {
            left = mid;
            f_left = f_mid;
        }

        iter += 1;
    }

    let converged = iter < max_iter;

    Ok(RootBracket {
        a: left,
        b: right,
        fa: f_left,
        fb: f_right,
        iterations: iter,
        converged,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: Tolerance = Tolerance {
        abs: 1e-10,
        rel: 1e-10,
    };
    const MAX_ITER: usize = 100;

    // ── 整数域测试 ──

    #[test]
    fn test_int_bisection_basic() {
        // f(i) = i - 5, root at i=5
        let f = |i: i64| (i - 5) as f64;
        let r = solve_bracket_int(0, 10, &f, 20).unwrap();
        assert!(r.has_root());
        assert!(r.hi - r.lo <= 1);
        // bracket should contain 5
        assert!(r.lo <= 5 && 5 <= r.hi);
    }

    #[test]
    fn test_int_bisection_negative() {
        let f = |i: i64| (i + 3) as f64;
        let r = solve_bracket_int(-10, 0, &f, 20).unwrap();
        assert!(r.lo <= -3 && -3 <= r.hi);
    }

    #[test]
    fn test_int_bisection_large_range() {
        let f = |i: i64| (i - 500) as f64;
        let r = solve_bracket_int(0, 1000, &f, 20).unwrap();
        assert!(r.lo <= 500 && 500 <= r.hi);
    }

    #[test]
    fn test_int_bisection_best_index() {
        let f = |i: i64| (i - 7) as f64;
        let r = solve_bracket_int(0, 10, &f, 20).unwrap();
        // f(6) = -1, f(7) = 0 -> best is 7
        assert_eq!(r.best_index(), 7);
    }

    #[test]
    fn test_int_bisection_not_bracketed() {
        let f = |i: i64| (i * i) as f64; // always >= 0
        let r = solve_bracket_int(-5, 5, &f, 20);
        assert!(r.is_err());
    }

    #[test]
    fn test_int_bisection_invalid_interval() {
        let f = |i: i64| i as f64;
        let r = solve_bracket_int(5, 3, &f, 20);
        assert!(r.is_err());
    }

    // ── 连续域测试 (原有) ──

    #[test]
    fn test_max_iterations_predict() {
        let n1 = max_iterations(10 - 1, 1, 2);
        assert_eq!(n1, 4);
        let n2 = max_iterations(9.2, 1.01, 2);
        assert_eq!(n2, 4);
    }

    #[test]
    fn test_invalid_interval() {
        let f = |x: f64| x * x - 1.0;

        // a > b
        let result = solve_bracket(2.0, 1.0, &f, Tolerance::new(1e-8, 1e-8), 100);
        assert!(matches!(result, Err(RootFindingError::InvalidInterval)));
    }

    #[test]
    fn test_function_not_bracketed() {
        let f = |x: f64| x * x + 1.0; // 没有实数根，f(x) > 0 always

        let result = solve_bracket(0.0, 2.0, &f, Tolerance::new(1e-8, 1e-8), 100);
        assert!(matches!(
            result,
            Err(RootFindingError::FunctionNotBracketed)
        ));
    }

    #[test]
    fn test_convergence_failure() {
        let f = |x: f64| x.powi(3) - 2.0 * x - 5.0;

        // 设置 max_iter = 3，强制无法收敛
        let result = solve_bracket(2.0, 3.0, &f, Tolerance::new(1e-12, 1e-12), 3);
        let root_bracket = result.unwrap();
        assert!(!root_bracket.is_converged());
        assert!(root_bracket.has_root());
    }

    #[test]
    fn test_cos_x_equals_x() {
        let f = |x: f64| x.cos() - x;
        let result = solve_bracket(0.0, 1.0, &f, TOL, MAX_ITER).unwrap();
        assert!(result.converged);
        assert!((result.b - result.a) < 2.0 * TOL.to_f64(1.0));
        assert!(result.fa * result.fb <= 0.0); // 异号
    }

    #[test]
    fn test_cubic_polynomial() {
        let f = |x: f64| x.powi(3) - x - 2.0;
        let result = solve_bracket(1.0, 2.0, &f, TOL, MAX_ITER).unwrap();
        assert!(result.converged);
        assert!((result.b - result.a) < 2.0 * TOL.to_f64(1.0));
        assert!(result.fa * result.fb <= 0.0);
    }

    #[test]
    fn test_oscillatory_function() {
        // sin(1/x)，避开 x=0，在 [0.1, 1.0] 内有多个根，我们选一个能收敛的区间
        let f = |x: f64| (1.0 / x).sin();
        let result = solve_bracket(0.1, 1.0, &f, TOL, MAX_ITER).unwrap();
        assert!(result.converged);
        assert!((result.b - result.a) < 2.0 * TOL.to_f64(0.9));
        assert!(result.fa * result.fb <= 0.0);
    }

    #[test]
    fn test_steep_function_near_root() {
        // (x - 1)^2 * log10(x)
        let f = |x: f64| (x - 1.0).powi(2) * x.log10();
        let result = solve_bracket(0.5, 2.0, &f, TOL, MAX_ITER).unwrap();
        assert!(result.converged);
        assert!((result.b - result.a) < 2.0 * TOL.to_f64(1.5));
        assert!(result.fa * result.fb <= 0.0);
    }

    #[test]
    fn test_root_at_left_endpoint() {
        let f = |x: f64| x; // f(0) = 0
        let result = solve_bracket(0.0, 1.0, &f, TOL, MAX_ITER).unwrap();
        assert!(result.converged);
        assert!(result.a == 0.0 && result.b == 0.0);
        assert!(result.fa == 0.0 && result.fb == 0.0);
    }

    #[test]
    fn test_root_at_right_endpoint() {
        let f = |x: f64| x - 2.0; // f(2) = 0
        let result = solve_bracket(1.0, 2.0, &f, TOL, MAX_ITER).unwrap();
        assert!(result.converged);
        assert!(result.a == 2.0 && result.b == 2.0);
        assert!(result.fa == 0.0 && result.fb == 0.0);
    }

    #[test]
    fn test_small_initial_interval() {
        let f = |x: f64| x.cos() - x;
        let result = solve_bracket(0.7390851332, 0.7390851333, &f, TOL, MAX_ITER).unwrap();
        assert!(result.converged);
        assert!((result.b - result.a) < 2.0 * TOL.to_f64(1.0));
        assert!(result.fa * result.fb <= 0.0);
    }

    #[test]
    fn test_unbracketed_function_should_fail() {
        let f = |x: f64| x.powi(2); // f(x) = x^2 >= 0
        let result = solve_bracket(-1.0, 1.0, &f, TOL, MAX_ITER);
        assert!(matches!(
            result,
            Err(RootFindingError::FunctionNotBracketed)
        ));
    }

    #[test]
    fn test_multiple_roots_in_interval() {
        let f = |x: f64| (x - 1.0) * (x - 2.0) * (x - 3.0); // 有三个根：1, 2, 3
        let result = solve_bracket(0.5, 3.1, &f, TOL, MAX_ITER).unwrap();
        assert!(result.converged);
        assert!((result.b - result.a) < 2.0 * TOL.to_f64(1.0));
        assert!(result.fa * result.fb <= 0.0);
    }

    #[test]
    fn test_discontinuous_function() {
        let f = |x: f64| if x < 0.0 { -1.0 } else { 1.0 }; // sign(x)
        let result = solve_bracket(-1.0, 1.0, &f, TOL, MAX_ITER).unwrap();
        assert!(result.converged);
        assert!((result.b - result.a) < 2.0 * TOL.to_f64(1.0));
        assert!(result.fa * result.fb <= 0.0); // TOMS748 可以处理跳跃间断点
    }

    #[test]
    fn test_high_precision_requirement() {
        let f = |x: f64| x.cos() - x;
        let tol = Tolerance::new(1e-15, 1e-15);
        let result = solve_bracket(0.0, 1.0, &f, tol, MAX_ITER).unwrap();
        assert!(result.converged);
        assert!((result.b - result.a) < 2.0 * tol.to_f64(1.0));
    }

    #[test]
    fn test_low_max_iterations() {
        let f = |x: f64| x.cos() - x;
        let result = solve_bracket(0.0, 1.0, &f, TOL, 3).unwrap(); // 迭代次数限制
        assert!(!result.converged); // 不应收敛
        assert!(result.iterations == 3);
        assert!(result.fa * result.fb <= 0.0); // 但仍包含根
    }
}
