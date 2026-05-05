use num_traits::{Num, ToPrimitive};

use crate::{
    rootfinding::{RootBracket, RootFindingError},
    Tolerance,
};

// ── 通用二分核心 ──────────────────────────────────────────────

/// 通用二分搜索核心，同时服务于连续域（f64）和离散域（i64）。
///
/// 参数化中点计算和收敛判定，避免重复实现。
fn bisect_core<T, F, M, C>(
    lo: T,
    hi: T,
    f: &F,
    midpoint: M,
    is_done: C,
    max_iter: usize,
) -> Result<(T, T, f64, f64, usize), RootFindingError>
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

    Ok((left, right, f_left, f_right, iterations))
}

// ── 整数域二分搜索 ────────────────────────────────────────────

/// 整数域二分搜索结果
#[derive(Debug, Clone, Copy)]
pub struct IntBracket {
    pub lo: i64,
    pub hi: i64,
    pub flo: f64,
    pub fhi: f64,
    pub iterations: usize,
}

impl IntBracket {
    pub fn has_root(&self) -> bool {
        self.flo * self.fhi <= 0.0
    }

    pub fn best_index(&self) -> i64 {
        if self.flo.abs() < self.fhi.abs() {
            self.lo
        } else {
            self.hi
        }
    }

    pub fn best_value(&self) -> f64 {
        if self.flo.abs() < self.fhi.abs() {
            self.flo
        } else {
            self.fhi
        }
    }
}

/// 整数域二分法查找函数根所在区间。
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

    let (left, right, f_left, f_right, iterations) =
        bisect_core(lo, hi, f, |a, b| (a + b) / 2, |a, b| b - a <= 1, max_iter)?;

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

/// 二分搜索法查找函数根（连续域）。
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

    // 检查端点是否为精确根
    let fa = f(a);
    let fb = f(b);
    if fa == 0.0 {
        return Ok(RootBracket {
            a,
            b: a,
            fa,
            fb: fa,
            iterations: 0,
            converged: true,
        });
    }
    if fb == 0.0 {
        return Ok(RootBracket {
            a: b,
            b,
            fa: fb,
            fb,
            iterations: 0,
            converged: true,
        });
    }

    let (left, right, f_left, f_right, iter) = bisect_core(
        a,
        b,
        f,
        |x, y| (x + y) / 2.0,
        |x, y| tol.approx_eq(x, y),
        max_iter,
    )?;

    // 检查中点上是否恰好命中根
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
        let f = |i: i64| (i - 5) as f64;
        let r = solve_bracket_int(0, 10, &f, 20).unwrap();
        assert!(r.has_root());
        assert!(r.hi - r.lo <= 1);
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
        assert_eq!(r.best_index(), 7);
    }

    #[test]
    fn test_int_bisection_not_bracketed() {
        let f = |i: i64| (i * i) as f64;
        let r = solve_bracket_int(-5, 5, &f, 20);
        assert!(r.is_err());
    }

    #[test]
    fn test_int_bisection_invalid_interval() {
        let f = |i: i64| i as f64;
        let r = solve_bracket_int(5, 3, &f, 20);
        assert!(r.is_err());
    }

    // ── 连续域测试 ──

    #[test]
    fn test_max_iterations_predict() {
        assert_eq!(max_iterations(10 - 1, 1, 2), 4);
        assert_eq!(max_iterations(9.2, 1.01, 2), 4);
    }

    #[test]
    fn test_invalid_interval() {
        let f = |x: f64| x * x - 1.0;
        let result = solve_bracket(2.0, 1.0, &f, Tolerance::new(1e-8, 1e-8), 100);
        assert!(matches!(result, Err(RootFindingError::InvalidInterval)));
    }

    #[test]
    fn test_function_not_bracketed() {
        let f = |x: f64| x * x + 1.0;
        let result = solve_bracket(0.0, 2.0, &f, Tolerance::new(1e-8, 1e-8), 100);
        assert!(matches!(
            result,
            Err(RootFindingError::FunctionNotBracketed)
        ));
    }

    #[test]
    fn test_convergence_failure() {
        let f = |x: f64| x.powi(3) - 2.0 * x - 5.0;
        let result = solve_bracket(2.0, 3.0, &f, Tolerance::new(1e-12, 1e-12), 3);
        let rb = result.unwrap();
        assert!(!rb.is_converged());
        assert!(rb.has_root());
    }

    #[test]
    fn test_cos_x_equals_x() {
        let f = |x: f64| x.cos() - x;
        let result = solve_bracket(0.0, 1.0, &f, TOL, MAX_ITER).unwrap();
        assert!(result.converged);
        assert!((result.b - result.a) < 2.0 * TOL.to_f64(1.0));
        assert!(result.fa * result.fb <= 0.0);
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
        let f = |x: f64| (1.0 / x).sin();
        let result = solve_bracket(0.1, 1.0, &f, TOL, MAX_ITER).unwrap();
        assert!(result.converged);
        assert!((result.b - result.a) < 2.0 * TOL.to_f64(0.9));
        assert!(result.fa * result.fb <= 0.0);
    }

    #[test]
    fn test_steep_function_near_root() {
        let f = |x: f64| (x - 1.0).powi(2) * x.log10();
        let result = solve_bracket(0.5, 2.0, &f, TOL, MAX_ITER).unwrap();
        assert!(result.converged);
        assert!((result.b - result.a) < 2.0 * TOL.to_f64(1.5));
        assert!(result.fa * result.fb <= 0.0);
    }

    #[test]
    fn test_root_at_left_endpoint() {
        let f = |x: f64| x;
        let result = solve_bracket(0.0, 1.0, &f, TOL, MAX_ITER).unwrap();
        assert!(result.converged);
        assert!(result.a == 0.0 && result.b == 0.0);
        assert!(result.fa == 0.0 && result.fb == 0.0);
    }

    #[test]
    fn test_root_at_right_endpoint() {
        let f = |x: f64| x - 2.0;
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
        let f = |x: f64| x.powi(2);
        let result = solve_bracket(-1.0, 1.0, &f, TOL, MAX_ITER);
        assert!(matches!(
            result,
            Err(RootFindingError::FunctionNotBracketed)
        ));
    }

    #[test]
    fn test_multiple_roots_in_interval() {
        let f = |x: f64| (x - 1.0) * (x - 2.0) * (x - 3.0);
        let result = solve_bracket(0.5, 3.1, &f, TOL, MAX_ITER).unwrap();
        assert!(result.converged);
        assert!((result.b - result.a) < 2.0 * TOL.to_f64(1.0));
        assert!(result.fa * result.fb <= 0.0);
    }

    #[test]
    fn test_discontinuous_function() {
        let f = |x: f64| if x < 0.0 { -1.0 } else { 1.0 };
        let result = solve_bracket(-1.0, 1.0, &f, TOL, MAX_ITER).unwrap();
        assert!(result.converged);
        assert!((result.b - result.a) < 2.0 * TOL.to_f64(1.0));
        assert!(result.fa * result.fb <= 0.0);
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
        let result = solve_bracket(0.0, 1.0, &f, TOL, 3).unwrap();
        assert!(!result.converged);
        assert!(result.iterations == 3);
        assert!(result.fa * result.fb <= 0.0);
    }
}
