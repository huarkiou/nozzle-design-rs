use crate::rootfinding::{RootBracket, RootFindingError};

/// TOMS748 一维方程求根算法，返回包含根的区间信息结构体
///
/// # 参数
/// - `a`, `b`: 初始区间，要求函数f在区间端点的函数值符号相反
/// - `f`: 连续函数
/// - `abs_tol`: 绝对误差容限
/// - `rel_tol`: 相对误差容限
/// - `max_iter`: 最大迭代次数
pub fn solve_bracket<F>(
    a: f64,
    b: f64,
    f: &F,
    abs_tol: f64,
    rel_tol: f64,
    max_iter: usize,
) -> Result<RootBracket, RootFindingError>
where
    F: Fn(f64) -> f64,
{
    if a > b {
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

    // 保证 fa <= 0 <= fb
    if f_left > 0.0 {
        std::mem::swap(&mut left, &mut right);
        std::mem::swap(&mut f_left, &mut f_right);
    }

    let mut c = left;
    let mut fc = f_left;
    let mut iterations = 0;

    for _ in 0..max_iter {
        iterations += 1;

        let m = 0.5 * (left + right);
        let tol = rel_tol * (right - left).abs() + abs_tol;

        if (right - left).abs() < 2.0 * tol {
            return Ok(RootBracket {
                a: left,
                b: right,
                fa: f_left,
                fb: f_right,
                iterations,
                converged: true,
            });
        }

        // IQI 插值尝试
        if f_left != f_right && f_left != fc && f_right != fc {
            let q = f_left / f_right;
            let r = fc / f_right;
            let s = fc / f_left;

            let p1 = r * (r - q) * (left - m) + (1.0 - q) * r * (c - m);
            let q_val = (r - 1.0) * (q - 1.0) * (s - 1.0);

            if q_val != 0.0 {
                let s_iqi = m + p1 / q_val;

                if (s_iqi - left) * (s_iqi - right) < 0.0 {
                    let fs = f(s_iqi);
                    if fs == 0.0 {
                        return Ok(RootBracket {
                            a: s_iqi,
                            b: s_iqi,
                            fa: 0.0,
                            fb: 0.0,
                            iterations,
                            converged: true,
                        });
                    }
                    if fs.signum() == f_left.signum() {
                        left = s_iqi;
                        f_left = fs;
                    } else {
                        right = s_iqi;
                        f_right = fs;
                    }
                    c = m;
                    fc = f(c);
                    continue;
                }
            }
        }

        // 割线法尝试
        let s_sec = right - f_right * (right - left) / (f_right - f_left);
        if (s_sec - left) * (s_sec - right) < 0.0 {
            let fs = f(s_sec);
            if fs == 0.0 {
                return Ok(RootBracket {
                    a: s_sec,
                    b: s_sec,
                    fa: 0.0,
                    fb: 0.0,
                    iterations,
                    converged: true,
                });
            }
            if fs.signum() == f_left.signum() {
                left = s_sec;
                f_left = fs;
            } else {
                right = s_sec;
                f_right = fs;
            }
        } else {
            // 回退到二分法
            let s_bi = 0.5 * (left + right);
            let fs = f(s_bi);
            if fs == 0.0 {
                return Ok(RootBracket {
                    a: s_bi,
                    b: s_bi,
                    fa: 0.0,
                    fb: 0.0,
                    iterations,
                    converged: true,
                });
            }
            if fs.signum() == f_left.signum() {
                left = s_bi;
                f_left = fs;
            } else {
                right = s_bi;
                f_right = fs;
            }
        }
    }

    // 达到最大迭代次数，仍返回当前区间
    Ok(RootBracket {
        a: left,
        b: right,
        fa: f_left,
        fb: f_right,
        iterations,
        converged: false,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    const ABS_TOL: f64 = 1e-10;
    const REL_TOL: f64 = 1e-10;
    const MAX_ITER: usize = 100;

    #[test]
    fn test_invalid_interval() {
        let f = |x: f64| x * x - 1.0;

        // a > b
        let result = solve_bracket(2.0, 1.0, &f, 1e-8, 1e-8, 100);
        assert!(matches!(result, Err(RootFindingError::InvalidInterval)));
    }

    #[test]
    fn test_function_not_bracketed() {
        let f = |x: f64| x * x + 1.0; // 没有实数根，f(x) > 0 always

        let result = solve_bracket(0.0, 2.0, &f, 1e-8, 1e-8, 100);
        assert!(matches!(
            result,
            Err(RootFindingError::FunctionNotBracketed)
        ));
    }

    #[test]
    fn test_convergence_failure() {
        let f = |x: f64| x.powi(3) - 2.0 * x - 5.0;

        // 设置 max_iter = 3，强制无法收敛
        let result = solve_bracket(2.0, 3.0, &f, 1e-12, 1e-12, 3);
        let root_bracket = result.unwrap();
        assert!(!root_bracket.is_converged());
        assert!(root_bracket.has_root());
    }

    #[test]
    fn test_cos_x_equals_x() {
        let f = |x: f64| x.cos() - x;
        let result = solve_bracket(0.0, 1.0, &f, ABS_TOL, REL_TOL, MAX_ITER).unwrap();
        assert!(result.converged);
        assert!((result.b - result.a) < 2.0 * (REL_TOL * (1.0) + ABS_TOL));
        assert!(result.fa * result.fb <= 0.0); // 异号
    }

    #[test]
    fn test_cubic_polynomial() {
        let f = |x: f64| x.powi(3) - x - 2.0;
        let result = solve_bracket(1.0, 2.0, &f, ABS_TOL, REL_TOL, MAX_ITER).unwrap();
        assert!(result.converged);
        assert!((result.b - result.a) < 2.0 * (REL_TOL * (1.0) + ABS_TOL));
        assert!(result.fa * result.fb <= 0.0);
    }

    #[test]
    fn test_oscillatory_function() {
        // sin(1/x)，避开 x=0，在 [0.1, 1.0] 内有多个根，我们选一个能收敛的区间
        let f = |x: f64| (1.0 / x).sin();
        let result = solve_bracket(0.1, 1.0, &f, ABS_TOL, REL_TOL, MAX_ITER).unwrap();
        assert!(result.converged);
        assert!((result.b - result.a) < 2.0 * (REL_TOL * (0.9) + ABS_TOL));
        assert!(result.fa * result.fb <= 0.0);
    }

    #[test]
    fn test_steep_function_near_root() {
        // (x - 1)^2 * log10(x)
        let f = |x: f64| (x - 1.0).powi(2) * x.log10();
        let result = solve_bracket(0.5, 2.0, &f, ABS_TOL, REL_TOL, MAX_ITER).unwrap();
        assert!(result.converged);
        assert!((result.b - result.a) < 2.0 * (REL_TOL * (1.5) + ABS_TOL));
        assert!(result.fa * result.fb <= 0.0);
    }

    #[test]
    fn test_root_at_left_endpoint() {
        let f = |x: f64| x; // f(0) = 0
        let result = solve_bracket(0.0, 1.0, &f, ABS_TOL, REL_TOL, MAX_ITER).unwrap();
        assert!(result.converged);
        assert!(result.a == 0.0 && result.b == 0.0);
        assert!(result.fa == 0.0 && result.fb == 0.0);
    }

    #[test]
    fn test_root_at_right_endpoint() {
        let f = |x: f64| x - 2.0; // f(2) = 0
        let result = solve_bracket(1.0, 2.0, &f, ABS_TOL, REL_TOL, MAX_ITER).unwrap();
        assert!(result.converged);
        assert!(result.a == 2.0 && result.b == 2.0);
        assert!(result.fa == 0.0 && result.fb == 0.0);
    }

    #[test]
    fn test_small_initial_interval() {
        let f = |x: f64| x.cos() - x;
        let result =
            solve_bracket(0.7390851332, 0.7390851333, &f, ABS_TOL, REL_TOL, MAX_ITER).unwrap();
        assert!(result.converged);
        assert!((result.b - result.a) < 2.0 * (REL_TOL + ABS_TOL));
        assert!(result.fa * result.fb <= 0.0);
    }

    #[test]
    fn test_unbracketed_function_should_fail() {
        let f = |x: f64| x.powi(2); // f(x) = x^2 >= 0
        let result = solve_bracket(-1.0, 1.0, &f, ABS_TOL, REL_TOL, MAX_ITER);
        assert!(matches!(
            result,
            Err(RootFindingError::FunctionNotBracketed)
        ));
    }

    #[test]
    fn test_multiple_roots_in_interval() {
        let f = |x: f64| (x - 1.0) * (x - 2.0) * (x - 3.0); // 有三个根：1, 2, 3
        let result = solve_bracket(0.5, 3.1, &f, ABS_TOL, REL_TOL, MAX_ITER).unwrap();
        assert!(result.converged);
        assert!((result.b - result.a) < 2.0 * (REL_TOL + ABS_TOL));
        assert!(result.fa * result.fb <= 0.0);
    }

    #[test]
    fn test_discontinuous_function() {
        let f = |x: f64| if x < 0.0 { -1.0 } else { 1.0 }; // sign(x)
        let result = solve_bracket(-1.0, 1.0, &f, ABS_TOL, REL_TOL, MAX_ITER).unwrap();
        assert!(result.converged);
        assert!((result.b - result.a) < 2.0 * (REL_TOL + ABS_TOL));
        assert!(result.fa * result.fb <= 0.0); // TOMS748 可以处理跳跃间断点
    }

    #[test]
    fn test_high_precision_requirement() {
        let f = |x: f64| x.cos() - x;
        let abs_tol = 1e-15;
        let rel_tol = 1e-15;
        let result = solve_bracket(0.0, 1.0, &f, abs_tol, rel_tol, MAX_ITER).unwrap();
        assert!(result.converged);
        assert!((result.b - result.a) < 2.0 * (rel_tol + abs_tol));
    }

    #[test]
    fn test_low_max_iterations() {
        let f = |x: f64| x.cos() - x;
        let result = solve_bracket(0.0, 1.0, &f, ABS_TOL, REL_TOL, 3).unwrap(); // 迭代次数限制
        assert!(!result.converged); // 不应收敛
        assert!(result.iterations == 3);
        assert!(result.fa * result.fb <= 0.0); // 但仍包含根
    }
}
