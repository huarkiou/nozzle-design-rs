/// 2D Newton 法求解 f(x) = 0.
///
/// 使用有限差分 Jacobian，迭代至收敛或达到最大迭代次数。
///
/// # 参数
/// - `f`: 函数 f([x0, x1]) -> [f0, f1]
/// - `x`: 初始猜测 [x0, x1]，会被原地修改为最终解
/// - `eps`: 收敛容限
/// - `max_iter`: 最大迭代次数
///
/// # 返回值
/// `true` 表示收敛，`false` 表示未收敛或失败
pub fn newton_2d<F>(f: &F, x: &mut [f64; 2], eps: f64, max_iter: usize) -> bool
where
    F: Fn(&[f64; 2]) -> [f64; 2],
{
    let h = 1e-6;
    for _ in 0..max_iter {
        let fx = f(x);

        if fx[0].abs() < eps && fx[1].abs() < eps {
            return true;
        }

        // 有限差分 Jacobian
        let mut x_h = *x;
        x_h[0] += h;
        let fx_dv = f(&x_h);
        let df0_dx0 = (fx_dv[0] - fx[0]) / h;
        let df1_dx0 = (fx_dv[1] - fx[1]) / h;

        x_h = *x;
        x_h[1] += h;
        let fx_dt = f(&x_h);
        let df0_dx1 = (fx_dt[0] - fx[0]) / h;
        let df1_dx1 = (fx_dt[1] - fx[1]) / h;

        let det = df0_dx0 * df1_dx1 - df0_dx1 * df1_dx0;
        if det.abs() < 1e-15 {
            return false;
        }

        let dx0 = -(fx[0] * df1_dx1 - fx[1] * df0_dx1) / det;
        let dx1 = -(df0_dx0 * fx[1] - df1_dx0 * fx[0]) / det;

        x[0] += dx0;
        x[1] += dx1;

        // 数值稳定性
        if !x[0].is_finite() || !x[1].is_finite() {
            return false;
        }
    }

    false
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_newton_2d_simple() {
        // f0 = x0^2 + x1^2 - 4 = 0 (circle radius 2)
        // f1 = x0 - x1 = 0 (line)
        // Solution: x0 = x1 = sqrt(2) ≈ 1.414
        let f = |x: &[f64; 2]| -> [f64; 2] { [x[0].powi(2) + x[1].powi(2) - 4.0, x[0] - x[1]] };

        let mut x = [1.0_f64, 1.0_f64];
        let converged = newton_2d(&f, &mut x, 1e-10, 50);
        assert!(converged);
        let expected = 2.0_f64.sqrt();
        assert!((x[0] - expected).abs() < 1e-8);
        assert!((x[1] - expected).abs() < 1e-8);
    }

    #[test]
    fn test_newton_2d_no_solution() {
        // f0 = x0^2 + 1 = 0 (no real solution)
        // f1 = x1 - 1 = 0
        let f = |x: &[f64; 2]| -> [f64; 2] { [x[0].powi(2) + 1.0, x[1] - 1.0] };

        let mut x = [0.0_f64, 0.0_f64];
        let converged = newton_2d(&f, &mut x, 1e-10, 10);
        assert!(!converged);
    }
}
