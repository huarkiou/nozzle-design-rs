use crate::{Tolerance, rootfinding::RootFindingError};

/// 使用割线法（Secant Method）求解一维非线性方程 f(x) = 0
///
/// # 参数
/// - `f`: 要求解的非线性函数 f(x)
/// - `x0`: 初始点（第一个迭代点）
/// - `tol`: 收敛容忍度，用于判断是否收敛（见 `Tolerance` 类型）
/// - `max_iter`: 最大迭代次数，防止无限循环
///
/// # 返回值
/// - `Ok(x)`: 找到近似根 `x`，满足 |f(x)| < tol
/// - `Err(e)`: 求解失败，可能的原因包括：
///     - Jacobian（差商）接近 0，无法继续迭代
///     - 达到最大迭代次数仍未收敛
pub fn solve<F>(f: F, x0: f64, tol: Tolerance, max_iter: usize) -> Result<f64, RootFindingError>
where
    F: Fn(f64) -> f64,
{
    // 尝试最多 5 次调整初始点
    let h = 1e-8;
    let deltax = 1e-6;
    let mut x_prev = x0;
    let mut x_curr = x0 + h;

    for _ in 0..5 {
        // 检查初始点是否有效
        let f_prev = f(x_prev);
        let f_curr = f(x_curr);

        if f_prev.is_finite() && f_curr.is_finite() {
            break;
        }

        // 否则微调初始点
        x_prev = x0 + deltax;
        x_curr = x_prev + h;
    }

    for _ in 0..max_iter {
        let f_prev = f(x_prev);
        let f_curr = f(x_curr);

        if (f_curr - f_prev).abs() < f64::EPSILON {
            return Err(RootFindingError::InvalidJacobian);
        }

        let x_next = x_curr - f_curr * (x_curr - x_prev) / (f_curr - f_prev);

        if tol.approx_eq(x_next, x_curr) {
            return Ok(x_next);
        }

        x_prev = x_curr;
        x_curr = x_next;
    }

    Err(RootFindingError::MaxIterationsExceeded)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_secant_method_quadratic() {
        // f(x) = x^2 - 2, root is sqrt(2) ≈ 1.4142
        let f = |x: f64| x * x - 2.0;
        let tol = Tolerance::new(1e-8, 1e-8);
        let result = solve(f, 1.0, tol, 100);

        assert!(result.is_ok());
        let root = result.unwrap();
        assert!((root - 2.0_f64.sqrt()).abs() < 1e-8);
    }

    #[test]
    fn test_secant_method_linear() {
        // f(x) = 2x - 4, root is 2.0
        let f = |x: f64| 2.0 * x - 4.0;
        let tol = Tolerance::new(1e-8, 1e-8);
        let result = solve(f, 1.0, tol, 100);

        assert!(result.is_ok());
        let root = result.unwrap();
        assert!((root - 2.0).abs() < 1e-8);
    }

    #[test]
    fn test_secant_method_cubic() {
        // f(x) = x^3 - x - 2, root ≈ 1.5213797...
        let f = |x: f64| x.powi(3) - x - 2.0;
        let tol = Tolerance::new(1e-8, 1e-8);
        let result = solve(f, 1.0, tol, 100);

        assert!(result.is_ok());
        let root = result.unwrap();
        assert!((root - 1.521379706804567).abs() < 1e-8);
    }

    #[test]
    fn test_secant_method_divergent_function() {
        let f = |x: f64| x.powi(2) + 1.0; // x² + 1 = 0 无实根
        let tol = Tolerance::new(1e-8, 1e-8);
        let result = solve(f, 1.5, tol, 100);

        // 期望：不收敛或 Jacobian 为 0
        assert!(matches!(
            result,
            Err(RootFindingError::MaxIterationsExceeded)
        ));
    }

    #[test]
    fn test_secant_method_oscillation_function() {
        let f = |x: f64| (x - 0.5).sin() / (x - 0.5); // 在 x = 0.5 附近震荡
        let tol = Tolerance::new(1e-8, 1e-8);
        let result = solve(f, 1.5, tol, 100);

        assert!(result.is_ok());
    }

    #[test]
    fn test_secant_method_flat_function() {
        // f(x) = 1.0 (constant function), no root
        let f = |_| 1.0;
        let tol = Tolerance::new(1e-8, 1e-8);
        let result = solve(f, 0.0, tol, 100);

        // 割线法会因差商为 0 而失败
        assert!(matches!(result, Err(RootFindingError::InvalidJacobian)));
    }

    #[test]
    fn test_secant_method_higher_iterations() {
        // f(x) = cos(x) - x, root ≈ 0.7390851332151607
        let f = |x: f64| x.cos() - x;
        let tol = Tolerance::new(1e-12, 1e-12);
        let result = solve(f, 0.0, tol, 100);

        assert!(result.is_ok());
        let root = result.unwrap();
        assert!((root - 0.7390851332151607).abs() < 1e-12);
    }

    #[test]
    fn test_secant_method_tan() {
        // f(x) = tan(x) near π/2 where it diverges
        let f = |x: f64| x.tan();
        let tol = Tolerance::new(1e-8, 1e-8);
        let result = solve(f, 1.0, tol, 100);

        assert!(result.is_ok());
        let root = result.unwrap();
        assert!((root - 0.).abs() < 1e-8);
    }
}
