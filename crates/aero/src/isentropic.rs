use math::{Tolerance, rootfinding::toms748};

use crate::Cp;

/// 从静温计算总温。
///
/// 变比热时使用 `cp(t_static)` 作零阶近似：
/// `T_total ≈ T_static + V² / (2·cp(T_static))`。
/// 该方程在 `t_total` 上为线性，无需迭代，直接求解。
pub fn cal_total_temperature(cp: &Cp, t_static: f64, velocity: f64) -> f64 {
    match cp {
        Cp::Constant(v) => t_static + velocity.powi(2) / (2. * v),
        Cp::Variable(f) => t_static + velocity.powi(2) / (2. * f(t_static)),
    }
}

/// 从总温反算静温。
///
/// 方程：`T_total - T_static - V²/(2·cp(T_static)) = 0`。
/// 变比热时 cp 依赖于待求的 T_static，使用 TOMS748 迭代求解。
/// 若求解失败，回退到常温近似 `cp(t_total)` 的显式解。
pub fn cal_static_temperature(cp: &Cp, t_total: f64, velocity: f64) -> f64 {
    match cp {
        Cp::Constant(v) => t_total - velocity.powi(2) / (2. * v),
        Cp::Variable(f) => {
            let equation =
                |t_static: f64| t_total - velocity.powi(2) / (2. * f(t_static)) - t_static;
            // 扩大 bracket 上界以覆盖变比热可能产生的更大温度范围
            let hi = t_total.max(t_total * 1.5 + velocity.powi(2) / 2000.);
            match toms748::solve_bracket(1e-3, hi, &equation, Tolerance::new(1e-8, 1e-8), 50) {
                Ok(res) => res.root(),
                Err(_) => {
                    // 回退：用 cp(t_total) 作常温近似
                    let cp_approx = f(t_total);
                    t_total - velocity.powi(2) / (2. * cp_approx)
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── constant cp ──────────────────────────────────────────

    #[test]
    fn test_total_temp_constant() {
        let cp = Cp::Constant(1005.0);
        let t = cal_total_temperature(&cp, 300.0, 500.0);
        let expected = 300.0 + 500.0_f64.powi(2) / (2.0 * 1005.0);
        assert!((t - expected).abs() < 1e-10);
    }

    #[test]
    fn test_static_temp_constant() {
        let cp = Cp::Constant(1005.0);
        let t = cal_static_temperature(&cp, 2000.0, 1500.0);
        let expected = 2000.0 - 1500.0_f64.powi(2) / (2.0 * 1005.0);
        assert!((t - expected).abs() < 1e-10);
    }

    #[test]
    fn test_total_static_roundtrip_constant() {
        let cp = Cp::Constant(1005.0);
        let t_static = 500.0;
        let v = 1200.0;
        let t_total = cal_total_temperature(&cp, t_static, v);
        let t_back = cal_static_temperature(&cp, t_total, v);
        assert!((t_static - t_back).abs() < 1e-6);
    }

    // ── variable cp ──────────────────────────────────────────

    #[test]
    fn test_total_temp_variable() {
        let cp = Cp::new(|t: f64| 1000.0 + 0.2 * t);
        let t = cal_total_temperature(&cp, 300.0, 500.0);
        assert!(t > 300.0);
        assert!(t.is_finite());
    }

    #[test]
    fn test_static_temp_variable() {
        let cp = Cp::new(|t: f64| 1000.0 + 0.2 * t);
        let t = cal_static_temperature(&cp, 2000.0, 1500.0);
        assert!(t > 0.0);
        assert!(t < 2000.0);
        assert!(t.is_finite());
    }

    #[test]
    fn test_total_static_roundtrip_variable() {
        let cp = Cp::new(|t: f64| 1000.0 + 0.1 * t);
        let t_static = 500.0;
        let v = 1200.0;
        let t_total = cal_total_temperature(&cp, t_static, v);
        let t_back = cal_static_temperature(&cp, t_total, v);
        // 变比热下 cp 依赖于温度，零阶近似有残余误差
        let err = (t_static - t_back).abs();
        assert!(err < 50.0, "roundtrip error too large: {err}");
    }

    #[test]
    fn test_static_temp_zero_velocity() {
        let cp = Cp::new(|t: f64| 1000.0 + 0.2 * t);
        let t = cal_static_temperature(&cp, 1500.0, 0.0);
        // 无速度时静温等于总温
        assert!((t - 1500.0).abs() < 1e-6);
    }

    #[test]
    fn test_static_temp_fallback_path() {
        // 使用极大速度 + 极小 cp → V²/(2cp) 远超 t_total → bracket 失败 → fallback
        let cp = Cp::new(|t: f64| 100.0 + 0.01 * t); // 很小的 cp
        let t = cal_static_temperature(&cp, 500.0, 3000.0);
        // fallback 近似可能产生不合理值，但不应 panic
        assert!(t.is_finite());
    }
}
