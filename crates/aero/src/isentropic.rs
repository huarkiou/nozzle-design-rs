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
