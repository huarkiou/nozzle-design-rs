use crate::moc::{
    MocPoint,
    unitprocess::{Context, ExitLineFunc, GeneralConfig, UnitProcess},
};

/// 无旋特征线法的基本计算过程
pub struct Irrotational {
    pub conf: GeneralConfig,
}

impl Irrotational {
    /// 计算指定点处的左行特征线相关参数 `(L, Q, R, S)`。
    ///
    /// # 参数
    ///
    /// - `point`: `MocPoint`流场中任一点。
    ///
    /// # 返回值
    ///
    /// 返回一个四元组 `(L, Q, R, S)`：
    ///
    /// - `L`: 特征线斜率，由流动方向与马赫角决定。
    /// - `Q`: 流动速度与声速关系的中间变量。
    /// - `R`: 线性组合项，结合了速度和L。
    /// - `S`: 轴对称修正项，若配置中启用轴对称则非零。
    fn cal_lqrs_left(&self, point: &MocPoint) -> (f64, f64, f64, f64) {
        let (soundspeed, ma) = point.sound_speed_and_mach_number();
        let l = (point.flow_direction() + (1.0 / ma).asin()).tan();
        let q = point.u * point.u - soundspeed * soundspeed;
        let r = 2.0 * point.u * point.v - q * l;
        let s = if self.conf.axisym {
            soundspeed.powi(2) * point.v / point.y
        } else {
            0.0
        };
        (l, q, r, s)
    }

    /// 计算指定点处的右行特征线相关参数 `(L, Q, R, S)`。
    ///
    /// # 参数
    ///
    /// - `point`: `MocPoint`流场中任一点。
    ///
    /// # 返回值
    ///
    /// 返回一个四元组 `(L, Q, R, S)`：
    ///
    /// - `L`: 特征线斜率，由流动方向与马赫角决定。
    /// - `Q`: 流动速度与声速关系的中间变量。
    /// - `R`: 线性组合项，结合了速度和L。
    /// - `S`: 轴对称修正项，若配置中启用轴对称则非零。
    fn cal_lqrs_right(&self, point: &MocPoint) -> (f64, f64, f64, f64) {
        let (soundspeed, ma) = point.sound_speed_and_mach_number();
        let l = (point.flow_direction() - (1.0 / ma).asin()).tan();
        let q = point.u * point.u - soundspeed * soundspeed;
        let r = 2.0 * point.u * point.v - q * l;
        let s = if self.conf.axisym {
            soundspeed.powi(2) * point.v / point.y
        } else {
            0.0
        };
        (l, q, r, s)
    }

    fn cal_lqrs_left_modified(&self, point: &MocPoint, modi: f64) -> (f64, f64, f64, f64) {
        let (soundspeed, ma) = point.sound_speed_and_mach_number();
        let l = (point.flow_direction() + (1.0 / ma).asin()).tan();
        let q = point.u * point.u - soundspeed * soundspeed;
        let r = 2.0 * point.u * point.v - q * l;
        let s = if self.conf.axisym {
            soundspeed.powi(2) * modi
        } else {
            0.0
        };
        (l, q, r, s)
    }
}

#[inline]
fn tuple_average_3f64(t1: (f64, f64, f64), t2: (f64, f64, f64)) -> (f64, f64, f64) {
    ((t1.0 + t2.0) / 2., (t1.1 + t2.1) / 2., (t1.2 + t2.2) / 2.)
}

impl UnitProcess for Irrotational {
    fn interior_point(&self, context: Context) -> Option<MocPoint> {
        let p1 = &context.next[context.idx_next];
        let p2 = &context.prev[context.idx_prev];
        let (mut lp, mut qp, mut rp, mut sp);
        if p2.y.abs() < self.conf.tol.abs {
            (lp, qp, rp, sp) = self.cal_lqrs_left_modified(p2, p1.v / p1.y);
        } else {
            (lp, qp, rp, sp) = self.cal_lqrs_left(p2);
        }
        let (mut lm, mut qm, mut rm, mut sm) = self.cal_lqrs_right(p1);

        let mut pr = (p1 + p2) / 2.;
        for _ in 0..self.conf.n_corr {
            let pr_prev = pr.clone();

            // 计算 pr.x 和 pr.y
            pr.x = (p1.y - p2.y - lm * p1.x + lp * p2.x) / (lp - lm);
            pr.y = p1.y + lm * (pr.x - p1.x);

            // 计算 Tp, Tm
            let tp = sp * (pr.x - p2.x) + qp * p2.u + rp * p2.v;
            let tm = sm * (pr.x - p1.x) + qm * p1.u + rm * p1.v;

            // 解速度
            let denominator = qm * rp - qp * rm;
            let u = (tm * rp - tp * rm) / denominator;
            let v = (qm * tp - qp * tm) / denominator;

            pr.u = u;
            pr.v = v;
            // 特别处理轴线附近的点
            if self.conf.axisym && pr.y.abs() < self.conf.tol.abs {
                pr.v = 0.;
            }

            // 检查是否收敛
            if pr.is_position_converged_with(&pr_prev, self.conf.tol)
                && pr.is_velocity_converged_with(&pr_prev, self.conf.tol)
            {
                break;
            }

            // 更新 LQRS 使用新中点
            let mid1 = (&pr + p2) / 2.0;
            let mid2 = (&pr + p1) / 2.0;
            (lp, qp, rp, sp) = self.cal_lqrs_left(&mid1);
            (lm, qm, rm, sm) = self.cal_lqrs_right(&mid2);
        }

        let (tt, pt, rt) = tuple_average_3f64(
            p1.total_temperature_pressure_density(),
            p2.total_temperature_pressure_density(),
        );
        pr.set_temperature_pressure_density(tt, pt, rt);

        Some(pr)
    }

    fn symmetry_axis_point(&self, context: Context) -> Option<MocPoint> {
        todo!()
    }

    fn inverse_wall_point(&self, context: Context) -> Option<MocPoint> {
        todo!()
    }

    fn exit_characteristics_point(
        &self,
        context: Context,
        cal_u_v: ExitLineFunc,
    ) -> Option<MocPoint> {
        todo!()
    }

    fn exit_characteristics_point_fixed_dist(
        &self,
        context: Context,
        cal_u_v: ExitLineFunc,
        dist: f64,
    ) -> Option<MocPoint> {
        todo!()
    }

    fn transition_interior_point(&self, context: Context) -> Option<MocPoint> {
        todo!()
    }

    fn transition_wall_point(&self, context: Context) -> Option<MocPoint> {
        todo!()
    }

    fn last_point(&self, context: Context, cal_u_v: ExitLineFunc) -> Option<MocPoint> {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use math::Tolerance;

    use crate::{Material, moc::CharLine};

    use super::*;

    #[test]
    fn test_interior_point_1() {
        let config = GeneralConfig {
            axisym: true,
            tol: Tolerance::new(1e-5, 1e-5),
            n_corr: 20,
        };

        let unitprocess = Irrotational { conf: config };
        let mat = Arc::new(Material::from_rgas_gamma(320.0, 1.2));

        // 构造两个输入点
        let p1 = MocPoint::from_compatible(
            0.131460,
            0.040118,
            2473.4,
            812.8,
            34042.0,
            3000.0,
            0.086151,
            mat.clone(),
        );
        let p2 = MocPoint::from_compatible(
            0.135683,
            0.037123,
            2502.8,
            737.6,
            32781.0,
            3000.0,
            0.083482,
            mat.clone(),
        );

        let velocity = 2628.726210082568_f64;
        let theta = 0.3013949963150419_f64;
        let target = MocPoint::new(
            0.14113488139562955,
            0.040536826493196544,
            velocity * theta.cos(),
            velocity * theta.sin(),
            28742.423476934055,
            1200.4683626106616,
            0.07481934065705345,
            mat.clone(),
        );

        // 创建 CharLine
        let mut next_line = CharLine::new();
        next_line.push(p1.clone());

        let mut prev_line = CharLine::new();
        prev_line.push(p2.clone());

        let context = Context {
            prev: &prev_line,
            next: &next_line,
            idx_prev: 0,
            idx_next: 0,
        };

        let result_point = unitprocess
            .interior_point(context)
            .expect("interior point should be Some");

        assert!(
            result_point.is_converged_with(&target, unitprocess.conf.tol),
            "Interior point did not converge to expected value:\nresult:{:15}\ntarget:{:15}\n  diff:{}",
            result_point,
            target,
            &result_point - &target
        );
    }

    #[test]
    fn test_interior_point_2() {
        let config = GeneralConfig {
            axisym: true,
            tol: Tolerance::new(1e-5, 1e-5),
            n_corr: 20,
        };

        let unitprocess = Irrotational { conf: config };
        let mat = Arc::new(Material::from_rgas_gamma(287.042, 1.4));

        // 构造两个输入点
        let velocity = 1948.3337719140004_f64;
        let theta = 0.43988113776612270_f64;
        let p1 = MocPoint::from_compatible(
            5.9734147955750752,
            1.1257175564978437,
            velocity * theta.cos(),
            velocity * theta.sin(),
            65.123505884851653,
            2000.0,
            0.0010063378824976170,
            mat.clone(),
        );

        let velocity = 1955.5966975668214_f64;
        let theta = 0.47112060764605074_f64;
        let p2 = MocPoint::from_compatible(
            6.0,
            1.1247947107403191,
            velocity * theta.cos(),
            velocity * theta.sin(),
            57.598901602349798,
            2000.0,
            0.00076988609083849971,
            mat.clone(),
        );

        let velocity = 1949.2684506799624_f64;
        let theta = 0.43805199441312853_f64;
        let target = MocPoint::new(
            6.036241021295845,
            1.1473941269582562,
            velocity * theta.cos(),
            velocity * theta.sin(),
            74.5290738178469,
            108.96389835620812,
            0.0010021388327088687,
            mat.clone(),
        );

        // 创建 CharLine
        let mut next_line = CharLine::new();
        next_line.push(p1.clone());

        let mut prev_line = CharLine::new();
        prev_line.push(p2.clone());

        let context = Context {
            prev: &prev_line,
            next: &next_line,
            idx_prev: 0,
            idx_next: 0,
        };

        let result_point = unitprocess
            .interior_point(context)
            .expect("interior point should be Some");

        assert!(
            result_point.is_converged_with(&target, unitprocess.conf.tol),
            "Interior point did not converge to expected value:\nresult:{:15}\ntarget:{:15}\n  diff:{}",
            result_point,
            target,
            &result_point - &target
        );
    }
}
