use crate::moc::{
    CharLine, MocPoint,
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
    #[inline]
    fn cal_lqrs_left(&self, point: &MocPoint) -> (f64, f64, f64, f64) {
        let (soundspeed, ma) = point.sound_speed_and_mach_number();
        let l = (point.flow_direction() + (1.0 / ma).asin()).tan();
        let q = point.u * point.u - soundspeed * soundspeed;
        let r = 2.0 * point.u * point.v - q * l;
        let s = if self.conf.is_axisymmetric() {
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
    #[inline]
    fn cal_lqrs_right(&self, point: &MocPoint) -> (f64, f64, f64, f64) {
        let (soundspeed, ma) = point.sound_speed_and_mach_number();
        let l = (point.flow_direction() - (1.0 / ma).asin()).tan();
        let q = point.u * point.u - soundspeed * soundspeed;
        let r = 2.0 * point.u * point.v - q * l;
        let s = if self.conf.is_axisymmetric() {
            soundspeed.powi(2) * point.v / point.y
        } else {
            0.0
        };
        (l, q, r, s)
    }

    #[inline]
    fn cal_lqrs_left_modified(&self, point: &MocPoint, modi: f64) -> (f64, f64, f64, f64) {
        let (soundspeed, ma) = point.sound_speed_and_mach_number();
        let l = (point.flow_direction() + (1.0 / ma).asin()).tan();
        let q = point.u * point.u - soundspeed * soundspeed;
        let r = 2.0 * point.u * point.v - q * l;
        let s = if self.conf.is_axisymmetric() {
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
            if self.conf.is_axisymmetric() && pr.y.abs() < self.conf.tol.abs {
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
        let p1 = &context.next[context.idx_next];

        let (mut lm, mut qm, mut rm, mut sm) = self.cal_lqrs_right(p1);

        let mut pr = p1.clone();
        for _ in 0..self.conf.n_corr {
            let pr_prev = pr.clone();

            // 计算 pr.x 和 pr.y
            pr.x = p1.x - p1.y / lm;
            pr.y = 0.;

            // 计算 Tp, Tm
            let tm = sm * (pr.x - p1.x) + qm * p1.u + rm * p1.v;

            // 解速度
            let u = tm / qm;
            let v = 0.;

            pr.u = u;
            pr.v = v;

            // 检查是否收敛
            if pr.is_position_converged_with(&pr_prev, self.conf.tol)
                && pr.is_velocity_converged_with(&pr_prev, self.conf.tol)
            {
                break;
            }

            // 更新 LQRS 使用新中点
            let mid = (&pr + p1) / 2.0;
            (lm, qm, rm, sm) = self.cal_lqrs_right(&mid);
        }

        let (tt, pt, rt) = p1.total_temperature_pressure_density();
        pr.set_temperature_pressure_density(tt, pt, rt);

        Some(pr)
    }

    fn inverse_wall_point(&self, context: Context, wall_info: (f64, f64, f64)) -> Option<MocPoint> {
        let p1 = &context.prev[context.idx_prev];
        let p2 = &context.prev[context.idx_prev + 1];

        let lm = (p2.y - p1.y) / (p2.x - p1.x);

        let mut p3 = (p1 + p2) / 2.; // p3为p1与p2之间某一点，预估p3参数

        let velo = p3.velocity();
        let mut pr = MocPoint {
            x: wall_info.0,
            y: wall_info.1,
            u: velo * wall_info.2.cos(),
            v: velo * wall_info.2.sin(),
            ..p3.clone()
        };
        for _ in 0..self.conf.n_corr {
            let pr_prev = pr.clone();

            // 计算此轮迭代中的p3参数
            for _ in 0..self.conf.n_corr {
                let p3_prev = p3.clone();
                let p_tmp = (&p3 + &pr) / 2.;
                let (_, mach) = p_tmp.sound_speed_and_mach_number();
                let lp = (p_tmp.flow_direction() + (1. / mach).asin()).tan();
                p3.x = (p1.y - pr.y + lp * pr.x - lm * p1.x) / (lp - lm);
                p3.y = pr.y + lp * (p3.x - pr.x);
                p3.interpolate_along(p1, p2);
                if !p3.is_valid() {
                    return None;
                }
                if p3.is_converged_with(&p3_prev, self.conf.tol) {
                    break;
                }
            }

            // 沿着左行特征线p3-pr
            let (_, qp, rp, sp) = self.cal_lqrs_left(&((&pr + &p3) / 2.));
            let tp = sp * (pr.x - p3.x) + qp * p3.u + rp * p3.v;
            let l0 = pr.v / pr.u;
            pr.u = tp / (qp + l0 * rp);
            pr.v = pr.u * l0;

            // 检查是否收敛
            if pr.is_position_converged_with(&pr_prev, self.conf.tol)
                && pr.is_velocity_converged_with(&pr_prev, self.conf.tol)
            {
                break;
            }
        }

        let (tt, pt, rt) = tuple_average_3f64(
            p1.total_temperature_pressure_density(),
            p2.total_temperature_pressure_density(),
        );
        pr.set_temperature_pressure_density(tt, pt, rt);

        Some(pr)
    }

    fn exit_characteristics_point(
        &self,
        context: Context,
        cal_u_v: ExitLineFunc,
    ) -> Option<MocPoint> {
        let p1 = &context.next[context.idx_next];
        let p2 = &context.prev[context.idx_prev];
        let mut lp;
        if p2.y.abs() < self.conf.tol.abs {
            (lp, _, _, _) = self.cal_lqrs_left_modified(p2, p1.v / p1.y);
        } else {
            (lp, _, _, _) = self.cal_lqrs_left(p2);
        }

        let mut pr = (p1 + p2) / 2.;
        (pr.u, pr.v) = cal_u_v(pr.y, &p1);
        for _ in 0..self.conf.n_corr {
            let pr_prev = pr.clone();

            // 流线
            let l0 = (p2.flow_direction() + pr.flow_direction()).tan();
            // 计算 pr.x 和 pr.y
            pr.x = (p2.y - p1.y + lp * p1.x - l0 * p2.x) / (lp - l0);
            pr.y = p1.y + lp * (pr.x - p1.x);

            // 解速度
            (pr.u, pr.v) = cal_u_v(pr.y, &p1);

            // 检查是否收敛
            if pr.is_position_converged_with(&pr_prev, self.conf.tol)
                && pr.is_velocity_converged_with(&pr_prev, self.conf.tol)
            {
                break;
            }

            // 更新 LQRS 使用新中点
            let mid1 = (&pr + p1) / 2.0;
            (lp, _, _, _) = self.cal_lqrs_left(&mid1);
        }

        let (tt, pt, rt) = tuple_average_3f64(
            p1.total_temperature_pressure_density(),
            p2.total_temperature_pressure_density(),
        );
        pr.set_temperature_pressure_density(tt, pt, rt);

        Some(pr)
    }

    fn exit_characteristics_point_fixed_dist(
        &self,
        context: Context,
        cal_u_v: ExitLineFunc,
        dist: f64,
    ) -> Option<MocPoint> {
        // let p1 = &context.next[context.idx_next];
        let p2 = &context.prev[context.idx_prev];

        // let mut pr = (p1 + p2) / 2.;
        let mut pr = p2.clone();
        let mut theta_p = p2.flow_direction() + (1. / p2.mach_number()).asin();
        for _ in 0..self.conf.n_corr {
            let pr_prev = pr.clone();

            pr.x = p2.x + dist * theta_p.cos();
            pr.y = p2.y + dist * theta_p.sin();
            (pr.u, pr.v) = cal_u_v(pr.y, p2);

            // 检查是否收敛
            if pr.is_position_converged_with(&pr_prev, self.conf.tol)
                && pr.is_velocity_converged_with(&pr_prev, self.conf.tol)
            {
                break;
            }

            let pt = (&pr + p2) / 2.;
            theta_p = pt.flow_direction() + (1. / pt.mach_number()).asin();
        }

        // let (tt, pt, rt) = tuple_average_3f64(
        //     p1.total_temperature_pressure_density(),
        //     p2.total_temperature_pressure_density(),
        // );
        let (tt, pt, rt) = p2.total_temperature_pressure_density();
        pr.set_temperature_pressure_density(tt, pt, rt);

        Some(pr)
    }

    fn transition_interior_point(&self, context: Context) -> Option<MocPoint> {
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

            // 检查是否收敛
            if pr.is_position_converged_with(&pr_prev, self.conf.tol)
                && pr.is_velocity_converged_with(&pr_prev, self.conf.tol)
            {
                break;
            }

            // 更新 LQRS 使用新中点
            let mid1 = (&pr + p2) / 2.0;
            (lp, qp, rp, sp) = self.cal_lqrs_left(&mid1);
            let mid2 = (&pr + p1) / 2.0;
            (lm, qm, rm, sm) = self.cal_lqrs_right(&mid2);
        }

        if pr.x > p1.x {
            return None;
        }

        let (tt, pt, rt) = tuple_average_3f64(
            p1.total_temperature_pressure_density(),
            p2.total_temperature_pressure_density(),
        );
        pr.set_temperature_pressure_density(tt, pt, rt);

        Some(pr)
    }

    fn transition_wall_point(&self, context: Context) -> Option<MocPoint> {
        let p1 = &context.next[context.idx_next];
        let p2 = &context.prev[context.idx_prev];
        let p3 = &context.prev[context.idx_prev + 1];

        // 右行特征线pr-p1
        let (mut lm, mut qm, mut rm, mut sm) = self.cal_lqrs_right(p1);
        // 流线p2-pr
        let mut l0 = p2.v / p2.u;
        // p4为右行特征线p2-p3上的一点
        let mut p4 = (p2 + p3) / 2.;
        // 从p4发出的左行特征线p4-pr
        let (mut lp, mut qp, mut rp, mut sp) = self.cal_lqrs_left(&p4);

        let mut pr = (p2 + p1) / 2.;
        for _ in 0..self.conf.n_corr {
            let pr_prev = pr.clone();

            // 计算p4
            for _ in 0..self.conf.n_corr {
                let p4_prev = p4.clone();
                let l23 = (p3.y - p2.y) / (p3.x - p2.x);
                p4.x = (p2.y - pr.y + lp * pr.x - l23 * p2.x) / (lp - l23);
                p4.y = pr.y + lp * (p4.x - pr.x);
                p4.interpolate_along(p2, p3);
                if !p4.is_valid() {
                    return None;
                }
                if p4.is_converged_with(&p4_prev, self.conf.tol) {
                    break;
                }
            }

            // 计算 pr.x 和 pr.y
            pr.x = (p1.y - p2.y + l0 * p2.x - lm * p1.x) / (l0 - lm);
            pr.y = p2.y + l0 * (pr.x - p2.x);

            // 解速度
            let tp = sp * (pr.x - p1.x) + qp * p1.u + rp * p1.v;
            let tm = sm * (pr.x - p4.x) + qm * p4.u + rm * p4.v;
            let denominator = qm * rp - qp * rm;
            let u = (tm * rp - tp * rm) / denominator;
            let v = (qm * tp - qp * tm) / denominator;
            pr.u = u;
            pr.v = v;
            // 修正
            let velo = pr.velocity();
            if velo < self.conf.tol.abs {
                pr.u = velo;
                pr.v = 0.;
            }

            // 检查是否收敛
            if pr.is_position_converged_with(&pr_prev, self.conf.tol)
                && pr.is_velocity_converged_with(&pr_prev, self.conf.tol)
            {
                break;
            }

            // 左行特征线p4-pr
            let mid1 = (&p4 + &pr) / 2.0;
            (lp, qp, rp, sp) = self.cal_lqrs_left(&mid1);
            // 右行特征线pr-p1
            let mid2 = (&pr + p1) / 2.0;
            (lm, qm, rm, sm) = self.cal_lqrs_right(&mid2);

            // 流线p2-pr
            l0 = (p2.v / p2.u + pr.v / pr.u) / 2.;
        }

        let (tt, pt, rt) = tuple_average_3f64(
            p1.total_temperature_pressure_density(),
            p2.total_temperature_pressure_density(),
        );
        pr.set_temperature_pressure_density(tt, pt, rt);

        Some(pr)
    }

    fn last_point(
        &self,
        context: Context,
        cal_u_v: ExitLineFunc,
        mfr_need: f64,
    ) -> Option<MocPoint> {
        let p2 = &context.prev[context.idx_prev];
        let p3 = &context.prev[context.idx_prev + 1];

        // 预估p3到pr的距离
        let mut d3r = p3.distance_to(&p2) / 2.;

        let mut pr = p2.clone();
        for _ in 0..self.conf.n_corr {
            let pt = (p3 + &pr) / 2.;
            let theta = pt.flow_direction() + (1. / pt.mach_number()).asin();
            pr.x = p3.x + d3r * theta.cos();
            pr.y = p3.y + d3r * theta.sin();
            (pr.u, pr.v) = cal_u_v(pr.y, &p3);
            let mfr_cur = CharLine::mass_flow_rate(&vec![pr.clone(), p3.clone()], self.conf.axisym);
            if (mfr_need - mfr_cur).abs() < mfr_need * self.conf.tol.rel {
                break;
            }
            d3r *= mfr_need / mfr_cur;
        }

        Some(pr)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        Material,
        moc::{AreaType, CharLine},
    };
    use math::Tolerance;

    #[test]
    fn test_interior_point_1() {
        let config = GeneralConfig {
            axisym: AreaType::Axisymmetric,
            tol: Tolerance::new(1e-5, 1e-5),
            n_corr: 20,
        };

        let unitprocess = Irrotational { conf: config };
        let mat = Material::from_rgas_gamma(320.0, 1.2);

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
            axisym: AreaType::Axisymmetric,
            tol: Tolerance::new(1e-5, 1e-5),
            n_corr: 20,
        };

        let unitprocess = Irrotational { conf: config };
        let mat = Material::from_rgas_gamma(287.042, 1.4);

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

    #[test]
    fn test_symmetry_axis_point_1() {
        let config = GeneralConfig {
            axisym: AreaType::Axisymmetric,
            tol: Tolerance::new(1e-5, 1e-5),
            n_corr: 20,
        };

        let unitprocess = Irrotational { conf: config };
        let mat = Material::from_rgas_gamma(320.0, 1.2);

        // 构造两个输入点
        let p1 = MocPoint::from_compatible(
            0.079625,
            0.001290,
            2306.1,
            35.7,
            170250.,
            3000.0,
            0.32947,
            mat.clone(),
        );

        let velocity = 2332.5011075972866;
        let theta = 0.0_f64;
        let target = MocPoint::new(
            0.0832869605831241,
            0.,
            velocity * theta.cos(),
            velocity * theta.sin(),
            151235.51317795238,
            1583.1871310045397,
            0.2985071667131088,
            mat.clone(),
        );

        // 创建 CharLine
        let mut next_line = CharLine::new();
        next_line.push(p1.clone());

        let mut prev_line = CharLine::new();
        prev_line.push(p1.clone());

        let context = Context {
            prev: &prev_line,
            next: &next_line,
            idx_prev: 0,
            idx_next: 0,
        };

        let result_point = unitprocess
            .symmetry_axis_point(context)
            .expect("symmetry axis point should be Some");

        assert!(
            result_point.is_converged_with(&target, unitprocess.conf.tol),
            "Symmetry axis point did not converge to expected value:\nresult:{:15}\ntarget:{:15}\n  diff:{}",
            result_point,
            target,
            &result_point - &target
        );
    }

    #[test]
    fn test_inverse_wall_point_1() {
        let config = GeneralConfig {
            axisym: AreaType::Axisymmetric,
            tol: Tolerance::new(1e-5, 1e-5),
            n_corr: 20,
        };

        let unitprocess = Irrotational { conf: config };
        let mat = Material::from_rgas_gamma(320.0, 1.2);

        // 构造两个输入点
        let p1 = MocPoint::from_compatible(
            0.005085,
            0.026080,
            1577.5,
            702.3,
            1160400.,
            3000.0,
            1.6309,
            mat.clone(),
        );

        let p2 = MocPoint::from_compatible(
            0.005495,
            0.026020,
            1578.3,
            705.7,
            1154500.,
            3000.0,
            1.6240,
            mat.clone(),
        );

        let velocity = 1727.8264988128874;
        let theta = 0.4196581804279899_f64;
        let target = MocPoint::new(
            0.00529,
            0.02605,
            velocity * theta.cos(),
            velocity * theta.sin(),
            1157449.3598894787,
            2222.5561432291665,
            1.6274499981586739,
            mat.clone(),
        );

        // 创建 CharLine
        let mut next_line = CharLine::new();
        next_line.push(target.clone());

        let mut prev_line = CharLine::new();
        prev_line.push(p1.clone());
        prev_line.push(p2.clone());

        let context = Context {
            prev: &prev_line,
            next: &next_line,
            idx_prev: 0,
            idx_next: 0,
        };

        let result_point = unitprocess
            .inverse_wall_point(context, (target.x, target.y, target.flow_direction()))
            .expect("inverse wall point should be Some");

        assert!(
            result_point.is_converged_with(&target, unitprocess.conf.tol),
            "Inverse wall point did not converge to expected value:\nresult:{:15}\ntarget:{:15}\n  diff:{}",
            result_point,
            target,
            &result_point - &target
        );
    }

    #[test]
    fn test_exit_characteristics_point_1() {
        let config = GeneralConfig {
            axisym: AreaType::Axisymmetric,
            tol: Tolerance::new(1e-5, 1e-5),
            n_corr: 20,
        };

        let unitprocess = Irrotational { conf: config };
        let mat = Material::from_rgas_gamma(287.04, 1.4);

        // 构造两个输入点
        let mut velocity = 578.68844312083490;
        let mut theta = 0.0_f64;
        let p1 = MocPoint::from_compatible(
            1.8951087644550260,
            0.0000000000000000,
            velocity * theta.cos(),
            velocity * theta.sin(),
            175574.31375951759,
            300.0,
            4.5876520884696967,
            mat.clone(),
        );

        velocity = 578.42134504689045;
        theta = 0.00052976493774764544_f64;
        let p2 = MocPoint::from_compatible(
            1.8919379416868118,
            0.0013837169816496176,
            velocity * theta.cos(),
            velocity * theta.sin(),
            176286.08472429001,
            300.0,
            4.6009075183326198,
            mat.clone(),
        );

        let velocity = 578.6884431208349;
        let theta = 0.0_f64;
        let target = MocPoint::new(
            1.898286422924196,
            0.0013870801847318665,
            velocity * theta.cos(),
            velocity * theta.sin(),
            175575.22539912476,
            133.33317695810655,
            4.587658492100335,
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
            .exit_characteristics_point(context, Box::new(|_, p: &MocPoint| (p.velocity(), 0.)))
            .expect("exit characteristics point should be Some");

        assert!(
            result_point.is_converged_with(&target, unitprocess.conf.tol),
            "Exit characteristics point did not converge to expected value:\nresult:{:15}\ntarget:{:15}\n  diff:{}",
            result_point,
            target,
            &result_point - &target
        );
    }

    #[test]
    fn test_exit_characteristics_point_fixed_dist_1() {
        let config = GeneralConfig {
            axisym: AreaType::Axisymmetric,
            tol: Tolerance::new(1e-5, 1e-5),
            n_corr: 20,
        };

        let unitprocess = Irrotational { conf: config };
        let mat = Material::from_rgas_gamma(287.04, 1.4);

        // 构造两个输入点
        let velocity = 578.42134504689045;
        let theta = 0.00052976493774764544_f64;
        let p1 = MocPoint::from_compatible(
            1.8919379416868118,
            0.0013837169816496176,
            velocity * theta.cos(),
            velocity * theta.sin(),
            176286.08472429001,
            300.,
            4.6009075183326198,
            mat.clone(),
        );

        let velocity = 578.68844312083490;
        let theta = 0.0_f64;
        let p2 = MocPoint::from_compatible(
            1.8951087644550260,
            0.0000000000000000,
            velocity * theta.cos(),
            velocity * theta.sin(),
            175574.31375951759,
            300.0,
            4.5876520884696967,
            mat.clone(),
        );

        let velocity = 578.6884431208349;
        let theta = 0.0_f64;
        let target = MocPoint::new(
            1.898279448151052,
            0.0013840356254392625,
            velocity * theta.cos(),
            velocity * theta.sin(),
            175575.22539912476,
            133.33317695810655,
            4.587658492100335,
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
            .exit_characteristics_point_fixed_dist(
                context,
                Box::new(|_, p: &MocPoint| (p.velocity(), 0.0)),
                p1.distance_to(&p2),
            )
            .expect("exit characteristics point should be Some");

        assert!(
            result_point.is_converged_with(&target, unitprocess.conf.tol),
            "Exit characteristics point did not converge to expected value:\nresult:{:15}\ntarget:{:15}\n  diff:{}",
            result_point,
            target,
            &result_point - &target
        );
    }

    #[test]
    fn test_transition_interior_point_1() {
        let config = GeneralConfig {
            axisym: AreaType::Axisymmetric,
            tol: Tolerance::new(1e-5, 1e-5),
            n_corr: 20,
        };

        let unitprocess = Irrotational { conf: config };
        let mat = Material::from_rgas_gamma(287.04, 1.4);

        // 构造两个输入点
        let velocity = 578.68844312083490_f64;
        let theta = 0.0_f64;
        let p1 = MocPoint::from_compatible(
            1.8982869690824169,
            0.0013870804740672891,
            velocity * theta.cos(),
            velocity * theta.sin(),
            175574.31375951762,
            300.0,
            4.5876310941627727,
            mat.clone(),
        );

        let velocity = 578.42134504689045_f64;
        let theta = 0.00052976493774764544_f64;
        let p2 = MocPoint::from_compatible(
            1.8919379416868118,
            0.0013837169816496176,
            velocity * theta.cos(),
            velocity * theta.sin(),
            176286.08472429001,
            300.0,
            4.6009075183326198,
            mat.clone(),
        );

        let velocity = 578.5103632704744_f64;
        let theta = 0.0005287809353698835_f64;
        let target = MocPoint::new(
            1.8951118904357467,
            0.002772289355895323,
            velocity * theta.cos(),
            velocity * theta.sin(),
            176048.36890983424,
            133.43573797015043,
            4.596475229646061,
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
            .transition_interior_point(context)
            .expect("transition interior point should be Some");

        assert!(
            result_point.is_converged_with(&target, unitprocess.conf.tol),
            "transition interior point did not converge to expected value:\nresult:{:15}\ntarget:{:15}\n  diff:{}",
            result_point,
            target,
            &result_point - &target
        );
    }

    #[test]
    fn test_transition_wall_point_1() {
        let config = GeneralConfig {
            axisym: AreaType::Axisymmetric,
            tol: Tolerance::new(1e-5, 1e-5),
            n_corr: 20,
        };

        let unitprocess = Irrotational { conf: config };
        let mat = Material::from_rgas_gamma(287.04, 1.4);

        // 构造两个输入点
        let velocity = 421.13870346605336_f64;
        let theta = 0.13644122164345995_f64;
        let p1 = MocPoint::from_compatible(
            0.0071288894325802748,
            0.99499809824891272,
            velocity * theta.cos(),
            velocity * theta.sin(),
            886014.89349445840,
            300.0,
            14.578529717965795,
            mat.clone(),
        );

        let velocity = 420.96812614461288_f64;
        let theta = 0.13493747221790942_f64;
        let p2 = MocPoint::from_compatible(
            0.0000000000000000,
            1.0000000000000000,
            velocity * theta.cos(),
            velocity * theta.sin(),
            887060.13126201590,
            300.0,
            14.590840930401765,
            mat.clone(),
        );

        let velocity = 420.91142292916430_f64;
        let theta = 0.13589472516792001_f64;
        let p3 = MocPoint::from_compatible(
            0.0070286901822888574,
            0.99487139839190986,
            velocity * theta.cos(),
            velocity * theta.sin(),
            887410.75364098849,
            300.0,
            14.594930905654611,
            mat.clone(),
        );

        let (velocity, theta): (f64, f64) = (420.6220130296241, 0.13621169438131384);
        let target = MocPoint::new(
            0.00021939206260824715,
            1.0000299275858477,
            velocity * theta.cos(),
            velocity * theta.sin(),
            889188.1360271955,
            211.94712641090678,
            14.61581981911365,
            mat.clone(),
        );

        // 创建 CharLine
        let mut next_line = CharLine::new();
        next_line.push(p1.clone());

        let mut prev_line = CharLine::new();
        prev_line.push(p2.clone());
        prev_line.push(p3.clone());

        let context = Context {
            prev: &prev_line,
            next: &next_line,
            idx_prev: 0,
            idx_next: 0,
        };

        let result_point = unitprocess
            .transition_wall_point(context)
            .expect("transition wall point should be Some");

        assert!(
            result_point.is_converged_with(&target, unitprocess.conf.tol),
            "transition wall point did not converge to expected value:\nresult:{:15}\ntarget:{:15}\n  diff:{}",
            result_point,
            target,
            &result_point - &target
        );
    }

    #[test]
    fn test_transition_wall_point_2() {
        let config = GeneralConfig {
            axisym: AreaType::Axisymmetric,
            tol: Tolerance::new(1e-5, 1e-5),
            n_corr: 20,
        };

        let unitprocess = Irrotational { conf: config };
        let mat = Material::from_rgas_gamma(287.04, 1.2874718741354083);

        // 构造两个输入点
        let (velocity, theta): (f64, f64) = (2035.7726107195051, 0.31325662671349330);
        let p1 = MocPoint::from_compatible(
            0.11604615315469557,
            0.43042968648297669,
            velocity * theta.cos(),
            velocity * theta.sin(),
            42576.927088002609,
            2868.0,
            0.11810927136810002,
            mat.clone(),
        );

        let (velocity, theta): (f64, f64) = (2035.9849779642398, 0.31583904158230208);
        let p2 = MocPoint::from_compatible(
            0.10087027099609050,
            0.42723582495402546,
            velocity * theta.cos(),
            velocity * theta.sin(),
            42525.885489298969,
            2868.0,
            0.11799928133601660,
            mat.clone(),
        );

        let (velocity, theta): (f64, f64) = (2032.8683375746696, 0.31146882264229447);
        let p3 = MocPoint::from_compatible(
            0.10095815521696822,
            0.42723335376744537,
            velocity * theta.cos(),
            velocity * theta.sin(),
            43279.197749120111,
            2868.0,
            0.11961963449605677,
            mat.clone(),
        );

        let (velocity, theta): (f64, f64) = (2035.0914547801467, 0.3138191082330746);
        let target = MocPoint::new(
            0.11109947528489121,
            0.4305670752761383,
            velocity * theta.cos(),
            velocity * theta.sin(),
            42740.89366076534,
            1257.157592612592,
            0.11846240650806635,
            mat.clone(),
        );

        // 创建 CharLine
        let mut next_line = CharLine::new();
        next_line.push(p1.clone());

        let mut prev_line = CharLine::new();
        prev_line.push(p2.clone());
        prev_line.push(p3.clone());

        let context = Context {
            prev: &prev_line,
            next: &next_line,
            idx_prev: 0,
            idx_next: 0,
        };

        let result_point = unitprocess
            .transition_wall_point(context)
            .expect("transition wall point should be Some");

        assert!(
            result_point.is_converged_with(&target, unitprocess.conf.tol),
            "transition wall point did not converge to expected value:\nresult:{:15}\ntarget:{:15}\n  diff:{}",
            result_point,
            target,
            &result_point - &target
        );
    }

    #[test]
    fn test_transition_wall_point_3() {
        let config = GeneralConfig {
            axisym: AreaType::Axisymmetric,
            tol: Tolerance::new(1e-5, 1e-5),
            n_corr: 20,
        };

        let unitprocess = Irrotational { conf: config };
        let mat = Material::from_rgas_gamma(320., 1.4);

        // 构造两个输入点
        let (velocity, theta): (f64, f64) = (1789.6295849830256, 0.32746227178513421);
        let p1 = MocPoint::from_compatible(
            0.35843974472954165,
            1.1199968743597424,
            velocity * theta.cos(),
            velocity * theta.sin(),
            9898.1499933631203,
            2000.0,
            0.054248073729845178,
            mat.clone(),
        );

        let (velocity, theta): (f64, f64) = (1794.7436828628317, 0.34103154223032950);
        let p2 = MocPoint::from_compatible(
            0.31381485362272077,
            1.0981396888218800,
            velocity * theta.cos(),
            velocity * theta.sin(),
            9409.7996711261931,
            2000.0,
            0.052322544523846583,
            mat.clone(),
        );

        let (velocity, theta): (f64, f64) = (1789.4319903359690, 0.33055297823245189);
        let p3 = MocPoint::from_compatible(
            0.31495905666087648,
            1.0981979811985894,
            velocity * theta.cos(),
            velocity * theta.sin(),
            9917.3456869649326,
            2000.0,
            0.054323198622279717,
            mat.clone(),
        );

        let (velocity, theta): (f64, f64) = (1794.1731949180319, 0.33416148703456666);
        let target = MocPoint::new(
            0.378701656042577,
            1.120917525312198,
            velocity * theta.cos(),
            velocity * theta.sin(),
            9463.47216247024,
            562.9207797489387,
            0.05253554418231826,
            mat.clone(),
        );

        // 创建 CharLine
        let mut next_line = CharLine::new();
        next_line.push(p1.clone());

        let mut prev_line = CharLine::new();
        prev_line.push(p2.clone());
        prev_line.push(p3.clone());

        let context = Context {
            prev: &prev_line,
            next: &next_line,
            idx_prev: 0,
            idx_next: 0,
        };

        let result_point = unitprocess
            .transition_wall_point(context)
            .expect("transition wall point should be Some");

        assert!(
            result_point.is_converged_with(&target, unitprocess.conf.tol),
            "transition wall point did not converge to expected value:\nresult:{:15}\ntarget:{:15}\n  diff:{}",
            result_point,
            target,
            &result_point - &target
        );
    }

    #[test]
    fn test_last_point_1() {
        let config = GeneralConfig {
            axisym: AreaType::Axisymmetric,
            tol: Tolerance::new(1e-5, 1e-5),
            n_corr: 20,
        };

        let unitprocess = Irrotational { conf: config };
        let mat = Material::from_rgas_gamma(287.04, 1.4);

        // 构造两个输入点
        let (velocity, theta): (f64, f64) = (578.67913564023843, 0.0);
        let p2 = MocPoint::from_compatible(
            5.6047812037733848,
            1.6203058834240001,
            velocity * theta.cos(),
            velocity * theta.sin(),
            175606.89227499434,
            300.0,
            4.5882134340057954,
            mat.clone(),
        );

        let (velocity, theta): (f64, f64) = (578.68844312083490, 0.0);
        let p3 = MocPoint::from_compatible(
            5.6062409066524017,
            1.6196688176592677,
            velocity * theta.cos(),
            velocity * theta.sin(),
            175574.31375951803,
            300.0,
            4.5875638563519372,
            mat.clone(),
        );

        let (velocity, theta): (f64, f64) = (578.6884431208349, 0.);
        let target = MocPoint::new(
            5.6077005803442805,
            1.6203058783093023,
            velocity * theta.cos(),
            velocity * theta.sin(),
            175606.89227499434,
            133.33853817022347,
            4.588213434005795,
            mat.clone(),
        );

        // 创建 CharLine
        let next_line = CharLine::new();

        let mut prev_line = CharLine::new();
        prev_line.push(p2.clone());
        prev_line.push(p3.clone());

        let context = Context {
            prev: &prev_line,
            next: &next_line,
            idx_prev: 0,
            idx_next: 0,
        };

        let result_point = unitprocess
            .last_point(
                context,
                Box::new(|_, p: &MocPoint| (p.velocity(), 0.0)),
                17.215908632646901,
            )
            .expect("last point should be Some");

        assert!(
            result_point.is_converged_with(&target, unitprocess.conf.tol),
            "last point did not converge to expected value:\nresult:{:15}\ntarget:{:15}\n  diff:{}",
            result_point,
            target,
            &result_point - &target
        );
    }
}
