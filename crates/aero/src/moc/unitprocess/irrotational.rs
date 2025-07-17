use crate::moc::MocPoint;
use crate::moc::unitprocess::GeneralConfig;
use crate::moc::unitprocess::UnitProcess;

/// 无旋特征线法的基本计算过程
pub struct Irrotational {
    conf: GeneralConfig,
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
    fn cal_lqrs_left(&self, point: MocPoint) -> (f64, f64, f64, f64) {
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
    fn cal_lqrs_right(&self, point: MocPoint) -> (f64, f64, f64, f64) {
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
}

impl UnitProcess for Irrotational {
    fn interior_point(&self, context: super::Context) -> Option<MocPoint> {
        todo!()
    }

    fn symmetry_axis_point(&self, context: super::Context) -> Option<MocPoint> {
        todo!()
    }

    fn inverse_wall_point(&self, context: super::Context) -> Option<MocPoint> {
        todo!()
    }

    fn exit_characteristics_point(
        &self,
        context: super::Context,
        cal_u_v: super::ExitLineFunc,
    ) -> Option<MocPoint> {
        todo!()
    }

    fn exit_characteristics_point_fixed_dist(
        &self,
        context: super::Context,
        cal_u_v: super::ExitLineFunc,
        dist: f64,
    ) -> Option<MocPoint> {
        todo!()
    }

    fn transition_interior_point(&self, context: super::Context) -> Option<MocPoint> {
        todo!()
    }

    fn transition_wall_point(&self, context: super::Context) -> Option<MocPoint> {
        todo!()
    }

    fn last_point(
        &self,
        context: super::Context,
        cal_u_v: super::ExitLineFunc,
    ) -> Option<MocPoint> {
        todo!()
    }
}
