use crate::moc::unitprocess::UnitProcess;
use crate::moc::unitprocess::UnitprocessConfig;

// /// 有旋特征线法的基本计算过程
pub struct Rotational {
    pub conf: UnitprocessConfig,
}

#[allow(unused)]
impl UnitProcess for Rotational {
    fn interior_point(&self, context: super::Context) -> Option<crate::moc::MocPoint> {
        todo!()
    }

    fn symmetry_axis_point(&self, context: super::Context) -> Option<crate::moc::MocPoint> {
        todo!()
    }

    fn inverse_wall_point(
        &self,
        context: super::Context,
        wall_info: (f64, f64, f64),
    ) -> Option<crate::moc::MocPoint> {
        todo!()
    }

    fn exit_characteristics_point(
        &self,
        context: super::Context,
        cal_u_v: super::ExitLineFunc,
    ) -> Option<crate::moc::MocPoint> {
        todo!()
    }

    fn exit_characteristics_point_fixed_dist(
        &self,
        context: super::Context,
        cal_u_v: super::ExitLineFunc,
        dist: f64,
    ) -> Option<crate::moc::MocPoint> {
        todo!()
    }

    fn transition_interior_point(&self, context: super::Context) -> Option<crate::moc::MocPoint> {
        todo!()
    }

    fn transition_wall_point(&self, context: super::Context) -> Option<crate::moc::MocPoint> {
        todo!()
    }

    fn last_point(
        &self,
        context: super::Context,
        cal_u_v: super::ExitLineFunc,
        mfr_need: f64,
    ) -> Option<crate::moc::MocPoint> {
        todo!()
    }
}
