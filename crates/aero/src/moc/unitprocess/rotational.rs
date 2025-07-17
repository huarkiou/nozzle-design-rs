use crate::moc::unitprocess::GeneralConfig;
use crate::moc::unitprocess::UnitProcess;

/// 有旋特征线法的基本计算过程
pub struct Rotational {
    conf: GeneralConfig,
}

impl UnitProcess for Rotational {
    fn interior_point(&self, context: super::Context) -> Option<crate::moc::MocPoint> {
        todo!()
    }

    fn symmetry_axis_point(&self, context: super::Context) -> Option<crate::moc::MocPoint> {
        todo!()
    }

    fn inverse_wall_point(&self, context: super::Context) -> Option<crate::moc::MocPoint> {
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
    ) -> Option<crate::moc::MocPoint> {
        todo!()
    }
}
