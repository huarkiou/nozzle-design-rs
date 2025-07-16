use crate::moc::{CharLine, MocPoint};

#[derive(Clone, Copy)]
pub struct UPContext<'a> {
    pub prev: &'a CharLine,
    pub next: &'a CharLine,
    pub idx_prev: usize,
    pub idx_next: usize,
}

#[derive(Clone, Copy)]
pub struct GeneralConfig {
    pub axisym: bool,
    pub tol: f64,
    pub n_corr: i64,
}

pub type ExitLineFunc = Box<dyn Fn(f64, MocPoint) -> (f64, f64)>;

pub trait UnitProcess {
    fn interior_point(conf: GeneralConfig, context: UPContext) -> Option<MocPoint>;

    fn symmetry_axis_point(conf: GeneralConfig, context: UPContext) -> Option<MocPoint>;

    fn direct_walld_point(conf: GeneralConfig, context: UPContext) -> Option<MocPoint>;

    fn inverse_wall_point(conf: GeneralConfig, context: UPContext) -> Option<MocPoint>;

    fn free_pressure_boundary_point(conf: GeneralConfig, context: UPContext) -> Option<MocPoint>;

    fn exit_characteristics_point(
        conf: GeneralConfig,
        context: UPContext,
        cal_u_v: ExitLineFunc,
    ) -> Option<MocPoint>;

    fn exit_characteristics_point_fixed_dist(
        conf: GeneralConfig,
        context: UPContext,
        cal_u_v: ExitLineFunc,
    ) -> Option<MocPoint>;

    fn transition_interior_point(conf: GeneralConfig, context: UPContext) -> Option<MocPoint>;

    fn unknown_wall_point(
        conf: GeneralConfig,
        context: UPContext,
        mfr_need: f64,
    ) -> Option<MocPoint>;

    fn transition_wall_point(conf: GeneralConfig, context: UPContext) -> Option<MocPoint>;

    fn last_point(
        conf: GeneralConfig,
        context: UPContext,
        mfr_need: f64,
        cal_u_v: ExitLineFunc,
    ) -> Option<MocPoint>;
}
