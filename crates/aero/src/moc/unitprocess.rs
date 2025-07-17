use crate::moc::{CharLine, MocPoint};

#[derive(Clone, Copy)]
pub struct UPContext<'a> {
    /// 前一条特征线，不单调总体上从左上到右下
    pub prev: &'a CharLine,
    /// 后一条特征线，不单调总体上从左上到右下
    pub next: &'a CharLine,
    /// 前一条特征线初始点下标
    pub idx_prev: usize,
    /// 后一条特征线初始点下标
    pub idx_next: usize,
}

#[derive(Clone, Copy)]
pub struct GeneralConfig {
    /// 轴对称/二维平面
    pub axisym: bool,
    /// 容差
    pub tol: f64,
    /// 最大预估校正次数
    pub n_corr: i64,
}

/// 出口最后一条特征线边界条件
///
/// # 输入参数
/// - `y:f64`, `p:MocPoint`: 待求点y轴坐标, 参考点气流参数p
/// # 返回值
/// - `u:f64`, `v:f64`: 待求点气流在x和y方向的速度u, v
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
