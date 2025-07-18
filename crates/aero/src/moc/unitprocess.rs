mod irrotational;
mod rotational;

pub use irrotational::Irrotational;
use math::Tolerance;
// pub use rotational::Rotational;

use crate::moc::{CharLine, MocPoint};

#[derive(Clone, Copy)]
pub struct Context<'a> {
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
    pub tol: Tolerance,
    /// 最大预估校正次数
    pub n_corr: u16,
}

/// 出口最后一条特征线边界条件
///
/// # 输入参数
/// - `y:f64`, `p:MocPoint`: 待求点y轴坐标, 参考点气流参数p
/// # 返回值
/// - `u:f64`, `v:f64`: 待求点气流在x和y方向的速度u, v
pub type ExitLineFunc = Box<dyn Fn(f64, MocPoint) -> (f64, f64)>;

/// 特征线法基本计算过程
pub trait UnitProcess {
    /// 求解流场内点
    ///
    /// # 参数
    /// - `&self`: 包含计算过程中所需的控制参数
    /// - `context`: 参数上下文，要求p1发出的右行特征线和p2发出的左行特征线交于点pr
    ///     - p1为context.next\[idx_next\]
    ///     - p2为context.prev\[idx_prev\]
    ///     - pr为待求点
    fn interior_point(&self, context: Context) -> Option<MocPoint>;

    /// 求解轴点 假设轴线方向为x轴正方向
    ///
    /// # 参数
    /// - `&self`: 包含计算过程中所需的控制参数
    /// - `context`: 参数上下文，要求p1发出的右行特征线和p2发出的流线交于点pr，且p2在轴线上
    ///     - p1为context.next\[idx_next\]
    ///     - p2为context.prev\[idx_prev\]
    ///     - pr为待求点
    fn symmetry_axis_point(&self, context: Context) -> Option<MocPoint>;

    /// 求解逆推壁面点
    ///
    /// # 参数
    /// - `&self`: 包含计算过程中所需的控制参数
    /// - `context`: 参数上下文，要求p1和p2之间某一点发出的左行特征线通过点pr，且点p1和pr为壁面点
    ///     - p1为context.prev\[idx_prev\]
    ///     - p2为context.prev\[idx_prev+1\]
    ///     - pr为待求点，context.next\[idx_next\]应存有点pr的坐标和倾角(x,y,θ)
    fn inverse_wall_point(&self, context: Context) -> Option<MocPoint>;

    /// 求解出口左行特征线上的点
    ///
    /// # 参数
    /// - `&self`: 包含计算过程中所需的控制参数
    /// - `context`: 参数上下文，要求p1发出的右行特征线和p2发出的左行特征线交于点pr，p2和pr均在出口左行特征线上
    ///     - p1为context.next\[idx_next\]
    ///     - p2为context.prev\[idx_prev\]
    ///     - pr为待求点
    /// - `cal_u_v`: 出口左行特征线边界条件
    fn exit_characteristics_point(
        &self,
        context: Context,
        cal_u_v: ExitLineFunc,
    ) -> Option<MocPoint>;

    /// 求解出口左行特征线上的点
    ///
    /// # 参数
    /// - `&self`: 包含计算过程中所需的控制参数
    /// - `context`: 参数上下文，要求p1发出的右行特征线和p2发出的左行特征线交于点pr，p2和pr均在出口左行特征线上
    ///     - p1为context.next\[idx_next\]
    ///     - p2为context.prev\[idx_prev\]
    ///     - pr为待求点
    /// - `cal_u_v`: 出口左行特征线边界条件
    /// - `dist`: 出口特征线上两点p2和pr的距离
    fn exit_characteristics_point_fixed_dist(
        &self,
        context: Context,
        cal_u_v: ExitLineFunc,
        dist: f64,
    ) -> Option<MocPoint>;

    /// 求解转向区内点
    ///
    /// # 参数
    /// - `&self`: 包含计算过程中所需的控制参数
    /// - `context`: 参数上下文，要求p1在pr发出的右行特征线上，p2发出的左行特征线过点pr，p3发出的左行特征线过点p1
    ///     - p1为context.next\[idx_next\]
    ///     - p2为context.prev\[idx_prev\]
    ///     - p3为context.prev\[idx_prev+1\]
    ///     - pr为待求点
    fn transition_interior_point(&self, context: Context) -> Option<MocPoint>;

    /// 求解转向区壁面点
    ///
    /// # 参数
    /// - `&self`: 包含计算过程中所需的控制参数
    /// - `context`: 参数上下文，要求p1在pr发出的右行特征线上，p2和p3上一点发出的左行特征线过点pr
    ///     - p1为context.next\[idx_next\]
    ///     - p2为context.prev\[idx_prev\]
    ///     - p3为context.prev\[idx_prev+1\]
    ///     - pr为待求点
    fn transition_wall_point(&self, context: Context) -> Option<MocPoint>;

    /// 求解最后一点，该点应为上壁面出口点，且为出口左行特征线上最后一点
    ///
    /// # 参数
    /// - `&self`: 包含计算过程中所需的控制参数
    /// - `context`: 参数上下文，要求p3在p2发出的右行特征线上，p3发出的左行特征线过点pr，p2发出的流线过点pr
    ///     - p2为context.prev\[idx_prev\]，为壁面点
    ///     - p3为context.prev\[idx_prev+1\]，为出口特征线点
    ///     - pr为待求点
    /// - `cal_u_v`: 出口左行特征线边界条件
    fn last_point(&self, context: Context, cal_u_v: ExitLineFunc) -> Option<MocPoint>;

    // 下面这三个暂时用不到
    // fn direct_walld_point(&self, context: Context) -> Option<MocPoint>;
    // fn free_pressure_boundary_point(&self, context: Context) -> Option<MocPoint>;
    // fn unknown_wall_point(&self, context: Context, mfr_need: f64) -> Option<MocPoint>;
}
