use math::Tolerance;

use crate::moc::AreaType;

#[derive(Clone, Copy)]
pub struct UnitprocessConfig {
    /// 轴对称/二维平面
    pub axisym: AreaType,
    /// 容差
    pub tol: Tolerance,
    /// 最大预估校正次数
    pub n_corr: u16,
}

impl UnitprocessConfig {
    pub fn is_axisymmetric(&self) -> bool {
        match self.axisym {
            AreaType::Axisymmetric => true,
            _ => false,
        }
    }
}
