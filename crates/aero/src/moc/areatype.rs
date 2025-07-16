/// 描述面积微元的建模方式，用于流场积分计算。
#[derive(Debug, Clone, Copy)]
pub enum AreaType {
    /// 平面模型，使用指定的固定深度 `d`（单位：m）
    Planar(f64),

    /// 轴对称模型，根据 y 坐标自动计算环形面积（单位：m）
    Axisymmetric,
}
