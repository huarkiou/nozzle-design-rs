use nalgebra::Vector3;
use std::ops::{Deref, DerefMut};

use crate::moc::{AreaType, MocPoint};

/// 一条特征线上所有MocPoint的集合
#[derive(Clone)]
pub struct CharLine {
    data: Vec<MocPoint>,
}

impl Deref for CharLine {
    type Target = Vec<MocPoint>;

    fn deref(&self) -> &Self::Target {
        &self.data
    }
}

impl DerefMut for CharLine {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.data
    }
}

impl CharLine {
    pub fn new() -> Self {
        Self { data: Vec::new() }
    }

    pub fn with_capacity(n: usize) -> Self {
        Self {
            data: Vec::with_capacity(n),
        }
    }
}

impl CharLine {
    /// 计算通过给定特征线上的质量流量。
    ///
    /// # 参数
    /// - `line`: `MocPoint` 的切片，表示特征线上的点
    /// - `area_type`: 面积计算方式，详见 [`AreaType`]：
    ///   - [`AreaType::Planar(d)`]：使用指定的固定深度 `d`
    ///   - [`AreaType::Axisymmetric`]：使用 y 坐标计算环形面积（轴对称模型）
    ///
    /// # 返回值
    /// 返回质量流量（单位：kg/s），若输入点数不足 2 个则返回 `0.0`
    ///
    /// # 公式说明
    /// 质量流量计算公式为：
    /// ```text
    /// mfr = ∫ ρ * (V · dA)
    /// ```
    /// 其中 ρ 为密度，V 为速度向量，dA 为面积微元向量。
    ///
    pub fn mass_flow_rate(line: &[MocPoint], area_type: AreaType) -> f64 {
        let mut mfr = 0.0;
        for window in line.windows(2) {
            if let [p1, p2] = window {
                // 中心差分计算边长向量
                let dx = p1.x - p2.x;
                let dy = p1.y - p2.y;
                let dl = Vector3::new(dx, dy, 0.0);
                // dA = dl × n （n 是 z 方向单位向量）
                let da = dl.cross(&Vector3::new(0.0, 0.0, 1.0));
                // 如果是轴对称（rot==true），使用环形面积修正
                let current_depth = match area_type {
                    AreaType::Planar(d) => d,
                    AreaType::Axisymmetric => 2.0 * std::f64::consts::PI * (p1.y + p2.y) / 2.0,
                };
                // 计算平均速度向量
                let v_avg = Vector3::new((p1.u + p2.u) / 2.0, (p1.v + p2.v) / 2.0, 0.0);
                // 计算面积向量
                let da_scaled = da * current_depth;
                // 平均密度
                let rho_avg = (p1.rho + p2.rho) / 2.0;
                // 质量流量增量
                mfr += rho_avg * da_scaled.dot(&v_avg);
            }
        }
        mfr
    }

    /// 计算指定线上在x方向的总推力，包括压差力和动量力。
    ///
    /// # 参数
    /// - `line`: `MocPoint` 的切片，表示一段特征线
    /// - `p_ambient`: 环境压力（单位：Pa）
    /// - `area_type`: 面积计算方式，详见 [`AreaType`]：
    ///   - [`AreaType::Planar(d)`]：使用指定的固定深度 `d`
    ///   - [`AreaType::Axisymmetric`]：使用 y 坐标计算环形面积（轴对称模型）
    ///
    /// # 返回值
    /// 返回总推力（单位：N），若输入点数不足 2 个则返回 `0.0`
    ///
    /// # 推力组成
    /// 推力由两部分组成：
    /// 1. **压差力**：`(p_avg - p_ambient) * dA · direction`
    /// 2. **动量力**：`rho_avg * (dA · V_avg) * (V_avg · direction)`
    ///
    pub fn thrust(line: &[MocPoint], p_ambient: f64, area_type: AreaType) -> f64 {
        let mut ret = 0.0;
        let direction = Vector3::new(1.0, 0.0, 0.0).normalize();
        for window in line.windows(2) {
            if let [p1, p2] = window {
                // 环形修正
                let current_depth = match area_type {
                    AreaType::Planar(d) => d,
                    AreaType::Axisymmetric => 2.0 * std::f64::consts::PI * (p1.y + p2.y) / 2.0,
                };
                // 面积微元 dA = dl × n * depth
                let dl = Vector3::new(p2.x - p1.x, p2.y - p1.y, 0.0);
                let da = dl.cross(&Vector3::new(0.0, 0.0, 1.0)) * current_depth;
                // 平均压力
                let p_avg = (p1.p + p2.p) / 2.0;
                // 平均速度
                let v_avg = Vector3::new((p1.u + p2.u) / 2.0, (p1.v + p2.v) / 2.0, 0.0);
                // 平均密度
                let rho_avg = (p1.rho + p2.rho) / 2.0;
                // 推力分量
                let pressure_force = (p_avg - p_ambient) * da.dot(&direction);
                let momentum_force = rho_avg * da.dot(&v_avg) * v_avg.dot(&direction);
                ret += pressure_force + momentum_force;
            }
        }
        ret
    }

    /// 计算壁面上由于压力分布在指定方向上的合力。
    ///
    /// # 参数
    /// - `line`: `MocPoint` 的切片，表示一系列壁面点（用于计算压力合力）
    /// - `direction`: 方向向量 `[x, y, z]`，表示受力方向。即使不是单位向量也会被自动归一化
    /// - `p_ambient`: 环境压力（单位：Pa）
    /// - `area_type`: 面积计算方式，详见 [`AreaType`]：
    ///   - [`AreaType::Planar(d)`]：使用指定的固定深度 `d`
    ///   - [`AreaType::Axisymmetric`]：使用 y 坐标计算环形面积（轴对称模型）
    ///
    /// # 返回值
    /// 返回在指定方向上的净压力合力（单位：N），若输入点数不足 2 个则返回 `0.0`
    ///
    /// # 公式说明
    /// 压力合力计算公式为：
    /// ```text
    /// F = ∫ (p_avg - p_ambient) * (dA · direction)
    /// ```
    ///
    pub fn force(
        line: &[MocPoint],
        direction: [f64; 3],
        p_ambient: f64,
        area_type: AreaType,
    ) -> f64 {
        let mut ret = 0.0;
        let direc = Vector3::new(direction[0], direction[1], direction[2]).normalize();
        for window in line.windows(2) {
            if let [p1, p2] = window {
                let width = match area_type {
                    AreaType::Planar(d) => d,
                    AreaType::Axisymmetric => 2.0 * std::f64::consts::PI * (p1.y + p2.y) / 2.0,
                };
                // 边向量
                let dl = Vector3::new(p1.x - p2.x, p1.y - p2.y, 0.0);
                // 面积微元 dA = dl × n * width
                let da = dl.cross(&Vector3::new(0.0, 0.0, 1.0)) * width;
                // 平均压力
                let p_avg = (p1.p + p2.p) / 2.0;
                // 压差力
                let force_component = (p_avg - p_ambient) * da.dot(&direc);
                ret += force_component;
            }
        }
        ret
    }
}

#[cfg(test)]
mod test {
    use crate::Material;

    use super::*;

    #[test]
    fn test_mass_flow_rate_empty() {
        let line = &CharLine::new();
        assert_eq!(
            CharLine::mass_flow_rate(&line[0..], AreaType::Axisymmetric),
            0.
        );
    }

    #[test]
    fn test_mass_flow_rate_one() {
        let single_point = vec![MocPoint {
            x: 0.0,
            y: 0.0,
            u: 100.0,
            v: 0.0,
            p: 101325.0,
            t: 298.0,
            rho: 1.225,
            mat: Material::from_rgas_gamma(287.0, 1.4),
        }];
        let mfr_single = CharLine::mass_flow_rate(&single_point, AreaType::Axisymmetric);
        assert_eq!(mfr_single, 0.0);
    }

    #[test]
    fn test_mass_flow_rate_two() {
        let mut line = vec![
            MocPoint {
                x: 0.0,
                y: 0.0,
                u: 100.0,
                v: 100.0,
                p: 101325.0,
                t: 298.0,
                rho: 1.225,
                mat: Material::from_rgas_gamma(287.0, 1.4),
            },
            MocPoint {
                x: 1.0,
                y: 0.0,
                u: 100.0,
                v: 100.0,
                p: 101325.0,
                t: 298.0,
                rho: 1.225,
                mat: Material::from_rgas_gamma(287.0, 1.4),
            },
        ];

        assert_eq!(
            CharLine::mass_flow_rate(&line, AreaType::Planar(1.0)),
            122.50000000000001
        );
        line.reverse();
        assert_eq!(
            CharLine::mass_flow_rate(&line, AreaType::Planar(1.0)),
            -122.50000000000001
        );
    }

    #[test]
    fn test_thrust_and_force() {
        let line = vec![
            MocPoint {
                x: 0.0,
                y: 0.0,
                u: 100.0,
                v: 100.0,
                p: 101325.0,
                t: 298.0,
                rho: 1.225,
                mat: Material::from_rgas_gamma(287.0, 1.4),
            },
            MocPoint {
                x: 1.0,
                y: 0.0,
                u: 100.0,
                v: 100.0,
                p: 101325.0,
                t: 298.0,
                rho: 1.225,
                mat: Material::from_rgas_gamma(287.0, 1.4),
            },
        ];

        let ambient = 100000.0;

        assert_eq!(
            CharLine::thrust(&line, ambient, AreaType::Planar(1.0)),
            -12250.000000000002
        );

        assert_eq!(
            CharLine::force(&line, [1.0, 0.0, 0.0], ambient, AreaType::Planar(1.0)),
            0.0
        );

        assert_eq!(
            CharLine::force(&line, [0.0, 1.0, 0.0], ambient, AreaType::Planar(1.0)),
            1325.0
        );

        assert_eq!(
            CharLine::force(&line, [0.0, -1.0, 0.0], ambient, AreaType::Planar(1.0)),
            -1325.0
        );

        assert_eq!(
            CharLine::force(&line, [1.0, 0.0, 0.0], ambient, AreaType::Planar(1.0)),
            0.
        );

        assert_eq!(
            CharLine::force(&line, [1.0, -1.0, 0.0], ambient, AreaType::Planar(1.0)),
            -936.9164850721754
        );

        assert_eq!(
            CharLine::force(&line, [1.0, 1.0, 0.0], ambient, AreaType::Planar(1.0)),
            936.9164850721754
        );

        assert_eq!(
            CharLine::force(&line, [-1.0, -1.0, 0.0], ambient, AreaType::Planar(1.0)),
            -936.9164850721754
        );
    }
}
