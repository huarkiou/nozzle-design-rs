use crate::cp::CpSegment;
use crate::Cp;
use math::quadrature;
use serde::{Deserialize, Deserializer, Serialize};

// MaterialProperty 结构体定义
#[derive(Clone, Debug, Serialize)]
pub struct Material {
    /// 摩尔质量 (kg/kmol)
    /// - **注意**：非国际单位制，1 kmol = 1000 mol
    molecular_weight: f64,

    /// 定压比热容 Cp(J/(kg·K))，输入温度 T(K)
    cp: Cp,
}

// ── TOML 反序列化辅助结构体 ─────────────────────────────────
// 用于同时支持两种格式:
//   1) cp = 1004.675              → 常数比热容
//   2) cp_segments = [...]        → 分段多项式变比热容
//   3) cp = nan                   → 向后兼容，触发 NASA 9 空气模型（在 config.rs 中处理）

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct MaterialHelper {
    molecular_weight: f64,
    #[serde(default)]
    cp: Option<f64>,
    #[serde(default, rename = "cp_segments")]
    cp_segments: Vec<CpSegment>,
}

impl<'de> Deserialize<'de> for Material {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let helper = MaterialHelper::deserialize(deserializer)?;
        let has_cp = helper.cp.is_some();
        let has_segments = !helper.cp_segments.is_empty();

        let cp = match (has_cp, has_segments) {
            (true, false) => Cp::Constant(helper.cp.unwrap()),
            (false, true) => Cp::from_piecewise_segments(helper.cp_segments),
            (true, true) => {
                return Err(serde::de::Error::custom(
                    "不能同时指定 'cp' 和 'cp_segments'，请只保留其中一个",
                ));
            }
            (false, false) => {
                return Err(serde::de::Error::custom(
                    "必须指定 'cp'（常数比热容，可设为 nan 使用 NASA 9 空气模型）或 'cp_segments'（分段多项式变比热容）",
                ));
            }
        };

        Ok(Material {
            molecular_weight: helper.molecular_weight,
            cp,
        })
    }
}

impl PartialEq for Material {
    fn eq(&self, other: &Self) -> bool {
        self.molecular_weight == other.molecular_weight && self.cp == other.cp
    }
}

impl Material {
    pub fn borrow_cp(&self) -> &Cp {
        &self.cp
    }

    pub fn borrow_mw(&self) -> &f64 {
        &self.molecular_weight
    }
}

impl Material {
    // 常量定义为 associated constants
    pub const UNIVERSAL_GAS_CONSTANT: f64 = 8.31446261815324; // J/(mol·K) R = N_A * k
    pub const AVOGADRO_CONSTANT: f64 = 6.02214076e23; // mol⁻¹
    pub const BOLTZMANN_CONSTANT: f64 = 1.380649e-23; // J/K
    pub const SPEED_OF_LIGHT: f64 = 299792458.0; // m/s

    /// 用摩尔质量kg/kmol和Cp(T)函数构造
    /// # 参数
    /// - `molecular_weight`: 摩尔质量(单位kg/kmol) 注意此处不是国际单位制
    /// - `cp`: 定压比热容(单位J/(kg·K))关于温度T(单位K)的函数
    pub fn new(molecular_weight: f64, cp: impl Fn(f64) -> f64 + Send + Sync + 'static) -> Self {
        Self {
            molecular_weight,
            cp: (Cp::new(cp)),
        }
    }

    /// 用摩尔质量kg/kmol和Cp(T)函数构造
    /// # 参数
    /// - `molecular_weight`: 摩尔质量(单位kg/kmol) 注意此处不是国际单位制
    /// - `cp`: 定压比热容(单位J/(kg·K)) 常数
    pub fn from_mw_cp(molecular_weight: f64, cp: f64) -> Self {
        Self {
            molecular_weight,
            cp: (Cp::Constant(cp)),
        }
    }

    /// 用气体常数和比热比构造
    /// # 参数
    /// - `r_gas`: 气体常数
    /// - `gamma`: 比热比
    pub fn from_rgas_gamma(r_gas: f64, gamma: f64) -> Self {
        Self {
            molecular_weight: Self::UNIVERSAL_GAS_CONSTANT / r_gas * 1e3, // kg/kmol
            cp: Cp::Constant(gamma / (gamma - 1.0) * r_gas),
        }
    }

    // 获取 Cp 值（在给定温度下）
    pub fn cp(&self, temperature: f64) -> f64 {
        self.cp.eval(temperature)
    }

    // 比热比 gamma，默认使用 Cp 计算 Cv = Cp - Rgas
    pub fn gamma(&self, temperature: f64) -> f64 {
        let cp = self.cp(temperature);
        let cv = cp - self.rgas();
        cp / cv
    }

    // 气体常数 Rgas = R / Mw （单位：J/(kg·K)）
    pub fn rgas(&self) -> f64 {
        Self::UNIVERSAL_GAS_CONSTANT / self.molecular_weight * 1e3
    }

    // 焓值 h(T)，参考点默认为 T_ref = 0.0
    pub fn enthalpy(&self, temperature: f64, t_ref: f64) -> f64 {
        if temperature == t_ref {
            return 0.0;
        }
        quadrature::gauss_legendre::integrate::<20, _>(|t| self.cp(t), t_ref, temperature)
    }
}

// 静态方法模拟：工厂函数
impl Material {
    // 空气（常数 Cp）
    pub fn air_constant() -> Self {
        Self::from_mw_cp(28.968, 1004.675)
    }

    // 空气（分段多项式 Cp）
    pub fn air_piecewise_polynomial() -> Self {
        Self::new(28.968, move |t| {
            let poly1 = |x: f64| {
                1161.482
                    + (-2.368819
                        + (0.01485511
                            + (-5.034909e-5
                                + (9.928570e-8
                                    + (-1.111097e-10 + (6.540196e-14 - 1.573588e-17 * x) * x)
                                        * x)
                                    * x)
                                * x)
                            * x)
                        * x
            };
            let poly2 = |x: f64| {
                -7069.814
                    + (33.70605
                        + (-0.05812760
                            + (5.421615e-5
                                + (-2.936679e-8
                                    + (9.237533e-12 + (-1.565553e-15 + 1.112335e-19 * x) * x)
                                        * x)
                                    * x)
                                * x)
                            * x)
                        * x
            };
            if (100. ..=1000.).contains(&t) {
                poly1(t)
            } else if (1000. ..=3000.).contains(&t) {
                poly2(t)
            } else if t < 100. {
                poly1(100.)
            } else if t > 3000. {
                poly2(3000.)
            } else {
                poly1(100.)
            }
        })
    }

    // 空气(NASA 9系数模型）
    pub fn air_nasa9piecewise_polynomial() -> Self {
        Self::new(28.968, move |t| {
            let poly1 = |x: f64| {
                (2898903. / x - 56496.26) / x
                    + 1437.799
                    + (-1.653609 + (0.003062254 + (-2.279138e-06 + 6.272365e-10 * x) * x) * x) * x
            };
            let poly2 = |x: f64| {
                (6.932494e+07 / x - 361053.2) / x
                    + 1476.665
                    + (-0.06138349 + (2.027963e-05 + (-3.075525e-09 + 1.888054e-13 * x) * x) * x)
                        * x
            };
            if (200. ..=1000.).contains(&t) {
                poly1(t)
            } else if (1000. ..=6000.).contains(&t) {
                poly2(t)
            } else if t < 200. {
                poly1(200.)
            } else if t > 6000. {
                poly2(6000.)
            } else {
                poly1(200.)
            }
        })
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_material_new() {
        let m = Material::new(25.0, |_| 1006.1);
        assert_eq!(m.cp(100.), 1006.1);
    }

    #[test]
    fn test_air_constant() {
        let m = Material::air_constant();
        assert!((m.cp(300.0) - 1004.675).abs() < 1e-10);
    }

    #[test]
    fn test_air_piecewise_polynomial() {
        let m = Material::air_piecewise_polynomial();
        // 验证在有效温度范围内返回合理值
        assert!(m.cp(300.0) > 0.0);
        assert!(m.cp(1500.0) > 0.0);
    }

    #[test]
    fn test_air_nasa9() {
        let m = Material::air_nasa9piecewise_polynomial();
        assert!(m.cp(300.0) > 0.0);
        assert!(m.cp(3000.0) > 0.0);
    }

    // ── TOML 反序列化测试 ─────────────────────────────────

    #[test]
    fn test_deserialize_constant_cp() {
        let toml_str = r#"
molecular_weight = 28.968
cp = 1004.675
"#;
        let m: Material = toml::from_str(toml_str).unwrap();
        assert!((m.molecular_weight - 28.968).abs() < 1e-10);
        assert!(
            (m.cp(300.0) - 1004.675).abs() < 1e-10,
            "got {}",
            m.cp(300.0)
        );
        assert!((m.cp(8000.0) - 1004.675).abs() < 1e-10);
    }

    #[test]
    fn test_deserialize_piecewise_cp() {
        let toml_str = r#"
molecular_weight = 28.968

[[cp_segments]]
t_min = 200.0
t_max = 1000.0
coefficients = [1000.0, 0.2]

[[cp_segments]]
t_min = 1000.0
t_max = 3000.0
coefficients = [800.0, 0.4]
"#;
        let m: Material = toml::from_str(toml_str).unwrap();
        // 第一段: 1000 + 0.2*500 = 1100
        assert!((m.cp(500.0) - 1100.0).abs() < 1e-10);
        // 第二段: 800 + 0.4*2000 = 1600
        assert!((m.cp(2000.0) - 1600.0).abs() < 1e-10);
    }

    #[test]
    fn test_deserialize_piecewise_constant_polynomial() {
        // 零次多项式（常数）→ 等效于 cp = 500
        let toml_str = r#"
molecular_weight = 44.01

[[cp_segments]]
t_min = 100.0
t_max = 5000.0
coefficients = [846.0]
"#;
        let m: Material = toml::from_str(toml_str).unwrap();
        assert!((m.cp(300.0) - 846.0).abs() < 1e-10);
        assert!((m.cp(2000.0) - 846.0).abs() < 1e-10);
    }

    #[test]
    fn test_deserialize_high_degree_polynomial() {
        // 5次多项式: 1 + 2T + 3T² + 4T³ + 5T⁴
        let toml_str = r#"
molecular_weight = 18.0

[[cp_segments]]
t_min = 0.0
t_max = 1000.0
coefficients = [1.0, 2.0, 3.0, 4.0, 5.0]
"#;
        let m: Material = toml::from_str(toml_str).unwrap();
        let t: f64 = 2.0;
        let expected = 1.0 + 2.0 * t + 3.0 * t.powi(2) + 4.0 * t.powi(3) + 5.0 * t.powi(4);
        assert!((m.cp(t) - expected).abs() < 1e-10);
    }

    #[test]
    fn test_deserialize_cp_nan() {
        // cp = nan → Constant(f64::NAN)，后续 config.rs 解析会替换为 NASA 9
        let toml_str = r#"
molecular_weight = 28.968
cp = nan
"#;
        let m: Material = toml::from_str(toml_str).unwrap();
        assert!(m.cp(300.0).is_nan());
    }

    #[test]
    fn test_deserialize_both_cp_and_segments_error() {
        let toml_str = r#"
molecular_weight = 28.968
cp = 1004.0

[[cp_segments]]
t_min = 200.0
t_max = 1000.0
coefficients = [1.0, 2.0]
"#;
        let result: Result<Material, _> = toml::from_str(toml_str);
        assert!(result.is_err());
        let err = result.unwrap_err().to_string();
        assert!(
            err.contains("cp"),
            "expected error about cp conflict, got: {err}"
        );
    }

    #[test]
    fn test_deserialize_neither_cp_nor_segments_error() {
        let toml_str = r#"
molecular_weight = 28.968
"#;
        let result: Result<Material, _> = toml::from_str(toml_str);
        assert!(result.is_err());
        let err = result.unwrap_err().to_string();
        assert!(
            err.contains("cp") || err.contains("cp_segments"),
            "expected error about missing cp, got: {err}"
        );
    }

    #[test]
    fn test_deserialize_unknown_field_error() {
        let toml_str = r#"
molecular_weight = 28.968
cp = 1004.0
unknown_field = 42
"#;
        // deny_unknown_fields 会拒绝未识别的字段
        let result: Result<Material, _> = toml::from_str(toml_str);
        assert!(result.is_err());
    }
}
