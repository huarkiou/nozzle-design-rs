use serde::{Deserialize, Serialize};
use std::path::Path;
use std::{f64, fs};

use crate::Material;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Control {
    #[serde(rename = "irrotational")]
    pub irrotational: bool, // 无旋特征线法(true)或有旋特征线法(false)
    #[serde(rename = "axisymmetric")]
    pub axisymmetric: bool, // false代表二维平面问题，true代表二维轴对称问题
    #[serde(rename = "eps")]
    pub eps: f64, // 残差小于eps视为相等/收敛
    #[serde(rename = "n_correction_max")]
    pub n_correction_max: i32, // 特征线法基本过程计算中欧拉预估校正迭代的最大校正次数
    #[serde(rename = "n_inlet")]
    pub n_inlet: i32, // 入口初值线上特征线网格点数 《===若计算发散可增大此参数重新尝试
}

impl Default for Control {
    fn default() -> Self {
        Self {
            irrotational: true,
            axisymmetric: true,
            eps: 1e-7,
            n_correction_max: 40,
            n_inlet: 101,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Geometry {
    #[serde(rename = "height")]
    pub height_i: f64, // 进口高度(axisymmetric=false)或截面半径(axisymmetric=true)(m)
    #[serde(rename = "height_e")]
    pub height_e: f64, // 目标出口高度(axisymmetric=false)或截面半径(axisymmetric=true)(m) *若以最大推力为目标，对其约束则设置为nan*
    #[serde(rename = "length")]
    pub length: f64, // 喷管目标长度(m)
    #[serde(rename = "width")]
    pub width: f64, // 喷管的横向宽度(m)(仅当axisymmetric=false时有效)
}

impl Default for Geometry {
    fn default() -> Self {
        Self {
            height_i: 1.0,
            height_e: f64::NAN,
            length: 6.0,
            width: 1.0,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Inlet {
    #[serde(rename = "p_total")]
    pub p_total: f64, // 来流总压(Pa)
    #[serde(rename = "T_total")]
    pub temperature_total: f64, // 来流总温(K)
    #[serde(rename = "Ma")]
    pub ma: f64, // 来流马赫数
    #[serde(rename = "theta")]
    pub theta: f64, // 来流气流方向角(°) 《===不建议使用此参数
}

impl Default for Inlet {
    fn default() -> Self {
        Self {
            p_total: 800000.0,
            temperature_total: 2000.0,
            ma: 1.2,
            theta: 0.0,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Throat {
    #[serde(rename = "R_t")]
    pub radius_throat: f64, // 过渡圆弧半径(m)
    #[serde(rename = "theta_a")]
    pub theta_a: f64, // 初始膨胀角(°) *若为负数或nan则由程序自动迭代计算选取，这也会导致[Geometry].height_e失效*
}

impl Default for Throat {
    fn default() -> Self {
        Self {
            radius_throat: 0.0,
            theta_a: f64::NAN,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Outlet {
    #[serde(rename = "p_ambient")]
    pub p_ambient: f64, // 设计出口背压(Pa)
}

impl Default for Outlet {
    fn default() -> Self {
        Self { p_ambient: 7000. }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IO {
    output_prefix: String, // 输出文件名的前缀
}

impl Default for IO {
    fn default() -> Self {
        Self {
            output_prefix: Default::default(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NozzleConfig {
    #[serde(rename = "MOCControl")]
    pub control: Control,
    #[serde(rename = "Material")]
    pub material: Material,
    #[serde(rename = "Geometry")]
    pub geometry: Geometry,
    #[serde(rename = "Inlet")]
    pub inlet: Inlet,
    #[serde(rename = "Throat")]
    pub throat: Throat,
    #[serde(rename = "Outlet")]
    pub outlet: Outlet,
    #[serde(rename = "IO")]
    pub io: IO,
}

impl NozzleConfig {
    /// 从 TOML 文件加载配置
    pub fn from_toml_file<P: AsRef<Path>>(path: P) -> Result<Self, Box<dyn std::error::Error>> {
        let contents = fs::read_to_string(path)?;
        let config: NozzleConfig = toml::from_str(&contents)?;

        Ok(config)
    }

    /// 保存为 TOML 文件
    pub fn to_toml_file<P: AsRef<Path>>(&self, path: P) -> Result<(), Box<dyn std::error::Error>> {
        let toml_str = toml::to_string_pretty(self)?;
        fs::write(path, toml_str)?;
        Ok(())
    }

    /// 验证配置参数的有效性
    pub fn validate(&self) -> Result<(), String> {
        if self.throat.radius_throat < 0.0 {
            return Err("喉部半径不能为负".to_string());
        }
        if self.outlet.p_ambient <= 0.0 {
            return Err("背压不能为负".to_string());
        }
        if self.geometry.height_i <= 0.0 {
            return Err("进出口高度必须大于0".to_string());
        }
        if self.geometry.length <= 0.0 {
            return Err("喷管长度必须大于0".to_string());
        }
        if self.inlet.ma < 1.0 {
            return Err("马赫数不能小于1".to_string());
        }
        Ok(())
    }

    /// 构建最终的喷管设计对象
    pub fn build_design(self) -> Result<NozzleDesign, String> {
        self.validate()?;
        Ok(NozzleDesign {
            control: self.control,
            material: self.material,
            inlet: self.inlet,
            geometry: self.geometry,
            throat: self.throat,
            outlet: self.outlet,
            io: self.io,
        })
    }
}

// 内部计算用的设计结构体
#[derive(Debug, Clone)]
pub struct NozzleDesign {
    pub control: Control,
    pub material: Material,
    pub inlet: Inlet,
    pub geometry: Geometry,
    pub throat: Throat,
    pub outlet: Outlet,
    pub io: IO,
}

impl From<NozzleConfig> for NozzleDesign {
    fn from(config: NozzleConfig) -> Self {
        NozzleDesign {
            control: config.control,
            material: config.material,
            inlet: config.inlet,
            geometry: config.geometry,
            throat: config.throat,
            outlet: config.outlet,
            io: config.io,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new() {
        let config = NozzleConfig {
            control: Control::default(),
            material: Material::air_piecewise_polynomial(), // 使用理想空气
            inlet: Inlet::default(),
            geometry: Geometry::default(),
            throat: Throat::default(),
            outlet: Outlet::default(),
            io: IO::default(),
        };

        // 验证并构建设计
        config.validate().expect("config invalid");
        let design = config.clone().build_design().unwrap();

        // 保存为配置文件
        let _ = config.to_toml_file("generated_config.toml");

        println!("\n程序化创建的喷管设计参数:");
        println!("{:#?}", design);
    }
}
