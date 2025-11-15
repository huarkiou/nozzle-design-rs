use serde::{Deserialize, Serialize};
use std::fs;
use std::path::Path;

use crate::Material;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Control {
    #[serde(rename = "irrotational")]
    pub irrotational: bool, // 无旋特征线法(true)或有旋特征线法(false)
    #[serde(rename = "axisymmetric")]
    pub axisymmetric: bool, // false代表二维平面问题，true代表二维轴对称问题
    #[serde(rename = "eps")]
    pub eps: f64, // 相对误差小于eps视为收敛
    #[serde(rename = "n_correction_max")]
    pub n_correction_max: i32, // 特征线法基本过程计算中欧拉预估校正迭代的最大校正次数
    #[serde(rename = "n_inlet")]
    pub n_inlet: i32, // 入口初值线上特征线网格点数
}

impl Default for Control {
    fn default() -> Self {
        Self {
            irrotational: false,
            axisymmetric: false,
            eps: 1e-7,
            n_correction_max: 40,
            n_inlet: 101,
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
    pub theta: f64, // 来流气流方向角(rad)
}

impl Default for Inlet {
    fn default() -> Self {
        Self {
            p_total: 0.0,
            temperature_total: 0.0,
            ma: 0.0,
            theta: 0.0,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Geometry {
    #[serde(rename = "height_i")]
    pub height_i: f64, // 进口高度(axisymmetric=false)或截面半径(axisymmetric=true)/m
    #[serde(rename = "height_e")]
    pub height_e: f64, // 目标出口高度(axisymmetric=false)或截面半径(axisymmetric=true)/m
    #[serde(rename = "length")]
    pub length: f64, // 目标总长度/m
    #[serde(rename = "width")]
    pub width: f64, // 喷管宽度(二维平面问题中axisymmetric=false)
}

impl Default for Geometry {
    fn default() -> Self {
        Self {
            height_i: 0.0,
            height_e: 0.0,
            length: 0.0,
            width: 0.0,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NozzleConfig {
    pub control: Control,
    pub material: Material,
    pub inlet: Inlet,
    pub geometry: Geometry,
    #[serde(rename = "radius_throat")]
    pub radius_throat: f64, // 喉部过渡圆弧半径 / m
    #[serde(rename = "theta_a")]
    pub theta_a: f64, // 初始膨胀角 / rad
    #[serde(rename = "p_ambient")]
    pub p_ambient: f64, // 背压 / Pa
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
        if self.radius_throat < 0.0 {
            return Err("喉部半径不能为负".to_string());
        }
        if self.p_ambient < 0.0 {
            return Err("背压不能为负".to_string());
        }
        if self.geometry.height_i <= 0.0 || self.geometry.height_e <= 0.0 {
            return Err("进出口高度必须大于0".to_string());
        }
        if self.geometry.length <= 0.0 {
            return Err("喷管长度必须大于0".to_string());
        }
        if self.inlet.ma < 0.0 {
            return Err("马赫数不能为负".to_string());
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
            radius_throat: self.radius_throat,
            theta_a: self.theta_a,
            p_ambient: self.p_ambient,
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
    pub radius_throat: f64,
    pub theta_a: f64,
    pub p_ambient: f64,
}

impl From<NozzleConfig> for NozzleDesign {
    fn from(config: NozzleConfig) -> Self {
        NozzleDesign {
            control: config.control,
            material: config.material,
            inlet: config.inlet,
            geometry: config.geometry,
            radius_throat: config.radius_throat,
            theta_a: config.theta_a,
            p_ambient: config.p_ambient,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new() {
        let config = NozzleConfig {
            control: Control {
                irrotational: true,
                axisymmetric: true,
                eps: 1e-6,
                n_correction_max: 50,
                n_inlet: 201,
            },
            material: Material::air_piecewise_polynomial(), // 使用理想空气
            inlet: Inlet {
                p_total: 200000.0,        // 200 kPa
                temperature_total: 800.0, // 800 K
                ma: 0.1,                  // 低马赫数入口
                theta: 0.0,
            },
            geometry: Geometry {
                height_i: 0.02, // 2cm 进口
                height_e: 0.08, // 8cm 出口
                length: 0.3,    // 30cm 长度
                width: 0.1,     // 10cm 宽度（平面问题）
            },
            radius_throat: 0.01, // 1cm 喉部半径
            theta_a: 0.2,        // 0.2rad 初始膨胀角
            p_ambient: 101325.0, // 标准大气压
        };

        // 验证并构建设计
        config.validate().expect("config invalid");
        let design = config.build_design().unwrap();

        // 保存为配置文件
        // config.to_toml_file("generated_config.toml");

        println!("\n程序化创建的喷管设计参数:");
        println!("{:#?}", design);
    }
}
