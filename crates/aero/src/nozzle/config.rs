use serde::{Deserialize, Serialize};
use std::{
    f64::{self, consts::PI},
    fs,
    path::Path,
};

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
    pub n_correction_max: u16, // 特征线法基本过程计算中欧拉预估校正迭代的最大校正次数
    #[serde(rename = "n_inlet")]
    pub n_inlet: u16, // 入口初值线上特征线网格点数 《===若计算发散可增大此参数重新尝试
}

impl Control {
    pub fn validate(&self) -> Result<(), String> {
        if self.eps > 1e-3 || self.eps < 1e-14 {
            Err("容差(eps)大小不合理".to_string())
        } else if self.n_correction_max < 1 {
            Err("预估校正迭代次数不合理".to_string())
        } else if self.n_inlet < 2 {
            Err("入口网格密度不合理".to_string())
        } else {
            Ok(())
        }
    }
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

impl Geometry {
    pub fn validate(&self) -> Result<(), String> {
        if self.height_i <= 0. {
            Err("进口高度不合理".to_string())
        } else if self.height_e <= self.height_i {
            Err("出口高度不合理".to_string())
        } else if self.length <= 0. {
            Err("目标长度不合理".to_string())
        } else {
            Ok(())
        }
    }
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
    pub theta: f64, // 来流气流方向角(rad) 《===不建议使用此参数
}

impl Inlet {
    pub fn validate(&self) -> Result<(), String> {
        if self.p_total <= 0. {
            Err("总压应为正数".to_string())
        } else if self.temperature_total <= 0. {
            Err("总温应为正数".to_string())
        } else if self.ma < 1. {
            Err("入口气流马赫数应不小于1".to_string())
        } else if self.theta != 0. {
            Err("Inlet.theta功能还未实现，应当将其设置为0".to_string())
        } else {
            Ok(())
        }
    }
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
    pub theta_a: f64, // 初始膨胀角(rad) *若为负数或nan则由程序自动迭代计算选取，这也会导致[Geometry].height_e失效*
}

impl Throat {
    pub fn validate(&self) -> Result<(), String> {
        if self.radius_throat < 0. {
            Err("喉部过渡圆弧半径应为正数".to_string())
        } else if self.theta_a < 0. || self.theta_a >= PI {
            Err("初始膨胀角只能为0~89.99°或NaN".to_string())
        } else {
            Ok(())
        }
    }
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

impl Outlet {
    pub fn validate(&self) -> Result<(), String> {
        if self.p_ambient <= 0. {
            Err("出口背压应为正数".to_string())
        } else {
            Ok(())
        }
    }
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

impl IO {
    pub fn validate(&self) -> Result<(), String> {
        Ok(())
    }
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
        let mut config: NozzleConfig = toml::from_str(&contents)?;
        config.inlet.theta = config.inlet.theta.to_radians();
        config.throat.theta_a = config.throat.theta_a.to_radians();
        Ok(config)
    }

    /// 保存为 TOML 文件
    pub fn to_toml_file<P: AsRef<Path>>(&self, path: P) -> Result<(), Box<dyn std::error::Error>> {
        let toml_str = toml::to_string_pretty(&NozzleConfig {
            control: self.control.clone(),
            material: self.material.clone(),
            geometry: self.geometry.clone(),
            inlet: Inlet {
                theta: self.inlet.theta.to_degrees(),
                ..self.inlet
            },
            throat: Throat {
                theta_a: self.throat.theta_a.to_degrees(),
                ..self.throat
            },
            outlet: self.outlet.clone(),
            io: self.io.clone(),
        })?;
        if let Some(parent) = Path::new(path.as_ref()).parent() {
            fs::create_dir_all(parent)?;
        }
        fs::write(path, toml_str)?;
        Ok(())
    }

    /// 验证配置参数的有效性
    pub fn validate(&self) -> Result<(), String> {
        self.control.validate()?;
        self.geometry.validate()?;
        self.inlet.validate()?;
        self.throat.validate()?;
        self.outlet.validate()?;
        self.io.validate()?;
        if self.control.axisymmetric == false && self.geometry.width <= 0. {
            Err("在二维平面问题中第三维深度应为正数".to_string())
        } else {
            Ok(())
        }
    }

    /// 构建最终的喷管设计对象
    pub fn build_design(self) -> Result<NozzleConfig, String> {
        self.validate()?;
        Ok(NozzleConfig {
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

#[cfg(test)]
mod tests {
    use std::env;

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
        let _ = config.clone().build_design().unwrap();

        // 保存为配置文件
        let output_dir = env::temp_dir().join("generated_config.toml");

        println!("{:}", output_dir.to_str().unwrap_or_default());
        let _ = config.to_toml_file(output_dir).unwrap();
    }
}
