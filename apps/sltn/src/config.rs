use std::path::PathBuf;

use serde::Deserialize;

#[derive(Deserialize)]
#[serde(rename_all = "PascalCase")]
pub struct Config {
    #[serde(default)]
    pub control: ControlConfig,
    pub base_fluid_field: BaseFluidFieldConfig,
    pub inlet: ShapeConfig,
    pub outlet: ShapeConfig,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ControlConfig {
    #[serde(default = "default_n_theta")]
    pub n_theta: usize,
    #[serde(default = "default_n_axis")]
    pub n_axis: usize,
    #[serde(default)]
    pub monotonic: bool,
    #[serde(default)]
    pub weight_parameter_a: f64,
    #[serde(default)]
    pub export_obj: bool,
}

fn default_n_theta() -> usize {
    66
}
fn default_n_axis() -> usize {
    111
}

impl Default for ControlConfig {
    fn default() -> Self {
        Self {
            n_theta: 66,
            n_axis: 111,
            monotonic: false,
            weight_parameter_a: 0.0,
            export_obj: false,
        }
    }
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BaseFluidFieldConfig {
    pub axisymmetric: bool,
    pub datasource_inlet: PathBuf,
    pub datasource_outlet: PathBuf,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ShapeConfig {
    #[serde(default = "default_normalized")]
    pub normalized: bool,
    pub shape: String,
    #[serde(default)]
    pub center: Option<Vec<f64>>,
    #[serde(default)]
    pub radius: Option<f64>,
    #[serde(default)]
    pub a: Option<f64>,
    #[serde(default)]
    pub b: Option<f64>,
    #[serde(default)]
    pub alpha: Option<f64>,
    #[serde(default)]
    pub length: Option<f64>,
    #[serde(default)]
    pub width: Option<f64>,
    #[serde(default)]
    pub n: Option<f64>,
    #[serde(default)]
    pub datasource: Option<PathBuf>,
}

fn default_normalized() -> bool {
    true
}
