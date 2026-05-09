use std::{path::PathBuf, process};

use aero::moc::{CharLines, read_charlines_from_file_checked};
use aero::streamline_trace::{StreamlineConfig, StreamlineTrace};
use clap::Parser;
use geometry::obj::ObjModel;
use geometry::wallpoints::write_boundaries;
use geometry::{Circle, ClosedCurve, Ellipse, Rectangular, SuperEllipse, UserDefined};

/// 3D streamline tracing nozzle design.
///
/// Reads a MOC base flow field (characteristic lines) and
/// inlet/outlet cross-section shapes, traces streamlines through the
/// field, and exports the resulting 3D nozzle geometry.
#[derive(Parser)]
#[command(
    name = "streamlinetracenozzle",
    version,
    about = "3D streamline tracing nozzle design using MOC base flow field"
)]
struct Cli {
    /// TOML configuration file
    #[arg(default_value = "StreamlineTraceNozzle.toml")]
    configfile: PathBuf,
}

fn main() {
    if let Err(e) = run() {
        eprintln!("Error: {e}");
        process::exit(1);
    }
}

// ─────────────────────────────────────────────────────────
// Configuration structs (mirror the TOML sections)
// ─────────────────────────────────────────────────────────

#[derive(serde::Deserialize)]
#[serde(rename_all = "PascalCase")]
struct Config {
    #[serde(default)]
    control: ControlConfig,
    base_fluid_field: BaseFluidFieldConfig,
    inlet: ShapeConfig,
    outlet: ShapeConfig,
}

#[derive(serde::Deserialize)]
#[serde(deny_unknown_fields)]
struct ControlConfig {
    #[serde(default = "default_n_theta")]
    n_theta: usize,
    #[serde(default = "default_n_axis")]
    n_axis: usize,
    #[serde(default)]
    monotonic: bool,
    #[serde(default)]
    weight_parameter_a: f64,
    #[serde(default)]
    export_obj: bool,
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

#[derive(serde::Deserialize)]
#[serde(deny_unknown_fields)]
struct BaseFluidFieldConfig {
    axisymmetric: bool,
    datasource_inlet: PathBuf,
    datasource_outlet: PathBuf,
}

#[derive(serde::Deserialize)]
#[serde(deny_unknown_fields)]
struct ShapeConfig {
    #[serde(default = "default_normalized")]
    normalized: bool,
    shape: String,
    #[serde(default)]
    center: Option<Vec<f64>>,
    #[serde(default)]
    radius: Option<f64>,
    #[serde(default)]
    a: Option<f64>,
    #[serde(default)]
    b: Option<f64>,
    #[serde(default)]
    alpha: Option<f64>,
    #[serde(default)]
    length: Option<f64>,
    #[serde(default)]
    width: Option<f64>,
    #[serde(default)]
    n: Option<f64>,
    #[serde(default)]
    datasource: Option<PathBuf>,
}

fn default_normalized() -> bool {
    true
}

// ─────────────────────────────────────────────────────────
// Shape factory
// ─────────────────────────────────────────────────────────

fn build_shape(
    shape_config: &ShapeConfig,
    normalize_factor: f64,
) -> Result<Box<dyn ClosedCurve>, String> {
    let factor = if shape_config.normalized {
        1.0
    } else {
        normalize_factor
    };
    let center: (f64, f64) = shape_config
        .center
        .as_ref()
        .map(|c| {
            if c.len() == 2 {
                Ok((c[0] / factor, c[1] / factor))
            } else {
                Err(format!(
                    "center must have exactly 2 elements, got {}",
                    c.len()
                ))
            }
        })
        .transpose()
        .map(|r| r.unwrap_or((0.0, 0.0)))
        .unwrap_or((0.0, 0.0));

    match shape_config.shape.as_str() {
        "circle" => {
            let r = shape_config
                .radius
                .ok_or("circle shape requires 'radius'")?
                / factor;
            Ok(Box::new(Circle::new(r, center)))
        }
        "ellipse" => {
            let a = shape_config.a.ok_or("ellipse shape requires 'a'")? / factor;
            let b = shape_config.b.ok_or("ellipse shape requires 'b'")? / factor;
            let alpha = shape_config
                .alpha
                .map(|deg| deg.to_radians())
                .unwrap_or(0.0);
            Ok(Box::new(Ellipse::new(a, b, alpha, center)))
        }
        "rectangular" => {
            let length = shape_config
                .length
                .ok_or("rectangular shape requires 'length'")?
                / factor;
            let width = shape_config
                .width
                .ok_or("rectangular shape requires 'width'")?
                / factor;
            let alpha = shape_config
                .alpha
                .map(|deg| deg.to_radians())
                .unwrap_or(0.0);
            Ok(Box::new(Rectangular::new(length, width, alpha, center)))
        }
        "superellipse" => {
            let a = shape_config.a.ok_or("superellipse shape requires 'a'")? / factor;
            let b = shape_config.b.ok_or("superellipse shape requires 'b'")? / factor;
            let power = shape_config.n.unwrap_or(2.0);
            let alpha = shape_config
                .alpha
                .map(|deg| deg.to_radians())
                .unwrap_or(0.0);
            Ok(Box::new(SuperEllipse::new(a, b, power, alpha, center)))
        }
        "userdefined" => {
            let ds = shape_config
                .datasource
                .as_ref()
                .ok_or("userdefined shape requires 'datasource'")?;
            let points =
                geometry::wallpoints::read_points_2d(ds.to_str().unwrap_or("")).map_err(|e| {
                    format!(
                        "failed to read userdefined datasource '{}': {e}",
                        ds.display()
                    )
                })?;
            let pts: Vec<geometry::Point3d> = points
                .into_iter()
                .map(|(z, y)| geometry::Point3d::new(z / factor, y / factor, 0.0))
                .collect();
            let alpha = shape_config
                .alpha
                .map(|deg| deg.to_radians())
                .unwrap_or(0.0);
            Ok(Box::new(UserDefined::new(pts, alpha, Some(center))))
        }
        other => Err(format!(
            "unknown shape '{}'. Supported: circle, ellipse, rectangular, superellipse, userdefined",
            other
        )),
    }
}

// ─────────────────────────────────────────────────────────
// Main pipeline
// ─────────────────────────────────────────────────────────

fn run() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();

    // ── 1. Load TOML config ──
    let config_path = &cli.configfile;
    eprintln!("Reading configuration from '{}'...", config_path.display());
    let config_toml = std::fs::read_to_string(config_path).map_err(|e| {
        format!(
            "failed to read config file '{}': {e}",
            config_path.display()
        )
    })?;
    let config: Config = toml::from_str(&config_toml)?;

    // ── 2. Read base flow field ──
    eprintln!("Reading base fluid field data...");
    let lines_inlet: CharLines = read_charlines_from_file_checked(
        config.base_fluid_field.datasource_inlet.to_str().unwrap(),
    )
    .map_err(|e| format!("Failed to load inlet charlines: {e}"))?;

    let lines_outlet: CharLines = read_charlines_from_file_checked(
        config.base_fluid_field.datasource_outlet.to_str().unwrap(),
    )
    .map_err(|e| format!("Failed to load outlet charlines: {e}"))?;

    if lines_inlet.is_empty() || lines_outlet.is_empty() {
        return Err("Inlet or outlet charlines are empty".into());
    }

    // ── 3. Normalize factor = inlet height (y of first point in first line) ──
    let normalize_factor = lines_inlet
        .first()
        .and_then(|l| l.first())
        .map(|p| p.y)
        .unwrap_or(1.0);

    // ── 4. Build cross‑section shapes ──
    eprintln!("Reading inlet/outlet shape...");
    let inlet_shape =
        build_shape(&config.inlet, normalize_factor).map_err(|e| format!("Inlet shape: {e}"))?;
    let outlet_shape =
        build_shape(&config.outlet, normalize_factor).map_err(|e| format!("Outlet shape: {e}"))?;

    // ── 5. Build StreamlineTrace and run ──
    let trace_config = StreamlineConfig {
        axisymmetric: config.base_fluid_field.axisymmetric,
        datasource_inlet: lines_inlet,
        datasource_outlet: lines_outlet,
        inlet_shape,
        outlet_shape,
        monotonic: config.control.monotonic,
        weight_parameter_a: config.control.weight_parameter_a,
    };

    eprintln!(
        "Tracing streamlines (n_theta = {}, n_axis = {})...",
        config.control.n_theta, config.control.n_axis
    );
    let mut tracer = StreamlineTrace::new(trace_config);
    tracer.run(config.control.n_theta, config.control.n_axis)?;

    // ── 6. Write geometry files ──
    eprintln!("Writing geometry files...");

    let write_model =
        |model: &[Vec<geometry::Point3d>], name: &str| -> Result<(), Box<dyn std::error::Error>> {
            if !model.is_empty() {
                // Repeat the first boundary at the end to close the loop.
                let mut closed = model.to_vec();
                closed.push(model[0].clone());

                let dat_path = format!("{name}.dat");
                write_boundaries(&dat_path, &closed)?;
                eprintln!("  Written: {dat_path}");

                if config.control.export_obj {
                    let mut obj = ObjModel::from_wallpoints(&closed);
                    obj.auto_normal();
                    let obj_path = format!("{name}.obj");
                    obj.write_obj(&obj_path)?;
                    eprintln!("  Written: {obj_path}");
                }
            }
            Ok(())
        };

    write_model(&tracer.model, "model")?;
    write_model(&tracer.downstream, "downstream")?;
    write_model(&tracer.upstream, "upstream")?;

    eprintln!("Done.");
    Ok(())
}
