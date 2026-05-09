use std::{path::PathBuf, process};

use aero::moc::read_charlines_from_file_checked;
use aero::streamline_trace::{StreamlineConfig, StreamlineTrace};
use clap::{Parser, Subcommand};
use geometry::obj::ObjModel;
use geometry::wallpoints::write_boundaries;

mod config;
mod shape;

use config::Config;

/// 3D streamline tracing nozzle design.
#[derive(Parser)]
#[command(
    name = "sltn",
    version,
    about = "3D streamline tracing nozzle design using MOC base flow field"
)]
struct Cli {
    #[command(subcommand)]
    command: Option<Command>,

    /// TOML configuration file
    #[arg(default_value = "sltn.toml")]
    configfile: PathBuf,
}

#[derive(Subcommand)]
enum Command {
    /// Generate a default sltn.toml and exit
    Init,
}

fn main() {
    let cli = Cli::parse();
    if matches!(cli.command, Some(Command::Init)) {
        if let Err(e) = generate_default_config(&cli.configfile) {
            eprintln!("Error: {e}");
            process::exit(1);
        }
        return;
    }
    if let Err(e) = run(&cli) {
        eprintln!("Error: {e}");
        process::exit(1);
    }
}

fn run(cli: &Cli) -> Result<(), Box<dyn std::error::Error>> {
    let config_path = &cli.configfile;
    eprintln!("Reading configuration from '{}'...", config_path.display());
    let config_toml = std::fs::read_to_string(config_path).map_err(|e| {
        format!(
            "failed to read config file '{}': {e}",
            config_path.display()
        )
    })?;
    let config: Config = toml::from_str(&config_toml)?;

    eprintln!("Reading base fluid field data...");
    let lines_inlet = read_charlines_from_file_checked(
        config.base_fluid_field.datasource_inlet.to_str().unwrap(),
    )
    .map_err(|e| format!("Failed to load inlet charlines: {e}"))?;
    let lines_outlet = read_charlines_from_file_checked(
        config.base_fluid_field.datasource_outlet.to_str().unwrap(),
    )
    .map_err(|e| format!("Failed to load outlet charlines: {e}"))?;

    if lines_inlet.is_empty() || lines_outlet.is_empty() {
        return Err("Inlet or outlet charlines are empty".into());
    }

    let normalize_factor = lines_inlet
        .first()
        .and_then(|l| l.first())
        .map(|p| p.y)
        .unwrap_or(1.0);

    eprintln!("Reading inlet/outlet shape...");
    let inlet_shape = shape::build_shape(&config.inlet, normalize_factor)
        .map_err(|e| format!("Inlet shape: {e}"))?;
    let outlet_shape = shape::build_shape(&config.outlet, normalize_factor)
        .map_err(|e| format!("Outlet shape: {e}"))?;

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

    eprintln!("Writing geometry files...");
    let write_model =
        |model: &[Vec<geometry::Point3d>], name: &str| -> Result<(), Box<dyn std::error::Error>> {
            if !model.is_empty() {
                let mut closed = model.to_vec();
                closed.push(model[0].clone());
                write_boundaries(&format!("{name}.dat"), &closed)?;
                eprintln!("  Written: {name}.dat");
                if config.control.export_obj {
                    let mut obj = ObjModel::from_wallpoints(&closed);
                    obj.auto_normal();
                    obj.write_obj(&format!("{name}.obj"))?;
                    eprintln!("  Written: {name}.obj");
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

fn generate_default_config(path: &PathBuf) -> Result<(), Box<dyn std::error::Error>> {
    if path.exists() {
        return Err(format!("{} already exists", path.display()).into());
    }
    std::fs::write(path, DEFAULT_CONFIG)?;
    eprintln!("Default configuration written to {}", path.display());
    Ok(())
}

const DEFAULT_CONFIG: &str = include_str!("../configs/default.toml");
