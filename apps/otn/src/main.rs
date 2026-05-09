use std::{path::PathBuf, process};

use clap::{Parser, Subcommand};

use aero::nozzle::{ConstraintNozzle, NozzleConfig};

mod export;

/// Optimal Thrust Nozzle generator using Method of Characteristics.
#[derive(Parser)]
#[command(
    name = "otn",
    version,
    about = "Generate optimal thrust nozzle (OTN) using Method of Characteristics (MOC)."
)]
struct Cli {
    #[command(subcommand)]
    command: Option<Command>,

    /// TOML configuration file
    #[arg(default_value = "otn.toml")]
    configfile: PathBuf,
}

#[derive(Subcommand)]
enum Command {
    /// Generate a default otn.toml and exit
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
    let config = NozzleConfig::from_toml_file(&cli.configfile)?;
    config
        .validate()
        .map_err(|e| format!("Config validation failed: {e}"))?;
    let config = config
        .build()
        .map_err(|e| format!("Config build failed: {e}"))?;

    let mut nozzle = ConstraintNozzle::new_otn(config.clone());
    nozzle.run();
    let charlines = nozzle.get_assembly_charlines();

    if charlines.is_empty() {
        return Err("Nozzle computation produced no charlines".into());
    }

    let prefix = &config.io.output_prefix;
    let field_filename = format!("{prefix}field_data.txt");
    charlines.write_to_file(&field_filename, false)?;
    eprintln!("  Field data written to {field_filename}");

    export::export_boundary_geometry(&charlines, &config)?;
    eprintln!("  Boundary geometry written to {}geo_all.dat", prefix);

    export::print_summary(&charlines, &config, nozzle.theta_a());

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
