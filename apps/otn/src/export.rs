use std::{
    fmt::Write as FmtWrite,
    fs::File,
    io::{self, BufWriter, Write as IoWrite},
};

use math::geometry::Coord2d;

use aero::moc::{CharLine, CharLines, MocPoint};
use aero::nozzle::NozzleConfig;

// ─────────────────────────────────────────────────────────
// Boundary geometry export
// ─────────────────────────────────────────────────────────

pub fn export_boundary_geometry(
    charlines: &CharLines,
    config: &NozzleConfig,
) -> Result<(), Box<dyn std::error::Error>> {
    let prefix = &config.io.output_prefix;
    let filename = format!("{prefix}geo_all.dat");
    let file = File::create(&filename)?;
    let mut w = BufWriter::new(file);

    let wall_points: Vec<&MocPoint> = charlines.iter().filter_map(|line| line.first()).collect();

    if wall_points.len() < 3 {
        return Err("Need at least 3 charlines for boundary geometry".into());
    }

    let n = wall_points.len();
    let inlet_wall = wall_points[0];
    let y_t = inlet_wall.y;
    let exit_x = wall_points[n - 1].x;
    let inlet_x = inlet_wall.x;
    let nozzle_len = exit_x - inlet_x;
    let y_max = wall_points.iter().map(|p| p.y).fold(0.0_f64, f64::max);

    let downward_ratio = 1.0 / 5.0;
    let inlet_indent = 1.5;
    let bottom_extend = 5.0;
    let right_extra = 3.0;

    let y_bottom = -(y_max * downward_ratio);
    let x_left = -y_t - inlet_indent * y_t;
    let x_right = exit_x + right_extra * nozzle_len;
    let x_bottom_right = exit_x + bottom_extend * nozzle_len;

    // 1. sj_wall_u: upper wall
    write_segment(
        &mut w,
        wall_points[1..n - 1].iter().map(|p| (p.x, p.y, 0.0)),
    )?;
    w.write_all(b"\n")?;
    // 2. sj_wall_d: lower wall (axis)
    write_segment(&mut w, [(inlet_x, 0.0, 0.0), (exit_x, 0.0, 0.0)])?;
    w.write_all(b"\n")?;
    // 3. sj_inlet
    write_segment(
        &mut w,
        [
            (-y_t, 0.0, 0.0),
            (-y_t, y_t, 0.0),
            (inlet_wall.x, inlet_wall.y, 0.0),
        ],
    )?;
    w.write_all(b"\n")?;
    // 4. sj_throatu
    write_segment(
        &mut w,
        [(-y_t, y_t, 0.0), (inlet_wall.x, inlet_wall.y, 0.0)],
    )?;
    w.write_all(b"\n")?;
    // 5. sj_throatd
    write_segment(&mut w, [(-y_t, 0.0, 0.0), (inlet_wall.x, 0.0, 0.0)])?;
    w.write_all(b"\n")?;
    // 6. sj_walld_d
    write_segment(&mut w, [(inlet_x, 0.0, 0.0), (x_bottom_right, 0.0, 0.0)])?;
    w.write_all(b"\n")?;
    // 7. bottomleft_u
    write_segment(&mut w, [(x_left, y_t, 0.0), (x_left, y_max, 0.0)])?;
    w.write_all(b"\n")?;
    // 8. bottomleft_in
    write_segment(&mut w, [(x_left, y_t, 0.0), (-y_t, y_t, 0.0)])?;
    w.write_all(b"\n")?;
    // 9. bottom
    write_segment(
        &mut w,
        [
            (x_left, y_bottom, 0.0),
            (x_right, y_bottom, 0.0),
            (x_bottom_right, y_bottom, 0.0),
        ],
    )?;
    w.write_all(b"\n")?;
    // 10. right
    write_segment(&mut w, [(x_right, y_max, 0.0), (x_right, y_bottom, 0.0)])?;
    w.write_all(b"\n")?;
    // 11. top
    write_segment(&mut w, [(x_left, y_max, 0.0), (x_right, y_max, 0.0)])?;
    w.write_all(b"\n")?;
    // 12. topleft_in
    write_segment(&mut w, [(x_left, y_max, 0.0), (x_left, y_t, 0.0)])?;
    w.write_all(b"\n")?;
    // 13. topleft_d
    write_segment(
        &mut w,
        [
            (x_left, y_t, 0.0),
            (x_left, 0.0, 0.0),
            (x_left, y_bottom, 0.0),
        ],
    )?;
    w.write_all(b"\n")?;
    // 14. sj_walld_plus
    write_segment(&mut w, [(exit_x, 0.0, 0.0), (x_bottom_right, 0.0, 0.0)])?;
    w.write_all(b"\n")?;
    // 15. sj_outlet
    let exit_wall = wall_points[n - 1];
    write_segment(
        &mut w,
        [(exit_wall.x, exit_wall.y, 0.0), (exit_x, 0.0, 0.0)],
    )?;

    w.flush()?;
    Ok(())
}

fn write_segment(
    w: &mut impl IoWrite,
    pts: impl IntoIterator<Item = (f64, f64, f64)>,
) -> io::Result<()> {
    for (x, y, z) in pts {
        writeln!(w, "{:.9} {:.9} {:.9}", x, y, z)?;
    }
    Ok(())
}

// ─────────────────────────────────────────────────────────
// Summary
// ─────────────────────────────────────────────────────────

pub fn print_summary(charlines: &CharLines, config: &NozzleConfig, theta_a: f64) {
    let wall_points: Vec<&MocPoint> = charlines.iter().filter_map(|line| line.first()).collect();

    if wall_points.len() < 2 {
        eprintln!("  WARNING: Not enough wall points for summary");
        return;
    }

    let inlet_wall = wall_points[0];
    let exit_wall = wall_points[wall_points.len() - 1];

    let wall_coords: Vec<Coord2d> = wall_points.iter().map(|p| Coord2d::new(p.x, p.y)).collect();
    let wall_len: f64 = wall_coords
        .windows(2)
        .map(|w| w[0].distance_to(&w[1]))
        .sum();

    let exit_ma = exit_wall.mach_number();
    let inlet_ma = inlet_wall.mach_number();

    let area_type = if config.control.axisymmetric {
        aero::moc::AreaType::Axisymmetric
    } else {
        aero::moc::AreaType::Planar(config.geometry.width)
    };

    let exit_line = charlines.last().cloned().unwrap_or_default();
    let thrust = CharLine::thrust(&exit_line, config.outlet.p_ambient, area_type);
    let mass_flow = CharLine::mass_flow_rate(&exit_line, area_type);

    let exit_height = exit_wall.y;
    let inlet_height = inlet_wall.y;

    let mut summary = String::new();
    let _ = writeln!(summary);
    let _ = writeln!(summary, "====== OptimumNozzle Summary ======");
    let _ = writeln!(summary, "  theta_a      = {:.6} deg", theta_a.to_degrees());
    let _ = writeln!(summary, "  wall length  = {:.6} m", wall_len);
    let _ = writeln!(summary, "  inlet height = {:.6} m", inlet_height);
    let _ = writeln!(summary, "  exit height  = {:.6} m", exit_height);
    let _ = writeln!(summary, "  inlet Ma     = {:.6}", inlet_ma);
    let _ = writeln!(summary, "  exit Ma      = {:.6}", exit_ma);
    let _ = writeln!(summary, "  thrust       = {:.3} N", thrust);
    let _ = writeln!(summary, "  mass flow    = {:.6} kg/s", mass_flow);
    let _ = writeln!(summary, "====================================");

    eprintln!("{summary}");
}
