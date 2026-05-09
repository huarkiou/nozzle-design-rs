use crate::moc::{CharLine, CharLines, MocPoint};
use geometry::{Point3d, WallPoints};

use super::intersect::calculate_streamline_intersection;

const DBG_STR_LOG: bool = false;

/// Trace a single streamline from `start_point` through the given
/// characteristic lines.
///
/// # Arguments
///
/// * `start_point` – The starting location on the first characteristic
///   line (or the inlet/outlet shape). Only `(x, y)` are used.
/// * `datasource` – Ordered characteristic lines. The first entry must
///   contain the start point's x‑coordinate.
/// * `x_positive` – Tracing direction: `true` = increasing x
///   (downstream), `false` = decreasing x (upstream).
///
/// # Returns
///
/// `Some(WallPoints)` containing the streamline points in order (start →
/// end), or `None` if tracing fails at the very first step.
pub fn trace_streamline(
    start_point: &Point3d,
    datasource: &CharLines,
    x_positive: bool,
) -> Option<WallPoints> {
    // ── basic validity ────────────────────────────────────────────
    if datasource.is_empty() {
        return None;
    }

    // ── step 1: bootstrap from the first characteristic line ──────
    let first_line = &datasource[0];
    if first_line.len() < 2 {
        return None;
    }

    let mut p_prev = initial_moc_point(start_point, first_line)?;
    if DBG_STR_LOG {
        eprintln!(
            "  trace start: ({:.6}, {:.6}) -> MOC ({:.6}, {:.6}, θ={:.4})",
            start_point.x,
            start_point.y,
            p_prev.x,
            p_prev.y,
            p_prev.flow_direction()
        );
    }

    let mut points: WallPoints = Vec::with_capacity(datasource.len().max(2));
    points.push(Point3d::new(p_prev.x, p_prev.y, 0.0));

    // ── step 2: march through the remaining characteristic lines ──
    for (line_idx, line) in datasource.iter().enumerate().skip(1) {
        if line.len() < 2 {
            continue;
        }

        // Try each segment of the line; take the first valid intersection
        let mut found = false;
        for seg in 0..(line.len() - 1) {
            if let Some(next_pt) = calculate_streamline_intersection(&p_prev, line, seg, x_positive)
            {
                // Skip duplicate positions (can happen when consecutive
                // charlines overlap at the same axial station).
                if (next_pt.x - p_prev.x).abs() < 1e-10 && (next_pt.y - p_prev.y).abs() < 1e-10 {
                    continue;
                }
                if DBG_STR_LOG {
                    eprintln!(
                        "    line {} seg {} -> ({:.6}, {:.6})",
                        line_idx, seg, next_pt.x, next_pt.y
                    );
                }
                points.push(Point3d::new(next_pt.x, next_pt.y, 0.0));
                p_prev = next_pt;
                found = true;
                break;
            }
        }

        if !found {
            if DBG_STR_LOG {
                eprintln!(
                    "  WARNING: no intersection on line {} for p_prev=({:.4},{:.4})",
                    line_idx, p_prev.x, p_prev.y
                );
            }
            // Skip this characteristic line — it may be a partial line
            // that does not span the radial position of this streamline.
            // Continue to the next line.
            continue;
        }
    }

    if points.is_empty() {
        None
    } else {
        Some(points)
    }
}

/// Project `start_point` onto the first characteristic line and return
/// the associated `MocPoint` (with interpolated flow parameters).
///
/// The start point comes from a cross‑section shape (`ClosedCurve`).
/// `start_point.x` is the lateral/spatial coordinate (z in shape space),
/// `start_point.y` is the vertical coordinate (y in shape space).
/// The distance from the symmetry axis is `hypot(x, y)`.
fn initial_moc_point(start_point: &Point3d, line: &CharLine) -> Option<MocPoint> {
    if line.is_empty() {
        return None;
    }

    // The characteristic line x is the axial station across all points.
    let line_x = line[0].x;

    // True radial distance from the symmetry axis (not just y).
    let r_target = start_point.x.hypot(start_point.y);

    // Clamp to the flow field radial bounds.
    let r_min = line.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);
    let r_max = line.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max);
    let r_target = r_target.clamp(r_min, r_max);

    // Find the segment that brackets r_target in the y (radial) dimension
    // of the characteristic line.
    for w in line.windows(2) {
        let p1 = &w[0];
        let p2 = &w[1];

        let seg_min = p1.y.min(p2.y);
        let seg_max = p1.y.max(p2.y);

        if r_target >= seg_min - 1e-8 && r_target <= seg_max + 1e-8 {
            let mut tmp = MocPoint::default();
            tmp.x = line_x;
            tmp.y = r_target;
            let mut moc = tmp.interpolate_along(p1, p2);
            moc.x = line_x;
            moc.y = r_target;
            return Some(moc);
        }
    }

    // Fallback: use the nearest point on the line
    let nearest = line.iter().min_by(|a, b| {
        (a.y - r_target)
            .abs()
            .partial_cmp(&(b.y - r_target).abs())
            .unwrap_or(std::cmp::Ordering::Equal)
    })?;

    let mut moc = nearest.clone();
    moc.x = line_x;
    moc.y = r_target;
    Some(moc)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Material;

    fn make_pt(x: f64, y: f64, u: f64, v: f64) -> MocPoint {
        MocPoint::new(
            x,
            y,
            u,
            v,
            100000.0,
            300.0,
            1.2,
            Material::from_rgas_gamma(287.0, 1.4),
        )
    }

    #[test]
    fn test_trace_empty_fails() {
        let datasource = CharLines::new();
        let start = Point3d::new(0.0, 0.5, 0.0);
        assert!(trace_streamline(&start, &datasource, true).is_none());
    }

    #[test]
    fn test_trace_single_line() {
        let mut line = CharLine::new();
        line.push(make_pt(0.0, 0.0, 500.0, 0.0));
        line.push(make_pt(0.0, 2.0, 500.0, 0.0));
        let mut cl = CharLines::new();
        cl.push(line);

        let start = Point3d::new(0.0, 1.0, 0.0);
        let result = trace_streamline(&start, &cl, true);
        // With only one char line we only get the bootstrap point.
        assert!(result.is_some());
        let pts = result.unwrap();
        assert_eq!(pts.len(), 1);
        assert!((pts[0].x - 0.0).abs() < 1e-10);
        assert!((pts[0].y - 1.0).abs() < 1e-10);
    }
}
