use crate::moc::MocPoint;

/// Find where a streamline (emanating from previous point `p_prev`)
/// intersects the segment `line[i]..line[i+1]` of a characteristic line.
///
/// Uses an iterative predictor‑corrector scheme (max 30 iterations):
///
/// 1. Estimate the flow direction at the intersection as
///    `L14 = tan((θ_prev + θ_est) / 2)` (first iteration: `θ_est = θ_prev`).
/// 2. Compute the slope of the char‑line segment `L23`.
/// 3. Solve for the intersection `(x, y)` of the two lines.
/// 4. Interpolate the flow parameters at `(x, y)` from the segment
///    endpoints → obtain a new estimate of `θ_est`.
/// 5. Repeat until `(x, y)` converges (tol 1e-8) or the iteration limit
///    is reached.
///
/// # Arguments
///
/// * `p_prev` – The previous point on the streamline (must have valid
///   `u, v`).
/// * `line`  – The characteristic line (a slice of `MocPoint`s).
/// * `i`     – Index of the left endpoint of the target segment.
/// * `x_positive` – Tracing direction: `true` means the intersection's
///   x‑coordinate must be ≥ `p_prev.x` (downstream), `false` means it
///   must be ≤ `p_prev.x` (upstream).
///
/// # Returns
///
/// `Some(MocPoint)` at the converged intersection, or `None` if the
/// segment is degenerate, the lines are parallel, or the iteration
/// fails to converge.
pub fn calculate_streamline_intersection(
    p_prev: &MocPoint,
    line: &[MocPoint],
    i: usize,
    x_positive: bool,
) -> Option<MocPoint> {
    // ── guard: valid segment ──────────────────────────────────────
    if i + 1 >= line.len() {
        return None;
    }
    let p_a = &line[i];
    let p_b = &line[i + 1];

    let dx_line = p_b.x - p_a.x;
    let dy_line = p_b.y - p_a.y;

    // Degenerate segment (collapsed endpoints)
    if dx_line.abs() < 1e-14 && dy_line.abs() < 1e-14 {
        return None;
    }

    // ── initial slope estimate ────────────────────────────────────
    let theta_prev = p_prev.flow_direction();
    let mut l14_old = theta_prev.tan();

    let mut x_old = f64::NAN;
    let mut y_old = f64::NAN;
    let mut current_theta: f64;

    for _iter in 0..30 {
        // ── char‑line segment slope ───────────────────────────────
        let l23 = if dx_line.abs() < 1e-14 {
            // Vertical characteristic line
            f64::INFINITY
        } else {
            dy_line / dx_line
        };

        // The two lines are:
        //   streamline: y = l14 * (x - p_prev.x) + p_prev.y
        //   char seg:   y = l23 * (x - p_a.x) + p_a.y
        //
        // Solve: l14*(x - px) + py = l23*(x - ax) + ay
        //   x*(l14 - l23) = l14*px - py - l23*ax + ay
        //   x = (l14*px - py - l23*ax + ay) / (l14 - l23)

        let (x_new, y_new) =
            if (l14_old - l23).abs() < 1e-14 || l23.is_infinite() && l14_old.is_infinite() {
                // Near‑parallel → no unique intersection
                return None;
            } else if dx_line.abs() < 1e-14 {
                // Vertical characteristic line — x known
                let x_int = p_a.x;
                let y_int = l14_old * (x_int - p_prev.x) + p_prev.y;
                (x_int, y_int)
            } else {
                let x_int = (l14_old * p_prev.x - p_prev.y - l23 * p_a.x + p_a.y) / (l14_old - l23);
                let y_int = l14_old * (x_int - p_prev.x) + p_prev.y;
                (x_int, y_int)
            };

        // Guard against NaN
        if !x_new.is_finite() || !y_new.is_finite() {
            return None;
        }

        // ── convergence check ─────────────────────────────────────
        if _iter > 0 && (x_new - x_old).abs() < 1e-8 && (y_new - y_old).abs() < 1e-8 {
            // Converged — return final interpolated point
            let mut tmp = MocPoint::default();
            tmp.x = x_new;
            tmp.y = y_new;
            let mut result = tmp.interpolate_along(p_a, p_b);
            result.x = x_new;
            result.y = y_new;
            return Some(result);
        }

        x_old = x_new;
        y_old = y_new;

        // ── check intersection lies within segment bounds ──────────
        let eps = 1e-10;
        let (x_min, x_max) = if p_a.x <= p_b.x {
            (p_a.x - eps, p_b.x + eps)
        } else {
            (p_b.x - eps, p_a.x + eps)
        };
        let (y_min, y_max) = if p_a.y <= p_b.y {
            (p_a.y - eps, p_b.y + eps)
        } else {
            (p_b.y - eps, p_a.y + eps)
        };

        if x_new < x_min || x_new > x_max || y_new < y_min || y_new > y_max {
            return None;
        }

        // ── direction check ────────────────────────────────────────
        if x_positive && x_new < p_prev.x - eps {
            return None;
        }
        if !x_positive && x_new > p_prev.x + eps {
            return None;
        }

        // ── interpolate flow parameters at (x_new, y_new) ──────────
        let mut tmp = MocPoint::default();
        tmp.x = x_new;
        tmp.y = y_new;
        let mut tmp = tmp.interpolate_along(p_a, p_b);
        tmp.x = x_new;
        tmp.y = y_new;
        current_theta = tmp.flow_direction();

        // ── update slope for next iteration ─────────────────────────
        l14_old = ((theta_prev + current_theta) / 2.0).tan();
    }

    // Exceeded max iterations without converging
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Material;

    fn make_point(x: f64, y: f64, u: f64, v: f64) -> MocPoint {
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
    fn test_simple_intersection() {
        // Previous point: (0, 0.5), flow direction ≈ 0 (horizontal right)
        let p_prev = make_point(0.0, 0.5, 500.0, 0.0);

        // Characteristic line segment: (2, 0) → (2, 2), vertical
        let line = vec![
            make_point(2.0, 0.0, 400.0, 50.0),
            make_point(2.0, 2.0, 400.0, 50.0),
        ];

        let result = calculate_streamline_intersection(&p_prev, &line, 0, true);
        assert!(result.is_some());
        let pt = result.unwrap();
        assert!((pt.x - 2.0).abs() < 1e-6);
        assert!((pt.y - 0.5).abs() < 0.2); // close to prev y
    }

    #[test]
    fn test_opposite_direction_rejected() {
        // p_prev is at x=5, but we require x_positive=true
        let p_prev = make_point(5.0, 1.0, 500.0, 0.0);
        let line = vec![
            make_point(2.0, 0.0, 400.0, 50.0),
            make_point(2.0, 2.0, 400.0, 50.0),
        ];
        let result = calculate_streamline_intersection(&p_prev, &line, 0, true);
        assert!(result.is_none());
    }

    #[test]
    fn test_parallel_lines() {
        let p_prev = make_point(0.0, 0.0, 500.0, 0.0);
        // Characteristic line parallel to flow direction
        let line = vec![
            make_point(1.0, 0.0, 400.0, 50.0),
            make_point(2.0, 0.0, 400.0, 50.0),
        ];
        let result = calculate_streamline_intersection(&p_prev, &line, 0, true);
        // Should be None because lines are nearly parallel (flow horizontal, line horizontal)
        assert!(result.is_none());
    }
}
