use geometry::{Point3d, WallPoints};

/// Merge a downstream profile and an upstream profile into a single
/// wall profile using a weight function.
///
/// Uses x‑coordinate based merging (matching the C++ reference):
/// 1. Align both profiles to the same x‑coordinates via interpolation.
/// 2. Blend the y and z coordinates using `weight_func(norm_x)` where
///    `norm_x = (x - x_min) / (x_max - x_min)`.
/// 3. Optionally enforce monotonically increasing x.
pub fn merge_line(
    downstream: &WallPoints,
    upstream: &WallPoints,
    monotonic: bool,
    weight_func: &(dyn Fn(f64) -> f64 + Send + Sync),
) -> WallPoints {
    // If either profile is empty, return the other unchanged.
    if downstream.is_empty() {
        return upstream.to_vec();
    }
    if upstream.is_empty() {
        return downstream.to_vec();
    }

    // Get the global x range.
    let x_min = downstream
        .first()
        .map(|p| p.x)
        .unwrap_or(0.0)
        .min(upstream.first().map(|p| p.x).unwrap_or(0.0));
    let x_max = downstream
        .last()
        .map(|p| p.x)
        .unwrap_or(0.0)
        .max(upstream.last().map(|p| p.x).unwrap_or(0.0));
    let x_len = x_max - x_min;

    // Collect all unique x coordinates from both profiles
    let mut xs: Vec<f64> = Vec::with_capacity(downstream.len() + upstream.len());
    xs.extend(downstream.iter().map(|p| p.x));
    xs.extend(upstream.iter().map(|p| p.x));
    xs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    xs.dedup_by(|a, b| (*a - *b).abs() < 1e-8);

    let mut result: WallPoints = Vec::with_capacity(xs.len());

    for &x in &xs {
        let ds_pt = interpolate_at_x(downstream, x);
        let us_pt = interpolate_at_x(upstream, x);

        // Normalised position for weight function: 0 at inlet, 1 at outlet
        let w = if x_len > 1e-15 {
            weight_func((x - x_min) / x_len)
        } else {
            0.0
        };

        let y = (1.0 - w) * ds_pt.y + w * us_pt.y;
        let z = (1.0 - w) * ds_pt.z + w * us_pt.z;
        result.push(Point3d::new(x, y, z));
    }

    if monotonic {
        result = enforce_monotonic_x(result);
    }

    result
}

/// Linear interpolation at a given x between two bracketing profile points.
/// Uses binary search — the profile must be sorted by x (ascending).
fn interpolate_at_x(profile: &WallPoints, x: f64) -> Point3d {
    if profile.len() == 1 {
        return profile[0];
    }
    let idx = profile.partition_point(|p| p.x < x);
    if idx == 0 {
        return Point3d::new(x, profile[0].y, profile[0].z);
    }
    if idx >= profile.len() {
        let last = &profile[profile.len() - 1];
        return Point3d::new(x, last.y, last.z);
    }
    let a = &profile[idx - 1];
    let b = &profile[idx];
    let dx = b.x - a.x;
    let t = if dx.abs() > 1e-15 {
        (x - a.x) / dx
    } else {
        0.5
    };
    Point3d::new(x, a.y + t * (b.y - a.y), a.z + t * (b.z - a.z))
}

/// Remove points that violate strictly‑increasing‑x monotonicity.
fn enforce_monotonic_x(mut pts: WallPoints) -> WallPoints {
    if pts.len() < 2 {
        return pts;
    }
    let mut kept: WallPoints = Vec::with_capacity(pts.len());
    kept.push(pts[0]);
    for pt in pts.drain(1..) {
        if pt.x > kept.last().unwrap().x {
            kept.push(pt);
        }
    }
    kept
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pt(x: f64, y: f64) -> Point3d {
        Point3d::new(x, y, 0.0)
    }

    #[test]
    fn test_merge_linear_weight() {
        let downstream = vec![pt(0.0, 0.0), pt(1.0, 0.5), pt(2.0, 1.0)];
        let upstream = vec![pt(0.0, 1.0), pt(1.0, 1.5), pt(2.0, 2.0)];
        let weight = |x: f64| x; // linear

        let result = merge_line(&downstream, &upstream, false, &weight);
        // At x=0, weight=0 → use downstream → y ≈ 0
        assert!((result[0].y - 0.0).abs() < 1e-10);
        // At x=2, weight=1 → use upstream → y ≈ 2
        assert!((result.last().unwrap().y - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_merge_monotonic() {
        let downstream = vec![pt(0.0, 0.0), pt(0.5, 1.0), pt(0.3, 2.0)];
        let upstream = vec![pt(0.0, 0.0), pt(0.5, 1.0), pt(1.0, 2.0)];
        let weight = |_: f64| 0.5;

        let result = merge_line(&downstream, &upstream, true, &weight);
        for w in result.windows(2) {
            assert!(w[1].x > w[0].x, "x not monotonic");
        }
    }
}
