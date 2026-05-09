use geometry::{Point3d, WallPoints};

/// Merge a downstream profile and an upstream profile into a single
/// wall profile using a weight function.
///
/// # Arguments
///
/// * `downstream` – Profile obtained by tracing from the inlet downstream.
/// * `upstream`   – Profile obtained by tracing from the outlet upstream
///   (reversed so it also goes inlet → outlet).
/// * `monotonic`  – If `true`, enforce that the output profile has
///   strictly increasing x‑coordinates.  Non‑monotonic points are
///   dropped.
/// * `weight_func` – A function `w ∈ [0, 1] → [0, 1]` that controls
///   the blend ratio.  `w(0) = 0` (fully downstream), `w(1) = 1`
///   (fully upstream).  The argument to `w` is the normalised
///   arc‑length position along the profile.
///
/// # Returns
///
/// The merged `WallPoints`.
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

    // ── build cumulative arc lengths for both profiles ────────────
    let ds_arc = cumulative_arc_lengths(downstream);
    let us_arc = cumulative_arc_lengths(upstream);

    let total_ds = *ds_arc.last().unwrap();
    let total_us = *us_arc.last().unwrap();

    // Use the longer profile as the "master" for resampling
    let master = if total_ds >= total_us {
        downstream
    } else {
        upstream
    };
    let master_arc = if total_ds >= total_us {
        &ds_arc
    } else {
        &us_arc
    };
    let master_total = total_ds.max(total_us);

    let n = master.len();
    let mut result = Vec::with_capacity(n);

    for i in 0..n {
        // Normalised arc-length position ∈ [0, 1]
        let s = master_arc[i] / master_total;
        let w = weight_func(s);

        // Find corresponding point on each profile at the same
        // arc-length fraction.
        let ds_pt = resample_at_arc_fraction(downstream, &ds_arc, s);
        let us_pt = resample_at_arc_fraction(upstream, &us_arc, s);

        let x = (1.0 - w) * ds_pt.x + w * us_pt.x;
        let y = (1.0 - w) * ds_pt.y + w * us_pt.y;
        let z = (1.0 - w) * ds_pt.z + w * us_pt.z;

        result.push(Point3d::new(x, y, z));
    }

    if monotonic {
        result = enforce_monotonic_x(result);
    }

    result
}

/// Compute cumulative arc lengths for a profile.
fn cumulative_arc_lengths(profile: &WallPoints) -> Vec<f64> {
    let mut cum = Vec::with_capacity(profile.len());
    cum.push(0.0);
    for w in profile.windows(2) {
        let d = w[0].distance_to(&w[1]);
        cum.push(cum.last().unwrap() + d);
    }
    cum
}

/// Resample a profile at the given arc-length fraction `s ∈ [0, 1]`.
fn resample_at_arc_fraction(profile: &WallPoints, cum_arc: &[f64], s: f64) -> Point3d {
    if profile.len() == 1 {
        return profile[0];
    }
    let total = *cum_arc.last().unwrap();
    if total < 1e-15 {
        return profile[0];
    }

    let target = s * total;

    for i in 1..profile.len() {
        if target <= cum_arc[i] {
            let seg_len = cum_arc[i] - cum_arc[i - 1];
            let t = if seg_len > 1e-15 {
                (target - cum_arc[i - 1]) / seg_len
            } else {
                0.5 // midpoint of zero-length segment
            };
            let a = &profile[i - 1];
            let b = &profile[i];
            return Point3d::new(
                a.x + t * (b.x - a.x),
                a.y + t * (b.y - a.y),
                a.z + t * (b.z - a.z),
            );
        }
    }

    profile[profile.len() - 1]
}

/// Remove points that violate strictly‑increasing‑x monotonicity.
fn enforce_monotonic_x(mut pts: WallPoints) -> WallPoints {
    if pts.len() < 2 {
        return pts;
    }

    // Single forward pass: keep the first point, then only keep
    // subsequent points whose x is strictly greater than the last kept
    // point.
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
        assert_eq!(result.len(), 3);
        // At x=0, weight=0 → use downstream
        assert!((result[0].y - 0.0).abs() < 1e-10);
        // At x=2, weight=1 → use upstream
        assert!((result[2].y - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_merge_monotonic() {
        // Non-monotonic downstream
        let downstream = vec![pt(0.0, 0.0), pt(0.5, 1.0), pt(0.3, 2.0)];
        let upstream = vec![pt(0.0, 0.0), pt(0.5, 1.0), pt(1.0, 2.0)];
        let weight = |_: f64| 0.5;

        let result = merge_line(&downstream, &upstream, true, &weight);
        for w in result.windows(2) {
            assert!(w[1].x > w[0].x, "x not monotonic");
        }
    }
}
