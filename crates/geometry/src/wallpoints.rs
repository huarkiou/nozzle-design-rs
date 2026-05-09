//! Utilities for processing `WallPoints` (ordered collections of 3D wall points).
//!
//! Functions for interpolation, resampling, merging, polar-angle matching,
//! file I/O (UG NX `.dat` format), and basic transformations.

use crate::point::{Point3d, WallPoints};

// ---------------------------------------------------------------------------
// 1. Interpolation at a specific x coordinate
// ---------------------------------------------------------------------------

/// Linear interpolation at a specific x coordinate between two points.
///
/// Returns `None` if `x` is outside `[p1.x, p2.x]` (within a small
/// tolerance) or if the two points have the same x coordinate.
pub fn interpolate_point_x(p1: &Point3d, p2: &Point3d, x: f64) -> Option<Point3d> {
    if !x.is_finite() || !p1.x.is_finite() || !p2.x.is_finite() {
        return None;
    }

    let (x_min, x_max) = if p1.x <= p2.x {
        (p1.x, p2.x)
    } else {
        (p2.x, p1.x)
    };

    let tol = 1e-12;
    if x < x_min - tol || x > x_max + tol {
        return None;
    }

    let dx = p2.x - p1.x;
    if dx.abs() < f64::EPSILON {
        return None;
    }

    let t = (x - p1.x) / dx;
    Some(Point3d::new(
        x,
        p1.y + t * (p2.y - p1.y),
        p1.z + t * (p2.z - p1.z),
    ))
}

// ---------------------------------------------------------------------------
// 2. Merge two x-ordered arrays
// ---------------------------------------------------------------------------

/// Merge two x-ordered `WallPoints` arrays (ascending by x).
///
/// For x values present in both arrays the y/z values are averaged.
/// Duplicate x values (within `eps = 1e-8`) are removed.  The result is
/// sorted by x ascending.
///
/// This is typically used to merge downstream/upstream streamline profiles.
pub fn insert_into_xorderd_array(arr1: &[Point3d], arr2: &[Point3d]) -> WallPoints {
    let eps = 1e-8;

    if arr1.is_empty() {
        return arr2.to_vec();
    }
    if arr2.is_empty() {
        return arr1.to_vec();
    }

    let mut result = WallPoints::with_capacity(arr1.len() + arr2.len());
    let (mut i, mut j) = (0usize, 0usize);
    let n1 = arr1.len();
    let n2 = arr2.len();

    while i < n1 || j < n2 {
        if i >= n1 {
            result.push(arr2[j]);
            j += 1;
        } else if j >= n2 {
            result.push(arr1[i]);
            i += 1;
        } else {
            let dx = arr1[i].x - arr2[j].x;
            if dx.abs() < eps {
                let x = (arr1[i].x + arr2[j].x) * 0.5;
                let y = (arr1[i].y + arr2[j].y) * 0.5;
                let z = (arr1[i].z + arr2[j].z) * 0.5;
                result.push(Point3d::new(x, y, z));
                i += 1;
                j += 1;
            } else if dx < 0.0 {
                result.push(arr1[i]);
                i += 1;
            } else {
                result.push(arr2[j]);
                j += 1;
            }
        }
    }

    result
}

// ---------------------------------------------------------------------------
// 3. Resample a single boundary to `n` evenly-spaced points
// ---------------------------------------------------------------------------

/// Resample a single x-ordered `WallPoints` to exactly `n` points with
/// evenly-spaced x coordinates between `min_x` and `max_x`.
///
/// Uses linear interpolation between existing points.  If the boundary has
/// fewer than 2 points, the original slice is returned as-is.
pub fn refine_or_coarsen_boundary(boundary: &[Point3d], n: usize) -> WallPoints {
    if boundary.len() < 2 || n < 2 {
        return boundary.to_vec();
    }

    let mut min_x = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    for p in boundary {
        if p.x.is_finite() {
            min_x = min_x.min(p.x);
            max_x = max_x.max(p.x);
        }
    }

    if !min_x.is_finite() || !max_x.is_finite() {
        return boundary.to_vec();
    }

    let range = max_x - min_x;
    if range.abs() < f64::EPSILON {
        return boundary.to_vec();
    }

    let mut result = WallPoints::with_capacity(n);
    let n_f = (n - 1) as f64;

    for i in 0..n {
        let frac = i as f64 / n_f;
        let x = min_x + frac * range;

        let idx = boundary.partition_point(|p| p.x < x);

        let point = if idx == 0 {
            Point3d::new(x, boundary[0].y, boundary[0].z)
        } else if idx >= boundary.len() {
            let last = &boundary[boundary.len() - 1];
            Point3d::new(x, last.y, last.z)
        } else {
            let p0 = &boundary[idx - 1];
            let p1 = &boundary[idx];
            let dx = p1.x - p0.x;
            if dx.abs() < f64::EPSILON {
                Point3d::new(x, (p0.y + p1.y) * 0.5, (p0.z + p1.z) * 0.5)
            } else {
                let t = (x - p0.x) / dx;
                Point3d::new(x, p0.y + t * (p1.y - p0.y), p0.z + t * (p1.z - p0.z))
            }
        };
        result.push(point);
    }

    result
}

// ---------------------------------------------------------------------------
// 4. Resample multiple boundaries
// ---------------------------------------------------------------------------

/// Apply [`refine_or_coarsen_boundary`] to each boundary in the collection.
pub fn refine_or_coarsen_boundaries(boundaries: &[WallPoints], n: usize) -> Vec<WallPoints> {
    boundaries
        .iter()
        .map(|b| refine_or_coarsen_boundary(b, n))
        .collect()
}

// ---------------------------------------------------------------------------
// 5. Match points by polar angle
// ---------------------------------------------------------------------------

/// Match two sets of points by polar angle.
///
/// Given `base` (reference) and `target`, for each point in `base` the
/// function finds the point in `target` whose polar angle is closest.  If
/// no exact match exists within `eps = 1e-6`, the result is linearly
/// interpolated between the two nearest angles in `target`.
///
/// The output has the same length as `base` and contains the matched (or
/// interpolated) points from `target`.  Used to ensure 1:1 correspondence
/// between inlet and outlet cross-sections.
pub fn interpolate_points_theta(base: &[Point3d], target: &[Point3d]) -> WallPoints {
    if base.is_empty() || target.is_empty() {
        return WallPoints::new();
    }

    let eps = 1e-6;

    let mut indices: Vec<usize> = (0..target.len()).collect();
    indices.sort_by(|&a, &b| {
        target[a]
            .polar_angle()
            .partial_cmp(&target[b].polar_angle())
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let thetas: Vec<f64> = indices.iter().map(|&i| target[i].polar_angle()).collect();
    let n = thetas.len();

    let mut result = WallPoints::with_capacity(base.len());

    for bp in base {
        let theta = bp.polar_angle();
        if !theta.is_finite() {
            continue;
        }

        let pos = thetas.partition_point(|&t| t < theta - eps);

        if pos < n && (thetas[pos] - theta).abs() < eps {
            result.push(target[indices[pos]]);
        } else if pos > 0 && (thetas[pos - 1] - theta).abs() < eps {
            result.push(target[indices[pos - 1]]);
        } else {
            let idx_lo = if pos == 0 { n - 1 } else { pos - 1 };
            let idx_hi = pos % n;

            let p_lo = &target[indices[idx_lo]];
            let p_hi = &target[indices[idx_hi]];

            let t_lo = thetas[idx_lo];
            let t_hi = thetas[idx_hi];

            let (t_start, t_end) = if idx_lo == n - 1 && idx_hi == 0 {
                if theta >= t_lo {
                    (t_lo, t_hi + std::f64::consts::TAU)
                } else {
                    (t_lo - std::f64::consts::TAU, t_hi)
                }
            } else {
                (t_lo, t_hi)
            };

            let gap = t_end - t_start;
            let frac = if gap.abs() < f64::EPSILON {
                0.5
            } else {
                ((theta - t_start) / gap).clamp(0.0, 1.0)
            };

            let x = p_lo.x + frac * (p_hi.x - p_lo.x);
            let y = p_lo.y + frac * (p_hi.y - p_lo.y);
            let z = p_lo.z + frac * (p_hi.z - p_lo.z);

            result.push(Point3d::new(x, y, z));
        }
    }

    result
}

// ---------------------------------------------------------------------------
// 6. Merge two theta-ordered arrays
// ---------------------------------------------------------------------------

/// Merge two `WallPoints` arrays sorted by polar angle.
///
/// Duplicate theta values (within `eps = 1e-6`) are averaged.  Returns a
/// merged array sorted by theta ascending in `[0, 2π)`.
pub fn insert_into_thetaorderd_array(arr1: &[Point3d], arr2: &[Point3d]) -> WallPoints {
    let eps = 1e-6;

    let mut all: Vec<Point3d> = Vec::with_capacity(arr1.len() + arr2.len());
    all.extend(arr1.iter().copied());
    all.extend(arr2.iter().copied());

    all.sort_by(|a, b| {
        a.polar_angle()
            .partial_cmp(&b.polar_angle())
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let mut result = WallPoints::with_capacity(all.len());
    let mut i = 0;
    while i < all.len() {
        let mut sum_x = all[i].x;
        let mut sum_y = all[i].y;
        let mut sum_z = all[i].z;
        let mut count: usize = 1;
        let theta_ref = all[i].polar_angle();

        let mut j = i + 1;
        while j < all.len() && (all[j].polar_angle() - theta_ref).abs() < eps {
            sum_x += all[j].x;
            sum_y += all[j].y;
            sum_z += all[j].z;
            count += 1;
            j += 1;
        }

        if count == 1 {
            result.push(all[i]);
        } else {
            let inv = 1.0 / count as f64;
            result.push(Point3d::new(sum_x * inv, sum_y * inv, sum_z * inv));
        }
        i = j;
    }

    result
}

// ---------------------------------------------------------------------------
// 7. File I/O — UG NX .dat format
// ---------------------------------------------------------------------------

/// Write a single boundary to a UG NX `.dat` file.
///
/// Each point is written as `"x y z\n"` with 9 decimal places.
pub fn write_boundary(path: &str, boundary: &[Point3d]) -> std::io::Result<()> {
    let mut out = String::with_capacity(boundary.len() * 64);
    for p in boundary {
        use std::fmt::Write;
        let _ = writeln!(out, "{:.9} {:.9} {:.9}", p.x, p.y, p.z);
    }
    std::fs::write(path, out)
}

/// Write multiple boundaries to a UG NX `.dat` file.
///
/// Boundaries are separated by a `ROW` line.
pub fn write_boundaries(path: &str, boundaries: &[WallPoints]) -> std::io::Result<()> {
    let mut out = String::new();
    for (idx, boundary) in boundaries.iter().enumerate() {
        for p in boundary {
            use std::fmt::Write;
            let _ = writeln!(out, "{:.9} {:.9} {:.9}", p.x, p.y, p.z);
        }
        if idx + 1 < boundaries.len() {
            out.push_str("ROW\n");
        }
    }
    std::fs::write(path, out)
}

/// Read 2D points from a text file (one `x y` pair per line, space-separated).
pub fn read_points_2d(path: &str) -> std::io::Result<Vec<(f64, f64)>> {
    let content = std::fs::read_to_string(path)?;
    let mut points = Vec::new();
    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() || line == "ROW" {
            continue;
        }
        let mut parts = line.split_whitespace();
        let x: f64 = parts.next().and_then(|s| s.parse().ok()).ok_or_else(|| {
            std::io::Error::new(std::io::ErrorKind::InvalidData, "invalid x coordinate")
        })?;
        let y: f64 = parts.next().and_then(|s| s.parse().ok()).ok_or_else(|| {
            std::io::Error::new(std::io::ErrorKind::InvalidData, "invalid y coordinate")
        })?;
        points.push((x, y));
    }
    Ok(points)
}

/// Read multiple boundaries from a `.dat` file (blocks separated by `ROW`).
pub fn read_boundaries(path: &str) -> std::io::Result<Vec<WallPoints>> {
    let content = std::fs::read_to_string(path)?;
    let mut boundaries: Vec<WallPoints> = Vec::new();
    let mut current = WallPoints::new();

    for line in content.lines() {
        let line = line.trim();
        if line == "ROW" {
            if !current.is_empty() {
                boundaries.push(std::mem::take(&mut current));
            }
            continue;
        }
        if line.is_empty() {
            continue;
        }
        let mut parts = line.split_whitespace();
        let x: f64 = parts.next().and_then(|s| s.parse().ok()).ok_or_else(|| {
            std::io::Error::new(std::io::ErrorKind::InvalidData, "invalid x coordinate")
        })?;
        let y: f64 = parts.next().and_then(|s| s.parse().ok()).ok_or_else(|| {
            std::io::Error::new(std::io::ErrorKind::InvalidData, "invalid y coordinate")
        })?;
        let z: f64 = parts.next().and_then(|s| s.parse().ok()).unwrap_or(0.0);
        current.push(Point3d::new(x, y, z));
    }

    if !current.is_empty() {
        boundaries.push(current);
    }

    Ok(boundaries)
}

// ---------------------------------------------------------------------------
// 8. Utility functions
// ---------------------------------------------------------------------------

/// Compute the arithmetic mean center of a set of points (x, y only).
///
/// Returns `(0.0, 0.0)` for an empty slice.
pub fn cal_center_point(points: &[Point3d]) -> (f64, f64) {
    if points.is_empty() {
        return (0.0, 0.0);
    }
    let n = points.len() as f64;
    let sum_x: f64 = points.iter().map(|p| p.x).sum();
    let sum_y: f64 = points.iter().map(|p| p.y).sum();
    (sum_x / n, sum_y / n)
}

/// Rotate all points around the z-axis through the origin by `alpha` radians.
///
/// Rotation is counter-clockwise in the XY plane (right-hand rule with z
/// pointing up).
pub fn rotate_wallpoints(points: &mut [Point3d], alpha: f64) {
    let cos_a = alpha.cos();
    let sin_a = alpha.sin();
    for p in points.iter_mut() {
        let x_new = p.x * cos_a - p.y * sin_a;
        let y_new = p.x * sin_a + p.y * cos_a;
        p.x = x_new;
        p.y = y_new;
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn interpolate_x_basic() {
        let a = Point3d::new(0.0, 0.0, 0.0);
        let b = Point3d::new(2.0, 4.0, 6.0);
        let mid = interpolate_point_x(&a, &b, 1.0).unwrap();
        assert!((mid.x - 1.0).abs() < 1e-12);
        assert!((mid.y - 2.0).abs() < 1e-12);
        assert!((mid.z - 3.0).abs() < 1e-12);
    }

    #[test]
    fn interpolate_x_out_of_range() {
        let a = Point3d::new(0.0, 0.0, 0.0);
        let b = Point3d::new(2.0, 4.0, 6.0);
        assert!(interpolate_point_x(&a, &b, -0.1).is_none());
        assert!(interpolate_point_x(&a, &b, 2.1).is_none());
    }

    #[test]
    fn interpolate_x_same_x() {
        let a = Point3d::new(1.0, 0.0, 0.0);
        let b = Point3d::new(1.0, 4.0, 6.0);
        assert!(interpolate_point_x(&a, &b, 1.0).is_none());
    }

    #[test]
    fn interpolate_x_nan() {
        let a = Point3d::new(0.0, 0.0, 0.0);
        let b = Point3d::new(2.0, 4.0, 6.0);
        assert!(interpolate_point_x(&a, &b, f64::NAN).is_none());
    }

    #[test]
    fn interpolate_x_reversed() {
        let a = Point3d::new(2.0, 4.0, 6.0);
        let b = Point3d::new(0.0, 0.0, 0.0);
        let mid = interpolate_point_x(&a, &b, 1.0).unwrap();
        assert!((mid.x - 1.0).abs() < 1e-12);
        assert!((mid.y - 2.0).abs() < 1e-12);
    }

    #[test]
    fn merge_x_ordered_basic() {
        let arr1 = vec![
            Point3d::new(0.0, 0.0, 0.0),
            Point3d::new(2.0, 2.0, 2.0),
            Point3d::new(4.0, 4.0, 4.0),
        ];
        let arr2 = vec![Point3d::new(1.0, 10.0, 10.0), Point3d::new(3.0, 30.0, 30.0)];
        let merged = insert_into_xorderd_array(&arr1, &arr2);
        assert_eq!(merged.len(), 5);
        for w in merged.windows(2) {
            assert!(w[0].x <= w[1].x);
        }
    }

    #[test]
    fn merge_x_ordered_duplicate() {
        let arr1 = vec![Point3d::new(0.0, 0.0, 10.0), Point3d::new(2.0, 2.0, 20.0)];
        let arr2 = vec![Point3d::new(2.0, 8.0, 80.0), Point3d::new(4.0, 4.0, 40.0)];
        let merged = insert_into_xorderd_array(&arr1, &arr2);
        assert_eq!(merged.len(), 3);
        let p = &merged[1];
        assert!((p.x - 2.0).abs() < 1e-12);
        assert!((p.y - 5.0).abs() < 1e-12);
        assert!((p.z - 50.0).abs() < 1e-12);
    }

    #[test]
    fn merge_x_ordered_empty() {
        let arr1: WallPoints = vec![];
        let arr2 = vec![Point3d::new(1.0, 1.0, 1.0)];
        let merged = insert_into_xorderd_array(&arr1, &arr2);
        assert_eq!(merged.len(), 1);
        assert!((merged[0].x - 1.0).abs() < 1e-12);
    }

    #[test]
    fn refine_boundary_less_than_2() {
        let b = vec![Point3d::new(1.0, 2.0, 3.0)];
        let r = refine_or_coarsen_boundary(&b, 10);
        assert_eq!(r.len(), 1);
    }

    #[test]
    fn refine_boundary_exact_n() {
        let b = vec![
            Point3d::new(0.0, 0.0, 0.0),
            Point3d::new(10.0, 100.0, 200.0),
        ];
        let r = refine_or_coarsen_boundary(&b, 5);
        assert_eq!(r.len(), 5);
        assert!((r[0].x - 0.0).abs() < 1e-12);
        assert!((r[4].x - 10.0).abs() < 1e-12);
        assert!((r[2].x - 5.0).abs() < 1e-12);
        assert!((r[2].y - 50.0).abs() < 1e-12);
    }

    #[test]
    fn refine_boundary_same_x() {
        let b = vec![Point3d::new(5.0, 1.0, 2.0), Point3d::new(5.0, 3.0, 4.0)];
        let r = refine_or_coarsen_boundary(&b, 5);
        assert_eq!(r.len(), 2);
    }

    #[test]
    fn refine_boundaries_multiple() {
        let boundaries = vec![
            vec![Point3d::new(0.0, 0.0, 0.0), Point3d::new(2.0, 2.0, 2.0)],
            vec![Point3d::new(0.0, 0.0, 0.0), Point3d::new(4.0, 4.0, 4.0)],
        ];
        let r = refine_or_coarsen_boundaries(&boundaries, 3);
        assert_eq!(r.len(), 2);
        assert_eq!(r[0].len(), 3);
        assert_eq!(r[1].len(), 3);
    }

    #[test]
    fn theta_match_exact() {
        let base = vec![Point3d::new(1.0, 0.0, 0.0), Point3d::new(0.0, 1.0, 0.0)];
        let target = vec![Point3d::new(10.0, 0.0, 0.0), Point3d::new(0.0, 10.0, 0.0)];
        let r = interpolate_points_theta(&base, &target);
        assert_eq!(r.len(), 2);
        assert!((r[0].x - 10.0).abs() < 1e-12);
        assert!((r[1].y - 10.0).abs() < 1e-12);
    }

    #[test]
    fn theta_match_interpolate() {
        let base = vec![Point3d::new(1.0, 0.0, 0.0)];
        let target = vec![Point3d::new(10.0, 0.1, 0.0), Point3d::new(9.9, -0.1, 0.0)];
        let r = interpolate_points_theta(&base, &target);
        assert_eq!(r.len(), 1);
        assert!(r[0].x >= 9.9 && r[0].x <= 10.0);
    }

    #[test]
    fn theta_match_empty() {
        let base: WallPoints = vec![];
        let target = vec![Point3d::new(1.0, 0.0, 0.0)];
        assert!(interpolate_points_theta(&base, &target).is_empty());
    }

    #[test]
    fn merge_theta_ordered_basic() {
        let arr1 = vec![Point3d::new(1.0, 0.0, 0.0), Point3d::new(0.0, 1.0, 0.0)];
        let arr2 = vec![Point3d::new(-1.0, 0.0, 0.0), Point3d::new(0.0, -1.0, 0.0)];
        let merged = insert_into_thetaorderd_array(&arr1, &arr2);
        assert_eq!(merged.len(), 4);
        let thetas: Vec<f64> = merged.iter().map(|p| p.polar_angle()).collect();
        for w in thetas.windows(2) {
            assert!(w[0] <= w[1]);
        }
    }

    #[test]
    fn merge_theta_ordered_duplicate() {
        let arr1 = vec![Point3d::new(1.0, 0.0, 10.0)];
        let arr2 = vec![Point3d::new(2.0, 0.0, 20.0)];
        let merged = insert_into_thetaorderd_array(&arr1, &arr2);
        assert_eq!(merged.len(), 1);
        assert!((merged[0].x - 1.5).abs() < 1e-12);
        assert!((merged[0].z - 15.0).abs() < 1e-12);
    }

    #[test]
    fn write_and_read_boundary() {
        let boundary = vec![Point3d::new(1.0, 2.0, 3.0), Point3d::new(4.0, 5.0, 6.0)];
        let path = "test_boundary.dat";
        write_boundary(path, &boundary).unwrap();
        let read = read_boundaries(path).unwrap();
        assert_eq!(read.len(), 1);
        assert_eq!(read[0].len(), 2);
        assert!((read[0][0].x - 1.0).abs() < 1e-12);
        assert!((read[0][0].y - 2.0).abs() < 1e-12);
        assert!((read[0][0].z - 3.0).abs() < 1e-12);
        std::fs::remove_file(path).ok();
    }

    #[test]
    fn write_and_read_boundaries() {
        let boundaries = vec![
            vec![Point3d::new(1.0, 2.0, 3.0)],
            vec![Point3d::new(4.0, 5.0, 6.0)],
        ];
        let path = "test_boundaries.dat";
        write_boundaries(path, &boundaries).unwrap();
        let read = read_boundaries(path).unwrap();
        assert_eq!(read.len(), 2);
        assert_eq!(read[0].len(), 1);
        assert_eq!(read[1].len(), 1);
        std::fs::remove_file(path).ok();
    }

    #[test]
    fn read_points_2d_basic() {
        let path = "test_2d.dat";
        std::fs::write(path, "1.0 2.0\n3.0 4.0\n").unwrap();
        let pts = read_points_2d(path).unwrap();
        assert_eq!(pts.len(), 2);
        assert!((pts[0].0 - 1.0).abs() < 1e-12);
        assert!((pts[0].1 - 2.0).abs() < 1e-12);
        assert!((pts[1].0 - 3.0).abs() < 1e-12);
        assert!((pts[1].1 - 4.0).abs() < 1e-12);
        std::fs::remove_file(path).ok();
    }

    #[test]
    fn center_point_normal() {
        let pts = vec![Point3d::new(1.0, 2.0, 100.0), Point3d::new(3.0, 4.0, 200.0)];
        let (cx, cy) = cal_center_point(&pts);
        assert!((cx - 2.0).abs() < 1e-12);
        assert!((cy - 3.0).abs() < 1e-12);
    }

    #[test]
    fn center_point_empty() {
        let pts: WallPoints = vec![];
        let (cx, cy) = cal_center_point(&pts);
        assert_eq!(cx, 0.0);
        assert_eq!(cy, 0.0);
    }

    #[test]
    fn rotate_points() {
        let mut pts = vec![Point3d::new(1.0, 0.0, 5.0), Point3d::new(0.0, 1.0, 5.0)];
        rotate_wallpoints(&mut pts, std::f64::consts::FRAC_PI_2);
        assert!((pts[0].x - 0.0).abs() < 1e-12);
        assert!((pts[0].y - 1.0).abs() < 1e-12);
        assert!((pts[0].z - 5.0).abs() < 1e-12);
        assert!((pts[1].x + 1.0).abs() < 1e-12);
        assert!((pts[1].y - 0.0).abs() < 1e-12);
    }
}
