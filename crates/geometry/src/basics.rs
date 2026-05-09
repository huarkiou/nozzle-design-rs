use std::ops::Neg;

use crate::point::Point3d;

/// Sign function: returns -1, 0, or 1 (as type `T`) depending on whether
/// `x` is negative, zero, or positive.
#[inline]
pub fn sgn<T>(x: T) -> T
where
    T: PartialOrd + From<i8> + Neg<Output = T>,
{
    if x < T::from(0) {
        T::from(-1)
    } else if x > T::from(0) {
        T::from(1)
    } else {
        T::from(0)
    }
}

/// Compute the signed area of a simple polygon using the Gaussian
/// shoelace formula.  `points` must be in order; the z-coordinate is
/// ignored.  Returns a positive value for CCW-oriented polygons.
pub fn polygon_area(points: &[Point3d]) -> f64 {
    let n = points.len();
    if n < 3 {
        return 0.0;
    }

    let mut area = 0.0;
    for i in 0..n {
        let j = (i + 1) % n;
        area += points[i].x * points[j].y - points[j].x * points[i].y;
    }
    area.abs() * 0.5
}

/// Compute the area-weighted centroid `(x, y)` of a simple polygon.
/// The z coordinate is ignored.  Returns `(0.0, 0.0)` for degenerate
/// (empty or zero-area) polygons.
pub fn polygon_centroid(points: &[Point3d]) -> (f64, f64) {
    let n = points.len();
    if n < 3 {
        return (0.0, 0.0);
    }

    let mut signed_area = 0.0;
    let mut cx = 0.0;
    let mut cy = 0.0;

    for i in 0..n {
        let j = (i + 1) % n;
        let cross = points[i].x * points[j].y - points[j].x * points[i].y;
        signed_area += cross;
        cx += (points[i].x + points[j].x) * cross;
        cy += (points[i].y + points[j].y) * cross;
    }

    if signed_area.abs() < f64::EPSILON {
        return (0.0, 0.0);
    }

    let factor = 1.0 / (3.0 * signed_area);
    (cx * factor, cy * factor)
}

/// Check whether the ray
///
/// ```text
/// origin + t_ray * (cos(theta), sin(theta)),   t_ray >= 0
/// ```
///
/// intersects the line segment `p1 → p2`.
///
/// Returns `Some((t, intersection_point))` where `t ∈ [0, 1]` is the
/// parameter along the segment, or `None` if there is no intersection.
pub fn ray_segment_intersection(
    origin: &Point3d,
    theta: f64,
    p1: &Point3d,
    p2: &Point3d,
) -> Option<(f64, Point3d)> {
    let cos_t = theta.cos();
    let sin_t = theta.sin();

    let dx = p2.x - p1.x;
    let dy = p2.y - p1.y;

    // Determinant of the 2×2 system
    let det = dx * sin_t - dy * cos_t;

    // Near-zero determinant → ray is (almost) parallel to segment
    if det.abs() < f64::EPSILON {
        return None;
    }

    let ox = origin.x;
    let oy = origin.y;

    // Solve for t (segment parameter) via Cramer's rule
    let t = (cos_t * (p1.y - oy) - sin_t * (p1.x - ox)) / det;

    if t < 0.0 || t > 1.0 {
        return None;
    }

    // Solve for t_ray (ray parameter)
    let t_ray = ((p1.y - oy) * dx - (p1.x - ox) * dy) / det;

    if t_ray < 0.0 {
        return None;
    }

    Some((t, Point3d::new(p1.x + t * dx, p1.y + t * dy, 0.0)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sgn() {
        assert_eq!(sgn(-3i32), -1);
        assert_eq!(sgn(0i32), 0);
        assert_eq!(sgn(5i32), 1);

        assert_eq!(sgn(-2.5f64), -1.0);
        assert_eq!(sgn(0.0f64), 0.0);
        assert_eq!(sgn(3.14f64), 1.0);
    }

    #[test]
    fn test_polygon_area() {
        // Unit square
        let square = vec![
            Point3d::new(0.0, 0.0, 0.0),
            Point3d::new(1.0, 0.0, 0.0),
            Point3d::new(1.0, 1.0, 0.0),
            Point3d::new(0.0, 1.0, 0.0),
        ];
        assert!((polygon_area(&square) - 1.0).abs() < 1e-12);

        // Triangle
        let triangle = vec![
            Point3d::new(0.0, 0.0, 0.0),
            Point3d::new(2.0, 0.0, 0.0),
            Point3d::new(0.0, 3.0, 0.0),
        ];
        assert!((polygon_area(&triangle) - 3.0).abs() < 1e-12);

        // Degenerate
        assert!((polygon_area(&[]) - 0.0).abs() < 1e-12);
        assert!((polygon_area(&[Point3d::new(0.0, 0.0, 0.0)]) - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_polygon_centroid() {
        let square = vec![
            Point3d::new(0.0, 0.0, 0.0),
            Point3d::new(2.0, 0.0, 0.0),
            Point3d::new(2.0, 2.0, 0.0),
            Point3d::new(0.0, 2.0, 0.0),
        ];
        let (cx, cy) = polygon_centroid(&square);
        assert!((cx - 1.0).abs() < 1e-12);
        assert!((cy - 1.0).abs() < 1e-12);

        // Degenerate
        assert_eq!(polygon_centroid(&[]), (0.0, 0.0));
    }

    #[test]
    fn test_ray_segment_intersection() {
        let origin = Point3d::new(0.0, 0.0, 0.0);
        let p1 = Point3d::new(1.0, -1.0, 0.0);
        let p2 = Point3d::new(1.0, 1.0, 0.0);

        // Ray pointing directly at (1, 0) — should hit middle
        let result = ray_segment_intersection(&origin, 0.0, &p1, &p2);
        assert!(result.is_some());
        let (t, pt) = result.unwrap();
        assert!((t - 0.5).abs() < 1e-12);
        assert!((pt.x - 1.0).abs() < 1e-12);
        assert!((pt.y - 0.0).abs() < 1e-12);

        // Ray pointing away
        let result = ray_segment_intersection(&origin, std::f64::consts::PI, &p1, &p2);
        assert!(result.is_none());

        // Ray parallel to segment
        let result = ray_segment_intersection(&origin, std::f64::consts::FRAC_PI_2, &p1, &p2);
        assert!(result.is_none());
    }
}
