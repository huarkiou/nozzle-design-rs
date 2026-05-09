use std::f64::consts::PI;

use crate::basics::ray_segment_intersection;
use crate::closed_curve::ClosedCurve;
use crate::point::Point3d;

/// A user-defined cross-section defined by an ordered list of boundary
/// points forming a closed polygon, plus rotation `alpha` and center offset.
#[derive(Debug, Clone)]
pub struct UserDefined {
    pub points: Vec<Point3d>,
    pub alpha: f64,
    pub center: (f64, f64),
}

impl UserDefined {
    pub fn new(points: Vec<Point3d>, alpha: f64, center_hint: Option<(f64, f64)>) -> Self {
        let center = match center_hint {
            Some(c) => c,
            None => {
                if points.is_empty() {
                    (0.0, 0.0)
                } else {
                    let n = points.len() as f64;
                    let (sx, sy) = points
                        .iter()
                        .fold((0.0, 0.0), |(ax, ay), p| (ax + p.x, ay + p.y));
                    (sx / n, sy / n)
                }
            }
        };
        Self {
            points,
            alpha,
            center,
        }
    }
}

impl ClosedCurve for UserDefined {
    fn generate_point(&self, theta: f64) -> Point3d {
        let (x0, y0) = self.center;
        let sin_a = self.alpha.sin();
        let cos_a = self.alpha.cos();

        if self.points.is_empty() {
            return Point3d::new(x0, y0, 0.0);
        }

        let local_theta = theta - self.alpha;
        let origin = Point3d::new(x0, y0, 0.0);

        let mut best: Option<Point3d> = None;
        let mut best_dist = f64::INFINITY;

        for i in 0..self.points.len() {
            let j = (i + 1) % self.points.len();
            let p1 = &self.points[i];
            let p2 = &self.points[j];

            if let Some((_, pt)) = ray_segment_intersection(&origin, local_theta, p1, p2) {
                // Rotate local-space intersection back to world space
                let dx = pt.x - x0;
                let dy = pt.y - y0;
                let wx = x0 + dx * cos_a - dy * sin_a;
                let wy = y0 + dx * sin_a + dy * cos_a;
                let d2 = (wx - x0) * (wx - x0) + (wy - y0) * (wy - y0);

                if d2 < best_dist {
                    best_dist = d2;
                    best = Some(Point3d::new(wx, wy, 0.0));
                }
            }
        }

        best.unwrap_or_else(|| Point3d::new(x0, y0, 0.0))
    }

    fn generate_points(&self, n: usize) -> Vec<Point3d> {
        let (x0, y0) = self.center;
        let sin_a = self.alpha.sin();
        let cos_a = self.alpha.cos();

        let mut candidates: Vec<Point3d> = Vec::with_capacity(n + self.points.len());

        for i in 0..n {
            let theta = 2.0 * PI * i as f64 / n as f64;
            candidates.push(self.generate_point(theta));
        }

        for pt in &self.points {
            let dx = pt.x - x0;
            let dy = pt.y - y0;
            let wx = x0 + dx * cos_a - dy * sin_a;
            let wy = y0 + dx * sin_a + dy * cos_a;
            candidates.push(Point3d::new(wx, wy, 0.0));
        }

        candidates.sort_by(|a, b| {
            let ta = (a.y - y0).atan2(a.x - x0);
            let tb = (b.y - y0).atan2(b.x - x0);
            ta.partial_cmp(&tb).unwrap_or(std::cmp::Ordering::Equal)
        });

        candidates
    }

    fn center(&self) -> (f64, f64) {
        self.center
    }

    fn clone_box(&self) -> Box<dyn ClosedCurve> {
        Box::new(self.clone())
    }
}
