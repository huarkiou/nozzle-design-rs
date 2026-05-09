use crate::basics::ray_segment_intersection;
use crate::closed_curve::ClosedCurve;
use crate::point::Point3d;

/// A rectangular cross-section with length `l` (along first axis), width `w`
/// (along second axis), rotation `alpha` (radians), and center `(x0, y0)`.
#[derive(Debug, Clone)]
pub struct Rectangular {
    pub length: f64,
    pub width: f64,
    pub alpha: f64,
    pub center: (f64, f64),
}

impl Rectangular {
    /// Create a new `Rectangular`.
    ///
    /// Returns an error if `length ≤ 0` or `width ≤ 0`.
    pub fn new(
        length: f64,
        width: f64,
        alpha: f64,
        center: (f64, f64),
    ) -> Result<Self, &'static str> {
        if length <= 0.0 || width <= 0.0 {
            return Err("rectangular length and width must be > 0");
        }
        Ok(Self {
            length,
            width,
            alpha,
            center,
        })
    }

    fn corners(&self) -> [(f64, f64); 4] {
        let hx = self.length / 2.0;
        let hy = self.width / 2.0;
        let cos_a = self.alpha.cos();
        let sin_a = self.alpha.sin();
        let (x0, y0) = self.center;

        let rotate = |x: f64, y: f64| -> (f64, f64) {
            (x0 + x * cos_a - y * sin_a, y0 + x * sin_a + y * cos_a)
        };

        [
            rotate(-hx, hy),
            rotate(hx, hy),
            rotate(hx, -hy),
            rotate(-hx, -hy),
        ]
    }
}

impl ClosedCurve for Rectangular {
    fn generate_point(&self, theta: f64) -> Point3d {
        let corners = self.corners();
        let (x0, y0) = self.center;
        let origin = Point3d::new(x0, y0, 0.0);
        let mut best: Option<Point3d> = None;
        let mut best_dist = f64::INFINITY;

        for i in 0..4 {
            let j = (i + 1) % 4;
            let p1 = Point3d::new(corners[i].0, corners[i].1, 0.0);
            let p2 = Point3d::new(corners[j].0, corners[j].1, 0.0);
            if let Some((_, pt)) = ray_segment_intersection(&origin, theta, &p1, &p2) {
                let dx = pt.x - x0;
                let dy = pt.y - y0;
                let d2 = dx * dx + dy * dy;
                if d2 < best_dist {
                    best_dist = d2;
                    best = Some(pt);
                }
            }
        }

        best.unwrap_or_else(|| Point3d::new(x0, y0, 0.0))
    }

    fn generate_points(&self, n: usize) -> Vec<Point3d> {
        let n_adj = ((n.max(4) + 2) / 4 * 4).max(4);
        let per_edge = n_adj / 4;

        let corners = self.corners();
        let mut points = Vec::with_capacity(n_adj);

        for i in 0..4 {
            let j = (i + 1) % 4;
            let (ax, ay) = corners[i];
            let (bx, by) = corners[j];
            for k in 0..per_edge {
                let t = k as f64 / per_edge as f64;
                let x = ax + t * (bx - ax);
                let y = ay + t * (by - ay);
                points.push(Point3d::new(x, y, 0.0));
            }
        }

        points
    }

    fn center(&self) -> (f64, f64) {
        self.center
    }

    fn clone_box(&self) -> Box<dyn ClosedCurve> {
        Box::new(self.clone())
    }
}
