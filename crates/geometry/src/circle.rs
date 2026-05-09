use crate::closed_curve::ClosedCurve;
use crate::point::Point3d;

/// A circular cross-section.
#[derive(Debug, Clone)]
pub struct Circle {
    /// Radius of the circle.
    pub radius: f64,
    /// Center offset (x0, y0).
    pub center: (f64, f64),
}

impl Circle {
    /// Create a new `Circle` with the given radius and center offset.
    pub fn new(radius: f64, center: (f64, f64)) -> Self {
        Self { radius, center }
    }
}

impl ClosedCurve for Circle {
    /// Generate a point at polar angle `theta` using standard polar
    /// coordinates: `x = x0 + r·cos(θ)`, `y = y0 + r·sin(θ)`.
    fn generate_point(&self, theta: f64) -> Point3d {
        let (x0, y0) = self.center;
        Point3d::new(
            x0 + self.radius * theta.cos(),
            y0 + self.radius * theta.sin(),
            0.0,
        )
    }

    fn center(&self) -> (f64, f64) {
        self.center
    }

    fn clone_box(&self) -> Box<dyn ClosedCurve> {
        Box::new(self.clone())
    }
}
