use std::f64::consts::PI;

use crate::closed_curve::ClosedCurve;
use crate::point::Point3d;
use math::rootfinding::secant;
use math::Tolerance;

/// An elliptical cross-section with semi-major axis `a`, semi-minor axis
/// `b`, rotation `alpha` (radians), and center offset `(x0, y0)`.
#[derive(Debug, Clone)]
pub struct Ellipse {
    /// Semi-major axis length.
    pub a: f64,
    /// Semi-minor axis length.
    pub b: f64,
    /// Rotation angle in radians.
    pub alpha: f64,
    /// Center offset (x0, y0).
    pub center: (f64, f64),
}

impl Ellipse {
    /// Create a new `Ellipse`.
    pub fn new(a: f64, b: f64, alpha: f64, center: (f64, f64)) -> Self {
        Self {
            a,
            b,
            alpha,
            center,
        }
    }

    /// Tolerance for root-finding.
    fn tol() -> Tolerance {
        Tolerance::new(1e-8, 1e-8)
    }
}

impl ClosedCurve for Ellipse {
    /// Generate a point at polar angle `theta` using Newton's (secant) method
    /// to solve for the parameter `β` that yields the correct polar angle.
    fn generate_point(&self, theta: f64) -> Point3d {
        let (x0, y0) = self.center;
        let a = self.a;
        let b = self.b;
        let sin_a = self.alpha.sin();
        let cos_a = self.alpha.cos();

        // f(β) = atan2(Y(β), X(β)) - θ
        //   X(β) = a·cos(β)·cos(α) - b·sin(β)·sin(α)
        //   Y(β) = a·cos(β)·sin(α) + b·sin(β)·cos(α)
        let f = |beta: f64| -> f64 {
            let cos_b = beta.cos();
            let sin_b = beta.sin();
            let rx = a * cos_b;
            let ry = b * sin_b;
            let x = rx * cos_a - ry * sin_a;
            let y = rx * sin_a + ry * cos_a;
            y.atan2(x) - theta
        };

        // Use theta as initial guess
        let beta = secant::solve(f, theta, Self::tol(), 50).unwrap_or(theta);

        let cos_b = beta.cos();
        let sin_b = beta.sin();
        let rx = a * cos_b;
        let ry = b * sin_b;

        let x = x0 + rx * cos_a - ry * sin_a;
        let y = y0 + rx * sin_a + ry * cos_a;

        Point3d::new(x, y, 0.0)
    }

    /// Generate `n` points using the analytic parametric formula (no
    /// root-finding).
    ///
    /// Points are computed as:
    /// `x = x0 + a·cos(θ)·cos(α) - b·sin(θ)·sin(α)`
    /// `y = y0 + a·cos(θ)·sin(α) + b·sin(θ)·cos(α)`
    fn generate_points(&self, n: usize) -> Vec<Point3d> {
        let (x0, y0) = self.center;
        let a = self.a;
        let b = self.b;
        let sin_a = self.alpha.sin();
        let cos_a = self.alpha.cos();

        (0..n)
            .map(|i| {
                let theta = 2.0 * PI * i as f64 / n as f64;
                let cos_t = theta.cos();
                let sin_t = theta.sin();
                let x = x0 + a * cos_t * cos_a - b * sin_t * sin_a;
                let y = y0 + a * cos_t * sin_a + b * sin_t * cos_a;
                Point3d::new(x, y, 0.0)
            })
            .collect()
    }

    fn center(&self) -> (f64, f64) {
        self.center
    }

    fn clone_box(&self) -> Box<dyn ClosedCurve> {
        Box::new(self.clone())
    }
}
