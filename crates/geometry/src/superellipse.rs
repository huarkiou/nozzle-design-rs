use std::f64::consts::PI;

use crate::basics::sgn;
use crate::closed_curve::ClosedCurve;
use crate::point::Point3d;
use math::Tolerance;
use math::rootfinding::secant;

/// A super-ellipse (Lamé curve) cross-section.
///
/// Defined by: `|x/a|^n + |y/b|^n = 1` where `power` corresponds to `n` in
/// the standard notation.
///
/// - `power = 2` → standard ellipse
/// - `power > 2` → rounded rectangle (squircle)
/// - `power < 2` → diamond-like (astroid)
#[derive(Debug, Clone)]
pub struct SuperEllipse {
    /// Semi-axis a.
    pub a: f64,
    /// Semi-axis b.
    pub b: f64,
    /// Super-ellipse power parameter (`n` in the standard notation).
    pub power: f64,
    /// Rotation angle in radians.
    pub alpha: f64,
    /// Center offset (x0, y0).
    pub center: (f64, f64),
}

impl SuperEllipse {
    /// Create a new `SuperEllipse`.
    pub fn new(a: f64, b: f64, power: f64, alpha: f64, center: (f64, f64)) -> Self {
        Self {
            a,
            b,
            power,
            alpha,
            center,
        }
    }

    /// Tolerance for root-finding.
    fn tol() -> Tolerance {
        Tolerance::new(1e-8, 1e-8)
    }

    /// Compute the radial coordinates for the super-ellipse at parameter
    /// `beta`.
    ///
    /// `rx = a · sgn(cos(β)) · |cos(β)|^(2/power)`
    /// `ry = b · sgn(sin(β)) · |sin(β)|^(2/power)`
    #[inline]
    fn radii(&self, beta: f64) -> (f64, f64) {
        let cos_b = beta.cos();
        let sin_b = beta.sin();
        let exp = 2.0 / self.power;
        let rx = self.a * sgn(cos_b) * cos_b.abs().powf(exp);
        let ry = self.b * sgn(sin_b) * sin_b.abs().powf(exp);
        (rx, ry)
    }
}

impl ClosedCurve for SuperEllipse {
    /// Generate a point at polar angle `theta` using Newton's (secant) method
    /// to solve for the parameter `β` that yields the correct polar angle.
    fn generate_point(&self, theta: f64) -> Point3d {
        let (x0, y0) = self.center;
        let sin_a = self.alpha.sin();
        let cos_a = self.alpha.cos();

        // f(β) = atan2(Y(β), X(β)) - θ
        //   X(β) = rx·cos(α) - ry·sin(α)
        //   Y(β) = rx·sin(α) + ry·cos(α)
        let f = |beta: f64| -> f64 {
            let (rx, ry) = self.radii(beta);
            let x = rx * cos_a - ry * sin_a;
            let y = rx * sin_a + ry * cos_a;
            y.atan2(x) - theta
        };

        // Use theta as initial guess
        let beta = secant::solve(f, theta, Self::tol(), 50).unwrap_or(theta);

        let (rx, ry) = self.radii(beta);
        let x = x0 + rx * cos_a - ry * sin_a;
        let y = y0 + rx * sin_a + ry * cos_a;

        Point3d::new(x, y, 0.0)
    }

    /// Generate `n` points using the analytic super-ellipse formula (no
    /// root-finding).
    fn generate_points(&self, n: usize) -> Vec<Point3d> {
        let (x0, y0) = self.center;
        let sin_a = self.alpha.sin();
        let cos_a = self.alpha.cos();

        (0..n)
            .map(|i| {
                let theta = 2.0 * PI * i as f64 / n as f64;
                let (rx, ry) = self.radii(theta);
                let x = x0 + rx * cos_a - ry * sin_a;
                let y = y0 + rx * sin_a + ry * cos_a;
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
