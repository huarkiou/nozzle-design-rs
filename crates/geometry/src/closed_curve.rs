use crate::point::Point3d;

/// A closed planar curve (in the XY plane) that can produce uniformly
/// sampled boundary points by polar angle.
///
/// The trait is object-safe so that it can be stored behind `dyn`
/// references in generic nozzle geometry code.
pub trait ClosedCurve: Send + Sync {
    /// Return the single boundary point at the given polar angle `theta`
    /// (in radians, in the range `[0, 2π)`).
    fn generate_point(&self, theta: f64) -> Point3d;

    /// Return `n` uniformly-sampled boundary points around `[0, 2π)`.
    fn generate_points(&self, n: usize) -> Vec<Point3d> {
        (0..n)
            .map(|i| {
                let theta = 2.0 * std::f64::consts::PI * i as f64 / n as f64;
                self.generate_point(theta)
            })
            .collect()
    }

    /// Center offset `(x0, y0)` of the curve.
    fn center(&self) -> (f64, f64);

    /// Clone the curve into a new heap allocation.
    fn clone_box(&self) -> Box<dyn ClosedCurve>;
}

/// Implement `Clone` for `Box<dyn ClosedCurve>` so that
/// `StreamlineConfig` can derive `Clone`.
impl Clone for Box<dyn ClosedCurve> {
    fn clone(&self) -> Self {
        self.clone_box()
    }
}
