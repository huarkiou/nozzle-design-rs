mod intersect;
mod merge;
mod trace;
mod weight;

use crate::moc::CharLines;
use geometry::wallpoints::{interpolate_points_theta, refine_or_coarsen_boundaries};
use geometry::{ClosedCurve, Point3d, WallPoints};

/// A streamline paired with its polar angle for sorting.
#[derive(Debug, Clone)]
pub struct StreamlineWithAngle {
    pub streamline: WallPoints,
    /// Polar angle of the originating cross‑section point (radians).
    pub theta: f64,
}

/// Configuration for the streamline trace algorithm.
#[derive(Clone)]
pub struct StreamlineConfig {
    /// Whether the flow field is axisymmetric.
    pub axisymmetric: bool,
    /// Characteristic lines from the inlet region.
    pub datasource_inlet: CharLines,
    /// Characteristic lines from the outlet region.
    pub datasource_outlet: CharLines,
    /// Inlet cross‑section shape.
    pub inlet_shape: Box<dyn ClosedCurve>,
    /// Outlet cross‑section shape.
    pub outlet_shape: Box<dyn ClosedCurve>,
    /// Whether to enforce monotonic x in the merged profile.
    pub monotonic: bool,
    /// Weight‑function control parameter (see [`weight::build_weight_function`]).
    pub weight_parameter_a: f64,
}

impl std::fmt::Debug for StreamlineConfig {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("StreamlineConfig")
            .field("axisymmetric", &self.axisymmetric)
            .field(
                "datasource_inlet",
                &format_args!("{} lines", self.datasource_inlet.len()),
            )
            .field(
                "datasource_outlet",
                &format_args!("{} lines", self.datasource_outlet.len()),
            )
            .field("monotonic", &self.monotonic)
            .field("weight_parameter_a", &self.weight_parameter_a)
            .finish()
    }
}

/// Main streamline trace engine.
///
/// Given a MOC flow field (two sets of characteristic lines for the
/// inlet and outlet regions) and inlet/outlet shapes, traces 3D
/// streamlines and merges upstream/downstream profiles.
pub struct StreamlineTrace {
    config: StreamlineConfig,
    weight_func: Box<dyn Fn(f64) -> f64 + Send + Sync>,
    /// Final merged profiles, one per azimuthal position.
    pub model: Vec<WallPoints>,
    /// Downstream profiles (inlet → outlet), one per azimuthal position.
    pub downstream: Vec<WallPoints>,
    /// Upstream profiles (outlet → inlet), one per azimuthal position.
    pub upstream: Vec<WallPoints>,
}

impl StreamlineTrace {
    /// Create a new `StreamlineTrace` engine from the given
    /// configuration.
    pub fn new(config: StreamlineConfig) -> Self {
        let weight_func = weight::build_weight_function(config.weight_parameter_a);
        Self {
            config,
            weight_func,
            model: Vec::new(),
            downstream: Vec::new(),
            upstream: Vec::new(),
        }
    }

    /// Run the streamline trace.
    ///
    /// # Arguments
    ///
    /// * `n_theta` – Number of azimuthal (polar‑angle) samples.
    /// * `n_axis`  – Number of axial points in each output profile.
    pub fn run(&mut self, n_theta: usize, n_axis: usize) -> Result<(), String> {
        if n_theta < 2 {
            return Err("n_theta must be ≥ 2".into());
        }
        if n_axis < 2 {
            return Err("n_axis must be ≥ 2".into());
        }
        if self.config.datasource_inlet.is_empty() {
            return Err("datasource_inlet is empty".into());
        }
        if self.config.datasource_outlet.is_empty() {
            return Err("datasource_outlet is empty".into());
        }

        // ── Step 1: generate n_theta points on inlet & outlet shapes ──
        let inlet_points_raw = self.config.inlet_shape.generate_points(n_theta);
        let outlet_points_raw = self.config.outlet_shape.generate_points(n_theta);

        // ── Step 2: match by polar angle ────────────────────────────
        // Interpolate outlet points to the same theta values as inlet.
        let matched_outlet = interpolate_points_theta(&inlet_points_raw, &outlet_points_raw);

        // ── Step 3: downstream tracing (inlet → outlet) ────────────
        let mut downstream_profiles: Vec<WallPoints> = Vec::with_capacity(n_theta);
        for (_i, pt) in inlet_points_raw.iter().enumerate() {
            let r0 = pt.x.hypot(pt.y); // radial distance from symmetry axis
            match trace::trace_streamline(pt, &self.config.datasource_inlet, true) {
                Some(wp) => {
                    let wp3d = transform_to_3d(&wp, pt, r0);
                    downstream_profiles.push(wp3d);
                }
                None => {
                    return Err(format!(
                        "Downstream trace failed for inlet point ({:.4}, {:.4})",
                        pt.x, pt.y
                    ));
                }
            }
        }

        // ── Step 4: upstream tracing (outlet → inlet) ──────────────
        // Reverse the order of characteristic lines: outlet becomes
        // first, inlet becomes last.  Do NOT reverse individual line
        // point order nor negate velocities — the flow direction is
        // still from inlet→outlet; we rely on x_positive=false to
        // enforce backward marching.
        let mut reversed_outlet = CharLines::with_capacity(self.config.datasource_outlet.len());
        for line in self.config.datasource_outlet.iter().rev() {
            reversed_outlet.push(line.clone());
        }

        let mut upstream_profiles: Vec<WallPoints> = Vec::with_capacity(n_theta);
        for (_i, pt) in matched_outlet.iter().enumerate() {
            let r0 = pt.x.hypot(pt.y);
            match trace::trace_streamline(pt, &reversed_outlet, false) {
                Some(mut wp) => {
                    wp.reverse(); // inlet→outlet ordering
                    let wp3d = transform_to_3d(&wp, pt, r0);
                    upstream_profiles.push(wp3d);
                }
                None => {
                    return Err(format!(
                        "Upstream trace failed for outlet point ({:.4}, {:.4})",
                        pt.x, pt.y
                    ));
                }
            }
        }

        // ── Step 5: merge downstream + upstream ────────────────────
        let mut model_profiles: Vec<WallPoints> = Vec::with_capacity(n_theta);
        for (ds, us) in downstream_profiles.iter().zip(upstream_profiles.iter()) {
            let merged = merge::merge_line(ds, us, self.config.monotonic, &*self.weight_func);
            model_profiles.push(merged);
        }

        // ── Step 6: resample to n_axis points ──────────────────────
        let model_profiles = refine_or_coarsen_boundaries(&model_profiles, n_axis);

        // ── Step 7: sort by x‑coordinate of the first point ────────
        let mut sorted: Vec<WallPoints> = model_profiles;
        sorted.sort_by(|a, b| {
            let xa = a.first().map(|p| p.x).unwrap_or(f64::INFINITY);
            let xb = b.first().map(|p| p.x).unwrap_or(f64::INFINITY);
            xa.partial_cmp(&xb).unwrap_or(std::cmp::Ordering::Equal)
        });

        self.model = sorted;
        self.downstream = downstream_profiles;
        self.upstream = upstream_profiles;

        Ok(())
    }
}

/// Transform a 2D streamline (x = axial, y = radial from flow field)
/// into 3D coordinates using the start point's spatial position.
///
/// The start point `origin` (from the cross‑section shape) has
/// coordinates `(z_s, y_s)` = `(origin.x, origin.y)`.  Its distance
/// from the symmetry axis is `r0 = hypot(z_s, y_s)`.
///
/// For each flow‑field point `(x_f, r_f)` the 3D point is:
///
/// ```text
///   x_3d = x_f
///   y_3d = r_f * (z_s / r0)     (lateral / spatial)
///   z_3d = r_f * (y_s / r0)     (vertical)
/// ```
///
/// If `r0 ≈ 0` the point sits on the symmetry axis and the 3D
/// coordinates collapse to `(x_f, 0, 0)`.
fn transform_to_3d(streamline: &WallPoints, origin: &Point3d, r0: f64) -> WallPoints {
    if r0 < 1e-15 {
        return streamline
            .iter()
            .map(|p| Point3d::new(p.x, 0.0, 0.0))
            .collect();
    }
    let zs = origin.x; // spatial (shape z → 3D y)
    let ys = origin.y; // vertical  (shape y → 3D z)
    streamline
        .iter()
        .map(|p| {
            let rf = p.y; // radial position from flow field
            Point3d::new(p.x, rf * zs / r0, rf * ys / r0)
        })
        .collect()
}
