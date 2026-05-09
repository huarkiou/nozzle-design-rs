mod intersect;
mod merge;
mod trace;
mod weight;

use std::{
    fmt,
    sync::atomic::{AtomicUsize, Ordering},
};

use rayon::prelude::*;

use crate::moc::CharLines;
use geometry::wallpoints::{interpolate_points_theta, refine_or_coarsen_boundaries};
use geometry::{ClosedCurve, Point3d, WallPoints};

pub use trace::{CharLineSource, RevCharLines};

// ── error type ─────────────────────────────────────────────────

/// Errors that can occur during streamline tracing.
#[derive(Debug)]
pub enum StreamlineTraceError {
    /// Invalid parameter.
    InvalidParameter(String),
    /// One or more downstream traces failed.
    DownstreamTraceFailed {
        /// The inlet point (z, y) whose trace failed.
        point: (f64, f64),
        /// How many points failed.
        failed_count: usize,
        /// Total number of points.
        total: usize,
    },
    /// One or more upstream traces failed.
    UpstreamTraceFailed {
        /// The outlet point (z, y) whose trace failed.
        point: (f64, f64),
        failed_count: usize,
        total: usize,
    },
}

impl fmt::Display for StreamlineTraceError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidParameter(msg) => write!(f, "invalid parameter: {msg}"),
            Self::DownstreamTraceFailed {
                point,
                failed_count,
                total,
            } => write!(
                f,
                "downstream trace failed: {failed_count}/{total} streamlines — first at ({:.4}, {:.4})",
                point.0, point.1
            ),
            Self::UpstreamTraceFailed {
                point,
                failed_count,
                total,
            } => write!(
                f,
                "upstream trace failed: {failed_count}/{total} streamlines — first at ({:.4}, {:.4})",
                point.0, point.1
            ),
        }
    }
}

impl std::error::Error for StreamlineTraceError {}

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
    pub fn run(&mut self, n_theta: usize, n_axis: usize) -> Result<(), StreamlineTraceError> {
        if n_theta < 2 {
            return Err(StreamlineTraceError::InvalidParameter(
                "n_theta must be ≥ 2".into(),
            ));
        }
        if n_axis < 2 {
            return Err(StreamlineTraceError::InvalidParameter(
                "n_axis must be ≥ 2".into(),
            ));
        }
        if self.config.datasource_inlet.is_empty() {
            return Err(StreamlineTraceError::InvalidParameter(
                "datasource_inlet is empty".into(),
            ));
        }
        if self.config.datasource_outlet.is_empty() {
            return Err(StreamlineTraceError::InvalidParameter(
                "datasource_outlet is empty".into(),
            ));
        }

        // ── Step 1: generate n_theta points on inlet & outlet shapes ──
        let inlet_points_raw = self.config.inlet_shape.generate_points(n_theta);
        let outlet_points_raw = self.config.outlet_shape.generate_points(n_theta);

        // ── Step 2: match by polar angle ────────────────────────────
        let matched_outlet = interpolate_points_theta(&inlet_points_raw, &outlet_points_raw);

        // ── Step 3+4: parallel downstream & upstream tracing ─────────
        let rev_outlet = RevCharLines::new(&self.config.datasource_outlet);
        let axisymmetric = self.config.axisymmetric;
        let n_total = n_theta;

        // Collect failures instead of bailing on first error.
        let (downstream_results, upstream_results) = rayon::join(
            || -> Vec<Result<WallPoints, (f64, f64)>> {
                let counter = AtomicUsize::new(0);
                let results: Vec<_> = inlet_points_raw
                    .par_iter()
                    .map(|pt| {
                        let r0 = if axisymmetric { pt.x.hypot(pt.y) } else { pt.y };
                        let wp = trace::trace_streamline(
                            pt,
                            &self.config.datasource_inlet,
                            true,
                            axisymmetric,
                        )
                        .ok_or((pt.x, pt.y))?;
                        let n = counter.fetch_add(1, Ordering::Relaxed) + 1;
                        if n % 10 == 0 || n == n_total {
                            eprintln!("  downstream  {n}/{n_total}");
                        }
                        Ok(transform_to_3d(wp, pt, r0, axisymmetric))
                    })
                    .collect();
                eprintln!("  downstream  done");
                results
            },
            || -> Vec<Result<WallPoints, (f64, f64)>> {
                let counter = AtomicUsize::new(0);
                let results: Vec<_> = matched_outlet
                    .par_iter()
                    .map(|pt| {
                        let r0 = if axisymmetric { pt.x.hypot(pt.y) } else { pt.y };
                        let mut wp = trace::trace_streamline(pt, &rev_outlet, false, axisymmetric)
                            .ok_or((pt.x, pt.y))?;
                        wp.reverse();
                        let n = counter.fetch_add(1, Ordering::Relaxed) + 1;
                        if n % 10 == 0 || n == n_total {
                            eprintln!("  upstream    {n}/{n_total}");
                        }
                        Ok(transform_to_3d(wp, pt, r0, axisymmetric))
                    })
                    .collect();
                eprintln!("  upstream    done");
                results
            },
        );

        // Check downstream failures
        let ds_fails: Vec<_> = downstream_results.iter().filter(|r| r.is_err()).collect();
        if !ds_fails.is_empty() {
            let (z, y) = ds_fails[0].as_ref().unwrap_err();
            return Err(StreamlineTraceError::DownstreamTraceFailed {
                point: (*z, *y),
                failed_count: ds_fails.len(),
                total: downstream_results.len(),
            });
        }

        // Check upstream failures
        let us_fails: Vec<_> = upstream_results.iter().filter(|r| r.is_err()).collect();
        if !us_fails.is_empty() {
            let (z, y) = us_fails[0].as_ref().unwrap_err();
            return Err(StreamlineTraceError::UpstreamTraceFailed {
                point: (*z, *y),
                failed_count: us_fails.len(),
                total: upstream_results.len(),
            });
        }

        let downstream_profiles: Vec<WallPoints> =
            downstream_results.into_iter().map(|r| r.unwrap()).collect();
        let upstream_profiles: Vec<WallPoints> =
            upstream_results.into_iter().map(|r| r.unwrap()).collect();

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
/// The start point `origin` has coordinates `(z_s, y_s)`.
///
/// **Axisymmetric**: each flow-field radial position `r_f` is projected
/// onto the 3D sphere:
///   `y_3d = r_f × (z_s / r₀)`, `z_3d = r_f × (y_s / r₀)`
/// where `r₀ = hypot(z_s, y_s)`.
///
/// **Planar**: the flow is uniform in the z‑direction.  y follows the
/// flow, z stays constant:
///   `y_3d = r_f × (y_s / r₀)`, `z_3d = z_s`
///
/// If `r₀ ≈ 0` the point sits on the symmetry axis and the 3D
/// coordinates collapse to `(x_f, 0, 0)`.
fn transform_to_3d(
    mut streamline: WallPoints,
    origin: &Point3d,
    r0: f64,
    axisymmetric: bool,
) -> WallPoints {
    let zs = origin.x; // spatial (shape z)
    let ys = origin.y; // vertical  (shape y)

    if r0 < 1e-15 {
        for p in &mut streamline {
            p.y = 0.0;
            p.z = 0.0;
        }
        return streamline;
    }

    if axisymmetric {
        for p in &mut streamline {
            let rf = p.y;
            p.y = rf * zs / r0;
            p.z = rf * ys / r0;
        }
    } else {
        for p in &mut streamline {
            let rf = p.y;
            p.y = rf * ys / r0;
            p.z = zs;
        }
    }
    streamline
}
