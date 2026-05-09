use geometry::{Circle, ClosedCurve, Ellipse, Point3d, Rectangular, SuperEllipse, UserDefined};

use crate::config::ShapeConfig;

/// Build a cross‑section shape from TOML configuration.
pub fn build_shape(
    shape_config: &ShapeConfig,
    normalize_factor: f64,
) -> Result<Box<dyn ClosedCurve>, String> {
    let factor = if shape_config.normalized {
        1.0
    } else {
        normalize_factor
    };
    let center: (f64, f64) = shape_config
        .center
        .as_ref()
        .map(|c| {
            if c.len() == 2 {
                Ok((c[0] / factor, c[1] / factor))
            } else {
                Err(format!(
                    "center must have exactly 2 elements, got {}",
                    c.len()
                ))
            }
        })
        .transpose()
        .map(|r| r.unwrap_or((0.0, 0.0)))
        .unwrap_or((0.0, 0.0));

    match shape_config.shape.as_str() {
        "circle" => {
            let r = shape_config.radius.ok_or("circle shape requires 'radius'")? / factor;
            Ok(Box::new(Circle::new(r, center).map_err(|e| format!("Circle: {e}"))?))
        }
        "ellipse" => {
            let a = shape_config.a.ok_or("ellipse shape requires 'a'")? / factor;
            let b = shape_config.b.ok_or("ellipse shape requires 'b'")? / factor;
            let alpha = shape_config.alpha.map(|deg| deg.to_radians()).unwrap_or(0.0);
            Ok(Box::new(Ellipse::new(a, b, alpha, center).map_err(|e| format!("Ellipse: {e}"))?))
        }
        "rectangular" => {
            let length = shape_config.length.ok_or("rectangular shape requires 'length'")? / factor;
            let width = shape_config.width.ok_or("rectangular shape requires 'width'")? / factor;
            let alpha = shape_config.alpha.map(|deg| deg.to_radians()).unwrap_or(0.0);
            Ok(Box::new(Rectangular::new(length, width, alpha, center).map_err(|e| format!("Rectangular: {e}"))?))
        }
        "superellipse" => {
            let a = shape_config.a.ok_or("superellipse shape requires 'a'")? / factor;
            let b = shape_config.b.ok_or("superellipse shape requires 'b'")? / factor;
            let power = shape_config.n.unwrap_or(2.0);
            let alpha = shape_config.alpha.map(|deg| deg.to_radians()).unwrap_or(0.0);
            Ok(Box::new(SuperEllipse::new(a, b, power, alpha, center).map_err(|e| format!("SuperEllipse: {e}"))?))
        }
        "userdefined" => {
            let ds = shape_config.datasource.as_ref().ok_or("userdefined shape requires 'datasource'")?;
            let points = geometry::wallpoints::read_points_2d(ds.to_str().unwrap_or(""))
                .map_err(|e| format!("failed to read userdefined datasource '{}': {e}", ds.display()))?;
            let pts: Vec<Point3d> = points
                .into_iter()
                .map(|(z, y)| Point3d::new(z / factor, y / factor, 0.0))
                .collect();
            let alpha = shape_config.alpha.map(|deg| deg.to_radians()).unwrap_or(0.0);
            Ok(Box::new(UserDefined::new(pts, alpha, Some(center)).map_err(|e| format!("UserDefined: {e}"))?))
        }
        other => Err(format!(
            "unknown shape '{}'. Supported: circle, ellipse, rectangular, superellipse, userdefined",
            other
        )),
    }
}
