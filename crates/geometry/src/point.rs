use std::ops::{Add, Div, Mul, Sub};

/// A point in 3D space.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point3d {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Point3d {
    /// Create a new `Point3d` from coordinates.
    #[inline]
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// Euclidean distance to another point.
    #[inline]
    pub fn distance_to(&self, other: &Point3d) -> f64 {
        self.distance_squared_to(other).sqrt()
    }

    /// Squared Euclidean distance to another point.
    #[inline]
    pub fn distance_squared_to(&self, other: &Point3d) -> f64 {
        (self.x - other.x).powi(2) + (self.y - other.y).powi(2) + (self.z - other.z).powi(2)
    }

    /// Polar angle in the XY plane (atan2(y, x)), normalized to [0, 2π).
    #[inline]
    pub fn polar_angle(&self) -> f64 {
        let angle = self.y.atan2(self.x);
        if angle < 0.0 {
            angle + std::f64::consts::TAU
        } else {
            angle
        }
    }

    /// Rotate this point counter-clockwise by `alpha` radians around a
    /// z-axis line passing through `origin`.
    #[inline]
    pub fn rotate_around_z(&mut self, origin: &Point3d, alpha: f64) {
        let dx = self.x - origin.x;
        let dy = self.y - origin.y;
        let cos_a = alpha.cos();
        let sin_a = alpha.sin();
        self.x = origin.x + dx * cos_a - dy * sin_a;
        self.y = origin.y + dx * sin_a + dy * cos_a;
    }

    /// Linearly interpolate between `self` and `other` at the given x
    /// coordinate.  Assumes the two points have distinct x values.
    #[inline]
    pub fn interpolate_x(&self, other: &Point3d, x: f64) -> Point3d {
        let dx = other.x - self.x;
        if dx.abs() < f64::EPSILON {
            // Points share the same x; return the midpoint projection.
            return Point3d {
                x,
                y: (self.y + other.y) * 0.5,
                z: (self.z + other.z) * 0.5,
            };
        }
        let t = (x - self.x) / dx;
        Point3d {
            x,
            y: self.y + t * (other.y - self.y),
            z: self.z + t * (other.z - self.z),
        }
    }

    /// Linearly interpolate between `self` and `other` at the given y
    /// coordinate.  Assumes the two points have distinct y values.
    #[inline]
    pub fn interpolate_y(&self, other: &Point3d, y: f64) -> Point3d {
        let dy = other.y - self.y;
        if dy.abs() < f64::EPSILON {
            return Point3d {
                x: (self.x + other.x) * 0.5,
                y,
                z: (self.z + other.z) * 0.5,
            };
        }
        let t = (y - self.y) / dy;
        Point3d {
            x: self.x + t * (other.x - self.x),
            y,
            z: self.z + t * (other.z - self.z),
        }
    }

    /// Vector magnitude (length from origin).
    #[inline]
    pub fn length(&self) -> f64 {
        (self.x.powi(2) + self.y.powi(2) + self.z.powi(2)).sqrt()
    }
}

impl Default for Point3d {
    fn default() -> Self {
        Self {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        }
    }
}

// ---------------------------------------------------------------
// Arithmetic operators
// ---------------------------------------------------------------

impl Add for Point3d {
    type Output = Point3d;

    #[inline]
    fn add(self, rhs: Self) -> Self::Output {
        Self::Output {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
        }
    }
}

impl Sub for Point3d {
    type Output = Point3d;

    #[inline]
    fn sub(self, rhs: Self) -> Self::Output {
        Self::Output {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

impl Mul<f64> for Point3d {
    type Output = Point3d;

    #[inline]
    fn mul(self, rhs: f64) -> Self::Output {
        Self::Output {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}

impl Div<f64> for Point3d {
    type Output = Point3d;

    #[inline]
    fn div(self, rhs: f64) -> Self::Output {
        Self::Output {
            x: self.x / rhs,
            y: self.y / rhs,
            z: self.z / rhs,
        }
    }
}

// ---------------------------------------------------------------
// Serde support (behind feature flag)
// ---------------------------------------------------------------

#[cfg(feature = "serde")]
impl serde::Serialize for Point3d {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        use serde::ser::SerializeTuple;
        let mut tup = serializer.serialize_tuple(3)?;
        tup.serialize_element(&self.x)?;
        tup.serialize_element(&self.y)?;
        tup.serialize_element(&self.z)?;
        tup.end()
    }
}

#[cfg(feature = "serde")]
impl<'de> serde::Deserialize<'de> for Point3d {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        deserializer.deserialize_tuple(3, Point3dVisitor)
    }
}

#[cfg(feature = "serde")]
struct Point3dVisitor;

#[cfg(feature = "serde")]
impl<'de> serde::de::Visitor<'de> for Point3dVisitor {
    type Value = Point3d;

    fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
        formatter.write_str("a tuple of three f64 values")
    }

    fn visit_seq<A: serde::de::SeqAccess<'de>>(self, mut seq: A) -> Result<Point3d, A::Error> {
        let x = seq
            .next_element()?
            .ok_or_else(|| serde::de::Error::invalid_length(0, &self))?;
        let y = seq
            .next_element()?
            .ok_or_else(|| serde::de::Error::invalid_length(1, &self))?;
        let z = seq
            .next_element()?
            .ok_or_else(|| serde::de::Error::invalid_length(2, &self))?;
        Ok(Point3d { x, y, z })
    }
}

/// A collection of wall points (ordered set of 3D points).
pub type WallPoints = Vec<Point3d>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new_default() {
        let p = Point3d::new(1.0, 2.0, 3.0);
        assert_eq!(p.x, 1.0);
        assert_eq!(p.y, 2.0);
        assert_eq!(p.z, 3.0);

        let q = Point3d::default();
        assert_eq!(q, Point3d::new(0.0, 0.0, 0.0));
    }

    #[test]
    fn test_distance() {
        let p = Point3d::new(0.0, 0.0, 0.0);
        let q = Point3d::new(3.0, 4.0, 0.0);
        assert!((p.distance_to(&q) - 5.0).abs() < 1e-12);
        assert!((p.distance_squared_to(&q) - 25.0).abs() < 1e-12);
    }

    #[test]
    fn test_polar_angle() {
        let p = Point3d::new(1.0, 0.0, 0.0);
        assert!((p.polar_angle() - 0.0).abs() < 1e-12);

        let q = Point3d::new(0.0, 1.0, 0.0);
        assert!((q.polar_angle() - std::f64::consts::FRAC_PI_2).abs() < 1e-12);

        let r = Point3d::new(-1.0, 0.0, 0.0);
        assert!((r.polar_angle() - std::f64::consts::PI).abs() < 1e-12);

        let s = Point3d::new(0.0, -1.0, 0.0);
        assert!((s.polar_angle() - 3.0 * std::f64::consts::FRAC_PI_2).abs() < 1e-12);
    }

    #[test]
    fn test_rotate_around_z() {
        let mut p = Point3d::new(1.0, 0.0, 5.0);
        let origin = Point3d::new(0.0, 0.0, 5.0);
        p.rotate_around_z(&origin, std::f64::consts::FRAC_PI_2);
        assert!((p.x - 0.0).abs() < 1e-12);
        assert!((p.y - 1.0).abs() < 1e-12);
        assert!((p.z - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_interpolate_x() {
        let a = Point3d::new(0.0, 0.0, 0.0);
        let b = Point3d::new(2.0, 4.0, 6.0);
        let mid = a.interpolate_x(&b, 0.5);
        assert!((mid.x - 0.5).abs() < 1e-12);
        assert!((mid.y - 1.0).abs() < 1e-12);
        assert!((mid.z - 1.5).abs() < 1e-12);
    }

    #[test]
    fn test_arithmetic() {
        let a = Point3d::new(1.0, 2.0, 3.0);
        let b = Point3d::new(4.0, 5.0, 6.0);
        let sum = a + b;
        assert_eq!(sum, Point3d::new(5.0, 7.0, 9.0));
        let diff = a - b;
        assert_eq!(diff, Point3d::new(-3.0, -3.0, -3.0));
        let scaled = a * 2.0;
        assert_eq!(scaled, Point3d::new(2.0, 4.0, 6.0));
        let div = a / 2.0;
        assert_eq!(div, Point3d::new(0.5, 1.0, 1.5));
    }

    #[test]
    fn test_length() {
        let p = Point3d::new(3.0, 4.0, 0.0);
        assert!((p.length() - 5.0).abs() < 1e-12);
    }
}
