pub mod basics;
pub mod circle;
pub mod closed_curve;
pub mod ellipse;
pub mod obj;
pub mod point;
pub mod rectangular;
pub mod superellipse;
pub mod userdefined;
pub mod wallpoints;

// Convenience re-exports of the main public types.
pub use basics::{ray_segment_intersection, sgn};
pub use circle::Circle;
pub use closed_curve::ClosedCurve;
pub use ellipse::Ellipse;
pub use obj::ObjModel;
pub use point::Point3d;
pub use point::WallPoints;
pub use rectangular::Rectangular;
pub use superellipse::SuperEllipse;
pub use userdefined::UserDefined;
