pub mod bisection;
mod error;
pub mod newton2d;
mod rootbracket;
pub mod secant;
pub mod toms748;

pub use error::RootFindingError;
pub use newton2d::newton_2d;
pub use rootbracket::RootBracket;
