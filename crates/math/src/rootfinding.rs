pub mod bisection;
mod error;
mod rootbracket;
pub mod secant;
pub mod toms748;

pub use error::RootFindingError;
pub use rootbracket::RootBracket;
