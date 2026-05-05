pub mod bisection;
mod error;
pub mod newton2d;
pub mod secant;
pub mod toms748;

pub use bisection::{max_iterations, solve_bracket, BisectBracket};
pub use error::RootFindingError;
pub use newton2d::newton_2d;
