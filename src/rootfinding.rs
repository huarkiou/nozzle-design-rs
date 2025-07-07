pub mod error;
pub mod rootbracket;
pub mod toms748;

pub use error::RootFindingError;
pub use rootbracket::RootBracket;
pub use toms748::toms748_bracket;
