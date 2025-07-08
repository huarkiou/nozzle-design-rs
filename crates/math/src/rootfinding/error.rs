use std::fmt::Display;

#[derive(Debug)]
pub enum RootFindingError {
    FunctionNotBracketed,
    MaxIterationsExceeded,
}

impl Display for RootFindingError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RootFindingError::FunctionNotBracketed => write!(f, "Function Not Bracketed"),
            RootFindingError::MaxIterationsExceeded => write!(f, "Max Iterations Exceeded"),
        }
    }
}
