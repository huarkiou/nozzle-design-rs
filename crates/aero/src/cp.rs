use std::sync::Arc;

use serde::{Deserialize, Deserializer, Serialize};

// 热容类型枚举
#[derive(Clone)]
pub enum Cp {
    Constant(f64),
    Variable(Arc<dyn Fn(f64) -> f64>),
}

impl PartialEq for Cp {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::Constant(l0), Self::Constant(r0)) => l0 == r0,
            (Self::Variable(l0), Self::Variable(r0)) => Arc::ptr_eq(l0, r0),
            _ => false,
        }
    }
}

impl Cp {
    pub fn new(cp: impl Fn(f64) -> f64 + 'static) -> Self {
        Self::Variable(Arc::new(cp))
    }

    pub fn eval(&self, temperature: f64) -> f64 {
        match self {
            Cp::Constant(v) => *v,
            Cp::Variable(f) => f(temperature),
        }
    }
}

impl From<f64> for Cp {
    fn from(value: f64) -> Self {
        Cp::Constant(value)
    }
}

impl std::fmt::Debug for Cp {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Constant(cp_value) => f.debug_tuple("Constant").field(cp_value).finish(),
            Self::Variable(cp_func) => f
                .debug_tuple("Variable")
                .field(&"unknown")
                .field(&cp_func(273.15))
                .finish(),
        }
    }
}

impl Serialize for Cp {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        match self {
            Self::Constant(cp_value) => serializer.serialize_f64(*cp_value),
            Self::Variable(cp_func) => serializer.serialize_f64(cp_func(273.15)),
        }
    }
}

impl<'de> Deserialize<'de> for Cp {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let value = f64::deserialize(deserializer)?;
        Ok(Cp::Constant(value))
    }
}
