use std::sync::Arc;

use serde::{Deserialize, Deserializer, Serialize};

// 热容类型枚举
#[derive(Clone)]
pub enum Cp {
    Constant(f64),
    Variable(Arc<dyn Fn(f64) -> f64 + Send + Sync>),
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
    pub fn new(cp: impl Fn(f64) -> f64 + Send + Sync + 'static) -> Self {
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
            Self::Variable(_) => serializer.serialize_f64(f64::NAN),
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_constant_eval() {
        let cp = Cp::Constant(1005.0);
        assert_eq!(cp.eval(300.0), 1005.0);
        assert_eq!(cp.eval(0.0), 1005.0);
        assert_eq!(cp.eval(5000.0), 1005.0);
    }

    #[test]
    fn test_variable_eval() {
        let cp = Cp::new(|t: f64| 1000.0 + 0.2 * t);
        assert!((cp.eval(300.0) - 1060.0).abs() < 1e-10);
        assert!((cp.eval(1000.0) - 1200.0).abs() < 1e-10);
    }

    #[test]
    fn test_new_creates_variable() {
        let cp = Cp::new(|t: f64| t * 2.0);
        match cp {
            Cp::Variable(_) => {} // ok
            Cp::Constant(_) => panic!("expected Variable"),
        }
    }

    #[test]
    fn test_from_f64() {
        let cp: Cp = 1005.0_f64.into();
        match cp {
            Cp::Constant(v) => assert_eq!(v, 1005.0),
            Cp::Variable(_) => panic!("expected Constant"),
        }
    }

    #[test]
    fn test_partial_eq_constant() {
        assert_eq!(Cp::Constant(1005.0), Cp::Constant(1005.0));
        assert_ne!(Cp::Constant(1005.0), Cp::Constant(1006.0));
    }

    #[test]
    fn test_partial_eq_variable_ptr() {
        let f = |t: f64| t;
        let a = Cp::new(f);
        // Same closure → different Arc pointers → not equal
        let b = Cp::new(|t: f64| t);
        assert_ne!(a, b);
    }

    #[test]
    fn test_partial_eq_mixed() {
        assert_ne!(Cp::Constant(1005.0), Cp::new(|t: f64| t));
    }

    #[test]
    fn test_debug_constant() {
        let cp = Cp::Constant(1005.0);
        let s = format!("{:?}", cp);
        assert!(s.contains("Constant"));
        assert!(s.contains("1005"));
    }

    #[test]
    fn test_debug_variable() {
        let cp = Cp::new(|t: f64| 1000.0 + 0.2 * t);
        let s = format!("{:?}", cp);
        assert!(s.contains("Variable"));
    }

    #[test]
    fn test_serialize_constant() {
        let cp = Cp::Constant(1005.0);
        let json = serde_json::to_string(&cp).unwrap();
        assert_eq!(json, "1005.0");
    }

    #[test]
    fn test_serialize_variable() {
        // Variable serializes as nan in TOML (round-trip triggers nasa9 fallback)
        let cp = Cp::new(|t: f64| 1000.0 + 0.2 * t);
        #[derive(Serialize)]
        struct W {
            cp: Cp,
        }
        let toml_str = toml::to_string(&W { cp }).unwrap();
        assert!(toml_str.contains("nan"), "expected nan, got: {toml_str}");
    }

    #[test]
    fn test_deserialize_always_constant() {
        let cp: Cp = serde_json::from_str("1005.0").unwrap();
        match cp {
            Cp::Constant(v) => assert_eq!(v, 1005.0),
            Cp::Variable(_) => panic!("deserialization should produce Constant"),
        }
    }

    #[test]
    fn test_clone_constant() {
        let a = Cp::Constant(1005.0);
        let b = a.clone();
        assert_eq!(a, b);
    }

    #[test]
    fn test_clone_variable() {
        let a = Cp::new(|t: f64| t);
        let b = a.clone();
        // Same Arc pointer → equal
        assert_eq!(a, b);
    }
}
