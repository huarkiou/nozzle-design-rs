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

// ── 分段多项式定义 ──────────────────────────────────────────

/// 分段多项式的一个温度区间
///
/// 多项式形式: `Cp(T) = c[0] + c[1]·T + c[2]·T² + … + c[n]·Tⁿ`
///
/// `coefficients` 数组长度 = 多项式次数 + 1，支持任意次数。
/// 求值使用 Horner 方法以保证数值稳定性。
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct CpSegment {
    /// 温度下界 (K)
    pub t_min: f64,
    /// 温度上界 (K)
    pub t_max: f64,
    /// 多项式系数 [c0, c1, c2, …, cn]
    #[serde(default)]
    pub coefficients: Vec<f64>,
}

impl CpSegment {
    /// Horner 法求值: c0 + t*(c1 + t*(c2 + … + t*cn))
    pub fn eval(&self, t: f64) -> f64 {
        self.coefficients.iter().rfold(0.0, |acc, &c| acc * t + c)
    }
}

// ── Cp 构造方法 ─────────────────────────────────────────────

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

    /// 从分段多项式构建变比热容
    ///
    /// - 温度区间内按多项式求值
    /// - 低于最低区间时，取最低区间下界的值（外推恒定）
    /// - 高于最高区间时，取最高区间上界的值（外推恒定）
    ///
    /// # Panics
    /// 如果 `segments` 为空。
    pub fn from_piecewise_segments(segments: Vec<CpSegment>) -> Self {
        assert!(!segments.is_empty(), "cp_segments must not be empty");
        let segs = segments;
        Self::new(move |t| {
            for seg in &segs {
                if t >= seg.t_min && t <= seg.t_max {
                    return seg.eval(t);
                }
            }
            // 超出所有区间 — 取最近边界的值
            if t < segs[0].t_min {
                segs[0].eval(segs[0].t_min)
            } else {
                let last = segs.last().unwrap();
                last.eval(last.t_max)
            }
        })
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

    // ── 分段多项式测试 ────────────────────────────────────

    #[test]
    fn test_segment_eval_linear() {
        // Cp(T) = 1000 + 0.2·T
        let seg = CpSegment {
            t_min: 200.0,
            t_max: 1000.0,
            coefficients: vec![1000.0, 0.2],
        };
        assert!((seg.eval(300.0) - 1060.0).abs() < 1e-10);
        assert!((seg.eval(500.0) - 1100.0).abs() < 1e-10);
    }

    #[test]
    fn test_segment_eval_quadratic() {
        // Cp(T) = 1 + 2T + 3T²
        let seg = CpSegment {
            t_min: 100.0,
            t_max: 500.0,
            coefficients: vec![1.0, 2.0, 3.0],
        };
        let t = 10.0;
        let expected = 1.0 + 2.0 * t + 3.0 * t * t;
        assert!((seg.eval(t) - expected).abs() < 1e-10);
    }

    #[test]
    fn test_segment_eval_constant_polynomial() {
        // 常数多项式: Cp(T) = 500
        let seg = CpSegment {
            t_min: 0.0,
            t_max: 1000.0,
            coefficients: vec![500.0],
        };
        assert!((seg.eval(300.0) - 500.0).abs() < 1e-10);
    }

    #[test]
    fn test_piecewise_in_range() {
        let cp = Cp::from_piecewise_segments(vec![
            CpSegment {
                t_min: 200.0,
                t_max: 1000.0,
                coefficients: vec![1000.0, 0.2], // 1000 + 0.2·T
            },
            CpSegment {
                t_min: 1000.0,
                t_max: 3000.0,
                coefficients: vec![800.0, 0.4], // 800 + 0.4·T
            },
        ]);
        // 第一段内
        let v = cp.eval(500.0);
        assert!((v - 1100.0).abs() < 1e-10, "got {v}");
        // 第二段内
        let v = cp.eval(2000.0);
        assert!((v - 1600.0).abs() < 1e-10, "got {v}");
        // 段边界
        let v = cp.eval(1000.0);
        assert!((v - 1200.0).abs() < 1e-10, "got {v}");
    }

    #[test]
    fn test_piecewise_below_range() {
        let cp = Cp::from_piecewise_segments(vec![CpSegment {
            t_min: 200.0,
            t_max: 1000.0,
            coefficients: vec![0.0, 1.0], // T
        }]);
        // 低于最低区间 → 取 t=200 的值
        assert!((cp.eval(100.0) - 200.0).abs() < 1e-10);
        assert!((cp.eval(0.0) - 200.0).abs() < 1e-10);
    }

    #[test]
    fn test_piecewise_above_range() {
        let cp = Cp::from_piecewise_segments(vec![CpSegment {
            t_min: 200.0,
            t_max: 1000.0,
            coefficients: vec![0.0, 1.0], // T
        }]);
        // 高于最高区间 → 取 t=1000 的值
        assert!((cp.eval(1500.0) - 1000.0).abs() < 1e-10);
    }

    #[test]
    #[should_panic(expected = "cp_segments must not be empty")]
    fn test_piecewise_empty_panics() {
        Cp::from_piecewise_segments(vec![]);
    }

    #[test]
    fn test_segment_eval_many_coefficients() {
        // Cp(T) = 1 + 2T + 3T² + 4T³ + 5T⁴ (5次系数)
        let seg = CpSegment {
            t_min: 0.0,
            t_max: 1000.0,
            coefficients: vec![1.0, 2.0, 3.0, 4.0, 5.0],
        };
        let t: f64 = 2.0;
        let expected = 1.0 + 2.0 * t + 3.0 * t.powi(2) + 4.0 * t.powi(3) + 5.0 * t.powi(4);
        assert!((seg.eval(t) - expected).abs() < 1e-10);
    }
}
