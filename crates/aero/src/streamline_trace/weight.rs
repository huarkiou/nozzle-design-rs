/// Build the weight function `f(x)` where `f(0) = 0`, `f(1) = 1`.
///
/// The single parameter `a` controls the shape:
///
/// * `a = 0` → linear: `f(x) = x`
/// * `a > 0` → S‑shaped (arctan family):
///   `f(x) = atan((2x − 1)·a) / atan(a) / 2 + 0.5`
/// * `a < 0` → inverse S‑shaped (tan family):
///   `f(x) = (tan(atan(a)·(2x − 1)) / a + 1) / 2`
///
/// The returned closure is `Send + Sync` so that it can be stored
/// inside `StreamlineTrace`.
pub fn build_weight_function(a: f64) -> Box<dyn Fn(f64) -> f64 + Send + Sync> {
    if a == 0.0 {
        Box::new(move |x: f64| x)
    } else if a > 0.0 {
        let a_atan = a.atan();
        Box::new(move |x: f64| {
            let t = 2.0 * x - 1.0;
            ((t * a).atan() / a_atan) / 2.0 + 0.5
        })
    } else {
        // a < 0
        let a_atan = a.atan();
        Box::new(move |x: f64| {
            let t = 2.0 * x - 1.0;
            ((t * a_atan).tan() / a + 1.0) / 2.0
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx_eq(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-10
    }

    #[test]
    fn test_linear_weight() {
        let f = build_weight_function(0.0);
        assert!(approx_eq(f(0.0), 0.0));
        assert!(approx_eq(f(0.5), 0.5));
        assert!(approx_eq(f(1.0), 1.0));
    }

    #[test]
    fn test_positive_a() {
        let f = build_weight_function(1.0);
        // f(0) should be close to 0, f(1) close to 1, symmetric
        assert!(f(0.0) >= -0.01 && f(0.0) <= 0.01);
        assert!(f(1.0) >= 0.99 && f(1.0) <= 1.01);
        // monotonic increasing
        assert!(f(0.3) < f(0.7));
    }

    #[test]
    fn test_negative_a() {
        let f = build_weight_function(-1.0);
        assert!(f(0.0) >= -0.01 && f(0.0) <= 0.01);
        assert!(f(1.0) >= 0.99 && f(1.0) <= 1.01);
        assert!(f(0.3) < f(0.7));
    }
}
