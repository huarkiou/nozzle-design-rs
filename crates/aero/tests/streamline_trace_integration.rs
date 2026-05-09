//! Integration test for the full streamline-trace pipeline using a
//! minimal synthetic flow field.

use aero::moc::{CharLine, CharLines, MocPoint};
use aero::streamline_trace::{CharLineSource, StreamlineConfig, StreamlineTrace};
use geometry::{Circle, Ellipse};

/// Build a minimal 3-line axisymmetric flow field.
///
/// Inlet (x=0): 3 points from axis (y=0) to wall (y=1).
/// Mid   (x=3): 3 points.
/// Outlet(x=6): 3 points, wall at y=0.7 (nozzle expands).
fn make_test_field() -> CharLines {
    let mat = aero::Material::from_rgas_gamma(287.0, 1.4);

    let make_pt = |x: f64, y: f64, u: f64| -> MocPoint {
        MocPoint::new(x, y, u, 0.0, 100_000.0, 300.0, 1.2, mat.clone())
    };

    let mut lines = CharLines::new();

    // Inlet line (x=0, vertical)
    let mut inlet = CharLine::new();
    inlet.push(make_pt(0.0, 0.0, 500.0)); // axis
    inlet.push(make_pt(0.0, 0.5, 500.0));
    inlet.push(make_pt(0.0, 1.0, 500.0)); // wall
    lines.push(inlet);

    // Mid line
    let mut mid = CharLine::new();
    mid.push(make_pt(3.0, 0.0, 520.0));
    mid.push(make_pt(3.0, 0.4, 520.0));
    mid.push(make_pt(3.0, 0.85, 520.0));
    lines.push(mid);

    // Outlet line (x=6, vertical)
    let mut outlet = CharLine::new();
    outlet.push(make_pt(6.0, 0.0, 540.0));
    outlet.push(make_pt(6.0, 0.35, 540.0));
    outlet.push(make_pt(6.0, 0.70, 540.0));
    lines.push(outlet);

    lines
}

#[test]
#[ignore = "integration test, run with -- --ignored"]
fn test_full_pipeline_axisymmetric() {
    let field = make_test_field();

    let config = StreamlineConfig {
        axisymmetric: true,
        datasource_inlet: field.clone(),
        datasource_outlet: field,
        inlet_shape: Box::new(Circle::new(1.0, (0.0, 0.0)).unwrap()),
        outlet_shape: Box::new(Ellipse::new(0.7, 0.5, 0.0, (0.0, 0.0)).unwrap()),
        monotonic: false,
        weight_parameter_a: 0.0,
    };

    let mut tracer = StreamlineTrace::new(config);
    tracer.run(12, 20).expect("trace should succeed");

    // model: one profile per theta (12 profiles, each with ~20 points)
    assert_eq!(tracer.model.len(), 12, "expected 12 azimuthal profiles");

    for (i, profile) in tracer.model.iter().enumerate() {
        assert!(
            profile.len() >= 2,
            "profile {i}: expected ≥2 points, got {}",
            profile.len()
        );
        // x should be monotonically increasing
        for w in profile.windows(2) {
            assert!(
                w[0].x < w[1].x,
                "profile {i}: x not monotonic: {} ≥ {}",
                w[0].x,
                w[1].x
            );
        }
        // inlet x ≈ 0, outlet x ≈ 6
        assert!(
            (profile.first().unwrap().x - 0.0).abs() < 0.1,
            "profile {i}: first x = {}",
            profile.first().unwrap().x
        );
        assert!(
            (profile.last().unwrap().x - 6.0).abs() < 0.1,
            "profile {i}: last x = {}",
            profile.last().unwrap().x
        );
        // downstream profiles should exist
    }

    assert!(!tracer.downstream.is_empty());
    assert!(!tracer.upstream.is_empty());
}

#[test]
#[ignore = "integration test, run with -- --ignored"]
fn test_full_pipeline_planar() {
    let field = make_test_field();

    let config = StreamlineConfig {
        axisymmetric: false,
        datasource_inlet: field.clone(),
        datasource_outlet: field,
        inlet_shape: Box::new(Circle::new(1.0, (0.0, 0.0)).unwrap()),
        outlet_shape: Box::new(Circle::new(0.7, (0.0, 0.0)).unwrap()),
        monotonic: false,
        weight_parameter_a: 0.0,
    };

    let mut tracer = StreamlineTrace::new(config);
    tracer.run(8, 15).expect("planar trace should succeed");

    assert_eq!(tracer.model.len(), 8);
    for profile in &tracer.model {
        assert!(profile.len() >= 2);
        // Planar: z should be constant per streamline
        // (each streamline has the same z as its inlet point)
    }
}

#[test]
#[ignore = "integration test, run with -- --ignored"]
fn test_rev_char_lines() {
    let field = make_test_field();
    let rev = aero::streamline_trace::RevCharLines::new(&field);

    assert_eq!(rev.len(), field.len());
    // First element of reversed = last of original
    let fwd_first = field.get_line(0).unwrap();
    let rev_first = rev.get_line(0).unwrap();
    assert_eq!(fwd_first[0].x, 0.0); // inlet x
    assert_eq!(rev_first[0].x, 6.0); // outlet x (reversed)
}

#[test]
#[ignore = "integration test, run with -- --ignored"]
fn test_error_on_empty_datasource() {
    let config = StreamlineConfig {
        axisymmetric: true,
        datasource_inlet: CharLines::new(),
        datasource_outlet: CharLines::new(),
        inlet_shape: Box::new(Circle::new(1.0, (0.0, 0.0)).unwrap()),
        outlet_shape: Box::new(Circle::new(0.5, (0.0, 0.0)).unwrap()),
        monotonic: false,
        weight_parameter_a: 0.0,
    };
    let mut tracer = StreamlineTrace::new(config);
    let err = tracer.run(4, 5).unwrap_err();
    let msg = format!("{err}");
    assert!(msg.contains("empty"), "expected 'empty' in error: {msg}");
}
