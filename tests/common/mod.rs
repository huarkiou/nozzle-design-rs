//! 集成测试共享辅助函数。
//!
//! 提供 OTN 运行、SLTN 运行、结果验证等公共逻辑，
//! 避免各测试文件重复实现。
//!
//! 不同测试文件按需导入，未使用的函数允许 dead_code。

#![allow(dead_code)]

use aero::{
    moc::CharLines,
    nozzle::{ConstraintNozzle, Control, Geometry, Inlet, NozzleConfig, Outlet, Throat, IO},
    streamline_trace::{StreamlineConfig, StreamlineTrace},
    Material,
};
use geometry::ClosedCurve;

// ── OTN 运行 ──────────────────────────────────────────────────────────

/// 标准轴对称 OTN（自动选 θ_a，NASA9 变比热）。
pub fn run_otn_axisymmetric() -> CharLines {
    run_otn_axisymmetric_with_theta(f64::NAN)
}

/// 标准轴对称 OTN，指定初始膨胀角。
pub fn run_otn_axisymmetric_with_theta(theta_a: f64) -> CharLines {
    let config = NozzleConfig {
        control: Control {
            axisymmetric: true,
            ..Control::default()
        },
        material: Material::air_nasa9piecewise_polynomial(),
        inlet: Inlet::default(),
        geometry: Geometry::default(),
        throat: Throat {
            radius_throat: 0.0,
            theta_a,
        },
        outlet: Outlet::default(),
        io: IO::default(),
    };
    let mut nozzle = ConstraintNozzle::new_otn(config);
    nozzle.run();
    nozzle.get_assembly_charlines()
}

/// 标准平面 OTN（固定 θ_a=19°，NASA9 变比热）。
pub fn run_otn_planar() -> CharLines {
    let config = NozzleConfig {
        control: Control {
            axisymmetric: false,
            ..Control::default()
        },
        material: Material::air_nasa9piecewise_polynomial(),
        inlet: Inlet::default(),
        geometry: Geometry::default(),
        throat: Throat {
            radius_throat: 0.0,
            theta_a: 19.0_f64.to_radians(),
        },
        outlet: Outlet::default(),
        io: IO::default(),
    };
    let mut nozzle = ConstraintNozzle::new_otn(config);
    nozzle.run();
    nozzle.get_assembly_charlines()
}

/// 通用 OTN 运行（完全自定义配置，固定 θ_a 避免自动搜索）。
pub fn run_otn_custom(
    axisymmetric: bool,
    material: Material,
    height_i: f64,
    length: f64,
    ma: f64,
    p_total: f64,
    t_total: f64,
    theta_a: f64,
    p_ambient: f64,
) -> CharLines {
    let config = NozzleConfig {
        control: Control {
            axisymmetric,
            ..Control::default()
        },
        material,
        inlet: Inlet {
            p_total,
            temperature_total: t_total,
            ma,
            theta: 0.0,
        },
        geometry: Geometry {
            height_i,
            height_e: f64::NAN,
            length,
            width: 1.0,
        },
        throat: Throat {
            radius_throat: 0.0,
            theta_a,
        },
        outlet: Outlet { p_ambient },
        io: IO::default(),
    };
    let mut nozzle = ConstraintNozzle::new_otn(config);
    nozzle.run();
    nozzle.get_assembly_charlines()
}

// ── SLTN 运行 ─────────────────────────────────────────────────────────

/// 运行标准 SLTN 并返回 tracer。
pub fn run_sltn(
    charlines: &CharLines,
    inlet_shape: Box<dyn ClosedCurve>,
    outlet_shape: Box<dyn ClosedCurve>,
    axisymmetric: bool,
    n_theta: usize,
    n_axis: usize,
) -> StreamlineTrace {
    let config = StreamlineConfig {
        axisymmetric,
        datasource_inlet: charlines.clone(),
        datasource_outlet: charlines.clone(),
        inlet_shape,
        outlet_shape,
        monotonic: false,
        weight_parameter_a: 0.0,
    };
    let mut tracer = StreamlineTrace::new(config);
    tracer
        .run(n_theta, n_axis)
        .expect("SLTN trace should succeed");
    tracer
}

/// 运行带 monotonic / weight 选项的 SLTN。
pub fn run_sltn_with_options(
    charlines: &CharLines,
    inlet_shape: Box<dyn ClosedCurve>,
    outlet_shape: Box<dyn ClosedCurve>,
    axisymmetric: bool,
    n_theta: usize,
    n_axis: usize,
    monotonic: bool,
    weight_parameter_a: f64,
) -> StreamlineTrace {
    let config = StreamlineConfig {
        axisymmetric,
        datasource_inlet: charlines.clone(),
        datasource_outlet: charlines.clone(),
        inlet_shape,
        outlet_shape,
        monotonic,
        weight_parameter_a,
    };
    let mut tracer = StreamlineTrace::new(config);
    tracer
        .run(n_theta, n_axis)
        .expect("SLTN trace should succeed");
    tracer
}

// ── 验证 ──────────────────────────────────────────────────────────────

/// 验证 OTN 产出的特征线基本有效性。
pub fn assert_charlines_valid(charlines: &CharLines) {
    assert!(!charlines.is_empty(), "charlines should not be empty");
    for (li, line) in charlines.iter().enumerate() {
        for (pi, pt) in line.iter().enumerate() {
            assert!(pt.is_valid(), "invalid MOC point: line={li} pt={pi}");
            assert!(
                pt.x >= -1e-9 && pt.y >= -1e-9,
                "coordinate out of bounds: line={li} pt={pi}: x={}, y={}",
                pt.x,
                pt.y
            );
        }
    }
}

/// 验证 SLTN 模型的完整性。
pub fn assert_sltn_valid(tracer: &StreamlineTrace, expected_len: f64) {
    assert!(!tracer.model.is_empty(), "SLTN model should not be empty");
    assert!(!tracer.downstream.is_empty(), "downstream profiles empty");
    assert!(!tracer.upstream.is_empty(), "upstream profiles empty");

    for (i, profile) in tracer.model.iter().enumerate() {
        assert!(profile.len() >= 2, "profile {i}: too few points");

        for w in profile.windows(2) {
            assert!(
                w[0].x < w[1].x,
                "profile {i}: x not monotonic: {} >= {}",
                w[0].x,
                w[1].x
            );
        }

        let first_x = profile.first().unwrap().x;
        assert!(
            first_x.abs() < 0.5,
            "profile {i}: first x = {first_x}, expected ≈ 0"
        );

        let last_x = profile.last().unwrap().x;
        assert!(
            (last_x - expected_len).abs() < expected_len * 0.15,
            "profile {i}: last x = {last_x}, expected ≈ {expected_len}"
        );
    }
}
