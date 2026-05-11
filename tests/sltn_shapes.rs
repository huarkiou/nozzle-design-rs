/// SLTN 截面形状参数覆盖集成测试。
///
/// 覆盖所有支持的截面形状（圆形、椭圆、矩形、超椭圆）及其组合，
/// 验证不同形状在轴对称/平面基准流场下的三维流线追踪正确性。
///
/// 所有测试标记为 `#[ignore]`，使用 `cargo test -- --ignored` 运行。
mod common;

use common::{assert_sltn_valid, run_otn_axisymmetric, run_otn_planar, run_sltn};
use geometry::{Circle, Ellipse, Rectangular, SuperEllipse};

// ═══════════════════════════════════════════════════════════════════════
// 单形状匹配（进口=出口同类型）
// ═══════════════════════════════════════════════════════════════════════

/// 圆形→圆形（等比例缩小），轴对称。
#[test]
#[ignore = "SLTN shape test, run with -- --ignored"]
fn test_sltn_circle_to_circle_axisymmetric() {
    let cl = run_otn_axisymmetric();
    let tracer = run_sltn(
        &cl,
        Box::new(Circle::new(1.0, (0.0, 0.0)).unwrap()),
        Box::new(Circle::new(0.6, (0.0, 0.0)).unwrap()),
        true,
        12,
        20,
    );
    assert_sltn_valid(&tracer, 6.0);
}

/// 椭圆→椭圆，轴对称。
#[test]
#[ignore = "SLTN shape test, run with -- --ignored"]
fn test_sltn_ellipse_to_ellipse_axisymmetric() {
    let cl = run_otn_axisymmetric();
    let tracer = run_sltn(
        &cl,
        Box::new(Ellipse::new(0.9, 0.7, 0.0, (0.0, 0.0)).unwrap()),
        Box::new(Ellipse::new(0.6, 0.5, 0.0, (0.0, 0.0)).unwrap()),
        true,
        12,
        20,
    );
    assert_sltn_valid(&tracer, 6.0);
}

/// 矩形→矩形，轴对称。
#[test]
#[ignore = "SLTN shape test, run with -- --ignored"]
fn test_sltn_rectangular_same_axisymmetric() {
    let cl = run_otn_axisymmetric();
    let tracer = run_sltn(
        &cl,
        Box::new(Rectangular::new(1.6, 0.8, 0.0, (0.0, 0.0)).unwrap()),
        Box::new(Rectangular::new(1.2, 1.0, 0.0, (0.0, 0.0)).unwrap()),
        true,
        16,
        30,
    );
    assert_sltn_valid(&tracer, 6.0);
}

/// 超椭圆→超椭圆 (n=3)，轴对称。
#[test]
#[ignore = "SLTN shape test, run with -- --ignored"]
fn test_sltn_superellipse_same_axisymmetric() {
    let cl = run_otn_axisymmetric();
    let tracer = run_sltn(
        &cl,
        Box::new(SuperEllipse::new(0.85, 0.85, 3.0, 0.0, (0.0, 0.0)).unwrap()),
        Box::new(SuperEllipse::new(0.6, 0.7, 3.0, 0.0, (0.0, 0.0)).unwrap()),
        true,
        16,
        25,
    );
    assert_sltn_valid(&tracer, 6.0);
}

// ═══════════════════════════════════════════════════════════════════════
// 混合形状（进口≠出口不同类型）
// ═══════════════════════════════════════════════════════════════════════

/// 圆→矩形，轴对称。
#[test]
#[ignore = "SLTN shape test, run with -- --ignored"]
fn test_sltn_circle_to_rectangular() {
    let cl = run_otn_axisymmetric();
    let tracer = run_sltn(
        &cl,
        Box::new(Circle::new(1.0, (0.0, 0.0)).unwrap()),
        Box::new(Rectangular::new(1.4, 0.8, 0.0, (0.0, 0.0)).unwrap()),
        true,
        16,
        25,
    );
    assert_sltn_valid(&tracer, 6.0);
}

/// 圆→超椭圆，轴对称。
#[test]
#[ignore = "SLTN shape test, run with -- --ignored"]
fn test_sltn_circle_to_superellipse() {
    let cl = run_otn_axisymmetric();
    let tracer = run_sltn(
        &cl,
        Box::new(Circle::new(1.0, (0.0, 0.0)).unwrap()),
        Box::new(SuperEllipse::new(0.7, 0.6, 2.5, 0.0, (0.0, 0.0)).unwrap()),
        true,
        16,
        25,
    );
    assert_sltn_valid(&tracer, 6.0);
}

/// 矩形→椭圆，轴对称。
#[test]
#[ignore = "SLTN shape test, run with -- --ignored"]
fn test_sltn_rectangular_to_ellipse() {
    let cl = run_otn_axisymmetric();
    let tracer = run_sltn(
        &cl,
        Box::new(Rectangular::new(1.6, 0.8, 0.0, (0.0, 0.0)).unwrap()),
        Box::new(Ellipse::new(0.7, 0.5, 0.0, (0.0, 0.0)).unwrap()),
        true,
        16,
        25,
    );
    assert_sltn_valid(&tracer, 6.0);
}

/// 超椭圆→矩形，轴对称。
#[test]
#[ignore = "SLTN shape test, run with -- --ignored"]
fn test_sltn_superellipse_to_rectangular() {
    let cl = run_otn_axisymmetric();
    let tracer = run_sltn(
        &cl,
        Box::new(SuperEllipse::new(0.85, 0.85, 2.5, 0.0, (0.0, 0.0)).unwrap()),
        Box::new(Rectangular::new(1.4, 0.7, 0.0, (0.0, 0.0)).unwrap()),
        true,
        20,
        40,
    );
    assert_sltn_valid(&tracer, 6.0);
}

// ═══════════════════════════════════════════════════════════════════════
// 平面基准流场
// ═══════════════════════════════════════════════════════════════════════

/// 平面流场 + 圆形截面。
#[test]
#[ignore = "SLTN shape test, run with -- --ignored"]
fn test_sltn_circle_planar() {
    let cl = run_otn_planar();
    let tracer = run_sltn(
        &cl,
        Box::new(Circle::new(1.0, (0.0, 0.0)).unwrap()),
        Box::new(Circle::new(0.7, (0.0, 0.0)).unwrap()),
        false,
        8,
        15,
    );
    assert_sltn_valid(&tracer, 6.0);
}

/// 平面流场 + 超椭圆截面。
#[test]
#[ignore = "SLTN shape test, run with -- --ignored"]
fn test_sltn_superellipse_planar() {
    let cl = run_otn_planar();
    let tracer = run_sltn(
        &cl,
        Box::new(SuperEllipse::new(0.9, 0.9, 3.0, 0.0, (0.0, 0.0)).unwrap()),
        Box::new(SuperEllipse::new(0.6, 0.8, 3.0, 0.0, (0.0, 0.0)).unwrap()),
        false,
        16,
        30,
    );
    assert_sltn_valid(&tracer, 6.0);
}

/// 平面流场 + 矩形截面。
#[test]
#[ignore = "SLTN shape test, run with -- --ignored"]
fn test_sltn_rectangular_planar() {
    let cl = run_otn_planar();
    let tracer = run_sltn(
        &cl,
        Box::new(Rectangular::new(1.4, 0.8, 0.0, (0.0, 0.0)).unwrap()),
        Box::new(Rectangular::new(1.0, 1.0, 0.0, (0.0, 0.0)).unwrap()),
        false,
        16,
        30,
    );
    assert_sltn_valid(&tracer, 6.0);
}

// ═══════════════════════════════════════════════════════════════════════
// 特殊截面参数
// ═══════════════════════════════════════════════════════════════════════

/// 旋转矩形 (45°)。
#[test]
#[ignore = "SLTN shape test, run with -- --ignored"]
fn test_sltn_rotated_rectangular() {
    let cl = run_otn_axisymmetric();
    let tracer = run_sltn(
        &cl,
        Box::new(Rectangular::new(1.2, 0.8, 45.0_f64.to_radians(), (0.0, 0.0)).unwrap()),
        Box::new(Rectangular::new(0.8, 0.8, 0.0, (0.0, 0.0)).unwrap()),
        true,
        20,
        30,
    );
    assert_sltn_valid(&tracer, 6.0);
}

/// 旋转椭圆 (30°)。
#[test]
#[ignore = "SLTN shape test, run with -- --ignored"]
fn test_sltn_rotated_ellipse() {
    let cl = run_otn_axisymmetric();
    let tracer = run_sltn(
        &cl,
        Box::new(Ellipse::new(0.9, 0.6, 30.0_f64.to_radians(), (0.0, 0.0)).unwrap()),
        Box::new(Ellipse::new(0.6, 0.5, 0.0, (0.0, 0.0)).unwrap()),
        true,
        16,
        25,
    );
    assert_sltn_valid(&tracer, 6.0);
}

/// 偏心截面（进口中心 (0, 0.2)）。
#[test]
#[ignore = "SLTN shape test, run with -- --ignored"]
fn test_sltn_offset_center() {
    let cl = run_otn_axisymmetric();
    let tracer = run_sltn(
        &cl,
        Box::new(Circle::new(0.8, (0.0, 0.2)).unwrap()),
        Box::new(Rectangular::new(1.0, 0.6, 0.0, (0.0, 0.1)).unwrap()),
        true,
        16,
        25,
    );
    assert_sltn_valid(&tracer, 6.0);
}

// ═══════════════════════════════════════════════════════════════════════
// 分辨率变化
// ═══════════════════════════════════════════════════════════════════════

/// 低分辨率 (n_theta=8, n_axis=15)。
#[test]
#[ignore = "SLTN shape test, run with -- --ignored"]
fn test_sltn_low_resolution() {
    let cl = run_otn_axisymmetric();
    let tracer = run_sltn(
        &cl,
        Box::new(Circle::new(1.0, (0.0, 0.0)).unwrap()),
        Box::new(Ellipse::new(0.7, 0.5, 0.0, (0.0, 0.0)).unwrap()),
        true,
        8,
        15,
    );
    assert_sltn_valid(&tracer, 6.0);
}

/// 中分辨率 (n_theta=24, n_axis=50)。
#[test]
#[ignore = "SLTN shape test, run with -- --ignored"]
fn test_sltn_medium_resolution() {
    let cl = run_otn_axisymmetric();
    let tracer = run_sltn(
        &cl,
        Box::new(Circle::new(1.0, (0.0, 0.0)).unwrap()),
        Box::new(Ellipse::new(0.7, 0.5, 0.0, (0.0, 0.0)).unwrap()),
        true,
        24,
        50,
    );
    assert_sltn_valid(&tracer, 6.0);
}

/// 高分辨率 (n_theta=66, n_axis=111) — 默认配置。
#[test]
#[ignore = "SLTN shape test, run with -- --ignored"]
fn test_sltn_high_resolution() {
    let cl = run_otn_axisymmetric();
    let tracer = run_sltn(
        &cl,
        Box::new(Circle::new(1.0, (0.0, 0.0)).unwrap()),
        Box::new(Ellipse::new(0.7, 0.55, 0.0, (0.0, 0.0)).unwrap()),
        true,
        66,
        111,
    );
    assert_sltn_valid(&tracer, 6.0);
}
