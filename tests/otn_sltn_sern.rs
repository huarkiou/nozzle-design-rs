/// 项目级继承集成测试：OTN → SLTN 完整工作流（SERN 设计）。
///
/// 这些测试验证从最优推力喷管 (OTN) 生成基准流场，到三维流线追踪喷管 (SLTN)
/// 设计的完整流程。覆盖常见超燃冲压发动机 SERN (Single Expansion Ramp Nozzle)
/// 喷管设计配置。
///
/// SERN 设计要点：
/// - 超燃冲压发动机排气喷管通常具有非圆形截面（矩形、超椭圆等）
/// - 进口截面面积较小，出口截面面积大（高膨胀比）
/// - 支持轴对称与平面两种基准流场类型
///
/// 所有测试标记为 `#[ignore]`，使用 `cargo test -- --ignored` 运行。
mod common;

use common::{
    assert_sltn_valid, run_otn_axisymmetric, run_otn_planar, run_sltn, run_sltn_with_options,
};
use geometry::{Circle, Ellipse, Rectangular, SuperEllipse};

// ═══════════════════════════════════════════════════════════════════════
// SERN 设计：轴对称 OTN → 各种 SLTN 截面组合
// ═══════════════════════════════════════════════════════════════════════

/// 经典圆形→椭圆 SERN（轴对称）。
#[test]
#[ignore = "full OTN→SLTN SERN, run with -- --ignored"]
fn test_sern_circle_ellipse_axisymmetric() {
    let cl = run_otn_axisymmetric();
    let tracer = run_sltn(
        &cl,
        Box::new(Circle::new(0.9, (0.0, 0.0)).unwrap()),
        Box::new(Ellipse::new(0.7, 0.55, 0.0, (0.0, 0.0)).unwrap()),
        true,
        12,
        20,
    );
    assert_sltn_valid(&tracer, 6.0);
}

/// 圆形→圆形（收缩出口），验证逆向追踪。
#[test]
#[ignore = "full OTN→SLTN SERN, run with -- --ignored"]
fn test_sern_circle_to_circle_shrink() {
    let cl = run_otn_axisymmetric();
    let tracer = run_sltn(
        &cl,
        Box::new(Circle::new(1.0, (0.0, 0.0)).unwrap()),
        Box::new(Circle::new(0.6, (0.0, 0.0)).unwrap()),
        true,
        8,
        15,
    );
    assert_sltn_valid(&tracer, 6.0);
}

/// 矩形→矩形 SERN（典型二维超燃冲压喷管）。
#[test]
#[ignore = "full OTN→SLTN SERN, run with -- --ignored"]
fn test_sern_rectangular() {
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

/// 平面流场 + 超椭圆截面 SERN。
#[test]
#[ignore = "full OTN→SLTN SERN, run with -- --ignored"]
fn test_sern_superellipse_planar() {
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

/// 超椭圆→矩形混合 SERN。
#[test]
#[ignore = "full OTN→SLTN SERN, run with -- --ignored"]
fn test_sern_superellipse_to_rect() {
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

/// 旋转矩形 SERN（45°）。
#[test]
#[ignore = "full OTN→SLTN SERN, run with -- --ignored"]
fn test_sern_rotated_rect_planar() {
    let cl = run_otn_planar();
    let tracer = run_sltn(
        &cl,
        Box::new(Rectangular::new(1.2, 0.8, 45.0_f64.to_radians(), (0.0, 0.0)).unwrap()),
        Box::new(Rectangular::new(0.8, 0.8, 0.0, (0.0, 0.0)).unwrap()),
        false,
        20,
        30,
    );
    assert_sltn_valid(&tracer, 6.0);
}

/// 偏心截面 SERN。
#[test]
#[ignore = "full OTN→SLTN SERN, run with -- --ignored"]
fn test_sern_offset_centers() {
    let cl = run_otn_axisymmetric();
    let tracer = run_sltn(
        &cl,
        Box::new(Ellipse::new(0.8, 0.7, 0.0, (0.0, 0.2)).unwrap()),
        Box::new(Rectangular::new(1.0, 0.6, 0.0, (0.0, 0.1)).unwrap()),
        true,
        16,
        25,
    );
    assert_sltn_valid(&tracer, 6.0);
}

// ═══════════════════════════════════════════════════════════════════════
// 鲁棒性测试
// ═══════════════════════════════════════════════════════════════════════

/// monotonic=true 强制单调。
#[test]
#[ignore = "full OTN→SLTN SERN, run with -- --ignored"]
fn test_sern_monotonic() {
    let cl = run_otn_axisymmetric();
    let tracer = run_sltn_with_options(
        &cl,
        Box::new(Circle::new(1.0, (0.0, 0.0)).unwrap()),
        Box::new(Rectangular::new(1.4, 0.8, 0.0, (0.0, 0.0)).unwrap()),
        true,
        12,
        20,
        true,
        0.0,
    );
    for (i, profile) in tracer.model.iter().enumerate() {
        for w in profile.windows(2) {
            assert!(
                w[0].x < w[1].x,
                "profile {i}: x not monotonic with monotonic=true"
            );
        }
    }
}

/// 加权过渡函数 (a=2.0)。
#[test]
#[ignore = "full OTN→SLTN SERN, run with -- --ignored"]
fn test_sern_weighted() {
    let cl = run_otn_axisymmetric();
    let tracer = run_sltn_with_options(
        &cl,
        Box::new(Circle::new(0.9, (0.0, 0.0)).unwrap()),
        Box::new(Ellipse::new(0.6, 0.5, 0.0, (0.0, 0.0)).unwrap()),
        true,
        16,
        25,
        false,
        2.0,
    );
    assert_sltn_valid(&tracer, 6.0);
}
