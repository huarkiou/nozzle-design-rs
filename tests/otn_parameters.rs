/// OTN 喷管参数覆盖集成测试。
///
/// 覆盖特征线法控制参数、几何约束、材料模型、进口条件、喉部参数
/// 等独立维度的变化，验证各参数组合下 OTN 计算正常产出有效流场。
///
/// 所有测试标记为 `#[ignore]`，使用 `cargo test -- --ignored` 运行。
mod common;

use aero::{
    Material,
    nozzle::{Control, Geometry, IO, Inlet, NozzleConfig, Outlet, Throat},
};
use common::{assert_charlines_valid, run_otn_custom};

// ═══════════════════════════════════════════════════════════════════════
// 流动类型
// ═══════════════════════════════════════════════════════════════════════

/// 轴对称流动，默认参数。
#[test]
#[ignore = "OTN param test, run with -- --ignored"]
fn test_otn_axisymmetric_default() {
    let cl = run_otn_custom(
        true,
        Material::air_nasa9piecewise_polynomial(),
        1.0,
        6.0,
        1.2,
        800_000.0,
        2000.0,
        19.0_f64.to_radians(),
        7000.0,
    );
    assert_charlines_valid(&cl);
}

/// 平面流动，默认参数。
#[test]
#[ignore = "OTN param test, run with -- --ignored"]
fn test_otn_planar_default() {
    let cl = run_otn_custom(
        false,
        Material::air_nasa9piecewise_polynomial(),
        1.0,
        6.0,
        1.2,
        800_000.0,
        2000.0,
        19.0_f64.to_radians(),
        7000.0,
    );
    assert_charlines_valid(&cl);
}

// ═══════════════════════════════════════════════════════════════════════
// 材料模型
// ═══════════════════════════════════════════════════════════════════════

/// NASA 9 系数变比热空气模型。
#[test]
#[ignore = "OTN param test, run with -- --ignored"]
fn test_otn_material_nasa9() {
    let cl = run_otn_custom(
        true,
        Material::air_nasa9piecewise_polynomial(),
        1.0,
        6.0,
        1.2,
        800_000.0,
        2000.0,
        19.0_f64.to_radians(),
        7000.0,
    );
    assert_charlines_valid(&cl);
}

/// 常数比热容空气模型（γ=1.4）。
#[test]
#[ignore = "OTN param test, run with -- --ignored"]
fn test_otn_material_constant_cp() {
    let cl = run_otn_custom(
        true,
        Material::from_rgas_gamma(287.042, 1.4),
        1.0,
        6.0,
        1.2,
        800_000.0,
        2000.0,
        19.0_f64.to_radians(),
        7000.0,
    );
    assert_charlines_valid(&cl);
}

// ═══════════════════════════════════════════════════════════════════════
// 几何参数
// ═══════════════════════════════════════════════════════════════════════

/// 短喷管 (length=3)。
#[test]
#[ignore = "OTN param test, run with -- --ignored"]
fn test_otn_geometry_short() {
    let cl = run_otn_custom(
        true,
        Material::air_nasa9piecewise_polynomial(),
        1.0,
        3.0,
        1.2,
        800_000.0,
        2000.0,
        19.0_f64.to_radians(),
        7000.0,
    );
    assert_charlines_valid(&cl);
}

/// 长喷管 (length=10)。
#[test]
#[ignore = "OTN param test, run with -- --ignored"]
fn test_otn_geometry_long() {
    let cl = run_otn_custom(
        true,
        Material::air_nasa9piecewise_polynomial(),
        1.0,
        10.0,
        1.2,
        800_000.0,
        2000.0,
        19.0_f64.to_radians(),
        7000.0,
    );
    assert_charlines_valid(&cl);
}

/// 小进口高度 (height=0.5)。
#[test]
#[ignore = "OTN param test, run with -- --ignored"]
fn test_otn_geometry_small_inlet() {
    let cl = run_otn_custom(
        true,
        Material::air_nasa9piecewise_polynomial(),
        0.5,
        6.0,
        1.2,
        800_000.0,
        2000.0,
        19.0_f64.to_radians(),
        7000.0,
    );
    assert_charlines_valid(&cl);
}

// ═══════════════════════════════════════════════════════════════════════
// 进口条件
// ═══════════════════════════════════════════════════════════════════════

/// 低马赫数 (Ma=1.05)。
#[test]
#[ignore = "OTN param test, run with -- --ignored"]
fn test_otn_inlet_low_ma() {
    let cl = run_otn_custom(
        true,
        Material::air_nasa9piecewise_polynomial(),
        1.0,
        6.0,
        1.05,
        800_000.0,
        2000.0,
        19.0_f64.to_radians(),
        7000.0,
    );
    assert_charlines_valid(&cl);
}

/// 高马赫数 (Ma=3.0)。
#[test]
#[ignore = "OTN param test, run with -- --ignored"]
fn test_otn_inlet_high_ma() {
    let cl = run_otn_custom(
        true,
        Material::air_nasa9piecewise_polynomial(),
        1.0,
        6.0,
        3.0,
        800_000.0,
        2000.0,
        19.0_f64.to_radians(),
        7000.0,
    );
    assert_charlines_valid(&cl);
}

/// 高总压 (10 MPa)。
#[test]
#[ignore = "OTN param test, run with -- --ignored"]
fn test_otn_inlet_high_pressure() {
    let cl = run_otn_custom(
        true,
        Material::air_nasa9piecewise_polynomial(),
        1.0,
        6.0,
        1.2,
        10_000_000.0,
        2000.0,
        19.0_f64.to_radians(),
        7000.0,
    );
    assert_charlines_valid(&cl);
}

/// 高总温 (3500 K)。
#[test]
#[ignore = "OTN param test, run with -- --ignored"]
fn test_otn_inlet_high_temperature() {
    let cl = run_otn_custom(
        true,
        Material::air_nasa9piecewise_polynomial(),
        1.0,
        6.0,
        1.2,
        800_000.0,
        3500.0,
        19.0_f64.to_radians(),
        7000.0,
    );
    assert_charlines_valid(&cl);
}

// ═══════════════════════════════════════════════════════════════════════
// 喉部参数
// ═══════════════════════════════════════════════════════════════════════

/// 小膨胀角 (θ_a=5°)。
#[test]
#[ignore = "OTN param test, run with -- --ignored"]
fn test_otn_throat_small_theta() {
    let cl = run_otn_custom(
        true,
        Material::air_nasa9piecewise_polynomial(),
        1.0,
        6.0,
        1.2,
        800_000.0,
        2000.0,
        5.0_f64.to_radians(),
        7000.0,
    );
    assert_charlines_valid(&cl);
}

/// 大膨胀角 (θ_a=30°)。
#[test]
#[ignore = "OTN param test, run with -- --ignored"]
fn test_otn_throat_large_theta() {
    let cl = run_otn_custom(
        true,
        Material::air_nasa9piecewise_polynomial(),
        1.0,
        6.0,
        1.2,
        800_000.0,
        2000.0,
        30.0_f64.to_radians(),
        7000.0,
    );
    assert_charlines_valid(&cl);
}

/// 非零喉部过渡圆弧半径 (R_t=0.5)。
#[test]
#[ignore = "OTN param test, run with -- --ignored"]
fn test_otn_throat_nonzero_radius() {
    let config = NozzleConfig {
        control: Control {
            axisymmetric: true,
            ..Control::default()
        },
        material: Material::air_nasa9piecewise_polynomial(),
        inlet: Inlet::default(),
        geometry: Geometry::default(),
        throat: Throat {
            radius_throat: 0.5,
            theta_a: 19.0_f64.to_radians(),
        },
        outlet: Outlet::default(),
        io: IO::default(),
    };
    let mut nozzle = aero::nozzle::ConstraintNozzle::new_otn(config);
    nozzle.run();
    assert_charlines_valid(&nozzle.get_assembly_charlines());
}

// ═══════════════════════════════════════════════════════════════════════
// 出口条件
// ═══════════════════════════════════════════════════════════════════════

/// 低背压 (1000 Pa)。
#[test]
#[ignore = "OTN param test, run with -- --ignored"]
fn test_otn_outlet_low_backpressure() {
    let cl = run_otn_custom(
        true,
        Material::air_nasa9piecewise_polynomial(),
        1.0,
        6.0,
        1.2,
        800_000.0,
        2000.0,
        19.0_f64.to_radians(),
        1000.0,
    );
    assert_charlines_valid(&cl);
}

/// 高背压 (50000 Pa)。
#[test]
#[ignore = "OTN param test, run with -- --ignored"]
fn test_otn_outlet_high_backpressure() {
    let cl = run_otn_custom(
        true,
        Material::air_nasa9piecewise_polynomial(),
        1.0,
        6.0,
        1.2,
        800_000.0,
        2000.0,
        30.0_f64.to_radians(),
        50_000.0,
    );
    assert_charlines_valid(&cl);
}
