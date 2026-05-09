/// 喷管完整集成测试。
///
/// 这些测试运行完整的喷管特征线法计算流程，耗时较长。
/// 仅验证最终输出（点有效性、出口长度），不测试内部状态。
use aero::{
    nozzle::{ConstraintNozzle, Control, Geometry, Inlet, NozzleConfig, Outlet, Throat, IO},
    Material,
};

/// 自动选择初始膨胀角，验证完整喷管流程。
#[test]
#[ignore = "slow integration test, run with -- --ignored"]
fn test_full_nozzle_auto_theta() {
    let config = NozzleConfig {
        control: Control::default(),
        material: Material::air_piecewise_polynomial(),
        inlet: Inlet::default(),
        geometry: Geometry::default(),
        throat: Throat {
            radius_throat: 0.0,
            theta_a: f64::NAN, // 自动选择
        },
        outlet: Outlet::default(),
        io: IO::default(),
    };
    let target_len = config.geometry.length;
    let mut n = ConstraintNozzle::new_otn(config);
    n.run();
    let lines = n.get_assembly_charlines();

    assert!(!lines.is_empty(), "应产生流场数据");

    for (li, line) in lines.iter().enumerate() {
        for (pi, point) in line.iter().enumerate() {
            assert!(point.is_valid(), "无效点: line={li} pt={pi}",);
            assert!(
                point.x >= -1e-9 && point.y >= -1e-9,
                "坐标越界: line={li} pt={pi}: x={}, y={}",
                point.x,
                point.y
            );
        }
    }

    let mut max_x: f64 = 0.0;
    for line in lines.iter() {
        for point in line.iter() {
            max_x = max_x.max(point.x);
            assert!(
                point.x < target_len + 1.0,
                "x={} 超出长度 {target_len}",
                point.x
            );
        }
    }
    assert!(
        max_x >= target_len - 1e-4,
        "出口未到达目标长度: max_x={:.6}, target={:.6}",
        max_x,
        target_len
    );
}

/// 固定膨胀角，验证含膨胀段的完整喷管。
#[test]
#[ignore = "slow integration test, run with -- --ignored"]
fn test_full_nozzle_fixed_theta() {
    let config = NozzleConfig {
        control: Control::default(),
        material: Material::from_rgas_gamma(287.042, 1.4),
        inlet: Inlet::default(),
        geometry: Geometry::default(),
        throat: Throat {
            radius_throat: 0.0,
            theta_a: 19.0_f64.to_radians(),
        },
        outlet: Outlet::default(),
        io: IO::default(),
    };
    let target_len = config.geometry.length;
    let mut n = ConstraintNozzle::new_otn(config);
    n.run();
    let lines = n.get_assembly_charlines();

    assert!(!lines.is_empty(), "应产生流场数据");

    for (li, line) in lines.iter().enumerate() {
        for (pi, point) in line.iter().enumerate() {
            assert!(point.is_valid(), "无效点: line={li} pt={pi}");
            assert!(
                point.x >= -1e-9 && point.y >= -1e-9,
                "坐标越界: line={li} pt={pi}: x={}, y={}",
                point.x,
                point.y
            );
        }
    }

    for line in lines.iter() {
        for point in line.iter() {
            assert!(
                point.x < target_len + 1.0,
                "x={} 超出长度 {target_len}",
                point.x
            );
        }
    }
}
