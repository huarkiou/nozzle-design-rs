/// 项目级继承集成测试：OTN → SLTN 完整工作流。
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
/// 所有测试标记为 `#[ignore]`，因为完整 OTN 计算耗时较长。
/// 使用 `cargo test -- --ignored` 运行。
use aero::{
    Material,
    moc::CharLines,
    nozzle::{ConstraintNozzle, Control, Geometry, IO, Inlet, NozzleConfig, Outlet, Throat},
    streamline_trace::{StreamlineConfig, StreamlineTrace},
};
use geometry::{Circle, ClosedCurve, Ellipse, Rectangular, SuperEllipse};

// ── 辅助函数 ──────────────────────────────────────────────────────────

/// 运行标准轴对称 OTN 喷管（自动选择初始膨胀角），返回特征线流场数据。
///
/// 使用 NASA 9 系数变比热空气模型以获得更真实的物理模拟。
/// 喷管几何：进口高度 1 m，目标长度 6 m，出口高度不受约束（最大推力模式）。
fn run_otn_axisymmetric() -> CharLines {
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
            theta_a: f64::NAN, // 自动迭代选择
        },
        outlet: Outlet::default(),
        io: IO::default(),
    };
    let mut nozzle = ConstraintNozzle::new_otn(config);
    nozzle.run();
    nozzle.get_assembly_charlines()
}

/// 运行固定膨胀角的平面 OTN 喷管，返回特征线流场数据。
///
/// 用于平面 SERN 喷管设计，固定膨胀角减少计算不确定性。
fn run_otn_planar() -> CharLines {
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

/// 验证 SLTN 产出的模型满足基本正确性条件。
///
/// 检查项：
/// - 每个方位角剖面点数 ≥ 2
/// - 每个剖面的 x 坐标严格单调递增
/// - 进口 x ≈ 0，出口 x ≈ 目标长度
/// - 顺流/逆流剖面均非空
fn assert_sltn_valid(tracer: &StreamlineTrace, expected_len: f64) {
    assert!(!tracer.model.is_empty(), "SLTN model should not be empty");
    assert!(
        !tracer.downstream.is_empty(),
        "downstream profiles should not be empty"
    );
    assert!(
        !tracer.upstream.is_empty(),
        "upstream profiles should not be empty"
    );

    for (i, profile) in tracer.model.iter().enumerate() {
        assert!(
            profile.len() >= 2,
            "profile {i}: expected ≥2 points, got {}",
            profile.len()
        );

        // x 严格单调递增
        for w in profile.windows(2) {
            assert!(
                w[0].x < w[1].x,
                "profile {i}: x not monotonic: {} >= {}",
                w[0].x,
                w[1].x
            );
        }

        // 进口 x ≈ 0
        let first_x = profile.first().unwrap().x;
        assert!(
            first_x.abs() < 0.5,
            "profile {i}: first x = {first_x}, expected ≈ 0"
        );

        // 出口 x ≈ 目标长度
        let last_x = profile.last().unwrap().x;
        assert!(
            (last_x - expected_len).abs() < expected_len * 0.15,
            "profile {i}: last x = {last_x}, expected ≈ {expected_len}"
        );
    }
}

/// 运行 SLTN 并返回 tracer。
fn run_sltn(
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

// ═══════════════════════════════════════════════════════════════════════
// 测试用例：OTN → SLTN 完整继承工作流
// ═══════════════════════════════════════════════════════════════════════

// ── SERN 设计 1: 轴对称 OTN → 圆形进口 × 椭圆形出口 ──

/// 经典轴对称喷管 + 椭圆出口 SERN 设计。
///
/// 进口圆形 (r=0.9)，出口椭圆 (a=0.7, b=0.55)。
/// 模拟超燃冲压发动机中圆形进口过渡到扁椭圆出口的常见场景。
#[test]
#[ignore = "full OTN→SLTN integration, run with -- --ignored"]
fn test_otn_to_sltn_circle_ellipse_axisymmetric() {
    let charlines = run_otn_axisymmetric();
    assert!(!charlines.is_empty(), "OTN should produce charlines");

    let inlet = Circle::new(0.9, (0.0, 0.0)).expect("valid circle");
    let outlet = Ellipse::new(0.7, 0.55, 0.0, (0.0, 0.0)).expect("valid ellipse");

    let tracer = run_sltn(&charlines, Box::new(inlet), Box::new(outlet), true, 12, 20);
    assert_sltn_valid(&tracer, 6.0);
}

// ── SERN 设计 2: 轴对称 OTN → 圆形进口 × 圆形出口（等比例缩小） ──

/// 圆形进口到更小圆形出口。
/// 验证收缩出口截面的逆向流线追踪正常工作。
#[test]
#[ignore = "full OTN→SLTN integration, run with -- --ignored"]
fn test_otn_to_sltn_circle_to_circle_shrink() {
    let charlines = run_otn_axisymmetric();
    assert!(!charlines.is_empty(), "OTN should produce charlines");

    let inlet = Circle::new(1.0, (0.0, 0.0)).expect("valid circle");
    let outlet = Circle::new(0.6, (0.0, 0.0)).expect("valid circle");

    let tracer = run_sltn(&charlines, Box::new(inlet), Box::new(outlet), true, 8, 15);
    assert_sltn_valid(&tracer, 6.0);
}

// ── SERN 设计 3: 轴对称 OTN → 矩形进口 × 矩形出口 ──

/// 矩形截面 SERN 设计（典型的二维超燃冲压发动机排气喷管）。
///
/// 矩形进口 1.6×0.8，矩形出口 1.2×1.0。
/// 矩形截面是超燃冲压发动机最常用的喷管截面类型。
#[test]
#[ignore = "full OTN→SLTN integration, run with -- --ignored"]
fn test_otn_to_sltn_rectangular_sern() {
    let charlines = run_otn_axisymmetric();
    assert!(!charlines.is_empty(), "OTN should produce charlines");

    let inlet = Rectangular::new(1.6, 0.8, 0.0, (0.0, 0.0)).expect("valid rect");
    let outlet = Rectangular::new(1.2, 1.0, 0.0, (0.0, 0.0)).expect("valid rect");

    let tracer = run_sltn(&charlines, Box::new(inlet), Box::new(outlet), true, 16, 30);
    assert_sltn_valid(&tracer, 6.0);
}

// ── SERN 设计 4: 平面 OTN → 超椭圆截面 ──

/// 平面基准流场 + 超椭圆截面 SERN 设计。
///
/// 超椭圆 (n=3.0) 形状近似为圆角矩形，是超燃冲压发动机喷管中
/// 兼顾结构强度与气动性能的常用截面形状。
#[test]
#[ignore = "full OTN→SLTN integration, run with -- --ignored"]
fn test_otn_to_sltn_superellipse_planar() {
    let charlines = run_otn_planar();
    assert!(!charlines.is_empty(), "planar OTN should produce charlines");

    let inlet = SuperEllipse::new(0.9, 0.9, 3.0, 0.0, (0.0, 0.0)).expect("valid superellipse");
    let outlet = SuperEllipse::new(0.6, 0.8, 3.0, 0.0, (0.0, 0.0)).expect("valid superellipse");

    let tracer = run_sltn(&charlines, Box::new(inlet), Box::new(outlet), false, 16, 30);
    assert_sltn_valid(&tracer, 6.0);
}

// ── SERN 设计 5: 轴对称 OTN → 超椭圆进口 × 矩形出口 ──

/// 混合形状 SERN 设计。
///
/// 进口超椭圆 (n=2.5) 过渡到矩形出口，模拟超燃冲压发动机中
/// 进口圆形/超椭圆过渡到矩形出口截面。
#[test]
#[ignore = "full OTN→SLTN integration, run with -- --ignored"]
fn test_otn_to_sltn_superellipse_to_rect() {
    let charlines = run_otn_axisymmetric();
    assert!(!charlines.is_empty(), "OTN should produce charlines");

    let inlet = SuperEllipse::new(0.85, 0.85, 2.5, 0.0, (0.0, 0.0)).expect("valid superellipse");
    let outlet = Rectangular::new(1.4, 0.7, 0.0, (0.0, 0.0)).expect("valid rect");

    let tracer = run_sltn(&charlines, Box::new(inlet), Box::new(outlet), true, 20, 40);
    assert_sltn_valid(&tracer, 6.0);
}

// ── SERN 设计 6: 平面 OTN → 矩形（旋转45°） ──

/// 旋转矩形截面 SERN 设计。
///
/// 进口矩形旋转 45°，出口矩形不旋转。
/// 验证旋转截面在流线追踪中的正确性。
#[test]
#[ignore = "full OTN→SLTN integration, run with -- --ignored"]
fn test_otn_to_sltn_rotated_rect_planar() {
    let charlines = run_otn_planar();
    assert!(!charlines.is_empty(), "planar OTN should produce charlines");

    let inlet = Rectangular::new(1.2, 0.8, 45.0_f64.to_radians(), (0.0, 0.0)).expect("valid rect");
    let outlet = Rectangular::new(0.8, 0.8, 0.0, (0.0, 0.0)).expect("valid rect");

    let tracer = run_sltn(&charlines, Box::new(inlet), Box::new(outlet), false, 20, 30);
    assert_sltn_valid(&tracer, 6.0);
}

// ── SERN 设计 7: 轴对称 OTN → 偏心椭圆进口 × 偏心矩形出口 ──

/// 偏心截面 SERN 设计。
///
/// 进口椭圆中心偏移 (0, 0.2)，出口矩形中心偏移 (0, 0.1)。
/// 验证截面中心不重合时的流线追踪正确性。
#[test]
#[ignore = "full OTN→SLTN integration, run with -- --ignored"]
fn test_otn_to_sltn_offset_centers() {
    let charlines = run_otn_axisymmetric();
    assert!(!charlines.is_empty(), "OTN should produce charlines");

    let inlet = Ellipse::new(0.8, 0.7, 0.0, (0.0, 0.2)).expect("valid ellipse");
    let outlet = Rectangular::new(1.0, 0.6, 0.0, (0.0, 0.1)).expect("valid rect");

    let tracer = run_sltn(&charlines, Box::new(inlet), Box::new(outlet), true, 16, 25);
    assert_sltn_valid(&tracer, 6.0);
}

// ── 鲁棒性测试 ────────────────────────────────────────────────────────

/// 单调性标志开启的 SLTN 测试。
///
/// 当 monotonic=true 时，SLTN 会在融合顺流/逆流剖面时
/// 强制抹平 x 坐标非单调区域。
#[test]
#[ignore = "full OTN→SLTN integration, run with -- --ignored"]
fn test_otn_to_sltn_monotonic() {
    let charlines = run_otn_axisymmetric();
    assert!(!charlines.is_empty(), "OTN should produce charlines");

    let inlet = Circle::new(1.0, (0.0, 0.0)).expect("valid circle");
    let outlet = Rectangular::new(1.4, 0.8, 0.0, (0.0, 0.0)).expect("valid rect");

    let config = StreamlineConfig {
        axisymmetric: true,
        datasource_inlet: charlines.clone(),
        datasource_outlet: charlines.clone(),
        inlet_shape: Box::new(inlet),
        outlet_shape: Box::new(outlet),
        monotonic: true,
        weight_parameter_a: 0.0,
    };
    let mut tracer = StreamlineTrace::new(config);
    tracer.run(12, 20).expect("SLTN monotonic should succeed");

    // monotonic=true 时每个剖面必须严格单调
    for (i, profile) in tracer.model.iter().enumerate() {
        for w in profile.windows(2) {
            assert!(
                w[0].x < w[1].x,
                "profile {i}: x not monotonic with monotonic=true: {} >= {}",
                w[0].x,
                w[1].x
            );
        }
    }
}

/// 加权过渡函数测试。
///
/// 使用 weight_parameter_a=2.0 的非线性加权融合顺流/逆流剖面。
#[test]
#[ignore = "full OTN→SLTN integration, run with -- --ignored"]
fn test_otn_to_sltn_weighted() {
    let charlines = run_otn_axisymmetric();
    assert!(!charlines.is_empty(), "OTN should produce charlines");

    let inlet = Circle::new(0.9, (0.0, 0.0)).expect("valid circle");
    let outlet = Ellipse::new(0.6, 0.5, 0.0, (0.0, 0.0)).expect("valid ellipse");

    let config = StreamlineConfig {
        axisymmetric: true,
        datasource_inlet: charlines.clone(),
        datasource_outlet: charlines.clone(),
        inlet_shape: Box::new(inlet),
        outlet_shape: Box::new(outlet),
        monotonic: false,
        weight_parameter_a: 2.0,
    };
    let mut tracer = StreamlineTrace::new(config);
    tracer.run(16, 25).expect("SLTN weighted should succeed");
    assert_sltn_valid(&tracer, 6.0);
}

/// 高分辨率 SLTN 测试 — 验证高密度网格下的稳定性。
#[test]
#[ignore = "full OTN→SLTN integration, high-resolution, run with -- --ignored"]
fn test_otn_to_sltn_high_resolution() {
    let charlines = run_otn_axisymmetric();
    assert!(!charlines.is_empty(), "OTN should produce charlines");

    let inlet = Circle::new(1.0, (0.0, 0.0)).expect("valid circle");
    let outlet = Ellipse::new(0.7, 0.55, 0.0, (0.0, 0.0)).expect("valid ellipse");

    // n_theta=66, n_axis=111 为 SLTN 默认高分辨率配置
    let tracer = run_sltn(&charlines, Box::new(inlet), Box::new(outlet), true, 66, 111);
    assert_sltn_valid(&tracer, 6.0);
}
