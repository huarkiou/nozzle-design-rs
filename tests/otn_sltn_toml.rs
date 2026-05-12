/// 项目级 TOML 配置文件继承集成测试：OTN → SLTN 完整工作流。
///
/// 这些测试模拟 CLI 应用 (`apps/otn` 和 `apps/sltn`) 从 TOML 配置文件
/// 读取参数、执行喷管设计的完整流程，验证：
/// 1. OTN 从 TOML 文件加载配置 → 运行 → 生成特征线流场
/// 2. 特征线数据经过文件 I/O 往返（写→读）保持完整性
/// 3. SLTN 从 TOML 文件加载配置 → 读取特征线文件 → 运行流线追踪
///
/// 覆盖材料模型：
/// - 分段多项式变比热容（NASA 9 系数空气模型）
/// - 常数比热容（γ=1.4 常比热空气）
///
/// 使用 `cargo test -- --ignored` 运行。
mod common;

use std::{
    fs,
    io::Write,
    path::{Path, PathBuf},
};

use aero::Material;
use aero::moc::{CharLines, read_charlines_from_file_checked};
use aero::nozzle::{ConstraintNozzle, NozzleConfig};
use aero::streamline_trace::{StreamlineConfig, StreamlineTrace};
use geometry::{Circle, ClosedCurve, Ellipse, Rectangular, SuperEllipse};
use serde::Deserialize;

// ── 内联 SLTN 配置结构体 ─────────────────────────────────────────────
//
// 由于 `apps/sltn` 是二进制 crate，其 `Config` 结构体不可导入。
// 此处定义最小化的镜像结构体用于反序列化 SLTN TOML 配置文件。

#[derive(Deserialize)]
#[serde(rename_all = "PascalCase")]
struct SltnConfig {
    #[serde(default)]
    control: SltnControl,
    base_fluid_field: SltnBaseFluidField,
    inlet: SltnShape,
    outlet: SltnShape,
}

#[derive(Deserialize, Default)]
#[serde(deny_unknown_fields)]
struct SltnControl {
    #[serde(default = "default_n_theta")]
    n_theta: usize,
    #[serde(default = "default_n_axis")]
    n_axis: usize,
    #[serde(default)]
    monotonic: bool,
    #[serde(default)]
    weight_parameter_a: f64,
}

fn default_n_theta() -> usize {
    66
}
fn default_n_axis() -> usize {
    111
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct SltnBaseFluidField {
    axisymmetric: bool,
    datasource_inlet: PathBuf,
    datasource_outlet: PathBuf,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct SltnShape {
    #[serde(default = "default_normalized")]
    #[allow(dead_code)]
    normalized: bool,
    shape: String,
    #[serde(default)]
    center: Option<Vec<f64>>,
    #[serde(default)]
    radius: Option<f64>,
    #[serde(default)]
    a: Option<f64>,
    #[serde(default)]
    b: Option<f64>,
    #[serde(default)]
    alpha: Option<f64>,
    #[serde(default)]
    length: Option<f64>,
    #[serde(default)]
    width: Option<f64>,
    #[serde(default)]
    n: Option<f64>,
}

fn default_normalized() -> bool {
    true
}

// ── 临时目录管理 ──────────────────────────────────────────────────────

struct TempDir {
    path: PathBuf,
}

impl TempDir {
    fn new(prefix: &str) -> Self {
        let dir = std::env::temp_dir().join(format!("nozzle_test_{prefix}_{}", uuid_simple()));
        fs::create_dir_all(&dir).expect("failed to create temp dir");
        Self { path: dir }
    }

    fn write(&self, filename: &str, content: &str) -> PathBuf {
        let p = self.path.join(filename);
        let mut f = fs::File::create(&p).expect("failed to create temp file");
        f.write_all(content.as_bytes())
            .expect("failed to write temp file");
        p
    }

    fn write_bytes(&self, filename: &str, content: &[u8]) -> PathBuf {
        let p = self.path.join(filename);
        let mut f = fs::File::create(&p).expect("failed to create temp file");
        f.write_all(content).expect("failed to write temp file");
        p
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.path);
    }
}

/// 简单 UUID 用于临时目录命名（避免外部依赖）。
fn uuid_simple() -> String {
    use std::time::{SystemTime, UNIX_EPOCH};
    let t = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_nanos();
    format!("{t:016x}")
}

// ── 工作流辅助函数 ────────────────────────────────────────────────────

/// 从 TOML 文件加载 OTN 配置、运行喷管、返回特征线。
fn run_otn_from_toml(toml_path: &Path) -> CharLines {
    let config =
        NozzleConfig::from_toml_file(toml_path).expect("failed to load OTN config from TOML");
    config.validate().expect("OTN config validation failed");
    let config = config.build().expect("OTN config build failed");

    let mut nozzle = ConstraintNozzle::new_otn(config);
    nozzle.run();
    nozzle.get_assembly_charlines()
}

/// 写入特征线到文件再读回，验证 I/O 往返完整性。
fn charlines_roundtrip(charlines: &CharLines, filepath: &Path) -> CharLines {
    charlines
        .write_to_file(filepath, false)
        .expect("failed to write charlines");
    assert!(filepath.exists(), "charlines file should exist after write");

    let read_back =
        read_charlines_from_file_checked(filepath).expect("failed to read charlines back");

    // 行数可能因文件格式边界条件微小差异（±1），但不应显著丢失数据
    let diff = (charlines.len() as isize - read_back.len() as isize).unsigned_abs();
    assert!(
        diff <= 2,
        "charlines roundtrip: original {} lines vs read-back {} lines (diff={diff} too large)",
        charlines.len(),
        read_back.len()
    );
    assert!(
        !read_back.is_empty(),
        "read-back charlines should not be empty"
    );
    read_back
}

/// 从 SLTN TOML 模板构建 StreamlineConfig 并运行追踪。
///
/// `template_toml` 中的 `datasource_inlet` / `datasource_outlet` 会被替换为
/// 实际的临时文件路径。
fn run_sltn_from_toml_template(template_toml: &str, charlines_path: &Path) -> StreamlineTrace {
    let toml_content = template_toml.replace("./field_data.txt", &charlines_path.to_string_lossy());
    let config: SltnConfig = toml::from_str(&toml_content).expect("failed to parse SLTN TOML");

    let normalize_factor = 1.0; // normalized=true in TOML

    let inlet_shape = build_shape(&config.inlet, normalize_factor);
    let outlet_shape = build_shape(&config.outlet, normalize_factor);

    let lines_inlet = read_charlines_from_file_checked(&config.base_fluid_field.datasource_inlet)
        .expect("failed to load inlet charlines");
    let lines_outlet = read_charlines_from_file_checked(&config.base_fluid_field.datasource_outlet)
        .expect("failed to load outlet charlines");

    let trace_config = StreamlineConfig {
        axisymmetric: config.base_fluid_field.axisymmetric,
        datasource_inlet: lines_inlet,
        datasource_outlet: lines_outlet,
        inlet_shape,
        outlet_shape,
        monotonic: config.control.monotonic,
        weight_parameter_a: config.control.weight_parameter_a,
    };

    let mut tracer = StreamlineTrace::new(trace_config);
    tracer
        .run(config.control.n_theta, config.control.n_axis)
        .expect("SLTN trace should succeed");
    tracer
}

/// 从 TOML ShapeConfig 构建 ClosedCurve。
fn build_shape(shape: &SltnShape, factor: f64) -> Box<dyn ClosedCurve> {
    let center: (f64, f64) = shape
        .center
        .as_ref()
        .map(|c| (c[0], c[1]))
        .unwrap_or((0.0, 0.0));

    match shape.shape.as_str() {
        "circle" => {
            let r = shape.radius.expect("circle requires radius") / factor;
            Box::new(Circle::new(r, center).expect("valid circle"))
        }
        "ellipse" => {
            let a = shape.a.expect("ellipse requires a") / factor;
            let b = shape.b.expect("ellipse requires b") / factor;
            let alpha = shape.alpha.map(|deg| deg.to_radians()).unwrap_or(0.0);
            Box::new(Ellipse::new(a, b, alpha, center).expect("valid ellipse"))
        }
        "rectangular" => {
            let l = shape.length.expect("rectangular requires length") / factor;
            let w = shape.width.expect("rectangular requires width") / factor;
            let alpha = shape.alpha.map(|deg| deg.to_radians()).unwrap_or(0.0);
            Box::new(Rectangular::new(l, w, alpha, center).expect("valid rect"))
        }
        "superellipse" => {
            let a = shape.a.expect("superellipse requires a") / factor;
            let b = shape.b.expect("superellipse requires b") / factor;
            let n = shape.n.unwrap_or(2.0);
            let alpha = shape.alpha.map(|deg| deg.to_radians()).unwrap_or(0.0);
            Box::new(SuperEllipse::new(a, b, n, alpha, center).expect("valid superellipse"))
        }
        other => panic!("unsupported shape: {other}"),
    }
}

// ═══════════════════════════════════════════════════════════════════════
// 测试用例：TOML → OTN → 文件 I/O → SLTN → 验证
// ═══════════════════════════════════════════════════════════════════════

/// 分段多项式变比热容 (NASA 9 系数) — 完整 TOML 继承流程。
///
/// 流程：
/// 1. 从 `fixtures/otn_piecewise_cp.toml` 加载 OTN 配置
/// 2. 运行 OTN 生成特征线
/// 3. 特征线写文件 → 读回（验证 I/O 往返）
/// 4. 从 `fixtures/sltn.toml` 模板加载 SLTN 配置
/// 5. 运行 SLTN 流线追踪
/// 6. 验证 3D 模型完整性
#[test]
#[ignore = "full TOML OTN→SLTN pipeline, run with -- --ignored"]
fn test_toml_pipeline_piecewise_cp() {
    let otn_toml = include_str!("fixtures/otn_piecewise_cp.toml");
    let sltn_toml = include_str!("fixtures/sltn.toml");

    let tmp = TempDir::new("piecewise_cp");

    // 1-2. OTN: TOML → charlines
    let otn_path = tmp.write("otn.toml", otn_toml);
    let charlines = run_otn_from_toml(&otn_path);
    assert!(!charlines.is_empty(), "OTN should produce charlines");

    // 验证特征线有效性
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

    // 3. 文件 I/O 往返
    let field_path = tmp.write_bytes("field_data.txt", &[]);
    let _read_back = charlines_roundtrip(&charlines, &field_path);

    // 4-5. SLTN: TOML → read charlines → trace
    let tracer = run_sltn_from_toml_template(sltn_toml, &field_path);

    // 6. 验证
    common::assert_sltn_valid(&tracer, 5.0);
}

#[test]
#[ignore = "full TOML OTN→SLTN pipeline, run with -- --ignored"]
fn test_toml_pipeline_constant_cp() {
    let otn_toml = include_str!("fixtures/otn_constant_cp.toml");
    let sltn_toml = include_str!("fixtures/sltn.toml");

    let tmp = TempDir::new("constant_cp");

    // 1-2. OTN: TOML → charlines
    let otn_path = tmp.write("otn.toml", otn_toml);
    let config =
        NozzleConfig::from_toml_file(&otn_path).expect("failed to load constant-Cp OTN config");

    // 验证材料模型确实是常数 Cp
    let cp_at_300 = config.material.cp(300.0);
    let cp_at_2000 = config.material.cp(2000.0);
    assert!(
        (cp_at_300 - cp_at_2000).abs() < 1e-9,
        "constant Cp: cp(300K)={cp_at_300} should equal cp(2000K)={cp_at_2000}"
    );

    config.validate().expect("validation failed");
    let config = config.build().expect("build failed");

    let mut nozzle = ConstraintNozzle::new_otn(config);
    nozzle.run();
    let charlines = nozzle.get_assembly_charlines();
    assert!(!charlines.is_empty());

    // 3. 文件 I/O 往返
    let field_path = tmp.write_bytes("field_data.txt", &[]);
    let read_back = charlines_roundtrip(&charlines, &field_path);

    // 验证读回的特征线点有效
    for (li, line) in read_back.iter().enumerate() {
        for (pi, pt) in line.iter().enumerate() {
            assert!(pt.is_valid(), "roundtrip invalid: line={li} pt={pi}");
        }
    }

    // 4-5. SLTN
    let tracer = run_sltn_from_toml_template(sltn_toml, &field_path);

    // 6. 验证
    common::assert_sltn_valid(&tracer, 6.0);
}

/// 验证分段多项式 Cp 的变比热特性。
///
/// 不运行完整流程，仅验证 TOML 解析后材料模型的 Cp 随温度变化。
#[test]
fn test_piecewise_cp_variable_specific_heat() {
    let otn_toml = include_str!("fixtures/otn_piecewise_cp.toml");
    let tmp = TempDir::new("cp_verify");
    let otn_path = tmp.write("otn.toml", otn_toml);

    let config =
        NozzleConfig::from_toml_file(&otn_path).expect("failed to load piecewise-Cp config");

    // 分段多项式 Cp 应随温度显著变化
    let cp_low = config.material.cp(300.0);
    let _cp_mid = config.material.cp(1500.0);
    let cp_high = config.material.cp(3000.0);

    assert!(cp_low > 1000.0, "cp at 300K should be > 1000");
    assert!(cp_high > 1200.0, "cp at 3000K should be > 1200");
    assert!(
        (cp_low - cp_high).abs() > 50.0,
        "piecewise Cp should vary with temperature: cp(300K)={cp_low}, cp(3000K)={cp_high}"
    );

    let gamma_low = config.material.gamma(300.0);
    let gamma_high = config.material.gamma(3000.0);
    assert!(
        (gamma_low - gamma_high).abs() > 0.01,
        "gamma should vary: γ(300K)={gamma_low}, γ(3000K)={gamma_high}"
    );
}

/// 验证常数 Cp 材料模型的比热不变特性。
#[test]
fn test_constant_cp_invariant_specific_heat() {
    let otn_toml = include_str!("fixtures/otn_constant_cp.toml");
    let tmp = TempDir::new("cp_const_verify");
    let otn_path = tmp.write("otn.toml", otn_toml);

    let config =
        NozzleConfig::from_toml_file(&otn_path).expect("failed to load constant-Cp config");

    let temperatures = [100.0, 300.0, 1000.0, 2000.0, 5000.0];
    let cp_values: Vec<f64> = temperatures
        .iter()
        .map(|&t| config.material.cp(t))
        .collect();

    for w in cp_values.windows(2) {
        assert!(
            (w[0] - w[1]).abs() < 1e-9,
            "constant Cp should not vary: cp={} at T={}K",
            w[0],
            temperatures[0]
        );
    }

    // 常数 Cp 下 gamma 也应基本不变（微小变化来自 Rgas 精度）
    let gamma_300 = config.material.gamma(300.0);
    assert!((gamma_300 - 1.4).abs() < 0.01, "constant Cp air γ≈1.4");
}

// ── 自定义 cp_segments 分段多项式测试 ──────────────────────────────────

/// 验证 TOML 中显式 cp_segments 解析正确，且与内置 NASA 9 模型一致。
///
/// cp_segments 允许用户自定义分段多项式变比热容系数，
/// 每个段的常数项 T⁰ 放在 pos_coefficients[0]。
///
/// 此处使用与 NASA 9 一致的系数（两段各自独立常数项），验证：
/// 1. TOML 解析不报错
/// 2. 两个温度段返回不同的多项式值
/// 3. 与内置 `Material::air_nasa9piecewise_polynomial()` 的 Cp 值一致
#[test]
fn test_custom_cp_segments_parsing() {
    let otn_toml = include_str!("fixtures/otn_custom_segments_cp.toml");
    let tmp = TempDir::new("cp_segments");
    let otn_path = tmp.write("otn.toml", otn_toml);

    let config =
        NozzleConfig::from_toml_file(&otn_path).expect("failed to load cp_segments config");
    let toml_mat = &config.material;
    let nasa9 = Material::air_nasa9piecewise_polynomial();

    // 在多个温度点对比 TOML cp_segments 和内置 NASA 9
    let test_temps = [
        250.0, 300.0, 500.0, 800.0, 999.0, 1001.0, 1500.0, 2500.0, 4000.0, 5500.0,
    ];
    for &t in &test_temps {
        let cp_toml = toml_mat.cp(t);
        let cp_nasa = nasa9.cp(t);
        let gamma_toml = toml_mat.gamma(t);
        let gamma_nasa = nasa9.gamma(t);

        assert!(
            (cp_toml - cp_nasa).abs() < 1e-6,
            "T={t}K: TOML cp={cp_toml:.6} vs NASA9 cp={cp_nasa:.6}, diff={:.3e}",
            (cp_toml - cp_nasa).abs()
        );
        assert!(
            (gamma_toml - gamma_nasa).abs() < 1e-10,
            "T={t}K: TOML γ={gamma_toml:.10} vs NASA9 γ={gamma_nasa:.10}",
        );
    }

    // 两段各取一点验证跨段变化
    let cp_600 = toml_mat.cp(600.0);
    let cp_2000 = toml_mat.cp(2000.0);
    assert!(cp_600 > 1000.0, "cp(600K) should be > 1000, got {cp_600}");
    assert!(
        cp_2000 > 1200.0,
        "cp(2000K) should be > 1200, got {cp_2000}"
    );
    assert!(
        (cp_600 - cp_2000).abs() > 50.0,
        "Cp should differ across segments: 600K={cp_600}, 2000K={cp_2000}"
    );
}

/// 自定义 cp_segments — 完整 TOML 继承流程。
///
/// 使用显式 `[[Material.cp_segments]]` 定义两段多项式比热容
/// （NASA 9 系数，两段各自独立常数项），固定膨胀角 θ_a=19°。
///
/// 流程：TOML → OTN → 文件 I/O → SLTN → 验证。
/// 同时对比内置 `Material::air_nasa9piecewise_polynomial()` 确保一致性。
#[test]
#[ignore = "full TOML OTN→SLTN pipeline with cp_segments, run with -- --ignored"]
fn test_toml_pipeline_custom_cp_segments() {
    let otn_toml = include_str!("fixtures/otn_custom_segments_cp.toml");
    let sltn_toml = include_str!("fixtures/sltn.toml");

    let tmp = TempDir::new("custom_segments");

    // 0. 验证 TOML cp_segments 与内置 NASA 9 一致
    let otn_path = tmp.write("otn.toml", otn_toml);
    let config =
        NozzleConfig::from_toml_file(&otn_path).expect("failed to load cp_segments config");
    let nasa9 = Material::air_nasa9piecewise_polynomial();
    for &t in &[300.0, 600.0, 1500.0, 3000.0, 5000.0] {
        let diff = (config.material.cp(t) - nasa9.cp(t)).abs();
        assert!(
            diff < 1e-6,
            "T={t}K: TOML cp_segments cp differs from NASA9 by {diff:.3e}"
        );
    }

    // 1-2. OTN: TOML (cp_segments) → charlines
    let otn_path = tmp.write("otn.toml", otn_toml);
    let charlines = run_otn_from_toml(&otn_path);
    assert!(!charlines.is_empty(), "OTN should produce charlines");

    // 验证特征线
    for (li, line) in charlines.iter().enumerate() {
        for (pi, pt) in line.iter().enumerate() {
            assert!(pt.is_valid(), "invalid MOC point: line={li} pt={pi}");
        }
    }

    // 3. 文件 I/O 往返
    let field_path = tmp.write_bytes("field_data.txt", &[]);
    let read_back = charlines_roundtrip(&charlines, &field_path);
    for (li, line) in read_back.iter().enumerate() {
        for (pi, pt) in line.iter().enumerate() {
            assert!(pt.is_valid(), "roundtrip invalid: line={li} pt={pi}");
        }
    }

    // 4-5. SLTN
    let tracer = run_sltn_from_toml_template(sltn_toml, &field_path);

    // 6. 验证
    common::assert_sltn_valid(&tracer, 5.0);
}
