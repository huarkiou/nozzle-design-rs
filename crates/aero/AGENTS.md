# CRATE KNOWLEDGE BASE — aero
**Generated:** 2026-05-13

## OVERVIEW
气动热力学核心库 — MOC solver, nozzle sections, 3D streamline tracing, material models.
Largest crate: 27 .rs files, 3 submodules. Deps: math, geometry, nalgebra 0.34, rayon 1.10,
serde, toml 0.9. Dev-deps: criterion 0.7, serde_json 1.0.
## STRUCTURE
```
crates/aero/
├── src/
│   ├── lib.rs                 # 6 pub modules
│   ├── cp.rs                  # 比热容模型 (NASA 9 / constant / piecewise polynomial)
│   ├── isentropic.rs          # 等熵关系 (T_total ↔ T_static)
│   ├── material.rs            # 工质属性 (MW, cp, γ, R, enthalpy)
│   ├── moc.rs + moc/          # 特征线法核心
│   │   ├── unitprocess.rs+unitprocess/ # Irrotational (1399 ln) + Rotational (STUB) + Config
│   │   ├── mocpoint.rs        # 特征点 (MocPoint, 663 lines)
│   │   ├── charline.rs        # 单条特征线
│   │   ├── charlines.rs       # 特征线集合
│   │   └── areatype.rs        # 流动类型 (Axisymmetric / Planar)
│   ├── nozzle.rs + nozzle/    # 喷管各段
│   │   ├── config.rs          # OTN TOML 反序列化 (446 lines)
│   │   ├── constraint_nozzle.rs # 约束喷管总装 (918 lines)
│   │   ├── expansion_section.rs # 膨胀段 (theta_a)
│   │   ├── initial_line.rs    # 喉部初值线
│   │   ├── initial_section.rs # 初值问题
│   │   ├── section.rs         # Section trait
│   │   ├── transition_section.rs # 转向段 (二分搜索出口长度)
│   │   └── uniform_section.rs # 均一区
│   └── streamline_trace/      # 三维流线追踪
│       ├── mod.rs
│       ├── trace.rs           # 向前/向后追踪
│       ├── intersect.rs       # 流线-特征线交点
│       ├── merge.rs           # 加权过渡融合
│       └── weight.rs          # 加权函数
├── benches/                   # Criterion (4 files, harness=false)
├── tests/                     # Crate 级集成测试 (2 files)
└── Cargo.toml
```
## WHERE TO LOOK
| Task | File | Notes |
|------|------|-------|
| Irrotational MOC solver | `src/moc/unitprocess/irrotational.rs` | 1399 ln, primary production solver |
| MOC unit process config | `src/moc/unitprocess/config.rs` | Grid resolution, tolerance |
| MOC point type | `src/moc/mocpoint.rs` | 663 ln, core data structure |
| Nozzle initial line | `src/nozzle/initial_line.rs` | Starting characteristic from throat |
| Expansion section | `src/nozzle/expansion_section.rs` | Controlled by theta_a |
| Transition section | `src/nozzle/transition_section.rs` | Binary search for exit length |
| Constraint nozzle | `src/nozzle/constraint_nozzle.rs` | 918 ln, full nozzle assembly |
| OTN config | `src/nozzle/config.rs` | TOML deserialization + validation |
| Material model | `src/material.rs` | MW, cp, γ, R, enthalpy |
| Cp models | `src/cp.rs` | NASA 9, constant, piecewise polynomial |
| Isentropic relations | `src/isentropic.rs` | T_total/T_static conversions |
| Streamline trace | `src/streamline_trace/trace.rs` | Forward/backward along charlines |
| Wall point merge | `src/streamline_trace/merge.rs` | Weighted transition blending |
| Intersection solver | `src/streamline_trace/intersect.rs` | Ray-charline intersection |
## CONVENTIONS
- `#[inline(always)]` only in hot-path MOC solver methods; benchmark before/after
- `rayon` parallel iterators in MOC grid computation and streamline trace
- Criterion `harness = false`, nightly toolchain: `cargo bench -p aero`
- `<name>.rs` + `<name>/` module pattern for `moc`, `nozzle`

## ANTI-PATTERNS (THIS CRATE)
- **CRITICAL**: `rotational.rs` is a COMPLETE STUB. All 8 `UnitProcess` trait methods are
  `todo!()`. `#[allow(unused)]` suppresses dead-code warnings. DO NOT depend on `Rotational`.
- DO NOT refactor `constraint_nozzle.rs` without understanding the full assembly pipeline:
  InitialLine → InitialSection → ExpansionSection → TransitionSection → UniformSection
- DO NOT skip benchmarks when modifying MOC hot paths: `mocpoint.rs`, `irrotational.rs`,
  `charlines.rs` — profile with `cargo bench -p aero`
