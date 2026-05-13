# PROJECT KNOWLEDGE BASE

**Generated:** 2026-05-13
**Commit:** 381c03f
**Branch:** main

## OVERVIEW
超声速喷管扩张段设计工具集 (Supersonic nozzle design toolkit). Rust workspace using Method of Characteristics (MOC). 3 lib crates (math → geometry → aero) + 2 CLI apps (otn, sltn).

## STRUCTURE
```
nozzle-design-rs/
├── crates/
│   ├── math/          # 数值计算基础库 (root-finding, quadrature, polynomials, 2D coords)
│   ├── geometry/      # 三维几何形状定义与处理 (ClosedCurve, OBJ export, wallpoints)
│   └── aero/          # 气动热力学核心 (MOC solver, nozzle sections, streamline trace)
├── apps/
│   ├── otn/           # 最优推力喷管 CLI (Optimal Thrust Nozzle)
│   └── sltn/          # 三维流线追踪喷管 CLI (StreamLine Tracing Nozzle)
├── tests/             # 项目级集成测试 (全部 #[ignore], 运行需 -- --ignored)
├── tools/             # Python 可视化脚本 (matplotlib, numpy, scipy)
└── .github/workflows/ # CI/CD (BuildTest + Release)
```

## WHERE TO LOOK
| Task | Location | Notes |
|------|----------|-------|
| OTN design logic | `crates/aero/src/nozzle/` | initial_line, expansion_section, transition_section, etc. |
| MOC unit processes | `crates/aero/src/moc/unitprocess/` | Irrotational (1399 lines, primary) + Rotational (stub) |
| Root finding algorithms | `crates/math/src/rootfinding/` | bisection, secant, TOMS748, Newton-2D |
| Gauss-Legendre quadrature | `crates/math/src/quadrature/gauss_legendre.rs` | 844 lines |
| Shape definitions | `crates/geometry/src/{circle,ellipse,rectangular,superellipse}.rs` | ClosedCurve trait |
| Wallpoint processing | `crates/geometry/src/wallpoints.rs` | interpolation, resampling, merge, I/O |
| Streamline tracing | `crates/aero/src/streamline_trace/` | trace, intersect, merge, weight |
| Material / cp models | `crates/aero/src/{material,cp}.rs` | NASA 9-coefficient, constant cp, piecewise polynomial |
| Isentropic relations | `crates/aero/src/isentropic.rs` | T_total ↔ T_static |
| Config parsing (OTN) | `crates/aero/src/nozzle/config.rs` | serde + toml deserialization |
| Config parsing (SLTN) | `apps/sltn/src/config.rs` | serde + toml deserialization |
| Integration tests | `tests/` | common/ helpers, fixtures/, 4 test files |
| Benchmarks | `crates/aero/benches/` | criterion, harness=false, requires nightly |
| CI | `.github/workflows/buildtest.yml` | fmt → build → test |

## CONVENTIONS
- **Edition 2024** (requires Rust ≥1.85)
- **Default rustfmt** — no `rustfmt.toml`; CI enforces via `cargo fmt --check`
- **No clippy config** — no `clippy.toml`; default rules apply
- **Chinese documentation/comments** throughout the codebase
- **TOML-based configuration** for both CLI apps (serde + toml 0.9)
- **nalgebra 0.34** for linear algebra, **rayon 1.10** for parallelism
- **No `unsafe` blocks** — entire project is safe Rust
- **Release profile**: `opt-level = 3`, `lto = "thin"`
- **Workspace resolver**: `"3"` (tied to edition 2024)

## ANTI-PATTERNS (THIS PROJECT)
- **DO NOT** add `unsafe` blocks — project policy is safe Rust only
- **DO NOT** modify `crates/aero/src/moc/unitprocess/rotational.rs` without implementing the full `UnitProcess` trait — it is an intentional stub (8x `todo!()`)
- **DO NOT** suppress warnings with undocumented `#[allow(...)]` — use `#[expect(...)]` (Rust 1.85+) when intentional
- **DO NOT** add undocumented dependencies without updating all `Cargo.toml` files

## UNIQUE STYLES
- Module file + sibling directory pattern: `moc.rs` + `moc/` for sub-modules
- Integration tests use section-comment banners: `// ═══════ Flow Type ═══════`
- Test helpers centralized in `tests/common/mod.rs`
- TOML fixture files embedded via `include_str!()` at compile time
- `#[ignore]` on all integration tests (computationally expensive)
- Custom `TempDir` RAII in `tests/otn_sltn_toml.rs` for test isolation

## COMMANDS
```bash
# Build
cargo build                        # Debug build all workspace members
cargo build --release              # Release (thin LTO, opt-level=3)
cargo build -p otn --release       # Build single app

# Test
cargo test --workspace             # Unit + doc tests (fast)
cargo test -- --ignored            # Integration tests (slow)
cargo test --test otn_parameters -- --ignored  # Specific test file

# Format
cargo fmt --check                  # Required by CI before commit

# Benchmark (nightly toolchain required)
cargo bench -p aero

# Python tools
uv sync                            # Install Python dependencies
```

## NOTES
- **Rotational MOC is a STUB**: `rotational.rs` has all 8 trait methods as `todo!()`. Listed as "未实现" in README. The `#[allow(unused)]` on the impl suppresses the dead-code warning.
- **Root `Cargo.toml` ghost package**: Has `[package]` section with no `src/` — serves as workspace aggregator for dev-dependencies and release profile.
- **No `examples/` directory exists** — API usage must be learned from integration tests or app source.
- **No `rustfmt.toml`, `clippy.toml`, or `rust-toolchain.toml`** — all Rust defaults apply.
- **Python tools require ≥3.13** — managed via `uv`, separate from Rust build.
