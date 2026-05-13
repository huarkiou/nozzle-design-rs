# GEOMETRY — CRATE KNOWLEDGE BASE

## OVERVIEW
三维几何形状定义与处理 — ClosedCurve trait, 标准截面形状 (Circle/Ellipse/SuperEllipse/Rectangular/UserDefined), Point3d/WallPoints, OBJ 模型导出, wallpoint 插值/重采样/合并. Depends on math, nalgebra 0.34, optional serde (default enabled).

## WHERE TO LOOK
| Task | File | Notes |
|------|------|-------|
| ClosedCurve trait | `src/closed_curve.rs` | 抽象接口: `sample_by_angle(n) → Vec<Point3d>`, `area()`, `centroid()` |
| Circle shape | `src/circle.rs` | 实现 ClosedCurve: radius, center |
| Ellipse shape | `src/ellipse.rs` | 实现 ClosedCurve: semi-axes a/b, center, rotation alpha |
| SuperEllipse shape | `src/superellipse.rs` | 实现 ClosedCurve: a/b 半轴 + n 指数控制方-圆过渡 |
| Rectangular shape | `src/rectangular.rs` | 实现 ClosedCurve: width/height, corner radius, center |
| UserDefined shape | `src/userdefined.rs` | 实现 ClosedCurve: 自定义离散点, 插值采样 |
| Point3d / WallPoints | `src/point.rs` | 三维点类型 + 壁面点集 (serde-gated), Display 格式化 |
| WallPoint 处理 | `src/wallpoints.rs` | 672 行 — 插值, 重采样, 合并, 文件 I/O (DAT/OBJ) |
| OBJ 导出 | `src/obj.rs` | ObjModel: vertices, quads, export-to-file |
| 基础几何工具 | `src/basics.rs` | 射线-线段求交, 符号函数 |
| Crate root | `src/lib.rs` | 模块声明 + 便利重导出 (ClosedCurve, Point3d, WallPoints, all shapes) |

## CONVENTIONS
- `serde` 是 optional feature, 默认启用 (`default = ["serde"]`)
- `#[cfg(feature = "serde")]` 门控用于 `Point3d` 和 `WallPoints` 的 Serialize/Deserialize derive (in `point.rs`)
- 无 crate 级集成测试目录 (tests/); 单元测试内联 (`#[cfg(test)] mod tests`)
- ClosedCurve trait 的 `sample_by_angle` 返回的 Point3d 坐标均位于 z=0 平面

## ANTI-PATTERNS
- DO NOT 添加新闭合曲线形状而不实现完整的 `ClosedCurve` trait
- DO NOT 在 `ClosedCurve::sample_by_angle` 中假设点数大于等于 3 (调用方负责传入合理值)
- DO NOT 使用 `unsafe` 代码 — 项目策略是纯安全 Rust
