# CRATE KNOWLEDGE BASE — math

## OVERVIEW
数值计算基础库 — root-finding, quadrature, polynomials, 2D geometry primitives. Depends on nalgebra 0.34, num-traits 0.2.

## STRUCTURE
```
crates/math/src/
├── lib.rs              # re-exports five items (geometry, quadrature, rootfinding, Tolerance, polynomial)
├── tolerance.rs        # 容差抽象 (Tolerance trait)
├── rootfinding.rs + rootfinding/
│   ├── bisection.rs    # 二分法
│   ├── secant.rs       # 割线法
│   ├── toms748.rs      # TOMS 748 有界求根
│   └── newton2d.rs     # 牛顿二维求解器
├── quadrature.rs + quadrature/
│   └── gauss_legendre.rs    # 高斯-勒让德数值积分 (844 lines)
├── polynomial.rs + polynomial/
│   ├── base.rs              # 多项式基础类型
│   └── piecewisepolynomial.rs  # 分段多项式
└── geometry.rs + geometry/
    ├── coord2d.rs        # 二维坐标 (Coord2d)
    └── polyline.rs       # 二维折线 (Polyline2d)
```

## WHERE TO LOOK
| Task | File |
|------|------|
| Bisection root finding | `src/rootfinding/bisection.rs` |
| Secant root finding | `src/rootfinding/secant.rs` |
| TOMS748 root finding | `src/rootfinding/toms748.rs` |
| Newton-2D solver | `src/rootfinding/newton2d.rs` |
| Rootfinding error types | `src/rootfinding/error.rs` |
| Gauss-Legendre quadrature | `src/quadrature/gauss_legendre.rs` |
| Polynomial base types | `src/polynomial/base.rs` |
| Piecewise polynomial | `src/polynomial/piecewisepolynomial.rs` |
| 2D coordinate type | `src/geometry/coord2d.rs` |
| 2D polyline operations | `src/geometry/polyline.rs` |
| Tolerance abstraction | `src/tolerance.rs` |

## CONVENTIONS
- No crate-level integration tests (tested transitively via aero crate and workspace tests)
- All unit tests inline (`#[cfg(test)] mod tests`)
- `tolerance` module is private; `Tolerance` trait re-exported at crate root
- Uses `num-traits::Float` bounds for generic numeric parameters
- Uses `nalgebra::SVector<F, 2>` for 2D vector types

## ANTI-PATTERNS
- DO NOT add dependencies beyond num-traits and nalgebra (keep foundation lean)
- DO NOT add integration tests here — they belong in workspace `tests/`
