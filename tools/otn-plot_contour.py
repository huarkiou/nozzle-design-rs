"""OTN 喷管特征线法流场可视化。

对每个 fluid_field*.txt 文件生成两幅图:
1. 特征线散点图 — 渐变色区分线序, 叠加壁面线和出口线
2. 流场参数云图 — Delaunay 三角剖分 + 壁面/出口叠加

用法:
    python tools/otn-plot_contour.py                      # 默认所有 fluid_field*.txt
    python tools/otn-plot_contour.py file1.txt file2.txt  # 指定文件
    python tools/otn-plot_contour.py --field p            # 压力云图
    python tools/otn-plot_contour.py --no-show            # 不弹窗

字段: ma(马赫数,默认), p(压力), t(温度), rho(密度), v(速度)
"""

import argparse
import os
import pathlib

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.tri import Triangulation

plt.style.use("default")
plt.rcParams.update(
    {
        "font.size": 11,
        "axes.labelsize": 13,
        "axes.titlesize": 13,
        "legend.fontsize": 9,
        "figure.figsize": (11, 6),
        "figure.dpi": 120,
    }
)

# 文件列: x, y, V, theta, p, rho, T, Rg, gamma, Ma
COL_X, COL_Y = 0, 1
COL_V, COL_THETA = 2, 3
COL_P, COL_RHO = 4, 5
COL_T, COL_RG = 6, 7
COL_GAMMA, COL_MA = 8, 9

FIELD_META = {
    "ma": (COL_MA, "Mach Number", "plasma"),
    "p": (COL_P, "Pressure (Pa)", "coolwarm"),
    "t": (COL_T, "Temperature (K)", "hot"),
    "rho": (COL_RHO, "Density (kg/m^3)", "viridis"),
    "v": (COL_V, "Velocity (m/s)", "turbo"),
}


# ── 数据加载 ──────────────────────────────────────────────────


def load_fluid_data(filepath: str) -> list[np.ndarray]:
    """按空行分组读取流场数据, 返回每条特征线的 ndarray。"""
    if not os.path.isfile(filepath):
        raise FileNotFoundError(f"文件不存在: {filepath}")

    groups: list[np.ndarray] = []
    current: list[list[float]] = []

    with open(filepath, "r") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                if current:
                    arr = np.array(current)
                    arr = arr[np.isfinite(arr).all(axis=1)]
                    if len(arr) >= 1:
                        groups.append(arr)
                    current = []
            else:
                try:
                    current.append([float(v) for v in line.split(",")])
                except ValueError:
                    continue
        if current:
            arr = np.array(current)
            arr = arr[np.isfinite(arr).all(axis=1)]
            if len(arr) >= 1:
                groups.append(arr)

    if not groups:
        raise ValueError("未找到有效数据")
    return groups


# ── 绘制 ──────────────────────────────────────────────────────


def plot_charline_scatter(
    groups: list[np.ndarray], ax: plt.Axes, title: str = "Characteristic Lines"
):
    """特征线散点图: 渐变色 + 壁面/出口高亮。"""
    n = len(groups)
    colors = plt.get_cmap("viridis")(np.linspace(0.05, 0.95, max(n, 2)))

    for i, g in enumerate(groups):
        x, y = g[:, COL_X], g[:, COL_Y]
        if len(g) >= 2:
            ax.plot(x, y, "-", linewidth=0.3, color=colors[i], alpha=0.7)
        else:
            ax.plot(x, y, "o", markersize=2, color=colors[i], alpha=0.9)

    # 壁面线: 每条特征线第一点 (跳过 IVL 和 Exit)
    wall_pts = [g[0] for g in groups[1:-1]]
    wall = np.array(wall_pts)
    ax.plot(
        wall[:, COL_X],
        wall[:, COL_Y],
        "-",
        color="black",
        linewidth=1.2,
        label="Wall",
        zorder=10,
    )

    # 出口线: 最后一条 (所有点 x~=length)
    last = groups[-1]
    if np.all(np.abs(last[:, COL_X] - last[0, COL_X]) < 1e-3):
        ax.plot(
            last[:, COL_X],
            last[:, COL_Y],
            "-",
            color="red",
            linewidth=2.0,
            label="Exit Line",
            zorder=10,
        )

    ax.set_aspect("equal")
    ax.set_xlabel(r"$x$ (m)")
    ax.set_ylabel(r"$y$ (m)")
    ax.set_title(title)
    ax.legend(loc="upper right", framealpha=0.8)
    ax.grid(True, alpha=0.2)


def plot_contour_tri(
    x,
    y,
    values,
    ax,
    attr_name="Mach Number",
    cmap="plasma",
    n_levels=40,
    wall_line=None,
    exit_line=None,
):
    """Delaunay 三角剖分填充等值线图, 叠加壁面/出口。"""
    if len(x) < 3:
        ax.text(0.5, 0.5, "Too few points", transform=ax.transAxes, ha="center")
        return

    v_min, v_max = float(np.nanmin(values)), float(np.nanmax(values))
    if v_min == v_max:
        v_min -= max(abs(v_min) * 0.01, 1e-6)
        v_max += max(abs(v_max) * 0.01, 1e-6)
    levels = np.linspace(v_min, v_max, n_levels + 1)

    tri = _safe_triangulation(x, y)
    if tri is None:
        ax.text(0.5, 0.5, "Triangulation failed", transform=ax.transAxes, ha="center")
        return

    tcf = ax.tricontourf(tri, values, levels=levels, cmap=cmap, alpha=0.9)

    if wall_line is not None and len(wall_line) >= 2:
        ax.plot(
            wall_line[:, 0],
            wall_line[:, 1],
            "-",
            color="black",
            linewidth=1.2,
            zorder=10,
        )
    if exit_line is not None and len(exit_line) >= 2:
        ax.plot(
            exit_line[:, 0], exit_line[:, 1], "-", color="red", linewidth=1.5, zorder=10
        )

    ax.set_aspect("equal")
    ax.set_xlabel(r"$x$ (m)")
    ax.set_ylabel(r"$y$ (m)")
    ax.set_title(f"{attr_name} Contour")

    cbar = plt.colorbar(tcf, ax=ax, orientation="horizontal", pad=0.10, aspect=35)
    cbar.set_label(attr_name, fontsize=11)
    tick_step = max(1, n_levels // 8)
    cbar.set_ticks(list(levels[::tick_step]))


def _safe_triangulation(x, y):
    """安全的三角剖分, 退化时加微小扰动。"""
    try:
        return Triangulation(x, y)
    except Exception:
        pass
    try:
        rng = np.random.default_rng(42)
        return Triangulation(
            x + rng.uniform(-1e-8, 1e-8, len(x)), y + rng.uniform(-1e-8, 1e-8, len(y))
        )
    except Exception:
        return None


# ── 摘要 ──────────────────────────────────────────────────────


def print_summary(groups: list[np.ndarray]):
    """打印特征线统计。"""
    n = len(groups)
    ivl = groups[0]
    last = groups[-1]
    total_pts = sum(len(g) for g in groups)

    print(f"  {n} charlines, {total_pts} total points")
    print(
        f"  IVL: {len(ivl)} pts, "
        f"x={ivl[0, COL_X]:.3f}, y=[{ivl[-1, COL_Y]:.3f}..{ivl[0, COL_Y]:.3f}]"
    )

    if np.all(np.abs(last[:, COL_X] - last[0, COL_X]) < 1e-3):
        print(
            f"  Exit line: {len(last)} pts, "
            f"x={last[0, COL_X]:.3f}, "
            f"y=[{last[-1, COL_Y]:.3f}..{last[0, COL_Y]:.3f}]"
        )

    mid = groups[n // 2]
    print(
        f"  x: [{ivl[0, COL_X]:.3f}..{last[-1, COL_X]:.3f}], "
        f"y: [{last[-1, COL_Y]:.3f}..{mid[0, COL_Y]:.3f}]"
    )


# ── 单文件处理 ────────────────────────────────────────────────


def process_one_file(
    data_path: pathlib.Path, out_dir: pathlib.Path, field: str = "ma", show: bool = True
):
    stem = data_path.stem
    file_out_dir = out_dir / stem
    file_out_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n{'=' * 60}")
    print(f"读取: {data_path}")
    try:
        groups = load_fluid_data(str(data_path))
    except FileNotFoundError:
        print("  skip: file not found")
        return
    except ValueError as e:
        print(f"  skip: {e}")
        return

    print_summary(groups)

    all_pts = np.vstack(groups)
    x_all, y_all = all_pts[:, COL_X], all_pts[:, COL_Y]

    # 壁面线 / 出口线
    wall_line = np.array([g[0] for g in groups[1:-1]])
    last = groups[-1]
    is_exit = np.all(np.abs(last[:, COL_X] - last[0, COL_X]) < 1e-3)
    exit_line = last if is_exit else None

    # 场参数
    col_idx, attr_name, cmap = FIELD_META.get(field, FIELD_META["ma"])
    field_vals = all_pts[:, col_idx]
    print(f"  {attr_name}: [{field_vals.min():.4g}, {field_vals.max():.4g}]")

    # 图1: 特征线
    fig1, ax1 = plt.subplots()
    plot_charline_scatter(groups, ax1, title=f"MOC — {stem}")
    fig1.tight_layout()
    p1 = file_out_dir / f"{stem}_charline_scatter.png"
    fig1.savefig(p1, dpi=200, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig1)
    print(f"  scatter → {p1}")

    # 图2: 云图
    fig2, ax2 = plt.subplots()
    plot_contour_tri(
        x_all,
        y_all,
        field_vals,
        ax2,
        attr_name=attr_name,
        cmap=cmap,
        wall_line=wall_line,
        exit_line=exit_line,
    )
    fig2.tight_layout()
    suffix = {
        "ma": "ma",
        "p": "pressure",
        "t": "temperature",
        "rho": "density",
        "v": "velocity",
    }.get(field, field)
    p2 = file_out_dir / f"{stem}_{suffix}_contour.png"
    fig2.savefig(p2, dpi=200, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig2)
    print(f"  contour → {p2}")


# ── 主流程 ────────────────────────────────────────────────────


def main():
    parser = argparse.ArgumentParser(description="OTN 喷管流场可视化")
    parser.add_argument("files", nargs="*", help="fluid_field*.txt 文件")
    parser.add_argument(
        "--field",
        "-f",
        default="ma",
        choices=["ma", "p", "t", "rho", "v"],
        help="云图字段 (default: ma)",
    )
    parser.add_argument("--no-show", action="store_true", help="不弹出 matplotlib 窗口")
    parser.add_argument(
        "--out", "-o", default=None, help="输出目录 (default: target/tmp/output)"
    )
    args = parser.parse_args()

    proj_root = pathlib.Path(__file__).resolve().parent.parent
    tmp_dir = proj_root / "target" / "tmp"
    out_dir = pathlib.Path(args.out) if args.out else (tmp_dir / "output")

    if args.files:
        file_paths = [pathlib.Path(p).resolve() for p in args.files]
    else:
        file_paths = sorted(tmp_dir.glob("fluid_field*.txt"))

    if not file_paths:
        print("未找到 fluid_field*.txt。请先运行:")
        print("  cargo test -p aero -- test_new_and_run test_with_expansion")
        return

    print(f"处理 {len(file_paths)} 文件, field={args.field}:")
    for fp in file_paths:
        print(f"  - {fp}")

    for fp in file_paths:
        process_one_file(fp, out_dir, field=args.field, show=not args.no_show)

    print(f"\n完成 → {out_dir}")


if __name__ == "__main__":
    main()
