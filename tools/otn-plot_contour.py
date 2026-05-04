"""OTN 喷管特征线法流场可视化。

对每个 fluid_field*.txt 文件生成两幅图:
1. 特征线散点图 — 每条右行特征线用不同颜色连线+标记点,
   用于观察特征线走向、间距, 判断算法逻辑是否存在问题。
2. 马赫数云图 — 基于 Delaunay 三角剖分的填充等值线,
   用于观察结果数值是否正确。

用法:
    # 绘制 target/tmp/ 下所有 fluid_field*.txt 文件
    python tools/otn-plot_contour.py

    # 指定单个文件
    python tools/otn-plot_contour.py target/tmp/fluid_field.txt

    # 指定多个文件
    python tools/otn-plot_contour.py target/tmp/fluid_field.txt target/tmp/fluid_field_expansion.txt

    (需先从项目根目录执行 `cargo test -p aero -- test_new_and_run test_with_expansion`
     生成 target/tmp/fluid_field*.txt)
"""

import os
import pathlib
import sys

import matplotlib.pyplot as plt
import numpy as np
import scienceplots  # noqa: F401
from matplotlib.tri import Triangulation

plt.style.use(["science", "no-latex"])
plt.rcParams.update(
    {
        "font.size": 12,
        "axes.labelsize": 14,
        "axes.titlesize": 14,
        "legend.fontsize": 10,
        "figure.figsize": (10, 6),
    }
)


# ── 数据加载 ──────────────────────────────────────────────────


def load_fluid_data_grouped(filepath: str) -> list[np.ndarray]:
    """按空行分组读取流场数据, 保留每条特征线的独立结构。

    fluid_field*.txt 格式: 逗号分隔, 10 列, 空行分隔不同的特征线组。
    列: x, y, V, theta, p, rho, T, Rg, gamma, Ma

    Returns:
        groups: 每条特征线为一个 (n_points × 10) 的 ndarray。
                第一条为初值线 (IVL), 后续为右行特征线。
    """
    if not os.path.isfile(filepath):
        raise FileNotFoundError(f"数据文件不存在: {filepath}")

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
        raise ValueError("未找到有效数据组 (所有行均为空或含 NaN)")
    return groups


# ── 绘制 ──────────────────────────────────────────────────────


def plot_charline_scatter(
    groups: list[np.ndarray],
    ax: plt.Axes,
    title: str = "Characteristic Lines",
):
    """逐条特征线绘制连线+散点, 用渐变色区分不同线。"""
    n = len(groups)
    colors = plt.cm.plasma(np.linspace(0.05, 0.95, n))  # ty:ignore[unresolved-attribute]

    for i, line in enumerate(groups):
        x, y = line[:, 0], line[:, 1]
        if len(line) >= 2:
            ax.plot(
                x, y, "-o",
                markersize=2.0 if i == 0 else 1.2,
                linewidth=0.6 if i == 0 else 0.4,
                color=colors[i],
                alpha=0.8,
            )
        else:
            ax.plot(
                x, y, "o",
                markersize=3.0, color=colors[i], alpha=0.9,
            )

    if n > 0:
        ivl = groups[0]
        ax.plot(
            ivl[:, 0],
            ivl[:, 1],
            "-o",
            markersize=2.5,
            linewidth=1.0,
            color="black",
            alpha=0.9,
            label="Initial Value Line (IVL)",
        )

    ax.set_aspect("equal")
    ax.set_xlabel(r"$x$ (m)")
    ax.set_ylabel(r"$y$ (m)")
    ax.set_title(title)
    ax.legend(loc="upper right")
    ax.grid(True, alpha=0.3)


def plot_contour_tri(
    x: np.ndarray,
    y: np.ndarray,
    values: np.ndarray,
    ax: plt.Axes,
    attr_name: str = "Mach Number",
    n_levels: int = 30,
    cmap: str = "plasma",
):
    """基于 Delaunay 三角剖分的填充等值线图。"""
    if len(x) < 3:
        raise ValueError("点太少, 无法三角剖分")

    v_min, v_max = float(np.nanmin(values)), float(np.nanmax(values))
    if v_min == v_max:
        v_min -= max(abs(v_min) * 0.01, 1e-6)
        v_max += max(abs(v_max) * 0.01, 1e-6)
    levels = np.linspace(v_min, v_max, n_levels + 1)

    tri = Triangulation(x, y)
    tcf = ax.tricontourf(tri, values, levels=levels, cmap=cmap, alpha=0.9)

    ax.set_aspect("equal")
    ax.set_xlabel(r"$x$ (m)")
    ax.set_ylabel(r"$y$ (m)")
    ax.set_title(f"{attr_name} Contour (Delaunay Triangulation)")

    cbar = plt.colorbar(tcf, ax=ax, orientation="horizontal", pad=0.12, aspect=30)
    cbar.set_label(attr_name, fontsize=12)
    tick_step = max(1, n_levels // 6)
    cbar.set_ticks(levels[::tick_step])  # ty:ignore[invalid-argument-type]


def print_group_summary(groups: list[np.ndarray], label_width: int = 20):
    """打印每条特征线的统计信息。"""
    for i, g in enumerate(groups):
        lbl = "IVL (初值线)" if i == 0 else f"特征线 #{i}"
        print(
            f"  {lbl:<{label_width}} {len(g):>5} 个点, "
            f"x=[{g[:, 0].min():.4f}, {g[:, 0].max():.4f}], "
            f"y=[{g[:, 1].min():.4f}, {g[:, 1].max():.4f}]"
        )


def process_one_file(data_path: pathlib.Path, out_dir: pathlib.Path):
    """对单个数据文件生成散点图+云图。

    Args:
        data_path: fluid_field*.txt 的路径。
        out_dir: 输出目录的根 (会创建以文件名命名的子目录)。
    """
    stem = data_path.stem  # e.g. "fluid_field" or "fluid_field_expansion"
    file_out_dir = out_dir / stem
    file_out_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n{'=' * 60}")
    print(f"读取数据: {data_path}")
    try:
        groups = load_fluid_data_grouped(str(data_path))
    except FileNotFoundError:
        print("  ⚠ 文件不存在, 跳过")
        return
    except ValueError as e:
        print(f"  ⚠ 数据错误: {e}")
        return

    print(f"共 {len(groups)} 条特征线")
    print_group_summary(groups, label_width=15)

    all_pts = np.vstack(groups)
    x_all, y_all = all_pts[:, 0], all_pts[:, 1]
    Ma_all = all_pts[:, 9]

    print(f"总点数: {len(all_pts)}, Ma 范围: [{Ma_all.min():.4f}, {Ma_all.max():.4f}]")

    # ── 图1: 特征线散点图 ──
    fig1, ax1 = plt.subplots()
    plot_charline_scatter(groups, ax1, title=f"MOC Characteristic Lines — {stem}")
    fig1.tight_layout()
    scatter_path = file_out_dir / f"{stem}_charline_scatter.png"
    fig1.savefig(scatter_path, dpi=200, bbox_inches="tight")
    plt.show()
    plt.close(fig1)
    print(f"  散点图 → {scatter_path}")

    # ── 图2: 马赫数云图 ──
    fig2, ax2 = plt.subplots()
    plot_contour_tri(x_all, y_all, Ma_all, ax2, attr_name="Mach Number", cmap="plasma")
    fig2.tight_layout()
    contour_path = file_out_dir / f"{stem}_ma_contour.png"
    fig2.savefig(contour_path, dpi=200, bbox_inches="tight")
    plt.close(fig2)
    print(f"  云图   → {contour_path}")


# ── 主流程 ────────────────────────────────────────────────────


def main():
    proj_root = pathlib.Path(__file__).resolve().parent.parent
    tmp_dir = proj_root / "target" / "tmp"
    out_dir = tmp_dir / "output"

    # 收集要处理的文件
    if len(sys.argv) > 1:
        # 命令行指定文件
        file_paths = [pathlib.Path(p).resolve() for p in sys.argv[1:]]
    else:
        # 默认: 所有 fluid_field*.txt
        file_paths = sorted(tmp_dir.glob("fluid_field*.txt"))

    if not file_paths:
        print("未找到任何 fluid_field*.txt 文件。")
        print("请先运行: cargo test -p aero -- test_new_and_run test_with_expansion")
        return

    print(f"将处理 {len(file_paths)} 个文件:")
    for fp in file_paths:
        print(f"  - {fp}")

    for fp in file_paths:
        process_one_file(fp, out_dir)

    print(f"\n完成。输出目录: {out_dir}")


if __name__ == "__main__":
    main()
