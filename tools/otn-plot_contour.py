"""OTN 喷管特征线法流场可视化。

生成两幅图:
1. 特征线散点图 — 每条右行特征线用不同颜色连线+标记点,
   用于观察特征线走向、间距, 判断算法逻辑是否存在问题。
2. 马赫数云图 — 基于 Delaunay 三角剖分的填充等值线,
   用于观察结果数值是否正确。

用法:
    python tools/otn-plot_contour.py
    (需先从项目根目录执行 `cargo test -p aero -- test_new_and_run`
     生成 target/tmp/fluid_field.txt)
"""

import os
import pathlib

import matplotlib.pyplot as plt
import numpy as np
import scienceplots  # noqa: F401
from matplotlib.tri import Triangulation

plt.style.use(["science", "no-latex"])
plt.rcParams.update({
    "font.size": 12,
    "axes.labelsize": 14,
    "axes.titlesize": 14,
    "legend.fontsize": 10,
    "figure.figsize": (10, 6),
})


# ── 数据加载 ──────────────────────────────────────────────────

def load_fluid_data_grouped(filepath: str) -> list[np.ndarray]:
    """按空行分组读取流场数据, 保留每条特征线的独立结构。

    fluid_field.txt 格式: 逗号分隔, 10 列, 空行分隔不同的特征线组。
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
                    # 筛掉含 NaN/Inf 的行
                    arr = arr[np.isfinite(arr).all(axis=1)]
                    if len(arr) >= 2:
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
            if len(arr) >= 2:
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
    """逐条特征线绘制连线+散点, 用渐变色区分不同线。

    Args:
        groups: load_fluid_data_grouped 的输出, 每个 ndarray 对应一条特征线。
        ax: matplotlib Axes 对象。
        title: 图标题。
    """
    n = len(groups)
    # 用 tab10 循环色 + 深色到浅色 (wall→axis)
    colors = plt.cm.plasma(np.linspace(0.05, 0.95, n))

    for i, line in enumerate(groups):
        x, y = line[:, 0], line[:, 1]
        ax.plot(x, y, "-o",
                markersize=2.0 if i == 0 else 1.2,
                linewidth=0.6 if i == 0 else 0.4,
                color=colors[i],
                alpha=0.8)

    # 标注初值线
    if n > 0:
        ivl = groups[0]
        ax.plot(ivl[:, 0], ivl[:, 1], "-o",
                markersize=2.5, linewidth=1.0,
                color="black", alpha=0.9, label="Initial Value Line (IVL)")

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
    """基于 Delaunay 三角剖分的填充等值线图。

    Args:
        x, y: 所有散点的坐标 (一维数组)。
        values: 对应每个点的标量值 (一维数组)。
        ax: matplotlib Axes 对象。
        attr_name: colorbar 标签。
        n_levels: 等值线层数。
        cmap: colormap 名称。
    """
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

    cbar = plt.colorbar(tcf, ax=ax, orientation="horizontal",
                        pad=0.12, aspect=30)
    cbar.set_label(attr_name, fontsize=12)
    tick_step = max(1, n_levels // 6)
    cbar.set_ticks(levels[::tick_step])


# ── 主流程 ────────────────────────────────────────────────────

def main():
    proj_root = pathlib.Path(__file__).resolve().parent.parent
    data_path = proj_root / "target" / "tmp" / "fluid_field.txt"
    out_dir = proj_root / "target" / "tmp" / "output"
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"读取数据: {data_path}")
    try:
        groups = load_fluid_data_grouped(str(data_path))
    except FileNotFoundError:
        print(f"❌ 数据文件不存在. 请先运行: cargo test -p aero -- test_new_and_run")
        return
    except ValueError as e:
        print(f"❌ 数据错误: {e}")
        return

    print(f"共 {len(groups)} 条特征线 (第一条为初值线 IVL)")

    # 统计每条特征线的点数
    for i, g in enumerate(groups):
        label = "IVL (初值线)" if i == 0 else f"特征线 #{i}"
        print(f"  {label}: {len(g)} 个点, "
              f"x=[{g[:, 0].min():.4f}, {g[:, 0].max():.4f}], "
              f"y=[{g[:, 1].min():.4f}, {g[:, 1].max():.4f}]")

    # 摊平所有点用于轮廓图 (保留分组结构用于散点图)
    all_pts = np.vstack(groups)
    x_all, y_all = all_pts[:, 0], all_pts[:, 1]
    Ma_all = all_pts[:, 9]

    print(f"总点数: {len(all_pts)}, "
          f"Ma 范围: [{Ma_all.min():.4f}, {Ma_all.max():.4f}]")

    # ── 图1: 特征线散点图 ──
    fig1, ax1 = plt.subplots()
    plot_charline_scatter(groups, ax1,
                          title="MOC Characteristic Lines (Right-Running)")
    fig1.tight_layout()
    scatter_path = out_dir / "charline_scatter.png"
    fig1.savefig(scatter_path, dpi=200, bbox_inches="tight")
    print(f"散点图已保存: {scatter_path}")

    # ── 图2: 马赫数云图 ──
    fig2, ax2 = plt.subplots()
    plot_contour_tri(x_all, y_all, Ma_all, ax2,
                     attr_name="Mach Number", cmap="plasma")
    fig2.tight_layout()
    contour_path = out_dir / "ma_contour.png"
    fig2.savefig(contour_path, dpi=200, bbox_inches="tight")
    print(f"云图已保存: {contour_path}")

    plt.show()


if __name__ == "__main__":
    main()
