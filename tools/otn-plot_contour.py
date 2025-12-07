import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.interpolate
import scienceplots  # noqa: F401

# 启用 scienceplots 风格（更学术美观）
plt.style.use(["science", "no-latex"])  # 或 'ieee', 'grid' 等
# plt.rcParams['text.usetex'] = True  # 如需 LaTeX 渲染可开启

# 全局配置（避免重复设置）
plt.rcParams.update(
    {
        "font.size": 12,
        "axes.labelsize": 14,
        "axes.titlesize": 14,
        "legend.fontsize": 12,
        "figure.figsize": (8, 6),
    }
)


def load_fluid_data(filepath: str) -> tuple[np.ndarray, ...]:
    """加载流场数据，自动跳过含 NaN 的行，并返回命名字段。

    返回: (x, y, p, rho, T, Ma)
    """
    if not os.path.isfile(filepath):
        raise FileNotFoundError(f"数据文件不存在: {filepath}")

    try:
        data = np.loadtxt(filepath, delimiter=",")
    except Exception as e:
        raise ValueError(f"无法加载数据文件 {filepath}: {e}")

    # 去除含 NaN/Inf 的行
    mask = np.isfinite(data).all(axis=1)
    data = data[mask]
    if data.size == 0:
        raise ValueError("数据为空或全为无效值（NaN/Inf）")

    # 字段映射（建议用 dict 或 namedtuple，此处保持简洁）
    x, y = data[:, 0], data[:, 1]
    p, rho, T = data[:, 4], data[:, 5], data[:, 6]
    Ma = data[:, 9]

    return x, y, p, rho, T, Ma


def plot_scatter_and_contour(
    x: np.ndarray,
    y: np.ndarray,
    value: np.ndarray,
    attr_name: str = "Ma",
    resolution: int = 256,
    n_levels: int = 30,
    cmap="viridis",  # 更推荐 viridis/plasma 等感知均匀色图
    split_at_y0: bool = True,
    save_path: str | None = None,
):
    """绘制散点 + 等值线/填充云图（支持分段插值防伪影）"""
    if len(x) == 0:
        raise ValueError("输入数据为空")

    fig, ax = plt.subplots()
    ax.set_aspect("equal")

    # 散点图（可选：降低 alpha 避免过密遮挡）
    scatter = ax.scatter(x, y, c="gray", s=2, alpha=0.4, label="Grid Points")

    # 检查有效范围
    x_min, x_max = x.min(), x.max()
    y_min, y_max = y.min(), y.max()
    if x_min == x_max or y_min == y_max:
        raise ValueError("x 或 y 坐标无变化，无法绘图")

    # 自动计算 colorbar 范围（排除极端离群值可选）
    v_min, v_max = np.nanmin(value), np.nanmax(value)
    if v_min == v_max:
        v_min -= abs(v_min) * 0.01 or 1e-3
        v_max += abs(v_max) * 0.01 or 1e-3

    levels = np.linspace(v_min, v_max, n_levels + 1)

    # 分段插值：按 y=0 拆分上下流场（避免跨激波插值失真）
    if split_at_y0:
        # 找到 y 最接近 0 的索引（更鲁棒：取 sign 变化处）
        # 简化版：按 y 排序后找零点分割
        idx_sorted = np.argsort(y)
        y_sorted = y[idx_sorted]
        sign_change = np.where(np.diff(np.sign(y_sorted)) != 0)[0]
        if len(sign_change) > 0:
            split_idx = idx_sorted[sign_change[0] + 1]
        else:
            split_idx = len(x) // 2  # fallback
    else:
        split_idx = len(x)

    def _contour_segment(x_seg, y_seg, val_seg):
        if len(x_seg) < 4:  # 插值至少需几个点
            return None
        try:
            xi = np.linspace(x_seg.min(), x_seg.max(), resolution)
            yi = np.linspace(y_seg.min(), y_seg.max(), resolution)
            Xi, Yi = np.meshgrid(xi, yi)
            Zi = scipy.interpolate.griddata(
                (x_seg, y_seg), val_seg, (Xi, Yi), method="linear", fill_value=np.nan
            )
            return ax.contourf(Xi, Yi, Zi, levels=levels, cmap=cmap, alpha=0.85)
        except Exception as e:
            print(f"插值失败（段）: {e}")
            return None

    cset = None
    # 上半场
    if split_idx > 0:
        cset = _contour_segment(x[:split_idx], y[:split_idx], value[:split_idx])
    # 下半场
    if split_idx < len(x):
        cset2 = _contour_segment(x[split_idx:], y[split_idx:], value[split_idx:])
        if cset is None and cset2 is not None:
            cset = cset2

    if cset is None:
        raise RuntimeError("所有插值段均失败，请检查数据分布")

    # 坐标轴
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_xlabel(r"$x / y_t$")
    ax.set_ylabel(r"$y / y_t$")
    ax.set_title(f"{attr_name} Contour")

    # Colorbar
    cbar = fig.colorbar(cset, ax=ax, orientation="horizontal", pad=0.12, aspect=30)
    cbar.set_label(attr_name, fontsize=14)
    # 智能 tick：约 5~7 个主刻度
    tick_step = max(1, n_levels // 6)
    cbar.set_ticks(levels[::tick_step])

    # 可选：添加等值线（增强结构）
    # ax.contour(Xi, Yi, Zi, levels=levels[::3], colors='k', linewidths=0.3, alpha=0.5)

    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"图像已保存至: {save_path}")

    plt.show()


def main():
    # 构建路径（更健壮：使用 pathlib）
    import pathlib

    tmp_dir = pathlib.Path(__file__).parent.resolve() / "../target/tmp"
    filepath = tmp_dir / "fluid_field.txt"
    filepath = filepath.resolve()

    try:
        x, y, p, rho, T, Ma = load_fluid_data(str(filepath))
        plot_scatter_and_contour(
            x,
            y,
            Ma,
            attr_name="Mach Number",
            resolution=256,  # 512 可能过大，256 足够且快；可调
            cmap="plasma",
            save_path=str(tmp_dir / "output/ma_contour.png"),  # 可选保存
        )
    except Exception as e:
        print(f"❌ 运行出错: {e}")
        import traceback

        traceback.print_exc()


if __name__ == "__main__":
    main()
