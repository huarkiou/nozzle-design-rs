[![BuildTest](https://github.com/huarkiou/nozzle-design-rs/actions/workflows/buildtest.yml/badge.svg)](https://github.com/huarkiou/nozzle-design-rs/actions/workflows/buildtest.yml)

# nozzle-design-rs

基于特征线法 (Method of Characteristics, MOC) 的超声速喷管扩张段设计工具集，使用 Rust 编写，提供 **最大推力喷管 (OTN)** 生成与 **三维流线追踪喷管 (SLTN)** 设计两大核心功能。

---

## 目录

- [背景](#背景)
- [项目架构](#项目架构)
- [依赖与构建](#依赖与构建)
- [应用](#应用)
  - [otn — 最优推力喷管](#otn--最优推力喷管)
  - [sltn — 三维流线追踪喷管](#sltn--三维流线追踪喷管)
- [核心库](#核心库)
- [Python 工具](#python-工具)
- [配置参考](#配置参考)
- [运行测试](#运行测试)
- [参考文献](#参考文献)

---

## 背景

超声速喷管是火箭发动机、高超声速风洞、空气动力学实验台等设备的核心部件。喷管气动型面的设计直接决定了推力效率与流场品质。

本项目实现了基于特征线法 (MOC) 的喷管型面设计：

- **二维轴对称 / 平面流动**：通过特征线网格求解超声速流场，沿壁面追踪特征线，获得喷管型面坐标。
- **三维型面设计**：在二维轴对称基准流场的基础上，通过流线追踪技术，将二维结果拓展至三维空间，生成任意截面形状的三维喷管。

支持的主要物理模型：

| 模型 | 说明 |
|------|------|
| 无旋特征线法 | 假设流场无旋，适用于均匀来流 |
| 有旋特征线法 | 考虑流场涡量，适用于非均匀来流(未实现) |
| 常比热 / 变比热 | 支持分段多项式、NASA 9系数等变比热模型 |
| 轴对称 / 平面 | 支持轴对称流动与二维平面流动 |

---

## 项目架构

```
nozzle-design-rs/
├── crates/                    # Rust 核心库
│   ├── math/                  # 数值计算基础库
│   ├── geometry/              # 几何形状定义与处理
│   └── aero/                  # 气动热力学核心
├── apps/                      # 命令行应用
│   ├── otn/                   # 最优推力喷管 (Optimal Thrust Nozzle)
│   └── sltn/                  # 三维流线追踪喷管 (StreamLine Tracing Nozzle)
├── tests/                     # 项目级集成测试
│   ├── common/                # 测试公共模块
│   ├── otn_parameters.rs      # OTN 参数覆盖测试
│   ├── sltn_shapes.rs         # SLTN 截面形状覆盖测试
│   ├── otn_sltn_sern.rs       # OTN→SLTN 完整工作流（SERN 设计）
│   ├── otn_sltn_toml.rs       # TOML 配置文件继承流程测试
│   └── fixtures/              # TOML 配置模板
├── tools/                     # Python 可视化与数据后处理脚本
├── .github/workflows/         # CI/CD (GitHub Actions)
├── Cargo.toml                 # Rust workspace 配置
├── pyproject.toml             # Python 工具依赖
└── readme.md
```

### 模块依赖关系

```
apps/otn ───────┐
                ├── aero ──────┬── geometry
apps/sltn ──────┘              └── math
```

- `math`：最底层，提供数值求根（二分法、TOMS748、牛顿二维）、高斯-勒让德积分、分段多项式、二维坐标/折线等工具。
- `geometry`：二维/三维几何抽象层，提供闭合曲线（圆、椭圆、超椭圆、矩形、自定义）、三维点、壁面点集处理、OBJ 模型导出等。
- `aero`：气动热力学核心，实现特征线法 (MOC)、等熵关系、材料属性、流线追踪。
- `otn` / `sltn`：面向最终用户的命令行应用，解析 TOML 配置文件并调用 `aero` 库执行设计。

---

## 依赖与构建

### 系统要求

- **Rust** ≥ 1.85 (edition 2024)
- **Python** ≥ 3.13（仅用于 `tools/` 目录中的可视化脚本）

### 编译

```bash
# 克隆仓库
git clone https://github.com/huarkiou/nozzle-design-rs.git
cd nozzle-design-rs

# 编译所有程序（release 模式，启用 thin LTO）
cargo build --release

# 编译指定程序
cargo build -p otn --release
cargo build -p sltn --release
```

### 安装 Python 工具依赖

```bash
# 使用 uv（推荐）
uv sync
```

---

## 应用

### otn — 最大推力喷管

基于 **特征线法 (MOC)** 生成二维轴对称/平面最优推力喷管型面。

```bash
# 生成默认配置文件
cargo run -p otn --release -- init

# 运行喷管设计（默认读取 ./otn.toml）
cargo run -p otn --release

# 指定配置文件
cargo run -p otn --release -- otn.toml
```

**输出文件**：
- `<prefix>field_data.txt` — 完整特征线流场数据
- `<prefix>geo_all.dat` — 喷管壁面几何坐标（UG NX `.dat` 格式）

**设计流程**：

1. **初值线** (Initial Line)：从喉部出发，计算初始特征线
2. **初值问题** (Initial Section)：由初值线向下游推进
3. **膨胀段** (Expansion Section)：控制初始膨胀角 `theta_a`，形成膨胀波系
4. **转向段** (Transition Section)：通过二分搜索确定转向段起始点，使出口达到目标长度
5. **均一区** (Uniform Section)：填补转向段与膨胀段之间的缺口

当 `theta_a` 设置为 `NaN` 或负值时，程序会自动迭代搜索满足出口背压（或出口高度）约束的最优初始膨胀角。

### sltn — 三维流线追踪喷管

在 OTN 生成的二维轴对称基准流场基础上，通过 **流线追踪** 技术生成任意截面形状的三维喷管。

```bash
# 生成默认配置文件
cargo run -p sltn --release -- init

# 运行流线追踪（默认读取 ./sltn.toml）
cargo run -p sltn --release

# 指定配置文件
cargo run -p sltn --release -- sltn.toml
```

**输出文件**：
- `model.dat` / `model.obj` — 喷管主型面
- `downstream.dat` / `downstream.obj` — 下游延伸段
- `upstream.dat` / `upstream.obj` — 上游收缩段

**支持的入口/出口截面形状**：

| 形状 | 说明 |
|------|------|
| `Circle` | 圆形截面 |
| `Ellipse` | 椭圆截面 |
| `SuperEllipse` | 超椭圆截面 |
| `Rectangular` | 矩形截面 |
| `UserDefined` | 自定义离散点 |

---

## 核心库

### `math` — 数值计算基础库

| 模块 | 说明 |
|------|------|
| `rootfinding` | 二分法 (bisection)、TOMS748 (有界求根)、牛顿二维求解器 |
| `quadrature` | 高斯-勒让德数值积分 |
| `polynomial` | 多项式与分段多项式支持 |
| `geometry` | 二维坐标点 (`Coord2d`)、折线 (`Polyline`) |
| `tolerance` | 容差抽象 |

### `geometry` — 几何形状定义与处理

| 模块 | 说明 |
|------|------|
| `ClosedCurve` trait | 闭合曲线抽象接口，按极角生成均匀采样点 |
| `Circle` / `Ellipse` / `SuperEllipse` / `Rectangular` | 标准截面形状 |
| `UserDefined` | 自定义离散点截面 |
| `Point3d` / `WallPoints` | 三维坐标点与壁面点集处理 |
| `wallpoints` | 插值、重采样、合并、文件 I/O 等工具函数 |
| `obj` | OBJ 三维模型导出 |

### `aero` — 气动热力学核心

| 模块 | 说明 |
|------|------|
| `moc` | 特征线法核心：特征点 (`MocPoint`)、特征线 (`CharLine`/`CharLines`)、无旋/有旋基本过程 (`UnitProcess`) |
| `nozzle` | 喷管各段计算：初值线、初值问题、膨胀段、转向段、均一区、约束喷管 (`ConstraintNozzle`) |
| `streamline_trace` | 三维流线追踪：沿基准流场向前/向后追踪流线，通过加权函数控制型面过渡 |
| `material` | 工质属性：定压比热、比热比、气体常数、焓值，支持多种空气模型 |
| `isentropic` | 等熵关系：总温/静温互算 |

#### 喷管各段说明

| 段 | `nozzle` 模块 | 说明 |
|----|--------------|------|
| 初值问题 | `InitialLine` + `InitialSection` | 从喉部出发的初始特征线网络 |
| 膨胀段 | `ExpansionSection` | 由初始膨胀角 `theta_a` 控制的膨胀波系 |
| 转向段 | `TransitionSection` | 将超声速气流转向至平行出口方向 |
| 均一区 | `UniformSection` | 填补各段之间的流场缺口 |

---

## Python 工具

`tools/` 目录提供用于数据后处理与可视化的 Python 脚本：

| 脚本 | 说明 |
|------|------|
| `otn-plot_contour.py` | OTN 流场马赫数/压力等值线图 |
| `plot_mocpoints.py` | 特征线网格可视化 |
| `sltn-c2e-crosssection_plot.py` | SLTN 圆形→椭圆截面过渡图 |
| `sltn-c2r-crosssection_plot.py` | SLTN 圆形→矩形截面过渡图 |
| `sltn-c2t-crosssection_plot.py` | SLTN 圆形→三角形截面过渡图 |
| `sltn-weight_funcs_plot.py` | 流线追踪加权函数可视化 |
| `cal_ellipse_area.py` | 椭圆面积计算 |
| `fit_circle.py` | 圆拟合 |
| `fit_ellipse.py` | 椭圆拟合 |
| `fit_plane.py` | 平面拟合 |

---

## 配置参考

### OTN 配置文件（otn.toml）

```toml
[MOCControl]
irrotational = true          # 无旋特征线法 (true) 或有旋特征线法 (false)
axisymmetric = true          # 二维轴对称 (true) 或二维平面 (false)
eps = 1e-5                   # 残差收敛容差
n_correction_max = 20        # 欧拉预估校正最大迭代次数
n_inlet = 61                 # 入口初值线网格点数

[Material]
molecular_weight = 28.968    # 摩尔质量 (kg/kmol)
cp = nan                     # 定压比热容 J/(kg·K)，nan 使用 NASA 9 空气模型

[Geometry]
height = 1.0                 # 进口高度/半径 (m)
height_e = NaN               # 目标出口高度/半径 (m)，NaN 表示以最大推力为目标
length = 6.0                 # 喷管目标长度 (m)
width = 1.0                  # 横向宽度 (仅 axisymmetric=false 时有效)

[Inlet]
p_total = 800000.0           # 来流总压 (Pa)
T_total = 2000.0             # 来流总温 (K)
Ma = 1.2                     # 来流马赫数 (必须 ≥ 1)
theta = 0.0                  # 来流方向角 (rad)，当前仅支持 0

[Throat]
R_t = 0.0                    # 喉部过渡圆弧半径 (m)
theta_a = NaN                # 初始膨胀角 (rad)，NaN 时自动迭代

[Outlet]
p_ambient = 7000.0           # 设计出口背压 (Pa)

[IO]
output_prefix = ""           # 输出文件前缀
```

#### Material 配置方式

```toml
# 方式 1：常数比热容
[Material]
molecular_weight = 28.968
cp = 1004.675                 # 常数定压比热容 J/(kg·K)

# 方式 2：内置 NASA 9 系数变比热空气模型
[Material]
molecular_weight = 28.968
cp = nan                      # nan 触发内置 NASA 9 系数空气模型

# 方式 3：自定义分段多项式变比热容
[Material]
molecular_weight = 28.968
# pos_coefficients[0]=T⁰（常数项）, [1]=T¹, [2]=T², …
# neg_coefficients[0]=T⁻¹, [1]=T⁻², …
[[Material.cp_segments]]
t_min = 200.0
t_max = 1000.0
pos_coefficients = [1437.799, -1.653609, 0.003062254, -2.279138e-06, 6.272365e-10]
neg_coefficients = [-56496.26, 2898903.0]
[[Material.cp_segments]]
t_min = 1000.0
t_max = 6000.0
pos_coefficients = [1476.665, -0.06138349, 2.027963e-05, -3.075525e-09, 1.888054e-13]
neg_coefficients = [-361053.2, 69324940.0]
# 向后兼容：pos_coefficients 为空时回退到材料级 cp 作为常数项
```

### SLTN 配置文件（sltn.toml）

```toml
[Control]
n_theta = 66                 # 周向采样点数
n_axis = 111                 # 轴向采样点数
monotonic = false            # 强制单调性
weight_parameter_a = 0       # 加权参数 a
export_obj = false           # 是否导出 OBJ 模型

[BaseFluidField]
axisymmetric = true          # 基准流场类型
datasource_inlet = './field_data.txt'
datasource_outlet = './field_data.txt'

[Inlet]
normalized = true            # 坐标归一化
shape = 'circle'
center = [0, 0]
radius = 0.9

[Outlet]
normalized = true
shape = 'ellipse'
center = [0, 0]
a = 0.6
b = 0.5
alpha = 0                    # 旋转角 (°)
```

---

## 运行测试

```bash
# 运行所有单元测试
cargo test --workspace

# 运行指定 crate 的测试
cargo test -p aero
cargo test -p math
cargo test -p geometry

# 运行项目级集成测试（含 --ignored，耗时较长）
cargo test -- --ignored

# 运行特定分类的集成测试
cargo test --test otn_parameters -- --ignored      # OTN 参数覆盖
cargo test --test sltn_shapes -- --ignored          # SLTN 截面形状覆盖
cargo test --test otn_sltn_sern -- --ignored        # SERN 设计工作流
cargo test --test otn_sltn_toml -- --ignored        # TOML 配置文件流程

# 包含详细输出
cargo test --workspace -- --nocapture

# 运行 benchmark（需 nightly toolchain）
cargo bench -p aero
```

### 集成测试分类

| 测试文件 | 说明 | 测试数 |
|---------|------|--------|
| `otn_parameters.rs` | OTN 参数覆盖：流动类型、材料模型、几何、进口/出口条件 | 15 |
| `sltn_shapes.rs` | SLTN 截面覆盖：圆形/椭圆/矩形/超椭圆、混合形状、平面流场、分辨率 | 16 |
| `otn_sltn_sern.rs` | OTN→SLTN 完整工作流，覆盖典型超燃冲压发动机 SERN 设计 | 9 |
| `otn_sltn_toml.rs` | TOML 配置文件输入路径，覆盖三种材料模型 | 5 |

---

## 参考文献

1. Zucrow, M. J., & Hoffman, J. D. (1976). *Gas Dynamics* (Vol. 1 & 2). John Wiley & Sons.
2. Anderson, J. D. (2003). *Modern Compressible Flow: With Historical Perspective* (3rd ed.). McGraw-Hill.
3. Shapiro, A. H. (1953). *The Dynamics and Thermodynamics of Compressible Fluid Flow*. Ronald Press.
