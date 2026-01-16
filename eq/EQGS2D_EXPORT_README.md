# EQGS2D 2D数据导出与可视化工具

## 概述

这个工具可以将 EQ 程序中的 2D 等离子体平衡数据导出为 CSV 文件，包括：
- 极向磁通 PSIRZ(R,Z) 在 R-Z 网格上的分布
- 磁通梯度 ∂ψ/∂R 和 ∂ψ/∂Z
- 磁面坐标 RPS(θ,ψ) 和 ZPS(θ,ψ)
- 磁面等高线数据
- 分界面（separatrix）边界
- 关键平衡参数

## 与 1D 导出的关系

本工具是 EQGS1D 导出功能的扩展：
- **EQGS1D 导出** (`SE` 命令): 23 个 1D 剖面参数
- **EQGS2D 导出** (`CE` 命令): 10 个 2D 场数据和几何数据

两者互补，提供完整的等离子体平衡数据。

---

## 使用方法

### 步骤 1: 编译（如果尚未完成）

```bash
cd /Users/dengxiaoya/TASK/latest/task/eq
make clean
make eq
```

### 步骤 2: 运行 EQ 程序并导出 2D 数据

#### 方式 A: 使用准备好的输入文件

```bash
./eq <in/eq.ITER01_export2d.in
```

#### 方式 B: 交互式运行

```bash
./eq
```

然后在程序中输入命令：
```
0                # 显示模式选择
f                # 文件模式
eq.ITER01.gs     # 图形文件名
c                # 继续
p                # 参数设置
[输入参数 namelist]
r                # 运行计算
c                # 继续
[输入其他参数]
g                # 进入图形菜单
c                # 选择计算
ce               # ← 新命令: 导出 2D 数据到 CSV
x                # 退出图形菜单
s                # 保存
q                # 退出
```

**新命令 `CE`** (C + E) 会导出所有 2D 平衡数据到 10 个 CSV 文件。

### 步骤 3: 验证导出的文件

```bash
ls -lh eqgs2d_*.csv
```

应该看到 10 个 CSV 文件：
```
eqgs2d_01_PSIRZ_grid.csv              - PSIRZ(R,Z) 极向磁通网格 (19KB)
eqgs2d_02_DPSIDR_grid.csv             - ∂ψ/∂R 径向梯度网格 (19KB)
eqgs2d_03_DPSIDZ_grid.csv             - ∂ψ/∂Z 垂直梯度网格 (19KB)
eqgs2d_04_RPS_flux_surfaces.csv       - RPS(θ,ψ) 磁面 R 坐标 (54KB)
eqgs2d_05_ZPS_flux_surfaces.csv       - ZPS(θ,ψ) 磁面 Z 坐标 (54KB)
eqgs2d_06_flux_surface_contours.csv   - 详细磁面等高线 (244KB)
eqgs2d_07_separatrix.csv              - 分界面上下边界 (4.7KB)
eqgs2d_08_parameters.csv              - 关键参数表 (351B)
eqgs2d_09_RZ_grids.csv                - R 和 Z 网格向量 (1.3KB)
eqgs2d_10_PSIP_profile.csv            - PSIP 归一化磁通剖面 (1.1KB)
```

### 步骤 4: 使用 Python 读取和可视化数据

#### 基本用法

```python
from read_eqgs2d_csv import load_psirz_grid, load_flux_surfaces, load_parameters

# 读取 PSIRZ 网格
R_grid, Z_grid, PSIRZ = load_psirz_grid()
print(f"R 范围: {R_grid.min():.2f} - {R_grid.max():.2f} m")
print(f"Z 范围: {Z_grid.min():.2f} - {Z_grid.max():.2f} m")
print(f"PSIRZ 形状: {PSIRZ.shape}")

# 读取磁面数据
theta, R_surfaces, Z_surfaces = load_flux_surfaces()
print(f"θ 点数: {len(theta)}")
print(f"磁面数: {R_surfaces.shape[1]}")

# 读取关键参数
params = load_parameters()
print(f"磁轴位置: R = {params['RAXIS']['value']:.3f} m")
print(f"磁轴位置: Z = {params['ZAXIS']['value']:.3e} m")
```

#### 绘制 PSIRZ 等高线图

```python
from read_eqgs2d_csv import plot_psirz_contour
import matplotlib.pyplot as plt

# 绘制带分界面和磁轴的等高线图
fig, ax = plot_psirz_contour(levels=20, save_fig=True)
plt.show()
```

#### 绘制磁面

```python
from read_eqgs2d_csv import plot_flux_surfaces
import matplotlib.pyplot as plt

# 每隔 5 条磁面绘制一条
fig, ax = plot_flux_surfaces(plot_every=5, save_fig=True)
plt.show()
```

#### 绘制所有 2D 数据组合图

```python
from read_eqgs2d_csv import plot_all_2d_data
import matplotlib.pyplot as plt

# 2x2 子图: PSIRZ, ∂ψ/∂R, ∂ψ/∂Z, 磁面
fig, axes = plot_all_2d_data(save_fig=True)
plt.show()
```

#### 运行完整示例

```bash
python read_eqgs2d_csv.py
```

这将：
1. 打印关键参数
2. 绘制 PSIRZ 等高线图（保存为 `eqgs2d_psirz_contour.png`）
3. 绘制磁面图（保存为 `eqgs2d_flux_surfaces.png`）
4. 绘制 2x2 组合图（保存为 `eqgs2d_all_plots.png`）
5. 显示所有图形

---

## 导出的 10 个文件详解

### 1. `eqgs2d_01_PSIRZ_grid.csv` - 极向磁通网格

**格式**: Z × R 矩阵（CSV 表格）

```csv
Z_m / R_m,3.67E+00,3.78E+00,...,8.59E+00
-4.46E+00,-1.03E+02,-1.03E+02,...,-1.03E+02
-4.19E+00,-1.03E+02,-1.03E+02,...,-1.03E+02
...
```

- **第一行**: R 网格坐标（列标题）
- **第一列**: Z 网格坐标（行索引）
- **数据**: PSIRZ(R,Z) 极向磁通值（单位: Wb）

**数据维度**: 33 × 33 网格

### 2-3. `eqgs2d_02_DPSIDR_grid.csv` / `eqgs2d_03_DPSIDZ_grid.csv`

**格式**: 同 PSIRZ_grid.csv

- **DPSIDR**: ∂ψ/∂R 径向磁通梯度（单位: Wb/m）
- **DPSIDZ**: ∂ψ/∂Z 垂直磁通梯度（单位: Wb/m）

**用途**: 计算磁场分量、电流密度等

### 4-5. `eqgs2d_04_RPS_flux_surfaces.csv` / `eqgs2d_05_ZPS_flux_surfaces.csv`

**格式**: θ × ψ 矩阵

```csv
theta_deg / surface,1,2,3,...,50
0.00E+00,6.38E+00,6.38E+00,...,7.84E+00
5.63E+00,6.38E+00,6.38E+00,...,7.85E+00
...
```

- **第一行**: 磁面序号（1 到 NRMAX）
- **第一列**: θ 角度（0° 到 360°，NTHMAX 个点）
- **数据**:
  - **RPS**: 每条磁面上各 θ 点的 R 坐标（单位: m）
  - **ZPS**: 每条磁面上各 θ 点的 Z 坐标（单位: m）

**数据维度**: 64 × 50（64 个 θ 点，50 条磁面）

**用途**: 绘制完整的磁面形状、计算几何参数

### 6. `eqgs2d_06_flux_surface_contours.csv` - 详细等高线

**格式**: 长表格（每行一个点）

```csv
surface,point,theta_deg,R_m,Z_m
1,1,0.00E+00,6.38E+00,-8.58E-15
1,2,5.63E+00,6.38E+00,6.26E-02
...
50,64,3.54E+02,7.84E+00,-6.38E-01
```

- **surface**: 磁面序号（1-50）
- **point**: 点序号（1-64）
- **theta_deg**: θ 角度
- **R_m**: R 坐标
- **Z_m**: Z 坐标

**数据量**: 50 × 64 = 3200 行

**用途**: 方便在 pandas 中进行过滤、分组分析

### 7. `eqgs2d_07_separatrix.csv` - 分界面

**格式**: 表格

```csv
RSU_m,ZSU_m,RSW_m,ZSW_m
7.84E+00,1.69E+00,7.84E+00,-1.69E+00
7.85E+00,1.69E+00,7.85E+00,-1.69E+00
...
```

- **RSU_m, ZSU_m**: 分界面上半部分的 (R, Z) 坐标
- **RSW_m, ZSW_m**: 分界面下半部分的 (R, Z) 坐标

**数据量**: NSUMAX 个点（本例 65 个）

**用途**: 标识等离子体边界

### 8. `eqgs2d_08_parameters.csv` - 关键参数

**格式**: 参数表

```csv
parameter,value,unit,description
RAXIS,6.34E+00,m,axis R
ZAXIS,-1.16E-14,m,axis Z
PSITA,1.20E+02,Wb,tor flux
PSIPA,1.03E+02,Wb,pol flux
PSI0,-1.03E+02,Wb,boundary
NRGMAX,33,-,R pts
NZGMAX,33,-,Z pts
NTHMAX,64,-,theta pts
NRMAX,50,-,radial surfs
NRPMAX,47,-,inside sep
NSUMAX,65,-,sep pts
```

**包含**:
- 磁轴位置 (RAXIS, ZAXIS)
- 环向/极向磁通 (PSITA, PSIPA)
- 边界磁通 (PSI0)
- 网格尺寸 (NRGMAX, NZGMAX, NTHMAX, NRMAX)
- 分界面点数 (NSUMAX)

### 9. `eqgs2d_09_RZ_grids.csv` - 网格向量

**格式**: 两列表格

```csv
index,R_m,Z_m
1,3.67E+00,-4.46E+00
2,3.78E+00,-4.19E+00
...
33,8.59E+00,4.46E+00
```

- **R_m**: R 方向网格坐标
- **Z_m**: Z 方向网格坐标

**用途**: 快速获取网格坐标向量

### 10. `eqgs2d_10_PSIP_profile.csv` - PSIP 剖面

**格式**: 两列表格

```csv
index,PSIP_normalized
1,0.00E+00
2,6.84E-02
...
50,1.00E+00
```

- **PSIP_normalized**: 归一化极向磁通（0.0 到 1.0）

**用途**: 将磁面序号映射到归一化磁通

---

## Python 函数说明

### 数据加载函数

```python
# 加载 PSIRZ 网格
R_grid, Z_grid, PSIRZ = load_psirz_grid(filename='eqgs2d_01_PSIRZ_grid.csv')

# 加载梯度网格
R_grid, Z_grid, DPSIDR = load_derivative_grid('eqgs2d_02_DPSIDR_grid.csv')
R_grid, Z_grid, DPSIDZ = load_derivative_grid('eqgs2d_03_DPSIDZ_grid.csv')

# 加载磁面数据
theta, R_surfaces, Z_surfaces = load_flux_surfaces(
    r_filename='eqgs2d_04_RPS_flux_surfaces.csv',
    z_filename='eqgs2d_05_ZPS_flux_surfaces.csv'
)

# 加载分界面
sep_df = load_separatrix(filename='eqgs2d_07_separatrix.csv')

# 加载参数
params = load_parameters(filename='eqgs2d_08_parameters.csv')
```

### 绘图函数

```python
# 绘制 PSIRZ 等高线
fig, ax = plot_psirz_contour(
    R_grid=None,      # 如果为 None，自动加载
    Z_grid=None,
    PSIRZ=None,
    levels=20,        # 等高线数量
    figsize=(12, 10),
    save_fig=False    # 是否保存图像
)

# 绘制磁面
fig, ax = plot_flux_surfaces(
    theta=None,
    R_surfaces=None,
    Z_surfaces=None,
    plot_every=5,     # 每隔几条磁面绘制一条
    figsize=(12, 10),
    save_fig=False
)

# 绘制所有 2D 数据组合图
fig, axes = plot_all_2d_data(save_fig=False)
```

---

## 常见问题

### Q: 如何从 2D 网格中提取特定位置的值？

A: 使用 pandas 的索引功能或 NumPy 插值：

```python
from read_eqgs2d_csv import load_psirz_grid
from scipy.interpolate import RectBivariateSpline

R_grid, Z_grid, PSIRZ = load_psirz_grid()

# 创建插值函数
interp = RectBivariateSpline(Z_grid, R_grid, PSIRZ.T)

# 在任意 (R, Z) 点插值
R_point = 6.2
Z_point = 0.0
psi_value = interp(Z_point, R_point)[0, 0]
print(f"PSIRZ at (R={R_point}, Z={Z_point}) = {psi_value:.4f} Wb")
```

### Q: 如何计算磁面的几何参数（椭圆度、三角形度）？

A: 从 RPS 和 ZPS 数据中提取：

```python
theta, R_surfaces, Z_surfaces = load_flux_surfaces()

# 选择某条磁面（例如第 25 条，约 0.5 ψ）
i_surf = 24  # 索引从 0 开始
R = R_surfaces[:, i_surf]
Z = Z_surfaces[:, i_surf]

# 计算几何参数
R_max = R.max()
R_min = R.min()
Z_max = Z.max()
Z_min = Z.min()

R_major = (R_max + R_min) / 2  # 主半径
a = (R_max - R_min) / 2        # 小半径
b = (Z_max - Z_min) / 2        # 垂直半径
kappa = b / a                   # 椭圆度

print(f"主半径 R0 = {R_major:.3f} m")
print(f"小半径 a = {a:.3f} m")
print(f"椭圆度 κ = {kappa:.3f}")
```

### Q: 如何绘制特定 PSIRZ 等值线？

A: 使用 matplotlib 的 contour 函数：

```python
import matplotlib.pyplot as plt
import numpy as np

R_grid, Z_grid, PSIRZ = load_psirz_grid()
R, Z = np.meshgrid(R_grid, Z_grid)

fig, ax = plt.subplots(figsize=(10, 8))

# 绘制特定值的等值线
psi_levels = [-100, -80, -60, -40, -20, 0]  # 指定 PSIRZ 值
contour = ax.contour(R, Z, PSIRZ.T, levels=psi_levels, colors='blue')
ax.clabel(contour, inline=True, fontsize=10)

ax.set_xlabel('R (m)')
ax.set_ylabel('Z (m)')
ax.set_aspect('equal')
plt.show()
```

### Q: 数据和 GS 图像中的一致吗？

A: 完全一致。导出的数据直接来自 `eqcomq.inc` 中的公共变量，与 `EQGS2D` 绘图函数使用相同的数据源。

---

## 技术细节

### 数据来源

所有数据来自 `eqcomq.inc` 中定义的公共变量：

```fortran
COMMON /Q02/ PSIRZ(81,81),DPSIDR(81,81),DPSIDZ(81,81),
           & RG(81),ZG(81),NRGMAX,NZGMAX
COMMON /Q07/ RPS(129,51),ZPS(129,51),PSIP(51),NTHMAX,NRMAX,NRPMAX
COMMON /Q08/ RRPSIAX,ZZPSIAX,PSITAX,PSIMAX,PSI0,RMAXIS,ZMAXIS
COMMON /Q10/ RSU(129),ZSU(129),RSW(129),ZSW(129),NSUMAX
```

这些变量在 `EQCALQ` 子程序中计算。

### 与 EQGS2D 的关系

`EQGS2D_EXPORT` 子程序直接读取与 `EQGS2D` 相同的公共变量，因此导出的数据与 GS 图像中绘制的数据完全一致。

### 坐标系统

- **R**: 主半径方向（大半径，从中心轴向外）
- **Z**: 垂直方向（向上为正）
- **θ**: 极向角（从外中平面逆时针，0° 到 360°）
- **ψ**: 极向磁通（归一化后 0.0 到 1.0）

---

## 文件列表

### Fortran 源代码

- **`eqgs2d_export.f`**: 2D 数据导出子程序
  - `SUBROUTINE EQGS2D_EXPORT`: 导出 10 个 CSV 文件

- **`eqgout.f`** (已修改):
  - 添加 'CE' 命令处理（第 41-43 行）

- **`Makefile`** (已修改):
  - 添加 `eqgs2d_export.f` 到 SRCS 列表

### Python 脚本

- **`read_eqgs2d_csv.py`**: 数据读取和可视化工具
  - `load_psirz_grid()`: 读取 PSIRZ 网格
  - `load_derivative_grid()`: 读取梯度网格
  - `load_flux_surfaces()`: 读取磁面坐标
  - `load_separatrix()`: 读取分界面
  - `load_parameters()`: 读取关键参数
  - `plot_psirz_contour()`: 绘制等高线图
  - `plot_flux_surfaces()`: 绘制磁面
  - `plot_all_2d_data()`: 绘制组合图

### 输入文件

- **`in/eq.ITER01_export2d.in`**: 自动化输入文件（包含 'ce' 命令）

### 文档

- **`EQGS1D_EXPORT_README.md`**: 1D 剖面导出文档
- **`EQGS2D_EXPORT_README.md`**: 本文档（2D 数据导出）
- **`gs_图像分析.md`**: GS 图像分析报告

---

## 示例工作流

### 完整分析流程

```bash
# 1. 运行 eq 并导出数据
cd /Users/dengxiaoya/TASK/latest/task/eq
./eq <in/eq.ITER01_export2d.in

# 2. 检查导出的文件
ls -lh eqgs2d_*.csv

# 3. 运行 Python 可视化
python read_eqgs2d_csv.py

# 4. 查看生成的图像
open eqgs2d_psirz_contour.png
open eqgs2d_flux_surfaces.png
open eqgs2d_all_plots.png
```

### 自定义分析

```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from read_eqgs2d_csv import *

# 加载所有数据
R_grid, Z_grid, PSIRZ = load_psirz_grid()
theta, R_surfaces, Z_surfaces = load_flux_surfaces()
sep = load_separatrix()
params = load_parameters()

# 提取磁轴位置
R_axis = params['RAXIS']['value']
Z_axis = params['ZAXIS']['value']

# 计算磁轴处的场
from scipy.interpolate import RectBivariateSpline
interp = RectBivariateSpline(Z_grid, R_grid, PSIRZ.T)
psi_axis = interp(Z_axis, R_axis)[0, 0]
print(f"磁轴磁通: {psi_axis:.4f} Wb")

# 分析磁面形状
i_surf = 30  # 第 30 条磁面
R = R_surfaces[:, i_surf]
Z = Z_surfaces[:, i_surf]
kappa = (Z.max() - Z.min()) / (R.max() - R.min())
print(f"第 {i_surf+1} 条磁面的椭圆度: {kappa:.3f}")

# 自定义绘图
fig, ax = plt.subplots(figsize=(10, 8))
ax.contour(R_grid, Z_grid, PSIRZ.T, levels=30, cmap='viridis')
ax.plot(R, Z, 'r-', linewidth=2, label=f'Surface {i_surf+1}')
ax.plot(R_axis, Z_axis, 'k*', markersize=15, label='Magnetic Axis')
ax.set_xlabel('R (m)')
ax.set_ylabel('Z (m)')
ax.set_aspect('equal')
ax.legend()
ax.grid(True, alpha=0.3)
plt.show()
```

---

## 下一步

1. **深入分析**: 使用导出的数据进行物理分析（稳定性、输运等）
2. **与实验对比**: 将计算结果与实验测量对比
3. **参数扫描**: 批量运行不同参数，分析趋势
4. **3D 可视化**: 使用 mayavi 或 plotly 进行 3D 可视化

---

## 相关链接

- **1D 导出**: 参见 `EQGS1D_EXPORT_README.md`
- **GS 分析**: 参见 `gs_图像分析.md`
- **源代码**: `eqgs2d_export.f`, `eqgout.f`

---

**创建日期**: 2026-01-14
**版本**: 1.0
