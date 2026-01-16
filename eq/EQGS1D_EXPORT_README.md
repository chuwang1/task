# EQGS1D 数据导出与分析工具

## 概述

这个工具可以将 EQ 程序中的 23 个等离子体平衡剖面参数导出为 CSV 文件，并提供 Python 脚本进行数据读取和分析。

## 方案对比

### 我的原始 CSV (`eqgr1d_parameters.csv`)
- **创建方式**: 手动分析源代码 `eqgout.f` 中的 `EQGR1D` 调用
- **内容**: 包含每个图的元数据信息（位置、标题、参数名等）
- **用途**: 了解每个图形的绘制参数和位置

### Fortran 导出的 CSV (`eqgs1d_*.csv`)
- **创建方式**: 在 EQ 程序运行时通过 Fortran 代码导出
- **内容**: 包含实际的物理量数值数据
- **用途**: 数据分析、绘图、后处理

**两者是互补的**：前者是"元数据表"，后者是"实际数据"。

---

## 使用方法

### 步骤 1: 修改代码（已完成）

我已经在代码中添加了导出功能：

1. **新增文件**: `eqgs1d_export.f` - 数据导出子程序
2. **修改文件**: `eqgout.f` - 添加 'SE' 命令选项
3. **更新**: `Makefile` - 包含新的源文件

### 步骤 2: 重新编译（已完成）

```bash
cd /Users/dengxiaoya/TASK/latest/task/eq
make clean
make eq
```

### 步骤 3: 运行 EQ 程序并导出数据

有两种方式：

#### 方式 A: 使用准备好的输入文件

```bash
./eq <in/eq.ITER01_export.in
```

#### 方式 B: 交互式运行

```bash
./eq
```

然后在程序中按以下顺序输入命令：
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
se               # ← 新命令: 导出数据到 CSV
x                # 退出图形菜单
s                # 保存
q                # 退出
```

**新命令 `SE`** (S + E) 会导出所有 23 个剖面到 CSV 文件。

### 步骤 4: 验证导出的文件

```bash
ls -lh eqgs1d_*.csv
```

应该看到 23 个 CSV 文件：
```
eqgs1d_01_PPS.csv      - 压强剖面
eqgs1d_02_TTS.csv      - 环向磁通剖面
eqgs1d_03_RLEN.csv     - 磁力线长度
...
eqgs1d_23_AVEGV.csv    - ⟨|∇V|⟩ 平均值
```

### 步骤 5: 使用 Python 读取和分析数据

#### 基本用法

```python
from read_eqgs1d_csv import load_profile, load_all_profiles

# 读取单个剖面（例如压强）
pps = load_profile(1)
print(pps.head())

# 读取所有剖面
data = load_all_profiles()
```

#### 绘图

```python
from read_eqgs1d_csv import plot_profile, plot_all_profiles
import matplotlib.pyplot as plt

# 绘制单个剖面
plot_profile(1)  # 绘制压强剖面
plt.show()

# 绘制所有 23 个剖面
plot_all_profiles(save_fig=True)
plt.show()
```

#### 导出到 Excel

```python
from read_eqgs1d_csv import export_to_excel

export_to_excel(output_file='my_eqgs1d_data.xlsx')
```

#### 统计分析

```python
from read_eqgs1d_csv import summary_statistics

summary_statistics()
```

#### 运行完整示例

```bash
python read_eqgs1d_csv.py
```

这将：
1. 读取所有 CSV 文件
2. 显示统计摘要
3. 绘制示例图 (PPS)
4. 导出所有数据到 Excel
5. 绘制所有 23 个剖面的总图

---

## 导出的 23 个参数

| 序号 | 文件名 | 变量名 | 描述 | 单位 |
|------|--------|--------|------|------|
| 1 | eqgs1d_01_PPS.csv | PPS | 压强剖面 | MPa |
| 2 | eqgs1d_02_TTS.csv | TTS | 环向磁通剖面 | Wb |
| 3 | eqgs1d_03_RLEN.csv | RLEN | 磁力线长度 | m |
| 4 | eqgs1d_04_QPS.csv | QPS | 安全因子 q | - |
| 5 | eqgs1d_05_VPS.csv | VPS | 体积剖面 | m³ |
| 6 | eqgs1d_06_SPS.csv | SPS | 磁面面积剖面 | m² |
| 7 | eqgs1d_07_MINMAX.csv | RRMIN/MAX,ZZMIN/MAX,BBMIN/MAX | R,Z,B 最小/最大值 | m,m,T |
| 8 | eqgs1d_08_BB.csv | BBMIN,BBMAX | 磁场强度最小/最大值 | T |
| 9 | eqgs1d_09_AVERR2.csv | AVERR2 | ⟨R²⟩ 磁面平均 | m² |
| 10 | eqgs1d_10_AVEIR2.csv | AVEIR2 | ⟨1/R²⟩ 磁面平均 | 1/m² |
| 11 | eqgs1d_11_AVEBB2.csv | AVEBB2 | ⟨B²⟩ 磁面平均 | T² |
| 12 | eqgs1d_12_AVEIB2.csv | AVEIB2 | ⟨1/B²⟩ 磁面平均 | 1/T² |
| 13 | eqgs1d_13_AVEGV2.csv | AVEGV2 | ⟨\|∇V\|²⟩ 磁面平均 | m⁴ |
| 14 | eqgs1d_14_AVEGVR2.csv | AVEGVR2 | ⟨\|∇V\|²/R²⟩ 磁面平均 | m² |
| 15 | eqgs1d_15_AVEGP2.csv | AVEGP2 | ⟨\|∇P\|²⟩ 磁面平均 | Pa²/m² |
| 16 | eqgs1d_16_AVEJ.csv | AVEJPR,AVEJTR | 平均极向/环向电流密度 | MA/m² |
| 17 | eqgs1d_17_RRPSI.csv | RRPSI | 主半径剖面 | m |
| 18 | eqgs1d_18_RSPSI.csv | RSPSI | 次半径剖面 | m |
| 19 | eqgs1d_19_ELIPPSI.csv | ELIPPSI | 椭圆度剖面 | - |
| 20 | eqgs1d_20_TRIGPSI.csv | TRIGPSI | 三角形度剖面 | - |
| 21 | eqgs1d_21_DVDPSIP.csv | DVDPSIP | dV/dψp | m³/Wb |
| 22 | eqgs1d_22_DVDPSIT.csv | DVDPSIT | dV/dψt | m³/Wb |
| 23 | eqgs1d_23_AVEGV.csv | AVEGV | ⟨\|∇V\|⟩ 磁面平均 | m² |

所有文件的 X 轴都是**归一化极向磁通 (PSIP)**。

---

## CSV 文件格式

每个 CSV 文件的格式：
```csv
PSIP,<变量名>_<单位>
0.00000000E+00,0.63999999E+00
0.06844354E+00,0.63936937E+00
...
```

- 第一列：归一化极向磁通 PSIP
- 后续列：对应的物理量（带单位后缀）
- 科学计数法格式 (Fortran E 格式)

---

## 文件说明

### Fortran 代码

- **`eqgs1d_export.f`**: 核心导出子程序
  - `SUBROUTINE EQGS1D_EXPORT`: 遍历 23 个剖面并写入 CSV
  - 使用 Fortran OPEN/WRITE/CLOSE 操作文件
  - 与 `EQGS1D` 使用相同的数据源

- **`eqgout.f`** (修改):
  - 添加 'SE' 命令处理
  - 调用 `EQGS1D_EXPORT`

- **`Makefile`** (修改):
  - 将 `eqgs1d_export.f` 添加到源文件列表

### Python 脚本

- **`read_eqgs1d_csv.py`**: 数据读取和分析工具
  - `load_profile(n)`: 读取第 n 个剖面
  - `load_all_profiles()`: 读取所有 23 个剖面
  - `plot_profile(n)`: 绘制第 n 个剖面
  - `plot_all_profiles()`: 绘制所有剖面（6x4 网格）
  - `export_to_excel()`: 导出到 Excel 文件
  - `summary_statistics()`: 显示统计摘要

### 文档

- **`eqgr1d_parameters.csv`**: EQGR1D 调用参数的元数据表
- **`gs_图像分析.md`**: GS 图像详细分析报告
- **`EQGS1D_EXPORT_README.md`**: 本文档

---

## 常见问题

### Q: 如何只导出部分数据？
A: 编辑 `eqgs1d_export.f`，注释掉不需要的部分。

### Q: 如何更改输出格式？
A: 修改 `eqgs1d_export.f` 中的 `WRITE` 语句格式。

### Q: Python 脚本需要哪些依赖？
A: 需要 `pandas`, `numpy`, `matplotlib`, `openpyxl`。安装：
```bash
pip install pandas numpy matplotlib openpyxl
```

### Q: 能否从 eqdata 二进制文件直接读取？
A: 可以，但需要了解 Fortran 二进制格式。建议使用本工具导出 CSV。

### Q: 导出的数据和 GS 图像中的数据一致吗？
A: 完全一致。两者使用相同的变量（来自 `eqcomq.inc`）。

---

## 技术细节

### 数据来源

所有数据来自 `eqcomq.inc` 中定义的公共变量，这些变量在 `EQCALQ` 子程序中计算：

- `PSIP(NR)`: 归一化极向磁通
- `PPS(NR)`: 压强剖面
- `QPS(NR)`: 安全因子
- 等等...

### 与 EQGR1D 的关系

`EQGS1D_EXPORT` 与 `EQGS1D` 使用相同的数据准备方式：
```fortran
! X 轴数据
DO NR=1,NRMAX
   GX(NR)=GUCLIP(PSIP(NR))
ENDDO

! Y 轴数据
DO NR=1,NRMAX
   GY(NR,1)=GUCLIP(PPS(NR))*1.E-6
ENDDO
```

因此导出的 CSV 数据与 GS 图像中显示的数据完全一致。

---

## 下一步

1. **扩展功能**: 添加更多剖面（如果需要）
2. **数据分析**: 使用 Python 脚本进行深入分析
3. **集成到工作流**: 将导出步骤集成到自动化脚本中

---

## 作者

自动生成工具 - 2026-01-14
