# `in/tr.CFEDR.in` 参数解释

基于 `doc/trcomm.tex` 和 `trparm.f90`，以下是 `in/tr.CFEDR.in` 中使用的参数的详细分类解释。

## 通用配置 (General Configuration)
- **`modelg=3`**: 设置几何模型。`MODELG=3` 表示模拟使用外部平衡输入（调用 `TREQIN`），通常映射到几何中心。
- **`KNAMEQ='../eq/in/gfile_efit_output'`**: 指定要读取的平衡文件（G-file/EFIT 输出）的路径。
- **`NSMAX=4`**: 设置模拟的热粒子种类数量（例如：电子、氘、氚、氦/杂质）。

## 等离子体分布 (初始条件)
这些参数定义了密度 ($n$) 和温度 ($T$) 的初始径向分布。
- **`PN`**: 磁轴上的初始数密度 ($10^{20} \text{m}^{-3}$)，对应每种粒子。
  - `0.7D0` (种类 1 - 电子)
  - `0.315D0` (种类 2)
  - `0.315D0` (种类 3)
  - `0.035D0` (种类 4)
- **`PNS`**: 等离子体边界上的初始数密度 ($10^{20} \text{m}^{-3}$)。
- **`PT`**: 轴上的初始温度 (keV)，对应每种粒子。
- **`PTS`**: 边界上的初始温度 (keV)。
- **`PROFN2`**: 密度分布的形状因子。
  - 分布通常建模为 $n(\rho) = (1 - \rho^{\text{PROFN1}})^{\text{PROFN2}}$。
- **`MDLNF=1`**: 模型标志，可能用于启用特定的快粒子或杂质模型。

## 加热与电流驱动 (Heating and Current Drive)
各种加热方法的参数。通常缩写含义如下：
- **`TOT`**: 总功率 (Total Power, MW)。
- **`R0`**: 沉积中心位置 (归一化半径)。
- **`RW`**: 沉积宽度 (Deposition width)。
- **`CD`**: 电流驱动效率或标志 (Current Drive)。
- **`NPR`**: 分布功率因子 (Profile power factor)。
- **`ENG`**: 束能量 (Beam Energy, keV)。

### 中性束注入 (NBI)
- **`PNBR0=0.D0`**: NBI 沉积中心（在轴上）。
- **`PNBRW=1.D0`**: NBI 沉积分布宽度（较宽）。
- **`PNBENG=1000.D0`**: 束能量 (1 MeV，ITER/CFEDR NNBI 的典型值)。
- **`PNBRTG`**: 切向半径或几何因子 (Tangency radius)。
- **`PNBTOT`**: NBI 总功率。

### 离子回旋共振加热 (ICRF)
- **`PICCD`**, **`PICR0`**, **`PICRW`**, **`PICNPR`**, **`PICTOT`**: 分别对应 ICRF 的电流驱动、沉积位置、宽度、形状因子和总功率。

### 电子回旋共振加热 (ECRF)
- **`PECCD`**, **`PECR0`**, **`PECRW`**, **`PECNPR`**, **`PECTOT`**: 分别对应 ECRF 的电流驱动、沉积位置、宽度、形状因子和总功率。

### 低杂波电流驱动 (LHRF)
- **`PLHCD`**, **`PLHR0`**, **`PLHRW`**, **`PLHNPR`**: 分别对应 LHRF 的电流驱动、沉积位置、宽度和形状因子。

## 时间与模拟控制 (Time and Simulation Control)
- **`DT=0.01D0`**: 时间步长（秒）。
- **`NTSTEP=1000`**: 每个主要循环/输出的时间步数。
- **`NTMAX`**: 运行的最大时间步数。脚本中多次更新此值以分段推进模拟 (1500 -> 2000 -> 6500)。

## 等离子体电流 (Plasma Current)
- **`RIPS`**: 初始等离子体电流 (MA)。
- **`RIPE`**: 最终等离子体电流 (MA)。

## 交互命令 (Interactive Commands)
文件中还包含用于交互式菜单的控制字符：
- **`r`**: 运行 (Run) - 开始或继续模拟。
- **`c`**: 继续? (Continue) - 可能用于更改参数后继续。
- **`g`**: 绘图 (Graph) - 进入绘图菜单。
- **`x`**: 退出/保存? (Exit) - 退出当前菜单。
- **`q`**: 退出程序 (Quit)。
