# TASK/EQ 模块文档（源码导向）

## 1. 模块定位
`EQ` 是 TASK 系统中的 Grad-Shafranov 平衡求解模块，支持：
- 固定边界/自由边界平衡求解
- 从外部平衡文件读取（TASK/EQ 二进制、EQDSK）
- 平衡后处理（通量面量、1D/2D 数据、图形输出、指标导出）

主程序入口在 `eqmain.f`，交互流程由 `eqmenu.f` 控制。

## 2. 主要功能

### 2.1 平衡求解
- 固定边界主流程：`EQCALC`（`eqcalc.f`）
- 自由边界求解：`EQCALX`（`eqcalx.f`）
- 磁轴/边界处理：`EQAXIS` 等（`eqsub.f`）

### 2.2 平衡读取
- TASK/EQ 二进制：`MODELG=3/9` -> `EQRTSK`（`eqfile.f`）
- EQDSK(gfile)：`MODELG=5` -> `EQDSKR`（`eq-eqdsk.f`）

### 2.3 剖面与物理量
- 剖面函数计算：`EQPPSI/EQFPSI/EQQPSI/EQJPSI/EQTPSI/EQOPSI`（`eqfunc.f`）
- 通量坐标与面平均量：`EQCALQ`（`eqcalq.f`）
- 真空区与外推：`EQCALQV`（`eqcalq.f`）

### 2.4 输出与可视化
- 保存平衡文件：`EQSAVE`（`eqfile.f`）
- 图形输出：`EQGOUT`（`eqgout.f`）
- CSV 导出：
  - 1D 导出：`eqgs1d_export.f`
  - 2D 导出：`eqgs2d_export.f`
- Metric 输出：`EQMETRIC`（`eqfile.f`）

## 3. 输入体系

## 3.1 交互命令（`eqmenu.f`）
- `R`：运行求解（`EQCALC`）
- `C`：在当前解基础上继续迭代（快速改 `PP0,PP1,PP2,PJ0,PJ1,PJ2,RIP,HM`）
- `P`：输入 `&EQ` 参数
- `V`：查看当前参数
- `I`：从参数文件读 `&EQ`
- `G`：图形输出
- `S`：保存平衡文件（`KNAMEQ`）
- `L`：读取平衡（TASK/EQ）
- `K`：读取 EQDSK
- `F`：输出 metric
- `H`：读取 ripple 轮廓并处理

## 3.2 `&EQ` 参数（`eqinit.f` 的 `NAMELIST /EQ/`）
参数总入口在 `EQNLIN`，可分组理解：

1. 几何与全局约束
- `RR, RA, RB, RKAP, RDLT, BB, Q0, QA, RIP, FRBIN`

2. 模型开关
- `MODELG`：几何/文件类型（2=托卡马克求解，5=EQDSK读入，3/9=TASK/EQ读入）
- `MODELQ`：`q` 模型
- `MDLEQF`：剖面形式（解析或样条）
- `MDLEQA, MDLEQC, MDLEQX, MDLEQV`

3. 剖面参数
- 压强：`PP0,PP1,PP2,PROFP0,PROFP1,PROFP2`
- 电流：`PJ0,PJ1,PJ2,PROFJ0,PROFJ1,PROFJ2`
- F 函数：`FF0,FF1,FF2,PROFF0,PROFF1,PROFF2`
- 温度/转动：`PT* / PROFTP* / PTSEQ`, `PV* / PROFV*`
- 形状参数：`PROFR0,PROFR1,PROFR2,RHOITB`

4. 网格与收敛
- 网格：`NSGMAX,NTGMAX,NUGMAX,NRGMAX,NZGMAX,NPSMAX,NRMAX,NTHMAX,NSUMAX,NRVMAX,NTVMAX`
- 收敛：`EPSEQ,NLPMAX,EPSNW,DELNW,NLPNW`
- 计算域：`RGMIN,RGMAX,ZGMIN,ZGMAX,ZLIMP,ZLIMM`

5. 边界/线圈（自由边界）
- `PSIB, NPFCMAX, RIPFC, RPFC, ZPFC, WPFC`

## 3.3 外部剖面输入能力
- `eq` 交互模式本身不直接支持“CSV 文件路径 -> 直接读入剖面”。
- 代码层提供接口 `TREQEX`（`treqin.f`），可传入数组 `PRHO/HJRHO/VTRHO/TRHO` 并走样条剖面求解（会设置 `MDLEQF=7/8`）。
- 若需要“从文件读特定剖面”，建议新增一个驱动程序：文件 -> 数组 -> `TREQEX`。

## 4. 输出体系

## 4.1 主要文件输出
1. 平衡二进制（`KNAMEQ`）
- 由 `EQSAVE` 写出，包含网格、`PSIRZ`、1D 剖面、轴位置、几何参数、通量坐标量等。

2. 图形输出（`.gs`）
- 由图形菜单触发。

3. CSV 导出
- 2D：`eqgs2d_01_PSIRZ_grid.csv` 等
- 1D：`eqgs1d_*.csv`

4. Metric
- `eq_metric.dat`

## 4.2 典型导出字段（TASK/EQ 二进制记录）
- 标量：`RR, BB, RIP, RAXIS, ZAXIS, PSITA, PSIPA, PSI0, ...`
- 网格：`RG, ZG, PSIRZ`
- 剖面：`PSIPS, PPPS, TTPS, TEPS, OMPS`
- 面平均：`PSIPNV, PSIPV, PSITV, QPV, TTV`

## 5. 关键源代码文件索引

## 5.1 核心流程
- `eqmain.f`：主程序入口，初始化与菜单调度
- `eqmenu.f`：交互命令和流程分发
- `eqinit.f`：默认参数、`&EQ` 读入、参数检查/显示

## 5.2 求解与物理
- `eqcalc.f`：固定边界求解主流程
- `eqcalx.f`：自由边界求解
- `eqcalq.f`：后处理（q、面平均、真空区处理）
- `eqfunc.f`：剖面函数及导数
- `eqsub.f`：磁轴/边界子程序

## 5.3 文件读写与格式
- `eqfile.f`：保存/读取 TASK/EQ 平衡、metric 输出
- `eq-eqdsk.f`：EQDSK 读取
- `equread.f90`：读取支持模块

## 5.4 图形与导出
- `eqgout.f`, `eqg2d.f`, `eqg3d.f`, `eqgsub.f`
- `eqgs1d_export.f`, `eqgs2d_export.f`

## 5.5 接口与扩展
- `treqin.f`：TR <-> EQ 剖面/平衡接口（`TREQEX`）
- `equintf.f`, `eqintf.f`：接口层
- `eqcom*.inc`：全局公共变量（common blocks）

## 6. 典型工作流

### 6.1 计算平衡（托卡马克）
1. `MODELG=2`，在 `&EQ` 设置几何与剖面参数
2. 菜单 `R` 求解，`C` 局部调参继续
3. `S` 保存 `KNAMEQ`
4. `G`/CSV 导出并做对比

### 6.2 从 gfile 对比/接轨
1. `MODELG=5` 读入 EQDSK
2. 导出 `PSIRZ/q/p/F` 等指标
3. 与计算解比较：`q(psi), p(psi), F(psi), LCFS, axis`

## 7. 常用调参建议（面向复现实验平衡）
- 先锁定：`RR, RA, RKAP, RDLT, BB, RIP`
- 再调：`PP* / PROFP*`（压强）、`PJ* / PROFJ*`（电流）
- 自由边界下根据目标形状调 `MDLEQX` 与边界/线圈参数
- 提升分辨率时同步调整 `NRGMAX/NZGMAX` 与收敛参数

## 8. 备注
- 该文档基于当前源码结构（Fortran 固定格式）总结。
- 若后续引入“文件化剖面输入”驱动，建议在本文件追加一节“Profile IO API”。
