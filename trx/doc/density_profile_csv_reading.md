# 外部密度 CSV 读取功能说明

**日期**: 2026-03-10  
**修改来源**: 从 `~/TASK/latest/task/trx/` 移植

## 功能概述

新增从外部 CSV 文件读取多粒子种密度剖面的功能，替代原有的解析函数生成密度。涉及三个参数：

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `model_prof` | INTEGER | 0 | 密度剖面模型选择。设为 **12** 启用 CSV 多粒子种读取 |
| `knam_prof` | CHARACTER | `'prof.data'` | CSV 文件路径 |
| `model_nevolve` | INTEGER | 0 | 设为 **1** 时，每个时间步将密度重置为 CSV 值（密度不演化） |

## CSV 文件格式

文件第一行为表头，之后每行 5 列：

```
rho, ne, nD, nT, nHe
0.01, 1.86e+20, 8.0e+19, 8.0e+19, 5.6e+18
0.02, 1.85e+20, 7.9e+19, 7.9e+19, 5.5e+18
...
```

- **rho**: 归一化小半径 (0~1)
- **ne**: 电子密度 (m⁻³)
- **nD**: 氘密度 (m⁻³)
- **nT**: 氚密度 (m⁻³)
- **nHe**: 氦密度 (m⁻³)

> 注意：密度单位为 m⁻³（绝对值），不是 10²⁰ m⁻³。

## 使用示例

在输入文件中设置：

```fortran
&TR
   model_prof=12
   knam_prof='density_input0114_rho.csv'
   model_nevolve=1
/
```

## 工作原理

### 初始化阶段 (`model_prof=12`)

1. 读取 CSV 文件，跳过表头
2. 对 4 个粒子种（e, D, T, He）分别建立三次样条插值
3. 在每个径向网格点上用样条插值设置密度 `RN(nr,ns)`
4. 超过 4 种的粒子按电子密度等比例缩放
5. 边界密度 `PNS` 取自最外层网格点的密度值

### 时间演化阶段 (`model_nevolve=1`)

每个时间步执行 `tr_exec` 之后，调用 `tr_reset_density`：
- 从初始化时存储的样条数据重新插值，覆盖当前密度
- 效果：密度固定为 CSV 给定的剖面，不随输运方程演化

## 修改的文件

| 文件 | 修改内容 |
|------|---------|
| `trcomm.f90` | 添加 `model_nevolve` 变量声明 |
| `trinit.f90` | 添加 `model_nevolve=0` 默认值 |
| `trparm.f90` | 将 `model_nevolve` 加入 NAMELIST /TR/ |
| `trprof.f90` | 添加 CASE(12) CSV 读取逻辑、样条插值、边界密度设置、`tr_reset_density` 子程序 |
| `trloop.f90` | 在时间循环中调用 `tr_reset_density` |
| `trview.f90` | 添加 `model_nevolve` 参数显示 |

所有原文件已备份为 `.bak` 后缀。

## model_prof 取值汇总

| 值 | 说明 |
|----|------|
| 0 | 默认：用 `PROFN1/PROFN2` 解析函数生成密度 |
| 11 | 从文件读取单粒子种（电子）密度，其余按比例缩放 |
| **12** | **（新增）从 CSV 读取 4 个粒子种密度** |
| 41,42 | 从 total profile 文件读取 |
