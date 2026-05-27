# TR 输出 PRC 的 PyMak ONETWO 复算验证

本验证用于检查 Fortran TR 中 `MDLPR=2` 的 ONETWO PRC 实现是否与
PyMak 的 `calculate_cyclotron_radiation_onetwo` 一致。

验证不在 `/tmp` 下运行。运行目录固定为：

```text
validation/onetwo_prc/tr_output_3step/
```

## 运行命令

从 `trx.260120` 目录执行：

```bash
validation/onetwo_prc/run_tr_output_3step.sh
```

脚本会执行以下步骤：

1. 在 `validation/onetwo_prc/tr_output_3step/` 下准备本地运行目录。
2. 链接 TR 运行所需文件：`tr2`、`density_input0114_rho.csv`、
   `omfit_prl_for_tr.csv`、`chi_div_grad2_drdrho.dat`。
3. 使用 `tr_onetwo_3step.in` 完整运行 3 个时间步。
4. 用 PyMak ONETWO 公式从 TR 输出 CSV 的 `nE`、`TE` 重新计算 PRC。
5. 将对比结果写入 `validation/onetwo_prc/results/`。

## TR 输入设置

关键输入文件：

```text
validation/onetwo_prc/tr_output_3step/tr_onetwo_3step.in
```

关键设置：

```text
MDLPR  = 2
REFRAD = 0.8D0
ntmax  = 3
dt     = 0.0001
NRMAX  = 200
```

## 对比数据来源

对比脚本会自动寻找 `TRGRR1/TRGRR2` 直接输出的当前数组，而不是
`GVR` 历史缓存输出。本次运行使用：

```text
tr_data_017.csv  nE
tr_data_019.csv  TE
tr_data_022.csv  PRC
```

这些文件位于：

```text
validation/onetwo_prc/tr_output_3step/
```

该目录下的原始 TR 输出 CSV 和日志由脚本重新生成，并被 `.gitignore`
忽略；提交中保留可复现输入、脚本、汇总 CSV 和图像。

## 几何口径

`MODELG=5` 会读取 gfile，并更新 TR 实际使用的几何/磁场口径。
因此 PyMak 复算不能直接使用输入卡中的 `BB=6.0`、`RR=8.03`、`RA=2.72`。

本验证按 TR/eq 读入后的 EQDSK 口径使用：

```text
RR = 8.0304174650 m
RA = 2.7215121350 m
BB = -5.9772733123 T
```

ONETWO PRC 公式中使用 `abs(BB)`，所以磁场符号不影响结果。

## 结果

输出文件：

```text
validation/onetwo_prc/results/tr_output_3step_prc_comparison.csv
validation/onetwo_prc/results/tr_output_3step_prc_comparison.png
```

本次 3 步验证结果：

```text
points                 = 200
max_abs_err_mw_m3      = 3.0207979673e-08
max_rel_err_masked     = 1.1142070094e-06
mean_rel_err_masked    = 3.3564965195e-07
core_tr                = 3.3908670000e-02 MW/m^3
core_pymak             = 3.3908667637e-02 MW/m^3
core_rel_err           = 6.9696558758e-08
```

这里的 `max_rel_err_masked` 只统计 `abs(TR PRC) > 1e-8 MW/m^3` 的网格点，
避免边界极小值放大相对误差。

## 结论

在使用 TR 实际 gfile 几何后，Fortran `MDLPR=2` 的 ONETWO PRC 输出与
PyMak ONETWO 复算一致，误差达到约 `1e-6` 相对量级。此前若使用输入卡
原始 `BB=6.0` 复算，会出现约 `1%` 偏差；该偏差来自几何/磁场口径不一致，
不是 ONETWO 公式移植误差。
