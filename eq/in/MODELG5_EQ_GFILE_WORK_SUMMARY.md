# MODELG=5 纯 EQ 路线工作总结（阶段文档）

## 1. 目标
主目标是：

1. 从 `MODELG=5` 读取 gfile 的压强/电流相关剖面（`p(psi), p'(psi), F(psi), q(psi)` 等）。
2. 在 **纯 `eq` 模块**内重算磁平衡（不依赖外部 Python GS 解作为主路径）。
3. 将重算结果与 gfile 平衡做系统对比。
4. 当前重点：搞清楚并对齐 `eq` 中 `F(psi)` 的计算链路与约定。

---

## 2. 已完成工作（对话内累计）

## 2.1 对比体系从“看 PSIRZ 图”升级为多指标
已将对比指标扩展为：

1. `q(psi)`
2. `p(psi)`
3. `F(psi)`
4. LCFS
5. 磁轴

并恢复/补强了 overlay 图和 summary 总图。

## 2.2 gfile 剖面导出与自动对比链路
已建立从 gfile 导出剖面到自动对比的流程：

1. gfile -> CSV（`p, p', F, FF', q` 等）
2. 与 `eqdata` 自动对比并汇总到大图
3. 后续继续叠加 namelist 的 `tot` 电流做交叉核验

## 2.3 关键坐标约定排查：`psi_n` vs `rho_tor`
已确认并反复验证：

1. `rho_tor = sqrt(normalized toroidal flux)` 需由 `q(psi_p)` 积分构造
2. 直接把均匀 `psi_n` 当 `rho_tor` 会引入系统性错配
3. 该问题会同时污染 `jphi`、`p'`、`FF'` 的比较解释

## 2.4 GS 残差与边界诊断
做过：

1. 矩形 Dirichlet 与 LCFS Dirichlet 版本
2. LCFS 标注与几何对比
3. GS 残差图诊断

结论：即使源项固定，边界离散、几何重构与坐标映射误差仍会导致明显差异。

## 2.5 纯 `eq` 路线的核心改造方向已明确
按你指定路线，主线固定：

1. `MDLEQF=9`
2. `QPSL` 不再依赖 `EQQPV` 初始近似，改由 gfile 的 `QQPS(PSIPS)` 插值
3. 由该 `q(psi_p)` 重构 `psi_t(psi_p)` 映射
4. 同步更新 `PSIPNV/PSITV/QPV/TTV` 与 `UPPSI/UFPSI/UQPSI` 链路，保证一致性

该方向的本质是先修复 `q -> psi_t` 映射，再看 `p/F` 末端行为。

## 2.6 当前绘图工具已支持 gfile 全量叠加
`in/plot_eqipqp_all_variables_big.py` 已增强为：

1. 不只叠加 `q/p/F`
2. 还从 gfile 反算并叠加 `DPPSI/RSV/AVIR2/VPV/AVRR2/PSITV/PSIPNV/dvdpsip/aveir2`
3. 第一幅 `F(psi)` 图已改成 `logscale`

用于直接定位 `eq` 与 gfile 在 `EQIPQP` 输入量层面的差异。

---

## 3. 当前仍存在的问题

1. `q(psi)` 在外侧（约 `psi>0.8`）出现“近常数化”。
2. `p(psi)` 或 `dp/dpsi_n` 在外侧（约 `psi>0.7`）出现“提前塌到零”。
3. `F(psi)` 与 gfile 仍有幅值/末端行为差异（不仅是简单符号问题）。
4. LCFS/axis 虽接近但未完全重合，提示几何与边界仍有偏差来源。

---

## 4. 现在正在做的事：`eq` 模块里 `F(psi)` 到底怎么来

## 4.1 `F(psi)` 在 `eq` 内有两条路径

1. **参数化路径（`MDLEQF<5`）**
   - `EQFPSI` 通过 `EQFUNC` 用参数化公式构造 `F`
2. **样条路径（`MDLEQF>=5`，当前主线）**
   - `EQFPSI` 从 `UFPSI` 样条取 `F(psi_t)`
   - 关键映射是 `psi_p -> psi_t`（`EQPSITN`）

## 4.2 为什么 `psi_t` 映射是关键
`EQFPSI` 的导数换算因子为：

- $FDN = QPV / PSITA$

这意味着：

1. `QPV/PSITA` 若被错误映射污染
2. `DPPSI/DFPSI` 都会被一起拉偏
3. `EQIPQP` 反推 `FIPV` 时误差会沿径向积分累积放大

## 4.3 `EQIPQP` 与 `F` 的关系
`EQIPQP` 用 `q, dp/dpsi, RSV, AVIR2, VPV, AVRR2` 从边界向轴反推 `FIPV`。
因此 `F` 对不上的常见根因不是单点误差，而是整条链路定义不一致：

1. `q` 映射（`psi_p <-> psi_t`）不一致
2. 几何量使用了不同定义（`EQCALV` vs `EQCALQ` 量混用）
3. 单位约定不一致（`F` vs `2*pi*F`）
4. 外侧剖面被样条/边界处理钳平

---

## 5. 下一步执行计划（围绕 F）

1. 固定单一约定：全链路统一使用 `F = 2*pi*R*Bphi`（并在图上明确标注）。
2. 在同一次运行中导出并对齐以下中间量：`PSIPNV, PSITV, QPV, UFPSI输入, FIPV输出`。
3. 在 `psi_n=0.7~1.0` 外侧单独做数值审计：
   - 样条节点是否被截断
   - `QPV` 是否被钳平
   - `FDN=QPV/PSITA` 是否异常变平
4. 对 `EQIPQP` 做逐环带差分检查（每个 `NRV` 的 `AL/BL/FL/YM/FIPV`），定位 `F` 偏差从哪一层开始出现。

---

## 6. 现用命令（当前版本）

1. 画全变量大图（含 gfile 叠加）：
```bash
conda run -n Fusion python3 in/plot_eqipqp_all_variables_big.py --prefix eqdata_modelg5_rebuild_mdleqf9_qmap
```

2. 若需显式指定 gfile 与符号约定：
```bash
conda run -n Fusion python3 in/plot_eqipqp_all_variables_big.py \
  --prefix eqdata_modelg5_rebuild_mdleqf9_qmap \
  --gfile in/g260206.20000_teq_0114 \
  --gfile-profile-csv in/g260206.20000_teq_0114_profiles_from_gfile.csv \
  --gfile-f-sign -1 \
  --gfile-pprime-sign 1
```

