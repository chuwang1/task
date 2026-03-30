# EQ `MDLEQF=9` 路线与当前问题说明

## 1. 目标与结论（先看）

- 目标：在 `MODELG=5`（从 gfile 读入）下，用 `MDLEQF=9`（给定 `P,q` 的样条 profile）重建并求解平衡。
- 目前结论：
  - **几何量/安全因子 -> F 的后处理公式链路是可行的**（Python 侧已验证）。
  - **前向 EQ 迭代链路在 `EQLOOP` 内失稳**，表现为 `EQMAGS` 失败、`EQAXIS` 牛顿失败，最终 `PSI/DELPSI` 出现 NaN。
  - 失败不是单一由“初值是否来自 gfile”导致；但初始 $\psi(R,Z)$ 与目标平衡的偏差幅度会显著影响稳定性。

---

## 2. 当前 `MDLEQF=9` 路线（Fortran / m5）

入口：`modelg5_rebuild.f`

主流程：

1. `EQINIT` + 参数设置：`NLPMAX=400`, `NTVMAX=1024`, `NTGMAX=64`, `NPRINT=2`（`modelg5_rebuild.f:77-80`）
2. `MODELG=5` + `EQLOAD` 从 gfile 读入（`modelg5_rebuild.f:82-85`）
3. `EQMESH -> EQPSIN -> EQDEFB`（`modelg5_rebuild.f:104-109`）
4. `EQSET_MODELG5_TREQ`：用 gfile 的 `PPPS/TTPS/QQPS` 重建内部样条与 `PSITA` 映射（`modelg5_rebuild.f:117-118`）
5. `EQLOOP` 迭代求解（`modelg5_rebuild.f:126-129`）
6. 成功后 `EQTORZ/EQCALP/EQCALQ/EQSAVE`（`modelg5_rebuild.f:134-150`）

边界来源：

- 默认 `IUSELCFS=1`，优先走 gfile LCFS（`RSU/ZSU`）映射到 `RHOM/RHOG`（`eqcalc.f:189-195`）。
- 可通过 `USE_GFILE_LCFS=0` 强制回退 `EQFBND` 参数化边界（`eqcalc.f:197-198`，`run_m5_all.sh:10,18,20`）。

---

## 3. 已定位并修复/加入的关键改动

### 3.1 失败不再“假成功”

- `m5` 关键步骤失败时返回非零退出（`STOP 1`），避免后续脚本继续读旧结果（`modelg5_rebuild.f:145-147`）。
- `run_m5_all.sh` 在运行前删除旧前缀文件并检查新输出是否生成（`run_m5_all.sh:19-22`）。

### 3.2 轴点前处理修正

- `MODELG=5` 下不再在 `EQPSIN/EQPSIR` 强制 `RAXIS=RR, ZAXIS=0`（`eqcalc.f:62-64`, `eqcalc.f:121-123`）。

### 3.3 诊断增强

- 在 `modelg5_rebuild.f` 增加阶段轴点差日志：`after EQLOAD/EQPSIN/EQDEFB/EQSET_TREQ/EQLOOP/EQCALQ`（`modelg5_rebuild.f:97,112,115,119,128,139`）。
- 在 `EQLOOP` 增加 NaN/非法和 `SUM0/SUM1` 检查，返回 `IERR=901/902/903`（`eqcalc.f:464-483`）。
- 在 `EQDERV` 增加 `|grad psi|` 过小/NaN 告警（`eqsub.f:228-234`）。

### 3.4 软固定轴尝试

- `EQAXIS` 牛顿失败时回退上一轮轴，且限制每轮轴步长（当前 `DAXLIM=0.25 m`）（`eqsub.f:30`, `eqsub.f:41-49`）。

---

## 4. 当前失败签名（无 FPE trap 诊断运行）

在临时去掉 `-ffpe-trap` 后可得到完整轨迹：

1. `EQLOOP` 前阶段轴点与 gfile 一致。
2. `EQLOOP` 内出现多次 `XX EQMAGS: NOT ENOUGH N`。
3. 随后 `XX EQAXIS: NEWTN ERROR`。
4. 再后 `PSI/DELPSI` 出现 NaN，被新检查拦截为 `IERR=901`。

即：**先磁面追踪失败 -> 再轴求解失败 -> 再全场变量坏化（NaN）**。

---

## 5. 与 Python 求解器的关键差异（为什么 Python 更稳）

文件：`in/solve_gs_from_gfile.py`

- Python 走规则 `R-Z` 网格上的 Picard + Jacobi，不使用 `EQMAGS` 场线闭合追踪。
- LCFS 用多边形 Dirichlet + outside 冻结，数值保护较强（clip/ghost-cell/mask）。
- 在部分测试中，使用 `--init-psi-csv in/eqdata0114_eq_psirz_rz.csv`（初值不是 gfile 原始 `psirz`）仍可收敛；但当初值与目标差别过大时，Python 侧同样可能出现数值问题。
- 已加 `--recompute-axis` 做动态轴测试，说明 Fortran 主失稳点仍在 `EQLOOP` 闭环（`EQMAGS/EQAXIS` 链路）；但“初值偏差过大”会放大失稳风险，而非可忽略因素。

结论：Python 的稳定不代表 Fortran EQ 闭环（含轴牛顿+磁面追踪）必然稳定。

---

## 6. 目前最可能的根因组合

1. LCFS 单值化到 `RHOM(theta)` 的近边界误差（`NTGMAX=32` 时尤重；提高到 64 后虽改善但不根治）。
2. 外层 profile（尤其 q 尾部）与边界/网格耦合后使 `EQMAGS` 更易失效。
3. `EQLOOP` 对局部坏态容错不足，坏值可扩散到 `PSI/DELPSI`。

---

## 7. 建议下一步（按优先级）

1. 在 `EQCALV/EQMAGS` 失败时中止当前迭代更新，避免坏值写回主场。
2. 继续收紧轴回退策略（例如步长阈值从 `0.25m` 降到 `0.10m`，并在越界时直接拒绝更新）。
3. 对 q 尾部正则化做参数扫描（`PSIQREG0/QTAILFAC`），寻找稳定与物理保真折中点。
4. 做一次“首个失败磁面”定位输出（`NR, RINIT, PSIP`），确认是哪一圈先崩。

---

## 8. 相关脚本/文件索引

- 主程序：`modelg5_rebuild.f`
- 边界定义：`eqcalc.f`（`EQDEFB`, `EQSET_RHOB_FROM_RSU`, `EQFBND`）
- 迭代与RHS：`eqcalc.f`（`EQLOOP`, `EQRHSV`）
- 轴与磁面追踪：`eqsub.f`（`EQAXIS`, `EQMAGS`, `EQDERV`）
- 批处理：`run_m5_all.sh`
- Python 对照：`in/solve_gs_from_gfile.py`
