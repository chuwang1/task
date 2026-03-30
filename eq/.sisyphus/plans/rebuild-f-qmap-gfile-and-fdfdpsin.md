# Rebuild F 流程切换 q(gfile) 并新增 F*dF/dpsi_n 对比面板

## TL;DR

> **Quick Summary**: 将 `rebuild_f.sh` 的 `F_rebuilt` 重建默认切换到 `q(gfile)`，同时保留 `q(psi)` 可见性（双输出共存），并在最终大图中用最小侵入方式替换 1 个面板为 `F*dF/dpsi_n` 对比。  
> **Deliverables**:
> - `rebuild_f.sh` 使用 `--use-gfile-qp` 且绘图读取本次重建产物
> - `in/rebuild_f_from_q_eqdata.py` 输出双 q 列（`q` + `q_gfile`/`q_used`）且消除 q 列一致性风险
> - `in/plot_eqipqp_all_variables_big.py` 新增 `F*dF/dpsi_n` 对比面板（替换现有 4x4 中 1 格）
>
> **Estimated Effort**: Short  
> **Parallel Execution**: YES - 2 waves  
> **Critical Path**: T1 → T6 → T9 → F1-F4

---

## Context

### Original Request
用户要求：
1) 运行 `rebuild_f.sh` 时，在 `F_rebuilt` 重建阶段将 `q(psi)`输出切换为 `q(gfile)`；
2) 最终输出图 `in/eqdata_modelg5_rebuild_mdleqf9_qmap_eqipqp_all_variables_big.png` 中增加一幅 `F*dF/dpsi_n` 对比图。

### Interview Summary
**Key Discussions**:
- q 策略：用户确认“双输出并存”（保留 q(psi) 可见性 + 使用 q(gfile)）
- 测试策略：不引入新的自动化测试框架；采用脚本执行+产物断言+证据留存
- 图布局：替换现有 4x4 中 1 个面板（最小侵入）
- 计算策略：`F_rebuilt` 默认强制走 `q(gfile)`

**Research Findings**:
- `q(gfile)`链路已存在（`export_gfile_profiles.py` + `_load_gfile_qp`）
- 现有风险：`rebuild_f_from_q_eqdata.py` 中当使用 `q_use` 重建时，输出列仍可能写入原 `q`（不一致）
- `plot_eqipqp_all_variables_big.py` 当前是 4x4 面板，panel 13 可作为最小替换位
- `rebuild_f.sh` 目前存在 rebuild 产物与 plot 读取前缀可能不一致风险，需显式连线

### Metis Review
**Identified Gaps (addressed in this plan)**:
- G1: rebuild 与 plot 读取 CSV 不一致风险 → T1/T6 显式指定 `--rebuild-csv`
- G2: q 列语义不一致风险 → T2 明确 dual columns + source 标识
- G3: `F*dF/dpsi_n` 比较语义模糊 → 本计划锁定为 `F_rebuilt` vs `F_gfile` 的同轴比较
- G4: 面板替换后诊断信息丢失风险 → Guardrail：仅替换 panel13，不扩散到其他 panel

---

## Work Objectives

### Core Objective
确保 `F_rebuilt` 主流程默认基于 `q(gfile)`重建，且输出与可视化对齐，并在最终大图中加入可执行验证的 `F*dF/dpsi_n` 对比结果。

### Concrete Deliverables
- `rebuild_f.sh`：重建命令启用 `--use-gfile-qp`，绘图命令显式读取本次重建 CSV
- `in/rebuild_f_from_q_eqdata.py`：输出中保留原 `q` 并新增 `q_gfile`/`q_used`（或等价双列）
- `in/plot_eqipqp_all_variables_big.py`：panel13 改为 `F*dF/dpsi_n` 对比面板

### Definition of Done
- [ ] `bash rebuild_f.sh` 完成且返回码 0
- [ ] 产物 CSV 同时含 q(psi) 与 q(gfile) 可见列
- [ ] 目标 PNG 成功生成并包含新面板

### Must Have
- 默认 `F_rebuilt` 使用 q(gfile)
- q 双输出共存（兼容对比）
- 面板替换方式实现（4x4 不扩展）

### Must NOT Have (Guardrails)
- 不新增测试框架（pytest/vitest 等）
- 不修改 panel 0-12、14-15 的语义与布局
- 不改动 `--use-gfile-qp` Python 参数默认值（仍由 shell 流程决定启用）
- 不做与本需求无关的重构（函数签名、模块拆分、依赖变更）

---

## Verification Strategy (MANDATORY)

> **ZERO HUMAN INTERVENTION** — 所有验收由执行 agent 通过命令和文件断言完成。

### Test Decision
- **Infrastructure exists**: NO（无 formal 单测框架）
- **Automated tests**: None（不新增框架）
- **Framework**: N/A

### QA Policy
每个任务都必须含 agent-executed QA 场景与证据文件，证据写入：`.sisyphus/evidence/task-{N}-{slug}.{ext}`。

- **CLI/Pipeline**: Bash + exit code + artifact existence
- **Data Consistency**: Python one-liner 读取 CSV 校验列与有限值
- **Plot Output**: 文件存在性 + 修改时间 + 可选像素维度/面板文本检查

---

## Execution Strategy

### Parallel Execution Waves

Wave 1 (Start Immediately — 实现主改动，5 并发):
├── Task 1: `rebuild_f.sh` 启用 q(gfile)并显式连线 rebuild_csv [quick]
├── Task 2: `rebuild_f_from_q_eqdata.py` 双 q 列输出与 source 标识 [unspecified-high]
├── Task 3: `rebuild_f_from_q_eqdata.py` 明确 `F_direct` 与 q source 语义 [unspecified-high]
├── Task 4: `plot_eqipqp_all_variables_big.py` 增加 `F*dF/dpsi_n` 计算段 [quick]
└── Task 5: `plot_eqipqp_all_variables_big.py` 替换 panel13 绘制逻辑 [visual-engineering]

Wave 2 (After Wave 1 — 集成与验证，4 并发):
├── Task 6: 前缀/路径契约验证（rebuild→plot 同一产物）[quick]
├── Task 7: 回归验证（不带 `--use-gfile-qp` 直调脚本路径不破坏）[unspecified-high]
├── Task 8: 端到端运行 `rebuild_f.sh` + 产物断言 [deep]
└── Task 9: 证据归档与验收报告生成 [writing]

Wave FINAL (After ALL tasks — independent review, 4 parallel):
├── Task F1: Plan compliance audit (oracle)
├── Task F2: Code quality review (unspecified-high)
├── Task F3: Real manual QA execution by agent scripts (unspecified-high)
└── Task F4: Scope fidelity check (deep)

Critical Path: 1 → 6 → 8 → F1-F4

### Dependency Matrix
- **1**: none → 6,8
- **2**: none → 6,8
- **3**: none → 7,8
- **4**: none → 5,8
- **5**: 4 → 8
- **6**: 1,2 → 8,9
- **7**: 3 → 9
- **8**: 1,2,5,6 → 9,F1-F4
- **9**: 6,7,8 → F1-F4

### Agent Dispatch Summary
- **Wave 1 (5)**: T1 quick, T2 unspecified-high, T3 unspecified-high, T4 quick, T5 visual-engineering
- **Wave 2 (4)**: T6 quick, T7 unspecified-high, T8 deep, T9 writing
- **Final (4)**: F1 oracle, F2 unspecified-high, F3 unspecified-high, F4 deep

---

## TODOs

- [ ] 1. 更新 `rebuild_f.sh`：默认启用 q(gfile) 且 plot 读取本次重建 CSV

  **What to do**:
  - 在调用 `in/rebuild_f_from_q_eqdata.py` 时加入 `--use-gfile-qp`
  - 在调用 `in/plot_eqipqp_all_variables_big.py` 时传入显式 `--rebuild-csv` 指向本次输出
  - 保持现有 `--prefix` 与 `--out` 命名兼容

  **Must NOT do**:
  - 不改 shell 里无关命令
  - 不改变输出主 PNG 路径

  **Recommended Agent Profile**:
  - **Category**: `quick`（单脚本参数连线改动）
  - **Skills**: `[]`

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 1
  - **Blocks**: 6,8
  - **Blocked By**: None

  **References**:
  - `rebuild_f.sh:4-10` - 重建命令与绘图命令入口
  - `in/plot_eqipqp_all_variables_big.py:507` - `--rebuild-csv` 默认解析逻辑

  **Acceptance Criteria**:
  - [ ] `bash rebuild_f.sh` 日志出现 `using gfile q/p' from:`
  - [ ] plot 实际读取 CSV 为本次 rebuild 输出路径（命令参数可见）

  **QA Scenarios**:
  ```
  Scenario: rebuild_f.sh 启用 gfile q 的 happy path
    Tool: Bash
    Preconditions: gfile profile CSV 存在
    Steps:
      1. 运行 `bash rebuild_f.sh > .sisyphus/evidence/task-1-run.log 2>&1`
      2. 搜索日志关键字 `using gfile q/p' from:`
      3. 断言退出码为 0
    Expected Result: 日志出现关键字且脚本成功退出
    Evidence: .sisyphus/evidence/task-1-run.log

  Scenario: 缺失 gfile CSV 的错误路径
    Tool: Bash
    Preconditions: 临时传入不存在的 gfile csv 参数
    Steps:
      1. 运行等价 python 命令并指向不存在文件
      2. 断言非0退出且错误信息包含 missing/not found
    Expected Result: 失败可解释，非 silent failure
    Evidence: .sisyphus/evidence/task-1-missing-gfile-error.log
  ```

  **Commit**: YES (group 1)

- [ ] 2. 修复/增强 `in/rebuild_f_from_q_eqdata.py` 的 q 输出一致性（双输出）

  **What to do**:
  - 输出 CSV 中保留原始 q(psi) 列（兼容）
  - 同时新增 `q_gfile`（或 `q_used` + `q_source`）体现重建实际来源
  - 确保当 `--use-gfile-qp` 时，`F_rebuilt` 对应 q 列可追溯

  **Must NOT do**:
  - 不删除 legacy `q` 列
  - 不重命名已有关键输出列导致下游崩溃

  **Recommended Agent Profile**:
  - **Category**: `unspecified-high`
  - **Skills**: `[]`

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 1
  - **Blocks**: 6,8
  - **Blocked By**: None

  **References**:
  - `in/rebuild_f_from_q_eqdata.py:259,279,396-401,480` - q / q_use 加载与输出关键行
  - `in/export_gfile_profiles.py:146` - gfile q 列名约定

  **Acceptance Criteria**:
  - [ ] 输出 CSV 表头包含双 q 可见字段
  - [ ] `F_rebuilt` 非 NaN 有效点数 > 50

  **QA Scenarios**:
  ```
  Scenario: 双 q 列存在与可读
    Tool: Bash
    Preconditions: 已生成 rebuild CSV
    Steps:
      1. `python3 -c "import pandas as pd; d=pd.read_csv('in/eqdata_modelg5_rebuild_mdleqf9_qmap_samegeom_f_from_q_eqipqp_stable.csv'); print(d.columns.tolist())" > .sisyphus/evidence/task-2-columns.txt`
      2. 断言输出包含 `q` 与 `q_gfile`(或 q_used)
    Expected Result: 双输出列存在
    Evidence: .sisyphus/evidence/task-2-columns.txt

  Scenario: q(gfile)列缺失时的失败
    Tool: Bash
    Preconditions: 临时使用旧CSV（不含新列）执行校验脚本
    Steps:
      1. 运行列断言脚本指向旧CSV
      2. 断言脚本失败并提示缺列
    Expected Result: 校验脚本能明确拦截缺列
    Evidence: .sisyphus/evidence/task-2-missing-column-error.log
  ```

  **Commit**: YES (group 1)

- [ ] 3. 明确 `F_direct` 与 q source 语义一致性

  **What to do**:
  - 审核并固定 `F_direct` 是否继续使用原 q(psi) 或切到 q_used
  - 若保持原语义，增加字段/注释/元信息避免误读
  - 若切换，确保与 `F_rebuilt` 的 q source 一致

  **Must NOT do**:
  - 不在未标注情况下混用不同 q source

  **Recommended Agent Profile**:
  - **Category**: `unspecified-high`
  - **Skills**: `[]`

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 1
  - **Blocks**: 7,8
  - **Blocked By**: None

  **References**:
  - `in/rebuild_f_from_q_eqdata.py:266` - `F_direct` 当前依赖 q
  - `in/rebuild_f_from_q_eqdata.py:483-487` - `F_direct`/`dF_direct` 输出列

  **Acceptance Criteria**:
  - [ ] CSV 中 `F_direct` 的 q source 语义可机读（字段或列名可判断）

  **QA Scenarios**:
  ```
  Scenario: F_direct 语义可追溯 happy path
    Tool: Bash
    Preconditions: 新CSV已生成
    Steps:
      1. Python脚本读取 F_direct 与 q-source 指示列
      2. 打印语义标签到证据文件
      3. 断言标签不为空且仅一个来源
    Expected Result: F_direct source 明确无歧义
    Evidence: .sisyphus/evidence/task-3-fdirect-source.txt

  Scenario: source 指示缺失
    Tool: Bash
    Preconditions: 人为指向无 source 标识版本
    Steps:
      1. 运行语义断言脚本
      2. 断言失败并输出缺失原因
    Expected Result: 失败可解释
    Evidence: .sisyphus/evidence/task-3-source-missing-error.log
  ```

  **Commit**: YES (group 1)

- [ ] 4. 在 plot 脚本中新增 `F*dF/dpsi_n` 计算段

  **What to do**:
  - 基于现有 F 曲线（`F_rebuilt` 与 `F_gfile`）按 `psi_n` 计算导数
  - 使用安全梯度策略（边界点可控，不引入 NaN 扩散）

  **Must NOT do**:
  - 不新增外部依赖
  - 不改变已有数据加载流程

  **Recommended Agent Profile**:
  - **Category**: `quick`
  - **Skills**: `[]`

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 1
  - **Blocks**: 5,8
  - **Blocked By**: None

  **References**:
  - `in/plot_eqipqp_all_variables_big.py:586-611` - gfile/eq 数据已齐备的位置
  - `in/compare_psi_3way.py:217` - `f * gradient(f, psi_n)` 模式参考

  **Acceptance Criteria**:
  - [ ] 计算段在数据可用时生成两条可绘制数组
  - [ ] 数组长度与 x 轴一致

  **QA Scenarios**:
  ```
  Scenario: 计算数组生成成功
    Tool: Bash
    Preconditions: plot脚本可运行
    Steps:
      1. 运行plot脚本并输出日志
      2. 日志/断言脚本检查新数组长度>0且finite比例>95%
    Expected Result: 可绘制数据就绪
    Evidence: .sisyphus/evidence/task-4-compute.log

  Scenario: gfile数据缺失时容错
    Tool: Bash
    Preconditions: 传入不存在gfile profile路径
    Steps:
      1. 运行plot脚本
      2. 断言脚本失败信息可解释（而非崩溃栈不明）
    Expected Result: 错误处理明确
    Evidence: .sisyphus/evidence/task-4-missing-gfile-error.log
  ```

  **Commit**: YES (group 2)

- [ ] 5. 替换 panel13 为 `F*dF/dpsi_n` 对比面板

  **What to do**:
  - 保持 4x4 布局不变，仅替换 `axs[13]` 内容
  - 绘制 `F_rebuilt*dF_rebuilt/dpsi_n` 与 `F_gfile*dF_gfile/dpsi_n`
  - 设置标题、坐标、legend、grid 与现有风格一致

  **Must NOT do**:
  - 不改动其他 panel 编号和语义

  **Recommended Agent Profile**:
  - **Category**: `visual-engineering`
  - **Skills**: [`frontend-ui-ux`]

  **Parallelization**:
  - **Can Run In Parallel**: NO
  - **Parallel Group**: Wave 1 (after T4)
  - **Blocks**: 8
  - **Blocked By**: 4

  **References**:
  - `in/plot_eqipqp_all_variables_big.py:776-788` - 现 panel13 位置
  - `in/plot_eqipqp_all_variables_big.py:620-774` - 统一绘图样式参考

  **Acceptance Criteria**:
  - [ ] 生成 PNG 含新 panel13 标题 `F*dF/dpsi_n`（或约定同义标题）
  - [ ] 图像成功保存至目标路径

  **QA Scenarios**:
  ```
  Scenario: 新面板出图 happy path
    Tool: Bash
    Preconditions: 数据与脚本依赖完整
    Steps:
      1. 运行plot脚本生成目标PNG
      2. 校验文件存在且mtime为本次执行
      3. 若脚本支持日志标题输出，断言包含 panel13 新标题
    Expected Result: 新面板成功落图
    Evidence: .sisyphus/evidence/task-5-plot-success.log

  Scenario: x轴长度不一致错误场景
    Tool: Bash
    Preconditions: 构造长度不一致输入(旧CSV)
    Steps:
      1. 运行plot脚本
      2. 断言错误提示指向维度不匹配
    Expected Result: 失败可解释
    Evidence: .sisyphus/evidence/task-5-dim-mismatch-error.log
  ```

  **Commit**: YES (group 2)

- [ ] 6. 固化 rebuild→plot 的文件契约（前缀/路径一致性）

  **What to do**:
  - 确认 `rebuild_f.sh` 产物路径与 plot 读取路径一致
  - 若存在 `_samegeom` 与基础前缀差异，显式参数覆盖而非隐式默认

  **Must NOT do**:
  - 不改变文件命名策略导致历史脚本失效

  **Recommended Agent Profile**:
  - **Category**: `quick`
  - **Skills**: []

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 2
  - **Blocks**: 8,9
  - **Blocked By**: 1,2

  **References**:
  - `rebuild_f.sh:4-10` - 当前命令链路
  - `in/plot_eqipqp_all_variables_big.py:507,515` - rebuild_csv 默认拼接逻辑

  **Acceptance Criteria**:
  - [ ] plot 命令中 `--rebuild-csv` 指向刚生成文件
  - [ ] 无“读取旧 run 文件”行为

  **QA Scenarios**:
  ```
  Scenario: 路径契约正确
    Tool: Bash
    Preconditions: rebuild_f.sh 已更新
    Steps:
      1. 执行脚本并记录命令日志
      2. 解析日志中的 rebuild 输出路径与 plot 输入路径
      3. 断言两者相等
    Expected Result: 同一CSV链路
    Evidence: .sisyphus/evidence/task-6-path-contract.log

  Scenario: 默认路径回退导致错读
    Tool: Bash
    Preconditions: 去掉 --rebuild-csv 进行对照执行
    Steps:
      1. 执行对照命令
      2. 断言读取路径与本次产物不一致并触发告警
    Expected Result: 问题可复现且可被检测
    Evidence: .sisyphus/evidence/task-6-default-mismatch.log
  ```

  **Commit**: YES (group 3)

- [ ] 7. 做无 gfile flag 路径回归验证（脚本级）

  **What to do**:
  - 直接调用 `in/rebuild_f_from_q_eqdata.py`（不带 `--use-gfile-qp`）
  - 对比关键列（`F_rebuilt` 统计量/有限值数）与基线一致性

  **Must NOT do**:
  - 不要求人工看图

  **Recommended Agent Profile**:
  - **Category**: `unspecified-high`
  - **Skills**: []

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 2
  - **Blocks**: 9
  - **Blocked By**: 3

  **References**:
  - `in/rebuild_f_from_q_eqdata.py:202-233` - CLI 参数入口
  - 现有基线 CSV：`in/eqdata_modelg5_rebuild_mdleqf9_qmap_f_from_q_eqipqp.csv`

  **Acceptance Criteria**:
  - [ ] 不带 gfile flag 运行成功
  - [ ] 回归统计差异在预设阈值内

  **QA Scenarios**:
  ```
  Scenario: 默认路径回归通过
    Tool: Bash
    Preconditions: 基线CSV可读
    Steps:
      1. 运行无flag命令生成对照CSV
      2. Python脚本比较 F_rebuilt 的 finite count/均值
      3. 断言差异低于阈值
    Expected Result: 默认行为未破坏
    Evidence: .sisyphus/evidence/task-7-regression.txt

  Scenario: 基线文件缺失
    Tool: Bash
    Preconditions: 临时指向不存在基线
    Steps:
      1. 运行比较脚本
      2. 断言报错说明基线缺失
    Expected Result: 失败可解释
    Evidence: .sisyphus/evidence/task-7-missing-baseline-error.log
  ```

  **Commit**: NO

- [ ] 8. 端到端执行 `rebuild_f.sh` 并做产物断言

  **What to do**:
  - 跑完整 pipeline
  - 断言 CSV 列完整、F 列有效、PNG 新鲜生成
  - 断言目标图路径：`in/eqdata_modelg5_rebuild_mdleqf9_qmap_eqipqp_all_variables_big.png`

  **Must NOT do**:
  - 不以“手动看着像对”作为唯一验收

  **Recommended Agent Profile**:
  - **Category**: `deep`
  - **Skills**: []

  **Parallelization**:
  - **Can Run In Parallel**: NO
  - **Parallel Group**: Wave 2 (integration)
  - **Blocks**: 9,F1-F4
  - **Blocked By**: 1,2,5,6

  **References**:
  - `rebuild_f.sh`
  - 目标 PNG 路径与 rebuild CSV 路径约定

  **Acceptance Criteria**:
  - [ ] `bash rebuild_f.sh` 返回码 0
  - [ ] CSV 含双 q 输出字段
  - [ ] 目标 PNG 存在且 mtime 为当前执行窗口

  **QA Scenarios**:
  ```
  Scenario: 全链路 happy path
    Tool: Bash
    Preconditions: conda环境和输入文件就绪
    Steps:
      1. `bash rebuild_f.sh > .sisyphus/evidence/task-8-e2e.log 2>&1`
      2. Python脚本断言CSV列与F_rebuilt finite count
      3. 断言目标PNG存在且mtime<120秒
    Expected Result: 全链路通过
    Evidence: .sisyphus/evidence/task-8-e2e.log

  Scenario: 关键输入缺失
    Tool: Bash
    Preconditions: 临时断开某输入CSV路径
    Steps:
      1. 执行pipeline
      2. 断言非0退出且报错文件名明确
    Expected Result: 失败快速暴露
    Evidence: .sisyphus/evidence/task-8-input-missing-error.log
  ```

  **Commit**: NO

- [ ] 9. 证据归档与机器可读验收报告

  **What to do**:
  - 汇总每个任务证据到 `.sisyphus/evidence/`
  - 生成验收摘要（通过/失败、失败原因、复现命令）

  **Must NOT do**:
  - 不省略失败场景证据

  **Recommended Agent Profile**:
  - **Category**: `writing`
  - **Skills**: []

  **Parallelization**:
  - **Can Run In Parallel**: NO
  - **Parallel Group**: Wave 2 收尾
  - **Blocks**: F1-F4
  - **Blocked By**: 6,7,8

  **References**:
  - `.sisyphus/evidence/task-*`

  **Acceptance Criteria**:
  - [ ] 验收报告覆盖 1-9 全部任务
  - [ ] 每个任务至少 1 个 happy + 1 个 error 证据

  **QA Scenarios**:
  ```
  Scenario: 证据完整性检查
    Tool: Bash
    Preconditions: 前置任务已生成证据
    Steps:
      1. 扫描 .sisyphus/evidence/task-*.log|txt
      2. 统计每任务 happy/error 证据数量
      3. 输出报告
    Expected Result: 每任务证据齐全
    Evidence: .sisyphus/evidence/task-9-evidence-audit.txt

  Scenario: 缺失证据场景
    Tool: Bash
    Preconditions: 人为缺失一个任务证据
    Steps:
      1. 运行审计脚本
      2. 断言报告明确指出缺失任务
    Expected Result: 审计失败且可定位
    Evidence: .sisyphus/evidence/task-9-missing-evidence-error.txt
  ```

  **Commit**: NO

---

## Final Verification Wave (MANDATORY)

- [ ] F1. **Plan Compliance Audit** — `oracle`
  - 校验 Must Have / Must NOT Have 与任务证据一致性

- [ ] F2. **Code Quality Review** — `unspecified-high`
  - 运行静态检查与脚本检查，排查临时调试残留

- [ ] F3. **Real Manual QA (agent-executed)** — `unspecified-high`
  - 按任务 QA 场景逐条执行，归档最终证据

- [ ] F4. **Scope Fidelity Check** — `deep`
  - 核对变更只覆盖计划范围，无额外扩散

---

## Commit Strategy

- **Group 1**: `fix(rebuild): force gfile-q flow and dual q outputs`
  - Files: `rebuild_f.sh`, `in/rebuild_f_from_q_eqdata.py`
- **Group 2**: `feat(plot): add F*dF/dpsi_n panel by replacing panel13`
  - Files: `in/plot_eqipqp_all_variables_big.py`
- **Group 3**: `chore(pipeline): harden rebuild-to-plot csv contract`
  - Files: `rebuild_f.sh` (if needed extra alignment)

---

## Success Criteria

### Verification Commands
```bash
bash rebuild_f.sh
python3 -c "import pandas as pd; d=pd.read_csv('in/eqdata_modelg5_rebuild_mdleqf9_qmap_samegeom_f_from_q_eqipqp_stable.csv'); assert 'q' in d.columns; assert ('q_gfile' in d.columns) or ('q_used' in d.columns); import numpy as np; assert np.isfinite(d['F_rebuilt']).sum()>50"
python3 -c "import os,time; f='in/eqdata_modelg5_rebuild_mdleqf9_qmap_eqipqp_all_variables_big.png'; assert os.path.exists(f); assert time.time()-os.path.getmtime(f)<120"
```

### Final Checklist
- [ ] All "Must Have" present
- [ ] All "Must NOT Have" absent
- [ ] Script-driven QA scenarios passed with evidence
