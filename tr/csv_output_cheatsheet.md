# TR 模拟输出 CSV 变量完整速查表

此表基于 `tr.CFEDR.in` 标准执行序列生成。

## 1. 时间演化 (Time Evolution)
**X轴**: 时间 `t` (s) | **用途**: 查看全局参数随时间的变化

| 文件 | 变量 | 物理含义 | 单位 |
|:-----|:-----|:---------|:-----|
| 001 | TE, TD, \<TE\>, \<TD\> | 电子/D中心温度, 平均温度 | keV |
| 002 | IP, IOH, INB, IRF, IBS, IPARA | 总电流, 欧姆, NB驱动, RF驱动, 自举, 平行电流 | MA |
| 003 | PIN, POH, PNB, PRF, PNF, PEX | 总输入/欧姆/NB/RF/聚变/外部功率 | MW |
| 004 | POUT, PCX, PIE, PRL, PCON | 损失/电荷交换/离子-电子交换/辐射/热传导 | MW |
| 005 | QF | 聚变增益 Q值 | - |
| 006 | TEAV, TDAV, TTAV, TAAV | 各粒子平均温度 (E/D/T/He) | keV |
| 007 | TE0, TD0, TT0, TA0 | 各粒子中心温度 | keV |
| 008 | NE0, ND0, NT0, NA0, \<NE\> | 各粒子中心密度, 平均电子密度 | 10²⁰/m³ |
| 009 | WF, WB, WI, WE | 快粒子/束离子/离子/电子储能 | MJ |
| 010 | BETAa, BETA0, BETAN | 平均/中心/归一化 Beta | % |
| 011 | TAUE1, TAUE2, TAUE89 | 能量约束时间 (不同定义) | s |
| 012 | BETAPa, BETAP0 | 极向 Beta (平均/中心) | - |
| 013 | VLOOP | 环电压 | V |
| 014 | Q(0) | 轴心安全因子 | - |
| 015 | ALI | 内电感 li | - |
| 016 | RQ1 | q=1面半径 | m |

---

## 2. 径向剖面快照 (Radial Profiles)
**X轴**: 半径 `r` (m) | **用途**: 查看最终时刻的物理量分布

| 文件 | 变量 | 物理含义 | 单位 |
|:-----|:-----|:---------|:-----|
| 017 | NE, ND, NT, NA | 电子/D/T/He密度 | 10²⁰/m³ |
| 018 | NB, NF | 束离子/聚变产物密度 | 10²⁰/m³ |
| 019 | TE, TD, TT, TA | 各粒子温度 | keV |
| 020 | PE, PD, PT, PA | 各粒子压力 | MPa |
| 021 | POH, PNB, PNF, PR, PRF | 加热功率密度 (OH/NB/聚变/辐射/RF) | MW/m³ |
| 022 | POH, PRB, PRC, PRL, PCX, PIE | 损失功率密度 (轫致/回旋/线辐射/CX/i-e交换) | MW/m³ |
| 023 | AKE, AKNCE, AKDWE | 电子热扩散 (总/新经典/湍流) | m²/s |
| 024 | AKD, AKNCD, AKDWD | 离子热扩散 (总/新经典/湍流) | m²/s |
| 025 | QP | 安全因子 q | - |
| 026 | EZOH | 感应电场 | V/m |
| 027 | JTOT, JOH, JNB, JRF, JBS | 电流密度 (总/欧姆/NB/RF/自举) | MA/m² |
| 028 | LOG:ETA | 电阻率对数 | Ω·m |
| 029 | G, s, alpha | 度规因子/磁剪切/压力梯度参数 | - |
| 030 | LOG:AKE, AKI | 扩散系数对数 | m²/s |
| 031 | AKE, AKNCE, AKDWE | 电子热扩散详情 | m²/s |
| 032 | AKD, AKNCD, AKDWD | 离子热扩散详情 | m²/s |

---

## 3. 剖面时间演化 (Profile Evolution)
**X轴**: 半径 `r` (m) | **列**: 不同时刻快照 | **用途**: 观察剖面随时间的变化

| 文件 | 变量 | 物理含义 | 单位 |
|:-----|:-----|:---------|:-----|
| 033 | NE | 电子密度演化 | 10²⁰/m³ |
| 034 | ND | D密度演化 | 10²⁰/m³ |
| 035 | TE | 电子温度演化 | keV |
| 036 | TD | D温度演化 | keV |
| 037 | NT | T密度演化 | 10²⁰/m³ |
| 038 | NA | He (α粒子) 密度演化 | 10²⁰/m³ |
| 039 | TT | T温度演化 | keV |
| 040 | TA | He (α粒子) 温度演化 | keV |
| 041 | QP | 安全因子 q 剖面演化 | - |
| 042 | AJ | 总电流密度演化 | MA/m² |
| 043 | EZ | 电场演化 | V/m |
| 044 | AJOH | 欧姆电流演化 | MA/m² |
| 045 | AJNB+AJRF | 驱动电流演化 (NB+RF) | MA/m² |
| 046 | AJBS | 自举电流演化 | MA/m² |
| 047 | PIN | 输入功率密度演化 | MW/m³ |
| 048 | POH | 欧姆加热功率演化 | MW/m³ |

---

## 变量命名规则

| 前缀/后缀 | 含义 |
|:---------|:-----|
| E | 电子 (Electron) |
| D | 氘 (Deuterium) |
| T | 氚 (Tritium) |
| A | α粒子/氦 (Alpha/Helium) |
| B | 束离子 (Beam ion) |
| F | 聚变产物 (Fusion product) |
| 0 | 中心值 (On-axis) |
| AV/\<\> | 平均值 (Volume average) |
| NC | 新经典 (Neoclassical) |
| DW | 湍流/漂移波 (Drift-wave) |
| OH | 欧姆 (Ohmic) |
| NB | 中性束 (Neutral Beam) |
| RF | 射频 (Radio Frequency) |
| BS | 自举 (Bootstrap) |
