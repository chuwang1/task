# CDBM 输运模型技术手册

## 1. 概述

**电流扩散气球模 (Current Diffusivity Ballooning Mode, CDBM)** 是 TR 代码中用于托卡马克等离子体模拟的主要反常输运模型。该模型基于电流密度梯度和压力梯度驱动的微撕裂和气球不稳定性理论，计算湍流热扩散系数 $\chi_e$、$\chi_i$。

模型实现在 `trcoef.f90` 的子程序 `TRCFDW` 中，通过输入参数 `MDLKAI` 选择。

### MDLKAI 范围

| MDLKAI | 说明 |
|--------|------|
| 30–39 | 内联 CDBM，使用 $|\alpha|^{3/2}$ 标度 |
| 40–50 | 内联 CDBM，使用 $|\alpha|$ 标度 + 电子热速度修正 |
| 130–139 | 调用外部子程序 `tr_cdbm`（推荐用于生产计算） |

---

## 2. 物理背景

### 2.1 CDBM 核心公式

CDBM 模型的热扩散系数表达式为：

$$\chi = C_K \cdot F(s,\,\alpha,\,\kappa_q) \cdot |\alpha|^{3/2} \cdot \delta_e^2 \cdot \frac{V_A}{q R}$$

各符号含义：

| 符号 | 定义 | 单位 |
|------|------|------|
| $C_K$ | 数值系数（电子用 `CK0`，离子用 `CK1`） | 无量纲 |
| $F(s,\,\alpha,\,\kappa_q)$ | 结构函数（形状因子） | 无量纲 |
| $\alpha$ | 归一化压力梯度 | 无量纲 |
| $\delta_e$ | 无碰撞电子趋肤深度 | m |
| $V_A$ | Alfvén 速度 | m/s |
| $q$ | 安全因子 | 无量纲 |
| $R$ | 大半径 | m |

### 2.2 关键中间量

**Alfvén 速度：**

$$V_A = \frac{B}{\sqrt{\mu_0 \, n_e \, m_i}}$$

其中有效离子质量：

$$m_i = \frac{m_D \, n_D + m_T \, n_T + m_{He} \, n_{He}}{n_D + n_T + n_{He}}$$

对应 Fortran 代码 (`trcoef.f90:386-388`)：

```fortran
PNI = ANDX + ANT + ANA
AMI = (AMD*ANDX + AMT*ANT + AMA*ANA) / PNI
VA  = SQRT(BB**2 / (RMU0 * ANE * 1.D20 * AMI))
```

**无碰撞趋肤深度：**

$$\delta_e^2 = \frac{c^2}{\omega_{pe}^2} = \frac{c^2 \, m_e \, \varepsilon_0}{n_e \, e^2}$$

对应代码 (`trcoef.f90:611-612`)：

```fortran
WPE2   = ANE * 1.D20 * AEE**2 / (AME * EPS0)
DELTA2 = VC**2 / WPE2
```

**归一化压力梯度 (Shafranov $\alpha$)：**

$$\alpha = -\frac{2\,\mu_0\,q^2\,R}{B^2}\,\frac{dp}{dr}$$

其中 $dp/dr$ 是总压力（离子 + 电子 + 快粒子）梯度，单位为 $[10^{20}\,\mathrm{m}^{-3}\cdot\mathrm{keV}/\mathrm{m}]$，乘以 $10^{20}\cdot e_\mathrm{keV}$ 转换为 $[\mathrm{Pa/m}]$。

对应代码 (`trcoef.f90:392`)：

```fortran
ALPHA(NR) = -2.D0 * RMU0 * QL**2 * RR / BB**2 * (DPP * 1.D20 * RKEV)
```

**磁剪切：**

$$s = \frac{\rho}{q}\,\frac{dq}{dr}$$

使用 3 点数值导数从 $q$ 剖面计算。

对应代码 (`trcoef.f90:319-325, 390`)：

```fortran
! 内部点：中心差分
DQ = (QP(NR+1) - QP(NR-1)) / (2.D0 * DR)
! NR=1：前向差分
DQ = (4.D0*QP(2) - 3.D0*QP(1) - QP(3)) / (2.D0*DR)
! NR=NRMAX：后向差分
DQ = (3.D0*QP(NRMAX) - 4.D0*QP(NRMAX-1) + QP(NRMAX-2)) / (2.D0*DR)
! 磁剪切
S(NR) = RHOG(NR) / QL * DQ
```

**磁曲率：**

$$\kappa_q = -\varepsilon\left(1 - \frac{1}{q^2}\right)$$

其中 $\varepsilon = r/R$ 是局部逆纵横比。

对应代码 (`trcoef.f90:401`)：

```fortran
RKCV(NR) = -EPS * (1.D0 - 1.D0/(QL*QL))
```

---

## 3. 结构函数（形状因子）

CDBM 模型提供多个结构函数，由 `MDLKAI` 的末位数字选择。

### 3.1 TRCOFS($s$, $\alpha$, $\kappa_q$) — 标准形状因子

用于 `MDLKAI = 31, 32, 33, 34, 41, 42, 43, 44`。

**算法：**

定义偏移剪切 $\sigma$：

$$\sigma = \begin{cases} s - \alpha & \text{当 } \alpha \geq 0 \quad (\text{曲率分支使用 } \kappa_q > 0) \\ \alpha - s & \text{当 } \alpha < 0 \quad (\text{曲率分支使用 } \kappa_q < 0) \end{cases}$$

计算 $F_1$（气球模项）：

$$F_1 = \begin{cases} \displaystyle\frac{1 + 9\sqrt{2}\,\sigma^{5/2}}{\sqrt{2}\,(1 - 2\sigma + 3\sigma^2 + 2\sigma^3)} & \text{当 } \sigma \geq 0 \\[12pt] \displaystyle\frac{1}{\sqrt{2\,(1 - 2\sigma)\,(1 - 2\sigma + 3\sigma^2)}} & \text{当 } \sigma < 0 \end{cases}$$

计算 $F_2$（曲率驱动项）：

$$F_2 = \begin{cases} \displaystyle\frac{\kappa_q^{3/2}}{s^2} & \text{当 } \alpha \geq 0 \text{ 且 } \kappa_q > 0 \\[8pt] \displaystyle\frac{(-\kappa_q)^{3/2}}{s^2} & \text{当 } \alpha < 0 \text{ 且 } \kappa_q < 0 \\[8pt] 0 & \text{其他情况} \end{cases}$$

**最终结果：**

$$\mathrm{TRCOFS} = \max(F_1,\; F_2)$$

> **重要提示：** $F_2$ 中的 $1/s^2$ 项在 $s \to 0$ 时产生奇异性，这是磁轴附近"$\chi$ 尖峰"的主要来源。

对应代码 (`trcoef.f90:1907-1943`)：

```fortran
FUNCTION TRCOFS(S, ALPHA, RKCV)
  IF(ALPHA.GE.0.D0) THEN
     SA = S - ALPHA
     IF(SA.GE.0.D0) THEN
        FS1 = (1.D0 + 9.0D0*SQRT(2.D0)*SA**2.5D0) &
             /(SQRT(2.D0)*(1.D0 - 2.D0*SA + 3.D0*SA*SA + 2.0D0*SA**3))
     ELSE
        FS1 = 1.D0/SQRT(2.D0*(1.D0-2.D0*SA)*(1.D0-2.D0*SA+3.D0*SA*SA))
     ENDIF
     IF(RKCV.GT.0.D0) THEN
        FS2 = SQRT(RKCV)**3 / (S*S)
     ELSE
        FS2 = 0.D0
     ENDIF
  ELSE
     SA = ALPHA - S
     ! ... (对称结构，kappa_q < 0 时使用 sqrt(-kappa_q))
  ENDIF
  TRCOFS = MAX(FS1, FS2)
END FUNCTION
```

### 3.2 TRCOFSS($s$, $\alpha$) — 简化形状因子

用于 `MDLKAI = 35, 36, 37, 38`。

$$\sigma = s - \alpha, \qquad F = \frac{2\,\sigma^2}{1 + \frac{2}{9}\,|\sigma|^{5/2}}$$

该形式避免了 $1/s^2$ 奇异性。

对应代码 (`trcoef.f90:1991-2002`)：

```fortran
FUNCTION TRCOFSS(S, ALPHA)
  SA = S - ALPHA
  FS1 = 2.D0*SA**2 / (1.D0 + (2.D0/9.D0)*SQRT(ABS(SA))**5)
  TRCOFSS = FS1
END FUNCTION
```

### 3.3 TRCOFSX($s$, $\alpha$, $\kappa_q$, $\varepsilon_a$) — 扩展形状因子

用于 `MDLKAI = 39, 49`。

与 `TRCOFS` 类似，但增加了：

- **Shafranov 位移修正：** $\sigma = s - \left(1 - \frac{2\alpha}{1 + 6\alpha}\right)\alpha$
- **延伸度修正：** $F_1$ 乘以 $(1 + \kappa_q)^{5/2}$
- **曲率项修正：** 使用 $\displaystyle\frac{(\kappa_q/\varepsilon_a)^{3/2}}{s^2}$

### 3.4 TRCOFT($s$, $\alpha$, $\kappa_q$, $\varepsilon_a$) — 备选形状因子

用于 `MDLKAI = 50`。与 `TRCOFS` 相同，但曲率项使用：

$$F_2 = \frac{(\kappa_q/\varepsilon_a)^{3/2}}{s^2}$$

---

## 4. 内联 CDBM 模型 (MDLKAI = 30–50)

### 4.1 CASE(30:39) — 标准 CDBM

**通用公式：**

$$\chi_e = C_{K0} \cdot F(s,\,\alpha,\,\kappa_q) \cdot |\alpha|^{3/2} \cdot \delta_e^2 \cdot \frac{V_A}{q\,R}$$

$$\chi_i = C_{K1} \cdot F(s,\,\alpha,\,\kappa_q) \cdot |\alpha|^{3/2} \cdot \delta_e^2 \cdot \frac{V_A}{q\,R}$$

**子模型一览：**

| MDLKAI | 形状因子 $F$ | ExB 抑制 | $F$ 中含 $\alpha$ |
|--------|-------------|----------|-------------------|
| 30 | $\displaystyle\frac{1}{1.7 + \sqrt{6}\,s}$ | 否 | 否 |
| 31 | $\mathrm{TRCOFS}(s,\;\mathrm{CALF}\cdot\alpha,\;\kappa_q)$ | 否 | 是 |
| 32 | $\mathrm{TRCOFS}(s,\;\mathrm{CALF}\cdot\alpha,\;\kappa_q) \cdot f_{E\times B}$ | 是 | 是 |
| 33 | $\mathrm{TRCOFS}(s,\;0,\;\kappa_q)$ | 否 | 否 |
| 34 | $\mathrm{TRCOFS}(s,\;0,\;\kappa_q) \cdot f_{E\times B}$ | 是 | 否 |
| 35 | $\mathrm{TRCOFSS}(s,\;\mathrm{CALF}\cdot\alpha)$ | 否 | 是 |
| 36 | $\mathrm{TRCOFSS}(s,\;\mathrm{CALF}\cdot\alpha) \cdot f_{E\times B}$ | 是 | 是 |
| 37 | $\mathrm{TRCOFSS}(s,\;0)$ | 否 | 否 |
| 38 | $\mathrm{TRCOFSS}(s,\;0) \cdot f_{E\times B}$ | 是 | 否 |
| 39 | $\mathrm{TRCOFSX}(s,\;\mathrm{CALF}\cdot\alpha,\;\kappa_q,\;a/R)$ | 否 | 是 |

**CDBM05 延伸度修正**（当 `MDLCD05 != 0` 时，适用于 MDLKAI = 31, 32）：

$$F \;\longrightarrow\; F \cdot \left(\frac{2\sqrt{\kappa}}{1 + \kappa^2}\right)^{3/2}$$

其中 $\kappa = \mathrm{RKPRHO}(NR)$ 是来自平衡的局部延伸度。

对应代码 (`trcoef.f90:645-646`)：

```fortran
IF(MDLCD05.NE.0) &
     FS = FS * (2.D0*SQRT(RKPRHO(NR)) / (1.D0+RKPRHO(NR)**2))**1.5D0
```

### 4.2 CASE(40:50) — 漂移波气球模

使用 $|\alpha|$ 代替 $|\alpha|^{3/2}$，并增加电子热速度因子：

$$\chi = C_K \cdot F \cdot |\alpha| \cdot \delta_e^2 \cdot \frac{V_A}{q\,R} \cdot \frac{v_{Te}}{V_A}$$

其中电子热速度 $v_{Te} = \sqrt{T_e \cdot e_\mathrm{keV} / m_e}$。

MDLKAI = 50 还额外包含抗磁稳定化效应：

$$\chi_e = C_{K0} \cdot F \cdot |\alpha| \cdot \delta_e^2 \cdot \frac{V_A}{q\,R} \cdot \frac{1}{\Lambda(\lambda)\,(1 + \omega_{*T}^2)}$$

$$\chi_i = C_{K1} \cdot F \cdot |\alpha| \cdot \delta_e^2 \cdot \frac{V_A}{q\,R} \cdot \frac{1}{1 + \omega_{*T}^2}$$

其中 $\Lambda$ 是修正 Bessel 函数，$\omega_{*T}$ 是归一化抗磁频率。

---

## 5. 外部 CDBM 子程序 (MDLKAI = 130–139)

### 5.1 概述

对于 `MDLKAI = 130–139`，$\chi$ 计算委托给外部子程序 `tr_cdbm`：

```fortran
MODEL = MDLKAI - 130
CALL tr_cdbm(BB, RR, RS, RKAPL, QL, SHEARL, PNEL, RHONI, DPDRL,
             DVEXBDRL, CALF, CKAP, CEXB, MODEL, chi_cdbm)
AKDWEL = (CK0/12.0) * chi_cdbm
AKDWIL = (CK1/12.0) * chi_cdbm
```

**子模型：**

| MDLKAI | MODEL | 说明 |
|--------|-------|------|
| 130 | 0 | CDBM 原始版 |
| 131 | 1 | CDBM05（含延伸度修正） |
| 132 | 2 | CDBM + 弱 $E\times B$ 剪切 |
| 133 | 3 | CDBM05 + 弱 $E\times B$ 剪切 |
| 134 | 4 | CDBM + 强 $E\times B$ 剪切 |

### 5.2 $E\times B$ 剪切计算

在调用 `tr_cdbm` 之前，计算 $E\times B$ 剪切抑制因子：

$$S_L = s^2 + 0.01$$

$$\omega_{E1} = -\frac{q\,R}{S_L \cdot V_A}\,\frac{dV_{E\times B}}{dr}$$

$$C_{E\times B} = C_\mathrm{WEB} \cdot F_{E\times B}(|\omega_{E1}|,\;s,\;\alpha)$$

对应代码 (`trcoef.f90:908-911`)：

```fortran
SL  = (S(NR)**2 + 0.1D0**2)
WE1 = -QL*RR / (SL*VA) * DVE
RG1 = CWEB * FEXB(ABS(WE1), S(NR), ALPHA(NR))
cexb = RG1
```

### 5.3 `tr_cdbm` 的输入输出

**输入：**

| 变量 | 含义 | 单位 |
|------|------|------|
| $B$ (`BB`) | 环向磁场 | T |
| $R$ (`RR`) | 大半径 | m |
| $r_s$ (`RS`) | 局部小半径 $= a \cdot \rho$ | m |
| $\kappa$ (`RKAPL`) | 局部延伸度 | 无量纲 |
| $q$ (`QL`) | 局部安全因子 | 无量纲 |
| $s$ (`SHEARL`) | 磁剪切（可能经过平滑） | 无量纲 |
| $n_e$ (`PNEL`) | 电子密度 $= n_e \times 10^{20}$ | m$^{-3}$ |
| $\rho_i$ (`RHONI`) | 离子质量密度 | kg/m$^3$ |
| $dp/dr$ (`DPDRL`) | 压力梯度 | Pa/m |
| $dV_{E\times B}/dr$ (`DVEXBDRL`) | $E\times B$ 速度梯度 | 1/s |
| `CALF` | $\alpha$ 缩放因子 | 无量纲 |
| `CKAP` | 延伸度因子 ($= 1.0$) | 无量纲 |
| $C_{E\times B}$ (`CEXB`) | 来自 $F_{E\times B}$ 的抑制因子 | 无量纲 |
| `MODEL` | 子模型索引 (0–4) | — |

**输出：**

| 变量 | 含义 | 单位 |
|------|------|------|
| $\chi_\mathrm{cdbm}$ | 原始 CDBM 热扩散系数 | m$^2$/s |

原始 $\chi_\mathrm{cdbm}$ 随后经过缩放：

$$\chi_e = \frac{C_{K0}}{12}\,\chi_\mathrm{cdbm}, \qquad \chi_i = \frac{C_{K1}}{12}\,\chi_\mathrm{cdbm}$$

---

## 6. $E\times B$ 剪切抑制函数：FEXB

函数 $F_{E\times B}(x,\,s,\,\alpha)$ 返回一个抑制因子，取值范围 $[0,\,1]$：

$$f_{E\times B} = \exp\!\left(-\beta \cdot x^{\gamma}\right)$$

其中 $x = |\omega_{E\times B}|$ 是归一化 $E\times B$ 剪切率。

**$\beta$ 的计算：**

$$\alpha_L = \max(|\alpha|,\;0.001)$$

$$\beta = \frac{0.5}{\alpha_L^{0.602}} \cdot \frac{13.018 - 22.289\,s + 17.018\,s^2}{1 - 0.278\,s + 1.429\,s^2}$$

**$\gamma$ 的计算：**

$$\gamma = \begin{cases} \displaystyle\frac{1}{1.1\,\sqrt{1 - s - 2s^2 - 3s^3}} + 0.75 & \text{当 } s < 0 \\[12pt] \displaystyle\frac{1 - 0.5\,s}{1.1 - 2s + As^2 + 4s^3} + 0.75 & \text{当 } s \geq 0 \end{cases}$$

其中 $A = -\frac{10}{3}\,\alpha + \frac{16}{3}$。

> **注意：** $\beta$ 和 $\gamma$ 都依赖于 $s$ 的**符号**，而不仅仅是其大小。这意味着 $s > 0$ 和 $s < 0$ 即使绝对值相同，$E\times B$ 抑制效果也不同。

对应代码 (`trcoef.f90:2083-2120`)：

```fortran
FUNCTION FEXB(X, S, ALPHA)
  IF(ABS(ALPHA).LT.1.D-3) THEN
     ALPHAL = 1.D-3
  ELSE
     ALPHAL = ABS(ALPHA)
  ENDIF
  BETA = 0.5D0*ALPHAL**(-0.602D0) * (13.018D0 - 22.28915D0*S + 17.018D0*S**2) &
       / (1.D0 - 0.277584D0*S + 1.42913D0*S**2)
  ! ... gamma 计算 ...
  FEXB = EXP(-BETA * XG)
END FUNCTION
```

---

## 7. 磁剪切平滑与数值稳定化 (`model_cdbm_smooth`)

### 7.1 控制参数

```
model_cdbm_smooth = 0   ! 原始 CDBM（不做任何修改）
model_cdbm_smooth = 1   ! 5 点平滑 + 剪切下限 + chi 限幅（默认值）
```

### 7.2 磁剪切的 5 点 boxcar 平滑

在主径向循环之前预计算：

$$\hat{s}_\mathrm{HM}(j) = \begin{cases} s(1) & j = 1 \\ \frac{1}{3}\bigl[s(1) + s(2) + s(3)\bigr] & j = 2 \\ \frac{1}{5}\displaystyle\sum_{k=j-2}^{j+2} s(k) & j = 3, \ldots, N-2 \\ \frac{1}{3}\bigl[s(N-2) + s(N-1) + s(N)\bigr] & j = N-1 \\ s(N) & j = N \end{cases}$$

当 `model_cdbm_smooth = 0` 时：$\hat{s}_\mathrm{HM} = s$（不平滑）。

**物理动机：** 磁剪切 $s$ 由安全因子 $q$ 的导数计算而来，在网格上容易出现数值振荡。这些振荡通过形状因子 $\mathrm{TRCOFS}$ 中的 $1/s^2$ 项被放大，导致 $\chi$ 在 $r/a \approx 0.2$–$0.3$ 处出现尖峰。5 点平滑消除了网格尺度的噪声，同时保留了大尺度结构。

对应代码 (`trcoef.f90:150-163`)：

```fortran
IF(model_cdbm_smooth.GE.1) THEN
   S_HM(1) = S(1)
   S_HM(2) = (S(1)+S(2)+S(3))/3.D0
   DO NR=3,NRMAX-2
      S_HM(NR) = (S(NR-2)+S(NR-1)+S(NR)+S(NR+1)+S(NR+2))/5.D0
   ENDDO
   S_HM(NRMAX-1) = (S(NRMAX-2)+S(NRMAX-1)+S(NRMAX))/3.D0
   S_HM(NRMAX) = S(NRMAX)
ELSE
   S_HM(1:NRMAX) = S(1:NRMAX)
ENDIF
```

### 7.3 剪切下限

对于 `MDLKAI = 31`（当 `model_cdbm_smooth >= 1` 时）：

$$s_\mathrm{eff} = \max(s,\;0.5)$$

对于 `MDLKAI = 130:139`（当 `model_cdbm_smooth >= 1` 时）：

$$s_\mathrm{eff} = \max(\hat{s}_\mathrm{HM},\;0.5)$$

> **注意：** 强制 $s \geq 0.5$ 会禁止低剪切和负剪切区域的物理效应。在反剪切位形中，$s < 0$ 是有物理意义的。对于精确模拟，可能需要使用更小的正则化值。

### 7.4 $\chi$ 上下限（仅 MDLKAI = 31）

当 `model_cdbm_smooth >= 1` 时：

$$0.1 \leq \chi_e \leq 50 \quad [\mathrm{m}^2/\mathrm{s}]$$

$$0.1 \leq \chi_i \leq 50 \quad [\mathrm{m}^2/\mathrm{s}]$$

对应代码 (`trcoef.f90:649-655`)：

```fortran
IF(model_cdbm_smooth.GE.1) THEN
   AKDWEL = MIN(AKDWEL, 50.D0)
   AKDWIL = MIN(AKDWIL, 50.D0)
   AKDWEL = MAX(AKDWEL, 0.1D0)
   AKDWIL = MAX(AKDWIL, 0.1D0)
ENDIF
```

---

## 8. 数值实现细节

### 8.1 坐标系与网格

TR 代码使用归一化小半径 $\rho = r/a$ 作为径向坐标，其中 $a$ (`RA`) 为等离子体小半径。网格为**等间距**：

$$\Delta\rho = \frac{1}{N_\mathrm{RMAX}}, \qquad \Delta r = a \cdot \Delta\rho$$

代码中有两套径向网格：

- **整数网格** $\rho_M(NR) = (NR - \tfrac{1}{2})\,\Delta\rho$（密度、温度定义在此）
- **半网格** $\rho_G(NR) = NR \cdot \Delta\rho$（输运系数、梯度定义在此）

所有输运系数在**半网格** $\rho_G$ 上计算。密度和温度取相邻整数网格点的平均值：

$$n_e\big|_{\rho_G(NR)} = \frac{1}{2}\bigl[n_e(\rho_M(NR+1)) + n_e(\rho_M(NR))\bigr]$$

$$T_e\big|_{\rho_G(NR)} = \frac{1}{2}\bigl[T_e(\rho_M(NR+1)) + T_e(\rho_M(NR))\bigr]$$

在边界处 ($NR = N_\mathrm{RMAX}$)，使用边界值 $n_s^\mathrm{surf}$、$T_s^\mathrm{surf}$，梯度使用三点外推 (`DERIV3P`)。

### 8.2 导数的坐标约定

代码中的导数操作混合使用了 $\rho$ 坐标和物理半径 $r$ 坐标。关键转换因子为：

$$\mathrm{DRL} = \frac{1}{\Delta\rho \cdot a} = \frac{1}{\Delta r}$$

对应代码 (`trcoef.f90:169`)：

```fortran
DRL = 1.D0 / (DR * RA)
```

因此代码中的差分结果含义如下：

| 代码表达式 | 数学含义 | 坐标 |
|-----------|---------|------|
| `(RPP - RPM) * DRL` | $dp/dr$ | 对 $r$ |
| `(QP(NR+1) - QP(NR-1)) / (2*DR)` | $dq/d\rho$ | 对 $\rho$ |
| `DERIV3(NR, RHOG, VEXB, ...)` | $dV_{E\times B}/d\rho$ | 对 $\rho$ |

### 8.3 压力梯度

总压力包括电子、所有离子物种和快粒子（单位 $[10^{20}\,\mathrm{m}^{-3}\cdot\mathrm{keV}]$）：

$$P_{NR+1} = \sum_\mathrm{ions} n_s T_s + n_e T_e + P_\mathrm{beam} + P_\mathrm{add}$$

压力梯度对物理半径 $r$ 求导：

$$\frac{dp}{dr} = \frac{P_{NR+1} - P_{NR}}{\Delta\rho \cdot a} = (P_{NR+1} - P_{NR}) \cdot \mathrm{DRL}$$

### 8.4 安全因子导数与磁剪切

安全因子的导数**对 $\rho$ 求导**（注意不是对 $r$）。

内部点使用三点中心差分：

$$\frac{dq}{d\rho}\bigg|_{NR} = \frac{q_{NR+1} - q_{NR-1}}{2\,\Delta\rho}$$

边界使用前向/后向差分：

$$\frac{dq}{d\rho}\bigg|_{NR=1} = \frac{4q_2 - 3q_1 - q_3}{2\,\Delta\rho}$$

$$\frac{dq}{d\rho}\bigg|_{NR=N} = \frac{3q_N - 4q_{N-1} + q_{N-2}}{2\,\Delta\rho}$$

磁剪切的定义为：

$$s = \frac{\rho}{q}\,\frac{dq}{d\rho}$$

这与通常物理定义 $s = (r/q)\,dq/dr$ 是等价的，因为：

$$\frac{\rho}{q}\,\frac{dq}{d\rho} = \frac{r/a}{q}\cdot a\,\frac{dq}{dr} = \frac{r}{q}\,\frac{dq}{dr}$$

对应代码 (`trcoef.f90:390`)：

```fortran
S(NR) = RHOG(NR) / QL * DQ    ! DQ = dq/d_rho
```

### 8.5 $\alpha$ 中的压力梯度

归一化压力梯度 $\alpha$ 的公式中使用的是**对 $r$ 的导数**：

$$\alpha = -\frac{2\,\mu_0\,q^2\,R}{B^2}\,\frac{dp}{dr}$$

代码中 `DPP` 已经是 $dp/dr$（因为乘了 `DRL = 1/(DR*RA)`），所以：

```fortran
ALPHA(NR) = -2.D0 * RMU0 * QL**2 * RR / BB**2 * (DPP * 1.D20 * RKEV)
```

其中 `DPP * 1.D20 * RKEV` 将 $[10^{20}\,\mathrm{m}^{-3}\cdot\mathrm{keV/m}]$ 转换为 $[\mathrm{Pa/m}]$。

### 8.6 $E\times B$ 速度和剪切率

$$V_{E\times B}(NR) = -\frac{E_r(NR)}{B}$$

$E\times B$ 剪切率定义为：

$$\omega_{E\times B} = (s - 1)\,\frac{V_{E\times B}}{\rho} + \frac{dV_{E\times B}}{d\rho}$$

注意这里 `DVE = DERIV3(NR, RHOG, VEXB, ...)` 返回的是 $dV_{E\times B}/d\rho$（对 $\rho$ 的导数），而第一项中的 $V_{E\times B}/\rho$ 也是 $\rho$ 坐标下的量，所以整个表达式在 $\rho$ 坐标下自洽。

对应代码 (`trcoef.f90:146, 405-406`)：

```fortran
VEXB(NR) = -ER(NR) / BB
DVE = DERIV3(NR, RHOG, VEXB, NRMAX, 1)   ! dV_ExB/d_rho
WEXB(NR) = (S(NR)-1.D0) * VEXB(NR)/RHOG(NR) + DVE
```

### 8.7 输运系数的时间松弛

为了防止数值振荡，输运系数的变化受时间松弛约束：

$$\tau_K(NR) = \frac{q\,R}{v_{Ti}} \cdot \mathrm{MDLTC}$$

其中 `MDLTC` 是用户设定的松弛倍数（默认 0 表示不松弛）。

---

## 9. 输入参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `MDLKAI` | 31 | 输运模型选择器 |
| `CK0` | 12.0 | 电子 $\chi$ 乘数 |
| `CK1` | 12.0 | 离子 $\chi$ 乘数 |
| `CWEB` | 1.0 | $E\times B$ 抑制系数 |
| `CALF` | 1.0 | 形状因子中的 $\alpha$ 缩放因子 |
| `MDLCD05` | 0 | 0：原始 CDBM；非零：CDBM05（含延伸度修正） |
| `model_cdbm_smooth` | 1 | 0：不平滑；1：5 点平滑 + 剪切下限 + $\chi$ 限幅 |
| `MDLTC` | 0 | 输运系数时间松弛倍数（0 = 不松弛） |

---

## 10. 输出变量

### 10.1 输运系数

| 变量 | 含义 |
|------|------|
| `AKDW(NR, 1)` = `AKE` | 电子反常热扩散系数 $\chi_e$ [m$^2$/s] |
| `AKDW(NR, 2)` = `AKD` | 离子 (D) 反常热扩散系数 $\chi_i$ [m$^2$/s] |
| `AKDW(NR, 3)` | 离子 (T) 反常热扩散系数 [m$^2$/s] |
| `AKDW(NR, 4)` | 离子 (He) 反常热扩散系数 [m$^2$/s] |

### 10.2 诊断变量（CASE 30:39）

| 变量 | 内容 |
|------|------|
| `VGR1(NR, 1)` | 形状因子 $F$ |
| `VGR1(NR, 2)` | 磁剪切 $s$ |
| `VGR1(NR, 3)` | 归一化 $\alpha$ |
| `VGR2(NR, 1)` | 径向电场 $E_r$ |
| `VGR2(NR, 2)` | $E\times B$ 漂移速度 $V_{E\times B}$ |
| `VGR3(NR, 1)` | $E\times B$ 抑制因子 $f_{E\times B}$ |
| `VGR3(NR, 2)` | $|\omega_{E\times B}|$ |

### 10.3 诊断变量（CASE 130:139）

| 变量 | 内容 |
|------|------|
| `VGR1(NR, 2)` | 磁剪切 $s$ |
| `VGR1(NR, 3)` | 归一化 $\alpha$ |
| `VGR2(NR, 1)` | 径向电场 $E_r$ |
| `VGR2(NR, 2)` | $E\times B$ 漂移速度 $V_{E\times B}$ |
| `VGR2(NR, 3)` | $E\times B$ 剪切率 $\omega_{E\times B}$ |

### 10.4 CSV 输出

通过 `TRGR1D` $\to$ `TR_WRITE_CSV`，以下剖面输出到 CSV 文件：

- `AKE, AKNCE, AKDWE [m²/s] vs r` — 电子 $\chi$ 各分量
- `AKD, AKNCD, AKDWD [m²/s] vs r` — 离子 $\chi$ 各分量

---

## 11. 磁剪切对输运系数的影响路径

在 `MDLKAI = 134` 中，磁剪切 $s$ 通过**两条路径**影响最终输运系数：

### 路径 1：直接通过形状因子

$$s \;\longrightarrow\; s_\mathrm{eff} \;\longrightarrow\; \mathrm{tr\_cdbm}(\ldots) \;\longrightarrow\; \chi_\mathrm{cdbm} \;\longrightarrow\; \chi_e^{(\mathrm{AKDWE})}$$

形状因子 $F(s,\,\alpha,\,\kappa_q)$ 中的关键依赖：

- $\sigma = s - \alpha$：当 $s \to \alpha$ 时，$\sigma \to 0$，$F_1$ 可能跳变
- $1/s^2$：当 $|s| \to 0$ 时，$F_2$ 发散

### 路径 2：间接通过 $E\times B$ 抑制

$$s \;\longrightarrow\; S_L = s^2 + 0.01 \;\longrightarrow\; \omega_{E1} = -\frac{qR}{S_L \cdot V_A}\,\frac{dV_{E\times B}}{dr}$$

$$\longrightarrow\; F_{E\times B}(|\omega_{E1}|,\;s,\;\alpha) \;\longrightarrow\; C_{E\times B} \;\longrightarrow\; \mathrm{tr\_cdbm}(\ldots,\;C_{E\times B},\;\ldots) \;\longrightarrow\; \chi_\mathrm{cdbm}$$

$F_{E\times B}$ 本身也直接依赖 $s$ 的符号（$\beta$ 和 $\gamma$ 的计算公式在 $s < 0$ 和 $s \geq 0$ 时不同）。

### 边界区域 ($0.8 \leq \rho \leq 1.0$) 的震荡原因

1. $q(\rho)$ 在边界附近变化剧烈，数值导数放大噪声
2. $s$ 中的剪切项与径向因子共同继承并放大这些噪声
3. $\mathrm{TRCOFS}$ 中的 $1/s^2$ 进一步放大
4. $F_{E\times B}$ 中 $s$ 的符号切换可能导致公式分支跳变
5. 以上效应叠加，导致 $\chi_e^{(\mathrm{AKDWE})}$ 在边界区域震荡

---

## 12. 局部 CDBM 闭合与磁通坐标输运方程的关系

这一节专门回答一个容易混淆的问题：

> CDBM 部分很多中间量看起来带有“半径型导数”的结构，那么最终得到的输运系数是不是“相对于几何半径 $r$ 的扩散系数”？

答案是：

- **从输运方程使用上看**，TR 的主径向坐标是磁通归一变量
  $$\rho_{TR}=\sqrt{\psi_{t,n}}$$
  并不是简单的几何 $r/a$；
- **从局部 CDBM 实现上看**，代码内部又引入了一个半径型尺度因子 `RHOG`，在当前这类算例里数值上满足
  $$RHOG \approx A_{\rho,\mathrm{exp}}\,\rho_{TR}$$
  其中 $A_{\rho,\mathrm{exp}}\approx 3.6\sim 3.7$；
- 因而 CDBM 既不是“纯粹按几何半径 $r$ 定义”，也不是“完全按裸的 $\rho_{TR}$ 定义”，而是一个**磁通坐标方程 + 半径型局部闭合因子**的组合实现。

也就是说，TRX 里存在两个层次：

### 12.1 层次一：局部湍流闭合模型

局部 CDBM 闭合回答的问题是：

> 给定某一个磁面上的局部平衡量（压力梯度、磁剪切、曲率、$E\times B$ 剪切），应当得到多大的湍流热扩散率？

在这一层里，核心输入量都带有明显的“半径型局部闭合”意义，但这个“半径”在当前实现中并不总是严格等于几何半径 $r$。

#### 压力梯度

代码中：

```fortran
DRL = 1.D0 / (DR * RA)
DPP = (RPP - RPM) * DRL
```

代码里这里使用的是

$$DRL = \frac{1}{DR\cdot RA}$$

因此该项在实现上等价于一个“按 $RA\,\rho$ 缩放后的半径型梯度”：

$$DPP \approx \frac{\Delta p}{RA\,\Delta\rho}$$

进而

$$\alpha = -\frac{2\mu_0 q^2 R}{B^2}\frac{dp}{dr}$$

因此这里的 $\alpha$ 仍然保留了传统 tokamak 局部模型里“压力梯度驱动”的形式，但要注意其实现依赖于代码采用的等效半径标度，而不应机械理解成严格几何 $r$ 的导数。

#### 磁剪切

代码中先计算：

```fortran
DQ = (QP(NR+1)-QP(NR-1))/(2.D0*DR)
S(NR) = RHOG(NR)/QL * DQ
```

这里 `DQ` 是对 TR 主径向坐标的差分：

$$DQ \approx \frac{dq}{d\rho_{TR}}$$

但真正乘在前面的不是裸的 $\rho_{TR}$，而是 `RHOG(NR)`：

$$s = \frac{RHOG}{q}\frac{dq}{d\rho_{TR}}$$

在当前算例中，我们从输出 CSV 反推出：

$$RHOG \approx A_{\rho,\mathrm{exp}}\,\rho_{TR}, \qquad A_{\rho,\mathrm{exp}}\approx 3.6\sim 3.7$$

因此代码里的 shear 更接近：

$$s \approx \frac{A_{\rho,\mathrm{exp}}\,\rho_{TR}}{q}\frac{dq}{d\rho_{TR}}$$

而不是简单的

$$\frac{\rho_{TR}}{q}\frac{dq}{d\rho_{TR}}$$

更不能直接等同于

$$\frac{r}{q}\frac{dq}{dr}$$

除非进一步证明这个有效半径因子正好对应真实几何半径。

#### ExB 剪切

在 `MDLKAI = 130:139` 中，代码使用：

```fortran
SL  = (S(NR)**2 + 0.1D0**2)
WE1 = -QL*RR/(SL*VA) * DVE
```

这里的 `DVE` 由 `DERIV3(..., RHOG, VEXB, ...)` 计算，因此 ExB 剪切率本身也已经混入了 `RHOG` 这一有效半径标度，而不是单独按几何 $r$ 或裸的 $\rho_{TR}$ 处理。

### 12.2 层次二：磁通坐标下的 1D 输运方程

当局部模型给出热扩散率 $\chi$ 之后，TRX 并不是把它简单塞进

$$-n\chi \frac{dT}{dr}$$

这种几何半径形式，而是写进磁通坐标守恒方程：

$$-\frac{1}{V'}\frac{\partial}{\partial \rho}
\left[
V'\left<|\nabla\rho|^2\right> n\chi \frac{\partial T}{\partial \rho}
\right]$$

这里：

- $V' = dV/d\rho$
- $\langle |\nabla\rho|^2 \rangle$ 负责坐标与磁几何变换
- 自变量已经是磁通面标签 $\rho$

所以最终进入离散矩阵的是：

$$Q_{\mathrm{diff}} = V'\left<|\nabla\rho|^2\right> n\chi \left(-\frac{dT}{d\rho}\right)$$

而不是简单的：

$$Q = -n\chi \frac{dT}{dr}$$

### 12.3 两个层次之间的关系

因此，TRX 中的输运系数应这样理解：

1. **来源上**：$\chi$ 是由局部 CDBM 闭合关系产生，依赖于一组“半径型局部量”，其中最重要的是压力梯度和 shear；但这些量在代码里是通过 $\rho_{TR}$ 与 `RHOG` 的组合实现的；
2. **使用上**：$\chi$ 被放入磁通坐标 $\rho_{TR}$ 的输运方程中，与几何因子 $V'$、$\langle|\nabla\rho|\rangle$、$\langle|\nabla\rho|^2\rangle$ 共同构成真正的热流与散度项。

可以用下面这张“层次图”概括：

$$
\bigl(\text{pressure-gradient-like term},\; \text{shear term with } RHOG,\; \kappa_q,\; E\times B\text{ shear}\bigr)
\;\Longrightarrow\;
\chi_{\mathrm{cdbm}}(\rho)
\;\Longrightarrow\;
V'\langle|\nabla\rho|^2\rangle n\chi \frac{dT}{d\rho}
$$

### 12.4 这对不同磁平衡意味着什么

#### 轴对称、弱 shaping 情况

如果平衡接近圆截面或弱 shaping，并且同时满足：

$$\rho_{TR}=\sqrt{\psi_{t,n}} \approx r/a$$

且

$$RHOG \approx C\,\rho_{TR}$$

其中 $C$ 近似常数，那么把局部 CDBM 闭合近似理解成“半径型闭合”通常问题不大，误差主要体现在一个温和的尺度因子上。

#### 轴对称、强 shaping 情况

如果存在明显的：

- 延伸度
- 三角度
- Shafranov 位移
- 强边界 shaping

那么就不能再把 $\rho_{TR}$ 简单当成几何 $r/a$。此时必须通过：

$$V',\qquad \langle|\nabla\rho|\rangle,\qquad \langle|\nabla\rho|^2\rangle$$

这些几何量，把局部闭合出来的 $\chi$ 正确嵌入输运方程。

#### 真正三维非轴对称平衡

如果平衡是本质 3D 的（例如显著磁波纹、stellarator 几何），那么 TRX 这里这套 1D 磁面平均形式就不能直接认为严格适用。原因不是单个 $\chi$ 的单位问题，而是：

- 局部闭合本身是 tokamak/轴对称风格；
- 输运方程也是轴对称磁面平均形式；
- 三维几何下的输运张量、漂移、轨道平均都可能发生根本变化。

### 12.5 最终结论

最准确的一句话是：

> TRX 中的 CDBM 输运系数，**在使用上明确嵌入磁通坐标 $\rho_{TR}=\sqrt{\psi_{t,n}}$ 的 1D 磁面平均输运方程**；而**在局部闭合实现上，又通过 `RHOG \approx A_{\rho,exp}\rho_{TR}` 这样的有效半径因子引入半径型物理量**。因此它不是纯粹的几何 $r$-坐标闭合，也不是简单的裸 flux 坐标闭合，而是两者的工程化组合。

因此它既不能被简单视为“纯粹对 $r$ 定义的 $\chi$”，也不能脱离几何因子单独解释。

---

## 13. 参考文献

1. M. Yagi et al., "Nonlinear simulation of tearing mode and the resulting anomalous transport", Plasma Phys. Control. Fusion **39** (1997) A421.
2. A. Fukuyama et al., "Transport simulation with current diffusive ballooning mode", Plasma Phys. Control. Fusion **37** (1995) 611.
3. K. Itoh, S.-I. Itoh, A. Fukuyama, "Transport and Structural Formation in Plasmas", IOP Publishing (1999).

---

## 附录 A：物理常数

| 常量 | 符号 | 值 | 单位 |
|------|------|-----|------|
| 基本电荷 | $e$ | $1.602176487 \times 10^{-19}$ | C |
| 电子质量 | $m_e$ | $9.10938215 \times 10^{-31}$ | kg |
| 质子质量 | $m_p$ | $1.672621637 \times 10^{-27}$ | kg |
| 光速 | $c$ | $2.99792458 \times 10^{8}$ | m/s |
| 磁导率 | $\mu_0$ | $4\pi \times 10^{-7}$ | H/m |
| 介电常数 | $\varepsilon_0$ | $1/(c^2 \mu_0)$ | F/m |
| keV 转焦耳 | $e_\mathrm{keV}$ | $1.602 \times 10^{-16}$ | J |

## 附录 B：源文件索引

| 文件 | 内容 |
|------|------|
| `trcoef.f90:TRCFDW` | 主输运系数计算子程序（第 19–1135 行） |
| `trcoef.f90:TRCOFS` | 标准形状因子 $F(s,\,\alpha,\,\kappa_q)$（第 1907–1943 行） |
| `trcoef.f90:TRCOFSS` | 简化形状因子（第 1991–2002 行） |
| `trcoef.f90:TRCOFSX` | 扩展形状因子（含 Shafranov 位移修正）（第 1945–1989 行） |
| `trcoef.f90:TRCOFT` | 备选形状因子（MDLKAI=50 用）（第 2004–2040 行） |
| `trcoef.f90:FEXB` | $E\times B$ 剪切抑制函数（第 2083–2120 行） |
| `trcomm.f90` | 参数声明（CK0, CK1, CWEB, CALF 等） |
| `trinit.f90` | 参数默认值 |
| `trparm.f90` | Namelist 输入定义 |
