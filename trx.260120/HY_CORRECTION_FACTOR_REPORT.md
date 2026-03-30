# HY_CORRECTION_FACTOR 分析报告

## 1. 背景

在聚变等离子体中，DT（或 DD、D-He3）聚变反应产生的高能 alpha 粒子在热化过程中会将能量分配给背景电子和离子。这个分配比例由经典的快粒子慢化理论决定，核心是 **Stix 慢化函数** $H(y)$。

TR 代码中 `HY_CORRECTION_FACTOR` 是一个人为修正因子，用于调整 $H(y)$ 的输出值，从而改变离子加热占比。本报告首先阐述原始（未修正）的能量分配原理，然后说明修正的实现方式。

---

## 2. 原始版本的能量分配原理

### 2.1 Alpha 粒子的产生

DT 聚变反应：

$$D + T \rightarrow \alpha\,(3.5\,\mathrm{MeV}) + n\,(14.1\,\mathrm{MeV})$$

alpha 粒子以固定能量 $E_\alpha = 3.5\,\mathrm{MeV}$ 出生（D-He3 反应为 $3.6\,\mathrm{MeV}$），对应的出生速度为：

$$v_f = \sqrt{\frac{2E_\alpha}{m_\alpha}}$$

对应代码 (`trpnf.f90:162`)：

```fortran
VF = SQRT(2.D0*3.5D3*RKEV/AMA)
```

### 2.2 临界速度

快粒子在背景等离子体中减速时，受两种阻力：
- **电子阻力**：与 $v^3$ 成正比（高速时主导）
- **离子阻力**：与 $v^0$ 近似（低速时主导）

两者相等时的速度称为**临界速度** $v_c$，由下式定义：

$$v_c^3 = v_{c,D}^3 + v_{c,T}^3 + v_{c,He}^3$$

其中每个离子物种的贡献为：

$$v_{c,j}^3 = \frac{3\sqrt{\pi/2}\,m_e}{n_e}\left(\frac{T_e}{m_e}\right)^{3/2} \frac{n_j Z_j^2}{m_j}$$

对应代码 (`trpnf.f90:265-270`)：

```fortran
P1   = 3.D0*SQRT(0.5D0*PI)*AME/ANE*(ABS(TE)*RKEV/AME)**1.5D0
VCD3 = P1*RN(NR,NS_D)*PZ(NS_D)**2/AMD
VCT3 = P1*RN(NR,NS_T)*PZ(NS_T)**2/AMT
VCA3 = P1*RN(NR,NS_He4)*PZ(NS_He4)**2/AMA
VC3  = VCD3+VCT3+VCA3
VCR  = VC3**(1.D0/3.D0)
```

### 2.3 Stix H(y) 函数

**物理含义：** $H(y)$ 给出的是快粒子在从出生速度 $v_f$ 减速到热速度过程中，**传递给背景离子的能量占总传递能量的比例**。

$$y = \frac{v_f}{v_c}$$

$$H(y) = \frac{2}{y^2}\left[\frac{1}{6}\ln\frac{y^3+1}{(y+1)^3} + \frac{1}{\sqrt{3}}\left(\arctan\frac{2y-1}{\sqrt{3}}+\frac{\pi}{6}\right)\right]$$

对应代码 (`trpnf.f90:651-661`)：

```fortran
FUNCTION HY(V)
  USE TRCOMM, ONLY : PI,rkind
  IMPLICIT NONE
  REAL(rkind), INTENT(IN) :: V
  REAL(rkind) :: HY
  HY = 2.D0*(LOG((V**3+1.D0)/(V+1.D0)**3)/6.D0 &
       +(ATAN((2.D0*V-1.D0)/SQRT(3.D0))+PI/6.D0)/SQRT(3.D0))/V**2
  RETURN
END FUNCTION HY
```

**极限行为：**

| 条件 | $y = v_f/v_c$ | $H(y)$ | 物理含义 |
|------|--------------|--------|---------|
| $v_f \gg v_c$ | $y \gg 1$ | $H \to 0$ | alpha 主要被电子减速，能量几乎全给电子 |
| $v_f \ll v_c$ | $y \ll 1$ | $H \to 1$ | alpha 主要被离子减速，能量几乎全给离子 |
| $v_f = v_c$ | $y = 1$ | $H \approx 0.28$ | 约 28% 给离子 |

对于典型 CFEDR 参数（$T_e \sim 20\,\mathrm{keV}$，$n_e \sim 1.4\times 10^{20}\,\mathrm{m}^{-3}$），$v_f/v_c$ 通常在 1.5–3 范围，$H$ 大约 0.1–0.25。

### 2.4 慢化时间

Spitzer 慢化时间：

$$\tau_s = \frac{0.2\,A_\alpha\,T_e^{3/2}}{Z_\alpha^2\,n_e\,\ln\Lambda}$$

对应代码 (`trpnf.f90:274`)：

```fortran
TAUS = 0.2D0*PA(NS_He4)*ABS(TE)**1.5D0 &
     /(PZ(NS_He4)**2*ANE*COULOG(1,2,ANE,TE))
```

有效热化时间（用于后续功率计算）：

$$\tau_f = \frac{1}{2}\tau_s(1-H)$$

### 2.5 能量分配

设 alpha 粒子的聚变功率密度（待分配的总功率）为 $P_{\mathrm{fus,in}}$，则：

$$P_{\mathrm{fus,in}} = \frac{W_f \cdot e_\mathrm{keV} \cdot 10^{20}}{\tau_f}$$

其中 $W_f$ 是 alpha 粒子压力。

**分给电子的功率：**

$$P_e = (1 - H)\,P_{\mathrm{fus,in}}$$

**分给各离子物种的功率：**

$$P_j = \frac{v_{c,j}^3}{v_c^3}\,H\,P_{\mathrm{fus,in}}$$

其中 $j = D, T, He$。

也就是说，离子加热按各物种对临界速度的贡献来分配。

对应代码 (`trpnf.f90:285-288`)：

```fortran
PNFCL_NSNNFNR(NS_e,  NNF,NR) =     (1.D0-HYF)*PNFIN_NNFNR(NNF,NR)
PNFCL_NSNNFNR(NS_D,  NNF,NR) = (VCD3/VC3)*HYF*PNFIN_NNFNR(NNF,NR)
PNFCL_NSNNFNR(NS_T,  NNF,NR) = (VCT3/VC3)*HYF*PNFIN_NNFNR(NNF,NR)
PNFCL_NSNNFNR(NS_He4,NNF,NR) = (VCA3/VC3)*HYF*PNFIN_NNFNR(NNF,NR)
```

### 2.6 能量守恒验证

$$P_e + P_D + P_T + P_{He} = (1-H)\,P + \frac{v_{c,D}^3+v_{c,T}^3+v_{c,He}^3}{v_c^3}\,H\,P = (1-H+H)\,P = P$$

因此原始模型严格满足能量守恒。

---

## 3. 原始模型应用于三种聚变反应

代码中有三个子程序使用完全相同的 $H(y)$ 逻辑：

| 子程序 | 反应 | Alpha 出生能 | 代码位置 |
|--------|------|-------------|---------|
| `TRNFDT` | D-T | 3.5 MeV | `trpnf.f90:150-292` |
| `TRNFDD` | D-D | 3.5 MeV | `trpnf.f90:321-423` |
| `TRNFDHe3` | D-He3 | 3.6 MeV | `trpnf.f90:543-603` |

---

## 4. HY_CORRECTION_FACTOR 的人为修正

### 4.1 修正动机

在与 TRANSP 代码的对比中发现，TR 计算的离子加热功率比 TRANSP 偏低。

可能的原因包括：
- 经典 Stix 理论假设完全各向同性、Maxwellian 背景，不完全适用于真实等离子体
- 快粒子轨道效应（有限轨道宽度、被俘粒子比例）会修正能量分配
- TRANSP 使用 Monte Carlo 方法追踪快粒子，考虑了更多物理效应
- 两个代码在数值实现细节上的差异也可能贡献一部分

### 4.2 修正实现

**声明** (`trcomm.f90:82`)：

```fortran
REAL(rkind):: HY_CORRECTION_FACTOR
```

**默认值** (`trinit.f90:213`)：

```fortran
HY_CORRECTION_FACTOR = 1.0D0
```

**应用方式** — 在三个子程序中完全相同（例如 `trpnf.f90:273`）：

```fortran
HYF = HY(VF/VCR)
HYF = MIN(HYF * HY_CORRECTION_FACTOR, 1.D0-1.D-10)
```

### 4.3 修正的数学效果

设原始 $H(y)$ 的值为 $H_0$，修正后为：

$$H' = \min\!\left(H_0 \cdot f_\mathrm{corr},\; 1 - 10^{-10}\right)$$

其中 $f_\mathrm{corr}$ = `HY_CORRECTION_FACTOR`。

能量分配变为：

$$P_e' = (1 - H')\,P_{\mathrm{fus,in}}'$$

$$P_{\mathrm{ions}}' = H'\,P_{\mathrm{fus,in}}'$$

注意 $P_{\mathrm{fus,in}}'$ 本身也会变，因为 $\tau_f = \frac{1}{2}\tau_s(1-H')$。

### 4.4 修正对各物理量的影响

设 $f > 1$（如 $f = 1.16$）：

| 量 | 变化方向 | 原因 |
|----|---------|------|
| $H'$ | 增大 | 直接乘以 $f$ |
| 离子加热占比 | 增大 | $H'$ 更大 |
| 电子加热占比 | 减小 | $1 - H'$ 更小 |
| $\tau_f = \frac{1}{2}\tau_s(1-H')$ | 减小 | $(1-H')$ 变小 |
| $P_{\mathrm{fus,in}} = W_f/\tau_f$ | 增大 | $\tau_f$ 变小 |
| 离子温度 | 趋于升高 | 离子加热增强 |
| 电子温度 | 趋于降低 | 电子加热减弱 |

### 4.5 能量守恒

修正后仍然严格满足能量守恒：

$$P_e' + \sum_j P_j' = (1-H') + H' = 1$$

因为修正只改变了分配比例，没有创造或消灭能量。

### 4.6 上限保护

```fortran
HYF = MIN(HYF * HY_CORRECTION_FACTOR, 1.D0-1.D-10)
```

保证 $H' < 1$，避免：
- $1 - H' = 0$ 导致 $\tau_f = 0$，进而 $P_{\mathrm{fus,in}} \to \infty$
- 电子加热变为负值

### 4.7 推荐值

根据 TRANSP 对比：

$$f_\mathrm{corr} \approx 1.16$$

即离子加热占比需要增大约 16% 才能匹配 TRANSP 结果。

---

## 5. 使用方法

在输入文件 `&TR` namelist 中设置：

```fortran
HY_CORRECTION_FACTOR = 1.16   ! 匹配 TRANSP 的离子加热分配
```

默认值为 `1.0`，表示不做修正。

---

## 6. 总结

| 项目 | 原始版本 | 修正后 |
|------|---------|--------|
| 离子加热比例 | $H(v_f/v_c)$ | $\min(H \cdot f_\mathrm{corr},\;1-\epsilon)$ |
| 电子加热比例 | $1 - H$ | $1 - H'$ |
| 能量守恒 | 严格满足 | 严格满足 |
| 物理基础 | 经典 Stix 慢化理论 | 经典 + 经验修正因子 |
| 适用范围 | 所有 DT、DD、D-He3 反应 | 同上 |
| 默认行为 | $f = 1.0$（无修正） | 需用户显式设置 |

### 代码改动总结

| 文件 | 改动 |
|------|------|
| `trcomm.f90:82` | 声明 `REAL(rkind):: HY_CORRECTION_FACTOR` |
| `trinit.f90:213` | 默认值 `HY_CORRECTION_FACTOR = 1.0D0` |
| `trparm.f90:43` | 加入 `NAMELIST /TR/` |
| `trpnf.f90:273` | DT 反应：`HYF = MIN(HYF * HY_CORRECTION_FACTOR, 1.D0-1.D-10)` |
| `trpnf.f90:399` | DD 反应：同上 |
| `trpnf.f90:585` | D-He3 反应：同上 |
