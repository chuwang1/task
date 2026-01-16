# TASK/EQU 理论与实现说明文档

> 文档作者：自动生成  
> 更新日期：2026/01/15  
> 程序开发者：Masafumi Azumi (日本原子能研究所 JAERI)

---

## 目录

1. [程序概述](#1-程序概述)
2. [物理背景与数学基础](#2-物理背景与数学基础)
3. [Grad-Shafranov方程](#3-grad-shafranov方程)
4. [数值求解方法](#4-数值求解方法)
5. [自由边界平衡](#5-自由边界平衡)
6. [极向场线圈电流计算](#6-极向场线圈电流计算)
7. [程序结构与流程](#7-程序结构与流程)
8. [主要子程序说明](#8-主要子程序说明)
9. [输入输出说明](#9-输入输出说明)

---

## 1. 程序概述

**TASK/EQU** 是一个用于计算托卡马克等离子体磁流体力学(MHD)平衡的程序。程序求解轴对称环形等离子体的**自由边界平衡**问题，可以：

- 给定等离子体参数（电流、压强剖面等），计算平衡位形
- 自动调整极向场线圈电流以实现目标等离子体形状
- 计算分界面（separatrix）和X点位置
- 输出安全因子、压强等径向剖面

---

## 2. 物理背景与数学基础

### 2.1 托卡马克坐标系

采用柱坐标系 $(R, \phi, Z)$：
- $R$：大半径方向（距对称轴的距离）
- $\phi$：环向角
- $Z$：垂直方向

### 2.2 磁场表示

在轴对称条件下，磁场可以表示为：

$$\mathbf{B} = \nabla\phi \times \nabla\psi + F(\psi)\nabla\phi$$

其中：
- $\psi(R,Z)$：极向磁通函数
- $F(\psi) = RB_\phi$：极向电流函数

磁场分量：
$$B_R = -\frac{1}{R}\frac{\partial\psi}{\partial Z}$$
$$B_Z = \frac{1}{R}\frac{\partial\psi}{\partial R}$$
$$B_\phi = \frac{F(\psi)}{R}$$

极向磁场：
$$B_p = \frac{|\nabla\psi|}{R} = \sqrt{B_R^2 + B_Z^2}$$

### 2.3 MHD平衡条件

等离子体处于平衡态时，满足力平衡方程：

$$\mathbf{J} \times \mathbf{B} = \nabla p$$

其中 $\mathbf{J}$ 是电流密度，$p$ 是压强。

由 Maxwell 方程：
$$\mu_0 \mathbf{J} = \nabla \times \mathbf{B}$$

---

## 3. Grad-Shafranov方程

### 3.1 方程推导

将力平衡方程和 Maxwell 方程结合，可以推导出 **Grad-Shafranov (GS) 方程**：

$$\Delta^* \psi = -\mu_0 R^2 \frac{dp}{d\psi} - F\frac{dF}{d\psi}$$

其中 $\Delta^*$ 是修正的拉普拉斯算符：

$$\Delta^* \psi = R\frac{\partial}{\partial R}\left(\frac{1}{R}\frac{\partial\psi}{\partial R}\right) + \frac{\partial^2\psi}{\partial Z^2}$$

展开形式：
$$\Delta^* \psi = \frac{\partial^2\psi}{\partial R^2} - \frac{1}{R}\frac{\partial\psi}{\partial R} + \frac{\partial^2\psi}{\partial Z^2}$$

### 3.2 右端项的物理意义

GS方程右端项与环向电流密度相关：

$$j_\phi = -\frac{1}{\mu_0 R}\Delta^*\psi = R\frac{dp}{d\psi} + \frac{F}{\mu_0 R}\frac{dF}{d\psi}$$

- 第一项 $R \frac{dp}{d\psi}$：压强梯度驱动的抗磁电流
- 第二项 $\frac{F}{\mu_0 R}\frac{dF}{d\psi}$：极向电流梯度产生的电流

### 3.3 边界条件

- **磁轴**：$\psi$ 在磁轴处取极值，$\nabla\psi = 0$
- **等离子体边界**：$\psi = \psi_{\text{boundary}}$（通常归一化为0或常数）
- **远场边界**：由外部线圈和等离子体电流产生的磁通

---

## 4. 数值求解方法

### 4.1 计算网格

程序在 $(R, Z)$ 平面上建立均匀网格：

$$R_i = R_{\min} + (i-1)\Delta R, \quad i = 1, \ldots, N_R$$
$$Z_j = Z_{\min} + (j-1)\Delta Z, \quad j = 1, \ldots, N_Z$$

### 4.2 有限差分离散化

GS方程的差分格式：

$$\frac{\psi_{i+1,j} - 2\psi_{i,j} + \psi_{i-1,j}}{(\Delta R)^2} - \frac{1}{R_i}\frac{\psi_{i+1,j} - \psi_{i-1,j}}{2\Delta R} + \frac{\psi_{i,j+1} - 2\psi_{i,j} + \psi_{i,j-1}}{(\Delta Z)^2} = S_{i,j}$$

其中源项：
$$S_{i,j} = -\mu_0 R_i^2 \left.\frac{dp}{d\psi}\right|_{\psi_{i,j}} - F(\psi_{i,j})\left.\frac{dF}{d\psi}\right|_{\psi_{i,j}}$$

### 4.3 迭代求解

由于右端项依赖于 $\psi$，需要**迭代求解**：

1. 给定初始猜测 $\psi^{(0)}$
2. 计算源项 $S^{(n)}$ 基于 $\psi^{(n)}$
3. 求解线性方程组得到 $\psi^{(n+1)}$
4. 检查收敛性：$\|\psi^{(n+1)} - \psi^{(n)}\| < \epsilon$
5. 未收敛则返回步骤2

### 4.4 FCT (Flux Coordinate Transform) 方法

程序使用 FCT 方法加速收敛：

- 引入混合参数 $\alpha$（对应参数 `bavmax`, `bavmin`）
- 更新公式：$\psi^{(n+1)} = \alpha \tilde{\psi}^{(n+1)} + (1-\alpha)\psi^{(n)}$
- 根据收敛情况自适应调整 $\alpha$

---

## 5. 自由边界平衡

### 5.1 问题描述

自由边界平衡中，等离子体边界形状由**极向场线圈电流**和**等离子体电流**共同决定。

总磁通：
$$\psi(R,Z) = \psi_{\text{plasma}}(R,Z) + \psi_{\text{coil}}(R,Z)$$

### 5.2 线圈产生的磁通

单匝线圈在位置 $(R_0, Z_0)$ 产生的磁通：

$$\psi_{\text{coil}}(R,Z) = \frac{\mu_0 I}{2\pi} \sqrt{RR_0} \left[\left(\frac{2}{k} - k\right)K(k) - \frac{2}{k}E(k)\right]$$

其中：
- $k^2 = \frac{4RR_0}{(R+R_0)^2 + (Z-Z_0)^2}$
- $K(k)$, $E(k)$ 是第一类和第二类完全椭圆积分

程序中使用查表法加速椭圆积分计算（见 `flux0`, `flux`, `flxfun` 函数）。

### 5.3 等离子体边界条件

等离子体边界上设置**标记点** $(R_j^s, Z_j^s)$，$j = 1, \ldots, J^s$

边界条件：边界上磁通为常数
$$\psi(R_j^s, Z_j^s) = \psi_{\text{boundary}} = \text{const}$$

---

## 6. 极向场线圈电流计算

### 6.1 最小二乘法

线圈电流 $I_k^v$（$k = 1, \ldots, K^v$）通过最小化以下变分确定：

$$\delta W(I_k^v) = \sum_{j=1}^{J^s} w_j^s \left(\psi^p(\mathbf{r}_j^s) + \sum_{k=0}^{K^v} I_k^v \psi^k(\mathbf{r}_j^s) - \psi_j^m\right)^2 + \sum_{k=1}^{K^v} w_k^v (I_k^v - I_k^{v0})^2$$

其中：
- $w_j^s$：标记点权重
- $w_k^v$：线圈电流权重（对应参数 `CVACWG`）
- $\psi^p$：等离子体电流产生的磁通
- $\psi^k$：第 $k$ 个线圈单位电流产生的磁通
- $\psi_j^m$：标记点目标磁通（通常为0）
- $I_k^{v0}$：线圈电流初始值（对应参数 `CVACST`）

### 6.2 法方程

对 $I_k^v$ 求变分得到线性方程组：

$$\sum_{l=0}^{K^v} A_{kl} I_l^v = b_k$$

其中：
$$A_{kl} = \sum_{j=1}^{J^s} w_j^s \psi^k(\mathbf{r}_j^s) \psi^l(\mathbf{r}_j^s) + \delta_{kl} w_k^v$$
$$b_k = \sum_{j=1}^{J^s} w_j^s \left(\psi_j^m - \psi^p(\mathbf{r}_j^s)\right) \psi^k(\mathbf{r}_j^s) + w_k^v I_k^{v0}$$

### 6.3 固定电流处理

当 `IVAC(k) < 0` 时，第 $k$ 个线圈电流固定：
- 将该线圈产生的磁通计入等离子体磁通
- 从变分目标中排除该线圈

---

## 7. 程序结构与流程

### 7.1 主程序流程

```
主程序 (eqmain.f)
├── equ_init     ← 初始化
├── tr_init      ← 输运初始化
├── equ_parm     ← 读取参数文件 (equparm)
└── eqmenu       ← 主菜单/计算循环
    ├── eqsetu   ← 初始设置
    │   ├── eqgrd   ← 设置网格
    │   ├── eqflin  ← 读取线圈数据
    │   ├── eqchek  ← 检查输入
    │   └── eqview  ← 显示参数
    ├── eqequ    ← 初始平衡求解
    │   ├── eqpds0  ← 设置压强/电流剖面
    │   ├── eqrcu   ← 计算电流密度
    │   ├── eqbnd   ← 边界条件
    │   ├── eqpde   ← 求解PDE
    │   ├── eqadj   ← 调整线圈电流
    │   ├── eqaxi   ← 寻找磁轴
    │   ├── eqsep   ← 寻找分界面
    │   └── eqlin   ← 追踪边界
    └── eqfct    ← FCT平衡求解
        ├── eqode   ← 求解ODE
        └── eqchk   ← 收敛检查
```

### 7.2 迭代收敛控制

| 参数 | 说明 | 典型值 |
|------|------|--------|
| `MSETUP` | 初始平衡最大迭代次数 | 20-100 |
| `ESETUP` | 初始平衡收敛判据 | 1.0E-3 ~ 1.0E-7 |
| `IEQMAX` | FCT最大迭代次数 | 20 |
| `EEQMAX` | FCT收敛判据 | 1.0E-3 |
| `IADMAX` | 真空场调整迭代次数 | 1-50 |
| `BAVMAX` | FCT混合参数上限 | 0.8 |
| `BAVMIN` | FCT混合参数下限 | 0.2 |

---

## 8. 主要子程序说明

### 8.1 核心求解子程序

| 子程序 | 功能 | 文件 |
|--------|------|------|
| `eqpde` | 求解GS偏微分方程（FFT+三对角） | eqsub.f |
| `eqrcu` | 计算电流密度分布 $R \cdot j_\phi$ | eqsub.f |
| `eqode` | 求解磁面平均量ODE | eqfct.f |

### 8.2 边界与几何子程序

| 子程序 | 功能 | 文件 |
|--------|------|------|
| `eqbnd` | 计算边界磁通 | eqsub.f |
| `eqaxi` | 搜索磁轴位置 | eqsub.f |
| `eqsep` | 搜索X点和分界面 | eqsub.f |
| `eqlin` | 追踪等离子体边界 | eqsub.f |
| `eqtrc` | 磁力线追踪 | eqsub.f |

### 8.3 线圈电流计算

| 子程序 | 功能 | 文件 |
|--------|------|------|
| `eqadj` | 最小二乘法调整线圈电流 | eqsub.f |
| `flux` | 计算线圈磁通（椭圆积分） | eqsub.f |
| `flxfun` | 批量计算磁通 | eqsub.f |
| `flxtab` | 初始化椭圆积分查表 | eqsub.f |

### 8.4 剖面函数

| 子程序 | 功能 | 文件 |
|--------|------|------|
| `eqpds0` | 设置 $dp/d\psi$ 和 $dF/d\psi$ | eqsub.f |
| `eqpfds` | 自定义剖面（ICP=9时） | eqpfds.f |

---

## 9. 输入输出说明

### 9.1 输入文件

| 文件 | 说明 |
|------|------|
| `equparm` | 主参数文件（Namelist格式） |
| `coildata` | 极向场线圈几何数据 |

详细参数说明请参考 [equparm_coildata_guide.md](equparm_coildata_guide.md)

### 9.2 主要输出量

| 变量 | 说明 | 单位 |
|------|------|------|
| `psi(R,Z)` | 极向磁通分布 | Wb |
| `raxis, zaxis` | 磁轴位置 | m |
| `saxis` | 磁轴处磁通值 | Wb |
| `qaxis` | 轴心安全因子 | - |
| `qsurf` | 边界安全因子 | - |
| `qqv(n)` | 安全因子径向剖面 | - |
| `prv(n)` | 压强径向剖面 | Pa |
| `rsu, zsu` | 等离子体边界坐标 | m |
| `cvac(k)` | 线圈电流 | AT |

### 9.3 等离子体参数输出

| 变量 | 说明 |
|------|------|
| `beta` | 环向比压 $\beta_t = 2\mu_0 \langle p \rangle / B_0^2$ |
| `bets` | 极向比压 $\beta_p$ |
| `ell` | 椭圆度（elongation） |
| `trg` | 三角形变形度（triangularity） |
| `zzli` | 内感（internal inductance） |

---

## 附录 A: 椭圆积分近似公式

程序使用以下多项式近似计算椭圆积分：

**第一类完全椭圆积分** $K(k)$：
$$K(k) \approx (a_0 + a_1 m + a_2 m^2) + (b_0 + b_1 m + b_2 m^2) \ln(1/m)$$

**第二类完全椭圆积分** $E(k)$：
$$E(k) \approx (c_0 + c_1 m + c_2 m^2) + m(d_1 + d_2 m) \ln(1/m)$$

其中 $m = 1 - k^2$，系数：
- $a_0 = 1.3862944$, $a_1 = 0.1119723$, $a_2 = 0.0725296$
- $b_0 = 0.5$, $b_1 = 0.1213478$, $b_2 = 0.0288729$
- $c_0 = 1.0$, $c_1 = 0.4630151$, $c_2 = 0.1077812$
- $d_1 = 0.2452727$, $d_2 = 0.0412496$

---

## 附录 B: 剖面函数选项

### ICP(1) = 1（默认）
$$\frac{dp}{d\psi} \propto (1 - c_4)(1 - v^{c_2})^{c_3} + c_4$$
$$\frac{dF}{d\psi} \propto (1 - c_8)(1 - v^{c_6})^{c_7} + c_8$$

其中 $v = V/V_{\text{plasma}}$ 是归一化体积坐标。

### ICP(1) = 2
$$\frac{dp}{d\psi} \propto 1 - c_2 v^{c_3} - (1 - c_2 + c_4) v^{c_5}$$
$$\frac{dF}{d\psi} \propto 1 - c_6 v^{c_7} - (1 - c_6 + c_8) v^{c_9}$$

### ICP(1) = 11（中空电流剖面）
$$\frac{dp}{d\psi} \propto (1 - c_4)(1 - v^{c_2})^{c_3} + c_4$$
$$\frac{dF}{d\psi} \propto (1 - v^{c_6})^{c_7} \left(1 + c_8 \exp\left(-\left(\frac{v - c_9}{c_{10}}\right)^2\right)\right)$$

---

## 参考文献

1. V.D. Shafranov, "Plasma Equilibrium in a Magnetic Field", Reviews of Plasma Physics, Vol. 2 (1966)
2. H. Grad and H. Rubin, "Hydromagnetic Equilibria and Force-Free Fields", Proceedings of the Second United Nations International Conference on the Peaceful Uses of Atomic Energy, Vol. 31 (1958)
3. L.L. Lao et al., "Reconstruction of current profile parameters and plasma shapes in tokamaks", Nuclear Fusion 25 (1985) 1611
4. M. Azumi et al., 参考文档目录下的 `GradPaper3789.pdf`, `Azumi_equ.pdf`
