# EQIPQP 算法数学说明与变量总表

## 1. 物理背景：Grad-Shafranov (GS) 方程

托卡马克等离子体在轴对称（即圆环方向上物理量不变化）情况下的理想磁流体力学（MHD）平衡，由著名的 **Grad-Shafranov (GS) 方程**描述。它的完整形式是：

$$
\Delta^* \psi = - \mu_0 R^2 \frac{dp(\psi)}{d\psi} - F(\psi) \frac{dF(\psi)}{d\psi}
$$

其中：
- $\psi$ 代表**极向磁通**（为了推导方便，经常省去下角标 $p$，有时除以 $2\pi$）。
- $R$ 是距离圆柱对称中心轴的半径（大半径方向）。
- $p(\psi)$ 是等离子体的压强剖面（仅是磁通 $\psi$ 的函数）。
- $F(\psi) = R B_\phi$ 是极向电流函数（在此文档的 $F$ 定义中可能带 $2\pi$ 系数，本质是环向磁场 $B_\phi$ 的量度）。
- $\mu_0$ 是真空磁导率。
- $\Delta^*$ 是一个特殊的二阶偏微分算子（称为 Grad-Shafranov 算子），在圆柱坐标系 $(R, Z, \phi)$ 下，由于轴对称性（对 $\phi$ 导数为0），它的完整展开形式为：
  $$ \Delta^* \psi = R \frac{\partial}{\partial R} \left( \frac{1}{R} \frac{\partial \psi}{\partial R} \right) + \frac{\partial^2 \psi}{\partial Z^2} $$

GS 方程告诉我们：空间的磁通分布（左边 $\Delta^*$ 项代表的“扭曲度”）是由**压强梯度**和**极向电流函数梯度**共同支撑的。

## 2. 目标

`EQIPQP` 的目标是在已知压强剖面 $p(\psi_p)$ 和安全因子剖面 $q(\psi_p)$ 的条件下，重建极向电流函数
$F(\psi_p)=2\pi R B_\phi$（代码中常记为 `FIPV/FPSI/TTV`）。

它是一个“从边界向磁轴反推”的一维离散积分问题。它本质上是把上面那个三维（二维轴对称）的偏微分方程，通过“磁面平均”降维成一维的常微分方程，然后进行积分求解。

---

## 3. 通量与归一化定义

- 极向磁通：$\psi_p$
- 环向磁通：$\psi_t$
- 归一化极向磁通：$\psi_{pn}=\psi_p/\psi_{p,a}$
- 归一化环向磁通：$\psi_{tn}=\psi_t/\psi_{t,a}$

代码中对应：

- `PSIPV(NRV)`：$\psi_p$ 网格（`NRV` 网格点）
- `PSIPNV(NRV)`：$\psi_{pn}$
- `PSITV(NRV)`：$\psi_t$

---

## 3. 基本重建变量与物理含义

在开始推导之前，我们需要先定义几个核心变量。我们不仅列出公式，还解释它们背后的物理直觉。

定义中间变量（代表与电流、磁场能量相关的量）：

$$
Y(\psi_p)=\frac{1}{2}F(\psi_p)^2
$$

定义几何-物理耦合量：

$$
X=\frac{RSV}{q}
$$

$$
A=4\pi^2 B_0^2 R_0^2 \frac{X}{AVIR2\cdot VPV}
$$

$$
B=4\pi^2\mu_0 R_0^2 \frac{dp/d\psi_p}{AVIR2}
$$

$$
G=VPV\cdot AVRR2\cdot X
$$

其中 $R_0=$ `RR`，$B_0=$ `BB`，$\mu_0=$ `RMU0`。

### 💡 物理小剧场：为什么电流会和压强、几何形状扯上关系？

在托卡马克这个“甜甜圈”里，等离子体能乖乖待在里面不跑出来，是因为达到了**力的平衡**。这背后有两个核心机制：

1. **“气球抗压”机制（对应 $B$ 项）：**
   等离子体非常热，内部的压强很大，就像一个拼命想往外胀破的气球。为了把这个“气球”捏住，磁场必须产生一个向内的力（洛伦兹力）来抗衡它。
   大家都知道 **磁场力 = 电流 × 磁场**。因此，如果你知道等离子体内部压强变化有多快（这就是压强梯度 $dp/d\psi_p$），你就能算出需要多大的电流去配合磁场产生足够的力。
   👉 **$B$ 就是“压强驱动项”**。它直接告诉你：为了顶住等离子体的压强差，需要多少电流（磁场）的变化。

2. **“麻花变形”机制（对应 $A$ 和 $G$ 项）：**
   托卡马克里的磁感线不是直的，而是像拧麻花一样一圈圈绕着的。这根麻花拧了多少圈，由**安全因子 $q$** 决定。
   同时，这个“甜甜圈”并不是完美的圆柱体，它的横截面会被拉长、挤扁（也就是几何形状在变化）。当你把甜甜圈切成一层层的洋葱圈时，每一层的表面积、体积都在发生复杂的形变。
   磁感线在这些弯曲、变形的表面上穿行，为了维持它原有的“扭曲程度（$q$）”，磁场分布必须随着几何形状的改变而调整。
   👉 **$X = RSV/q$ 是什么？** 
   $RSV$ 是基于环向磁通定义的“等效小半径”。你可以把极向磁通 $\psi_p$ 想象成“竖着切一刀”看到的磁通，把环向磁通 $\psi_t$ 想象成“横着切一刀”看到的磁通。安全因子 $q$ 的定义本质上就是 $\Delta \psi_t / \Delta \psi_p$（绕一圈大圈对应绕几圈小圈）。
   在数学换算中，为了把关于极向磁通的推导，安全且不破坏单位地关联到真实的空间半径上，我们引入了 $X = RSV/q$。它就像是一个“带有扭曲度烙印的空间标尺”，把纯抽象的磁通变化转换成了带有物理距离（米）的实际变化。
   👉 **$G$ 是“几何特征容量”**，它把体积($VPV$)、截面形状($AVRR2$)和扭曲度($q$，藏在$X$里)打包成了一个综合指标。
   👉 **$A$ 是“几何变形的转换系数”**。当我们从洋葱的一层走到下一层时，几何形状 $G$ 发生了变化，这个变化必须由磁场（也就是电流 $F$）来买单。$A$ 就像一个汇率，告诉你“几何形状变了1块钱，电流需要变多少”。

### 3.1 `B` 项的数学推导（为什么是这个形状？）

为了让你不仅知其然，还知其所以然，我们来稍微推导一下 $B$ 项是怎么来的。这涉及到托卡马克中最核心的力学平衡方程——**Grad-Shafranov (GS) 方程**。

**第一步：最基本的力平衡**
在等离子体中，压强向外推的力（梯度）必须被磁场所产生的安培力（洛伦兹力）拉住：
$$ \nabla p = \mathbf{J} \times \mathbf{B} $$
其中 $p$ 是压强，$\mathbf{J}$ 是电流密度，$\mathbf{B}$ 是磁场。

**第二步：拆解磁场和电流**
在圆环形（托卡马克）几何中，磁场和电流可以用极向（绕截面转）和环向（绕大圈转）来分解。经过一系列严谨的矢量代数运算，物理学家们将三维的力平衡方程化简成了一个关于极向磁通 $\psi$ 的一维方程（即 GS 方程的变形）。
在这个方程中，极向电流函数 $F(\psi) = 2\pi R B_\phi$ （也常记为 $f$）满足以下关系：
$$ \frac{1}{\mu_0} F \frac{dF}{d\psi} = - \mu_0 R^2 \frac{dp}{d\psi} - \Delta^* \psi $$
（注：这里的 $\psi$ 有时指代极向磁通，不同文献差一个 $2\pi$，我们这里注重物理结构的推导，不纠结常数细节）。

**第三步：磁面平均（把 2D 变成 1D）**
上面的方程里有 $R^2$，但是在同一个磁通面（同一个洋葱圈层）上，$R$ 是在不断变化的（内侧 $R$ 小，外侧 $R$ 大）。为了变成一个只依赖于 $\psi$ 的一维方程，我们需要对整个磁面求平均，记作 $\langle \dots \rangle$。
我们对等式两边除以 $R^2$ 然后求平均：
$$ \frac{1}{\mu_0} F \frac{dF}{d\psi} \left\langle \frac{1}{R^2} \right\rangle = - \mu_0 \frac{dp}{d\psi} \langle 1 \rangle - \left\langle \frac{\Delta^* \psi}{R^2} \right\rangle $$
由于 $\langle 1 \rangle = 1$，等式左边刚好凑出了 $F \frac{dF}{d\psi}$。
根据微积分求导法则：$\frac{d}{d\psi}\left(\frac{1}{2}F^2\right) = F \frac{dF}{d\psi}$。
我们把 $\frac{1}{2}F^2$ 记作 $Y$，那么：
$$ \frac{1}{\mu_0} \frac{dY}{d\psi} \left\langle \frac{1}{R^2} \right\rangle = - \mu_0 \frac{dp}{d\psi} + \text{几何形变项(对应 A 和 G)} $$

**第四步：整理出压强项 $B$**
我们把等式左边除了 $dY/d\psi$ 之外的东西都除到右边去，专门看看压强这一项变成了什么样：
$$ \frac{dY}{d\psi} = - \mu_0^2 \frac{dp/d\psi}{\langle 1/R^2 \rangle} + \dots $$
为了适应代码里的单位制和归一化（比如大半径 $R_0$，引入 $4\pi^2$ 因子），代码中定义了 $AVIR2 \propto \langle 1/R^2 \rangle$。
稍作单位配平后，压强对 $dY/d\psi$ 的贡献系数就变成了代码里的 $B$ 项：
$$ B \propto \frac{dp/d\psi}{\langle 1/R^2 \rangle} \longrightarrow B = 4\pi^2\mu_0 R_0^2 \frac{dp/d\psi_p}{AVIR2} $$

**结论：**
你看，积分公式里的 $B_i \times (\psi_{i} - \psi_{i-1})$ 其实就是 $\frac{dY}{d\psi} \times d\psi$ 中的压强部分。$B$ 项分母上之所以会出现 $AVIR2$（即 $\langle 1/R^2 \rangle$），是因为我们在对 Grad-Shafranov 方程做磁面平均时，为了把 $F dF/d\psi$ 提取出来，强行在等式两边除以了 $R^2$！

### 3.2 `A` 和 `G` 项的数学推导（解密几何变形）

刚才我们在第三步留了一个“尾巴”：
$$ \frac{1}{\mu_0} \frac{dY}{d\psi} \left\langle \frac{1}{R^2} \right\rangle = - \mu_0 \frac{dp}{d\psi} - \left\langle \frac{\Delta^* \psi}{R^2} \right\rangle $$
等式最右边那个长得很可怕的项 $-\langle \Delta^* \psi / R^2 \rangle$，就是代码中 $A$ 和 $G$ 的真正来源。

**第五步：磁面几何的“降维打击”**
$\Delta^* \psi$ 是一个非常复杂的空间二阶偏导数算符，它描述了甜甜圈内部等磁通面在三维空间中的“弯曲、挤压和变形程度”。
幸运的是，在磁面坐标系的微积分（结合高斯散度定理）中，物理学家发现：**对这个复杂的二阶算子求整个磁面的平均，可以被降维变成一个“纯几何物理量”的一阶导数！**

具体的数学等式是这样的：
$$ -\left\langle \frac{\Delta^* \psi}{R^2} \right\rangle \propto \frac{1}{V'} \frac{d}{d\psi_p} \left( V' \left\langle \frac{|\nabla \psi_p|^2}{R^2} \right\rangle \right) $$
这里的 $V'$ 代表了磁面体积随磁通的变化率（就是代码里的 $VPV$）。
而 $\langle |\nabla \psi_p|^2 / R^2 \rangle$ 描述了截面轮廓的形变和磁场的梯度关系。

**第六步：打包成“几何综合容量” $G$**
为了把这个一阶导数算得准，代码把括号里的那一大坨东西打包起来，并且结合了描述磁场扭曲度的 $X = RSV/q$（用来转化单位和归一化）。

为什么非要引入 $X = RSV/q$？
在推导中，很多项是基于极向磁通 $\psi_p$ 求导的，而实际的物理形状最好用空间半径来衡量。$RSV$ 是一个等效的小半径，定义为 $RSV \propto \sqrt{\psi_t}$（与环向磁通相关）。而 $q$ 的定义本质上是 $d\psi_t / d\psi_p$。所以，把极向的变量和环向的变量转换来转换去时，$X = RSV/q$ 就成了一个绝佳的“桥梁”，把磁通的微小变化精准映射到了真实空间的几何距离上。

最终，这个大包裹被定义为 $G$：
$$ G = VPV \cdot AVRR2 \cdot X $$
这里的 $AVRR2$ 在代码中正比于截面的几何形变项 $\langle |\nabla \psi_p|^2 / R^2 \rangle$。
所以，$G$ 就是某一层磁通面的**“综合几何变形容量”**。那复杂的一阶导数部分，就变成了简简单单的 $dG/d\psi_p$。

**第七步：拼凑系数 $A$**
既然复杂的算子变成了 $dG/d\psi_p$，我们把公式写完整：
$$ \frac{dY}{d\psi} = \text{压强项} + \frac{\text{常数}}{\langle 1/R^2 \rangle V'} \frac{dG}{d\psi_p} $$
你看，为了平衡等式两边，在这个导数 $dG/d\psi_p$ 的前面，自然而然地多出了一堆系数。
分母上出现了 $\langle 1/R^2 \rangle$（也就是 $AVIR2$）和 $V'$（也就是 $VPV$）。为了把单位配平，并且引入大半径 $R_0$、参考磁场 $B_0$ 等常数，代码将前面这坨系数定义为了 $A$：
$$ A = 4\pi^2 B_0^2 R_0^2 \frac{X}{AVIR2 \cdot VPV} $$

**最终结论：**
原本三维空间里极其可怕的二阶几何变形项 $\langle \Delta^* \psi / R^2 \rangle d\psi_p$，经过“降维打击”后，完美地化简成了：
$$ A \times dG $$
这也就是为什么在离散代码里，我们只需要算一个乘法和减法：$A_i \times (G_i - G_{i-1})$，就能精准地等效三维空间里复杂的甜甜圈几何形变对电流的约束作用！

## 4. 离散反推公式（Fortran `EQIPQP`）

为了在计算机上计算，我们需要把连续的空间切成一圈一圈的“网格”。假设我们将等离子体从磁轴（中心）到边界分为 $N$ 层。下标 $N$ 代表最外层的边界，下标 $1$ 代表最中心的磁轴。我们已知边界的条件，要一层层往里推算。

### 4.1 边界条件（起点与安培环路定理）

最外层（等离子体边界）的极向电流函数 $F$ 主要由外部线圈产生的真空环向磁场决定，因此：

$$
F_N = 2\pi R_0 B_0
$$

**为什么是这个公式？这里的 $\mu_0$ 去哪了？**

根据高中物理的**安培环路定理**，如果我们绕着托卡马克的对称中心画一个半径为 $R$ 的圆，那么环向磁场 $B_\phi$ 和穿过这个圆的极向总电流 $I$ 满足：
$$ 2\pi R B_\phi = \mu_0 I $$

在等离子体物理中，为了计算磁场方便，科学家们并没有把 $F$ 定义为真实的电流安培数，而是**把等式左边直接定义为“极向电流函数” $F$**，即：
$$ F \equiv 2\pi R B_\phi $$
所以，$F$ 的物理本质其实是 $\mu_0 I$，它的单位是 **特斯拉·米 (T·m)** 而不是安培 (A)。$\mu_0$ 被“吸收”进了 $F$ 的定义中。

在等离子体的最外层（真空区域），没有等离子体电流的干扰，穿过中心孔的只有外部装置线圈提供的恒定电流。因此，在真空中 $2\pi R B_\phi$ 是一个处处相等的常数。我们取参考点大半径 $R_0$ 和该处的真空磁场 $B_0$，就得到了边界起跑线：$F_N = 2\pi R_0 B_0$。

为了方便积分计算，我们定义一个中间变量 $Y = \frac{1}{2}F^2$。那么在边界处：

$$
Y_N = \frac{1}{2}F_N^2
$$

### 4.2 从第 $i$ 层推导第 $i-1$ 层（向内走一步）

假设我们已经知道了第 $i$ 层的 $Y_i$，想要算出往里走一步的第 $i-1$ 层的 $Y_{i-1}$。
在微积分中，如果知道导数和其中一点的值，可以通过梯形面积公式估算另一点的值（也就是所谓的**梯形积分法**）。在物理上，这种递推关系可以写成：

$$
Y_{i-1} = Y_i + \text{压强变化的贡献} + \text{几何和磁场共同的贡献}
$$

具体写成公式就是：

$$
Y_{i-1}
=Y_i
+\underbrace{ \frac{B_i+B_{i-1}}{2}\,(\psi_{p,i}-\psi_{p,i-1}) }_{\text{这部分来自等离子体的压强差}}
+\underbrace{ \frac{A_i+A_{i-1}}{2}\,(G_i-G_{i-1}) }_{\text{这部分来自几何形状和磁场的耦合}}
$$

只要算出了 $Y_{i-1}$，就可以通过开根号还原出我们真正想要的 $F_{i-1}$：

$$
F_{i-1}=\sqrt{2Y_{i-1}}
$$

**这就是一段一段的“接力赛”：**
算出 $F_{N-1}$，然后再用 $F_{N-1}$ 算出 $F_{N-2}$，一路算回磁轴 $F_1$。

这对应代码中的变量为：

- `YP = 0.5D0*FIPV(NRV)**2`  （相当于上面公式里的 $Y_i$）
- `YM = YP + 0.5*(BLP+BLM)*(PSIPLP-PSIPLM) + 0.5*(ALP+ALM)*(FLP-FLM)`  （相当于算出了 $Y_{i-1}$）
- `FIPV(NRV-1)=SQRT(2*YM)`  （开根号得到 $F_{i-1}$）

### 4.3 轴端特殊处理（为什么最中心要不一样？）

当积分进行到最内层（$i=2$），即要算出最中心的 $F_1$ 时，有一个数学上的麻烦。
公式中有一个变量叫 $AVIR2$，它代表的是几何上 $\langle 1/R^2 \rangle$ 的平均值，而 $VPV$ 代表体积随磁通的变化率。在几何中心（磁轴）处，极向磁通 $\psi_p$ 的变化非常平缓，有些项（例如中心极小体积的计算）会出现“除以极小值”或者“ $0 \times \infty$ ”之类的奇异问题（Singularity）。

为了避免在中心点（$i=1$）直接除以 $VPV_1$ 导致程序崩溃或误差爆炸（Fortran 特例），程序做了一个小修正。

当 `NRV=2` 时，求平均贡献的 $ALM$ 项，巧妙地借用了外一点（$i=2$ 即 `NRV`）的参数 $VPV_2$，和 $X_2$（代码中的 `XP`）进行组合：

$$
ALM=4\pi^2 B_0^2R_0^2\frac{X_2}{AVIR2_{1}\,VPV_2}
$$

而对于一般不是中心的点，则是正常使用它自己（$i-1$ 也就是 `NRV-1`）这一层的参数：

$$
ALM=4\pi^2 B_0^2R_0^2\frac{X_{i-1}}{AVIR2_{i-1}\,VPV_{i-1}}
$$

这样，我们就在不丢失太多精度的情况下，成功避开了中心点的计算雷区，使得整条从外向内的推导线能够安稳降落到终点！

---

## 5. 变量总表（算法相关）

| 符号 | Fortran/Python变量 | 含义 | 单位 | 来源 |
|---|---|---|---|---|
| $\psi_p$ | `PSIPV` | 极向磁通网格值 | Wb | eqdata |
| $\psi_{pn}$ | `PSIPNV` | 归一化极向磁通 | 1 | eqdata |
| $p$ | `PPSI` / `PPPS` | 压强 | Pa（`PPPS`常见以 Pa 存） | 剖面插值 |
| $dp/d\psi_p$ | `DPPSI` | 压强对极向磁通导数 | Pa/Wb | `EQPPSI` |
| $q$ | `QPSI` / `QPV` | 安全因子 | 1 | `EQQPSI` / eqdata |
| $F$ | `FIPV` / `FPSI` / `TTV` | 极向电流函数 $2\pi R B_\phi$ | T·m | 反推/参考 |
| $Y$ | `YP`,`YM` | $\frac{1}{2}F^2$ | (T·m)$^2$ | 中间量 |
| $RSV$ | `RSV` | 等效半径（由 $\psi_t$ 定义） | m | `EQCALV` |
| $VPV$ | `VPV` | `EQCALV` 中通量面体积系数 | m$^3$ | `EQCALV` |
| $AVIR2$ | `AVIR2` | $\langle 1/R^2\rangle$ 相关无量纲化量（`EQCALV`定义） | 1 | `EQCALV` |
| $AVRR2$ | `AVRR2` | `EQCALV` 定义的几何平均量 | m$^2$ | `EQCALV` |
| $X$ | `XP`,`XM` | $RSV/q$ | m | 中间量 |
| $A$ | `ALP`,`ALM` | 几何-磁场耦合系数 | 1/Wb 量纲等效 | 中间量 |
| $B$ | `BLP`,`BLM` | 压强梯度驱动项系数 | (T·m)$^2$/Wb | 中间量 |
| $G$ | `FLP`,`FLM` | 几何组合量 $VPV\cdot AVRR2\cdot X$ | m$^6$（组合量） | 中间量 |
| $R_0$ | `RR` | 参考大半径 | m | 参数 |
| $B_0$ | `BB` | 参考环向场 | T | 参数 |
| $\mu_0$ | `RMU0` | 真空磁导率 | H/m | 常数 |
| $\pi$ | `PI` | 圆周率 | 1 | 常数 |
| 网格点数 | `NRVMAX` | 反推径向网格长度 | 1 | 配置 |
| 网格索引 | `NRV` | 径向索引 | 1 | 循环变量 |
| `QIPV` | `QIPV(NRV)` | `EQIPQP` 循环内实际使用的 `QPSI` 快照 | 1 | `EQIPQP` 内部记录 |
| `DPIPV` | `DPIPV(NRV)` | `EQIPQP` 循环内实际使用的 `DPPSI` 快照 | Pa/Wb | `EQIPQP` 内部记录 |

---

## 6. 直接代数式（`direct`）用于对照

脚本中还计算对照解：

$$
F_{\text{direct}}
=\frac{4\pi^2 q}{(dV/d\psi_p)\,\langle 1/R^2\rangle}
$$

其中：

- $dV/d\psi_p$ 来自 `DVDPSIP`
- $\langle 1/R^2\rangle$ 来自 `AVEIR2`

该式是局部代数关系，不包含 `EQIPQP` 的边界反推积分过程。

---

## 7. 输入与输出（当前脚本与导出链路）

脚本：`in/rebuild_f_from_q_eqdata.py`

### 输入

- `<prefix>_radial_profiles.csv`
- `<prefix>_1D_profiles.csv`
- `<prefix>_parameters.csv`
- `eqgs1d_21_DVDPSIP.csv`
- `eqgs1d_10_AVEIR2.csv`
- `eqgs1d_14_AVEGVR2.csv`
- `eqgs1d_24_EQIPQP_INPUTS.csv`（默认优先使用，除非显式 `--no-eqipqp-input`）

### 关于 `eqgs1d_24_EQIPQP_INPUTS.csv` 的“精确性”

当前实现已将下列量作为 **eqdata 快照** 存储并在导出时优先直接使用：

1. `RSV, VPV, AVIR2, AVRR2, FIPV`
2. `QIPV, DPIPV`（即 `EQIPQP` 循环内每个 `NRV` 实际调用 `EQQPSI/EQPPSI` 的结果）

因此在 `eqgs1d_export.f` 的 plot24 中：

1. 若 `IEQSNAP>=2`，导出的是迭代时的快照输入（不再后处理重算）。
2. 仅当快照不存在时，才回退到 `EQCALV + EQIPQP` 现场重算。

### 输出

- `<out-prefix>.csv`
- `<out-prefix>.png`

CSV 关键列：

- `F_direct`：代数式解
- `F_rebuilt`：选定方法（`direct` 或 `eqipqp`）重建值
- `F_ref`：eqdata 参考 `TTV`
- `dF_direct = F_direct - F_ref`
- `dF = F_rebuilt - F_ref`

---

## 8. 为什么会出现“同公式但结果不同”

`EQIPQP` 对输入定义高度敏感。若 $dp/d\psi_p$、$q$、`RSV/VPV/AVIR2/AVRR2` 的定义、单位、网格任何一项不一致，反推积分误差会沿径向累积并放大。

因此对齐原则是：

1. 同一平衡、同一网格；
2. 同一变量定义（尤其 `EQCALV` 与 `EQCALQ` 量不能直接混用）；
3. 足够导出精度（避免低精度 CSV 带来递推漂移）。

### 本轮定位到的典型错配源

1. 使用 `QPV` 代替 `QPSI`：`EQIPQP` 理论上应使用 `EQQPSI` 返回值，而不是任意外部 `QPV`。
2. `DPPSI` 不是同一链路值：若在导出后再独立重算 `DPPSI`，会与迭代中实际值不同，轴区最敏感。
3. 导出时重跑 `EQCALV`：会导致几何量 (`RSV/VPV/AVIR2/AVRR2`) 与迭代时网格/状态偏离。

### 现状态（已对齐后）

在同一 `prefix` 下，使用快照导出的 `eqgs1d_24_EQIPQP_INPUTS.csv` 进行 Python 重建：

1. `F_rebuilt` 与 Fortran `FIPV` 可达机器精度（RMSE ~ `1e-13`）。
2. `F_rebuilt` 与 `F_ref(TTV)` 也可达机器精度（RMSE ~ `1e-13`）。

---

## 9. 推荐的最小复现实验流程

1. 生成平衡并保存快照（`m5`）：

```bash
printf '0\nc\n' | ./m5 in/g260206.20000_teq_0114 eqdata_modelg5_rebuild_mdleqf9_qmap
```

2. 从该 `eqdata` 导出 plot24（优先使用快照）：

```bash
./eq < in/eq_export_qmap.in
```

3. Python 重建并核验：

```bash
conda run -n Fusion python3 in/rebuild_f_from_q_eqdata.py \
  --prefix eqdata_modelg5_rebuild_mdleqf9_qmap \
  --method eqipqp \
  --out-prefix in/eqdata_modelg5_rebuild_mdleqf9_qmap_f_from_q_eqipqp
```

若输出中 `F vs Fortran FIPV` 与 `F compare vs TTV` 的 RMSE 都在 `1e-12~1e-13` 量级，则链路已对齐。
