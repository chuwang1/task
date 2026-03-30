# `solve_gs_from_gfile.py` 数学推导与直观物理说明

本文档对应脚本：
`/Users/dengxiaoya/TASK/latest/task/eq/in/solve_gs_from_gfile.py`

目标：把脚本中的 Grad-Shafranov (GS) 求解过程用**严谨的数学推导**写清楚，同时穿插**高中物理级别的直观比喻**（“蹦床”、“拔河”、“电子表格”），并和代码实现一一对应。

---

## 1. 物理图像与控制方程

### 1.1 直观比喻：等离子体“蹦床”上的“拔河”
想象一块被紧紧绷在圆环形框架上的**蹦床**。
- **$\psi(R, Z)$（极向磁通）**：就好像蹦床在每个位置 $(R, Z)$ 的**高度或形变**。
- **方程左端（$\Delta^* \psi$）**：代表蹦床本身的**弹性张力**，由于表面张力的存在，它总是试图把形变拉平。
- **方程右端（RHS）**：等离子体内部的**压力**（$p'$）和**极向电流**（$FF'$），就像是在蹦床上放置的重物或者有人在下面拉拽的力。

平衡态（GS方程的解）就是：**蹦床的弹性恢复力**与**外加的压迫力**之间达到了“拔河”式的完美平衡。

### 1.2 严谨数学推导：为什么方程长这样？
托卡马克中的磁场处于轴对称（$\partial/\partial \phi = 0$）状态。根据磁场无散度 $\nabla \cdot \mathbf{B} = 0$，我们可以将磁场表示为极向部分和环向部分的叠加：
$$ \mathbf{B} = \frac{1}{R} \nabla \psi \times \hat{\phi} + \frac{F(\psi)}{R} \hat{\phi} $$
这里 $\psi$ 是极向磁通函数，而 $F(\psi) = R B_\phi$ 代表环向电流分布带来的函数。

根据安培环路定理 $\nabla \times \mathbf{B} = \mu_0 \mathbf{J}$，我们对 $\mathbf{B}$ 取旋度，提取其环向分量 $J_\phi$。
在柱坐标系 $(R, \phi, Z)$ 下，环向分量展开为：
$$ \mu_0 J_\phi = (\nabla \times \mathbf{B})_\phi = -\frac{1}{R} \left( \frac{\partial^2 \psi}{\partial R^2} - \frac{1}{R}\frac{\partial \psi}{\partial R} + \frac{\partial^2 \psi}{\partial Z^2} \right) \equiv -\frac{1}{R} \Delta^* \psi $$
这严谨地定义了特有的 Grad-Shafranov 椭圆算子 $\Delta^*$：
$$ \Delta^* \psi = R \frac{\partial}{\partial R} \left( \frac{1}{R} \frac{\partial \psi}{\partial R} \right) + \frac{\partial^2 \psi}{\partial Z^2} = \frac{\partial^2 \psi}{\partial R^2} - \frac{1}{R}\frac{\partial \psi}{\partial R} + \frac{\partial^2 \psi}{\partial Z^2} $$

另一方面，由理想MHD的力学平衡条件 $\nabla p = \mathbf{J} \times \mathbf{B}$，代入上述 $\mathbf{B}$ 和 $\mathbf{J}$，并在 $\hat{\phi}$ 方向展开投影，可以推导出环向电流密度的另一个物理表达式：
$$ J_\phi = R p'(\psi) + \frac{1}{\mu_0 R} F(\psi)F'(\psi) $$
将两个 $J_\phi$ 的表达式相等，就得到了著名的 **Grad-Shafranov (GS) 方程**：
$$ \Delta^* \psi = -\mu_0 R^2 p'(\psi) - F(\psi)F'(\psi) $$

脚本中的 `RHS` 指的就是方程右端（即上述比喻中的“外加力”项）：
$$ \mathrm{RHS}(R,Z) = -\mu_0 R^2 p'(\psi(R,Z)) - FF'(\psi(R,Z)) $$

**代码实现对应**：
```python
rhs = source_scale * (-MU0 * (R[:, None] ** 2) * pprime - ffprime)
```

---

## 2. 非线性迭代：Picard Iteration（拔河的动态调整）

### 2.1 直观比喻：寻找自洽的受力点
在我们的蹦床比喻中，有一个奇怪的设定：**重物的重量和分布并不是固定的，而是取决于蹦床当前的形变高度 $\psi$**。也就是说，形状改变了，受力也会跟着变。这就是非线性问题的难点。
为了找到最终平衡，我们必须“一步步试探”：
1. 先猜一个蹦床形状。
2. 根据这个形状放上对应分布的重物。
3. 算出在这些重物下蹦床的新形状。
4. 重复这个过程，直到蹦床形状不再变化。这就叫**外迭代 (Outer Loop)**。

### 2.2 严谨数学推导：Picard 迭代格式
GS 方程是一个强非线性椭圆偏微分方程。算子 $\Delta^*$ 是线性的，但右端项 $f(\psi, R) = -\mu_0 R^2 p'(\psi) - FF'(\psi)$ 依赖于未知量 $\psi$。
脚本采用 **Picard 迭代法**（一种定点迭代法）将非线性问题转化为一系列线性子问题求解。

设第 $k$ 步的近似解为 $\psi^{(k)}$，我们构造如下迭代格式来求解第 $k+1$ 步的 $\psi^{(k+1)}$：
$$ \Delta^* \psi^{(k+1)} = -\mu_0 R^2 p'(\psi^{(k)}) - FF'(\psi^{(k)}) $$
即方程右端固定为上一时刻的状态：
$$ \Delta^* \psi^{(k+1)} = \mathrm{RHS}^{(k)}(R, Z) $$

**算法流程（代码中的 `for k in range(max_outer)`）**：
1. 给出初始猜测 $\psi^{(0)}$。
2. 在第 $k$ 步，计算归一化磁通 $\psi_n^{(k)} = \frac{\psi^{(k)}-\psi_{axis}}{\psi_{bdy}-\psi_{axis}}$。
3. 利用 $\psi_n^{(k)}$ 对输入的源项剖面插值，得到空间上的 $p'$ 和 $FF'$，进而组装出 $\mathrm{RHS}^{(k)}$。
4. 此时方程变为带有已知常数右端的线性泊松类方程。求解此线性椭圆偏微分方程得到新的磁通分布 $\psi^{(k+1)}$。
5. 检查残差 $\|\psi^{(k+1)} - \psi^{(k)}\|$，若足够小则迭代收敛，否则令 $k = k+1$ 返回第2步。

---

## 3. 网格离散推导：将连续空间转化为“电子表格” (Finite Difference)

### 3.1 直观比喻：网格化的“电子表格”
计算机无法直接理解连续的曲面。我们需要把 $R-Z$ 平面划分为很多小网格，就像 Excel 的**二维电子表格**。
在这个电子表格中，某一个单元格的值 $\psi_{i,j}$，通过差分公式被设定为只与其**上下左右**四个相邻单元格的值，以及该点本身受到的外部力 $\mathrm{RHS}_{i,j}$ 相关联。通过公式不断自我刷新，整张表最终会达到稳定态。

### 3.2 严谨数学推导：中心差分法的泰勒展开
设我们有步长为 $\Delta R$ 和 $\Delta Z$ 的均匀网格：
$$ R_i = R_{min} + i\Delta R, \quad Z_j = Z_{min} + j\Delta Z $$

我们需要将微分算子 $\Delta^*$ 在点 $(i, j)$ 处进行有限差分离散。我们利用泰勒展开(Taylor Expansion)：
$$ \psi(R+\Delta R) = \psi(R) + \Delta R \frac{\partial \psi}{\partial R} + \frac{\Delta R^2}{2} \frac{\partial^2 \psi}{\partial R^2} + O(\Delta R^3) $$
$$ \psi(R-\Delta R) = \psi(R) - \Delta R \frac{\partial \psi}{\partial R} + \frac{\Delta R^2}{2} \frac{\partial^2 \psi}{\partial R^2} - O(\Delta R^3) $$

将上面两式相减并略去高阶项，得到一阶导数的**中心差分**：
$$ \frac{\partial \psi}{\partial R} \approx \frac{\psi_{i+1,j} - \psi_{i-1,j}}{2\Delta R} $$
将上面两式相加，得到二阶导数的**中心差分**：
$$ \frac{\partial^2 \psi}{\partial R^2} \approx \frac{\psi_{i+1,j} - 2\psi_{i,j} + \psi_{i-1,j}}{\Delta R^2} $$
同理可得 $Z$ 方向的离散：
$$ \frac{\partial^2 \psi}{\partial Z^2} \approx \frac{\psi_{i,j+1} - 2\psi_{i,j} + \psi_{i,j-1}}{\Delta Z^2} $$

将所有离散式代入方程 $\Delta^*\psi = \mathrm{RHS}$：
$$ \left[\frac{\psi_{i+1,j} - 2\psi_{i,j} + \psi_{i-1,j}}{\Delta R^2}\right] - \frac{1}{R_i} \left[\frac{\psi_{i+1,j} - \psi_{i-1,j}}{2\Delta R}\right] + \left[\frac{\psi_{i,j+1} - 2\psi_{i,j} + \psi_{i,j-1}}{\Delta Z^2}\right] = \mathrm{RHS}_{i,j} $$

重新排列组合，合并同类项，归纳出各自的空间系数：
*   左侧 $\psi_{i-1,j}$ 系数： $A_i = \frac{1}{\Delta R^2} + \frac{1}{2R_i\Delta R}$
*   右侧 $\psi_{i+1,j}$ 系数： $B_i = \frac{1}{\Delta R^2} - \frac{1}{2R_i\Delta R}$
*   上下 $\psi_{i,j-1}$ 及 $\psi_{i,j+1}$ 系数： $C = \frac{1}{\Delta Z^2}$
*   中心点 $\psi_{i,j}$ 系数： $-(\frac{2}{\Delta R^2} + \frac{2}{\Delta Z^2}) = -(A_i + B_i + 2C)$

我们把总权重记作 $D_i = A_i + B_i + 2C$。将包含周围邻点的项移到等号另一边，我们可以显式表达出中心点的**理论更新候选值 (Candidate)**：
$$ \psi_{i,j}^{cand} = \frac{A_i\psi_{i-1,j} + B_i\psi_{i+1,j} + C\psi_{i,j-1} + C\psi_{i,j+1} - \mathrm{RHS}_{i,j}}{D_i} $$

代码中的 `A`, `B`, `C`, `D` 矩阵以及 `cand` 变量计算直接严格对应了这个数学推导。

---

## 4. 线性子问题求解：加权 Jacobi (Weighted Jacobi)

在计算出上一步的候选值 $\psi_{i,j}^{cand}$ 后，为了求解上面庞大的联立方程组，脚本使用了一种迭代求解器——加权 Jacobi 方法（内层迭代 Inner loop）。

$$ \psi_{i,j}^{new} = (1-\omega)\psi_{i,j}^{old} + \omega \psi_{i,j}^{cand}, \quad 0 < \omega \le 1 $$

*   如果 $\omega = 1$（标准 Jacobi），网格点每次完全听从新算出的候选值。
*   如果 $\omega < 1$（**欠松弛 Under-relaxation**），网格点每次只向新候选值方向移动一小步。这在物理上相当于加入了“阻尼”，能有效避免电子表格在刷新过程中发生数值震荡，以略微牺牲速度的代价换取极高的稳定性。

每次迭代结束后，脚本评估 `step_change = max(|new - old|)`，如果 `step_change < tol_inner`，即判定这个线性子问题已经解出了稳定结果。

---

## 5. 边界条件与区域掩膜

如果“蹦床”没有坚固的外框固定，就无法产生张力。这里的“外框固定”就是**边界条件 (Boundary Conditions)**。脚本处理边界的方式有：

1. `rectangle` 模式：
   - 最简单的 Dirichlet 边界。整个矩形计算域四周被牢牢“钉死”，值被固定为输入文件（gfile）指定的背景通量。

2. `lcfs` (Ring Dirichlet) 模式：
   - 真实托卡马克的约束区域只在最后闭合磁通面 (LCFS) 内部。
   - 此模式把真正的边界强行定在多边形 LCFS 上，使其值为给定的 $\psi_{bdy}$。LCFS 之外的区域直接停止更新。

3. `lcfs_curve` (Ghost-cell / Cut-cell) 模式：
   - 数学更严谨。物理上的连续 LCFS 曲线几乎不可能完美落在离散的二维网格点上。
   - 假设网格点 $P$ 位于域内，其相邻点 $W$ 位于域外。真正的边界 $\psi_b$ 落在线段 $PW$ 之间。
   - 设边界点到 $P$ 的距离占比为 $1-\lambda$，到 $W$ 为 $\lambda$。利用线性插值近似：
     $$ \psi_b \approx (1-\lambda)\psi_P + \lambda \psi_W^{ghost} $$
   - 我们倒推出域外那个点如果存在的话，它必须具备什么值才能满足插值曲线条件（即 **幽灵点 Ghost Cell**）：
     $$ \psi_W^{ghost} = \frac{\psi_b - (1-\lambda)\psi_P}{\lambda} $$
   - 在进行空间差分（电子表格刷新）时，当探测到邻居在域外，就代入这个 $ghost$ 值，从而让等离子体边界完美契合复杂的几何曲线。

---

## 6. 一句话总结

`solve_gs_from_gfile.py` 的核心逻辑是：
为了求解复杂的等离子体磁流体平衡，脚本首先通过**泰勒展开**将 Grad-Shafranov 偏微分算子转化为网格上的**有限差分系统**（电子表格）；面对由于 $p'$ 和 $FF'$ 改变带来的非线性，使用 **Picard迭代**不断试探并更新源项（拔河的重量动态调整）；而在每次组装好线性子问题后，应用带有阻尼的**加权 Jacobi 方法**刷出平滑解；最终配合**幽灵点插值边界**，求出严谨自洽的极向磁通分布 $\psi(R,Z)$。