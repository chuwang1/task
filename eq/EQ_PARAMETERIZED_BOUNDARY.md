# EQ 参数化边界说明

## 1. 结论先行

在 `TASK/EQ` 中，最外闭合磁面（LCFS）在固定边界流程下是显式参数化的，核心几何参数是：

- `RR`：大半径（几何中心）
- `RA`：小半径（水平尺度）
- `RKAP`：elongation（拉伸）
- `RDLT`：triangularity（三角形度）

其中 `BB`、`RIP` 主要控制磁场/电流和通量量级，`RB` 主要控制壁面或计算域范围，不直接作为 LCFS 参数化方程的形状参数。

---

## 2. 参数化方程

等价的边界参数形式可写成（$\tau$ 为参数角）：

$$
R(\tau)=RR+RA\cos(\tau+RDLT\sin\tau),\qquad
Z(\tau)=RA\cdot RKAP\sin\tau
$$

该形式对应到源码中的实现路径：

1. 通过 `EQFBND` 求 $THDASH=\tau$（Brent 根求解）  
   `eqcalc.f` 中：
   - `EQFBND = ZBRF*COS(X+RDLT*SIN(X)) - RKAP*SIN(X)`  
   - 见 `EQFBND` 与 `EQDEFB`：`eqcalc.f` 行 170-271 附近

2. 得到边界半径 `RHOM/RHOG`  
   $$
   \rho_b = RA\sqrt{\cos^2(\tau+RDLT\sin\tau)+RKAP^2\sin^2\tau}
   $$
   - 见 `eqcalc.f` 行 187-193

3. 映射到几何坐标  
   $$
   R=RR+\rho_b\cos\theta,\qquad Z=\rho_b\sin\theta
   $$
   在边界上（$\sigma=1$）与上面的参数式等价。
   - 见 `eqcalc.f` 行 216-230

---

## 3. 固定边界与自由边界中的角色

## 3.1 固定边界（常用 `MODELG=2`, 菜单 `R` -> `EQCALC`）

固定边界主流程：

1. `EQMESH`
2. `EQPSIN/EQPSIR`（初值）
3. `EQDEFB`（由 `RR,RA,RKAP,RDLT` 定义边界）
4. `EQLOOP`（迭代求解）
5. `EQTORZ`、`EQCALP`

入口见 `eqcalc.f` 行 7-23。  
因此在不走自由边界时，这四个形状参数仍是核心输入。

## 3.2 自由边界（`EQCALX`, `MDLEQX`）

`MDLEQX=1/2` 时，会用 Newton 过程调 `PSIB`，使解出的几何量逼近目标 `RR,RA,RKAP,RDLT`：

- 残差定义：$RSLT(1:4)=(RRR-RR,\;RRA-RA,\;RRK-RKAP,\;RRD-RDLT)$
- 见 `eqcalx.f` 行 195-200

即：自由边界下，这四个参数是“目标形状约束”；固定边界下，它们是“直接参数化边界”。

---

## 4. 各参数对形状的几何影响

- 增大 `RR`：整体向外平移（`R` 方向）
- 增大 `RA`：整体放大截面（同时放大 `R`、`Z` 尺度）
- 增大 `RKAP`：上下拉长（`Z` 向伸长）
- 增大 `RDLT`：上/下顶点向内偏，形成 D 形截面

---

## 5. 与本仓库 Python 脚本的一致性

脚本：

- `in/calc_lcfs_shape.py`

该脚本直接使用等价参数方程采样生成 LCFS 点，便于快速画图与导出，不改变 `eq` 求解器逻辑。

示例（你给的参数）：

- $RR=8.03$
- $RA=2.72$
- $RKAP=1.89$
- $RDLT=0.59$

运行命令（Fusion 环境）：

```bash
conda run -n Fusion python in/calc_lcfs_shape.py \
  --RR 8.03 --RA 2.72 --RKAP 1.89 --RDLT 0.59 \
  --csv in/lcfs_from_user_params.csv \
  --plot in/lcfs_from_user_params.png
```

对应输出：

- `in/lcfs_from_user_params.csv`
- `in/lcfs_from_user_params.png`

---

## 6. 参数化边界不包含的内容

下列物理量不由边界参数化方程直接决定：

- 压强/电流剖面（`PP*`, `PJ*`, `FF*` 及其 `PROF*`）
- 全局电流和磁场量级（`RIP`, `BB`）
- 真空区外推与壁内外网格细节（`RB`, `FRBIN`, `MDLEQV` 等）

它们会改变平衡解和通量分布，但不直接改写上面的 LCFS 参数化几何公式。
