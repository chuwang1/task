# ONETWO 回旋辐射 PRC 的 Fortran 接入与验证方案

## 1. 目标

本文档描述如何将 PyMak 中的 ONETWO-style 回旋/同步辐射解析模型移植到 TR Fortran 代码中，并在正式接入输运演化前，先用一个独立 Fortran 小模块与 PyMak 进行同输入验证。

核心原则是先验证公式和单位，再接入 `PRC` 源项：

1. 用 Fortran 实现一个纯函数式 ONETWO PRC 计算模块。
2. 输入与 PyMak 完全一致的 `ne, Te, B, a, R0, wall_reflection`。
3. 对比 Fortran 与 `PyMak/cyclotron_radiation.py` 的 `calculate_cyclotron_radiation_onetwo` 输出。
4. 确认误差达到浮点舍入级别后，再把该模块接到 TR 的 `MDLPR` 分支中。

ONETWO 解析模型本身不需要 CYTRAN 的频率网格、偏振模式、体积单元和边界面积数组。它只依赖局部温度密度剖面、整体几何半径、中心磁场和壁反射率。

## 2. 现有 TR 中 PRC 的位置

当前 Fortran 代码中，回旋辐射损失 `PRC` 的计算路径是：

- `trcalc.f90` 中调用 `TR_CYTRAN`。
- `trcytran.f90` 中把 `RN, RT, ABB1RHO, PVOLRHOG, PSURRHOG` 等量传给 `CYTRAN`。
- `TR_CYTRAN` 得到 `psync_rm [keV/m^3/s]` 后转换为 `PRC [W/m^3]`。
- `tradat.f90` 的 `TRLOSS` 中把 `PRC` 加入 `PRSUM`。

目前关键逻辑是：

```fortran
IF(MDLPR.GT.0) CALL TR_CYTRAN
```

以及：

```fortran
SELECT CASE(MDLPR)
CASE(0)
   PRSUM(NR)=PRL(NR)+PRB(NR)
CASE(1,2)
   PRSUM(NR)=PRL(NR)+PRB(NR)+PRC(NR)
END SELECT
```

这说明 `TRLOSS` 已经允许 `MDLPR=2` 包含 `PRC`。因此建议保留现有 `MDLPR=1` 作为 CYTRAN，新增 `MDLPR=2` 作为 ONETWO。

推荐最终分支为：

```fortran
SELECT CASE(MDLPR)
CASE(1)
   CALL TR_CYTRAN
CASE(2)
   CALL TR_ONETWO_PRC
END SELECT
```

## 3. PyMak ONETWO 模型输入输出

PyMak 中对应函数为：

```python
calculate_cyclotron_radiation_onetwo(
    geometry,
    ne,
    Te,
    bt_T=None,
    wall_reflection=0.8,
    minor_radius_m=None,
    major_radius_m=None,
)
```

主要输入单位：

| 量 | PyMak 输入 | Fortran/TR 对应量 | 单位 |
| --- | --- | --- | --- |
| `ne` | `ne` | `RN(NR,1)` | `1e20 / m^3` |
| `Te` | `Te` | `RT(NR,1)` | `keV` |
| `bt_T` | `bt_T` | `BB`，或后续评估 `ABB1RHO(NR)` | `T` |
| `minor_radius_m` | `a` | `RA` | `m` |
| `major_radius_m` | `R0` | `RR` | `m` |
| `wall_reflection` | `wall_reflection` | 建议新增 `REFRAD` | 无量纲 |

PyMak 输出：

| 量 | 单位 |
| --- | --- |
| `power_density_total_MW_m3` | `MW/m^3` |
| `power_density_total_keV_m3_s` | `keV/m^3/s` |
| `phi_bar` | 无量纲 |

TR 内部 `PRC` 单位是 `W/m^3`。因此 Fortran 接入 TR 时必须保存：

```fortran
PRC(NR) = qsync_w_m3
```

不要把 `MW/m^3` 直接写入 `PRC`。CSV 和图形输出阶段会自行乘 `1.D-6`。

## 4. ONETWO 解析公式

PyMak 当前实现按 ONETWO 源码的 Gaussian-cgs 单位执行。移植时建议 Fortran 也直接使用相同 cgs 常数，避免 SI 与 cgs 混合造成隐藏系数错误。

常数：

```text
e_cgs      = 4.803204712570263e-10 statC
m_e_cgs    = 9.1093837015e-28 g
c_cgs      = 2.99792458e10 cm/s
kev_to_erg = 1.602176634e-9 erg/keV
```

输入转换：

```text
B_gauss  = abs(B_T) * 1e4
a_cm     = a_m * 100
R0_cm    = R0_m * 100
ne_cm3   = ne_1e20_m3 * 1e14
Te_erg   = Te_keV * kev_to_erg
mec2_erg = m_e_cgs * c_cgs^2
```

频率：

```text
omega_ce = e_cgs * B_gauss / (m_e_cgs * c_cgs)
omega_pe = sqrt(4*pi*ne_cm3*e_cgs^2/m_e_cgs)
```

几何修正：

```text
rm  = a_cm / R0_cm
chi = rm * sqrt(mec2_erg / Te_erg)
```

透明度因子：

```text
phi_bar = 60 * (Te_erg / mec2_erg)^1.5
        * sqrt(c_cgs*omega_ce/(a_cm*omega_pe^2)
               * (1 - wall_reflection)
               * (1 + chi))
```

单粒子加速度项：

```text
vdotsq = omega_ce^2 * 2 * Te_erg / m_e_cgs
```

功率密度：

```text
qsync_kev_cm3_s = ne_cm3/1.5 * e_cgs^2/c_cgs^3 * vdotsq * phi_bar / kev_to_erg
qsync_kev_m3_s  = qsync_kev_cm3_s * 1e6
qsync_w_m3      = qsync_kev_m3_s * 1.602176634e-16
qsync_mw_m3     = qsync_w_m3 * 1e-6
```

这里 `1.602176634e-16` 是 `keV -> J` 的转换因子。

## 5. 为什么先写独立 Fortran 小模块

ONETWO 是局域解析模型，它不像 CYTRAN 那样依赖 TR 的完整几何网格和辐射输运上下文。因此在正式接入前，可以把问题拆成一个非常小的、可独立验证的数值模块。

这样做有三个好处：

1. 把公式移植错误和 TR 演化耦合错误分离。
2. 可以用同一组温度密度数组直接对比 PyMak 和 Fortran。
3. 验证通过后，接入 TR 只剩下输入映射和 `PRC` 单位转换问题。

建议新增一个纯计算模块，例如：

```text
tr_onetwo_prc_kernel.f90
```

模块只暴露两个接口：

```fortran
MODULE tr_onetwo_prc_kernel
  USE bpsd_kinds
  IMPLICIT NONE
CONTAINS

  PURE SUBROUTINE onetwo_prc_profile( &
       nrmax, ne_1e20_m3, te_kev, bt_t, a_m, r0_m, wall_reflection, &
       prc_w_m3, phi_bar)

  PURE REAL(KIND=rkind) FUNCTION onetwo_prc_point( &
       ne_1e20_m3, te_kev, bt_t, a_m, r0_m, wall_reflection, phi_out)

END MODULE tr_onetwo_prc_kernel
```

其中：

- `onetwo_prc_point` 只算单个径向点。
- `onetwo_prc_profile` 对数组循环，调用 `onetwo_prc_point`。
- `prc_w_m3` 输出 TR 内部单位。
- `phi_bar` 可选输出，用于和 PyMak 诊断量对比。

模块内部不要 `USE trcomm`，不要直接读取 `RN, RT, RR, RA, BB`。这样它可以独立编译和单元测试。

## 6. 独立验证输入格式

建议先用一个简单 CSV 作为 Fortran 和 PyMak 的共同输入，例如：

```csv
rho,ne_1e20_m3,te_kev
0.00,1.20,25.0
0.25,1.05,22.0
0.50,0.80,16.0
0.75,0.45,8.0
1.00,0.15,1.0
```

几何和壁反射率使用单独参数：

```text
bt_T            = 6.0
minor_radius_m  = 2.72
major_radius_m  = 8.03
wall_reflection = 0.8
```

Fortran 验证程序输出：

```csv
rho,ne_1e20_m3,te_kev,prc_w_m3,prc_mw_m3,phi_bar
```

PyMak 验证脚本读取同一输入，调用：

```python
calculate_cyclotron_radiation_onetwo(
    geometry,
    ne=ne,
    Te=te,
    bt_T=6.0,
    minor_radius_m=2.72,
    major_radius_m=8.03,
    wall_reflection=0.8,
)
```

然后输出同样字段。

## 7. 验证判据

推荐比较：

```text
abs_err = abs(prc_fortran - prc_pymak)
rel_err = abs_err / max(abs(prc_pymak), eps)
```

通过标准：

```text
max_rel_err(PRC)     < 1e-12 到 1e-10
max_rel_err(phi_bar) < 1e-12 到 1e-10
```

如果使用不同编译器优化或输出位数较少，可以放宽到 `1e-9`。若误差达到百分比量级，通常说明以下之一有问题：

- `ne` 的 `1e20/m^3 -> cm^-3` 转换错了。
- `keV -> erg` 或 `keV -> J` 转换错了。
- `B[T] -> Gauss` 转换错了。
- `wall_reflection` 被误当成了 `SYNCABS`。
- `PRC` 被写成了 `MW/m^3` 而非 `W/m^3`。

## 8. 壁反射率参数建议

不要直接复用 `SYNCABS` 作为 ONETWO 的 `wall_reflection`。

当前 CYTRAN 参数含义是：

```text
SYNCABS  = fraction of cyclotron radiation absorption by walls
SYNCSELF = fraction of x/o mode reflected as same mode
```

CYTRAN 中壁反射强度近似为：

```text
1 - SYNCABS
```

而 ONETWO 公式直接使用：

```text
wall_reflection
```

并通过：

```text
sqrt(1 - wall_reflection)
```

控制净逃逸损失。

如果把 `wall_reflection` 和 `SYNCABS` 混用，容易出现物理语义反转。例如当前输入文件中 `SYNCABS=0.0`，如果令 `wall_reflection=1-SYNCABS=1.0`，ONETWO 净损失会趋近于 0。

因此建议新增 namelist 参数：

```fortran
REFRAD = 0.8D0
```

含义：

```text
REFRAD = wall reflection coefficient for ONETWO PRC
```

`REFRAD` 只用于 `MDLPR=2`。`SYNCABS/SYNCSELF` 继续只服务 CYTRAN。

## 9. 正式接入 TR 的步骤

验证小模块通过后，再做正式接入。

### 9.1 增加 TR 包装子程序

新增包装层：

```text
tr_onetwo_prc.f90
```

包装层可以 `USE trcomm`，负责把 TR 内部变量映射到纯计算模块：

```fortran
SUBROUTINE TR_ONETWO_PRC
  USE trcomm, ONLY: nrmax, rn, rt, prc, rr, ra, bb, refrad
  USE tr_onetwo_prc_kernel, ONLY: onetwo_prc_profile

  CALL onetwo_prc_profile( &
       nrmax, rn(1:nrmax,1), rt(1:nrmax,1), bb, ra, rr, refrad, &
       prc(1:nrmax), phi_dummy)
END SUBROUTINE TR_ONETWO_PRC
```

初版建议使用 `BB` 作为磁场输入，与 PyMak 的 `bt_T`/`geometry.bcentr` 对齐。后续若希望引入径向磁场效应，可再评估改为 `ABB1RHO(NR)`，但这会偏离当前 PyMak ONETWO 输入定义。

### 9.2 修改 `trcalc.f90`

把：

```fortran
IF(MDLPR.GT.0) CALL TR_CYTRAN
```

改为：

```fortran
SELECT CASE(MDLPR)
CASE(1)
   CALL TR_CYTRAN
CASE(2)
   CALL TR_ONETWO_PRC
END SELECT
```

### 9.3 增加参数

需要修改：

- `trcomm.f90`：增加 `REAL(rkind):: REFRAD`
- `trinit.f90`：设置默认值 `REFRAD=0.8D0`
- `trparm.f90`：加入 namelist `/TR/`
- `trinit.f90` 注释：说明 `MDLPR=2` 是 ONETWO

建议注释：

```text
MDLPR:
  0: PRSUM = PRB + PRL
  1: PRSUM = PRB + PRL + PRC, PRC from CYTRAN
  2: PRSUM = PRB + PRL + PRC, PRC from ONETWO analytic model

REFRAD:
  wall reflection coefficient used by ONETWO PRC
```

### 9.4 修改 Makefile

若新增独立文件，需要把源文件加入 `SRCS`：

```make
tr_onetwo_prc_kernel.f90 tr_onetwo_prc.f90
```

并添加必要依赖关系。

为了减少 Makefile 改动，也可以先把 `TR_ONETWO_PRC` 放进现有 `trcytran.f90` 模块。但长期看，ONETWO 和 CYTRAN 物理模型不同，单独文件更清楚。

## 10. TR 级别验证

正式接入后，推荐做两层验证。

第一层：冻结剖面验证。

- 在同一时刻读取 `RN(:,1), RT(:,1)`。
- 用 Fortran `TR_ONETWO_PRC` 生成 `PRC`。
- 用 PyMak 读取同一剖面生成 `PRC`。
- 比较 `PRC(NR)`。

第二层：完整输运运行验证。

- 复制输入文件，设置 `MDLPR=2, REFRAD=0.8D0`。
- 运行 `./tr2 < input`。
- 对比输出 CSV 中的 `PRC` 剖面。
- 与 `MDLPR=1` 的 CYTRAN 结果比较量级和趋势，但不要求二者相等，因为物理模型不同。

注意：当前 `TR_WRITE_CSV` 使用 `STATUS='REPLACE'` 写 `tr_data_001.csv` 等文件。完整运行前应在临时目录运行，或先备份已有 CSV。

## 11. 数值保护

Fortran 实现中建议加入以下保护：

```text
if ne <= 0      -> PRC = 0
if Te <= 0      -> PRC = 0
if abs(B) <= 0  -> PRC = 0
if a <= 0       -> PRC = 0
if R0 <= 0      -> PRC = 0
wall_reflection = min(max(wall_reflection, 0), 1)
```

对分母使用小正数保护：

```text
max(Te_erg, 1e-30)
max(a_cm * omega_pe^2, 1e-30)
max(R0_cm, 1e-30)
```

但不要用保护掩盖正常输入错误。验证阶段若发现 `Te<=0` 或 `ne<=0`，应输出诊断。

## 12. 推荐实施顺序

1. 写 `tr_onetwo_prc_kernel.f90`，只实现公式。
2. 写一个独立 test driver，读取小 CSV，输出 Fortran PRC。
3. 写 PyMak 对照脚本，读取同一 CSV，输出 PyMak PRC。
4. 比较 `PRC` 和 `phi_bar`，直到逐点一致。
5. 增加 TR 包装子程序 `TR_ONETWO_PRC`。
6. 用 `MDLPR=2` 接入 `trcalc.f90`。
7. 增加 `REFRAD` namelist 参数。
8. 编译并跑完整 TR 验证。

这一路径把移植风险分成两个独立问题：公式数值一致性，以及 TR 演化耦合一致性。前者通过小模块解决，后者再通过 `MDLPR=2` 的完整运行验证。
