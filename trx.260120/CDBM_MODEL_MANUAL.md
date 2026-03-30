# CDBM Transport Model Technical Manual

## 1. Overview

The **Current Diffusivity Ballooning Mode (CDBM)** model is the primary anomalous transport model used in the TR code for tokamak plasma simulations. It calculates turbulent thermal diffusivities (`chi_e`, `chi_i`) based on the theory that micro-tearing and ballooning instabilities, driven by current density gradients and pressure gradients, dominate the anomalous transport.

The model is implemented in the subroutine `TRCFDW` in `trcoef.f90`, selected by the input parameter `MDLKAI`.

### MDLKAI ranges

| MDLKAI | Description |
|--------|-------------|
| 30–39 | Inline CDBM with `|alpha|^{3/2}` scaling |
| 40–50 | Inline CDBM with `|alpha|` scaling and electron thermal velocity correction |
| 130–139 | External `tr_cdbm` subroutine (recommended for production) |

---

## 2. Physical Background

### 2.1 Core CDBM Formula

The thermal diffusivity in CDBM is expressed as:

```
chi = C_K * F(s, alpha, kappa_q) * |alpha|^{3/2} * delta_e^2 * V_A / (q * R)
```

where:

| Symbol | Definition | Unit |
|--------|-----------|------|
| `C_K` | Numerical coefficient (`CK0` for electrons, `CK1` for ions) | — |
| `F(s, alpha, kappa_q)` | Structure function (form factor) | — |
| `alpha` | Normalized pressure gradient | — |
| `delta_e` | Collisionless electron skin depth | m |
| `V_A` | Alfvén velocity | m/s |
| `q` | Safety factor | — |
| `R` | Major radius | m |

### 2.2 Key Intermediate Quantities

**Alfvén velocity:**

```
V_A = B / sqrt(mu_0 * n_e * m_i)
```

where `m_i = (m_D * n_D + m_T * n_T + m_He * n_He) / (n_D + n_T + n_He)` is the effective ion mass.

**Collisionless skin depth:**

```
delta_e^2 = c^2 / omega_{pe}^2 = c^2 * m_e * epsilon_0 / (n_e * e^2)
```

**Normalized pressure gradient (Shafranov alpha):**

```
alpha = -2 * mu_0 * q^2 * R / B^2 * dp/dr
```

where `dp/dr` is the total (ion + electron + fast particle) pressure gradient in units of `[10^{20} m^{-3} * keV / m]`, converted to `[Pa/m]` by multiplying `10^{20} * e_{keV}`.

**Magnetic shear:**

```
s = (rho / q) * dq/dr
```

computed from the `q` profile using a 3-point numerical derivative.

**Magnetic curvature:**

```
kappa_q = -epsilon * (1 - 1/q^2)
```

where `epsilon = r/R` is the local inverse aspect ratio.

---

## 3. Structure Functions

The CDBM model provides multiple structure (form) functions. The choice is controlled by the last digit of `MDLKAI`.

### 3.1 TRCOFS(s, alpha, kappa_q) — Standard form factor

Used by `MDLKAI = 31, 32, 33, 34, 41, 42, 43, 44`.

**Algorithm:**

Define the shifted shear:

```
If alpha >= 0:
    sigma = s - alpha
    (curvature branch uses kappa_q > 0)
Else:
    sigma = alpha - s
    (curvature branch uses kappa_q < 0)
```

Compute `F_1`:

```
If sigma >= 0:
    F_1 = (1 + 9*sqrt(2) * sigma^{5/2}) / (sqrt(2) * (1 - 2*sigma + 3*sigma^2 + 2*sigma^3))
Else:
    F_1 = 1 / sqrt(2 * (1 - 2*sigma) * (1 - 2*sigma + 3*sigma^2))
```

Compute `F_2` (curvature-driven term):

```
If alpha >= 0 and kappa_q > 0:
    F_2 = sqrt(kappa_q)^3 / s^2
Else if alpha < 0 and kappa_q < 0:
    F_2 = sqrt(-kappa_q)^3 / s^2
Else:
    F_2 = 0
```

**Result:**

```
TRCOFS = max(F_1, F_2)
```

> **Note:** The `1/s^2` term in `F_2` causes a singularity as `s -> 0`, which is the primary source of the "chi spike" near the magnetic axis.

### 3.2 TRCOFSS(s, alpha) — Simplified form factor

Used by `MDLKAI = 35, 36, 37, 38`.

```
sigma = s - alpha
F = 2 * sigma^2 / (1 + (2/9) * |sigma|^{5/2})
```

This form avoids the `1/s^2` singularity.

### 3.3 TRCOFSX(s, alpha, kappa_q, epsilon_a) — Extended form factor

Used by `MDLKAI = 39, 49`.

Like `TRCOFS` but includes:

- Shafranov shift correction: `sigma = s - (1 - 2*alpha/(1 + 6*alpha)) * alpha`
- Elongation modification: `F_1` multiplied by `(1 + kappa_q)^{5/2}`
- Curvature term uses `sqrt(kappa_q / epsilon_a)^3 / s^2`

### 3.4 TRCOFT(s, alpha, kappa_q, epsilon_a) — Alternative form factor

Used by `MDLKAI = 50`.

Identical to `TRCOFS` except the curvature term uses:

```
F_2 = sqrt(kappa_q / epsilon_a)^3 / s^2
```

---

## 4. Inline CDBM Models (MDLKAI = 30–50)

### 4.1 CASE(30:39) — Standard CDBM

**General formula:**

```
chi_e = CK0 * F(s, alpha, kappa_q) * |alpha|^{3/2} * delta_e^2 * V_A / (q * R)
chi_i = CK1 * F(s, alpha, kappa_q) * |alpha|^{3/2} * delta_e^2 * V_A / (q * R)
```

**Sub-models:**

| MDLKAI | Form factor `F` | ExB suppression | alpha in F |
|--------|-----------------|-----------------|------------|
| 30 | `1/(1.7 + sqrt(6)*s)` | No | — |
| 31 | `TRCOFS(s, CALF*alpha, kappa_q)` | No | Yes |
| 32 | `TRCOFS(s, CALF*alpha, kappa_q) * f_{ExB}` | Yes | Yes |
| 33 | `TRCOFS(s, 0, kappa_q)` | No | No |
| 34 | `TRCOFS(s, 0, kappa_q) * f_{ExB}` | Yes | No |
| 35 | `TRCOFSS(s, CALF*alpha)` | No | Yes |
| 36 | `TRCOFSS(s, CALF*alpha) * f_{ExB}` | Yes | Yes |
| 37 | `TRCOFSS(s, 0)` | No | No |
| 38 | `TRCOFSS(s, 0) * f_{ExB}` | Yes | No |
| 39 | `TRCOFSX(s, CALF*alpha, kappa_q, a/R)` | No | Yes |

**CDBM05 elongation correction** (when `MDLCD05 != 0`, applicable to MDLKAI = 31, 32):

```
F -> F * (2*sqrt(kappa) / (1 + kappa^2))^{3/2}
```

where `kappa = RKPRHO(NR)` is the local elongation from the equilibrium.

### 4.2 CASE(40:50) — Drift-wave ballooning

Uses `|alpha|` instead of `|alpha|^{3/2}`, plus an electron thermal velocity factor:

```
chi = CK * F * |alpha| * delta_e^2 * V_A / (q * R) * (v_{Te} / V_A)
```

where `v_{Te} = sqrt(T_e * e_{keV} / m_e)`.

MDLKAI = 50 additionally includes diamagnetic stabilization via the modified Bessel function `Lambda(lambda)`:

```
chi_e = CK0 * F * |alpha| * delta_e^2 * V_A / (q*R) / (Lambda * (1 + omega_{*T}^2))
chi_i = CK1 * F * |alpha| * delta_e^2 * V_A / (q*R) / (1 + omega_{*T}^2)
```

---

## 5. External CDBM Subroutine (MDLKAI = 130–139)

### 5.1 Overview

For `MDLKAI = 130–139`, the chi calculation is delegated to the external subroutine `tr_cdbm`:

```fortran
MODEL = MDLKAI - 130
CALL tr_cdbm(BB, RR, RS, RKAPL, QL, SHEARL, PNEL, RHONI, DPDRL,
             DVEXBDRL, CALF, CKAP, CEXB, MODEL, chi_cdbm)
AKDWEL = (CK0/12.0) * chi_cdbm
AKDWIL = (CK1/12.0) * chi_cdbm
```

**Sub-models:**

| MDLKAI | MODEL | Description |
|--------|-------|-------------|
| 130 | 0 | CDBM original |
| 131 | 1 | CDBM05 (elongation correction) |
| 132 | 2 | CDBM + weak ExB shear |
| 133 | 3 | CDBM05 + weak ExB shear |
| 134 | 4 | CDBM + strong ExB shear |

### 5.2 ExB shear calculation

Before calling `tr_cdbm`, the ExB shear suppression factor is computed:

```fortran
SL    = s^2 + 0.1^2
WE1   = -q * R / (SL * V_A) * dV_ExB/dr
CEXB  = CWEB * FEXB(|WE1|, s, alpha)
```

### 5.3 Input/Output of `tr_cdbm`

**Inputs:**

| Variable | Meaning | Unit |
|----------|---------|------|
| `BB` | Toroidal magnetic field | T |
| `RR` | Major radius | m |
| `RS` | Local minor radius `= a * rho` | m |
| `RKAPL` | Local elongation | — |
| `QL` | Local safety factor | — |
| `SHEARL` | Magnetic shear (possibly smoothed) | — |
| `PNEL` | Electron density `= n_e * 10^{20}` | m^{-3} |
| `RHONI` | Ion mass density | kg/m^3 |
| `DPDRL` | Pressure gradient | Pa/m |
| `DVEXBDRL` | ExB velocity gradient | 1/s |
| `CALF` | alpha scaling factor | — |
| `CKAP` | Elongation factor (= 1.0) | — |
| `CEXB` | ExB suppression factor from `FEXB` | — |
| `MODEL` | Sub-model index (0–4) | — |

**Output:**

| Variable | Meaning | Unit |
|----------|---------|------|
| `chi_cdbm` | Raw CDBM thermal diffusivity | m^2/s |

The raw `chi_cdbm` is then scaled: `chi_e = (CK0/12) * chi_cdbm`, `chi_i = (CK1/12) * chi_cdbm`.

---

## 6. ExB Shear Suppression: FEXB

The function `FEXB(x, s, alpha)` returns a suppression factor `f_{ExB} in [0, 1]`:

```
f_{ExB} = exp(-beta * x^gamma)
```

where `x = |omega_{ExB}|` is the normalized ExB shearing rate.

**beta:**

```
alpha_L = max(|alpha|, 0.001)
beta = 0.5 * alpha_L^{-0.602} * (13.018 - 22.289*s + 17.018*s^2) / (1 - 0.278*s + 1.429*s^2)
```

**gamma:**

```
If s < 0:
    gamma = 1 / (1.1 * sqrt(1 - s - 2*s^2 - 3*s^3)) + 0.75
Else:
    A = -10/3 * alpha + 16/3
    gamma = (1 - 0.5*s) / (1.1 - 2*s + A*s^2 + 4*s^3) + 0.75
```

> **Note:** Both `beta` and `gamma` depend on the **sign** of `s`, not just its magnitude.

---

## 7. Smoothing and Stabilization (`model_cdbm_smooth`)

### 7.1 Parameter

```
model_cdbm_smooth = 0   ! Original CDBM (no modification)
model_cdbm_smooth = 1   ! 5-point smoothing + shear floor + chi limits (default)
```

### 7.2 Five-point boxcar smoothing of magnetic shear

Applied before the main radial loop:

```
S_HM(1)          = S(1)
S_HM(2)          = (S(1) + S(2) + S(3)) / 3
S_HM(NR)         = (S(NR-2) + S(NR-1) + S(NR) + S(NR+1) + S(NR+2)) / 5    for NR = 3, ..., NRMAX-2
S_HM(NRMAX-1)    = (S(NRMAX-2) + S(NRMAX-1) + S(NRMAX)) / 3
S_HM(NRMAX)      = S(NRMAX)
```

When `model_cdbm_smooth = 0`: `S_HM = S` (no smoothing).

### 7.3 Shear floor

For `MDLKAI = 31` (when `model_cdbm_smooth >= 1`):

```
SHEARL = max(s, 0.5)
```

For `MDLKAI = 130:139` (when `model_cdbm_smooth >= 1`):

```
SHEARL = max(S_HM, 0.5)
```

### 7.4 Chi limits (MDLKAI = 31 only)

When `model_cdbm_smooth >= 1`:

```
0.1 <= chi_e <= 50    [m^2/s]
0.1 <= chi_i <= 50    [m^2/s]
```

---

## 8. Numerical Implementation Details

### 8.1 Grid quantities

All transport coefficients are computed on the **half-grid** (`RG`, between mesh points `RM`). Densities and temperatures are averaged:

```
n_e  = 0.5 * (RN(NR+1, 1) + RN(NR, 1))
T_e  = 0.5 * (RT(NR+1, 1) + RT(NR, 1))
```

At the boundary (`NR = NRMAX`), surface values `PNS`, `PTS` are used, and gradients use three-point extrapolation (`DERIV3P`).

### 8.2 Pressure gradient

Total pressure includes electrons, all ion species, and fast particles:

```
P_{NR+1} = sum_{ions} n_s * T_s + n_e * T_e + P_{beam} + P_{add}
P_{NR}   = (same at NR)
dp/dr    = (P_{NR+1} - P_{NR}) / (dr * a)
```

### 8.3 Safety factor derivative

Three-point centered difference for interior points:

```
dq/dr = (q_{NR+1} - q_{NR-1}) / (2 * dr)
```

Forward/backward difference at boundaries:

```
NR=1:     dq/dr = (4*q_2 - 3*q_1 - q_3) / (2*dr)
NR=NRMAX: dq/dr = (3*q_N - 4*q_{N-1} + q_{N-2}) / (2*dr)
```

### 8.4 ExB velocity and shear

```
V_ExB(NR) = -E_r(NR) / B
omega_ExB(NR) = (s - 1) * V_ExB / rho + dV_ExB/dr
```

---

## 9. Input Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `MDLKAI` | 31 | Transport model selector |
| `CK0` | 12.0 | Electron chi multiplier |
| `CK1` | 12.0 | Ion chi multiplier |
| `CWEB` | 1.0 | ExB suppression coefficient |
| `CALF` | 1.0 | Alpha scaling in form factor |
| `MDLCD05` | 0 | 0: original CDBM, nonzero: CDBM05 with elongation |
| `model_cdbm_smooth` | 1 | 0: no smoothing, 1: 5-point smoothing + limits |

---

## 10. Output Variables

### 10.1 Transport coefficients

| Variable | Meaning |
|----------|---------|
| `AKDW(NR, 1)` = `AKE` | Electron anomalous thermal diffusivity chi_e [m^2/s] |
| `AKDW(NR, 2)` = `AKD` | Ion (D) anomalous thermal diffusivity chi_i [m^2/s] |
| `AKDW(NR, 3)` | Ion (T) anomalous thermal diffusivity [m^2/s] |
| `AKDW(NR, 4)` | Ion (He) anomalous thermal diffusivity [m^2/s] |

### 10.2 Diagnostic variables (for CASE 30:39)

| Variable | Content |
|----------|---------|
| `VGR1(NR, 1)` | Form factor F |
| `VGR1(NR, 2)` | Magnetic shear s |
| `VGR1(NR, 3)` | Normalized alpha |
| `VGR2(NR, 1)` | Radial electric field E_r |
| `VGR2(NR, 2)` | ExB drift velocity V_ExB |
| `VGR3(NR, 1)` | ExB suppression factor f_ExB |
| `VGR3(NR, 2)` | |omega_ExB| |

### 10.3 CSV output

Via `TRGR1D` -> `TR_WRITE_CSV`, the following profiles are written to CSV files:

- `AKE, AKNCE, AKDWE [m^2/s] vs r` — electron chi components
- `AKD, AKNCD, AKDWD [m^2/s] vs r` — ion chi components

---

## 11. References

1. M. Yagi et al., "Nonlinear simulation of tearing mode and the resulting anomalous transport", Plasma Phys. Control. Fusion **39** (1997) A421.
2. A. Fukuyama et al., "Transport simulation with current diffusive ballooning mode", Plasma Phys. Control. Fusion **37** (1995) 611.
3. K. Itoh, S.-I. Itoh, A. Fukuyama, "Transport and Structural Formation in Plasmas", IOP Publishing (1999).

---

## Appendix A: Physical Constants

| Constant | Symbol | Value | Unit |
|----------|--------|-------|------|
| Elementary charge | `e` | 1.602176487 × 10^{-19} | C |
| Electron mass | `m_e` | 9.10938215 × 10^{-31} | kg |
| Proton mass | `m_p` | 1.672621637 × 10^{-27} | kg |
| Speed of light | `c` | 2.99792458 × 10^8 | m/s |
| Permeability | `mu_0` | 4π × 10^{-7} | H/m |
| Permittivity | `eps_0` | `1/(c^2 * mu_0)` | F/m |
| keV to Joule | `e_{keV}` | 1.602 × 10^{-16} | J |

## Appendix B: Source File Map

| File | Content |
|------|---------|
| `trcoef.f90:TRCFDW` | Main transport coefficient routine (lines 19–1135) |
| `trcoef.f90:TRCOFS` | Standard form factor F(s, alpha, kappa_q) |
| `trcoef.f90:TRCOFSS` | Simplified form factor |
| `trcoef.f90:TRCOFSX` | Extended form factor with Shafranov shift |
| `trcoef.f90:TRCOFT` | Alternative form factor (MDLKAI=50) |
| `trcoef.f90:FEXB` | ExB shear suppression function |
| `trcomm.f90` | Parameter declarations (CK0, CK1, CWEB, CALF, ...) |
| `trinit.f90` | Default parameter values |
| `trparm.f90` | Namelist definitions for input |
