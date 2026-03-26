# CDBM Smoothing Improvement (`model_cdbm_smooth`)

## Problem

The CDBM (Current Diffusivity Ballooning Mode) transport model in TR produces a sharp spike in the thermal diffusivity chi at r/a = 0.2-0.3. This spike originates from oscillations in the magnetic shear profile S(r) near the magnetic axis, which are amplified through the CDBM form factor function `TRCOFS(S, alpha, curv)`.

When the magnetic shear crosses zero or becomes negative (as it often does near the axis where dq/dr changes sign), the form factor exhibits a discontinuous jump. This jump propagates into chi and can cause:

- Numerical instability in time-stepping
- Unphysically large or small transport coefficients
- Oscillatory temperature profiles in the core

## Solution

Three coordinated fixes are applied, controlled by a single parameter `model_cdbm_smooth`:

### Fix 1: 5-point boxcar smoothing of magnetic shear

Before the transport coefficient loop, a smoothed shear array `S_HM` is pre-computed:

```
S_HM(1)          = S(1)
S_HM(2)          = (S(1) + S(2) + S(3)) / 3
S_HM(3:NRMAX-2)  = 5-point average of S
S_HM(NRMAX-1)    = (S(NRMAX-2) + S(NRMAX-1) + S(NRMAX)) / 3
S_HM(NRMAX)      = S(NRMAX)
```

This removes the grid-scale oscillation in S without altering its large-scale structure.

**Affected code**: `trcoef.f90`, between the VEXB loop and the main NR loop.

### Fix 2: Minimum shear threshold

A floor value of 0.5 is applied to the shear before it enters the CDBM form factor:

```fortran
SHEARL = MAX(S(NR), 0.5)       ! case(31)
SHEARL = MAX(S_HM(NR), 0.5)    ! case(130:139)
```

This prevents the form factor from diverging when S approaches zero or goes negative.

**Affected code**: `trcoef.f90`, `case(31)` and `case(130:139)`.

### Fix 3: Chi upper and lower limits (case 31 only)

For `case(31)` (MDLKAI=31), the computed chi is clamped to a physically reasonable range:

```fortran
AKDWEL = MIN(AKDWEL, 50.0)   ! max chi_e = 50 m^2/s
AKDWIL = MIN(AKDWIL, 50.0)   ! max chi_i = 50 m^2/s
AKDWEL = MAX(AKDWEL, 0.1)    ! min chi_e = 0.1 m^2/s
AKDWIL = MAX(AKDWIL, 0.1)    ! min chi_i = 0.1 m^2/s
```

**Affected code**: `trcoef.f90`, `case(31)`.

## Usage

The parameter `model_cdbm_smooth` is set in the input namelist `&TR`:

```
model_cdbm_smooth = 1   ! Enable all three fixes (default)
model_cdbm_smooth = 0   ! Disable: use original CDBM behavior
```

### Behavior summary

| Fix | `model_cdbm_smooth=0` | `model_cdbm_smooth=1` |
|-----|----------------------|----------------------|
| S_HM smoothing | S_HM = S (pass-through) | 5-point boxcar average |
| Shear minimum | No floor | `MAX(S or S_HM, 0.5)` |
| Chi limits (case 31) | No clamp | 0.1 <= chi <= 50 m^2/s |

## Affected MDLKAI values

- **MDLKAI=31**: CDBM with s-alpha and curvature (`case(31)`)
- **MDLKAI=130-139**: CDBM with ExB shear and elongation (`case(130:139)`)

Other MDLKAI values are not affected by this parameter.

## Files modified

| File | Change |
|------|--------|
| `trcomm.f90` | Declare `INTEGER:: model_cdbm_smooth` |
| `trinit.f90` | Default `model_cdbm_smooth = 1` |
| `trparm.f90` | Add to `NAMELIST /TR/` |
| `trcoef.f90` | Conditional smoothing, shear floor, and chi clamp |
