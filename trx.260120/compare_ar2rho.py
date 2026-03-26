#!/usr/bin/env python3
"""
Compare AR2RHO from TASK/TR with EXPRO_ave_grad_r from OMFIT.

TASK/TR:  AR2RHO = <|grad rho|^2>   (rho = rho_tor, dimensionless)
OMFIT:    EXPRO_ave_grad_r = <|grad r|>  (r = rmin, in meters)

Relationship:
  grad(rho) = grad(r) * (drho/dr) 
  <|grad rho|^2> = <|grad r|^2> * (drho/dr)^2

But EXPRO_ave_grad_r is <|grad r|>, not <|grad r|^2>.
<|grad r|>^2 != <|grad r|^2> in general (Jensen's inequality).

Also note: GACODE's 'r' in grad_r may be r/ARHO_EXP (normalized), not physical rmin.
We need to check carefully.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# ============================================================
# 1. Read TASK/TR AR2RHO data
# ============================================================
tr_lines = []
with open('/Users/dengxiaoya/TASK/CFEDR/git/task/trx.260120/tr_data_129.csv', 'r') as f:
    for line in f:
        line = line.strip()
        if line.startswith('Title') or line.startswith('X,') or line == '':
            continue
        parts = line.split(',')
        if len(parts) == 2:
            try:
                tr_lines.append((float(parts[0]), float(parts[1])))
            except ValueError:
                continue

tr_data = np.array(tr_lines)
tr_rho = tr_data[:, 0]
tr_ar2rho = tr_data[:, 1]  # <|grad rho|^2>

print("TR AR2RHO: {} points, range [{:.6f}, {:.6f}]".format(
    len(tr_rho), tr_ar2rho.min(), tr_ar2rho.max()))

# ============================================================
# 2. Read OMFIT data
# ============================================================
extra = pd.read_csv('/Users/dengxiaoya/CFEDRSW/OMFIT_out/outputs_new/input_profiles_extra.csv')
allp = pd.read_csv('/Users/dengxiaoya/CFEDRSW/OMFIT_out/outputs_new/all_profiles.csv')

omfit_rho = allp['rho'].values
omfit_rmin = allp['rmin'].values
expro_ave_grad_r = extra['EXPRO_ave_grad_r'].values   # <|grad r|>
expro_grad_r0 = extra['EXPRO_grad_r0'].values         # |grad r| at theta=0 (outboard midplane)

ARHO_EXP = 3.6146142

# ============================================================
# 3. Compute drho/dr_min numerically
# ============================================================
drho_drmin = np.gradient(omfit_rho, omfit_rmin)

# ============================================================
# 4. Try various conversions
# ============================================================

# Method A: <|grad r|>^2 * (drho/dr)^2
#   This is (<|grad r|>)^2 * (drho/dr)^2, but note <|grad r|>^2 != <|grad r|^2>
omfit_A = expro_ave_grad_r**2 * drho_drmin**2

# Method B: Same as A but using ARHO_EXP for normalization
#   If EXPRO_ave_grad_r is actually <|grad(r/ARHO_EXP)|> (dimensionless grad of normalized r)
#   then: <|grad rho|^2> ~ <|grad(r/ARHO_EXP)|^2> * (drho/d(r/ARHO_EXP))^2
#   = <|grad(r/ARHO_EXP)|>^2 * (ARHO_EXP * drho/dr)^2
drho_dr_norm = ARHO_EXP * drho_drmin  # drho/d(r/ARHO_EXP)
omfit_B = expro_ave_grad_r**2 * drho_dr_norm**2

# Method C: Maybe EXPRO_ave_grad_r is already <|grad(r/a)|> where a=rmin[-1]?
a_geo = omfit_rmin[-1]
drho_dr_a = a_geo * drho_drmin  # drho/d(r/a)  
omfit_C = expro_ave_grad_r**2 * drho_dr_a**2

# Method D: Direct - just look at the ratio
tr_interp = interp1d(tr_rho, tr_ar2rho, kind='cubic', fill_value='extrapolate')
tr_on_omfit = tr_interp(omfit_rho)

# What factor makes EXPRO_ave_grad_r^2 match AR2RHO?
mask = (omfit_rho > 0.1) & (omfit_rho < 0.9)
raw_ratio = tr_on_omfit[mask] / (expro_ave_grad_r[mask]**2)
print("\nTR_AR2RHO / EXPRO_ave_grad_r^2 : mean={:.6f}, std={:.6f}".format(
    np.mean(raw_ratio), np.std(raw_ratio)))
print("  (if constant, this would be (drho/dr)^2 or similar)")

ratio_A = tr_on_omfit[mask] / omfit_A[mask]
print("Method A (ave_grad_r^2 * (drho/drmin)^2):  ratio={:.6f}, std={:.6f}".format(
    np.mean(ratio_A), np.std(ratio_A)))

ratio_B = tr_on_omfit[mask] / omfit_B[mask]
print("Method B (ave_grad_r^2 * (ARHO*drho/drmin)^2): ratio={:.6f}, std={:.6f}".format(
    np.mean(ratio_B), np.std(ratio_B)))

ratio_C = tr_on_omfit[mask] / omfit_C[mask]
print("Method C (ave_grad_r^2 * (a*drho/drmin)^2):  ratio={:.6f}, std={:.6f}".format(
    np.mean(ratio_C), np.std(ratio_C)))

# Method E: the grad in GACODE is grad of (r/a) not physical r
# EXPRO_ave_grad_r = <|grad(rmin/a)|> where a is the normalization
# Then <|grad rho|^2> ~ <|grad(rmin/a)|>^2 * (drho/d(rmin/a))^2
# drho/d(rmin/a) = a * drho/drmin
# Let's try with different 'a' values
for a_try in [a_geo, ARHO_EXP, 1.0]:
    drho_da = a_try * drho_drmin
    test = expro_ave_grad_r**2 * drho_da**2
    r = tr_on_omfit[mask] / test[mask]
    print("  a={:.4f}: ratio={:.6f}, std={:.6f}".format(a_try, np.mean(r), np.std(r)))

# ============================================================
# 5. Plot
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Panel 1: Raw data
ax = axes[0, 0]
ax.plot(tr_rho, tr_ar2rho, 'b-o', markersize=3, label=r'TR: AR2RHO = $\langle|\nabla\rho|^2\rangle$')
ax.plot(omfit_rho, expro_ave_grad_r**2, 'r-', linewidth=2, 
        label=r'OMFIT: EXPRO\_ave\_grad\_r$^2$')
ax.set_xlabel(r'$\rho_{tor}$')
ax.set_ylabel('Value')
ax.set_title('Raw comparison (different units)')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# Panel 2: Best conversion
# Use the formula: <|grad rho|^2> = <|grad r|>^2 * (drho/dr)^2
# where r = rmin (physical), and we accept <|grad r|>^2 ~ <|grad r|^2>
ax = axes[0, 1]
ax.plot(tr_rho, tr_ar2rho, 'b-o', markersize=3, 
        label=r'TR: $\langle|\nabla\rho|^2\rangle$')
ax.plot(omfit_rho, omfit_A, 'r-', linewidth=2,
        label=r'$\langle|\nabla r|\rangle^2 \times (d\rho/dr)^2$')
ax.set_xlabel(r'$\rho_{tor}$')
ax.set_ylabel('Value')
ax.set_title(r'With conversion: $r = r_{min}$ (physical)')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# Panel 3: ratio analysis  
ax = axes[1, 0]
mask_plot = omfit_rho > 0.02
ax.plot(omfit_rho[mask_plot], tr_on_omfit[mask_plot] / omfit_A[mask_plot], 'r-', linewidth=2,
        label='Method A (physical r)')
ax.plot(omfit_rho[mask_plot], tr_on_omfit[mask_plot] / omfit_B[mask_plot], 'g--', linewidth=2,
        label='Method B (r * ARHO_EXP)')
ax.axhline(y=1.0, color='k', linestyle='-', linewidth=0.5)
ax.set_xlabel(r'$\rho_{tor}$')
ax.set_ylabel('TR / OMFIT ratio')
ax.set_title('Ratio analysis')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
ax.set_ylim(0, 2)

# Panel 4: relative difference for best method
ax = axes[1, 1]
# Find which method gives ratio closest to 1
best = omfit_A  # will update below
best_label = 'A'
for meth, label, data in [('A', 'phys r', omfit_A), ('B', 'ARHO r', omfit_B), ('C', 'a_geo r', omfit_C)]:
    r = np.mean(tr_on_omfit[mask] / data[mask])
    if abs(r - 1.0) < abs(np.mean(tr_on_omfit[mask] / best[mask]) - 1.0):
        best = data
        best_label = label

rel_diff = (best[mask_plot] - tr_on_omfit[mask_plot]) / tr_on_omfit[mask_plot] * 100
ax.plot(omfit_rho[mask_plot], rel_diff, 'r-', linewidth=2,
        label='Best method: {}'.format(best_label))
ax.axhline(y=0, color='k', linestyle='-', linewidth=0.5)
ax.set_xlabel(r'$\rho_{tor}$')
ax.set_ylabel('Relative difference [%]')
ax.set_title('(OMFIT - TR) / TR * 100%')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/Users/dengxiaoya/TASK/CFEDR/git/task/trx.260120/compare_ar2rho.png', dpi=150)
plt.close()
print("\nPlot saved to: compare_ar2rho.png")

# ============================================================
# 6. Comparison table
# ============================================================
omfit_A_interp = interp1d(omfit_rho, omfit_A, kind='cubic', fill_value='extrapolate')
print("\n" + "="*80)
print("{:>6s}  {:>12s}  {:>12s}  {:>12s}  {:>9s}".format(
    'rho', 'TR AR2RHO', 'OMFIT conv', 'ave_grad_r^2', 'diff%'))
print("="*80)
for r in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95]:
    tr_val = float(tr_interp(r))
    om_val = float(omfit_A_interp(r))
    agr2 = float(interp1d(omfit_rho, expro_ave_grad_r**2, kind='cubic')(r))
    diff = (om_val - tr_val) / tr_val * 100
    print("{:6.2f}  {:12.6f}  {:12.6f}  {:12.6f}  {:8.2f}%".format(
        r, tr_val, om_val, agr2, diff))
print("="*80)
