#!/usr/bin/env python3
"""
Compare AR1RHO from TASK/TR with EXPRO_ave_grad_r from OMFIT.

TASK/TR:  AR1RHO = <|grad rho|>    (rho = rho_tor, dimensionless)
OMFIT:    EXPRO_ave_grad_r = <|grad r|>  (r = rmin)

Since rho = f(r) is a flux function (constant on each flux surface),
drho/dr is constant on each surface, so:
  <|grad rho|> = <|grad r|> * |drho/dr|
  => AR1RHO = EXPRO_ave_grad_r * |drho/dr_min|

where drho/dr_min is computed numerically from the rho(rmin) mapping.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# ============================================================
# 1. Read TASK/TR AR1RHO
# ============================================================
tr_lines = []
with open('/Users/dengxiaoya/TASK/CFEDR/git/task/trx.260120/tr_data_134.csv', 'r') as f:
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
tr_ar1rho = tr_data[:, 1]  # <|grad rho|>

print("TR AR1RHO: {} points, range [{:.6f}, {:.6f}]".format(
    len(tr_rho), tr_ar1rho.min(), tr_ar1rho.max()))

# ============================================================
# 2. Read OMFIT data
# ============================================================
extra = pd.read_csv('/Users/dengxiaoya/CFEDRSW/OMFIT_out/outputs_new/input_profiles_extra.csv')
allp = pd.read_csv('/Users/dengxiaoya/CFEDRSW/OMFIT_out/outputs_new/all_profiles.csv')

omfit_rho = allp['rho'].values
omfit_rmin = allp['rmin'].values
expro_ave_grad_r = extra['EXPRO_ave_grad_r'].values
ARHO_EXP = 3.6146142

# ============================================================
# 3. Compute drho/dr_min numerically
# ============================================================
drho_drmin = np.gradient(omfit_rho, omfit_rmin)  # [1/m]

# ============================================================
# 4. Convert: <|grad rho|> = <|grad r|> * |drho/dr|
# ============================================================
omfit_ar1rho = expro_ave_grad_r * np.abs(drho_drmin)

print("OMFIT converted: range [{:.6f}, {:.6f}]".format(
    omfit_ar1rho[1:].min(), omfit_ar1rho[1:].max()))

# ============================================================
# 5. Also try using EXPRO_drdrho with ARHO_EXP correction
#    drho/dr_min = 1 / (dr_min/drho) = 1 / (EXPRO_drdrho * ARHO_EXP)
# ============================================================
drho_drmin_expro = 1.0 / (extra['EXPRO_drdrho'].values * ARHO_EXP)
omfit_ar1rho_expro = expro_ave_grad_r * np.abs(drho_drmin_expro)

# ============================================================
# 6. Interpolate and compare
# ============================================================
tr_interp = interp1d(tr_rho, tr_ar1rho, kind='cubic', fill_value='extrapolate')
tr_on_omfit = tr_interp(omfit_rho)

mask = (omfit_rho > 0.05) & (omfit_rho < 0.95)
ratio_num = omfit_ar1rho[mask] / tr_on_omfit[mask]
ratio_expro = omfit_ar1rho_expro[mask] / tr_on_omfit[mask]

print("\nMethod A (numerical drho/drmin):          mean ratio={:.6f}, std={:.6f}".format(
    np.mean(ratio_num), np.std(ratio_num)))
print("Method B (1/(EXPRO_drdrho*ARHO_EXP)):     mean ratio={:.6f}, std={:.6f}".format(
    np.mean(ratio_expro), np.std(ratio_expro)))

# ============================================================
# 7. Plot
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Panel 1: Direct comparison
ax = axes[0, 0]
ax.plot(tr_rho, tr_ar1rho, 'b-o', markersize=3, linewidth=1.5,
        label=r'TR: AR1RHO = $\langle|\nabla\rho|\rangle$')
ax.plot(omfit_rho, omfit_ar1rho, 'r-', linewidth=2,
        label=r'OMFIT: $\langle|\nabla r|\rangle \times |d\rho/dr|$ (numerical)')
ax.plot(omfit_rho, omfit_ar1rho_expro, 'g--', linewidth=2,
        label=r'OMFIT: via EXPRO\_drdrho $\times$ ARHO\_EXP')
ax.set_xlabel(r'$\rho_{tor}$')
ax.set_ylabel(r'$\langle|\nabla\rho|\rangle$')
ax.set_title(r'$\langle|\nabla\rho|\rangle$ comparison')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# Panel 2: Raw EXPRO_ave_grad_r vs AR1RHO
ax = axes[0, 1]
ax.plot(tr_rho, tr_ar1rho, 'b-o', markersize=3, linewidth=1.5,
        label=r'TR: $\langle|\nabla\rho|\rangle$')
ax.plot(omfit_rho, expro_ave_grad_r, 'r-', linewidth=2,
        label=r'OMFIT: $\langle|\nabla r_{min}|\rangle$ (raw)')
ax.set_xlabel(r'$\rho_{tor}$')
ax.set_ylabel('Value')
ax.set_title('Raw (unconverted) — different coordinates')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# Panel 3: Relative difference
ax = axes[1, 0]
mask_plot = omfit_rho > 0.02
rel_diff_num = (omfit_ar1rho[mask_plot] - tr_on_omfit[mask_plot]) / tr_on_omfit[mask_plot] * 100
rel_diff_exp = (omfit_ar1rho_expro[mask_plot] - tr_on_omfit[mask_plot]) / tr_on_omfit[mask_plot] * 100

ax.plot(omfit_rho[mask_plot], rel_diff_num, 'r-', linewidth=2,
        label='Method A: numerical drho/dr')
ax.plot(omfit_rho[mask_plot], rel_diff_exp, 'g--', linewidth=2,
        label='Method B: via EXPRO_drdrho')
ax.axhline(y=0, color='k', linestyle='-', linewidth=0.5)
ax.set_xlabel(r'$\rho_{tor}$')
ax.set_ylabel('Relative difference [%]')
ax.set_title('(OMFIT - TR) / TR * 100%')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# Panel 4: Jacobian comparison
ax = axes[1, 1]
ax.plot(omfit_rho[1:], np.abs(drho_drmin[1:]), 'b-', linewidth=2,
        label=r'$|d\rho/dr_{min}|$ (numerical)')
ax.plot(omfit_rho[1:], np.abs(drho_drmin_expro[1:]), 'r--', linewidth=2,
        label=r'$1/(EXPRO\_drdrho \times ARHO\_EXP)$')
ax.set_xlabel(r'$\rho_{tor}$')
ax.set_ylabel(r'$|d\rho/dr_{min}|$ [1/m]')
ax.set_title('Jacobian: two methods')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/Users/dengxiaoya/TASK/CFEDR/git/task/trx.260120/compare_ar1rho.png', dpi=150)
plt.close()
print("\nPlot saved to: compare_ar1rho.png")

# ============================================================
# 8. Comparison table
# ============================================================
omfit_interp = interp1d(omfit_rho, omfit_ar1rho, kind='cubic', fill_value='extrapolate')
omfit_interp_exp = interp1d(omfit_rho, omfit_ar1rho_expro, kind='cubic', fill_value='extrapolate')

print("\n" + "="*95)
print("{:>6s}  {:>12s}  {:>12s}  {:>12s}  {:>9s}  {:>9s}".format(
    'rho', 'TR AR1RHO', 'OMFIT(num)', 'OMFIT(expro)', 'diff_n%', 'diff_e%'))
print("="*95)
for r in [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95]:
    tr_val = float(tr_interp(r))
    om_n = float(omfit_interp(r))
    om_e = float(omfit_interp_exp(r))
    dn = (om_n - tr_val) / tr_val * 100
    de = (om_e - tr_val) / tr_val * 100
    print("{:6.2f}  {:12.6f}  {:12.6f}  {:12.6f}  {:8.2f}%  {:8.2f}%".format(
        r, tr_val, om_n, om_e, dn, de))
print("="*95)
print("""
Conversion formula:
  <|grad rho|> = <|grad r_min|> * |drho/dr_min|
  
  Method A: drho/dr_min from np.gradient(rho, rmin) — numerical
  Method B: drho/dr_min = 1/(EXPRO_drdrho * ARHO_EXP)
  Both methods are equivalent (ARHO_EXP = {:.4f}).
""".format(ARHO_EXP))
