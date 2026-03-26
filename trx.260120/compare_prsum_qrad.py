#!/usr/bin/env python3
"""
Compare TR PRSUM with OMFIT qrad to identify boundary discrepancies
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Read TR PRSUM data
tr_data = []
with open('tr_data_022.csv', 'r') as f:
    lines = f.readlines()
    for line in lines[2:]:
        if line.strip():
            values = [float(x) for x in line.split(',')]
            tr_data.append(values)
tr_data = np.array(tr_data)
r_a_tr = tr_data[:, 0]
PRSUM_tr = tr_data[:, 1]

# Read OMFIT qrad data
omfit_data = pd.read_csv('omfit_qrad_for_tr.csv')
r_a_omfit = omfit_data['r_a'].values
qrad_omfit = omfit_data['qrad_MW_m3'].values

# Interpolate OMFIT to TR grid
qrad_interp = np.interp(r_a_tr, r_a_omfit, qrad_omfit)

# Compare
print("=" * 70)
print("Comparison: TR PRSUM vs OMFIT qrad (interpolated to TR grid)")
print("=" * 70)
print(f"{'r/a':>8} {'PRSUM(TR)':>12} {'qrad(OMFIT)':>12} {'Diff':>10} {'Ratio':>8}")
print("-" * 70)

for i in range(len(r_a_tr)):
    diff = PRSUM_tr[i] - qrad_interp[i]
    ratio = PRSUM_tr[i] / qrad_interp[i] if qrad_interp[i] > 0.001 else 0

    # Highlight large differences
    flag = " ***" if abs(diff) > 0.02 else ""
    print(f"{r_a_tr[i]:>8.4f} {PRSUM_tr[i]:>12.6f} {qrad_interp[i]:>12.6f} {diff:>10.6f} {ratio:>8.2f}{flag}")

# Plot comparison
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Plot 1: Full profile comparison
ax1 = axes[0, 0]
ax1.plot(r_a_tr, PRSUM_tr, 'b-o', markersize=4, label='TR PRSUM')
ax1.plot(r_a_omfit, qrad_omfit, 'r-', linewidth=2, label='OMFIT qrad')
ax1.set_xlabel('r/a')
ax1.set_ylabel('Radiation [MW/m³]')
ax1.set_title('Full Profile Comparison')
ax1.legend()
ax1.grid(True, alpha=0.3)
ax1.set_xlim([0, 1])

# Plot 2: Edge region (r/a > 0.85)
ax2 = axes[0, 1]
mask_tr = r_a_tr > 0.85
mask_omfit = r_a_omfit > 0.85
ax2.plot(r_a_tr[mask_tr], PRSUM_tr[mask_tr], 'b-o', markersize=6, label='TR PRSUM')
ax2.plot(r_a_omfit[mask_omfit], qrad_omfit[mask_omfit], 'r-', linewidth=2, label='OMFIT qrad')
ax2.set_xlabel('r/a')
ax2.set_ylabel('Radiation [MW/m³]')
ax2.set_title('Edge Region (r/a > 0.85)')
ax2.legend()
ax2.grid(True, alpha=0.3)

# Plot 3: Difference
ax3 = axes[1, 0]
ax3.plot(r_a_tr, PRSUM_tr - qrad_interp, 'g-o', markersize=4)
ax3.axhline(y=0, color='k', linestyle='--')
ax3.set_xlabel('r/a')
ax3.set_ylabel('PRSUM - qrad [MW/m³]')
ax3.set_title('Difference (TR - OMFIT)')
ax3.grid(True, alpha=0.3)
ax3.set_xlim([0, 1])

# Plot 4: Ratio
ax4 = axes[1, 1]
ratio = np.where(qrad_interp > 0.01, PRSUM_tr / qrad_interp, np.nan)
ax4.plot(r_a_tr, ratio, 'm-o', markersize=4)
ax4.axhline(y=1, color='k', linestyle='--')
ax4.set_xlabel('r/a')
ax4.set_ylabel('PRSUM / qrad')
ax4.set_title('Ratio (TR / OMFIT)')
ax4.grid(True, alpha=0.3)
ax4.set_xlim([0, 1])
ax4.set_ylim([0, 6])

plt.tight_layout()
plt.savefig('prsum_qrad_comparison.png', dpi=150)
print(f"\nSaved: prsum_qrad_comparison.png")

# Check the last few points in OMFIT CSV
print("\n" + "=" * 70)
print("OMFIT qrad values near edge (last 10 points):")
print("=" * 70)
for i in range(-10, 0):
    print(f"r/a = {r_a_omfit[i]:.6f}, qrad = {qrad_omfit[i]:.6f} MW/m³")
