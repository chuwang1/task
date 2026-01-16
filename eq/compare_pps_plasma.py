#!/usr/bin/env python3
"""
Compare PPS from different sources - plasma region only.
"""
import matplotlib
matplotlib.use('Agg')
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Read data
eqgs1d = pd.read_csv('eqgs1d_01_PPS.csv')
eqgs1d.columns = eqgs1d.columns.str.strip()
eqdata = pd.read_csv('eqdata_1D_profiles.csv')

# Filter eqgs1d to plasma region only (where PPS > 0)
plasma_idx = eqgs1d['PPS_MPa'] > 0.001
eqgs1d_plasma = eqgs1d[plasma_idx]

psipa_plasma = eqgs1d_plasma['PSIP'].max()
psips_max = eqdata['PSIPS'].max()

print(f"eqgs1d plasma boundary PSIP: {psipa_plasma:.4f}")
print(f"eqdata PSIPS max: {psips_max:.4f}")
print(f"Difference: {abs(psipa_plasma - psips_max):.4f} ({abs(psipa_plasma - psips_max)/psips_max*100:.1f}%)")

# Create figure
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Plot 1: Only plasma region with raw axes
ax1 = axes[0]
ax1.plot(eqgs1d_plasma['PSIP'], eqgs1d_plasma['PPS_MPa'], 'b-o', markersize=4, label=f'eqgs1d: PPS (PSIP_max={psipa_plasma:.1f})')
ax1.plot(eqdata['PSIPS'], eqdata['PPPS']*1e-6, 'r-x', markersize=3, label=f'eqdata: PPPS (PSIPS_max={psips_max:.1f})')
ax1.set_xlabel('Poloidal Flux [Wb]')
ax1.set_ylabel('Pressure [MPa]')
ax1.set_title('Plasma Region Only (Raw X-axis)')
ax1.legend()
ax1.grid(True)

# Plot 2: Normalized x-axis (plasma only)
ax2 = axes[1]
psip_norm = eqgs1d_plasma['PSIP'] / psipa_plasma
psips_norm = eqdata['PSIPS'] / psips_max
ax2.plot(psip_norm, eqgs1d_plasma['PPS_MPa'], 'b-o', markersize=4, label='eqgs1d: PPS')
ax2.plot(psips_norm, eqdata['PPPS']*1e-6, 'r-x', markersize=3, label='eqdata: PPPS')
ax2.set_xlabel('Normalized Poloidal Flux (0=axis, 1=plasma edge)')
ax2.set_ylabel('Pressure [MPa]')
ax2.set_title('Normalized X-axis (Plasma Region)')
ax2.legend()
ax2.grid(True)

plt.tight_layout()
plt.savefig('pps_comparison_plasma.png', dpi=150)
print("\nSaved: pps_comparison_plasma.png")
