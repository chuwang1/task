#!/usr/bin/env python3
"""
Compare PPS from different sources to verify if differences are due to x-axis.
"""
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Read eqgs1d PPS (PSIP as x-axis, PPS in MPa)
eqgs1d = pd.read_csv('eqgs1d_01_PPS.csv')
eqgs1d.columns = eqgs1d.columns.str.strip()

# Read eqdata PPS (PSIPS as x-axis, PPPS in Pa)
eqdata = pd.read_csv('eqdata_1D_profiles.csv')

# Get PSIPA for normalization (max value of PSIP)
psipa = eqgs1d['PSIP'].max()
print(f"PSIPA (max PSIP from eqgs1d): {psipa:.4f}")
print(f"Max PSIPS from eqdata: {eqdata['PSIPS'].max():.4f}")

# Create figure with 2 subplots
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Plot 1: Raw data with original x-axes
ax1 = axes[0]
ax1.plot(eqgs1d['PSIP'], eqgs1d['PPS_MPa'], 'b-o', markersize=3, label='eqgs1d: PSIP vs PPS')
ax1.plot(eqdata['PSIPS'], eqdata['PPPS']*1e-6, 'r-x', markersize=3, label='eqdata: PSIPS vs PPPS')
ax1.set_xlabel('PSIP / PSIPS (different scales!)')
ax1.set_ylabel('Pressure [MPa]')
ax1.set_title('Raw Data (Different X-axes)')
ax1.legend()
ax1.grid(True)

# Plot 2: Normalized x-axis comparison
ax2 = axes[1]
# Normalize PSIP to [0,1]
psip_norm = eqgs1d['PSIP'] / psipa
# Normalize PSIPS to [0,1]  
psips_norm = eqdata['PSIPS'] / eqdata['PSIPS'].max()
ax2.plot(psip_norm, eqgs1d['PPS_MPa'], 'b-o', markersize=3, label='eqgs1d: PPS(PSIP_norm)')
ax2.plot(psips_norm, eqdata['PPPS']*1e-6, 'r-x', markersize=3, label='eqdata: PPPS(PSIPS_norm)')
ax2.set_xlabel('Normalized Poloidal Flux (0=axis, 1=edge)')
ax2.set_ylabel('Pressure [MPa]')
ax2.set_title('Normalized X-axis Comparison')
ax2.legend()
ax2.grid(True)

plt.tight_layout()
plt.savefig('pps_comparison_axes.png', dpi=150)
print("Saved: pps_comparison_axes.png")

# Print some statistics
print(f"\neqgs1d PPS: {len(eqgs1d)} points, P(0)={eqgs1d['PPS_MPa'].iloc[0]:.4f} MPa, P(edge)={eqgs1d['PPS_MPa'].iloc[-1]:.4f} MPa")
print(f"eqdata PPPS: {len(eqdata)} points, P(0)={eqdata['PPPS'].iloc[0]*1e-6:.4f} MPa, P(edge)={eqdata['PPPS'].iloc[-1]*1e-6:.4f} MPa")
plt.show()
