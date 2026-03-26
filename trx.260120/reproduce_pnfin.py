#!/usr/bin/env python3
"""
Reproduce PNFIN calculation from density and temperature
PNFIN is the fusion power input (alpha particle heating power)

Key formula:
  PNF = <σv> × nD × nT × Eα × 10^40 [W/m³]

where:
  <σv> = D-T fusion reactivity [m³/s]
  nD, nT = D, T densities [10^20/m³]
  Eα = 3.5 MeV (alpha particle energy)

In steady state: PNFIN ≈ PNF
"""
import numpy as np
import matplotlib.pyplot as plt

# Physical constants
RKEV = 1.6021773e-16    # keV to J

def read_tr_csv(filename):
    """Read TR output CSV file"""
    data = []
    with open(filename, 'r') as f:
        lines = f.readlines()
        for line in lines[2:]:
            if line.strip():
                values = [float(x) for x in line.split(',')]
                data.append(values)
    return np.array(data)

def sigmav_dt(Ti_keV):
    """
    D-T fusion reactivity <σv> [m³/s]
    Based on NRL Plasma Formulary / Bosch-Hale parametrization

    Input: Ti in keV
    Output: <σv> in m³/s
    """
    # Simplified Bosch-Hale formula for D-T
    # Valid for Ti = 0.2 - 100 keV
    Ti = np.clip(Ti_keV, 0.2, 100.0)

    # Coefficients for D-T reaction
    BG = 34.3827  # Gamow constant
    mrc2 = 1124656  # keV (reduced mass × c²)

    C1 = 1.17302e-9
    C2 = 1.51361e-2
    C3 = 7.51886e-2
    C4 = 4.60643e-3
    C5 = 1.35000e-2
    C6 = -1.06750e-4
    C7 = 1.36600e-5

    theta = Ti / (1.0 - (Ti*(C2 + Ti*(C4 + Ti*C6))) / (1.0 + Ti*(C3 + Ti*(C5 + Ti*C7))))
    xi = (BG**2 / (4.0 * theta))**(1.0/3.0)

    sigmav = C1 * theta * np.sqrt(xi / (mrc2 * Ti**3)) * np.exp(-3.0 * xi)

    return sigmav * 1e-6  # Convert from cm³/s to m³/s

def sigmav_dt_simple(Ti_keV):
    """
    Simplified D-T fusion reactivity from TASK code (SIGMAM function)
    Based on trpnf.f90
    """
    Ti = np.abs(Ti_keV)

    # From libsigma or similar - using approximation
    # σv ≈ 3.7e-18 × Ti^(-2/3) × exp(-20/Ti^(1/3)) / H(Ti)
    # where H = Ti/37 + 5.45/(3 + Ti*(1 + (Ti/37.5)^2.8))

    H = Ti/37.0 + 5.45/(3.0 + Ti*(1.0 + (Ti/37.5)**2.8))
    arg = -20.0/Ti**(1.0/3.0)

    # Avoid overflow
    sigmav = np.where(arg >= -100,
                      3.7e-18 * Ti**(-2.0/3.0) * np.exp(arg) / H,
                      0.0)

    return sigmav * 1e-6  # Convert from cm³/s to m³/s

# Read CSV data
print("Reading CSV data...")
data_n = read_tr_csv('tr_data_017.csv')   # n(NS) vs r
data_T = read_tr_csv('tr_data_019.csv')   # T(NS) vs r
data_pnfcl = read_tr_csv('tr_data_028.csv')  # PNFIN, PNFCL vs r (original)

# Extract data
r_a_n = data_n[:, 0]
ne = data_n[:, 1]   # [10^20/m^3]
nD = data_n[:, 2]
nT = data_n[:, 3]
nA = data_n[:, 4]

r_a_T = data_T[:, 0]
Te = data_T[:, 1]   # [keV]
TD = data_T[:, 2]
TT = data_T[:, 3]
TA = data_T[:, 4]

# Original PNFIN from CSV
r_a_orig = data_pnfcl[:, 0]
PNFIN_orig = data_pnfcl[:, 1]

# Interpolate to same grid
from scipy.interpolate import interp1d
f_nD = interp1d(r_a_n, nD, fill_value='extrapolate')
f_nT = interp1d(r_a_n, nT, fill_value='extrapolate')
f_TD = interp1d(r_a_T, TD, fill_value='extrapolate')
f_TT = interp1d(r_a_T, TT, fill_value='extrapolate')

r_a = r_a_orig
nD_interp = f_nD(r_a)
nT_interp = f_nT(r_a)
TD_interp = f_TD(r_a)
TT_interp = f_TT(r_a)

# Calculate effective ion temperature for DT reaction
# Ti_eff = (3*TD + 2*TT) / 5  (from trpnf.f90 SIGMAM function)
Ti_eff = (3.0 * TD_interp + 2.0 * TT_interp) / 5.0

print(f"\nTemperature range: {Ti_eff.min():.1f} - {Ti_eff.max():.1f} keV")
print(f"Density nD range: {nD_interp.min():.3f} - {nD_interp.max():.3f} [10^20/m³]")
print(f"Density nT range: {nT_interp.min():.3f} - {nT_interp.max():.3f} [10^20/m³]")

# Calculate fusion reactivity
sigmav = sigmav_dt(Ti_eff)  # [m³/s]
sigmav_simple = sigmav_dt_simple(Ti_eff)

print(f"\n<σv> range: {sigmav.min():.2e} - {sigmav.max():.2e} m³/s")

# Calculate fusion power density
# PNF = σv × nD × nT × Eα × 10^40
# Units: [m³/s] × [10^20/m³] × [10^20/m³] × [J] × 10^40 = W/m³
Ealpha = 3.5e3 * RKEV  # 3.5 MeV in Joules

# nD and nT are in 10^20/m³, so nD*nT*1e40 gives /m^6
PNF_calc = sigmav * nD_interp * nT_interp * Ealpha * 1e40 / 1e6  # MW/m³
PNF_simple = sigmav_simple * nD_interp * nT_interp * Ealpha * 1e40 / 1e6  # MW/m³

# In steady state, PNFIN ≈ PNF (alpha particle power released to plasma)
PNFIN_calc = PNF_calc

print(f"\nCalculated PNFIN range: {PNFIN_calc.min():.3f} - {PNFIN_calc.max():.3f} MW/m³")
print(f"Original PNFIN range:   {PNFIN_orig.min():.3f} - {PNFIN_orig.max():.3f} MW/m³")

# Create plots
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Plot 1: PNFIN comparison
ax1 = axes[0, 0]
ax1.plot(r_a, PNFIN_orig, 'k-', linewidth=2, label='CSV (Fortran)')
ax1.plot(r_a, PNFIN_calc, 'r--', linewidth=2, label='Calculated (Bosch-Hale)')
ax1.plot(r_a, PNF_simple, 'b:', linewidth=2, label='Calculated (Simple)')
ax1.set_xlabel('r/a', fontsize=12)
ax1.set_ylabel('PNFIN [MW/m³]', fontsize=12)
ax1.set_title('PNFIN: CSV vs Calculated from nD, nT, Ti', fontsize=14)
ax1.legend()
ax1.grid(True, alpha=0.3)
ax1.set_xlim([0, 1])

# Plot 2: Ratio
ax2 = axes[0, 1]
mask = PNFIN_orig > 0.01
ratio = np.zeros_like(PNFIN_orig)
ratio[mask] = PNFIN_calc[mask] / PNFIN_orig[mask]
ax2.plot(r_a[mask], ratio[mask], 'b-', linewidth=2)
ax2.axhline(y=1.0, color='k', linestyle='--', linewidth=1)
ax2.set_xlabel('r/a', fontsize=12)
ax2.set_ylabel('Ratio (Calc/CSV)', fontsize=12)
ax2.set_title('Ratio of Calculated to CSV PNFIN', fontsize=14)
ax2.grid(True, alpha=0.3)
ax2.set_xlim([0, 1])
ax2.set_ylim([0.5, 1.5])

# Plot 3: Input profiles (nD, nT)
ax3 = axes[1, 0]
ax3.plot(r_a, nD_interp, 'b-', linewidth=2, label='nD')
ax3.plot(r_a, nT_interp, 'r-', linewidth=2, label='nT')
ax3.set_xlabel('r/a', fontsize=12)
ax3.set_ylabel('Density [10²⁰/m³]', fontsize=12)
ax3.set_title('Input: Density Profiles', fontsize=14)
ax3.legend()
ax3.grid(True, alpha=0.3)
ax3.set_xlim([0, 1])

# Plot 4: Temperature and σv
ax4 = axes[1, 1]
ax4_twin = ax4.twinx()
ax4.plot(r_a, Ti_eff, 'r-', linewidth=2, label='Ti_eff')
ax4_twin.semilogy(r_a, sigmav, 'b-', linewidth=2, label='<σv>')
ax4.set_xlabel('r/a', fontsize=12)
ax4.set_ylabel('Ti [keV]', fontsize=12, color='red')
ax4_twin.set_ylabel('<σv> [m³/s]', fontsize=12, color='blue')
ax4.set_title('Input: Temperature and Fusion Reactivity', fontsize=14)
ax4.tick_params(axis='y', labelcolor='red')
ax4_twin.tick_params(axis='y', labelcolor='blue')
ax4.grid(True, alpha=0.3)
ax4.set_xlim([0, 1])

plt.tight_layout()
plt.savefig('pnfin_reproduction.png', dpi=150)
plt.show()

# Print summary
print("\n" + "="*60)
print("PNFIN Calculation Summary")
print("="*60)
print("\nPNFIN is controlled by:")
print("  1. D-T fusion reactivity <σv>(Ti)")
print("  2. Deuterium density nD")
print("  3. Tritium density nT")
print("\nFormula: PNFIN ≈ PNF = <σv> × nD × nT × Eα")
print(f"\nwhere Eα = 3.5 MeV = {Ealpha:.3e} J")

# Calculate mean relative error
if np.sum(mask) > 0:
    rel_err = np.mean(np.abs(PNFIN_calc[mask] - PNFIN_orig[mask]) / PNFIN_orig[mask]) * 100
    print(f"\nMean relative error: {rel_err:.1f}%")

    # Find scaling factor for better match
    scale = np.mean(PNFIN_orig[mask] / PNFIN_calc[mask])
    print(f"Best scaling factor: {scale:.3f}")

print("\n" + "="*60)
print("\nKey finding: PNFIN ∝ <σv>(Ti) × nD × nT")
print("  - Higher Ti → Higher <σv> → Higher PNFIN")
print("  - Higher nD, nT → Higher PNFIN")
print("  - PNFIN peaks at plasma center where Ti and n are highest")
print("="*60)
print("\nSaved: pnfin_reproduction.png")
