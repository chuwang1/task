#!/usr/bin/env python3
"""
Reproduce PNFCL calculation from CSV parameters
Based on trpnf.f90 TRNFDT subroutine
"""
import numpy as np
import matplotlib.pyplot as plt

# Physical constants (from trcomm.f90)
PI = np.pi
AME = 9.1093897e-31     # electron mass [kg]
AMP = 1.6726231e-27     # proton mass [kg]
AMD = 2.0 * AMP         # deuterium mass
AMT = 3.0 * AMP         # tritium mass
AMA = 4.0 * AMP         # alpha (He4) mass
RKEV = 1.6021773e-16    # keV to J
VC = 2.99792458e8       # speed of light [m/s]

# Charge numbers
Z_e = -1
Z_D = 1
Z_T = 1
Z_He4 = 2

# Mass numbers (for COULOG calculation)
PA_He4 = 4.0

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

def HY(V):
    """
    Energy partition function HY(V)
    From trpnb.f90 line 651-661
    """
    return 2.0 * (np.log((V**3 + 1.0) / (V + 1.0)**3) / 6.0 +
                  (np.arctan((2.0*V - 1.0) / np.sqrt(3.0)) + PI/6.0) / np.sqrt(3.0)) / V**2

def COULOG(ne_20, Te_keV):
    """
    Coulomb logarithm (simplified)
    From trcoef.f90
    """
    # ln(Lambda) ≈ 16.1 - 1.15*log10(ne) + 2.30*log10(Te)
    # ne in 10^20/m^3, Te in keV
    ne_cm3 = ne_20 * 1e14  # convert to cm^-3
    return 16.1 - 1.15*np.log10(ne_cm3) + 2.30*np.log10(Te_keV*1000)  # Te in eV

def calculate_pnfcl(r_a, ne, nD, nT, nA, Te, TD, TT, TA, WF):
    """
    Calculate PNFCL based on trpnf.f90 TRNFDT subroutine

    Inputs (all arrays):
        ne, nD, nT, nA: densities [10^20/m^3]
        Te, TD, TT, TA: temperatures [keV]
        WF: fast ion energy density [keV * 10^20/m^3] (needs to be derived from WF[MJ])

    Returns:
        PNFIN, PNFCL_e, PNFCL_D, PNFCL_T, PNFCL_He4 [MW/m^3]
    """
    nr = len(r_a)

    # Alpha particle velocity (3.5 MeV)
    VF = np.sqrt(2.0 * 3.5e3 * RKEV / AMA)  # [m/s]

    # Initialize output arrays
    PNFIN = np.zeros(nr)
    PNFCL_e = np.zeros(nr)
    PNFCL_D = np.zeros(nr)
    PNFCL_T = np.zeros(nr)
    PNFCL_He4 = np.zeros(nr)
    HYF_arr = np.zeros(nr)
    TAUF_arr = np.zeros(nr)

    for i in range(nr):
        ANE = ne[i]  # [10^20/m^3]
        TE = Te[i]   # [keV]

        if ANE < 1e-10 or TE < 0.1 or WF[i] < 1e-10:
            continue

        # P1 parameter (trpnf.f90 line 265)
        # P1 = 3*sqrt(0.5*pi)*me/ne * (Te*RKEV/me)^1.5
        P1 = 3.0 * np.sqrt(0.5 * PI) * AME / (ANE * 1e20) * (np.abs(TE) * RKEV / AME)**1.5

        # Velocity parameters for each species (trpnf.f90 lines 266-269)
        VCD3 = P1 * (nD[i] * 1e20) * Z_D**2 / AMD
        VCT3 = P1 * (nT[i] * 1e20) * Z_T**2 / AMT
        VCA3 = P1 * (nA[i] * 1e20) * Z_He4**2 / AMA
        VC3 = VCD3 + VCT3 + VCA3
        VCR = VC3**(1.0/3.0)

        # Energy partition factor HYF (trpnf.f90 line 271)
        V_ratio = VF / VCR
        HYF = HY(V_ratio)
        HYF_arr[i] = HYF

        # Coulomb logarithm
        ln_lambda = COULOG(ANE, TE)

        # Slowing down time TAUS (trpnf.f90 lines 272-273)
        # TAUS = 0.2 * A_He4 * Te^1.5 / (Z_He4^2 * ne * COULOG)
        TAUS = 0.2 * PA_He4 * np.abs(TE)**1.5 / (Z_He4**2 * ANE * ln_lambda)

        # Fast ion thermalization time (trpnf.f90 line 274)
        TAUF = 0.5 * TAUS * (1.0 - HYF)
        TAUF_arr[i] = TAUF

        # PNFIN calculation (trpnf.f90 line 282)
        # WF is in keV * 10^20/m^3 units, PNFIN should be in W/m^3
        # PNFIN = WF * RKEV * 1e20 / TAUF
        PNFIN[i] = WF[i] * RKEV * 1e20 / TAUF / 1e6  # Convert to MW/m^3

        # PNFCL distribution (trpnf.f90 lines 283-286)
        PNFCL_e[i] = (1.0 - HYF) * PNFIN[i]
        PNFCL_D[i] = (VCD3 / VC3) * HYF * PNFIN[i]
        PNFCL_T[i] = (VCT3 / VC3) * HYF * PNFIN[i]
        PNFCL_He4[i] = (VCA3 / VC3) * HYF * PNFIN[i]

    return PNFIN, PNFCL_e, PNFCL_D, PNFCL_T, PNFCL_He4, HYF_arr, TAUF_arr

# Read CSV data
print("Reading CSV data...")
data_n = read_tr_csv('tr_data_017.csv')   # n(NS) vs r
data_T = read_tr_csv('tr_data_019.csv')   # T(NS) vs r
data_WF = read_tr_csv('tr_data_026.csv')  # WB, WF vs r
data_pnfcl = read_tr_csv('tr_data_028.csv')  # PNFIN, PNFCL vs r (original)

# Extract data
r_a = data_n[:, 0]
ne = data_n[:, 1]   # [10^20/m^3]
nD = data_n[:, 2]
nT = data_n[:, 3]
nA = data_n[:, 4]

Te = data_T[:, 1]   # [keV]
TD = data_T[:, 2]
TT = data_T[:, 3]
TA = data_T[:, 4]

# WF from CSV is in different units - needs interpretation
# In tr_data_026.csv, WF appears to be energy density related
# The code uses RW which is WF in keV * 10^20/m^3
# Let's use the WF column and assume it needs scaling
WF_raw = data_WF[:, 2]  # Raw WF from CSV

# Original PNFCL from CSV
r_a_orig = data_pnfcl[:, 0]
PNFIN_orig = data_pnfcl[:, 1]
PNFCL_e_orig = data_pnfcl[:, 2]
PNFCL_D_orig = data_pnfcl[:, 3]
PNFCL_T_orig = data_pnfcl[:, 4]
PNFCL_He4_orig = data_pnfcl[:, 5]

# Interpolate to same grid
from scipy.interpolate import interp1d
if len(r_a) != len(r_a_orig):
    f_ne = interp1d(r_a, ne, fill_value='extrapolate')
    f_nD = interp1d(r_a, nD, fill_value='extrapolate')
    f_nT = interp1d(r_a, nT, fill_value='extrapolate')
    f_nA = interp1d(r_a, nA, fill_value='extrapolate')
    f_Te = interp1d(r_a, Te, fill_value='extrapolate')
    f_TD = interp1d(r_a, TD, fill_value='extrapolate')
    f_TT = interp1d(r_a, TT, fill_value='extrapolate')
    f_TA = interp1d(r_a, TA, fill_value='extrapolate')
    f_WF = interp1d(r_a, WF_raw, fill_value='extrapolate')

    r_a = r_a_orig
    ne = f_ne(r_a)
    nD = f_nD(r_a)
    nT = f_nT(r_a)
    nA = f_nA(r_a)
    Te = f_Te(r_a)
    TD = f_TD(r_a)
    TT = f_TT(r_a)
    TA = f_TA(r_a)
    WF_raw = f_WF(r_a)

# WF needs proper unit conversion
# From trpnf.f90: RW(NR,NNBMAX+NNF) is the fast ion energy density
# Let's estimate WF from PNFIN and TAUF relationship
# PNFIN = WF * RKEV * 1e20 / TAUF
# We need to find the right scaling factor

print("\nCalculating PNFCL from parameters...")

# First, let's estimate WF that would give the correct PNFIN
# We'll iterate to find the scaling factor
def find_wf_scaling(scale):
    WF = WF_raw * scale
    PNFIN_calc, _, _, _, _, _, _ = calculate_pnfcl(r_a, ne, nD, nT, nA, Te, TD, TT, TA, WF)
    # Compare at center
    idx = 5  # near center
    if PNFIN_orig[idx] > 0 and PNFIN_calc[idx] > 0:
        return np.abs(PNFIN_calc[idx] / PNFIN_orig[idx] - 1.0)
    return 1e10

from scipy.optimize import minimize_scalar
result = minimize_scalar(find_wf_scaling, bounds=(0.001, 1000), method='bounded')
wf_scale = result.x
print(f"WF scaling factor found: {wf_scale:.4f}")

# Calculate with optimal scaling
WF = WF_raw * wf_scale
PNFIN_calc, PNFCL_e_calc, PNFCL_D_calc, PNFCL_T_calc, PNFCL_He4_calc, HYF_arr, TAUF_arr = \
    calculate_pnfcl(r_a, ne, nD, nT, nA, Te, TD, TT, TA, WF)

# Create comparison plots
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Plot 1: PNFIN comparison
ax1 = axes[0, 0]
ax1.plot(r_a, PNFIN_orig, 'k-', linewidth=2, label='CSV (original)')
ax1.plot(r_a, PNFIN_calc, 'r--', linewidth=2, label='Calculated')
ax1.set_xlabel('r/a', fontsize=12)
ax1.set_ylabel('PNFIN [MW/m³]', fontsize=12)
ax1.set_title('PNFIN Comparison', fontsize=14)
ax1.legend()
ax1.grid(True, alpha=0.3)
ax1.set_xlim([0, 1])

# Plot 2: PNFCL_e comparison
ax2 = axes[0, 1]
ax2.plot(r_a, PNFCL_e_orig, 'k-', linewidth=2, label='CSV (original)')
ax2.plot(r_a, PNFCL_e_calc, 'r--', linewidth=2, label='Calculated')
ax2.set_xlabel('r/a', fontsize=12)
ax2.set_ylabel('PNFCL_e [MW/m³]', fontsize=12)
ax2.set_title('PNFCL to Electrons Comparison', fontsize=14)
ax2.legend()
ax2.grid(True, alpha=0.3)
ax2.set_xlim([0, 1])

# Plot 3: PNFCL to ions comparison
ax3 = axes[1, 0]
PNFCL_i_orig = PNFCL_D_orig + PNFCL_T_orig + PNFCL_He4_orig
PNFCL_i_calc = PNFCL_D_calc + PNFCL_T_calc + PNFCL_He4_calc
ax3.plot(r_a, PNFCL_i_orig, 'k-', linewidth=2, label='CSV (original)')
ax3.plot(r_a, PNFCL_i_calc, 'r--', linewidth=2, label='Calculated')
ax3.set_xlabel('r/a', fontsize=12)
ax3.set_ylabel('PNFCL_i [MW/m³]', fontsize=12)
ax3.set_title('PNFCL to Ions Comparison', fontsize=14)
ax3.legend()
ax3.grid(True, alpha=0.3)
ax3.set_xlim([0, 1])

# Plot 4: HYF and TAUF profiles
ax4 = axes[1, 1]
ax4_twin = ax4.twinx()
ax4.plot(r_a, HYF_arr, 'b-', linewidth=2, label='HYF (ion fraction)')
ax4_twin.plot(r_a, TAUF_arr, 'g-', linewidth=2, label='TAUF [s]')
ax4.set_xlabel('r/a', fontsize=12)
ax4.set_ylabel('HYF', fontsize=12, color='blue')
ax4_twin.set_ylabel('TAUF [s]', fontsize=12, color='green')
ax4.set_title('Energy Partition HYF and Thermalization Time', fontsize=14)
ax4.tick_params(axis='y', labelcolor='blue')
ax4_twin.tick_params(axis='y', labelcolor='green')
ax4.grid(True, alpha=0.3)
ax4.set_xlim([0, 1])

plt.tight_layout()
plt.savefig('pnfcl_comparison.png', dpi=150)
plt.show()

# Print comparison statistics
print("\n" + "="*60)
print("PNFCL Calculation Comparison")
print("="*60)

# Calculate relative errors
mask = PNFIN_orig > 0.01
if np.sum(mask) > 0:
    rel_err_pnfin = np.mean(np.abs(PNFIN_calc[mask] - PNFIN_orig[mask]) / PNFIN_orig[mask]) * 100
    rel_err_e = np.mean(np.abs(PNFCL_e_calc[mask] - PNFCL_e_orig[mask]) / PNFCL_e_orig[mask]) * 100
    rel_err_i = np.mean(np.abs(PNFCL_i_calc[mask] - PNFCL_i_orig[mask]) / PNFCL_i_orig[mask]) * 100

    print(f"Mean relative error (r/a < 0.9):")
    print(f"  PNFIN:   {rel_err_pnfin:.1f}%")
    print(f"  PNFCL_e: {rel_err_e:.1f}%")
    print(f"  PNFCL_i: {rel_err_i:.1f}%")

print(f"\nHYF (ion heating fraction) range: {HYF_arr[mask].min():.3f} - {HYF_arr[mask].max():.3f}")
print(f"TAUF (thermalization time) range: {TAUF_arr[mask].min():.3f} - {TAUF_arr[mask].max():.3f} s")
print(f"\nElectron fraction (1-HYF): {(1-HYF_arr[5])*100:.1f}%")
print(f"Ion fraction (HYF): {HYF_arr[5]*100:.1f}%")
print("="*60)
print("\nSaved: pnfcl_comparison.png")
