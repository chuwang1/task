#!/usr/bin/env python3
"""
Analyze missing terms in temperature reproduction.

From trexec.f90, the temperature equation source term is:
    D(NEQ,NR) = (PIN(NR,NSSN)/(RKEV*1.D20))*DV53
             + (-VI(NEQ,NEQ-1,2,NSW) + C83*DI(NEQ,NEQ-1,2,NSW))*RNV
             + (-VI(NEQ,NEQ,2,NSW) + C83*DI(NEQ,NEQ,2,NSW))*RPV

Where:
- PIN: net power input (heating - radiation - exchange) [W/m³]
- VI, DI: convection velocity and diffusion coefficient
- RNV: density n [10^20/m³]
- RPV: pressure n*T [10^20/m³ * keV]
- C83 = 8/3
- DV53 = DVRHO^(5/3), DV23 = DVRHO^(2/3)

The VI/DI terms represent:
- Heat convection (pinch): AVK (heat pinch) + AV*CC (particle pinch × 1.5)
- Heat diffusion: AK (thermal diffusivity, chi)
- Cross-diffusion: (AD*CC - AK)*RTV (particle diffusion effect on heat)

This script estimates the contribution of these missing terms.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
import os

# CFEDR parameters
R0 = 8.03        # Major radius [m]
A_MINOR = 2.68   # Minor radius [m]
KAPPA = 1.89     # Elongation
KAPPA_S = np.sqrt(KAPPA)

TR_DIR = '/Users/dengxiaoya/TASK/latest/task/trx'


def find_csv_by_title(title_keyword, search_dir=TR_DIR):
    """Find TR CSV file by title keyword."""
    for f in sorted(glob.glob(os.path.join(search_dir, 'tr_data_*.csv'))):
        with open(f, 'r') as fh:
            title = fh.readline().strip()
        if title_keyword in title:
            return f
    return None


def read_tr_csv(filename):
    """Read TR output CSV file."""
    with open(filename, 'r') as f:
        title = f.readline().strip()
    df = pd.read_csv(filename, skiprows=1, index_col=False)
    df = df.loc[:, ~df.columns.str.contains('^Unnamed')]
    df.columns = df.columns.str.strip()
    x = df.iloc[:, 0].values
    data = {col.strip(): df[col].values for col in df.columns[1:]}
    return x, data, title


def main():
    print("="*70)
    print("ANALYSIS OF MISSING TERMS IN TEMPERATURE EQUATION")
    print("="*70)

    # Load TR data
    csv_files = {
        'density': find_csv_by_title('n(NS)'),
        'temperature': find_csv_by_title('T(NS)'),
        'chi_e': find_csv_by_title('AKE,AKNCE'),
        'chi_i': find_csv_by_title('AKD,AKNCD'),
        'pin': find_csv_by_title('PIN [MW'),
    }

    # Read data
    r, dens, _ = read_tr_csv(csv_files['density'])
    _, temp, _ = read_tr_csv(csv_files['temperature'])
    _, chi_e, _ = read_tr_csv(csv_files['chi_e'])
    _, chi_i, _ = read_tr_csv(csv_files['chi_i'])
    _, pin, _ = read_tr_csv(csv_files['pin'])

    ne = dens.get('nE', dens.get('NE', np.zeros_like(r)))
    Te = temp.get('TE', temp.get('Te', np.zeros_like(r)))
    Ti = temp.get('TD', temp.get('Ti', np.zeros_like(r)))

    AKE = chi_e.get('AKE', np.zeros_like(r))
    AKNCE = chi_e.get('AKNCE', np.zeros_like(r))

    PIN_e = pin.get('PIN_1', np.zeros_like(r))

    nr = len(r)
    a = A_MINOR
    r_m = r * a
    dr = r_m[1] - r_m[0] if nr > 1 else a/50

    print("\n1. Data loaded:")
    print(f"   r/a range: [{r[0]:.4f}, {r[-1]:.4f}], {nr} points")
    print(f"   Te_center = {Te[0]:.2f} keV")
    print(f"   ne_center = {ne[0]:.4f} x10^20/m^3")
    print(f"   AKE_center = {AKE[0]:.4f} m^2/s")
    print(f"   PIN_e_center = {PIN_e[0]:.4f} MW/m^3")

    # Geometry factors from trmetric.f90
    print("\n2. Geometry factors (CFEDR κ={:.2f}):".format(KAPPA))
    DVRHO = 4 * np.pi**2 * KAPPA * a**2 * R0 * r  # = 4π²κa²R₀ρ
    DV53 = DVRHO**(5/3)
    DV23 = DVRHO**(2/3)

    # At midpoint
    idx = nr // 2
    print(f"   At r/a = {r[idx]:.2f}:")
    print(f"     DVRHO = {DVRHO[idx]:.4e}")
    print(f"     DV53  = {DV53[idx]:.4e}")
    print(f"     DV23  = {DV23[idx]:.4e}")

    # Calculate the PIN term in the source
    RKEV = 1.602176634e-16  # keV to Joules
    PIN_term = (PIN_e * 1e6 / (RKEV * 1e20)) * DV53

    print("\n3. PIN term analysis:")
    print(f"   PIN_e_center = {PIN_e[0]:.4f} MW/m³ = {PIN_e[0]*1e6:.4e} W/m³")
    print(f"   PIN_term_center = {PIN_term[0]:.4e} (in matrix units)")

    # Estimate convection/diffusion contribution
    # From trexec.f90:
    # VI = FA * (AVK + AV*CC) * DV23
    # DI = FB * AK * DV23
    #
    # Typical values:
    # - AVK ~ 0 or small (heat pinch velocity) [m/s]
    # - AV ~ small (particle pinch) [m/s]
    # - CC = 1.5
    # - AK = chi (thermal diffusivity) [m²/s]

    # FA and FB are geometric factors:
    # FB = DVRHO * AR2RHO / DR² = 4π²R₀ρ / DR² (κ cancels)
    AR2RHO = 1 / (KAPPA_S * a)**2
    FB = DVRHO * AR2RHO / dr**2

    print("\n4. Diffusion term analysis (DI = FB * AK * DV23):")
    print(f"   AR2RHO = {AR2RHO:.6e}")
    print(f"   FB (at r/a=0.5) = {FB[idx]:.4e}")
    print(f"   AKE (at r/a=0.5) = {AKE[idx]:.4f} m²/s")
    DI_term = FB * AKE * DV23
    print(f"   DI_term (at r/a=0.5) = {DI_term[idx]:.4e}")

    # C83 = 8/3
    C83 = 8.0/3.0

    # RPV = n * T (pressure)
    RPV = ne * Te  # [10^20/m³ * keV]

    # The additional source term: C83 * DI * RPV
    # This represents heat flux due to diffusion
    additional_term = C83 * DI_term * RPV

    print("\n5. Additional diffusion source term (C83 * DI * RPV):")
    print(f"   C83 = {C83:.4f}")
    print(f"   RPV (at r/a=0.5) = {RPV[idx]:.4e} [10^20/m³ * keV]")
    print(f"   Additional term (at r/a=0.5) = {additional_term[idx]:.4e}")

    # Compare magnitudes
    print("\n6. Comparison of term magnitudes at r/a=0.5:")
    print(f"   PIN term:        {PIN_term[idx]:.4e}")
    print(f"   Additional term: {additional_term[idx]:.4e}")
    ratio = additional_term[idx] / PIN_term[idx] if PIN_term[idx] != 0 else 0
    print(f"   Ratio:           {ratio:.4f} ({ratio*100:.1f}%)")

    # Temperature gradient effect
    # The term DD(NEQ,NEQ-1) = FB*(AD*CC - AK)*RTV*DV23
    # represents cross-coupling between density and temperature
    print("\n7. Cross-coupling term analysis:")
    print("   DD(NEQ,NEQ-1) = FB * (AD*CC - AK) * RTV * DV23")
    print("   This term couples density gradient to temperature equation")
    print("   If AD ≈ 0 (no particle diffusion), then:")
    print("     DD(NEQ,NEQ-1) ≈ -FB * AK * RTV * DV23")

    # Estimate temperature gradient
    dTdr = np.gradient(Te, r_m)  # keV/m
    print(f"   |dTe/dr| at r/a=0.5 = {abs(dTdr[idx]):.4f} keV/m")

    RTV = 0.5 * (Te[:-1] + Te[1:])  # Average temperature
    RTV = np.append(RTV, RTV[-1])

    cross_term = FB * AKE * RTV * DV23  # Simplified (assuming AD≈0)
    print(f"   Cross-coupling magnitude: {cross_term[idx]:.4e}")

    # Summary
    print("\n" + "="*70)
    print("SUMMARY: Possible missing contributions")
    print("="*70)
    print("""
The temperature equation source term in TR includes:

1. PIN term (net power input):
   - Heating: NBI, fusion, RF, ohmic
   - Losses: radiation, charge exchange
   - Exchange: ion-electron coupling
   → This is what we use in Python reproduction

2. Diffusion source term (C83 * DI * RPV):
   - Represents heat flux contribution to energy balance
   - Magnitude: ~{:.1f}% of PIN term at r/a=0.5
   → NOT included in our Python calculation

3. Convection term (-VI * RPV):
   - Heat pinch (AVK) and particle pinch (AV×1.5)
   - Typically small unless strong pinch exists
   → NOT included, but probably small

4. Cross-coupling term (FB*(AD*CC-AK)*RTV*DV23):
   - Couples density gradient to temperature
   - Can be significant if AD ≠ 0
   → NOT included

The missing terms could account for the remaining ~3.8% error.
""".format(ratio*100))

    # Plot analysis
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    ax = axes[0, 0]
    ax.plot(r, PIN_term, 'b-', lw=2, label='PIN term')
    ax.plot(r, additional_term, 'r--', lw=2, label='C83*DI*RPV')
    ax.set_xlabel('r/a')
    ax.set_ylabel('Source term magnitude')
    ax.set_title('Comparison of Source Terms')
    ax.legend()
    ax.grid(True)
    ax.set_yscale('symlog')

    ax = axes[0, 1]
    ratio_profile = additional_term / (np.abs(PIN_term) + 1e-30) * 100
    ax.plot(r, ratio_profile, 'b-', lw=2)
    ax.axhline(0, color='k', ls='--', lw=1)
    ax.set_xlabel('r/a')
    ax.set_ylabel('Ratio [%]')
    ax.set_title('Additional Term / PIN Term')
    ax.grid(True)
    ax.set_ylim([-50, 50])

    ax = axes[1, 0]
    ax.plot(r, FB * DV23, 'b-', lw=2, label='FB × DV23')
    ax.set_xlabel('r/a')
    ax.set_ylabel('Coefficient')
    ax.set_title('Diffusion Coefficient Factor')
    ax.legend()
    ax.grid(True)

    ax = axes[1, 1]
    ax.plot(r, dTdr, 'b-', lw=2, label='dTe/dr')
    ax.set_xlabel('r/a')
    ax.set_ylabel('dT/dr [keV/m]')
    ax.set_title('Temperature Gradient')
    ax.legend()
    ax.grid(True)

    plt.tight_layout()
    plt.savefig('missing_terms_analysis.png', dpi=150)
    print("\nPlot saved: missing_terms_analysis.png")


if __name__ == '__main__':
    main()
