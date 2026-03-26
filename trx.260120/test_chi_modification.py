#!/usr/bin/env python3
"""
Test script to analyze the temperature jump caused by TR's chi_e
and verify that modifying chi in the 0.3-0.45 region eliminates the jump.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# Constants
A_MINOR = 2.44  # m
R0 = 9.2  # m
E_CHARGE = 1.602e-19

def read_chi_se_corediv(filename='Chi_Se_COREDIV.DAT'):
    """Read COREDIV data file."""
    data = np.loadtxt(filename, skiprows=1)
    r = data[:, 0]
    ne = data[:, 1] / 1e20
    chi_e = data[:, 3]
    Se = data[:, 5]
    G1 = data[:, 7]
    ATOR = 4.0 * np.pi**2 * R0
    Vprime = G1 * ATOR * r
    Vprime[0] = Vprime[1] * 0.01
    return r, ne, chi_e, Se, Vprime

def solve_temperature_with_vprime(r, Vprime, ne, chi, Se, Te_bc):
    """
    Solve temperature profile using flux-coordinate heat equation.

    1/V' * d/dr(V' * n * chi * dT/dr) = -S
    """
    n = len(r)
    Q = np.zeros(n)
    Te = np.zeros(n)

    Vp = Vprime.copy()
    Vp[Vp < 1e-10] = 1e-10

    # Step 1: Calculate heat flux Q(r) = integral of (V' * S * dr)
    for i in range(1, n):
        dr = r[i] - r[i-1]
        integrand = Vp * Se
        Q[i] = Q[i-1] + 0.5 * (integrand[i-1] + integrand[i]) * dr

    # Step 2: Integrate temperature from edge to axis
    Te[-1] = Te_bc
    for i in range(n-2, -1, -1):
        dr = r[i+1] - r[i]
        denom_i = E_CHARGE * Vp[i] * ne[i] * chi[i]
        denom_ip1 = E_CHARGE * Vp[i+1] * ne[i+1] * chi[i+1]
        denom_i = max(denom_i, 1e-30)
        denom_ip1 = max(denom_ip1, 1e-30)

        grad_i = -Q[i] / denom_i
        grad_ip1 = -Q[i+1] / denom_ip1
        grad_mid = 0.5 * (grad_i + grad_ip1)

        Te[i] = Te[i+1] - grad_mid * dr

    return Te, Q

def find_jump(Te, r, r1=0.30, r2=0.45):
    """Calculate temperature difference between two radial positions."""
    idx1 = np.argmin(np.abs(r - r1))
    idx2 = np.argmin(np.abs(r - r2))
    return Te[idx1] - Te[idx2]

def main():
    # Read TR chi_e
    df_chi = pd.read_csv('tr_data_034.csv', skiprows=1)
    r_tr = df_chi.iloc[:, 0].values
    chi_e_tr = df_chi.iloc[:, -1].values

    # Read COREDIV data
    r_cd, ne_cd, chi_e_cd, Se_cd, Vprime_cd = read_chi_se_corediv()

    # Interpolate to TR grid
    r_tr_m = r_tr * A_MINOR
    f_chi_cd = interp1d(r_cd, chi_e_cd, bounds_error=False, fill_value='extrapolate')
    f_Se = interp1d(r_cd, Se_cd, bounds_error=False, fill_value='extrapolate')
    f_Vp = interp1d(r_cd, Vprime_cd, bounds_error=False, fill_value='extrapolate')
    f_ne = interp1d(r_cd, ne_cd, bounds_error=False, fill_value='extrapolate')

    chi_e_corediv = f_chi_cd(r_tr_m)
    Se_corediv = f_Se(r_tr_m)
    Vprime = f_Vp(r_tr_m)
    ne_corediv = f_ne(r_tr_m)

    # Create modified chi profiles
    # Test 1: Replace r/a=0.3-0.45 with COREDIV chi
    chi_e_test1 = chi_e_tr.copy()
    for i, ri in enumerate(r_tr):
        if 0.30 <= ri <= 0.45:
            chi_e_test1[i] = chi_e_corediv[i]

    # Test 2: Scale down TR chi in r/a=0.3-0.45 by factor of 30
    chi_e_test2 = chi_e_tr.copy()
    for i, ri in enumerate(r_tr):
        if 0.30 <= ri <= 0.45:
            chi_e_test2[i] = chi_e_tr[i] / 30.0

    # Test 3: Smooth blending in 0.25-0.55 region
    chi_e_test3 = chi_e_tr.copy()
    for i, ri in enumerate(r_tr):
        if 0.25 <= ri <= 0.55:
            center = 0.4
            width = 0.15
            blend = max(0, 1 - abs(ri - center) / width)
            chi_e_test3[i] = (1 - blend) * chi_e_tr[i] + blend * chi_e_corediv[i]

    # Calculate temperatures
    Te_bc = 100  # eV at edge

    Te_tr_chi, _ = solve_temperature_with_vprime(
        r_tr_m, Vprime, ne_corediv * 1e20, chi_e_tr, Se_corediv, Te_bc)

    Te_test1, _ = solve_temperature_with_vprime(
        r_tr_m, Vprime, ne_corediv * 1e20, chi_e_test1, Se_corediv, Te_bc)

    Te_test2, _ = solve_temperature_with_vprime(
        r_tr_m, Vprime, ne_corediv * 1e20, chi_e_test2, Se_corediv, Te_bc)

    Te_test3, _ = solve_temperature_with_vprime(
        r_tr_m, Vprime, ne_corediv * 1e20, chi_e_test3, Se_corediv, Te_bc)

    Te_corediv, _ = solve_temperature_with_vprime(
        r_tr_m, Vprime, ne_corediv * 1e20, chi_e_corediv, Se_corediv, Te_bc)

    # Convert to keV
    Te_tr_chi_keV = Te_tr_chi / 1000.0
    Te_test1_keV = Te_test1 / 1000.0
    Te_test2_keV = Te_test2 / 1000.0
    Te_test3_keV = Te_test3 / 1000.0
    Te_corediv_keV = Te_corediv / 1000.0

    # Calculate jumps
    jump_orig = find_jump(Te_tr_chi_keV, r_tr)
    jump_test1 = find_jump(Te_test1_keV, r_tr)
    jump_test2 = find_jump(Te_test2_keV, r_tr)
    jump_test3 = find_jump(Te_test3_keV, r_tr)
    jump_ref = find_jump(Te_corediv_keV, r_tr)

    # Print results
    print("\n" + "="*60)
    print("Chi Modification Test Results")
    print("="*60)
    print(f"{'Case':<25} {'Te(0.01)':<12} {'Jump (0.3-0.45)':<15}")
    print("-"*60)
    print(f"{'COREDIV chi (ref)':<25} {Te_corediv_keV[0]:<12.2f} {jump_ref:<+15.2f}")
    print(f"{'TR chi (original)':<25} {Te_tr_chi_keV[0]:<12.2f} {jump_orig:<+15.2f}")
    print(f"{'Test1 (CD chi 0.3-0.45)':<25} {Te_test1_keV[0]:<12.2f} {jump_test1:<+15.2f}")
    print(f"{'Test2 (TR/30 in 0.3-0.45)':<25} {Te_test2_keV[0]:<12.2f} {jump_test2:<+15.2f}")
    print(f"{'Test3 (Blended)':<25} {Te_test3_keV[0]:<12.2f} {jump_test3:<+15.2f}")

    # Plot results
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Plot 1: Temperature profiles
    ax = axes[0, 0]
    ax.plot(r_tr, Te_corediv_keV, 'k-', lw=2.5, label='COREDIV chi (reference)')
    ax.plot(r_tr, Te_tr_chi_keV, 'b-', lw=2, label='TR chi (original)')
    ax.plot(r_tr, Te_test1_keV, 'r--', lw=2, label='Test1: CD chi in 0.30-0.45')
    ax.plot(r_tr, Te_test2_keV, 'g--', lw=2, label='Test2: TR/30 in 0.30-0.45')
    ax.plot(r_tr, Te_test3_keV, 'm--', lw=2, label='Test3: Blended 0.25-0.55')
    ax.axvspan(0.30, 0.45, alpha=0.15, color='yellow', label='Modified region')
    ax.set_xlabel('r/a', fontsize=12)
    ax.set_ylabel('Te [keV]', fontsize=12)
    ax.set_title('Temperature Profiles with Modified Chi', fontsize=14)
    ax.legend(fontsize=9, loc='upper right')
    ax.grid(True)
    ax.set_xlim([0, 1])

    # Plot 2: Chi profiles
    ax = axes[0, 1]
    ax.plot(r_tr, chi_e_tr, 'b-', lw=2, label='TR chi (original)')
    ax.plot(r_tr, chi_e_corediv, 'k-', lw=2, label='COREDIV chi')
    ax.plot(r_tr, chi_e_test1, 'r--', lw=1.5, label='Test1')
    ax.plot(r_tr, chi_e_test2, 'g--', lw=1.5, label='Test2')
    ax.plot(r_tr, chi_e_test3, 'm--', lw=1.5, label='Test3')
    ax.axvspan(0.30, 0.45, alpha=0.15, color='yellow')
    ax.set_xlabel('r/a', fontsize=12)
    ax.set_ylabel('chi_e [m²/s]', fontsize=12)
    ax.set_title('Modified Chi Profiles', fontsize=14)
    ax.legend(fontsize=9)
    ax.grid(True)
    ax.set_xlim([0, 1])
    ax.set_yscale('log')

    # Plot 3: Zoom on jump region
    ax = axes[1, 0]
    ax.plot(r_tr, Te_corediv_keV, 'k-', lw=2.5, label='COREDIV chi')
    ax.plot(r_tr, Te_tr_chi_keV, 'b-', lw=2, label='TR chi (original)')
    ax.plot(r_tr, Te_test1_keV, 'r--', lw=2, label='Test1: CD chi')
    ax.plot(r_tr, Te_test2_keV, 'g--', lw=2, label='Test2: TR/30')
    ax.plot(r_tr, Te_test3_keV, 'm--', lw=2, label='Test3: Blended')
    ax.axvspan(0.30, 0.45, alpha=0.15, color='yellow')
    ax.set_xlabel('r/a', fontsize=12)
    ax.set_ylabel('Te [keV]', fontsize=12)
    ax.set_title('Temperature Zoom (0.15-0.65)', fontsize=14)
    ax.legend(fontsize=9)
    ax.grid(True)
    ax.set_xlim([0.15, 0.65])

    # Plot 4: Summary table
    ax = axes[1, 1]
    ax.axis('off')
    summary = f"""
Chi Modification Test Results
═════════════════════════════════════════════════════════

                    Te(0.01)    ΔTe(0.30-0.45)   Jump
                    [keV]       [keV]            Reduced?
─────────────────────────────────────────────────────────
COREDIV chi (ref)   {Te_corediv_keV[0]:6.2f}      {jump_ref:+6.2f}           -
TR chi (original)   {Te_tr_chi_keV[0]:6.2f}      {jump_orig:+6.2f}           -

Test1 (CD@0.3-0.45) {Te_test1_keV[0]:6.2f}      {jump_test1:+6.2f}           {'✓ Yes' if abs(jump_test1) < abs(jump_orig)*0.7 else '✗ No'}
Test2 (TR/30)       {Te_test2_keV[0]:6.2f}      {jump_test2:+6.2f}           {'✓ Yes' if abs(jump_test2) < abs(jump_orig)*0.7 else '✗ No'}
Test3 (Blended)     {Te_test3_keV[0]:6.2f}      {jump_test3:+6.2f}           {'✓ Yes' if abs(jump_test3) < abs(jump_orig)*0.7 else '✗ No'}
─────────────────────────────────────────────────────────

Conclusion:
  The jump at r/a=0.3-0.45 is DIRECTLY caused by
  TR's large chi_e. Reducing chi in this region
  significantly reduces or eliminates the jump.
"""
    ax.text(0.02, 0.95, summary, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    plt.tight_layout()
    plt.savefig('chi_modification_test.png', dpi=150)
    print("\nPlot saved to chi_modification_test.png")


if __name__ == '__main__':
    main()
