#!/usr/bin/env python3
"""
Generate Ar cooling function from COREDIV atomic data files.

This script reads the COREDIV atomic data files (ARTABLE1, ARPROM1, ARJON)
and generates:
1. A cooling function table L_Ar(Te)
2. Polynomial fit coefficients for Fortran implementation
3. Comparison plots

Usage:
    python generate_ar_cooling.py

Output:
    - ar_cooling_table.dat: Cooling function table
    - ar_cooling_coefficients.txt: Polynomial fit coefficients
    - ar_cooling_function.png: Comparison plot
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os

# Path to atomic data
ATOMDATA_PATH = 'ZAG/DATA/ATOMDATA'


def read_artable1(filepath):
    """
    Read ionization/recombination coefficients from ARTABLE1.

    Format: For each ionization state L (1-18):
        ALFA(K,L), BETA(K,L), GAMMA(K,L) for K=1-500

    Returns:
        alfa: shape (500, 18) - ionization rate coefficients
        beta: shape (500, 18) - recombination rate coefficients
    """
    alfa = np.zeros((500, 18))
    beta = np.zeros((500, 18))

    with open(filepath, 'r') as f:
        data = []
        for line in f:
            data.extend([float(x) for x in line.split()])

    # Reshape: 18 states * 500 temps * 3 coefficients
    data = np.array(data).reshape(18, 500, 3)

    for L in range(18):
        alfa[:, L] = data[L, :, 0]
        beta[:, L] = data[L, :, 1]

    return alfa, beta


def read_arprom1(filepath):
    """
    Read line radiation coefficients from ARPROM1.

    Format: SLIN(L, K) for L=1-19 (ionization states), K=1-500 (temperature)

    Returns:
        slin: shape (19, 500) - line radiation coefficients [W*m^3]
    """
    with open(filepath, 'r') as f:
        data = []
        for line in f:
            data.extend([float(x) for x in line.split()])

    # Total: 19 * 500 = 9500 values
    slin = np.array(data).reshape(19, 500)

    return slin


def read_arjon(filepath):
    """
    Read ionization potentials from ARJON.

    Returns:
        potjon: shape (18,) - ionization potentials [eV]
    """
    with open(filepath, 'r') as f:
        data = []
        for line in f:
            data.extend([float(x) for x in line.split()])

    return np.array(data[:18])


def temp_to_index(Te_eV):
    """Convert temperature [eV] to COREDIV index (1-500)."""
    K = 499.0 / 5.0 * np.log10(max(Te_eV, 1.0)) + 1.0
    return max(1, min(int(K), 499))


def index_to_temp(K):
    """Convert COREDIV index (1-500) to temperature [eV]."""
    return 10.0 ** ((K - 1) * 5.0 / 499.0)


def compute_coronal_equilibrium(Te_eV, alfa, beta):
    """
    Compute ionization state distribution in coronal equilibrium.

    In coronal equilibrium:
        n_z * S_z = n_{z+1} * alpha_{z+1}
        where S_z = ne * alfa_z (ionization rate)
              alpha_z = ne * beta_z (recombination rate)

    Args:
        Te_eV: electron temperature [eV]
        alfa: ionization rate coefficients (500, 18)
        beta: recombination rate coefficients (500, 18)

    Returns:
        f: fractional abundance for each ionization state (19,)
           f[0] = neutral, f[1] = Ar+, ..., f[18] = Ar18+
    """
    K = temp_to_index(Te_eV)

    f = np.zeros(19)
    f[0] = 1.0

    for z in range(18):
        # f[z+1] / f[z] = alfa[z] / beta[z+1]
        # Note: alfa index starts from 0 (for Ar0 -> Ar+)
        ratio = alfa[K, z] / max(beta[K, z], 1e-30)
        f[z + 1] = f[z] * ratio

    # Normalize
    f = f / np.sum(f)

    return f


def compute_cooling_function(Te_eV, f_z, slin):
    """
    Compute total cooling function L_z(Te).

    L_z = sum_z (f_z * SLIN_z)

    Args:
        Te_eV: electron temperature [eV]
        f_z: fractional abundance (19,)
        slin: line radiation coefficients (19, 500)

    Returns:
        L_total: total cooling function [W*m^3]
    """
    K = temp_to_index(Te_eV)

    L_total = 0.0
    for z in range(1, 19):  # Only ionized states contribute
        L_total += f_z[z] * slin[z, K]

    return L_total


def generate_cooling_table(alfa, beta, slin, n_points=200):
    """
    Generate cooling function table over temperature range.

    Args:
        alfa, beta, slin: atomic data arrays
        n_points: number of temperature points

    Returns:
        Te_keV: temperature array [keV]
        L_Ar: cooling function array [W*m^3]
        f_z_all: ionization state distributions at each temperature
    """
    # Temperature range: 1 eV to 100 keV
    Te_eV = np.logspace(0, 5, n_points)
    Te_keV = Te_eV / 1000.0

    L_Ar = np.zeros(n_points)
    f_z_all = np.zeros((n_points, 19))

    for i, T in enumerate(Te_eV):
        f_z = compute_coronal_equilibrium(T, alfa, beta)
        f_z_all[i, :] = f_z
        L_Ar[i] = compute_cooling_function(T, f_z, slin)

    return Te_keV, L_Ar, f_z_all


def fit_cooling_function(Te_keV, L_Ar):
    """
    Fit cooling function with piecewise polynomials in log-log space.

    Following the TR code convention:
        log10(L) = A0 + A1*x + A2*x^2 + A3*x^3 + A4*x^4 + A5*x^5
        where x = log10(Te[keV])

    Temperature ranges (matching TRRPC/TRRPFE):
        Range 1: Te <= 0.003 keV (3 eV)
        Range 2: 0.003 < Te <= 0.02 keV (20 eV)
        Range 3: 0.02 < Te <= 0.2 keV (200 eV)
        Range 4: 0.2 < Te <= 2 keV (2000 eV)
        Range 5: 2 < Te <= 20 keV
        Range 6: Te > 20 keV

    Returns:
        coeffs: dictionary with polynomial coefficients for each range
    """
    # Define temperature ranges
    ranges = [
        (0.001, 0.003, 'Range1'),
        (0.003, 0.02, 'Range2'),
        (0.02, 0.2, 'Range3'),
        (0.2, 2.0, 'Range4'),
        (2.0, 20.0, 'Range5'),
        (20.0, 100.0, 'Range6'),
    ]

    coeffs = {}

    for Te_min, Te_max, name in ranges:
        mask = (Te_keV >= Te_min) & (Te_keV <= Te_max)
        if np.sum(mask) < 6:
            continue

        x = np.log10(Te_keV[mask])
        y = np.log10(L_Ar[mask] + 1e-40)

        # Fit polynomial
        try:
            p = np.polyfit(x, y, 5)
            coeffs[name] = {
                'Te_min': Te_min,
                'Te_max': Te_max,
                'A5': p[0], 'A4': p[1], 'A3': p[2],
                'A2': p[3], 'A1': p[4], 'A0': p[5]
            }
        except:
            print(f"Warning: Could not fit {name}")

    return coeffs


def generate_fortran_code(coeffs):
    """Generate Fortran code for TRRPAR function."""

    code = '''
!     ***********************************************************
!           RADIATION POWER - Ar (Argon)
!           Generated from COREDIV atomic data
!     ***********************************************************

      FUNCTION TRRPAR(TE)

      USE TRCOMM,ONLY: rkind
      IMPLICIT NONE
      REAL(rkind):: TE, TEL, ARG, TRRPAR

'''

    # Add coefficients
    for i, (name, c) in enumerate(coeffs.items(), 1):
        code += f'''      ! {name}: {c["Te_min"]} <= Te <= {c["Te_max"]} keV
      REAL(rkind):: AAR{i}0={c["A0"]:15.8E}, AAR{i}1={c["A1"]:15.8E}, AAR{i}2={c["A2"]:15.8E}, &
                    AAR{i}3={c["A3"]:15.8E}, AAR{i}4={c["A4"]:15.8E}, AAR{i}5={c["A5"]:15.8E}

'''

    code += '''
      IF(TE.LE.0.D0) THEN
         TRRPAR = 0.D0
      ELSE
         TEL = LOG10(TE)
'''

    # Add conditional statements
    for i, (name, c) in enumerate(coeffs.items(), 1):
        if i == 1:
            code += f'''         IF(TE.LE.{c["Te_max"]}D0) THEN
            IF(TE.LT.{c["Te_min"]}D0) TEL = LOG10({c["Te_min"]}D0)
            ARG = AAR{i}0 + AAR{i}1*TEL + AAR{i}2*TEL**2 + &
                  AAR{i}3*TEL**3 + AAR{i}4*TEL**4 + AAR{i}5*TEL**5
'''
        elif i < len(coeffs):
            code += f'''         ELSEIF(TE.LE.{c["Te_max"]}D0) THEN
            ARG = AAR{i}0 + AAR{i}1*TEL + AAR{i}2*TEL**2 + &
                  AAR{i}3*TEL**3 + AAR{i}4*TEL**4 + AAR{i}5*TEL**5
'''
        else:
            code += f'''         ELSE
            ARG = AAR{i}0 + AAR{i}1*TEL + AAR{i}2*TEL**2 + &
                  AAR{i}3*TEL**3 + AAR{i}4*TEL**4 + AAR{i}5*TEL**5
'''

    code += '''         END IF
         TRRPAR = 10.D0**(ARG - 13.D0)
      END IF

      RETURN
      END FUNCTION TRRPAR
'''

    return code


def plot_results(Te_keV, L_Ar, f_z_all, coeffs, potjon):
    """Create diagnostic plots."""

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # 1. Cooling function
    ax = axes[0, 0]
    ax.loglog(Te_keV * 1000, L_Ar, 'b-', linewidth=2, label='L_Ar (from data)')

    # Plot fitted function
    for name, c in coeffs.items():
        mask = (Te_keV >= c['Te_min']) & (Te_keV <= c['Te_max'])
        Te_fit = Te_keV[mask]
        x = np.log10(Te_fit)
        L_fit = 10 ** (c['A0'] + c['A1']*x + c['A2']*x**2 +
                       c['A3']*x**3 + c['A4']*x**4 + c['A5']*x**5)
        ax.loglog(Te_fit * 1000, L_fit, 'r--', linewidth=1)

    ax.set_xlabel('Te [eV]')
    ax.set_ylabel('L_Ar [W m³]')
    ax.set_title('Ar Cooling Function (Coronal Equilibrium)')
    ax.legend()
    ax.grid(True, which='both', alpha=0.3)
    ax.set_xlim([1, 1e5])

    # 2. Ionization state distribution
    ax = axes[0, 1]
    Te_plot = [10, 100, 1000, 10000, 50000]  # eV
    colors = plt.cm.viridis(np.linspace(0, 1, len(Te_plot)))

    for T, color in zip(Te_plot, colors):
        idx = np.argmin(np.abs(Te_keV * 1000 - T))
        ax.bar(np.arange(19), f_z_all[idx, :], alpha=0.5,
               label=f'Te = {T} eV', color=color)

    ax.set_xlabel('Ionization state z')
    ax.set_ylabel('Fractional abundance')
    ax.set_title('Ar Ionization State Distribution')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 3. Mean charge state
    ax = axes[1, 0]
    z_mean = np.sum(f_z_all * np.arange(19), axis=1)
    ax.semilogx(Te_keV * 1000, z_mean, 'b-', linewidth=2)
    ax.set_xlabel('Te [eV]')
    ax.set_ylabel('Mean charge state <Z>')
    ax.set_title('Ar Mean Ionization State')
    ax.grid(True, alpha=0.3)
    ax.set_xlim([1, 1e5])
    ax.set_ylim([0, 18])

    # 4. Ionization potentials
    ax = axes[1, 1]
    ax.bar(np.arange(1, 19), potjon, color='steelblue')
    ax.set_xlabel('Ionization state z')
    ax.set_ylabel('Ionization potential [eV]')
    ax.set_title('Ar Ionization Potentials')
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('ar_cooling_function.png', dpi=150)
    print("Saved: ar_cooling_function.png")


def main():
    print("=" * 60)
    print("Generating Ar Cooling Function from COREDIV Atomic Data")
    print("=" * 60)

    # Read atomic data
    print("\n1. Reading atomic data files...")

    artable1_path = os.path.join(ATOMDATA_PATH, 'ARTABLE1')
    arprom1_path = os.path.join(ATOMDATA_PATH, 'ARPROM1')
    arjon_path = os.path.join(ATOMDATA_PATH, 'ARJON')

    if not os.path.exists(artable1_path):
        print(f"Error: {artable1_path} not found")
        return

    alfa, beta = read_artable1(artable1_path)
    print(f"   ARTABLE1: alfa shape = {alfa.shape}, beta shape = {beta.shape}")

    slin = read_arprom1(arprom1_path)
    print(f"   ARPROM1: slin shape = {slin.shape}")

    potjon = read_arjon(arjon_path)
    print(f"   ARJON: {len(potjon)} ionization potentials")

    # Generate cooling table
    print("\n2. Computing coronal equilibrium cooling function...")
    Te_keV, L_Ar, f_z_all = generate_cooling_table(alfa, beta, slin)
    print(f"   Temperature range: {Te_keV[0]*1000:.1f} eV to {Te_keV[-1]*1000:.1f} eV")
    print(f"   L_Ar range: {L_Ar.min():.2e} to {L_Ar.max():.2e} W*m³")

    # Save cooling table
    print("\n3. Saving cooling function table...")
    np.savetxt('ar_cooling_table.dat',
               np.column_stack([Te_keV, L_Ar]),
               header='Te[keV]  L_Ar[W*m^3]',
               fmt='%.6e')
    print("   Saved: ar_cooling_table.dat")

    # Fit polynomial coefficients
    print("\n4. Fitting polynomial coefficients...")
    coeffs = fit_cooling_function(Te_keV, L_Ar)

    for name, c in coeffs.items():
        print(f"   {name} ({c['Te_min']:.3f} - {c['Te_max']:.3f} keV):")
        print(f"      A0 = {c['A0']:.6e}")

    # Save coefficients
    with open('ar_cooling_coefficients.txt', 'w') as f:
        f.write("# Ar Cooling Function Polynomial Coefficients\n")
        f.write("# log10(L) = A0 + A1*x + A2*x^2 + A3*x^3 + A4*x^4 + A5*x^5\n")
        f.write("# where x = log10(Te[keV])\n\n")

        for name, c in coeffs.items():
            f.write(f"# {name}: {c['Te_min']} <= Te <= {c['Te_max']} keV\n")
            f.write(f"A0 = {c['A0']:.10e}\n")
            f.write(f"A1 = {c['A1']:.10e}\n")
            f.write(f"A2 = {c['A2']:.10e}\n")
            f.write(f"A3 = {c['A3']:.10e}\n")
            f.write(f"A4 = {c['A4']:.10e}\n")
            f.write(f"A5 = {c['A5']:.10e}\n\n")

    print("   Saved: ar_cooling_coefficients.txt")

    # Generate Fortran code
    print("\n5. Generating Fortran code...")
    fortran_code = generate_fortran_code(coeffs)

    with open('trrpar_function.f90', 'w') as f:
        f.write(fortran_code)

    print("   Saved: trrpar_function.f90")

    # Create plots
    print("\n6. Creating diagnostic plots...")
    plot_results(Te_keV, L_Ar, f_z_all, coeffs, potjon)

    # Summary
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    print(f"Peak cooling at Te ~ {Te_keV[np.argmax(L_Ar)]*1000:.0f} eV")
    print(f"L_Ar(peak) = {L_Ar.max():.2e} W*m³")
    print("\nTo use in TR code:")
    print("1. Copy trrpar_function.f90 content to tradat.f90")
    print("2. Add PLAR = ANE * ANAR(NR) * TRRPAR(TE) * 1.D40 in TRLOSS")
    print("3. Add PRL(NR) = PLFE + PLC + PLAR")


if __name__ == '__main__':
    main()
