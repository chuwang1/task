#!/usr/bin/env python3
"""
Analyze energy balance from TR output CSV files.

This script reads density, source terms, and transport coefficients
from TR output files and analyzes the energy balance at a given time point.

Available data from CSV files:
- Density: tr_data_017 (nE, nD, nT, nA vs r)
- Temperature: tr_data_019 (TE, TD, TT, TA vs r)
- Power input: tr_data_021 (POH, PNB, PNF, PRSUM, PRF vs r)
- Power loss: tr_data_022 (PRSUM, PRB, PRC, PRL, PCX, PIE vs r)
- Transport coefficients: tr_data_023/024 (AKE, AKD vs r)
- Heat input by species: tr_data_037 (PIN_1-4 vs r)

Energy equation (simplified 1D):
  d(nT)/dt = -1/r * d/dr(r * n * chi * dT/dr) + S

In steady state:
  1/r * d/dr(r * n * chi * dT/dr) = S / (1.5 * kB)

"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import integrate
from scipy.interpolate import interp1d

# Physical parameters for CFEDR
R0 = 7.8   # Major radius [m]
A_MINOR = 2.44  # Minor radius [m]

def read_tr_csv(filename):
    """Read TR output CSV file, return x (r/a) and data columns."""
    with open(filename, 'r') as f:
        title = f.readline().strip()

    df = pd.read_csv(filename, skiprows=1, index_col=False)
    df = df.loc[:, ~df.columns.str.contains('^Unnamed')]

    x = df.iloc[:, 0].values
    data = {col: df[col].values for col in df.columns[1:]}

    return x, data, title


def read_omfit_profiles(filename='omfit_profiles_for_tr.csv'):
    """Read OMFIT profiles exported CSV file."""
    if not os.path.exists(filename):
        return None
    df = pd.read_csv(filename)
    return df


def integrate_power(r, p, R0=R0, a=A_MINOR):
    """
    Integrate power density over plasma volume.

    Args:
        r: normalized radius (r/a)
        p: power density [MW/m³]
        R0: major radius [m]
        a: minor radius [m]

    Returns:
        Total integrated power [MW]
    """
    r_real = r * a  # Convert to real radius [m]
    # Volume element: dV = 2*pi*R0 * 2*pi*r*dr (torus)
    dV = 2 * np.pi * R0 * 2 * np.pi * r_real
    return np.trapz(p * dV, r_real)


def compare_with_omfit(r_tr, data, omfit_file='omfit_profiles_for_tr.csv'):
    """
    Compare TR results with OMFIT profiles.

    Args:
        r_tr: TR radial grid (r/a)
        data: TR data dictionary
        omfit_file: path to OMFIT profiles CSV
    """
    omfit = read_omfit_profiles(omfit_file)
    if omfit is None:
        print(f"\nWarning: OMFIT file {omfit_file} not found, skipping comparison")
        return

    print("\n" + "="*70)
    print("COMPARISON: TR Simulation vs OMFIT (CFEDR)")
    print("="*70)

    # Get OMFIT radial grid
    r_omfit = omfit['r/a'].values

    # Interpolate OMFIT data to TR grid
    def interp_to_tr(omfit_data):
        f = interp1d(r_omfit, omfit_data, kind='linear', fill_value='extrapolate')
        return f(r_tr)

    # Get TR power data
    if 'power_in' in data and 'power_loss' in data:
        PNF_tr = data['power_in']['values'].get('PNF', np.zeros_like(r_tr))
        PRF_tr = data['power_in']['values'].get('PRF', np.zeros_like(r_tr))
        POH_tr = data['power_in']['values'].get('POH', np.zeros_like(r_tr))
        PRSUM_tr = data['power_loss']['values'].get('PRSUM', np.zeros_like(r_tr))
    else:
        print("Error: Power data not found in TR results")
        return

    # Get OMFIT power data (already in MW/m³)
    qfuse_omfit = omfit['qfuse'].values if 'qfuse' in omfit.columns else np.zeros_like(r_omfit)
    qrfe_omfit = omfit['qrfe'].values if 'qrfe' in omfit.columns else np.zeros_like(r_omfit)
    qrad_omfit = omfit['qrad'].values if 'qrad' in omfit.columns else np.zeros_like(r_omfit)

    # Interpolate OMFIT to TR grid for comparison
    qfuse_interp = interp_to_tr(qfuse_omfit)
    qrfe_interp = interp_to_tr(qrfe_omfit)
    qrad_interp = interp_to_tr(qrad_omfit)

    # Calculate integrated powers
    P_fuse_tr = integrate_power(r_tr, PNF_tr)
    P_rf_tr = integrate_power(r_tr, PRF_tr)
    P_oh_tr = integrate_power(r_tr, POH_tr)
    P_rad_tr = integrate_power(r_tr, PRSUM_tr)

    P_fuse_omfit = integrate_power(r_omfit, qfuse_omfit)
    P_rf_omfit = integrate_power(r_omfit, qrfe_omfit)
    P_rad_omfit = integrate_power(r_omfit, qrad_omfit)

    # Print comparison table
    print("\n--- POWER DENSITY [MW/m³] ---")
    print(f"{'Source':<15} {'TR center':>12} {'OMFIT center':>12} {'TR avg':>10} {'OMFIT avg':>10}")
    print("-" * 60)
    print(f"{'Fusion':<15} {PNF_tr[0]:>12.4f} {qfuse_interp[0]:>12.4f} {np.mean(PNF_tr):>10.4f} {np.mean(qfuse_omfit):>10.4f}")
    print(f"{'RF':<15} {PRF_tr[0]:>12.4f} {qrfe_interp[0]:>12.4f} {np.mean(PRF_tr):>10.4f} {np.mean(qrfe_omfit):>10.4f}")
    print(f"{'Radiation':<15} {PRSUM_tr[0]:>12.4f} {qrad_interp[0]:>12.4f} {np.mean(PRSUM_tr):>10.4f} {np.mean(qrad_omfit):>10.4f}")

    print("\n--- INTEGRATED POWER [MW] ---")
    print(f"{'Source':<15} {'TR':>12} {'OMFIT':>12} {'Ratio (TR/OMFIT)':>18}")
    print("-" * 60)
    ratio_fuse = P_fuse_tr / P_fuse_omfit if P_fuse_omfit > 0 else float('inf')
    ratio_rf = P_rf_tr / P_rf_omfit if P_rf_omfit > 0 else float('inf')
    ratio_rad = P_rad_tr / P_rad_omfit if P_rad_omfit > 0 else float('inf')
    print(f"{'Fusion':<15} {P_fuse_tr:>12.1f} {P_fuse_omfit:>12.1f} {ratio_fuse:>18.2f}")
    print(f"{'RF':<15} {P_rf_tr:>12.1f} {P_rf_omfit:>12.1f} {ratio_rf:>18.2f}")
    print(f"{'Ohmic':<15} {P_oh_tr:>12.1f} {'-':>12} {'-':>18}")
    print(f"{'Radiation':<15} {P_rad_tr:>12.1f} {P_rad_omfit:>12.1f} {ratio_rad:>18.2f}")

    # Power balance
    P_in_tr = P_fuse_tr + P_rf_tr + P_oh_tr
    P_in_omfit = P_fuse_omfit + P_rf_omfit

    print("\n--- POWER BALANCE [MW] ---")
    print(f"{'Item':<20} {'TR':>12} {'OMFIT':>12}")
    print("-" * 50)
    print(f"{'Total Power IN':<20} {P_in_tr:>12.1f} {P_in_omfit:>12.1f}")
    print(f"{'Radiation Loss':<20} {P_rad_tr:>12.1f} {P_rad_omfit:>12.1f}")
    print(f"{'Net (IN - Rad)':<20} {P_in_tr - P_rad_tr:>12.1f} {P_in_omfit - P_rad_omfit:>12.1f}")

    # Analysis
    print("\n--- ANALYSIS ---")
    if P_fuse_tr > P_fuse_omfit * 1.5:
        print(f"! Fusion power too HIGH: TR={P_fuse_tr:.0f} MW vs OMFIT={P_fuse_omfit:.0f} MW")
        print(f"  Likely cause: Temperature too high (fusion ~ T^4)")
    elif P_fuse_tr < P_fuse_omfit * 0.5:
        print(f"! Fusion power too LOW: TR={P_fuse_tr:.0f} MW vs OMFIT={P_fuse_omfit:.0f} MW")
        print(f"  Likely cause: Temperature too low")
    else:
        print(f"✓ Fusion power reasonable: TR={P_fuse_tr:.0f} MW vs OMFIT={P_fuse_omfit:.0f} MW")

    if P_rad_tr < P_rad_omfit * 0.5:
        print(f"! Radiation too LOW: TR={P_rad_tr:.0f} MW vs OMFIT={P_rad_omfit:.0f} MW")
        print(f"  Missing: impurity line radiation (Ar), synchrotron")
    elif P_rad_tr > P_rad_omfit * 1.5:
        print(f"! Radiation too HIGH: TR={P_rad_tr:.0f} MW vs OMFIT={P_rad_omfit:.0f} MW")
    else:
        print(f"✓ Radiation reasonable: TR={P_rad_tr:.0f} MW vs OMFIT={P_rad_omfit:.0f} MW")

    # Plot comparison
    plot_power_comparison(r_tr, PNF_tr, PRF_tr, PRSUM_tr,
                          r_omfit, qfuse_omfit, qrfe_omfit, qrad_omfit)

    return {
        'P_fuse_tr': P_fuse_tr, 'P_fuse_omfit': P_fuse_omfit,
        'P_rf_tr': P_rf_tr, 'P_rf_omfit': P_rf_omfit,
        'P_rad_tr': P_rad_tr, 'P_rad_omfit': P_rad_omfit
    }


def plot_power_comparison(r_tr, PNF_tr, PRF_tr, PRSUM_tr,
                          r_omfit, qfuse_omfit, qrfe_omfit, qrad_omfit):
    """Plot power profile comparison between TR and OMFIT."""
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))

    # Fusion power
    ax = axes[0]
    ax.plot(r_tr, PNF_tr, 'b-', linewidth=2, label='TR (PNF)')
    ax.plot(r_omfit, qfuse_omfit, 'r--', linewidth=2, label='OMFIT (qfuse)')
    ax.set_xlabel('r/a')
    ax.set_ylabel('P [MW/m³]')
    ax.set_title('Fusion Heating')
    ax.legend()
    ax.grid(True)
    ax.set_xlim([0, 1])

    # RF power
    ax = axes[1]
    ax.plot(r_tr, PRF_tr, 'b-', linewidth=2, label='TR (PRF)')
    ax.plot(r_omfit, qrfe_omfit, 'r--', linewidth=2, label='OMFIT (qrfe)')
    ax.set_xlabel('r/a')
    ax.set_ylabel('P [MW/m³]')
    ax.set_title('RF Heating')
    ax.legend()
    ax.grid(True)
    ax.set_xlim([0, 1])

    # Radiation
    ax = axes[2]
    ax.plot(r_tr, PRSUM_tr, 'b-', linewidth=2, label='TR (PRSUM)')
    ax.plot(r_omfit, qrad_omfit, 'r--', linewidth=2, label='OMFIT (qrad)')
    ax.set_xlabel('r/a')
    ax.set_ylabel('P [MW/m³]')
    ax.set_title('Radiation Loss')
    ax.legend()
    ax.grid(True)
    ax.set_xlim([0, 1])

    plt.tight_layout()
    plt.savefig('power_comparison_tr_omfit.png', dpi=150)
    print(f"\nPower comparison plot saved to power_comparison_tr_omfit.png")


def analyze_energy_balance(time_index=None, compare_omfit=True):
    """
    Analyze energy balance at a specific time point.

    For radial profiles (files 017-044), data is at final time.
    For time-space evolution (files 045+), can select time index.
    """

    print("="*60)
    print("TR Energy Balance Analysis")
    print("="*60)

    # Read radial profiles (final time snapshot)
    files_to_read = {
        'density': 'tr_data_017.csv',      # nE, nD, nT, nA
        'temperature': 'tr_data_019.csv',  # TE, TD, TT, TA
        'power_in': 'tr_data_021.csv',     # POH, PNB, PNF, PRSUM, PRF
        'power_loss': 'tr_data_022.csv',   # PRSUM, PRB, PRC, PRL, PCX, PIE
        'chi_e': 'tr_data_023.csv',        # AKE, AKNCE, AKDWE
        'chi_i': 'tr_data_024.csv',        # AKD, AKNCD, AKDWD
        'pin_species': 'tr_data_037.csv',  # PIN_1-4
    }

    data = {}
    for key, filename in files_to_read.items():
        if os.path.exists(filename):
            r, vals, title = read_tr_csv(filename)
            data[key] = {'r': r, 'values': vals, 'title': title}
            print(f"Read {filename}: {list(vals.keys())}")
        else:
            print(f"Warning: {filename} not found")

    if 'density' not in data or 'temperature' not in data:
        print("Error: Required files not found")
        return

    r = data['density']['r']
    dr = r[1] - r[0] if len(r) > 1 else 0.02

    # Extract key quantities
    ne = data['density']['values'].get('nE', np.zeros_like(r))
    nD = data['density']['values'].get('nD', np.zeros_like(r))

    Te = data['temperature']['values'].get('TE', np.zeros_like(r))
    TD = data['temperature']['values'].get('TD', np.zeros_like(r))

    print("\n" + "="*60)
    print("Profile Summary")
    print("="*60)

    print(f"\nRadial grid: {len(r)} points, r/a = {r[0]:.3f} to {r[-1]:.3f}")
    print(f"\nElectron density ne [10^20/m^3]:")
    print(f"  Center: {ne[0]:.3f}, Edge: {ne[-1]:.3f}, Average: {np.mean(ne):.3f}")
    print(f"\nElectron temperature Te [keV]:")
    print(f"  Center: {Te[0]:.3f}, Edge: {Te[-1]:.3f}, Average: {np.mean(Te):.3f}")
    print(f"\nIon temperature TD [keV]:")
    print(f"  Center: {TD[0]:.3f}, Edge: {TD[-1]:.3f}, Average: {np.mean(TD):.3f}")

    # Transport coefficients
    if 'chi_e' in data:
        chi_e = data['chi_e']['values'].get('AKE', np.zeros_like(r))
        chi_e_nc = data['chi_e']['values'].get('AKNCE', np.zeros_like(r))
        chi_e_an = data['chi_e']['values'].get('AKDWE', np.zeros_like(r))
        print(f"\nElectron thermal diffusivity chi_e [m^2/s]:")
        print(f"  Total: {np.mean(chi_e):.3f} (avg), NC: {np.mean(chi_e_nc):.3f}, Anom: {np.mean(chi_e_an):.3f}")

    if 'chi_i' in data:
        chi_i = data['chi_i']['values'].get('AKD', np.zeros_like(r))
        chi_i_nc = data['chi_i']['values'].get('AKNCD', np.zeros_like(r))
        chi_i_an = data['chi_i']['values'].get('AKDWD', np.zeros_like(r))
        print(f"\nIon thermal diffusivity chi_i [m^2/s]:")
        print(f"  Total: {np.mean(chi_i):.3f} (avg), NC: {np.mean(chi_i_nc):.3f}, Anom: {np.mean(chi_i_an):.3f}")

    # Power balance
    if 'power_in' in data:
        vals = data['power_in']['values']
        POH = vals.get('POH', np.zeros_like(r))
        PNB = vals.get('PNB', np.zeros_like(r))
        PNF = vals.get('PNF', np.zeros_like(r))
        PRF = vals.get('PRF', np.zeros_like(r))

        print(f"\n" + "="*60)
        print("Power Density [MW/m^3]")
        print("="*60)
        print(f"Ohmic (POH):    center={POH[0]:.4f}, avg={np.mean(POH):.4f}")
        print(f"NBI (PNB):      center={PNB[0]:.4f}, avg={np.mean(PNB):.4f}")
        print(f"Fusion (PNF):   center={PNF[0]:.4f}, avg={np.mean(PNF):.4f}")
        print(f"RF (PRF):       center={PRF[0]:.4f}, avg={np.mean(PRF):.4f}")

    if 'power_loss' in data:
        vals = data['power_loss']['values']
        PRSUM = vals.get('PRSUM', np.zeros_like(r))
        PRB = vals.get('PRB', np.zeros_like(r))
        PRC = vals.get('PRC', np.zeros_like(r))
        PRL = vals.get('PRL', np.zeros_like(r))
        PCX = vals.get('PCX', np.zeros_like(r))
        PIE = vals.get('PIE', np.zeros_like(r))

        print(f"\nBremsstrahlung (PRB): center={PRB[0]:.4f}, avg={np.mean(PRB):.4f}")
        print(f"Synchrotron (PRC):    center={PRC[0]:.4f}, avg={np.mean(PRC):.4f}")
        print(f"Line rad (PRL):       center={PRL[0]:.4f}, avg={np.mean(PRL):.4f}")
        print(f"Charge exch (PCX):    center={PCX[0]:.4f}, avg={np.mean(PCX):.4f}")
        print(f"e-i exchange (PIE):   center={PIE[0]:.4f}, avg={np.mean(PIE):.4f}")

    # Estimate temperature from steady-state heat equation
    print(f"\n" + "="*60)
    print("Steady-State Temperature Estimate")
    print("="*60)

    T_calc = None
    if 'chi_e' in data and 'power_in' in data and 'power_loss' in data:
        # Total heating = Ohmic + RF + Fusion + NBI
        S_heating = POH + PRF + PNF + PNB
        # Total radiation loss
        S_loss = PRSUM
        T_calc = estimate_temperature_from_balance(r, ne, chi_e, S_heating, S_loss, Te, data)

    # Plot comparison
    plot_energy_balance(r, data, T_calc)

    # Compare with OMFIT if requested
    if compare_omfit:
        compare_with_omfit(r, data)

    return data

def solve_temperature_profile(r, n, chi, S_net, T_edge, a=A_MINOR):
    """
    Solve for temperature profile from steady-state heat equation.

    Steady-state 1D cylindrical heat equation:
    1/r * d/dr(r * n * chi * dT/dr) = -S_net

    where S_net is the net volumetric heating rate [MW/m³].

    Method:
    1. Integrate S_net * r from 0 to r to get cumulative heat
    2. Heat flux: q(r) = (1/r) * integral(S_net * r * dr)
    3. dT/dr = -q / (n * chi)
    4. Integrate from edge (T_edge) inward to get T(r)

    Args:
        r: normalized radius (r/a), must be increasing from 0 to 1
        n: density profile [10^20/m^3]
        chi: thermal diffusivity profile [m^2/s]
        S_net: net heating rate profile [MW/m³] (positive = heating)
        T_edge: boundary temperature at r/a=1 [keV]
        a: minor radius [m]

    Returns:
        T_calc: calculated temperature profile [keV]
    """
    # Convert r/a to real radius [m]
    r_real = r * a

    # Unit conversion for S_net:
    # S [MW/m³] = 1e6 W/m³ = 1e6 J/s/m³
    # n [10^20/m³] = 1e20 /m³
    # chi [m²/s]
    # T [keV], 1 keV = 1.602e-16 J
    #
    # Heat equation: 1/r * d/dr(r * n * chi * dT/dr) = -S / (1.5)
    # where 1.5 is the factor from (3/2)nT energy density
    #
    # Converting S [MW/m³] to [keV * 10^20/m³ / s]:
    # 1 MW/m³ = 1e6 J/s/m³ = 1e6 / 1.602e-16 keV/s/m³ = 6.24e21 keV/s/m³
    # Per 10^20 particles: 6.24e21 / 1e20 = 62.4 keV/s per (10^20/m³)

    S_conv = S_net * 62.4 / 1.5  # Convert to [keV * 10^20/m³ / s] / 1.5

    # Step 1: Compute cumulative integral of S * r * dr
    # This gives the total heat flowing through radius r
    integrand = S_conv * r_real
    cumulative_heat = np.zeros_like(r)
    for i in range(1, len(r)):
        cumulative_heat[i] = np.trapz(integrand[:i+1], r_real[:i+1])

    # Step 2: Heat flux q(r) = cumulative_heat / r
    # At r=0, use L'Hopital's rule: q(0) = 0
    q = np.zeros_like(r)
    q[1:] = cumulative_heat[1:] / r_real[1:]

    # Step 3: dT/dr = -q / (n * chi)
    # Avoid division by zero
    n_safe = np.maximum(n, 1e-10)
    chi_safe = np.maximum(chi, 1e-10)

    dTdr = -q / (n_safe * chi_safe)

    # Step 4: Integrate from edge inward to get T(r)
    # T(r) = T_edge + integral from r to 1 of dT/dr' * dr'
    # which equals T_edge - integral from r to 1 of q/(n*chi) dr'
    T_calc = np.zeros_like(r)
    T_calc[-1] = T_edge

    # Integrate backwards from edge to center
    for i in range(len(r) - 2, -1, -1):
        # T(i) = T(i+1) - dT/dr * dr
        dr = r_real[i + 1] - r_real[i]
        T_calc[i] = T_calc[i + 1] - dTdr[i + 1] * dr

    return T_calc


def estimate_temperature_from_balance(r, n, chi, S_heating, S_loss, T_actual, data=None):
    """
    Estimate temperature from heat equation in steady state.

    Solves the heat equation with net source term and compares to actual.

    Args:
        r: normalized radius (r/a)
        n: density profile [10^20/m³]
        chi: thermal diffusivity profile [m²/s]
        S_heating: total heating power density [MW/m³]
        S_loss: total loss power density [MW/m³] (radiation, etc.)
        T_actual: actual temperature profile from simulation [keV]
        data: full data dictionary (optional, for accessing additional fields)

    Returns:
        T_calc: calculated temperature profile [keV]
    """
    # Net source = heating - losses
    S_net = S_heating - S_loss

    # Get boundary temperature from actual profile
    T_edge = T_actual[-1]

    print(f"\nSolving steady-state heat equation:")
    print(f"  S_heating (center) = {S_heating[0]:.4f} MW/m³")
    print(f"  S_loss (center)    = {S_loss[0]:.4f} MW/m³")
    print(f"  S_net (center)     = {S_net[0]:.4f} MW/m³")
    print(f"  T_edge (BC)        = {T_edge:.3f} keV")

    # Solve for temperature profile
    T_calc = solve_temperature_profile(r, n, chi, S_net, T_edge)

    # Print comparison
    print(f"\n  Temperature comparison:")
    print(f"    {'Location':<15} {'Actual [keV]':>15} {'Calculated [keV]':>18} {'Ratio':>10}")
    print(f"    {'-'*60}")
    print(f"    {'Center (r=0)':<15} {T_actual[0]:>15.2f} {T_calc[0]:>18.2f} {T_calc[0]/T_actual[0] if T_actual[0] > 0 else 0:>10.2f}")

    # Find r/a = 0.5
    idx_half = np.argmin(np.abs(r - 0.5))
    print(f"    {'Mid (r/a=0.5)':<15} {T_actual[idx_half]:>15.2f} {T_calc[idx_half]:>18.2f} {T_calc[idx_half]/T_actual[idx_half] if T_actual[idx_half] > 0 else 0:>10.2f}")
    print(f"    {'Edge (r/a=1)':<15} {T_actual[-1]:>15.2f} {T_calc[-1]:>18.2f} {'1.00':>10}")

    # Calculate integrated power
    P_heating = integrate_power(r, S_heating)
    P_loss = integrate_power(r, S_loss)
    P_net = integrate_power(r, S_net)

    print(f"\n  Integrated powers:")
    print(f"    Heating: {P_heating:.1f} MW")
    print(f"    Loss:    {P_loss:.1f} MW")
    print(f"    Net:     {P_net:.1f} MW")

    # If calculated T is much lower than actual, suggests missing heating
    # If calculated T is much higher than actual, suggests missing losses
    if T_calc[0] < T_actual[0] * 0.5:
        print(f"\n  ! Calculated T too LOW -> Missing heating sources or chi too high")
    elif T_calc[0] > T_actual[0] * 1.5:
        print(f"\n  ! Calculated T too HIGH -> Missing losses (radiation?) or chi too low")
    else:
        print(f"\n  ✓ Temperature profile reasonably matched")

    return T_calc

def plot_energy_balance(r, data, T_calc=None):
    """Plot energy balance profiles.

    Args:
        r: radial grid (r/a)
        data: dictionary containing all profile data
        T_calc: calculated temperature profile from heat equation (optional)
    """

    fig, axes = plt.subplots(2, 3, figsize=(14, 8))

    # 1. Density profiles
    ax = axes[0, 0]
    if 'density' in data:
        for key, vals in data['density']['values'].items():
            ax.plot(r, vals, label=key)
    ax.set_xlabel('r/a')
    ax.set_ylabel('n [10^20/m^3]')
    ax.set_title('Density Profiles')
    ax.legend()
    ax.grid(True)

    # 2. Temperature profiles (with calculated T overlay)
    ax = axes[0, 1]
    if 'temperature' in data:
        for key, vals in data['temperature']['values'].items():
            ax.plot(r, vals, '-', linewidth=2, label=f'{key} (sim)')
    # Overlay calculated temperature
    if T_calc is not None:
        ax.plot(r, T_calc, 'k--', linewidth=2, label='Te (calculated)')
    ax.set_xlabel('r/a')
    ax.set_ylabel('T [keV]')
    ax.set_title('Temperature Profiles')
    ax.legend()
    ax.grid(True)

    # 3. Transport coefficients
    ax = axes[0, 2]
    if 'chi_e' in data:
        for key, vals in data['chi_e']['values'].items():
            ax.plot(r, vals, label=f'e: {key}')
    if 'chi_i' in data:
        for key, vals in data['chi_i']['values'].items():
            ax.plot(r, vals, '--', label=f'i: {key}')
    ax.set_xlabel('r/a')
    ax.set_ylabel('chi [m^2/s]')
    ax.set_title('Transport Coefficients')
    ax.legend()
    ax.grid(True)

    # 4. Power input
    ax = axes[1, 0]
    if 'power_in' in data:
        for key, vals in data['power_in']['values'].items():
            if not key.startswith('Y'):
                ax.plot(r, vals, label=key)
    ax.set_xlabel('r/a')
    ax.set_ylabel('P [MW/m^3]')
    ax.set_title('Power Input')
    ax.legend()
    ax.grid(True)

    # 5. Power loss
    ax = axes[1, 1]
    if 'power_loss' in data:
        for key, vals in data['power_loss']['values'].items():
            ax.plot(r, vals, label=key)
    ax.set_xlabel('r/a')
    ax.set_ylabel('P [MW/m^3]')
    ax.set_title('Power Loss')
    ax.legend()
    ax.grid(True)

    # 6. Net power balance
    ax = axes[1, 2]
    if 'power_in' in data and 'power_loss' in data:
        P_in = data['power_in']['values']
        P_out = data['power_loss']['values']

        total_in = sum(v for k, v in P_in.items() if not k.startswith('Y') and k != 'PRSUM')
        total_out = P_out.get('PRSUM', np.zeros_like(r)) + P_out.get('PCX', np.zeros_like(r))

        ax.plot(r, total_in, 'b-', label='P_in (POH+PNB+PNF+PRF)')
        ax.plot(r, total_out, 'r-', label='P_out (rad+CX)')
        ax.plot(r, total_in - total_out, 'g--', label='Net (P_in - P_out)')
        ax.axhline(y=0, color='k', linestyle=':')
    ax.set_xlabel('r/a')
    ax.set_ylabel('P [MW/m^3]')
    ax.set_title('Power Balance')
    ax.legend()
    ax.grid(True)

    plt.tight_layout()
    plt.savefig('energy_balance_analysis.png', dpi=150)
    print(f"\nPlot saved to energy_balance_analysis.png")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Analyze TR energy balance')
    parser.add_argument('-t', '--time-index', type=int, default=None,
                        help='Time index for evolution files (default: final)')
    parser.add_argument('--no-omfit', action='store_true',
                        help='Skip OMFIT comparison')
    parser.add_argument('--omfit-file', type=str, default='omfit_profiles_for_tr.csv',
                        help='Path to OMFIT profiles CSV file')

    args = parser.parse_args()

    analyze_energy_balance(time_index=args.time_index, compare_omfit=not args.no_omfit)
