#!/usr/bin/env python3
"""
Calculate temperature profiles from steady-state heat equation and compare with TR simulation.

This script solves the 1D cylindrical heat equation:
    1/r * d/dr(r * n * chi * dT/dr) = -S_net

For both electrons and ions, using the source terms and transport coefficients
from TR simulation output.

Usage:
    python3 calculate_temperature_profile.py
    python3 calculate_temperature_profile.py --save-csv  # Also save to CSV
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import interp1d

# Physical parameters for CFEDR
R0 = 7.8      # Major radius [m]
A_MINOR = 2.44  # Minor radius [m]


def read_tr_csv(filename):
    """Read TR output CSV file, return r/a and data columns."""
    with open(filename, 'r') as f:
        title = f.readline().strip()

    df = pd.read_csv(filename, skiprows=1, index_col=False)
    df = df.loc[:, ~df.columns.str.contains('^Unnamed')]

    r = df.iloc[:, 0].values
    data = {col: df[col].values for col in df.columns[1:]}

    return r, data, title


def integrate_power(r, p, R0=R0, a=A_MINOR):
    """Integrate power density over plasma volume."""
    r_real = r * a
    dV = 2 * np.pi * R0 * 2 * np.pi * r_real
    return np.trapz(p * dV, r_real)


def solve_temperature_profile(r, n, chi, S_net, T_edge, a=A_MINOR, chi_factor=0.5):
    """
    Solve for temperature profile from steady-state heat equation.

    Steady-state 1D cylindrical heat equation:
    1/r * d/dr(r * n * chi * dT/dr) = -S_net

    Args:
        r: normalized radius (r/a)
        n: density profile [10^20/m^3]
        chi: thermal diffusivity profile [m^2/s]
        S_net: net heating rate profile [MW/m³]
        T_edge: boundary temperature at r/a=1 [keV]
        a: minor radius [m]
        chi_factor: correction factor for chi (default 0.5)

            The chi values from TR output need to be scaled by ~0.5 to match
            the steady-state calculation. This is because TR uses volume-weighted
            coordinates with geometric factors:
            - Source term: multiplied by DVRHO^(5/3)
            - Diffusion term: multiplied by DVRHO^(2/3) * FB
            - FB = DVRHO * AR2RHO / DR^2 where AR2RHO = 1/(κ*a)²

            This volume-weighted discretization effectively reduces the
            apparent thermal diffusivity by approximately factor of 2 compared
            to simple cylindrical geometry.

    Returns:
        T_calc: calculated temperature profile [keV]
        q: heat flux profile
        dTdr: temperature gradient profile
    """
    r_real = r * a

    # Apply chi correction factor
    chi_eff = chi * chi_factor

    # Unit conversion:
    # S [MW/m³] -> [keV * 10^20/m³ / s]
    # 1 MW/m³ = 1e6 J/s/m³ = 1e6 / 1.602e-16 keV/s/m³ = 6.24e21 keV/s/m³
    # Per 10^20: 6.24e21 / 1e20 = 62.4 keV/s per (10^20/m³)
    # Divide by 1.5 for (3/2)nT energy density factor
    S_conv = S_net * 62.4 / 1.5

    # Step 1: Cumulative integral of S * r * dr
    integrand = S_conv * r_real
    cumulative_heat = np.zeros_like(r)
    for i in range(1, len(r)):
        cumulative_heat[i] = np.trapz(integrand[:i+1], r_real[:i+1])

    # Step 2: Heat flux q(r) = cumulative_heat / r
    q = np.zeros_like(r)
    q[1:] = cumulative_heat[1:] / r_real[1:]

    # Step 3: dT/dr = -q / (n * chi_eff)
    n_safe = np.maximum(n, 1e-10)
    chi_safe = np.maximum(chi_eff, 1e-10)
    dTdr = -q / (n_safe * chi_safe)

    # Step 4: Integrate from edge inward
    T_calc = np.zeros_like(r)
    T_calc[-1] = T_edge

    for i in range(len(r) - 2, -1, -1):
        dr = r_real[i + 1] - r_real[i]
        T_calc[i] = T_calc[i + 1] - dTdr[i + 1] * dr

    return T_calc, q, dTdr


def load_tr_data():
    """Load all necessary TR output data."""
    data = {}

    # Density profiles
    r, dens, _ = read_tr_csv('tr_data_017.csv')
    data['r'] = r
    data['ne'] = dens.get('nE', np.zeros_like(r))
    data['nD'] = dens.get('nD', np.zeros_like(r))
    data['nT'] = dens.get('nT', np.zeros_like(r))
    data['nA'] = dens.get('nA', np.zeros_like(r))
    # Total ion density
    data['ni'] = data['nD'] + data['nT'] + data['nA']

    # Temperature profiles
    _, temp, _ = read_tr_csv('tr_data_019.csv')
    data['Te'] = temp.get('TE', np.zeros_like(r))
    data['TD'] = temp.get('TD', np.zeros_like(r))
    data['TT'] = temp.get('TT', np.zeros_like(r))
    data['TA'] = temp.get('TA', np.zeros_like(r))
    # Use D temperature as representative ion temperature
    data['Ti'] = data['TD']

    # Power input profiles
    _, power_in, _ = read_tr_csv('tr_data_021.csv')
    data['POH'] = power_in.get('POH', np.zeros_like(r))
    data['PNB'] = power_in.get('PNB', np.zeros_like(r))
    data['PNF'] = power_in.get('PNF', np.zeros_like(r))
    data['PRF'] = power_in.get('PRF', np.zeros_like(r))

    # Power loss profiles
    _, power_loss, _ = read_tr_csv('tr_data_022.csv')
    data['PRSUM'] = power_loss.get('PRSUM', np.zeros_like(r))
    data['PRB'] = power_loss.get('PRB', np.zeros_like(r))
    data['PRC'] = power_loss.get('PRC', np.zeros_like(r))
    data['PRL'] = power_loss.get('PRL', np.zeros_like(r))
    data['PCX'] = power_loss.get('PCX', np.zeros_like(r))
    data['PIE'] = power_loss.get('PIE', np.zeros_like(r))

    # Transport coefficients
    _, chi_e, _ = read_tr_csv('tr_data_023.csv')
    data['chi_e'] = chi_e.get('AKE', np.zeros_like(r))

    _, chi_i, _ = read_tr_csv('tr_data_024.csv')
    data['chi_i'] = chi_i.get('AKD', np.zeros_like(r))

    # Species-specific power input (PIN_NS)
    if os.path.exists('tr_data_037.csv'):
        _, pin_ns, _ = read_tr_csv('tr_data_037.csv')
        data['PIN_e'] = pin_ns.get('PIN_1', np.zeros_like(r))    # electron
        data['PIN_D'] = pin_ns.get('PIN_2', np.zeros_like(r))    # D
        data['PIN_T'] = pin_ns.get('PIN_3', np.zeros_like(r))    # T
        data['PIN_He4'] = pin_ns.get('PIN_4', np.zeros_like(r))  # He4 (alpha)

    return data


def calculate_temperatures(data, use_pin_ns=True, chi_factor=0.5):
    """
    Calculate electron and ion temperature profiles from heat equation.

    Args:
        data: dictionary with TR simulation data
        use_pin_ns: if True, use species-specific PIN data from tr_data_037.csv
                    if False, use assumed 20/80 split of PNF
        chi_factor: correction factor for chi (default 0.5)

            TR uses volume-weighted coordinates that result in effective
            chi values about 2x larger than what's needed in simple
            cylindrical geometry. Use chi_factor=0.5 to match TR results.
            Use chi_factor=1.0 to use chi values as-is.

    Electron heat balance (using PIN values from TR):
        S_e_net = PIN_e  (NOTE: PIN already has losses subtracted in TR)

    Ion heat balance:
        S_i_net = PIN_D + PIN_T + PIN_He4 (already net values)
    """
    r = data['r']

    print("="*70)
    print("TEMPERATURE CALCULATION FROM HEAT EQUATION")
    print("="*70)

    if use_pin_ns and 'PIN_e' in data:
        # Use actual species-specific power from TR
        # IMPORTANT: PIN values from TR already include loss subtraction!
        # PIN_e = heating - PRSUM - PIE (Fortran: trcalc.f90)
        # PIN_D = heating - PCX*fraction
        # PIN_T = heating - PCX*fraction
        # PIN_He4 = heating only
        # So PIN is already the NET source term for the temperature equation

        S_e_net = data['PIN_e']  # Already net (losses subtracted in Fortran)
        S_i_net = data['PIN_D'] + data['PIN_T'] + data['PIN_He4']  # Already net

        print("\n--- Using actual species-specific power (PIN_NS) ---")
        print("  NOTE: PIN values already include loss subtraction in TR Fortran code")
        print(f"  PIN_e (electron net): {data['PIN_e'][0]:.4f} MW/m³")
        print(f"  PIN_D (deuterium net): {data['PIN_D'][0]:.4f} MW/m³")
        print(f"  PIN_T (tritium net):   {data['PIN_T'][0]:.4f} MW/m³")
        print(f"  PIN_He4 (alpha):       {data['PIN_He4'][0]:.4f} MW/m³")

        # Integrated powers for fraction calculation
        P_e = integrate_power(r, S_e_net)
        P_i = integrate_power(r, S_i_net)
        P_total = P_e + P_i
        print(f"\n  Power split: electrons {P_e/P_total*100:.1f}%, ions {P_i/P_total*100:.1f}%")

        # For display purposes
        S_e_heating = data['PIN_e']
        S_i_heating = S_i_net
        S_e_loss = np.zeros_like(r)  # Losses already in PIN
        S_i_loss = np.zeros_like(r)

    else:
        # Use assumed split (20% electrons, 80% ions)
        alpha_e_frac = 0.2
        alpha_i_frac = 0.8

        S_e_heating = (data['POH'] +
                       data['PRF'] +
                       alpha_e_frac * data['PNF'])

        S_i_heating = (data['PNB'] +
                       alpha_i_frac * data['PNF'] +
                       data['PIE'])

        print("\n--- Using assumed power split (20% e, 80% i) ---")
        print(f"  PNF_e (20%): {alpha_e_frac * data['PNF'][0]:.4f} MW/m³")
        print(f"  PNF_i (80%): {alpha_i_frac * data['PNF'][0]:.4f} MW/m³")

        # Electron losses
        S_e_loss = data['PRSUM'] + data['PIE']  # Radiation + e-i exchange

        # Net electron source
        S_e_net = S_e_heating - S_e_loss

        # Ion losses
        S_i_loss = data['PCX']  # Charge exchange

        # Net ion source
        S_i_net = S_i_heating - S_i_loss

    # Electron temperature
    print("\n--- ELECTRON TEMPERATURE ---")
    print(f"  chi_factor = {chi_factor} (TR uses volume-weighted coords)")
    print(f"  S_e_heating (center): {S_e_heating[0]:.4f} MW/m³")
    print(f"  S_e_loss (center):    {S_e_loss[0]:.4f} MW/m³")
    print(f"    PRSUM: {data['PRSUM'][0]:.4f}, PIE: {data['PIE'][0]:.4f}")
    print(f"  S_e_net (center):     {S_e_net[0]:.4f} MW/m³")

    Te_edge = data['Te'][-1]
    Te_calc, q_e, dTedr = solve_temperature_profile(
        r, data['ne'], data['chi_e'], S_e_net, Te_edge, chi_factor=chi_factor
    )

    # Ion temperature
    print("\n--- ION TEMPERATURE ---")
    print(f"  chi_factor = {chi_factor}")
    print(f"  S_i_heating (center): {S_i_heating[0]:.4f} MW/m³")
    print(f"  S_i_loss (center):    {S_i_loss[0]:.4f} MW/m³")
    print(f"    PCX: {data['PCX'][0]:.4f}")
    print(f"  S_i_net (center):     {S_i_net[0]:.4f} MW/m³")

    Ti_edge = data['Ti'][-1]
    Ti_calc, q_i, dTidr = solve_temperature_profile(
        r, data['ni'], data['chi_i'], S_i_net, Ti_edge, chi_factor=chi_factor
    )

    # Store results
    results = {
        'Te_calc': Te_calc,
        'Ti_calc': Ti_calc,
        'q_e': q_e,
        'q_i': q_i,
        'S_e_net': S_e_net,
        'S_i_net': S_i_net,
        'S_e_heating': S_e_heating,
        'S_i_heating': S_i_heating,
        'S_e_loss': S_e_loss,
        'S_i_loss': S_i_loss,
    }

    return results


def print_comparison(data, results):
    """Print temperature comparison table."""
    r = data['r']

    print("\n" + "="*70)
    print("TEMPERATURE COMPARISON: Calculated vs Simulation")
    print("="*70)

    # Find indices for key positions
    idx_center = 0
    idx_half = np.argmin(np.abs(r - 0.5))
    idx_edge = -1

    print(f"\n{'Location':<15} {'Te_sim':>10} {'Te_calc':>10} {'Ratio':>8} | "
          f"{'Ti_sim':>10} {'Ti_calc':>10} {'Ratio':>8}")
    print("-" * 85)

    locations = [
        ('Center (r=0)', idx_center),
        ('Mid (r/a=0.5)', idx_half),
        ('Edge (r/a=1)', idx_edge),
    ]

    for name, idx in locations:
        Te_sim = data['Te'][idx]
        Te_calc = results['Te_calc'][idx]
        Te_ratio = Te_calc / Te_sim if Te_sim > 0 else 0

        Ti_sim = data['Ti'][idx]
        Ti_calc = results['Ti_calc'][idx]
        Ti_ratio = Ti_calc / Ti_sim if Ti_sim > 0 else 0

        print(f"{name:<15} {Te_sim:>10.2f} {Te_calc:>10.2f} {Te_ratio:>8.2f} | "
              f"{Ti_sim:>10.2f} {Ti_calc:>10.2f} {Ti_ratio:>8.2f}")

    # Integrated powers
    print("\n--- INTEGRATED POWERS [MW] ---")
    print(f"{'Source':<20} {'Electron':>12} {'Ion':>12}")
    print("-" * 45)

    P_e_heat = integrate_power(r, results['S_e_heating'])
    P_i_heat = integrate_power(r, results['S_i_heating'])
    P_e_loss = integrate_power(r, results['S_e_loss'])
    P_i_loss = integrate_power(r, results['S_i_loss'])
    P_e_net = integrate_power(r, results['S_e_net'])
    P_i_net = integrate_power(r, results['S_i_net'])

    print(f"{'Heating':<20} {P_e_heat:>12.1f} {P_i_heat:>12.1f}")
    print(f"{'Losses':<20} {P_e_loss:>12.1f} {P_i_loss:>12.1f}")
    print(f"{'Net (conducted out)':<20} {P_e_net:>12.1f} {P_i_net:>12.1f}")


def plot_temperature_comparison(data, results, output_file='temperature_comparison.png'):
    """Plot comprehensive temperature comparison."""
    r = data['r']

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))

    # 1. Electron temperature comparison
    ax = axes[0, 0]
    ax.plot(r, data['Te'], 'b-', linewidth=2, label='Te (simulation)')
    ax.plot(r, results['Te_calc'], 'b--', linewidth=2, label='Te (calculated)')
    ax.set_xlabel('r/a')
    ax.set_ylabel('Te [keV]')
    ax.set_title('Electron Temperature')
    ax.legend()
    ax.grid(True)
    ax.set_xlim([0, 1])

    # 2. Ion temperature comparison
    ax = axes[0, 1]
    ax.plot(r, data['Ti'], 'r-', linewidth=2, label='Ti (simulation)')
    ax.plot(r, results['Ti_calc'], 'r--', linewidth=2, label='Ti (calculated)')
    ax.set_xlabel('r/a')
    ax.set_ylabel('Ti [keV]')
    ax.set_title('Ion Temperature')
    ax.legend()
    ax.grid(True)
    ax.set_xlim([0, 1])

    # 3. Both temperatures overlaid
    ax = axes[0, 2]
    ax.plot(r, data['Te'], 'b-', linewidth=2, label='Te (sim)')
    ax.plot(r, results['Te_calc'], 'b--', linewidth=2, label='Te (calc)')
    ax.plot(r, data['Ti'], 'r-', linewidth=2, label='Ti (sim)')
    ax.plot(r, results['Ti_calc'], 'r--', linewidth=2, label='Ti (calc)')
    ax.set_xlabel('r/a')
    ax.set_ylabel('T [keV]')
    ax.set_title('Temperature Comparison')
    ax.legend()
    ax.grid(True)
    ax.set_xlim([0, 1])

    # 4. Electron thermal diffusivity (log scale)
    ax = axes[1, 0]
    ax.semilogy(r, data['chi_e'], 'b-', linewidth=2, label='χ_e')
    ax.semilogy(r, data['chi_i'], 'r-', linewidth=2, label='χ_i')
    ax.set_xlabel('r/a')
    ax.set_ylabel('χ [m²/s]')
    ax.set_title('Thermal Diffusivity (log scale)')
    ax.legend()
    ax.grid(True)
    ax.set_xlim([0, 1])

    # 5. Net heating source
    ax = axes[1, 1]
    ax.plot(r, results['S_e_net'], 'b-', linewidth=2, label='S_e_net')
    ax.plot(r, results['S_i_net'], 'r-', linewidth=2, label='S_i_net')
    ax.axhline(y=0, color='k', linestyle='--', linewidth=0.5)
    ax.set_xlabel('r/a')
    ax.set_ylabel('S [MW/m³]')
    ax.set_title('Net Heating Source')
    ax.legend()
    ax.grid(True)
    ax.set_xlim([0, 1])

    # 6. Temperature ratio (calculated/simulation)
    ax = axes[1, 2]
    Te_ratio = np.where(data['Te'] > 0.1, results['Te_calc'] / data['Te'], np.nan)
    Ti_ratio = np.where(data['Ti'] > 0.1, results['Ti_calc'] / data['Ti'], np.nan)
    ax.plot(r, Te_ratio, 'b-', linewidth=2, label='Te_calc/Te_sim')
    ax.plot(r, Ti_ratio, 'r-', linewidth=2, label='Ti_calc/Ti_sim')
    ax.axhline(y=1.0, color='k', linestyle='--', linewidth=1, label='Perfect match')
    ax.set_xlabel('r/a')
    ax.set_ylabel('Ratio')
    ax.set_title('Temperature Ratio (Calc/Sim)')
    ax.legend()
    ax.grid(True)
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 2])

    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    print(f"\nPlot saved to {output_file}")


def plot_heat_balance_detail(data, results, output_file='heat_balance_detail.png'):
    """Plot detailed heat balance for electrons and ions."""
    r = data['r']

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # 1. Electron heating sources
    ax = axes[0, 0]
    if 'PIN_e' in data:
        ax.plot(r, data['PIN_e'], 'b-', linewidth=2, label='PIN_e (from TR)')
    else:
        ax.plot(r, data['POH'], 'g-', linewidth=2, label='POH (Ohmic)')
        ax.plot(r, data['PRF'], 'c-', linewidth=2, label='PRF (RF)')
        ax.plot(r, 0.2 * data['PNF'], 'm-', linewidth=2, label='0.2*PNF (Fusion→e)')
    ax.plot(r, results['S_e_heating'], 'b--', linewidth=2, label='Total heating')
    ax.set_xlabel('r/a')
    ax.set_ylabel('P [MW/m³]')
    ax.set_title('Electron Heating Sources')
    ax.legend()
    ax.grid(True)
    ax.set_xlim([0, 1])

    # 2. Electron losses
    ax = axes[0, 1]
    ax.plot(r, data['PRSUM'], 'r-', linewidth=2, label='PRSUM (Radiation)')
    ax.plot(r, data['PIE'], 'orange', linewidth=2, label='PIE (e→i exchange)')
    ax.plot(r, results['S_e_loss'], 'r--', linewidth=2, label='Total loss')
    ax.set_xlabel('r/a')
    ax.set_ylabel('P [MW/m³]')
    ax.set_title('Electron Losses')
    ax.legend()
    ax.grid(True)
    ax.set_xlim([0, 1])

    # 3. Ion heating sources
    ax = axes[1, 0]
    if 'PIN_D' in data:
        ax.plot(r, data['PIN_D'], 'g-', linewidth=2, label='PIN_D')
        ax.plot(r, data['PIN_T'], 'm-', linewidth=2, label='PIN_T')
        ax.plot(r, data['PIN_He4'], 'c-', linewidth=2, label='PIN_He4')
    else:
        ax.plot(r, data['PNB'], 'g-', linewidth=2, label='PNB (NBI)')
        ax.plot(r, 0.8 * data['PNF'], 'm-', linewidth=2, label='0.8*PNF (Fusion→i)')
    ax.plot(r, data['PIE'], 'orange', linewidth=2, label='PIE (e→i exchange)')
    ax.plot(r, results['S_i_heating'], 'r--', linewidth=2, label='Total heating')
    ax.set_xlabel('r/a')
    ax.set_ylabel('P [MW/m³]')
    ax.set_title('Ion Heating Sources')
    ax.legend()
    ax.grid(True)
    ax.set_xlim([0, 1])

    # 4. Ion losses
    ax = axes[1, 1]
    ax.plot(r, data['PCX'], 'purple', linewidth=2, label='PCX (Charge exchange)')
    ax.plot(r, results['S_i_loss'], 'r--', linewidth=2, label='Total loss')
    ax.set_xlabel('r/a')
    ax.set_ylabel('P [MW/m³]')
    ax.set_title('Ion Losses')
    ax.legend()
    ax.grid(True)
    ax.set_xlim([0, 1])

    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    print(f"Heat balance detail plot saved to {output_file}")


def save_to_csv(data, results, output_file='temperature_calculation.csv'):
    """Save calculation results to CSV."""
    df = pd.DataFrame({
        'r/a': data['r'],
        'Te_sim': data['Te'],
        'Te_calc': results['Te_calc'],
        'Ti_sim': data['Ti'],
        'Ti_calc': results['Ti_calc'],
        'chi_e': data['chi_e'],
        'chi_i': data['chi_i'],
        'S_e_net': results['S_e_net'],
        'S_i_net': results['S_i_net'],
        'ne': data['ne'],
        'ni': data['ni'],
    })

    df.to_csv(output_file, index=False, float_format='%.6e')
    print(f"Results saved to {output_file}")


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description='Calculate temperature profiles from heat equation')
    parser.add_argument('--save-csv', action='store_true',
                        help='Save results to CSV file')
    parser.add_argument('--output', type=str, default='temperature_comparison.png',
                        help='Output plot filename')

    args = parser.parse_args()

    # Check if required files exist
    required_files = [
        'tr_data_017.csv', 'tr_data_019.csv', 'tr_data_021.csv',
        'tr_data_022.csv', 'tr_data_023.csv', 'tr_data_024.csv'
    ]

    for f in required_files:
        if not os.path.exists(f):
            print(f"Error: Required file {f} not found")
            return

    # Load data
    print("Loading TR simulation data...")
    data = load_tr_data()

    # Calculate temperatures
    results = calculate_temperatures(data)

    # Print comparison
    print_comparison(data, results)

    # Plot results
    plot_temperature_comparison(data, results, args.output)
    plot_heat_balance_detail(data, results)

    # Save to CSV if requested
    if args.save_csv:
        save_to_csv(data, results)

    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)


if __name__ == '__main__':
    main()
