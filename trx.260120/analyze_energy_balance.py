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

def read_tr_csv(filename):
    """Read TR output CSV file, return x (r/a) and data columns."""
    with open(filename, 'r') as f:
        title = f.readline().strip()

    df = pd.read_csv(filename, skiprows=1, index_col=False)
    df = df.loc[:, ~df.columns.str.contains('^Unnamed')]

    x = df.iloc[:, 0].values
    data = {col: df[col].values for col in df.columns[1:]}

    return x, data, title


def read_evolution_csv(filename, time_index=-1):
    """
    Read time-space evolution CSV file at a specific time index.

    Args:
        filename: path to evolution CSV file (e.g., tr_data_047.csv)
        time_index: which time point to read (0=first, -1=last, or specific index)

    Returns:
        x: radial positions (r/a)
        values: array of values at the selected time
        time_indices: list of available time indices
        title: plot title
    """
    with open(filename, 'r') as f:
        title = f.readline().strip()

    df = pd.read_csv(filename, skiprows=1, index_col=False)
    df = df.loc[:, ~df.columns.str.contains('^Unnamed')]

    x = df.iloc[:, 0].values
    y_cols = list(df.columns[1:])
    n_times = len(y_cols)

    # Handle negative index
    if time_index < 0:
        time_index = n_times + time_index

    time_index = max(0, min(time_index, n_times - 1))

    selected_col = y_cols[time_index]
    values = df[selected_col].values

    return x, values, list(range(n_times)), title, time_index


def get_time_from_csv(time_csv='tr_data_001.csv'):
    """
    Get time values from a time-evolution CSV file.
    Returns array of time values corresponding to each time index.
    """
    if not os.path.exists(time_csv):
        return None

    df = pd.read_csv(time_csv, skiprows=1, index_col=False)
    df = df.loc[:, ~df.columns.str.contains('^Unnamed')]

    # First column is time (X)
    return df.iloc[:, 0].values


def read_evolution_data(time_index):
    """
    Read data from evolution files at a specific time index.

    Evolution files mapping:
    - tr_data_045: NE (electron density)
    - tr_data_046: ND (D density)
    - tr_data_047: TE (electron temperature)
    - tr_data_048: TD (D temperature)
    - tr_data_049: NT (T density)
    - tr_data_050: NA (alpha density)
    - tr_data_051: TT (T temperature)
    - tr_data_052: TA (alpha temperature)
    - tr_data_053: QP (safety factor)
    - tr_data_059: PIN (power input)
    - tr_data_060: POH (ohmic power)
    - tr_data_064: AKD (ion thermal diffusivity)
    """

    evolution_files = {
        'NE': 'tr_data_045.csv',
        'ND': 'tr_data_046.csv',
        'TE': 'tr_data_047.csv',
        'TD': 'tr_data_048.csv',
        'NT': 'tr_data_049.csv',
        'NA': 'tr_data_050.csv',
        'TT': 'tr_data_051.csv',
        'TA': 'tr_data_052.csv',
        'QP': 'tr_data_053.csv',
        'AJ': 'tr_data_054.csv',
        'PIN': 'tr_data_059.csv',
        'POH': 'tr_data_060.csv',
        'AKD': 'tr_data_064.csv',
    }

    data = {}
    r = None

    for var_name, filename in evolution_files.items():
        if os.path.exists(filename):
            x, values, time_indices, title, actual_idx = read_evolution_csv(filename, time_index)
            if r is None:
                r = x
            print(f"Read {filename} at time index {actual_idx}: {var_name}")

            # Group into density and temperature categories
            if var_name in ['NE', 'ND', 'NT', 'NA']:
                if 'density' not in data:
                    data['density'] = {'r': r, 'values': {}, 'title': 'Density'}
                # Map to standard names
                name_map = {'NE': 'nE', 'ND': 'nD', 'NT': 'nT', 'NA': 'nA'}
                data['density']['values'][name_map[var_name]] = values

            elif var_name in ['TE', 'TD', 'TT', 'TA']:
                if 'temperature' not in data:
                    data['temperature'] = {'r': r, 'values': {}, 'title': 'Temperature'}
                data['temperature']['values'][var_name] = values

            elif var_name == 'AKD':
                if 'chi_i' not in data:
                    data['chi_i'] = {'r': r, 'values': {}, 'title': 'Transport'}
                data['chi_i']['values']['AKD'] = values

            elif var_name == 'PIN':
                if 'power_in' not in data:
                    data['power_in'] = {'r': r, 'values': {}, 'title': 'Power Input'}
                data['power_in']['values']['PIN'] = values

            elif var_name == 'POH':
                if 'power_in' not in data:
                    data['power_in'] = {'r': r, 'values': {}, 'title': 'Power Input'}
                data['power_in']['values']['POH'] = values

    return data

def analyze_energy_balance(time_index=None):
    """
    Analyze energy balance at a specific time point.

    Args:
        time_index: which time point to analyze
                    None or -1: final time (use radial profile files 017-044)
                    0: first time point
                    1,2,3,...: specific time index (uses evolution files 045-064)

    For radial profiles (files 017-044), data is at final time.
    For time-space evolution (files 045+), can select time index.
    """

    # Get available time points
    times = get_time_from_csv('tr_data_001.csv')
    n_times = len(times) if times is not None else 0

    print("="*60)
    print("TR Energy Balance Analysis")
    print("="*60)

    if times is not None:
        print(f"\nAvailable time points: {n_times} (t = {times[0]:.3f} to {times[-1]:.3f} s)")

    # Determine if we use evolution files or snapshot files
    use_evolution = (time_index is not None and time_index >= 0 and time_index < n_times - 1)

    if use_evolution:
        actual_time = times[time_index] if times is not None else time_index
        print(f"Analyzing time index {time_index} (t = {actual_time:.4f} s)")

        # Read from evolution files
        data = read_evolution_data(time_index)
    else:
        if time_index is None or time_index == -1:
            print("Analyzing FINAL time point (using snapshot files)")
        else:
            print(f"Time index {time_index} -> using final time snapshot")

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

    if 'chi_e' in data and 'power_in' in data:
        estimate_temperature_from_balance(r, ne, chi_e, POH + PRF, Te)

    # Plot comparison
    plot_energy_balance(r, data)

    return data

def estimate_temperature_from_balance(r, n, chi, S, T_actual):
    """
    Estimate temperature from heat equation in steady state.

    Simplified 1D cylindrical heat equation:
    1/r * d/dr(r * n * chi * dT/dr) = -S / (1.5 * kB)

    Integrate twice to get T(r) given boundary condition T(a) = T_edge
    """

    # Convert units: S [MW/m^3], n [10^20/m^3], chi [m^2/s]
    # Heat capacity: 1.5 * n * kB, kB in keV = 1

    # For simple estimate, assume chi and n are roughly constant
    # Then: d^2T/dr^2 + (1/r)*dT/dr ≈ -S / (1.5 * n * chi)

    n_avg = np.mean(n[n > 0]) if np.any(n > 0) else 1.0
    chi_avg = np.mean(chi[chi > 0]) if np.any(chi > 0) else 1.0
    S_center = S[0] if len(S) > 0 else 0

    # Rough estimate of central temperature rise
    # T_center - T_edge ≈ S * a^2 / (4 * 1.5 * n * chi)
    # where a is minor radius (r/a = 1 at edge)

    a_norm = 1.0  # normalized radius

    # S in MW/m^3 = 1e6 W/m^3
    # n in 10^20/m^3
    # chi in m^2/s
    # T in keV, 1 keV = 1.6e-16 J

    # Conversion: S [MW/m^3] -> [keV / (10^20/m^3) / s]
    # 1 MW/m^3 = 1e6 J/s/m^3 = 1e6 / 1.6e-16 keV/s/m^3 = 6.25e21 keV/s/m^3
    # Per particle: 6.25e21 / 1e20 = 62.5 keV/s per (10^20/m^3)

    S_keV = S_center * 6.25  # approximate conversion factor for display

    dT_estimate = S_keV * a_norm**2 / (4 * 1.5 * chi_avg) if chi_avg > 0 else 0

    print(f"\nSimplified steady-state analysis:")
    print(f"  n_avg = {n_avg:.3f} [10^20/m^3]")
    print(f"  chi_avg = {chi_avg:.3f} [m^2/s]")
    print(f"  S_center = {S_center:.4f} [MW/m^3]")
    print(f"  Estimated dT ≈ {dT_estimate:.2f} keV (very rough)")
    print(f"  Actual T_center = {T_actual[0]:.3f} keV, T_edge = {T_actual[-1]:.3f} keV")
    print(f"  Actual dT = {T_actual[0] - T_actual[-1]:.3f} keV")

    # Check power balance
    # P_cond = n * chi * dT/dr integrated over surface
    dTdr = np.gradient(T_actual, r)
    q_cond = -n * chi * dTdr  # heat flux [keV * 10^20/m^3 * m^2/s / m] = [keV * 10^20 / m^2 / s]

    print(f"\nHeat flux at edge: q = {q_cond[-1]:.3f} [keV * 10^20/m^2/s]")

def plot_energy_balance(r, data):
    """Plot energy balance profiles."""

    fig, axes = plt.subplots(2, 3, figsize=(14, 8))

    # 1. Density profiles
    ax = axes[0, 0]
    if 'density' in data:
        for key, vals in data['density']['values'].items():
            # Handle potential dimension mismatch
            if len(vals) == len(r):
                ax.plot(r, vals, label=key)
            else:
                ax.plot(vals, label=key)  # Use index as x
    ax.set_xlabel('r/a')
    ax.set_ylabel('n [10^20/m^3]')
    ax.set_title('Density Profiles')
    ax.legend()
    ax.grid(True)

    # 2. Temperature profiles
    ax = axes[0, 1]
    if 'temperature' in data:
        for key, vals in data['temperature']['values'].items():
            if len(vals) == len(r):
                ax.plot(r, vals, label=key)
            else:
                ax.plot(vals, label=key)
    ax.set_xlabel('r/a')
    ax.set_ylabel('T [keV]')
    ax.set_title('Temperature Profiles')
    ax.legend()
    ax.grid(True)

    # 3. Transport coefficients
    ax = axes[0, 2]
    if 'chi_e' in data:
        for key, vals in data['chi_e']['values'].items():
            if len(vals) == len(r):
                ax.plot(r, vals, label=f'e: {key}')
            else:
                ax.plot(vals, label=f'e: {key}')
    if 'chi_i' in data:
        for key, vals in data['chi_i']['values'].items():
            if len(vals) == len(r):
                ax.plot(r, vals, '--', label=f'i: {key}')
            else:
                ax.plot(vals, '--', label=f'i: {key}')
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
                if len(vals) == len(r):
                    ax.plot(r, vals, label=key)
                else:
                    ax.plot(vals, label=key)
    ax.set_xlabel('r/a')
    ax.set_ylabel('P [MW/m^3]')
    ax.set_title('Power Input')
    ax.legend()
    ax.grid(True)

    # 5. Power loss
    ax = axes[1, 1]
    if 'power_loss' in data:
        for key, vals in data['power_loss']['values'].items():
            if len(vals) == len(r):
                ax.plot(r, vals, label=key)
            else:
                ax.plot(vals, label=key)
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

    parser = argparse.ArgumentParser(
        description='Analyze TR energy balance at a specific time point',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  python3 analyze_energy_balance.py              # Analyze final time (default)
  python3 analyze_energy_balance.py -t 0         # Analyze first time point
  python3 analyze_energy_balance.py -t 10        # Analyze time index 10
  python3 analyze_energy_balance.py -t -1        # Analyze final time
  python3 analyze_energy_balance.py --list-times # Show available time points
''')
    parser.add_argument('-t', '--time-index', type=int, default=None,
                        help='Time index (0=first, -1=final, or specific index)')
    parser.add_argument('--list-times', action='store_true',
                        help='List available time points and exit')

    args = parser.parse_args()

    if args.list_times:
        times = get_time_from_csv('tr_data_001.csv')
        if times is not None:
            print(f"Available time points: {len(times)}")
            print(f"Index range: 0 to {len(times)-1}")
            print(f"Time range: {times[0]:.4f} to {times[-1]:.4f} s")
            print("\nSample time points:")
            indices = np.linspace(0, len(times)-1, min(10, len(times)), dtype=int)
            for i in indices:
                print(f"  Index {i:3d}: t = {times[i]:.4f} s")
        else:
            print("Could not read time data from tr_data_001.csv")
    else:
        analyze_energy_balance(time_index=args.time_index)
