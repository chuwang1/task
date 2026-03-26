#!/usr/bin/env python3
"""
Calculate Zeff from TR CSV output and compare with the original TR Zeff.

Reads:
  - n(NS) CSV: thermal species densities (nE, nD, nT, nA=thermal He4)
  - NB,NF CSV: fast ion densities (NB=beam, NF=fast alpha)
  - ZEFF CSV: TR-calculated Zeff profile
  - PZC,PZFE CSV: carbon and iron charge states (function of Te)

Compares:
  - Zeff_thermal: calculated from thermal species only (no impurities)
  - Zeff_with_fast: including fast alpha contribution
  - Zeff_TR: original TR output
  - Impurity contribution: back-calculated from the difference
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import glob
import os


def find_csv_by_title(title_keyword, search_dir='.'):
    """Find TR CSV file by searching the title line for a keyword."""
    candidates = []
    for f in sorted(glob.glob(os.path.join(search_dir, 'tr_data_*.csv'))):
        with open(f, 'r') as fh:
            title = fh.readline().strip()
        if title_keyword in title:
            candidates.append(f)
    if not candidates:
        return None
    return os.path.basename(candidates[-1])


def read_tr_csv(filename, skip_title=True):
    """Read TR output CSV file"""
    df = pd.read_csv(filename, skiprows=1 if skip_title else 0)
    df.columns = df.columns.str.strip()
    return df


def main():
    print("="*70)
    print("Zeff Calculation from TR CSV Output")
    print("="*70)

    # Find CSV files by title
    file_map = {
        'dens':  find_csv_by_title('@n(NS)'),           # nE, nD, nT, nA
        'fast':  find_csv_by_title('@NB,NF'),            # fast ion densities
        'zeff':  find_csv_by_title('@ZEFF  vs r@'),      # TR Zeff
        'pz':    find_csv_by_title('@PZC,PZFE  vs r@'),  # charge states
    }

    print("\nCSV file mapping:")
    for key, fname in file_map.items():
        print(f"  {key:6s} -> {fname}")

    missing = [k for k, v in file_map.items() if v is None]
    if missing:
        print(f"\nERROR: Could not find CSV files for: {missing}")
        return

    # Read data
    df_dens = read_tr_csv(file_map['dens'])
    df_fast = read_tr_csv(file_map['fast'])
    df_zeff = read_tr_csv(file_map['zeff'])
    df_pz   = read_tr_csv(file_map['pz'])

    print(f"\nDensity columns: {list(df_dens.columns)}")
    print(f"Fast ion columns: {list(df_fast.columns)}")
    print(f"Zeff columns: {list(df_zeff.columns)}")
    print(f"PZ columns: {list(df_pz.columns)}")

    # Extract radial coordinate
    rho = df_dens['X'].values

    # Thermal species densities [10^20/m^3]
    ne  = df_dens['nE'].values   # electron density
    nD  = df_dens['nD'].values   # deuterium
    nT  = df_dens['nT'].values   # tritium
    nA  = df_dens['nA'].values   # thermal He4 (alpha)

    # Fast ion densities [10^20/m^3]
    fast_cols = [c for c in df_fast.columns if c != 'X']
    print(f"\nFast ion species ({len(fast_cols)}): {fast_cols}")

    # Typically: NB_1 = beam ion, NF_1 = fast alpha
    # NFMAX = NNBMAX + NNFMAX, last column(s) are fast alphas
    n_fast_total = np.zeros_like(rho)
    n_fast_alpha = np.zeros_like(rho)
    for i, col in enumerate(fast_cols):
        n_fast_total += df_fast[col].values

    # The last fast ion species is the fast alpha (from DT fusion)
    if len(fast_cols) >= 1:
        n_fast_alpha = df_fast[fast_cols[-1]].values
        print(f"  Fast alpha column: {fast_cols[-1]}")
    if len(fast_cols) >= 2:
        n_beam = df_fast[fast_cols[0]].values
        print(f"  Beam ion column: {fast_cols[0]}")

    # TR Zeff
    zeff_tr = df_zeff['ZEFF'].values
    rho_zeff = df_zeff['X'].values

    # Carbon and Iron charge states
    pzc  = df_pz['PZC'].values    # Z_C(Te)
    pzfe = df_pz['PZFE'].values   # Z_Fe(Te)
    rho_pz = df_pz['X'].values

    print(f"\n--- Density profiles at r/a ~ 0.5 ---")
    idx_mid = np.argmin(np.abs(rho - 0.5))
    print(f"  ne  = {ne[idx_mid]:.4f} x10^20/m^3")
    print(f"  nD  = {nD[idx_mid]:.4f}")
    print(f"  nT  = {nT[idx_mid]:.4f}")
    print(f"  nA (thermal) = {nA[idx_mid]:.6f}")
    print(f"  nF (fast alpha) = {n_fast_alpha[idx_mid]:.6f}")
    print(f"  nF/ne = {n_fast_alpha[idx_mid]/ne[idx_mid]*100:.2f}%")

    # === Calculate Zeff ===

    # 1. Zeff from thermal main ions only (no impurities)
    #    Z_D = Z_T = 1, Z_He4 = 2
    zeff_main = (1**2 * nD + 1**2 * nT + 2**2 * nA) / ne

    # 2. Zeff including fast alpha (Z=2)
    #    Need to also add fast alpha electrons to ne for consistency
    ne_total = ne  # ne already includes fast alpha electrons in TR
    zeff_main_fast = (1**2 * nD + 1**2 * nT + 2**2 * nA + 2**2 * n_fast_alpha) / ne_total

    # 3. Back-calculate impurity contribution from TR Zeff
    #    Zeff_TR = (sum Z_s^2 * n_s + Z_C^2 * n_C + Z_Fe^2 * n_Fe) / ne
    #    impurity_term = Zeff_TR * ne - sum(Z_s^2 * n_s)
    impurity_term = zeff_tr * ne - (1**2 * nD + 1**2 * nT + 2**2 * nA)
    # This equals PZC^2 * ANC + PZFE^2 * ANFE

    # 4. Zeff with fast alpha + impurities
    zeff_full = (1**2 * nD + 1**2 * nT + 2**2 * nA + 2**2 * n_fast_alpha + impurity_term) / ne_total

    print(f"\n{'='*70}")
    print(f"Zeff Comparison")
    print(f"{'='*70}")
    print(f"{'Location':<20} {'Zeff_TR':>10} {'Main only':>10} {'+ fast α':>10} {'diff':>10}")
    print(f"{'-'*60}")

    for label, idx in [("Center (r/a~0)", 0), ("r/a = 0.3", np.argmin(np.abs(rho-0.3))),
                        ("r/a = 0.5", idx_mid), ("r/a = 0.7", np.argmin(np.abs(rho-0.7))),
                        ("Edge (r/a~1)", -1)]:
        diff = zeff_full[idx] - zeff_tr[idx]
        print(f"  {label:<18} {zeff_tr[idx]:10.4f} {zeff_main[idx]:10.4f} "
              f"{zeff_main_fast[idx]:10.4f} {diff:+10.4f}")

    print(f"\n--- Fast alpha impact on Zeff ---")
    delta_zeff = zeff_main_fast - zeff_main
    print(f"  Max delta Zeff from fast alpha: {np.max(np.abs(delta_zeff)):.6f}")
    print(f"  Mean delta Zeff from fast alpha: {np.mean(np.abs(delta_zeff)):.6f}")

    print(f"\n--- Impurity contribution ---")
    print(f"  Impurity Zeff contribution (center): {impurity_term[0]/ne[0]:.4f}")
    print(f"  Impurity Zeff contribution (mid):    {impurity_term[idx_mid]/ne[idx_mid]:.4f}")
    print(f"  PZC (center) = {pzc[0]:.2f}, PZFE (center) = {pzfe[0]:.2f}")

    # === Plotting ===
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))

    # Panel 1: Zeff comparison
    ax = axes[0, 0]
    ax.plot(rho, zeff_tr, 'b-', lw=2, label='TR output')
    ax.plot(rho, zeff_main, 'g--', lw=1.5, label='Main ions only')
    ax.plot(rho, zeff_main_fast, 'r-.', lw=1.5, label='Main + fast α')
    ax.plot(rho, zeff_full, 'm:', lw=2, label='Full (+ impurity)')
    ax.set_xlabel('r/a')
    ax.set_ylabel('Zeff')
    ax.set_title('Zeff Comparison')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 1)

    # Panel 2: Density profiles
    ax = axes[0, 1]
    ax.plot(rho, ne, 'b-', lw=2, label='ne')
    ax.plot(rho, nD, 'g--', lw=1.5, label='nD')
    ax.plot(rho, nT, 'r--', lw=1.5, label='nT')
    ax.plot(rho, nA, 'purple', lw=1.5, ls='-.', label='nα (thermal)')
    ax.plot(rho, n_fast_alpha, 'orange', lw=2, label='nα (fast)')
    ax.set_xlabel('r/a')
    ax.set_ylabel('n [10²⁰/m³]')
    ax.set_title('Density Profiles')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 1)

    # Panel 3: Fast alpha impact
    ax = axes[1, 0]
    ax.plot(rho, delta_zeff, 'r-', lw=2, label='ΔZeff (fast α)')
    ax.plot(rho, impurity_term / ne, 'b--', lw=1.5, label='Impurity contribution')
    ax.set_xlabel('r/a')
    ax.set_ylabel('ΔZeff')
    ax.set_title('Zeff Contributions')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 1)

    # Panel 4: Charge states
    ax = axes[1, 1]
    ax.plot(rho_pz, pzc, 'b-', lw=2, label='Z_C (carbon)')
    ax.plot(rho_pz, pzfe, 'r-', lw=2, label='Z_Fe (iron)')
    ax.set_xlabel('r/a')
    ax.set_ylabel('Charge state Z')
    ax.set_title('Impurity Charge States')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 1)

    plt.tight_layout()
    plt.savefig('zeff_comparison.png', dpi=150)
    print(f"\nPlot saved: zeff_comparison.png")

    # Save data
    df_out = pd.DataFrame({
        'r_a': rho,
        'Zeff_TR': zeff_tr,
        'Zeff_main': zeff_main,
        'Zeff_main_fast': zeff_main_fast,
        'Zeff_full': zeff_full,
        'ne': ne, 'nD': nD, 'nT': nT,
        'nA_thermal': nA,
        'nA_fast': n_fast_alpha,
        'impurity_term': impurity_term,
    })
    df_out.to_csv('zeff_comparison.csv', index=False)
    print(f"Data saved: zeff_comparison.csv")


if __name__ == "__main__":
    main()
