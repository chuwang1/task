#!/usr/bin/env python3
"""
Reproduce TR chi profile at final time step using CDBM model (MDLKAI=132)

This script reads TR output CSV files and recalculates chi to verify
the transport coefficient calculation.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d

# Physical constants (same as Fortran code)
AEE = 1.602176487e-19   # elementary charge [C]
AME = 9.10938215e-31    # electron mass [kg]
AMP = 1.672621637e-27   # proton mass [kg]
VC = 2.99792458e8       # speed of light [m/s]
RMU0 = 4.0e-7 * np.pi   # permeability [H/m]
EPS0 = 1.0 / (VC**2 * RMU0)  # permittivity [F/m]
RKEV = 1.0e3 * AEE      # keV to J

# =============================================================================
# Machine parameters from tr.CFEDR0114.in
# =============================================================================
RR = 8.03   # Major radius [m]
RA = 2.72   # Minor radius [m]
BB = 6.0    # Toroidal field [T]
CK0 = 9.8   # Chi coefficient for electrons
CK1 = 9.8   # Chi coefficient for ions
CSPRS = 0.5 # Edge suppression factor
NREDGE_FRAC = 0.93  # Edge region starts at r/a = 0.93


def trcofs(shear, alpha, curv):
    """
    Form factor in CDBM model (identical to Fortran trcofs function)
    """
    shear = np.atleast_1d(shear)
    alpha = np.atleast_1d(alpha)
    curv = np.atleast_1d(curv)

    fs = np.zeros_like(shear)

    for i in range(len(shear)):
        if alpha[i] >= 0:
            sa = shear[i] - alpha[i]
            if sa >= 0:
                denom = np.sqrt(2.0) * (1.0 - 2.0*sa + 3.0*sa*sa + 2.0*sa**3)
                if abs(denom) > 1e-10:
                    fs1 = (1.0 + 9.0 * np.sqrt(2.0) * sa**2.5) / denom
                else:
                    fs1 = 0.0
            else:
                denom = 2.0 * (1.0 - 2.0*sa) * (1.0 - 2.0*sa + 3.0*sa*sa)
                if denom > 0:
                    fs1 = 1.0 / np.sqrt(denom)
                else:
                    fs1 = 0.0

            if curv[i] > 0 and abs(shear[i]) > 0.01:
                fs2 = np.sqrt(curv[i])**3 / (shear[i]**2)
            else:
                fs2 = 0.0
        else:
            # alpha < 0
            sa = alpha[i] - shear[i]
            if sa >= 0:
                denom = np.sqrt(2.0) * (1.0 - 2.0*sa + 3.0*sa*sa + 2.0*sa**3)
                if abs(denom) > 1e-10:
                    fs1 = (1.0 + 9.0 * np.sqrt(2.0) * sa**2.5) / denom
                else:
                    fs1 = 0.0
            else:
                denom = 2.0 * (1.0 - 2.0*sa) * (1.0 - 2.0*sa + 3.0*sa*sa)
                if denom > 0:
                    fs1 = 1.0 / np.sqrt(denom)
                else:
                    fs1 = 0.0

            if curv[i] < 0 and abs(shear[i]) > 0.01:
                fs2 = np.sqrt(-curv[i])**3 / (shear[i]**2)
            else:
                fs2 = 0.0

        fs[i] = max(fs1, fs2)

    return fs


def calculate_cdbm_chi_from_alpha(BB, RR, rs, rkap, qp, shear, ne, alpha, rhoni,
                                   dvexbdr=None, calf=1.0, ckap=1.0, cexb=1.0, model=2):
    """
    Calculate CDBM thermal diffusivity using TR's alpha directly

    model: 0=CDBM original, 1=CDBM05, 2=CDBM+weak ExB, 3=CDBM05+weak ExB
    """
    ckcdbm = 12.0  # Fixed numerical factor

    # TR code applies minimum shear of 0.5 (trcoef.f90 line 878)
    shear_eff = np.maximum(shear, 0.5)

    # Alfven velocity
    va = np.sqrt(BB**2 / (RMU0 * rhoni))

    # Square of plasma frequency
    wpe2 = ne * AEE**2 / (AME * EPS0)

    # Square of collisionless skin depth
    delta2 = VC**2 / wpe2

    # Magnetic curvature
    curv = -(rs / RR) * (1.0 - 1.0 / (qp**2))

    # Rotational shear for ExB calculation
    shearl = np.sqrt(shear_eff**2 + 0.1**2)

    # ExB shear effect
    if dvexbdr is None:
        dvexbdr = np.zeros_like(rs)
    wexb = -qp * RR / (shearl * va) * dvexbdr

    # Model-dependent factors
    # fk: elongation factor
    if model % 2 == 0:
        fk = 1.0
    else:
        # CDBM05 with elongation
        fk = (2.0 * np.sqrt(rkap) / (1.0 + rkap**2))**1.5

    # fe: ExB shear reduction factor
    model_type = (model // 2) % 3
    if model_type == 0:
        fe = 1.0
    elif model_type == 1:
        # Weak ExB shear (model=2,3)
        fe = 1.0 / (1.0 + cexb * wexb**2)
    else:
        fe = 1.0

    # Form factor using TR's alpha directly
    fs = trcofs(shear_eff, calf * alpha, ckap * curv)

    # CDBM chi
    chi_cdbm = ckcdbm * fs * fk * fe * np.abs(alpha)**1.5 * delta2 * va / (qp * RR)

    return chi_cdbm, fs, fe, shear_eff, curv, va, delta2


def read_tr_csv(filename, skip_title=True):
    """Read TR output CSV file"""
    df = pd.read_csv(filename, skiprows=1 if skip_title else 0)
    df.columns = df.columns.str.strip()
    return df


def main():
    print("="*70)
    print("Reproducing TR chi profile using CDBM model (MDLKAI=132)")
    print("="*70)

    # Read TR output files
    try:
        # Read profiles at final time step
        df_q = read_tr_csv('tr_data_053.csv')      # QP (safety factor)
        df_s = read_tr_csv('tr_data_061.csv')      # s (magnetic shear)
        df_alpha = read_tr_csv('tr_data_063.csv')  # alpha
        df_ne = read_tr_csv('tr_data_045.csv')     # NE
        df_Te = read_tr_csv('tr_data_047.csv')     # TE
        df_Ti = read_tr_csv('tr_data_048.csv')     # TD
        df_chi = read_tr_csv('tr_data_064.csv')    # AKD (ion chi from TR)

        # Also read ion densities for accurate rhoni calculation
        # tr_data_017.csv contains all densities: nE, nD, nT, nA (alpha)
        df_dens = read_tr_csv('tr_data_017.csv')

        # Get r/a coordinates from different grids
        rho_ne = df_ne['X'].values   # Mid-points (50 points: 0.01, 0.03, ..., 0.99)
        rho_q = df_q['X'].values     # Grid points (51 points: 0.00, 0.02, ..., 1.00)
        rho_chi = df_chi['X'].values # Chi grid (49 points for gradient calculation)

        # Get final time step values (last column with data)
        q_cols = [c for c in df_q.columns if c.startswith('QP_')]
        s_cols = [c for c in df_s.columns if c.startswith('s_')]
        alpha_cols = [c for c in df_alpha.columns if c.startswith('alpha_')]
        ne_cols = [c for c in df_ne.columns if c.startswith('NE_')]
        Te_cols = [c for c in df_Te.columns if c.startswith('TE_')]
        Ti_cols = [c for c in df_Ti.columns if c.startswith('TD_')]
        chi_cols = [c for c in df_chi.columns if c.startswith('AKD_')]

        # Sort by time step number and get last
        q_raw = df_q[sorted(q_cols, key=lambda x: int(x.split('_')[1]))[-1]].values
        s_raw = df_s[sorted(s_cols, key=lambda x: int(x.split('_')[1]))[-1]].values
        alpha_raw = df_alpha[sorted(alpha_cols, key=lambda x: int(x.split('_')[1]))[-1]].values
        ne_raw = df_ne[sorted(ne_cols, key=lambda x: int(x.split('_')[1]))[-1]].values * 1e20
        Te_raw = df_Te[sorted(Te_cols, key=lambda x: int(x.split('_')[1]))[-1]].values
        Ti_raw = df_Ti[sorted(Ti_cols, key=lambda x: int(x.split('_')[1]))[-1]].values
        chi_tr = df_chi[sorted(chi_cols, key=lambda x: int(x.split('_')[1]))[-1]].values

        # Read ion densities from tr_data_017.csv (nD, nT, nA columns)
        nD_raw = df_dens['nD'].values * 1e20  # D density [m^-3]
        nT_raw = df_dens['nT'].values * 1e20  # T density [m^-3]
        nHe_raw = df_dens['nA'].values * 1e20  # He/alpha density [m^-3]
        rho_dens = df_dens['X'].values

        # Use chi grid as the common grid (since we want to compare with chi_tr)
        rho = rho_chi.copy()

        # Interpolate all variables to the chi grid
        q = interp1d(rho_q, q_raw, kind='linear', fill_value='extrapolate')(rho)
        s = interp1d(rho_ne, s_raw, kind='linear', fill_value='extrapolate')(rho)
        alpha_tr = interp1d(rho_ne, alpha_raw, kind='linear', fill_value='extrapolate')(rho)
        ne = interp1d(rho_ne, ne_raw, kind='linear', fill_value='extrapolate')(rho)
        Te = interp1d(rho_ne, Te_raw, kind='linear', fill_value='extrapolate')(rho)
        Ti = interp1d(rho_ne, Ti_raw, kind='linear', fill_value='extrapolate')(rho)
        nD = interp1d(rho_dens, nD_raw, kind='linear', fill_value='extrapolate')(rho)
        nT = interp1d(rho_dens, nT_raw, kind='linear', fill_value='extrapolate')(rho)
        nHe = interp1d(rho_dens, nHe_raw, kind='linear', fill_value='extrapolate')(rho)

        print(f"\nRead TR profiles (final time step):")
        print(f"  r/a points: {len(rho)}")
        print(f"  q range: {q.min():.3f} - {q.max():.3f}")
        print(f"  s range: {s.min():.3f} - {s.max():.3f}")
        print(f"  alpha range: {alpha_tr.min():.3f} - {alpha_tr.max():.3f}")
        print(f"  ne range: {ne.min()/1e20:.3f} - {ne.max()/1e20:.3f} [10^20 m^-3]")
        print(f"  Te range: {Te.min():.2f} - {Te.max():.2f} keV")
        print(f"  chi_TR range: {chi_tr.min():.3f} - {chi_tr.max():.3f} m^2/s")

    except Exception as e:
        print(f"Error reading TR files: {e}")
        return

    # Calculate derived quantities
    rs = rho * RA  # Minor radius [m]
    rkap = 1.89    # Elongation from input file

    # Debug: Print array shapes
    print(f"\nArray shapes:")
    print(f"  rho: {rho.shape}, rs: {rs.shape}")
    print(f"  ne: {ne.shape}, Te: {Te.shape}, Ti: {Ti.shape}")
    print(f"  q: {q.shape}, s: {s.shape}")
    print(f"  nD: {nD.shape}, nT: {nT.shape}, nHe: {nHe.shape}")

    # Ion mass density using actual TR densities
    # rhoni = sum(mass_i * n_i) for all ion species
    rhoni = (2.0 * AMP * nD) + (3.0 * AMP * nT) + (4.0 * AMP * nHe)

    print(f"\nDerived quantities:")
    print(f"  rhoni range: {rhoni.min():.3e} - {rhoni.max():.3e} kg/m^3")

    # Calculate CDBM chi using TR's alpha directly
    chi_cdbm, fs, fe, shear_eff, curv, va, delta2 = calculate_cdbm_chi_from_alpha(
        BB, RR, rs, rkap, q, s, ne, alpha_tr, rhoni, model=2
    )

    # Apply CK factor: chi = (CK/12) * chi_cdbm
    chi_calc = (CK1 / 12.0) * chi_cdbm

    # Apply edge model (MDLEDGE=1)
    # TR uses NREDGE=NINT(0.93*NRMAX), NR>=NREDGE triggers edge model
    # Find index where r/a >= 0.93
    nredge = np.searchsorted(rho, NREDGE_FRAC)
    chi_calc_edge = chi_calc.copy()
    chi_calc_edge[nredge:] = CSPRS * chi_calc_edge[nredge:]

    print(f"  Edge model starts at index {nredge}, r/a = {rho[nredge]:.4f}")

    print(f"\nCalculation parameters:")
    print(f"  RR = {RR} m, RA = {RA} m, BB = {BB} T")
    print(f"  CK0 = CK1 = {CK1}")
    print(f"  Edge model: CSPRS = {CSPRS} for r/a > {NREDGE_FRAC}")
    print(f"  Shear minimum: 0.5 (as in TR code)")

    print(f"\nCalculated chi range: {chi_calc.min():.3f} - {chi_calc.max():.3f} m^2/s")
    print(f"With edge model: {chi_calc_edge.min():.3f} - {chi_calc_edge.max():.3f} m^2/s")

    # Compare
    diff = chi_calc_edge - chi_tr
    rel_diff = np.abs(diff) / (chi_tr + 1e-10) * 100

    print(f"\nComparison with TR output:")
    print(f"  Max absolute difference: {np.max(np.abs(diff)):.4f} m^2/s")
    print(f"  Mean relative difference: {np.mean(rel_diff[rho < 0.9]):.1f}%")

    # Plot comparison
    fig, axes = plt.subplots(2, 3, figsize=(14, 9))

    # Chi comparison
    ax1 = axes[0, 0]
    ax1.plot(rho, chi_tr, 'b-', lw=2, label='TR output (AKD)')
    ax1.plot(rho, chi_calc_edge, 'r--', lw=2, label='Python calc')
    ax1.axvline(NREDGE_FRAC, color='gray', ls=':', label=f'Edge (r/a={NREDGE_FRAC})')
    ax1.set_xlabel('r/a')
    ax1.set_ylabel(r'$\chi$ [m$^2$/s]')
    ax1.set_title('Chi Comparison (Final Time)')
    ax1.legend()
    ax1.grid(True)
    ax1.set_xlim(0, 1)

    # Chi log scale
    ax2 = axes[0, 1]
    ax2.semilogy(rho, chi_tr, 'b-', lw=2, label='TR output')
    ax2.semilogy(rho, chi_calc_edge, 'r--', lw=2, label='Python calc')
    ax2.set_xlabel('r/a')
    ax2.set_ylabel(r'$\chi$ [m$^2$/s]')
    ax2.set_title('Chi (Log Scale)')
    ax2.legend()
    ax2.grid(True)
    ax2.set_xlim(0, 1)

    # Alpha and curvature
    ax3 = axes[0, 2]
    ax3.plot(rho, alpha_tr, 'b-', lw=2, label=r'$\alpha$ (TR)')
    ax3.plot(rho, curv, 'r--', lw=2, label='curv')
    ax3.set_xlabel('r/a')
    ax3.set_ylabel(r'$\alpha$, curv')
    ax3.set_title('Alpha and Magnetic Curvature')
    ax3.legend()
    ax3.grid(True)
    ax3.set_xlim(0, 1)

    # Magnetic shear
    ax4 = axes[1, 0]
    ax4.plot(rho, s, 'b-', lw=2, label='TR shear (s)')
    ax4.plot(rho, shear_eff, 'r--', lw=2, label='Effective shear (min=0.5)')
    ax4.axhline(0.5, color='gray', ls=':', label='s_min = 0.5')
    ax4.set_xlabel('r/a')
    ax4.set_ylabel('s = (r/q)(dq/dr)')
    ax4.set_title('Magnetic Shear')
    ax4.legend()
    ax4.grid(True)
    ax4.set_xlim(0, 1)

    # Form factor
    ax5 = axes[1, 1]
    ax5.plot(rho, fs, 'b-', lw=2)
    ax5.set_xlabel('r/a')
    ax5.set_ylabel(r'$f_s$')
    ax5.set_title('CDBM Form Factor')
    ax5.grid(True)
    ax5.set_xlim(0, 1)

    # Safety factor
    ax6 = axes[1, 2]
    ax6.plot(rho, q, 'b-', lw=2)
    ax6.set_xlabel('r/a')
    ax6.set_ylabel('q')
    ax6.set_title('Safety Factor')
    ax6.grid(True)
    ax6.set_xlim(0, 1)

    plt.tight_layout()
    plt.savefig('reproduce_chi_profile.png', dpi=150)
    print(f"\nPlot saved: reproduce_chi_profile.png")

    # Save comparison to CSV
    df_out = pd.DataFrame({
        'r_a': rho,
        'chi_TR': chi_tr,
        'chi_calc': chi_calc_edge,
        'chi_diff': diff,
        'rel_diff_pct': rel_diff,
        'alpha': alpha_tr,
        's': s,
        's_eff': shear_eff,
        'fs': fs,
        'curv': curv,
        'q': q,
        'va': va,
        'delta2': delta2
    })
    df_out.to_csv('chi_comparison.csv', index=False)
    print(f"Data saved: chi_comparison.csv")


if __name__ == "__main__":
    main()
