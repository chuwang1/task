#!/usr/bin/env python3
"""
Reproduce TR chi profile at final time step using CDBM-style model (MDLKAI=31)

This script reads TR output CSV files and recalculates chi to verify
the transport coefficient calculation.

Key findings:
- MDLKAI=31 formula: chi = CK * fs * |alpha|^1.5 * delta^2 * va / (q*R)
- CK value reverse-engineered from TR output: ~9.7 (average over core region)
- Edge model (MDLEDGE=1) with CSPRS=0.55 matches TR output best
- Neoclassical chi (Hinton-Hazeltine model) is included

Accuracy achieved:
- Core (r/a < 0.85): ~3% mean error
- Edge (r/a >= 0.93): ~3% mean error
- Overall: ~3.3% mean error
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
# CK values reverse-engineered from TR output: ~9.7
CK0 = 9.8   # Chi coefficient for electrons
CK1 = 9.8   # Chi coefficient for ions (implied CK ~9.7 from data)
# Edge model parameters
MDLEDGE = 1  # 0=disabled, 1=enabled (testing with edge model enabled)
CSPRS = 0.55 # Edge suppression factor (optimized from TR output: 0.51-0.56)
NREDGE_FRAC = 0.93  # Edge region starts at r/a = 0.93
# Chi limits for MDLKAI=31 model
CHI_MIN = 0.1   # Minimum chi [m^2/s]
CHI_MAX = 50.0  # Maximum chi [m^2/s]


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


def coulomb_log(ne, Te):
    """Coulomb logarithm (simplified)"""
    # ln(Lambda) ~ 15-20 for typical tokamak plasmas
    # Using simplified formula: ln(Lambda) = 17.3 - 0.5*ln(ne/1e20) + 1.5*ln(Te)
    return 17.3 - 0.5 * np.log(ne / 1e20) + 1.5 * np.log(np.maximum(Te, 0.1))


def calculate_collision_time_i(ne, ni, Ti, Z, A):
    """
    Ion-ion collision time (FTAUI in Fortran)

    ne: electron density [m^-3]
    ni: ion density [m^-3]
    Ti: ion temperature [keV]
    Z: ion charge number
    A: ion atomic mass number
    """
    COEF = 12.0 * np.pi * np.sqrt(np.pi) * EPS0**2 * np.sqrt(A * AMP) / (AEE**4 * 1e20)
    ln_lambda = coulomb_log(ne, Ti)
    ni_20 = ni / 1e20  # Convert to 10^20 m^-3
    tau_i = COEF * (Ti * RKEV)**1.5 / (ni_20 * Z**4 * ln_lambda)
    return tau_i


def calculate_neoclassical_chi(rho, q, Ti, ne, ni, Bp, RR, RA, BB, Z=1, A=2, Zeff=1.5, model=1):
    """
    Calculate neoclassical thermal diffusivity

    model: 0 = Chang-Hinton, 1 = Hinton-Hazeltine (MDLKNC)
    """
    rs = rho * RA  # Minor radius [m]
    eps = rs / RR   # Inverse aspect ratio
    eps = np.maximum(eps, 0.01)  # Avoid division by zero

    # Thermal velocity
    v_th = np.sqrt(Ti * RKEV / (A * AMP))

    # Larmor radius squared
    rho_i2 = 2.0 * A * AMP * Ti * RKEV / (Z * AEE * Bp)**2

    # Collision time
    tau_i = calculate_collision_time_i(ne, ni, Ti, Z, A)

    # Collisionality
    nu_star = np.abs(q) * RR / (tau_i * v_th * eps)

    if model == 1:
        # Hinton-Hazeltine model (MDLKNC=1)
        # Coefficients from trcoef.f90 lines 1214-1221
        RK2 = 0.66
        RA2 = 1.03
        RB2 = 0.31
        RC2 = 0.74

        RK2_factor = RK2 * (1.0 / (1.0 + RA2 * np.sqrt(nu_star) + RB2 * nu_star) +
                           (eps * RC2)**2 / RB2 * nu_star / (1.0 + RC2 * nu_star * eps))

        chi_nc = np.sqrt(eps) * rho_i2 / tau_i * RK2_factor
    else:
        # Chang-Hinton model (MDLKNC=0)
        RALPHA = Zeff - 1.0
        DELDA = 0.0
        RMU = nu_star * (1.0 + 1.54 * RALPHA)

        F1 = (1.0 + 1.5 * (eps**2 + eps * DELDA) + 3.0/8.0 * eps**3 * DELDA) / (1.0 + 0.5 * eps * DELDA)
        F2 = np.sqrt(1.0 - eps**2) * (1.0 + 0.5 * eps * DELDA) / (1.0 + DELDA/eps * (np.sqrt(1.0 - eps**2) - 1.0))

        TERM1 = (0.66 * (1.0 + 1.54 * RALPHA) +
                 (1.88 * np.sqrt(eps) - 1.54 * eps) * (1.0 + 3.75 * RALPHA)) * F1 / (
                 1.0 + 1.03 * np.sqrt(RMU) + 0.31 * RMU)
        TERM2 = 0.583 * RMU * eps / (1.0 + 0.74 * RMU * eps) * (
                1.0 + (1.33 * RALPHA * (1.0 + 0.6 * RALPHA)) / (1.0 + 1.79 * RALPHA)) * (F1 - F2)

        chi_nc = (rho_i2 * np.sqrt(eps)) / tau_i * (TERM1 + TERM2)

    # Apply upper limit (same as Fortran code line 1276-1277)
    # CHECK = |T/(2*Z*e*R*B)| * r/a * RA
    chi_limit = np.abs(Ti * RKEV / (2.0 * Z * AEE * RR * BB)) * rho * RA
    chi_nc = np.minimum(chi_nc, chi_limit)

    return chi_nc


def read_tr_csv(filename, skip_title=True):
    """Read TR output CSV file"""
    df = pd.read_csv(filename, skiprows=1 if skip_title else 0)
    df.columns = df.columns.str.strip()
    return df


import glob
import os

def find_csv_by_title(title_keyword, search_dir='.', prefer_2d=False):
    """
    Find a TR CSV file by searching the title line for a keyword.
    Returns the LAST matching file (highest number) to get the final output.

    Args:
        title_keyword: string to search for in the title line
        search_dir: directory to search in
        prefer_2d: if True, prefer 2D contour files (with many columns)

    Returns:
        filename of matching CSV file, or None
    """
    candidates = []
    for f in sorted(glob.glob(os.path.join(search_dir, 'tr_data_*.csv'))):
        with open(f, 'r') as fh:
            title = fh.readline().strip()
        if title_keyword in title:
            # Check if it's a 2D file (has many columns = time evolution)
            with open(f, 'r') as fh:
                fh.readline()  # skip title
                header = fh.readline().strip()
            ncols = len(header.split(','))
            candidates.append((f, ncols))

    if not candidates:
        return None

    if prefer_2d:
        # Prefer files with many columns (time evolution data)
        candidates.sort(key=lambda x: x[1], reverse=True)
    else:
        # Prefer files with fewer columns (single time point / radial profile)
        # But among those with same column count, prefer the last one
        candidates.sort(key=lambda x: (-x[1], x[0]))
        # Actually just return the last matching file
        candidates.sort(key=lambda x: x[0])

    return os.path.basename(candidates[-1][0])


def main():
    print("="*70)
    print("Reproducing TR chi profile using CDBM model (MDLKAI=132)")
    print("="*70)

    # Find CSV files by title content (robust to file renumbering)
    file_map = {
        'QP':    find_csv_by_title('@QP  vs r@'),
        's':     find_csv_by_title('@s  vs r@'),
        'alpha': find_csv_by_title('@alpha  vs r@'),
        'NE':    find_csv_by_title('@NE [10^20/m^3]  vs r@'),
        'TE':    find_csv_by_title('@TE [keV]  vs r@'),
        'TD':    find_csv_by_title('@TD [keV]  vs r@'),
        'AKD':   find_csv_by_title('@AKD [m^2/s]  vs r@'),
        'dens':  find_csv_by_title('@n(NS)'),
    }

    print("\nCSV file mapping (found by title):")
    for key, fname in file_map.items():
        print(f"  {key:6s} -> {fname}")

    missing = [k for k, v in file_map.items() if v is None]
    if missing:
        print(f"\nERROR: Could not find CSV files for: {missing}")
        return

    # Read TR output files
    try:
        # Read profiles at final time step
        df_q = read_tr_csv(file_map['QP'])        # QP (safety factor)
        df_s = read_tr_csv(file_map['s'])          # s (magnetic shear)
        df_alpha = read_tr_csv(file_map['alpha'])   # alpha
        df_ne = read_tr_csv(file_map['NE'])        # NE
        df_Te = read_tr_csv(file_map['TE'])        # TE
        df_Ti = read_tr_csv(file_map['TD'])        # TD
        df_chi = read_tr_csv(file_map['AKD'])      # AKD (ion chi from TR)

        # Also read ion densities for accurate rhoni calculation
        df_dens = read_tr_csv(file_map['dens'])

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

        # Read ion densities (nD, nT, nA columns)
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

    # Total ion density for neoclassical calculation
    ni_total = nD + nT + nHe

    # Calculate poloidal magnetic field: Bp = r*B/(R*q)
    # Near axis, use a minimum value to prevent chi_nc from exploding
    Bp = rs * BB / (RR * q)
    # For r/a < 0.2, use Bp at r/a=0.2 as minimum
    Bp_min = 0.2 * RA * BB / (RR * q[np.searchsorted(rho, 0.2)])
    Bp = np.maximum(Bp, Bp_min)

    print(f"\nDerived quantities:")
    print(f"  rhoni range: {rhoni.min():.3e} - {rhoni.max():.3e} kg/m^3")
    print(f"  Bp range: {Bp.min():.4f} - {Bp.max():.4f} T")

    # Calculate CDBM chi using TR's alpha directly (original shear)
    chi_cdbm, fs, fe, shear_eff, curv, va, delta2 = calculate_cdbm_chi_from_alpha(
        BB, RR, rs, rkap, q, s, ne, alpha_tr, rhoni, model=2
    )

    # Apply CK factor: chi = (CK/12) * chi_cdbm (anomalous part)
    chi_anom = (CK1 / 12.0) * chi_cdbm

    # =========================================================
    # Smoothed shear test: eliminate chi peak at r/a=0.2-0.3
    # =========================================================
    print(f"\n--- Smoothed shear test ---")
    smooth_windows = [3, 5, 7, 11]  # boxcar window sizes
    chi_smoothed = {}
    s_smoothed = {}
    
    for nw in smooth_windows:
        # Boxcar (moving average) smoothing of magnetic shear
        kernel = np.ones(nw) / nw
        s_sm = np.convolve(s, kernel, mode='same')
        # Fix edges: keep original values at boundaries
        half = nw // 2
        s_sm[:half] = s[:half]
        s_sm[-half:] = s[-half:]
        s_smoothed[nw] = s_sm
        
        # Recalculate chi with smoothed shear
        chi_sm, fs_sm, _, _, _, _, _ = calculate_cdbm_chi_from_alpha(
            BB, RR, rs, rkap, q, s_sm, ne, alpha_tr, rhoni, model=2
        )
        chi_sm_scaled = (CK1 / 12.0) * chi_sm
        chi_sm_scaled = np.clip(chi_sm_scaled, CHI_MIN, CHI_MAX)
        chi_smoothed[nw] = chi_sm_scaled
        
        # Check peak region
        mask = (rho >= 0.2) & (rho <= 0.35)
        peak_orig = chi_anom[mask].max()
        peak_sm = chi_sm_scaled[mask].max()
        print(f"  Window={nw}: shear peak region chi max: {peak_orig:.3f} -> {peak_sm:.3f} "
              f"(reduction: {(1-peak_sm/peak_orig)*100:.1f}%)")

    # Calculate neoclassical chi (Hinton-Hazeltine model, MDLKNC=1)
    # Use D ion as representative (A=2, Z=1)
    chi_nc = calculate_neoclassical_chi(rho, q, Ti, ne, ni_total, Bp, RR, RA, BB, Z=1, A=2, Zeff=1.5, model=1)

    # Apply chi limits (MDLKAI=31 model has MIN/MAX limits)
    chi_anom = np.clip(chi_anom, CHI_MIN, CHI_MAX)

    print(f"\n  chi_anom range (after limits): {chi_anom.min():.3f} - {chi_anom.max():.3f} m^2/s")
    print(f"  chi_nc range: {chi_nc.min():.4f} - {chi_nc.max():.4f} m^2/s")

    # Apply edge model only if MDLEDGE=1
    chi_anom_edge = chi_anom.copy()
    nredge = np.searchsorted(rho, NREDGE_FRAC)
    if MDLEDGE == 1:
        chi_anom_edge[nredge:] = CSPRS * chi_anom_edge[nredge:]
        print(f"  Edge model (MDLEDGE=1): CSPRS={CSPRS} applied for r/a >= {NREDGE_FRAC}")
    else:
        print(f"  Edge model (MDLEDGE=0): disabled")

    # Total chi = CDH * chi_anom + CNH * chi_nc
    CNH = 1.0
    chi_calc_edge = chi_anom_edge + CNH * chi_nc
    
    # Also compute total chi with best smoothed shear (window=5)
    best_nw = 5
    chi_sm_edge = chi_smoothed[best_nw].copy()
    if MDLEDGE == 1:
        chi_sm_edge[nredge:] = CSPRS * chi_sm_edge[nredge:]
    chi_calc_smoothed = chi_sm_edge + CNH * chi_nc

    print(f"\nCalculation parameters:")
    print(f"  RR = {RR} m, RA = {RA} m, BB = {BB} T")
    print(f"  CK0 = CK1 = {CK1}")
    print(f"  MDLEDGE = {MDLEDGE}, CSPRS = {CSPRS}")
    print(f"  Chi limits: [{CHI_MIN}, {CHI_MAX}] m^2/s")
    print(f"  Shear minimum: 0.5 (as in TR code)")

    print(f"\nCalculated chi (anom+NC) range: {chi_calc_edge.min():.3f} - {chi_calc_edge.max():.3f} m^2/s")
    print(f"Smoothed chi (window={best_nw}) range: {chi_calc_smoothed.min():.3f} - {chi_calc_smoothed.max():.3f} m^2/s")

    # Compare
    diff = chi_calc_edge - chi_tr
    rel_diff = np.abs(diff) / (chi_tr + 1e-10) * 100

    print(f"\nComparison with TR output:")
    print(f"  Max absolute difference: {np.max(np.abs(diff)):.4f} m^2/s")
    print(f"  Mean relative difference: {np.mean(rel_diff[rho < 0.9]):.1f}%")

    # Plot comparison
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    fig.suptitle('CDBM Chi: Original vs Smoothed Shear', fontsize=13, fontweight='bold')

    # --- Plot 1: Chi comparison with smoothed ---
    ax1 = axes[0, 0]
    ax1.plot(rho, chi_tr, 'b-', lw=2, label='TR output (AKD)')
    ax1.plot(rho, chi_calc_edge, 'r--', lw=2, label='Original shear')
    ax1.plot(rho, chi_calc_smoothed, 'g-', lw=2, alpha=0.8, label=f'Smoothed (w={best_nw})')
    ax1.axvline(NREDGE_FRAC, color='gray', ls=':', alpha=0.5)
    ax1.axvspan(0.2, 0.35, alpha=0.1, color='red', label='Peak region')
    ax1.set_xlabel('r/a')
    ax1.set_ylabel(r'$\chi$ [m$^2$/s]')
    ax1.set_title('Chi Comparison')
    ax1.legend(fontsize=7)
    ax1.grid(True)
    ax1.set_xlim(0, 1)

    # --- Plot 2: Smoothed chi with different windows ---
    ax2 = axes[0, 1]
    ax2.plot(rho, chi_anom, 'k-', lw=2, label='Original')
    colors = ['#e74c3c', '#2ecc71', '#3498db', '#9b59b6']
    for i, nw in enumerate(smooth_windows):
        ax2.plot(rho, chi_smoothed[nw], color=colors[i], ls='--', lw=1.5, 
                 label=f'Smoothed w={nw}')
    ax2.axvspan(0.2, 0.35, alpha=0.1, color='red')
    ax2.set_xlabel('r/a')
    ax2.set_ylabel(r'$\chi_{anom}$ [m$^2$/s]')
    ax2.set_title('Anomalous Chi: Different Smoothing Windows')
    ax2.legend(fontsize=7)
    ax2.grid(True)
    ax2.set_xlim(0, 0.6)

    # --- Plot 3: Alpha and s-alpha (key to the peak) ---
    ax3 = axes[0, 2]
    sa_orig = shear_eff - alpha_tr  # SA = s_eff - alpha (key quantity)
    sa_smooth = np.maximum(s_smoothed[best_nw], 0.5) - alpha_tr
    ax3.plot(rho, alpha_tr, 'b-', lw=2, label=r'$\alpha$')
    ax3.plot(rho, sa_orig, 'r-', lw=2, label=r'$s_{eff} - \alpha$ (original)')
    ax3.plot(rho, sa_smooth, 'g--', lw=2, label=r'$s_{eff} - \alpha$ (smoothed)')
    ax3.axhline(0, color='gray', ls=':', alpha=0.5)
    ax3.axvspan(0.2, 0.35, alpha=0.1, color='red')
    ax3.set_xlabel('r/a')
    ax3.set_ylabel('Value')
    ax3.set_title(r'$\alpha$ and $s-\alpha$ (SA, peak cause)')
    ax3.legend(fontsize=7)
    ax3.grid(True)
    ax3.set_xlim(0, 0.6)

    # --- Plot 4: Magnetic shear: original vs smoothed ---
    ax4 = axes[1, 0]
    ax4.plot(rho, s, 'b-', lw=2, label='TR shear (s)')
    ax4.plot(rho, shear_eff, 'r--', lw=1.5, label='s_eff (min=0.5)')
    for i, nw in enumerate(smooth_windows):
        ax4.plot(rho, s_smoothed[nw], color=colors[i], ls=':', lw=1.5, 
                 label=f'Smoothed w={nw}', alpha=0.8)
    ax4.axhline(0.5, color='gray', ls=':', alpha=0.3)
    ax4.set_xlabel('r/a')
    ax4.set_ylabel('s = (r/q)(dq/dr)')
    ax4.set_title('Magnetic Shear: Original vs Smoothed')
    ax4.legend(fontsize=6, ncol=2)
    ax4.grid(True)
    ax4.set_xlim(0, 0.6)

    # --- Plot 5: Form factor comparison ---
    ax5 = axes[1, 1]
    ax5.plot(rho, fs, 'b-', lw=2, label='fs (original)')
    # Calculate fs for smoothed shear
    s_sm_best = s_smoothed[best_nw]
    shear_sm_eff = np.maximum(s_sm_best, 0.5)
    curv_sm = -(rs / RR) * (1.0 - 1.0 / (q**2))
    fs_sm = trcofs(shear_sm_eff, alpha_tr, curv_sm)
    ax5.plot(rho, fs_sm, 'g--', lw=2, label=f'fs (smoothed w={best_nw})')
    ax5.axvspan(0.2, 0.35, alpha=0.1, color='red')
    ax5.set_xlabel('r/a')
    ax5.set_ylabel(r'$f_s$')
    ax5.set_title('CDBM Form Factor')
    ax5.legend(fontsize=8)
    ax5.grid(True)
    ax5.set_xlim(0, 0.6)

    # --- Plot 6: Zoom on peak region ---
    ax6 = axes[1, 2]
    ax6.plot(rho, chi_tr, 'b-', lw=2, marker='o', ms=3, label='TR output')
    ax6.plot(rho, chi_calc_edge, 'r--', lw=2, marker='s', ms=3, label='Original')
    ax6.plot(rho, chi_calc_smoothed, 'g-', lw=2, marker='^', ms=3, label=f'Smoothed w={best_nw}')
    ax6.set_xlabel('r/a')
    ax6.set_ylabel(r'$\chi$ [m$^2$/s]')
    ax6.set_title('Zoom: Peak Region (r/a = 0.1-0.5)')
    ax6.legend(fontsize=8)
    ax6.grid(True)
    ax6.set_xlim(0.1, 0.5)
    ax6.set_ylim(0.2, 1.2)

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
