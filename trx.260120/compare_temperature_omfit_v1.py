#!/usr/bin/env python3
"""
Compare temperature calculations using TR vs OMFIT radiation data.

This script:
1. Reads OMFIT target temperatures and power profiles
2. Reads TR simulation results (chi, density, temperature)
3. Calculates temperature using:
   a) TR's PIN source terms (original)
   b) OMFIT heating/radiation profiles
4. Compares calculated temperatures with both TR and OMFIT targets

The goal is to see if using OMFIT's radiation profile improves temperature agreement.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import f90nml
import os
from scipy.interpolate import interp1d
from scipy.linalg import solve_banded

# Physical parameters for CFEDR
R0 = 7.8        # Major radius [m]
A_MINOR = 2.44  # Minor radius [m]
KAPPA = 1.89    # Elongation

def read_tr_csv(filename):
    """Read TR output CSV file."""
    with open(filename, 'r') as f:
        title = f.readline().strip()
    df = pd.read_csv(filename, skiprows=1, index_col=False)
    df = df.loc[:, ~df.columns.str.contains('^Unnamed')]
    x = df.iloc[:, 0].values
    data = {col: df[col].values for col in df.columns[1:]}
    return x, data, title


def read_omfit_profiles(namelist_file='profiles_CFEDR.namelist'):
    """
    Read OMFIT profiles including temperature, heating, and radiation.

    Returns:
        dict: Profiles on rho grid with keys:
            - rho: sqrt(toroidal flux) grid
            - Te, Ti: temperature [keV]
            - qrad: radiation power density [MW/m³]
            - qfuse, qfusi: fusion heating to e/i [MW/m³]
            - qrfe: RF heating to electrons [MW/m³]
            - pow_ei: e-i exchange power [MW]
    """
    nml = f90nml.read(namelist_file)

    profiles = {}

    # Power profiles from ONETWO (contains rho grid)
    if 'powers_particle_flux_onetwo' in nml:
        pwo = nml['powers_particle_flux_onetwo']
        profiles['rho'] = np.array(pwo['rho'])
        # Heating rates are in W/m³, convert to MW/m³
        if 'qrad' in pwo:
            profiles['qrad'] = np.array(pwo['qrad']) / 1e6
        if 'qfuse' in pwo:
            profiles['qfuse'] = np.array(pwo['qfuse']) / 1e6
        if 'qfusi' in pwo:
            profiles['qfusi'] = np.array(pwo['qfusi']) / 1e6
        if 'qrfe' in pwo:
            profiles['qrfe'] = np.array(pwo['qrfe']) / 1e6
        if 'qrfi' in pwo:
            profiles['qrfi'] = np.array(pwo['qrfi']) / 1e6
        if 'qbeame' in pwo:
            profiles['qbeame'] = np.array(pwo['qbeame']) / 1e6
        if 'qbeami' in pwo:
            profiles['qbeami'] = np.array(pwo['qbeami']) / 1e6
        if 'qdelt' in pwo:
            profiles['qdelt'] = np.array(pwo['qdelt']) / 1e6  # e-i exchange rate
        # Power flux profiles (MW)
        if 'pow_ei' in pwo:
            profiles['pow_ei'] = np.array(pwo['pow_ei'])

    # Electron temperature (has 201 points, same as power profiles)
    if 'electron' in nml:
        profiles['Te'] = np.array(nml['electron']['temperature'])
        profiles['ne'] = np.array(nml['electron']['density']) / 1e20  # to 10^20/m³

    # Ion temperature (use D as main ion)
    if 'ions_1' in nml:
        profiles['Ti'] = np.array(nml['ions_1']['temperature'])
        profiles['ni'] = np.array(nml['ions_1']['density']) / 1e20

    return profiles


def load_rho_to_ra_mapping(csv_file='~/CFEDRSW/rho_rmin.csv'):
    """Load rho to r/a mapping."""
    csv_file = os.path.expanduser(csv_file)
    if not os.path.exists(csv_file):
        # Try local file
        if os.path.exists('rho_rmin.csv'):
            csv_file = 'rho_rmin.csv'
        else:
            return None
    df = pd.read_csv(csv_file)
    return df['rho'].values, df['rmin'].values


def interpolate_to_tr_grid(omfit_data, rho_omfit, r_tr, rho_rmin=None):
    """
    Interpolate OMFIT data to TR r/a grid.

    Args:
        omfit_data: data on rho grid
        rho_omfit: OMFIT rho grid
        r_tr: TR r/a grid
        rho_rmin: tuple of (rho, rmin) for coordinate mapping

    Returns:
        data interpolated to r_tr grid
    """
    if rho_rmin is not None:
        rho_map, rmin_map = rho_rmin
        a = np.max(rmin_map)
        r_over_a_map = rmin_map / a
        # Convert rho_omfit to r/a
        r_over_a_omfit = np.interp(rho_omfit, rho_map, r_over_a_map)
    else:
        # Approximate: assume rho ≈ r/a (valid for circular plasmas)
        r_over_a_omfit = rho_omfit

    # Interpolate to TR grid
    f = interp1d(r_over_a_omfit, omfit_data, kind='linear',
                 bounds_error=False, fill_value='extrapolate')
    return f(r_tr)


def solve_temperature_matrix(r, n, chi, S_net, T_edge, chi_factor=1.0):
    """
    Solve steady-state heat equation using finite volume method.

    1/r * d/dr(r * n * chi * dT/dr) = -S_net

    Args:
        r: normalized radius (r/a)
        n: density [10^20/m³]
        chi: thermal diffusivity [m²/s]
        S_net: net source [MW/m³]
        T_edge: boundary temperature [keV]
        chi_factor: scaling factor for chi

    Returns:
        T: temperature profile [keV]
    """
    nr = len(r)
    a = A_MINOR
    r_real = r * a
    dr = r_real[1] - r_real[0] if nr > 1 else a / 50

    # Apply chi factor with minimum value to avoid singular matrix
    chi_eff = np.maximum(chi * chi_factor, 0.01)

    # Half-point radii
    r_half = np.zeros(nr + 1)
    r_half[1:-1] = 0.5 * (r_real[:-1] + r_real[1:])
    r_half[0] = 0.5 * r_real[0]  # Small non-zero value at center
    r_half[-1] = r_real[-1] + 0.5 * dr

    # Half-point values
    n_half = np.zeros(nr + 1)
    chi_half = np.zeros(nr + 1)
    n_half[1:-1] = 0.5 * (n[:-1] + n[1:])
    chi_half[1:-1] = 0.5 * (chi_eff[:-1] + chi_eff[1:])
    n_half[0] = n[0]
    chi_half[0] = chi_eff[0]
    n_half[-1] = n[-1]
    chi_half[-1] = chi_eff[-1]

    # Ensure positive values
    n_half = np.maximum(n_half, 0.01)
    chi_half = np.maximum(chi_half, 0.01)

    # Diffusion coefficient at half-points
    D_half = r_half * n_half * chi_half / dr

    # Source term conversion: MW/m³ -> keV·10^20/m³/s
    S_conv = S_net * 62.4 / 1.5

    # Volume element
    dV = np.maximum(r_real * dr, 1e-6)

    # Build tridiagonal matrix
    A_coef = np.zeros(nr)
    B_coef = np.zeros(nr)
    C_coef = np.zeros(nr)
    D_rhs = np.zeros(nr)

    for i in range(1, nr-1):
        A_coef[i] = D_half[i]
        C_coef[i] = D_half[i+1]
        B_coef[i] = -(A_coef[i] + C_coef[i])
        D_rhs[i] = -S_conv[i] * dV[i]

    # Center: dT/dr = 0 (use small regularization)
    C_coef[0] = D_half[1]
    B_coef[0] = -D_half[1] - 1e-10
    D_rhs[0] = -S_conv[0] * dV[0]

    # Edge: Dirichlet
    B_coef[-1] = 1.0
    A_coef[-1] = 0.0
    D_rhs[-1] = T_edge

    # Solve
    ab = np.zeros((3, nr))
    ab[0, 1:] = C_coef[:-1]
    ab[1, :] = B_coef
    ab[2, :-1] = A_coef[1:]

    try:
        T = solve_banded((1, 1), ab, D_rhs)
    except np.linalg.LinAlgError:
        # Fallback: simple integration method
        T = np.full(nr, T_edge)
        for i in range(nr-2, -1, -1):
            denom = n_half[i+1] * chi_half[i+1] * r_half[i+1] / dr + 1e-10
            dT = S_conv[i] * dV[i] / denom
            T[i] = T[i+1] + dT * dr

    return T


def integrate_power(r, p, R0=R0, a=A_MINOR):
    """Integrate power density over plasma volume [MW]."""
    r_real = r * a
    dV = 2 * np.pi * R0 * 2 * np.pi * r_real
    return np.trapz(p * dV, r_real)


def read_chi_se_corediv(filename='Chi_Se_COREDIV.DAT', a=A_MINOR):
    """
    Read Chi_Se_COREDIV.DAT file from OMFIT/COREDIV.

    File format (space-separated, 1 header line):
        r[m]  ne[m^-3]  ni[m^-3]  chi_e[m²/s]  chi_i[m²/s]  Se[W/m³]  Si[W/m³]  G1  G2  G3

    Returns:
        dict with keys: r_over_a, ne, ni, chi_e, chi_i, Se, Si
    """
    data = []
    with open(filename, 'r') as f:
        f.readline()  # Skip header
        for line in f:
            parts = line.split()
            if len(parts) >= 7:
                data.append([float(x) for x in parts[:7]])

    data = np.array(data)

    return {
        'r': data[:, 0],                    # r in meters
        'r_over_a': data[:, 0] / a,         # r/a
        'ne': data[:, 1] / 1e20,            # convert to 10^20/m³
        'ni': data[:, 2] / 1e20,            # convert to 10^20/m³
        'chi_e': data[:, 3],                # m²/s
        'chi_i': data[:, 4],                # m²/s
        'Se': data[:, 5] / 1e6,             # convert W/m³ to MW/m³
        'Si': data[:, 6] / 1e6,             # convert W/m³ to MW/m³
    }


def main():
    print("="*70)
    print("TEMPERATURE COMPARISON: TR vs OMFIT RADIATION")
    print("="*70)

    # Load TR data
    print("\n1. Loading TR simulation data...")
    r_tr, dens, _ = read_tr_csv('tr_data_017.csv')
    _, temp, _ = read_tr_csv('tr_data_019.csv')
    _, chi_e_data, _ = read_tr_csv('tr_data_023.csv')
    _, chi_i_data, _ = read_tr_csv('tr_data_024.csv')
    _, pin_data, _ = read_tr_csv('tr_data_037.csv')
    _, power_loss, _ = read_tr_csv('tr_data_021.csv')

    # Column names may be lowercase or uppercase
    ne_tr = dens.get('nE', dens.get('NE', np.zeros_like(r_tr)))
    ni_tr = (dens.get('nD', dens.get('ND', np.zeros_like(r_tr))) +
             dens.get('nT', dens.get('NT', np.zeros_like(r_tr))))
    Te_tr = temp.get('TE', temp.get('Te', np.zeros_like(r_tr)))
    Ti_tr = temp.get('TD', temp.get('Ti', np.zeros_like(r_tr)))
    chi_e_tr = chi_e_data.get('AKE', np.zeros_like(r_tr))
    chi_i_tr = chi_i_data.get('AKD', np.zeros_like(r_tr))
    PIN_e_tr = pin_data.get('PIN_1', np.zeros_like(r_tr))
    PIN_i_tr = (pin_data.get('PIN_2', np.zeros_like(r_tr)) +
                pin_data.get('PIN_3', np.zeros_like(r_tr)) +
                pin_data.get('PIN_4', np.zeros_like(r_tr)))
    PRSUM_tr = power_loss.get('PRSUM', np.zeros_like(r_tr))
    PIE_tr = power_loss.get('PIE', np.zeros_like(r_tr))

    print(f"   TR grid: {len(r_tr)} points, r/a = [{r_tr[0]:.4f}, {r_tr[-1]:.4f}]")
    print(f"   Te_center (TR sim) = {Te_tr[0]:.2f} keV")
    print(f"   Ti_center (TR sim) = {Ti_tr[0]:.2f} keV")

    # Load OMFIT data
    print("\n2. Loading OMFIT profiles...")
    omfit = read_omfit_profiles('profiles_CFEDR.namelist')
    rho_omfit = omfit['rho']

    # Load coordinate mapping
    rho_rmin = load_rho_to_ra_mapping()

    print(f"   OMFIT grid: {len(rho_omfit)} points")
    print(f"   Te_center (OMFIT) = {omfit['Te'][0]:.2f} keV")
    print(f"   Ti_center (OMFIT) = {omfit['Ti'][0]:.2f} keV")

    # Interpolate OMFIT to TR grid
    print("\n3. Interpolating OMFIT to TR grid...")
    Te_omfit = interpolate_to_tr_grid(omfit['Te'], rho_omfit, r_tr, rho_rmin)
    Ti_omfit = interpolate_to_tr_grid(omfit['Ti'], rho_omfit, r_tr, rho_rmin)
    qrad_omfit = interpolate_to_tr_grid(omfit.get('qrad', np.zeros_like(rho_omfit)),
                                         rho_omfit, r_tr, rho_rmin)
    qfuse_omfit = interpolate_to_tr_grid(omfit.get('qfuse', np.zeros_like(rho_omfit)),
                                          rho_omfit, r_tr, rho_rmin)
    qfusi_omfit = interpolate_to_tr_grid(omfit.get('qfusi', np.zeros_like(rho_omfit)),
                                          rho_omfit, r_tr, rho_rmin)
    qrfe_omfit = interpolate_to_tr_grid(omfit.get('qrfe', np.zeros_like(rho_omfit)),
                                         rho_omfit, r_tr, rho_rmin)
    qdelt_omfit = interpolate_to_tr_grid(omfit.get('qdelt', np.zeros_like(rho_omfit)),
                                          rho_omfit, r_tr, rho_rmin)

    # OMFIT net source terms
    # Electron: fusion + RF - radiation - e-i exchange
    S_e_omfit = qfuse_omfit + qrfe_omfit - qrad_omfit - qdelt_omfit
    # Ion: fusion + e-i exchange
    S_i_omfit = qfusi_omfit + qdelt_omfit

    print(f"   qrad_center (OMFIT) = {qrad_omfit[0]:.4f} MW/m³")
    print(f"   PRSUM_center (TR)   = {PRSUM_tr[0]:.4f} MW/m³")
    print(f"   Radiation ratio: OMFIT/TR = {qrad_omfit[0]/PRSUM_tr[0]:.2f}" if PRSUM_tr[0] > 0 else "")

    # Integrated powers
    P_rad_omfit = integrate_power(r_tr, qrad_omfit)
    P_rad_tr = integrate_power(r_tr, PRSUM_tr)
    print(f"\n   Integrated radiation: OMFIT={P_rad_omfit:.1f} MW, TR={P_rad_tr:.1f} MW")

    # Load COREDIV data (Chi_Se_COREDIV.DAT)
    print("\n3b. Loading COREDIV data (Chi_Se_COREDIV.DAT)...")
    corediv = None
    try:
        corediv = read_chi_se_corediv('Chi_Se_COREDIV.DAT')
        print(f"   COREDIV grid: {len(corediv['r'])} points, r = [0, {corediv['r'][-1]:.2f}] m")
        print(f"   ne_center = {corediv['ne'][0]:.4f} x10^20/m³")
        print(f"   chi_e_center = {corediv['chi_e'][0]:.4f} m²/s")
        print(f"   Se_center = {corediv['Se'][0]:.4f} MW/m³")

        # Interpolate COREDIV to TR grid
        f_ne_cd = interp1d(corediv['r_over_a'], corediv['ne'], bounds_error=False, fill_value='extrapolate')
        f_ni_cd = interp1d(corediv['r_over_a'], corediv['ni'], bounds_error=False, fill_value='extrapolate')
        f_chi_e_cd = interp1d(corediv['r_over_a'], corediv['chi_e'], bounds_error=False, fill_value='extrapolate')
        f_chi_i_cd = interp1d(corediv['r_over_a'], corediv['chi_i'], bounds_error=False, fill_value='extrapolate')
        f_Se_cd = interp1d(corediv['r_over_a'], corediv['Se'], bounds_error=False, fill_value='extrapolate')
        f_Si_cd = interp1d(corediv['r_over_a'], corediv['Si'], bounds_error=False, fill_value='extrapolate')

        ne_corediv = f_ne_cd(r_tr)
        ni_corediv = f_ni_cd(r_tr)
        chi_e_corediv = f_chi_e_cd(r_tr)
        chi_i_corediv = f_chi_i_cd(r_tr)
        Se_corediv = f_Se_cd(r_tr)
        Si_corediv = f_Si_cd(r_tr)
        print(f"   Interpolated to TR grid ({len(r_tr)} points)")
    except Exception as e:
        print(f"   Warning: Could not load COREDIV data: {e}")

    # Calculate temperatures with different source terms
    print("\n4. Calculating temperatures...")
    Te_edge = Te_tr[-1]
    Ti_edge = Ti_tr[-1]

    # Find optimal chi factors for TR case
    from scipy.optimize import minimize_scalar

    def find_chi_factor(n, chi, S_net, T_edge, T_target):
        def objective(factor):
            T = solve_temperature_matrix(r_tr, n, chi, S_net, T_edge, chi_factor=factor)
            return (T[0] - T_target[0])**2
        result = minimize_scalar(objective, bounds=(0.1, 2.0), method='bounded')
        return result.x

    # Case 1: TR source (PIN), optimal chi
    chi_factor_e_tr = find_chi_factor(ne_tr, chi_e_tr, PIN_e_tr, Te_edge, Te_tr)
    chi_factor_i_tr = find_chi_factor(ni_tr, chi_i_tr, PIN_i_tr, Ti_edge, Ti_tr)

    Te_calc_tr = solve_temperature_matrix(r_tr, ne_tr, chi_e_tr, PIN_e_tr, Te_edge,
                                           chi_factor=chi_factor_e_tr)
    Ti_calc_tr = solve_temperature_matrix(r_tr, ni_tr, chi_i_tr, PIN_i_tr, Ti_edge,
                                           chi_factor=chi_factor_i_tr)

    print(f"\n   Case 1: TR source (PIN), optimal chi")
    print(f"   chi_factor_e = {chi_factor_e_tr:.4f}, chi_factor_i = {chi_factor_i_tr:.4f}")
    print(f"   Te_center = {Te_calc_tr[0]:.2f} keV (target TR: {Te_tr[0]:.2f})")
    print(f"   Ti_center = {Ti_calc_tr[0]:.2f} keV (target TR: {Ti_tr[0]:.2f})")

    # Case 2: OMFIT source, same chi factors as TR
    Te_calc_omfit_same = solve_temperature_matrix(r_tr, ne_tr, chi_e_tr, S_e_omfit, Te_edge,
                                                   chi_factor=chi_factor_e_tr)
    Ti_calc_omfit_same = solve_temperature_matrix(r_tr, ni_tr, chi_i_tr, S_i_omfit, Ti_edge,
                                                   chi_factor=chi_factor_i_tr)

    print(f"\n   Case 2: OMFIT source, TR chi factors")
    print(f"   Te_center = {Te_calc_omfit_same[0]:.2f} keV (target OMFIT: {Te_omfit[0]:.2f})")
    print(f"   Ti_center = {Ti_calc_omfit_same[0]:.2f} keV (target OMFIT: {Ti_omfit[0]:.2f})")

    # Case 3: OMFIT source, optimal chi to match OMFIT temperature
    chi_factor_e_omfit = find_chi_factor(ne_tr, chi_e_tr, S_e_omfit, Te_edge, Te_omfit)
    chi_factor_i_omfit = find_chi_factor(ni_tr, chi_i_tr, S_i_omfit, Ti_edge, Ti_omfit)

    Te_calc_omfit_opt = solve_temperature_matrix(r_tr, ne_tr, chi_e_tr, S_e_omfit, Te_edge,
                                                  chi_factor=chi_factor_e_omfit)
    Ti_calc_omfit_opt = solve_temperature_matrix(r_tr, ni_tr, chi_i_tr, S_i_omfit, Ti_edge,
                                                  chi_factor=chi_factor_i_omfit)

    print(f"\n   Case 3: OMFIT source, optimal chi for OMFIT target")
    print(f"   chi_factor_e = {chi_factor_e_omfit:.4f}, chi_factor_i = {chi_factor_i_omfit:.4f}")
    print(f"   Te_center = {Te_calc_omfit_opt[0]:.2f} keV (target: {Te_omfit[0]:.2f})")
    print(f"   Ti_center = {Ti_calc_omfit_opt[0]:.2f} keV (target: {Ti_omfit[0]:.2f})")

    # Case 4: TR heating + OMFIT radiation
    # Reconstruct electron source: PIN_e already has TR losses subtracted
    # Add back TR losses and subtract OMFIT radiation
    # PIN_e = heating - PRSUM - PIE, so heating = PIN_e + PRSUM + PIE
    S_e_hybrid = PIN_e_tr + PRSUM_tr - qrad_omfit  # Use OMFIT radiation instead
    S_i_hybrid = PIN_i_tr  # Keep TR ion source

    chi_factor_e_hybrid = find_chi_factor(ne_tr, chi_e_tr, S_e_hybrid, Te_edge, Te_omfit)
    chi_factor_i_hybrid = find_chi_factor(ni_tr, chi_i_tr, S_i_hybrid, Ti_edge, Ti_omfit)

    Te_calc_hybrid = solve_temperature_matrix(r_tr, ne_tr, chi_e_tr, S_e_hybrid, Te_edge,
                                               chi_factor=chi_factor_e_hybrid)
    Ti_calc_hybrid = solve_temperature_matrix(r_tr, ni_tr, chi_i_tr, S_i_hybrid, Ti_edge,
                                               chi_factor=chi_factor_i_hybrid)

    print(f"\n   Case 4: TR heating + OMFIT radiation (hybrid)")
    print(f"   chi_factor_e = {chi_factor_e_hybrid:.4f}, chi_factor_i = {chi_factor_i_hybrid:.4f}")
    print(f"   Te_center = {Te_calc_hybrid[0]:.2f} keV (target OMFIT: {Te_omfit[0]:.2f})")
    print(f"   Ti_center = {Ti_calc_hybrid[0]:.2f} keV (target OMFIT: {Ti_omfit[0]:.2f})")

    # Case 5: COREDIV data (Se, ne, chi all from COREDIV file)
    Te_calc_corediv = None
    Ti_calc_corediv = None
    chi_factor_e_corediv = None
    chi_factor_i_corediv = None
    if corediv is not None:
        # Use COREDIV source (Se, Si), density (ne, ni), and chi directly
        # chi_factor=1.0 means use COREDIV chi as-is
        Te_calc_corediv_direct = solve_temperature_matrix(
            r_tr, ne_corediv, chi_e_corediv, Se_corediv, Te_edge, chi_factor=1.0)
        Ti_calc_corediv_direct = solve_temperature_matrix(
            r_tr, ni_corediv, chi_i_corediv, Si_corediv, Ti_edge, chi_factor=1.0)

        print(f"\n   Case 5a: COREDIV (Se, ne, chi) directly, chi_factor=1.0")
        print(f"   Te_center = {Te_calc_corediv_direct[0]:.2f} keV (target OMFIT: {Te_omfit[0]:.2f})")
        print(f"   Ti_center = {Ti_calc_corediv_direct[0]:.2f} keV (target OMFIT: {Ti_omfit[0]:.2f})")

        # Find optimal chi factor for COREDIV
        chi_factor_e_corediv = find_chi_factor(ne_corediv, chi_e_corediv, Se_corediv, Te_edge, Te_omfit)
        chi_factor_i_corediv = find_chi_factor(ni_corediv, chi_i_corediv, Si_corediv, Ti_edge, Ti_omfit)
        Te_calc_corediv = solve_temperature_matrix(
            r_tr, ne_corediv, chi_e_corediv, Se_corediv, Te_edge, chi_factor=chi_factor_e_corediv)
        Ti_calc_corediv = solve_temperature_matrix(
            r_tr, ni_corediv, chi_i_corediv, Si_corediv, Ti_edge, chi_factor=chi_factor_i_corediv)

        print(f"\n   Case 5b: COREDIV (Se, ne, chi), optimal chi for OMFIT target")
        print(f"   chi_factor_e = {chi_factor_e_corediv:.4f}, chi_factor_i = {chi_factor_i_corediv:.4f}")
        print(f"   Te_center = {Te_calc_corediv[0]:.2f} keV (target: {Te_omfit[0]:.2f})")
        print(f"   Ti_center = {Ti_calc_corediv[0]:.2f} keV (target: {Ti_omfit[0]:.2f})")

    # Plot comparison
    print("\n5. Plotting results...")
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))

    # Electron temperature comparison
    ax = axes[0, 0]
    ax.plot(r_tr, Te_tr, 'b-', lw=2, label='TR simulation')
    ax.plot(r_tr, Te_omfit, 'r-', lw=2, label='OMFIT target')
    ax.plot(r_tr, Te_calc_tr, 'g--', lw=1.5, label='Calc (TR source)')
    ax.plot(r_tr, Te_calc_omfit_same, 'm--', lw=1.5, label='Calc (OMFIT source)')
    if Te_calc_corediv is not None:
        ax.plot(r_tr, Te_calc_corediv, 'c-.', lw=2, label='Calc (COREDIV)')
    ax.set_xlabel('r/a')
    ax.set_ylabel('Te [keV]')
    ax.set_title('Electron Temperature')
    ax.legend(fontsize=8)
    ax.grid(True)
    ax.set_xlim([0, 1])

    # Ion temperature comparison
    ax = axes[0, 1]
    ax.plot(r_tr, Ti_tr, 'b-', lw=2, label='TR simulation')
    ax.plot(r_tr, Ti_omfit, 'r-', lw=2, label='OMFIT target')
    ax.plot(r_tr, Ti_calc_tr, 'g--', lw=1.5, label='Calc (TR source)')
    ax.plot(r_tr, Ti_calc_omfit_same, 'm--', lw=1.5, label='Calc (OMFIT source)')
    if Ti_calc_corediv is not None:
        ax.plot(r_tr, Ti_calc_corediv, 'c-.', lw=2, label='Calc (COREDIV)')
    ax.set_xlabel('r/a')
    ax.set_ylabel('Ti [keV]')
    ax.set_title('Ion Temperature')
    ax.legend(fontsize=8)
    ax.grid(True)
    ax.set_xlim([0, 1])

    # Radiation comparison
    ax = axes[0, 2]
    ax.plot(r_tr, PRSUM_tr, 'b-', lw=2, label='TR (PRSUM)')
    ax.plot(r_tr, qrad_omfit, 'r-', lw=2, label='OMFIT (qrad)')
    ax.set_xlabel('r/a')
    ax.set_ylabel('P_rad [MW/m³]')
    ax.set_title('Radiation Profile')
    ax.legend()
    ax.grid(True)
    ax.set_xlim([0, 1])

    # Electron source comparison
    ax = axes[1, 0]
    ax.plot(r_tr, PIN_e_tr, 'b-', lw=2, label='TR (PIN_e)')
    ax.plot(r_tr, S_e_omfit, 'r-', lw=2, label='OMFIT (net e)')
    ax.plot(r_tr, S_e_hybrid, 'g--', lw=1.5, label='Hybrid')
    if corediv is not None:
        ax.plot(r_tr, Se_corediv, 'c-.', lw=2, label='COREDIV (Se)')
    ax.set_xlabel('r/a')
    ax.set_ylabel('S_e [MW/m³]')
    ax.set_title('Electron Net Source')
    ax.legend(fontsize=8)
    ax.grid(True)
    ax.set_xlim([0, 1])

    # Ion source comparison
    ax = axes[1, 1]
    ax.plot(r_tr, PIN_i_tr, 'b-', lw=2, label='TR (PIN_i)')
    ax.plot(r_tr, S_i_omfit, 'r-', lw=2, label='OMFIT (net i)')
    if corediv is not None:
        ax.plot(r_tr, Si_corediv, 'c-.', lw=2, label='COREDIV (Si)')
    ax.set_xlabel('r/a')
    ax.set_ylabel('S_i [MW/m³]')
    ax.set_title('Ion Net Source')
    ax.legend(fontsize=8)
    ax.grid(True)
    ax.set_xlim([0, 1])

    # Chi factor comparison
    ax = axes[1, 2]
    if chi_factor_e_corediv is not None:
        cases = ['TR\n→TR', 'OMFIT\n→OMFIT', 'Hybrid\n→OMFIT', 'COREDIV\n→OMFIT']
        chi_e_vals = [chi_factor_e_tr, chi_factor_e_omfit, chi_factor_e_hybrid, chi_factor_e_corediv]
        chi_i_vals = [chi_factor_i_tr, chi_factor_i_omfit, chi_factor_i_hybrid, chi_factor_i_corediv]
    else:
        cases = ['TR source\n→TR target', 'OMFIT source\n→OMFIT target', 'Hybrid\n→OMFIT target']
        chi_e_vals = [chi_factor_e_tr, chi_factor_e_omfit, chi_factor_e_hybrid]
        chi_i_vals = [chi_factor_i_tr, chi_factor_i_omfit, chi_factor_i_hybrid]
    x = np.arange(len(cases))
    width = 0.35
    ax.bar(x - width/2, chi_e_vals, width, label='chi_e factor', color='blue', alpha=0.7)
    ax.bar(x + width/2, chi_i_vals, width, label='chi_i factor', color='red', alpha=0.7)
    ax.set_ylabel('Chi factor')
    ax.set_title('Optimal Chi Factors')
    ax.set_xticks(x)
    ax.set_xticklabels(cases, fontsize=8)
    ax.legend()
    ax.axhline(y=0.5, color='k', ls='--', lw=1, label='0.5 reference')
    ax.set_ylim([0, 1])
    ax.grid(True, axis='y')

    plt.tight_layout()
    plt.savefig('temperature_comparison_omfit.png', dpi=150)
    print("   Plot saved to temperature_comparison_omfit.png")

    # Summary table
    print("\n" + "="*70)
    print("SUMMARY TABLE")
    print("="*70)
    print(f"{'Case':<35} {'Te_center':>10} {'Ti_center':>10} {'chi_e':>8} {'chi_i':>8}")
    print("-"*70)
    print(f"{'TR simulation (reference)':<35} {Te_tr[0]:>10.2f} {Ti_tr[0]:>10.2f} {'-':>8} {'-':>8}")
    print(f"{'OMFIT target':<35} {Te_omfit[0]:>10.2f} {Ti_omfit[0]:>10.2f} {'-':>8} {'-':>8}")
    print("-"*70)
    print(f"{'TR source, opt chi → TR':<35} {Te_calc_tr[0]:>10.2f} {Ti_calc_tr[0]:>10.2f} {chi_factor_e_tr:>8.4f} {chi_factor_i_tr:>8.4f}")
    print(f"{'OMFIT source, TR chi factors':<35} {Te_calc_omfit_same[0]:>10.2f} {Ti_calc_omfit_same[0]:>10.2f} {chi_factor_e_tr:>8.4f} {chi_factor_i_tr:>8.4f}")
    print(f"{'OMFIT source, opt chi → OMFIT':<35} {Te_calc_omfit_opt[0]:>10.2f} {Ti_calc_omfit_opt[0]:>10.2f} {chi_factor_e_omfit:>8.4f} {chi_factor_i_omfit:>8.4f}")
    print(f"{'Hybrid (TR heat + OMFIT rad)':<35} {Te_calc_hybrid[0]:>10.2f} {Ti_calc_hybrid[0]:>10.2f} {chi_factor_e_hybrid:>8.4f} {chi_factor_i_hybrid:>8.4f}")

    print("\n" + "="*70)
    print("ANALYSIS")
    print("="*70)

    # Compare temperatures
    Te_diff_tr_omfit = Te_tr[0] - Te_omfit[0]
    Ti_diff_tr_omfit = Ti_tr[0] - Ti_omfit[0]
    print(f"\nTemperature difference (TR - OMFIT):")
    print(f"   Te: {Te_diff_tr_omfit:+.2f} keV ({Te_diff_tr_omfit/Te_omfit[0]*100:+.1f}%)")
    print(f"   Ti: {Ti_diff_tr_omfit:+.2f} keV ({Ti_diff_tr_omfit/Ti_omfit[0]*100:+.1f}%)")

    # Compare radiation
    P_rad_diff = P_rad_tr - P_rad_omfit
    print(f"\nRadiation difference (TR - OMFIT):")
    print(f"   P_rad: {P_rad_diff:+.1f} MW ({P_rad_diff/P_rad_omfit*100:+.1f}%)" if P_rad_omfit > 0 else "")

    # Conclusions
    print(f"\nConclusions:")
    if Te_tr[0] > Te_omfit[0]:
        if P_rad_tr < P_rad_omfit:
            print("   - TR temperature higher than OMFIT, radiation lower")
            print("   - Using OMFIT radiation would INCREASE losses → LOWER temperature")
            print("   - This should bring Te closer to OMFIT target")
        else:
            print("   - TR temperature higher and radiation also higher")
            print("   - Other heating sources differ between TR and OMFIT")
    else:
        print("   - TR temperature lower than OMFIT")
        print("   - May need to increase heating or reduce radiation")

    return {
        'r': r_tr,
        'Te_tr': Te_tr, 'Ti_tr': Ti_tr,
        'Te_omfit': Te_omfit, 'Ti_omfit': Ti_omfit,
        'Te_calc_tr': Te_calc_tr, 'Ti_calc_tr': Ti_calc_tr,
        'Te_calc_omfit': Te_calc_omfit_opt, 'Ti_calc_omfit': Ti_calc_omfit_opt,
        'chi_factor_e_tr': chi_factor_e_tr, 'chi_factor_i_tr': chi_factor_i_tr,
        'chi_factor_e_omfit': chi_factor_e_omfit, 'chi_factor_i_omfit': chi_factor_i_omfit
    }


if __name__ == '__main__':
    main()
