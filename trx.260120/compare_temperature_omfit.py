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
R0 = 8.03        # Major radius [m]
A_MINOR = 2.68  # Minor radius [m]
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
            profiles['qdelt'] = -np.array(pwo['qdelt']) / 1e6  # e-i exchange rate
        if 'qdelt_i' in pwo:
            profiles['qdelt_i'] = np.array(pwo['qdelt_i']) / 1e6
        if 'qcx' in pwo:
            profiles['qcx'] = np.array(pwo['qcx']) / 1e6
        if 'qioni' in pwo:
            profiles['qioni'] = np.array(pwo['qioni']) / 1e6
        if 'qe2d' in pwo:
            profiles['qe2d'] = np.array(pwo['qe2d']) / 1e6
        if 'qmag' in pwo:
            profiles['qmag'] = np.array(pwo['qmag']) / 1e6
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

    # Read optional species densities (e.g., Ar) from ions_* blocks
    for key, ion_block in nml.items():
        if not str(key).lower().startswith('ions_'):
            continue
        if not isinstance(ion_block, dict):
            continue
        name = str(ion_block.get('name', '')).strip().lower()
        if 'density' not in ion_block:
            continue
        density_1e20 = np.array(ion_block['density']) / 1e20
        if name == 'ar':
            profiles['nAr'] = density_1e20
        elif name in ('he', 'he4', 'alpha'):
            profiles['nHe'] = density_1e20

    return profiles


def load_rho_to_ra_mapping(csv_file='~/CFEDRSW/OMFIT_out/rho_rmin.csv'):
    """Load rho to r/a mapping."""
    csv_file = os.path.expanduser(csv_file)
    if not os.path.exists(csv_file):
        # Try local file
        if os.path.exists('rho_rmin.csv'):
            csv_file = 'rho_rmin.csv'
        else:
            return None
    df = pd.read_csv(csv_file)
    return df['rho'].values, df['rho'].values


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


def read_chi_se_corediv(filename='Chi_Se_COREDIV.DAT', a=A_MINOR, R0=R0):
    """
    Read Chi_Se_COREDIV.DAT file from OMFIT/COREDIV.

    File format (space-separated, 1 header line):
        r[m]  ne[m^-3]  ni[m^-3]  chi_e[m²/s]  chi_i[m²/s]  Se[W/m³]  Si[W/m³]  G1  G2  G3

    G1 = Vprime / (4π²R * r), so Vprime = G1 * 4π²R * r
    G2 = <|∇r|²> (flux surface average)

    Returns:
        dict with keys: r, r_over_a, ne, ni, chi_e, chi_i, Se, Si, Vprime, G1, G2
    """
    data = []
    with open(filename, 'r') as f:
        f.readline()  # Skip header
        for line in f:
            parts = line.split()
            if len(parts) >= 10:
                data.append([float(x) for x in parts[:10]])
            elif len(parts) >= 7:
                # Old format without G1, G2, G3
                data.append([float(x) for x in parts[:7]] + [0, 0, 0])

    data = np.array(data)
    r = data[:, 0]
    G1 = data[:, 7]
    G2 = data[:, 8]

    # Calculate Vprime from G1: Vprime = G1 * 4π²R * r
    ATOR = 4.0 * np.pi**2 * R0
    Vprime = G1 * ATOR * r
    # Handle axis (r=0): use small value to avoid division by zero
    if Vprime[0] < 1e-10:
        Vprime[0] = Vprime[1] * 0.01

    return {
        'r': r,                             # r in meters
        'r_over_a': r / a,                  # r/a
        'ne': data[:, 1] / 1e20,            # convert to 10^20/m³
        'ni': data[:, 2] / 1e20,            # convert to 10^20/m³
        'chi_e': data[:, 3],                # m²/s
        'chi_i': data[:, 4],                # m²/s
        'Se': data[:, 5],                   # W/m³ (keep in W for proper solver)
        'Si': data[:, 6],                   # W/m³
        'Vprime': Vprime,                   # dV/dr [m²]
        'G1': G1,                           # Vprime/(4π²R*r)
        'G2': G2,                           # <|∇r|²>
    }


def solve_temperature_with_vprime(r, Vprime, ne, chi, Se, Te_bc, from_edge=True):
    """
    Solve steady-state heat equation using the proper flux-coordinate method.

    Heat equation in flux coordinates:
        1/V' * d/dr(V' * n * chi * dT/dr) = -S

    Algorithm (from ReproduceT.ipynb):
        1. Q(r) = integral from 0 to r of (V' * S * dr)  [heat flux, W]
        2. dT/dr = -Q / (e * V' * n * chi)              [temperature gradient]
        3. T(r) = integrate dT/dr from boundary

    Args:
        r: radius [m]
        Vprime: dV/dr [m²]
        ne: density [m⁻³]
        chi: thermal diffusivity [m²/s]
        Se: source [W/m³]
        Te_bc: boundary temperature [eV]
        from_edge: if True, integrate from edge; if False, from axis

    Returns:
        Te: temperature [eV]
        Q: heat flux [W]
    """
    E_CHARGE = 1.602176634e-19  # J/eV

    r = np.asarray(r)
    Vp = np.asarray(Vprime)
    ne = np.asarray(ne)
    chi = np.asarray(chi)
    Se = np.asarray(Se)

    N = len(r)

    # Ensure positive values to avoid numerical issues
    Vp = np.maximum(Vp, 1e-10)
    ne = np.maximum(ne, 1e10)  # At least 1e10 m^-3
    chi = np.maximum(chi, 0.01)

    # Step 1: Integrate source to get heat flux Q(r)
    # Q(r) = integral from 0 to r of (V' * S * dr)
    integrand = Vp * Se
    Q = np.zeros(N)
    for i in range(1, N):
        dr = r[i] - r[i-1]
        Q[i] = Q[i-1] + 0.5 * (integrand[i-1] + integrand[i]) * dr

    # Step 2: Temperature gradient from Fourier's law
    # dT/dr = -Q / (e * V' * n * chi)
    denom = E_CHARGE * Vp * ne * chi
    denom = np.maximum(denom, 1e-30)  # Avoid division by zero
    dTdr = -Q / denom

    # Step 3: Integrate to get temperature
    Te = np.zeros(N)

    if from_edge:
        # Integrate from edge to axis
        Te[-1] = Te_bc
        for i in range(N-2, -1, -1):
            dr = r[i+1] - r[i]
            grad_mid = 0.5 * (dTdr[i] + dTdr[i+1])
            Te[i] = Te[i+1] - grad_mid * dr
    else:
        # Integrate from axis to edge
        Te[0] = Te_bc
        for i in range(N-1):
            dr = r[i+1] - r[i]
            grad_mid = 0.5 * (dTdr[i] + dTdr[i+1])
            Te[i+1] = Te[i] + grad_mid * dr

    return Te, Q


def find_csv_by_title(title_keyword, search_dir='.'):
    """Find a TR CSV file by searching for a keyword in its title line."""
    import glob
    for f in sorted(glob.glob(os.path.join(search_dir, 'tr_data_*.csv'))):
        with open(f, 'r') as fh:
            title = fh.readline().strip()
        if title_keyword in title:
            return f
    return None


def main():
    print("="*70)
    print("TEMPERATURE COMPARISON: TR vs OMFIT RADIATION")
    print("="*70)

    # Auto-find CSV files by title content
    csv_dens     = find_csv_by_title('n(NS)')
    csv_temp     = find_csv_by_title('T(NS)')
    csv_chi_e    = find_csv_by_title('AKE,AKNCE')
    csv_chi_i    = find_csv_by_title('AKD,AKNCD')
    csv_pin      = find_csv_by_title('@PIN [MW')
    csv_power    = find_csv_by_title('POH,PNB,PNF')
    csv_rad      = find_csv_by_title('PRSUM,PRB,PRC')
    csv_pnb      = find_csv_by_title('PNBIN,PNBCL')
    csv_pnf      = find_csv_by_title('PNFIN,PNFCL')

    csv_map = {
        'n(NS)':        csv_dens,
        'T(NS)':        csv_temp,
        'AKE':          csv_chi_e,
        'AKD':          csv_chi_i,
        'PIN':          csv_pin,
        'POH/PNB/PNF':  csv_power,
        'PRSUM/PRB':    csv_rad,
        'PNBIN/PNBCL':  csv_pnb,
        'PNFIN/PNFCL':  csv_pnf,
    }

    missing = [k for k, v in csv_map.items() if v is None]
    if missing:
        raise FileNotFoundError(f"Could not find CSV files for: {missing}")

    print("\n   Auto-detected CSV files:")
    for k, v in csv_map.items():
        print(f"     {k:<16s} -> {os.path.basename(v)}")

    # Load TR data
    print("\n1. Loading TR simulation data...")
    r_tr, dens, _ = read_tr_csv(csv_dens)
    _, temp, _ = read_tr_csv(csv_temp)
    _, chi_e_data, _ = read_tr_csv(csv_chi_e)
    _, chi_i_data, _ = read_tr_csv(csv_chi_i)
    _, pin_data, _ = read_tr_csv(csv_pin)
    _, power_src, _ = read_tr_csv(csv_power)
    _, rad_exch, _ = read_tr_csv(csv_rad)
    _, pnb_split, _ = read_tr_csv(csv_pnb)
    _, pnf_split, _ = read_tr_csv(csv_pnf)

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
    PRSUM_tr = rad_exch.get('PRSUM', power_src.get('-PRSUM', np.zeros_like(r_tr)))
    PRB_tr = rad_exch.get('PRB', np.zeros_like(r_tr))
    PRC_tr = rad_exch.get('PRC', np.zeros_like(r_tr))
    PRL_tr = rad_exch.get('PRL', np.zeros_like(r_tr))
    PIE_tr = rad_exch.get('PIE', np.zeros_like(r_tr))
    PCX_tr = rad_exch.get('PCX', np.zeros_like(r_tr))
    QEI_tr = rad_exch.get('QEI', np.zeros_like(r_tr))

    # TR component-wise net source (statefile-like with explicit QEI diagnostic)
    POH_tr = power_src.get('POH', np.zeros_like(r_tr))
    # PRF column: try 'PRF' first, then 'PRF(1:NS)' for backward compat
    PRF_e_tr = power_src.get('PRF', power_src.get('PRF(1:NS)', np.zeros_like(r_tr)))
    PRF_D_tr = power_src.get('Y6', np.zeros_like(r_tr))
    PRF_T_tr = power_src.get('Y7', np.zeros_like(r_tr))
    PRF_He_tr = power_src.get('Y8', np.zeros_like(r_tr))

    # PNBCL column: try 'PNBCL' first, then 'PNBCL(1:NS)'
    PNBCL_e_tr = pnb_split.get('PNBCL', pnb_split.get('PNBCL(1:NS)', np.zeros_like(r_tr)))
    PNBCL_D_tr = pnb_split.get('Y3', np.zeros_like(r_tr))
    PNBCL_T_tr = pnb_split.get('Y4', np.zeros_like(r_tr))
    PNBCL_He_tr = pnb_split.get('Y5', np.zeros_like(r_tr))

    # PNFCL column: try 'PNFCL' first, then 'PNFCL(1:NS)'
    PNFCL_e_tr = pnf_split.get('PNFCL', pnf_split.get('PNFCL(1:NS)', np.zeros_like(r_tr)))
    PNFCL_D_tr = pnf_split.get('Y3', np.zeros_like(r_tr))
    PNFCL_T_tr = pnf_split.get('Y4', np.zeros_like(r_tr))
    PNFCL_He_tr = pnf_split.get('Y5', np.zeros_like(r_tr))

    Se_components_tr = PNBCL_e_tr + PNFCL_e_tr + PRF_e_tr + POH_tr - PRSUM_tr - PIE_tr
    Si_components_tr = (
        PNBCL_D_tr + PNBCL_T_tr + PNBCL_He_tr +
        PNFCL_D_tr + PNFCL_T_tr + PNFCL_He_tr +
        PRF_D_tr + PRF_T_tr + PRF_He_tr - PCX_tr
    )

    # Requested net-source replacement:
    # Se_net ≈ (PIN_e related components) + QEI
    # Si_net ≈ (PIN_i related components) - QEI
    Se_net_tr = Se_components_tr + QEI_tr
    Si_net_tr = Si_components_tr - QEI_tr

    print(f"   Net-source replacement check (max |PIN-Sum|):")
    print(f"     e: {np.max(np.abs(PIN_e_tr - Se_components_tr)):.4e} MW/m³")
    print(f"     i: {np.max(np.abs(PIN_i_tr - Si_components_tr)):.4e} MW/m³")

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
    ne_omfit = interpolate_to_tr_grid(omfit.get('ne', np.zeros_like(rho_omfit)),
                                       rho_omfit, r_tr, rho_rmin)
    ni_omfit = interpolate_to_tr_grid(omfit.get('ni', np.zeros_like(rho_omfit)),
                                       rho_omfit, r_tr, rho_rmin)
    nAr_omfit = interpolate_to_tr_grid(omfit.get('nAr', np.zeros_like(rho_omfit)),
                                        rho_omfit, r_tr, rho_rmin)
    nHe_omfit = interpolate_to_tr_grid(omfit.get('nHe', np.zeros_like(rho_omfit)),
                                        rho_omfit, r_tr, rho_rmin)
    qrad_omfit = interpolate_to_tr_grid(omfit.get('qrad', np.zeros_like(rho_omfit)),
                                         rho_omfit, r_tr, rho_rmin)
    qfuse_omfit = interpolate_to_tr_grid(omfit.get('qfuse', np.zeros_like(rho_omfit)),
                                          rho_omfit, r_tr, rho_rmin)
    qfusi_omfit = interpolate_to_tr_grid(omfit.get('qfusi', np.zeros_like(rho_omfit)),
                                          rho_omfit, r_tr, rho_rmin)
    qrfe_omfit = interpolate_to_tr_grid(omfit.get('qrfe', np.zeros_like(rho_omfit)),
                                         rho_omfit, r_tr, rho_rmin)
    qrfi_omfit = interpolate_to_tr_grid(omfit.get('qrfi', np.zeros_like(rho_omfit)),
                                         rho_omfit, r_tr, rho_rmin)
    qbeame_omfit = interpolate_to_tr_grid(omfit.get('qbeame', np.zeros_like(rho_omfit)),
                                           rho_omfit, r_tr, rho_rmin)
    qbeami_omfit = interpolate_to_tr_grid(omfit.get('qbeami', np.zeros_like(rho_omfit)),
                                           rho_omfit, r_tr, rho_rmin)
    qe2d_omfit = interpolate_to_tr_grid(omfit.get('qe2d', np.zeros_like(rho_omfit)),
                                         rho_omfit, r_tr, rho_rmin)
    qioni_omfit = interpolate_to_tr_grid(omfit.get('qioni', np.zeros_like(rho_omfit)),
                                          rho_omfit, r_tr, rho_rmin)
    qmag_omfit = interpolate_to_tr_grid(omfit.get('qmag', np.zeros_like(rho_omfit)),
                                         rho_omfit, r_tr, rho_rmin)
    qcx_omfit = interpolate_to_tr_grid(omfit.get('qcx', np.zeros_like(rho_omfit)),
                                        rho_omfit, r_tr, rho_rmin)
    qdelt_omfit = interpolate_to_tr_grid(omfit.get('qdelt', np.zeros_like(rho_omfit)),
                                          rho_omfit, r_tr, rho_rmin)
    qdelt_i_omfit = interpolate_to_tr_grid(omfit.get('qdelt_i', np.zeros_like(rho_omfit)),
                                            rho_omfit, r_tr, rho_rmin)

    # Save OMFIT profiles used in output plots on TR grid
    if rho_rmin is not None:
        rho_map, rmin_map = rho_rmin
        a_map = np.max(rmin_map)
        r_over_a_map = rmin_map / a_map
        rho_on_tr = np.interp(r_tr, r_over_a_map, rho_map)
    else:
        rho_on_tr = r_tr.copy()

    # df_omfit_plot = pd.DataFrame({
    #     'r_m': r_tr * A_MINOR,
    #     'r_over_a': r_tr,
    #     'rho': rho_on_tr,
    #     'ne_omfit_1e20_m3': ne_omfit,
    #     'nAr_omfit_1e20_m3': nAr_omfit,
    #     'nHe_omfit_1e20_m3': nHe_omfit,
    #     'Te_omfit_keV': Te_omfit,
    # })
    df_omfit_plot = pd.DataFrame({
        'rho': rho_on_tr,
        'ne_omfit_1e20_m3': ne_omfit,
        'nAr_omfit_1e20_m3': nAr_omfit,
        'nHe_omfit_1e20_m3': nHe_omfit,
        'Te_omfit_keV': Te_tr,
    })
    omfit_plot_csv = 'omfit_r_rho_ne_te_from_temperature_plot.csv'
    df_omfit_plot.to_csv(omfit_plot_csv, index=False)
    print(f"   Saved OMFIT (r, rho, ne, Te) to {omfit_plot_csv}")

    # OMFIT net source terms (statefile-sum convention, see compare_Se_statefile.py)
    S_e_omfit = qfuse_omfit + qrfe_omfit + qbeame_omfit + qe2d_omfit + qmag_omfit + qdelt_omfit - qrad_omfit
    S_i_omfit = qfusi_omfit + qrfi_omfit + qbeami_omfit + qioni_omfit + qdelt_i_omfit - qcx_omfit

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
        corediv = read_chi_se_corediv('Chi_Se_COREDIV_rho.DAT')
        print(f"   COREDIV grid: {len(corediv['r'])} points, r = [0, {corediv['r'][-1]:.2f}] m")
        print(f"   ne_center = {corediv['ne'][0]:.4f} x10^20/m³")
        print(f"   chi_e_center = {corediv['chi_e'][0]:.4f} m²/s")
        print(f"   Se_center = {corediv['Se'][0]/1e6:.4f} MW/m³")
        print(f"   G1_center = {corediv['G1'][1]:.4f} (Vprime normalization)")
        print(f"   Vprime_center = {corediv['Vprime'][1]:.2f} m²")

        # Interpolate COREDIV to TR grid (use r in meters, not r/a)
        r_cd = corediv['r']
        r_tr_m = r_tr * A_MINOR  # Convert r/a to meters

        f_ne_cd = interp1d(r_cd, corediv['ne'], bounds_error=False, fill_value='extrapolate')
        f_ni_cd = interp1d(r_cd, corediv['ni'], bounds_error=False, fill_value='extrapolate')
        f_chi_e_cd = interp1d(r_cd, corediv['chi_e'], bounds_error=False, fill_value='extrapolate')
        f_chi_i_cd = interp1d(r_cd, corediv['chi_i'], bounds_error=False, fill_value='extrapolate')
        f_Se_cd = interp1d(r_cd, corediv['Se'], bounds_error=False, fill_value='extrapolate')
        f_Si_cd = interp1d(r_cd, corediv['Si'], bounds_error=False, fill_value='extrapolate')
        f_Vprime_cd = interp1d(r_cd, corediv['Vprime'], bounds_error=False, fill_value='extrapolate')

        ne_corediv = f_ne_cd(r_tr_m)      # 10^20/m³
        ni_corediv = f_ni_cd(r_tr_m)      # 10^20/m³
        chi_e_corediv = f_chi_e_cd(r_tr_m)  # m²/s
        chi_i_corediv = f_chi_i_cd(r_tr_m)  # m²/s
        Se_corediv = f_Se_cd(r_tr_m)      # W/m³
        Si_corediv = f_Si_cd(r_tr_m)      # W/m³
        Vprime_corediv = f_Vprime_cd(r_tr_m)  # m²
        # Handle axis
        if Vprime_corediv[0] < 1e-10:
            Vprime_corediv[0] = Vprime_corediv[1] * 0.01
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

    # Case 1: TR net source (components ± QEI), optimal chi
    chi_factor_e_tr = find_chi_factor(ne_tr, chi_e_tr, Se_net_tr, Te_edge, Te_tr)
    chi_factor_i_tr = find_chi_factor(ni_tr, chi_i_tr, Si_net_tr, Ti_edge, Ti_tr)

    Te_calc_tr = solve_temperature_matrix(r_tr, ne_tr, chi_e_corediv, Se_net_tr, Te_edge,
                                            chi_factor=chi_factor_e_tr)
    Ti_calc_tr = solve_temperature_matrix(r_tr, ni_tr, chi_e_corediv, Si_net_tr, Ti_edge,
                                            chi_factor=chi_factor_i_tr)

    print(f"\n   Case 1: TR net source (PIN-related ± QEI), optimal chi")
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
    # Se_net_tr = heating - PRSUM - PIE + QEI  -> heating = Se_net_tr + PRSUM + PIE - QEI
    # Replace TR radiation with OMFIT radiation while keeping other terms:
    # S_e_hybrid = heating - qrad_omfit - PIE + QEI = Se_net_tr + PRSUM - qrad_omfit
    S_e_hybrid = Se_net_tr + PRSUM_tr - qrad_omfit
    S_i_hybrid = Si_net_tr  # Keep TR ion net source

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

    # Case 5: COREDIV data with proper Vprime solver (from ReproduceT.ipynb)
    Te_calc_corediv = None
    Ti_calc_corediv = None
    chi_factor_e_corediv = None
    chi_factor_i_corediv = None
    Te_calc_corediv_direct = None
    Ti_calc_corediv_direct = None
    Te_calc_corediv_trchi = None
    Ti_calc_corediv_trchi = None

    if corediv is not None:
        # Use the proper flux-coordinate solver with Vprime
        # Convert units: ne from 10^20/m³ to m^-3, Te_edge from keV to eV
        Te_edge_eV = Te_edge * 1000.0  # keV -> eV
        Ti_edge_eV = Ti_edge * 1000.0

        # Solve with chi_factor=1.0 (no adjustment)
        Te_direct_eV, Qe_direct = solve_temperature_with_vprime(
            r_tr_m, Vprime_corediv, ne_corediv * 1e20, chi_e_corediv,
            Se_corediv, 100, from_edge=True)
        Ti_direct_eV, Qi_direct = solve_temperature_with_vprime(
            r_tr_m, Vprime_corediv, ni_corediv * 1e20, chi_i_corediv,
            Si_corediv, 100, from_edge=True)

        Te_calc_corediv_direct = Te_direct_eV / 1000.0  # eV -> keV
        Ti_calc_corediv_direct = Ti_direct_eV / 1000.0

        print(f"\n   Case 5a: COREDIV with Vprime solver (chi_factor=1.0)")
        print(f"   Te_center = {Te_calc_corediv_direct[0]:.2f} keV (target OMFIT: {Te_omfit[0]:.2f})")
        print(f"   Ti_center = {Ti_calc_corediv_direct[0]:.2f} keV (target OMFIT: {Ti_omfit[0]:.2f})")
        print(f"   Q_e(edge) = {Qe_direct[-1]/1e6:.2f} MW, Q_i(edge) = {Qi_direct[-1]/1e6:.2f} MW")

        # Also try with chi scaling to match OMFIT temperature
        def find_chi_factor_vprime(chi, n, Se, Vprime, T_edge_eV, T_target_keV):
            def objective(factor):
                Te_eV, _ = solve_temperature_with_vprime(
                    r_tr_m, Vprime, n * 1e20, chi * factor, Se, T_edge_eV, from_edge=True)
                return (Te_eV[0]/1000.0 - T_target_keV[0])**2
            result = minimize_scalar(objective, bounds=(0.1, 5.0), method='bounded')
            return result.x

        chi_factor_e_corediv = find_chi_factor_vprime(
            chi_e_corediv, ne_corediv, Se_corediv, Vprime_corediv, Te_edge_eV, Te_omfit)
        chi_factor_i_corediv = find_chi_factor_vprime(
            chi_i_corediv, ni_corediv, Si_corediv, Vprime_corediv, Ti_edge_eV, Ti_omfit)

        Te_scaled_eV, _ = solve_temperature_with_vprime(
            r_tr_m, Vprime_corediv, ne_corediv * 1e20, chi_e_corediv * chi_factor_e_corediv,
            Se_corediv, Te_edge_eV, from_edge=True)
        Ti_scaled_eV, _ = solve_temperature_with_vprime(
            r_tr_m, Vprime_corediv, ni_corediv * 1e20, chi_i_corediv * chi_factor_i_corediv,
            Si_corediv, Ti_edge_eV, from_edge=True)

        Te_calc_corediv = Te_scaled_eV / 1000.0
        Ti_calc_corediv = Ti_scaled_eV / 1000.0

        print(f"\n   Case 5b: COREDIV with Vprime, optimal chi for OMFIT target")
        print(f"   chi_factor_e = {chi_factor_e_corediv:.4f}, chi_factor_i = {chi_factor_i_corediv:.4f}")
        print(f"   Te_center = {Te_calc_corediv[0]:.2f} keV (target: {Te_omfit[0]:.2f})")
        print(f"   Ti_center = {Ti_calc_corediv[0]:.2f} keV (target: {Ti_omfit[0]:.2f})")

        # Case 5c: COREDIV source + Vprime, but TR chi (no adjustment)
        Te_corediv_trchi_eV, Qe_trchi = solve_temperature_with_vprime(
            r_tr_m, Vprime_corediv, ne_corediv * 1e20, chi_e_tr,
            Se_corediv, 100, from_edge=True)
        Ti_corediv_trchi_eV, Qi_trchi = solve_temperature_with_vprime(
            r_tr_m, Vprime_corediv, ni_corediv * 1e20, chi_i_tr,
            Si_corediv, 100, from_edge=True)

        Te_calc_corediv_trchi = Te_corediv_trchi_eV / 1000.0  # eV -> keV
        Ti_calc_corediv_trchi = Ti_corediv_trchi_eV / 1000.0

        print(f"\n   Case 5c: COREDIV source + Vprime, TR chi (no adjustment)")
        print(f"   Te_center = {Te_calc_corediv_trchi[0]:.2f} keV (target OMFIT: {Te_omfit[0]:.2f})")
        print(f"   Ti_center = {Ti_calc_corediv_trchi[0]:.2f} keV (target OMFIT: {Ti_omfit[0]:.2f})")
        print(f"   Q_e(edge) = {Qe_trchi[-1]/1e6:.2f} MW, Q_i(edge) = {Qi_trchi[-1]/1e6:.2f} MW")

    # Plot comparison
    print("\n5. Plotting results...")
    fig, axes = plt.subplots(2, 4, figsize=(20, 10))

    # Electron temperature comparison
    ax = axes[0, 0]
    ax.plot(r_tr, Te_tr, 'b-', lw=2, label='TR simulation')
    ax.plot(r_tr, Te_omfit, 'r-', lw=2, label='OMFIT target')
    # ax.plot(r_tr, Te_calc_tr, 'g--', lw=1.5, label='Calc (TR source+OMFIT source')
    # ax.plot(r_tr, Te_calc_omfit_same, 'm--', lw=1.5, label='Calc (OMFIT source)')
    # if Te_calc_corediv_direct is not None:
        # ax.plot(r_tr, Te_calc_corediv_direct, 'g-', lw=2, label='OMFIT Source+Transport')
    # if Te_calc_corediv_trchi is not None:
        # ax.plot(r_tr, Te_calc_corediv_trchi, 'orange', lw=2, ls='-.', label='OMFIT Source + TR Transport')
    ax.set_xlabel('rho')
    ax.set_ylabel('Te [keV]')
    ax.set_title('Electron Temperature')
    ax.legend(fontsize=7, loc='upper right')
    ax.grid(True)
    ax.set_xlim([0, 1])

    # Ion temperature comparison
    ax = axes[0, 1]
    ax.plot(r_tr, Ti_tr, 'b-', lw=2, label='TR simulation')
    ax.plot(r_tr, Ti_omfit, 'r-', lw=2, label='OMFIT target')
    # ax.plot(r_tr, Ti_calc_tr, 'g--', lw=1.5, label='Calc (TR source)')
    # ax.plot(r_tr, Ti_calc_omfit_same, 'm--', lw=1.5, label='OMFIT')
    # if Ti_calc_corediv_direct is not None:
        # ax.plot(r_tr, Ti_calc_corediv_direct, 'g-', lw=2, label='OMFIT Source+Transport')
    # if Ti_calc_corediv_trchi is not None:
        # ax.plot(r_tr, Ti_calc_corediv_trchi, 'orange', lw=2, ls='-.', label='OMFIT Source + TR Transport')
    ax.set_xlabel('rho')
    ax.set_ylabel('Ti [keV]')
    ax.set_title('Ion Temperature')
    ax.legend(fontsize=8)
    ax.grid(True)
    ax.set_xlim([0, 1])

    # Radiation comparison
    ax = axes[0, 2]
    ax.plot(r_tr, PRSUM_tr, color='k', lw=2, label='TR PRSUM total')
    ax.plot(r_tr, PRB_tr, color='tab:blue', ls='--', lw=1.8, label='TR PRB')
    ax.plot(r_tr, PRC_tr, color='tab:green', ls='-.', lw=1.8, label='TR PRC')
    ax.plot(r_tr, PRL_tr, color='tab:purple', ls=':', lw=2.0, label='TR PRL')
    ax.plot(r_tr, qrad_omfit, 'r-', lw=2, label='OMFIT (qrad)')
    ax.set_xlabel('rho')
    ax.set_ylabel('P_rad [MW/m³]')
    ax.set_title('Radiation Profile (TR split vs OMFIT)')
    ax.legend(fontsize=8)
    ax.grid(True)
    ax.set_xlim([0, 1])

    # Electron source comparison
    ax = axes[1, 0]
    ax.plot(r_tr, Se_net_tr, 'b-', lw=2, label='TR (Se_net ≈ PIN_e-components + QEI)')
    ax.plot(r_tr, PIN_e_tr, 'b:', lw=1.5, label='TR PIN_e (reference)')
    # ax.plot(r_tr, S_e_omfit, 'r-', lw=2, label='OMFIT (net e)')
    # ax.plot(r_tr, S_e_hybrid, 'g--', lw=1.5, label='Hybrid')
    sources=pd.read_csv('rho_rmin_Si_Se_qe_qi_sum.csv')
    if corediv is not None:
        # ax.plot(r_tr, Se_corediv/1e6, 'r-.', lw=2, label='OMFIT (Se)')
        ax.plot(sources['rho'], sources['qe_sum']/1e6,'r-.', lw=2, label='OMFIT (Se)')
    ax.set_xlabel('rho')
    ax.set_ylabel('S_e [MW/m³]')
    ax.set_title('Electron Net Source')
    ax.legend(fontsize=8)
    ax.grid(True)
    ax.set_xlim([0, 1])

    # Ion source comparison
    ax = axes[1, 1]
    ax.plot(r_tr, Si_net_tr, 'b-', lw=2, label='TR (Si_net ≈ PIN_i-components - QEI)')
    ax.plot(r_tr, PIN_i_tr, 'b:', lw=1.5, label='TR PIN_i (reference)')
    # ax.plot(r_tr, S_i_omfit, 'r-', lw=2, label='OMFIT (net i)')
    if corediv is not None:
        # ax.plot(r_tr, Si_corediv/1e6, 'r-.', lw=2, label='OMFIT (Si)')
        ax.plot(sources['rho'], sources['qi_sum']/1e6, 'r-.', lw=2, label='OMFIT (Si)')
    ax.set_xlabel('rho')
    ax.set_ylabel('S_i [MW/m³]')
    ax.set_title('Ion Net Source')
    ax.legend(fontsize=8)
    ax.grid(True)
    ax.set_xlim([0, 1])
    ax.set_ylim([-0.2, 1])

    # Density comparison (TR vs OMFIT)
    ax = axes[1, 2]
    ax.plot(r_tr, ne_tr, 'b-', lw=2, label='TR n_e')
    ax.plot(r_tr, ne_omfit, 'b--', lw=2, label='OMFIT n_e')
    ax.plot(r_tr, ni_tr, 'r-', lw=2, label='TR n_i')
    ax.plot(r_tr, 2*ni_omfit, 'r--', lw=2, label='OMFIT n_i')
    ax.set_xlabel('rho')
    ax.set_ylabel('n [10²⁰/m³]')
    ax.set_title('Density Profile (TR vs OMFIT)')
    ax.legend(fontsize=8)
    ax.grid(True)
    ax.set_xlim([0, 1])

    # Electron transport coefficient comparison
    ax = axes[0, 3]
    ax.plot(r_tr, chi_e_tr, 'b-', lw=2, label='TR chi_e')
    print(chi_e_tr)
    if corediv is not None:
        ax.plot(r_tr, chi_e_corediv, 'r-', lw=2, label='OMFIT chi_e')
    ax.set_xlabel('rho')
    ax.set_ylabel('chi_e [m²/s]')
    ax.set_title('Electron Transport Coefficient')
    ax.legend()
    ax.grid(True)
    # ax.set_xlim([0, 1])
    ax.set_yscale('log')

    # Ion transport coefficient comparison
    ax = axes[1, 3]
    ax.plot(r_tr, chi_i_tr, 'b-', lw=2, label='TR chi_i')
    if corediv is not None:
        ax.plot(r_tr, chi_i_corediv, 'r-', lw=2, label=' chi_i')
    ax.set_xlabel('rho')
    ax.set_ylabel('chi_i [m²/s]')
    ax.set_title('Ion Transport Coefficient')
    ax.legend()
    ax.grid(True)
    ax.set_xlim([0, 1])

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
    if Te_calc_corediv_direct is not None:
        print(f"{'COREDIV Se+chi (Vprime)':<35} {Te_calc_corediv_direct[0]:>10.2f} {Ti_calc_corediv_direct[0]:>10.2f} {'CD':>8} {'CD':>8}")
    if Te_calc_corediv_trchi is not None:
        print(f"{'COREDIV Se + TR chi (Vprime)':<35} {Te_calc_corediv_trchi[0]:>10.2f} {Ti_calc_corediv_trchi[0]:>10.2f} {'TR':>8} {'TR':>8}")
    if Te_calc_corediv is not None:
        print(f"{'COREDIV+Vprime (opt chi)':<35} {Te_calc_corediv[0]:>10.2f} {Ti_calc_corediv[0]:>10.2f} {chi_factor_e_corediv:>8.4f} {chi_factor_i_corediv:>8.4f}")

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
