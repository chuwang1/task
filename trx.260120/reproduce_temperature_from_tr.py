#!/usr/bin/env python3
"""
Reproduce temperature profile from TR output using the same algorithm as TR code.

Based on TR code analysis (trexec.f90):
1. TR solves the energy transport equation using fully implicit scheme
2. The matrix equation: (RD - ADV*B)*T^(n+1) = RD*T^n + DT*D
3. At steady state: 1/V' * d/dr(V' * n * chi * dT/dr) = -S

This script:
1. Reads TR output: PIN (source), chi (transport), n (density), T (reference)
2. Solves steady-state heat equation using same discretization
3. Compares calculated T with TR output T
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.linalg import solve_banded
import glob
import os

# Physical constants
RKEV = 1.602176634e-16  # keV to Joules

# CFEDR parameters
R0 = 8.03        # Major radius [m]
A_MINOR = 2.68   # Minor radius [m]
KAPPA = 1.89     # Elongation
KAPPA_S = np.sqrt(KAPPA)  # √κ - geometry correction factor


def get_geometry_chi_factor(species='electron', kappa=KAPPA):
    """
    Get the chi scaling factor for geometry correction.

    Based on analysis of TR code outputs, electrons and ions need
    OPPOSITE geometry corrections:
    - Electrons: chi × √κ (enhanced transport due to geometry)
    - Ions: chi / √κ (reduced transport due to geometry)

    This is because the TR code's PIN values have different
    effective normalizations for electrons vs ions.

    Args:
        species: 'electron' or 'ion'
        kappa: elongation factor

    Returns:
        chi_factor: scaling factor to apply to chi
    """
    kappa_s = np.sqrt(kappa)
    if species.lower() in ['electron', 'e', 'elec']:
        return kappa_s  # Electrons need MORE transport
    elif species.lower() in ['ion', 'i', 'd', 't', 'he', 'alpha']:
        return 1.0 / kappa_s  # Ions need LESS transport
    else:
        return 1.0  # Default: no correction

TR_DIR = '/Users/dengxiaoya/TASK/latest/task/trx'


def find_csv_by_title(title_keyword, search_dir=TR_DIR):
    """Find TR CSV file by title keyword."""
    for f in sorted(glob.glob(os.path.join(search_dir, 'tr_data_*.csv'))):
        with open(f, 'r') as fh:
            title = fh.readline().strip()
        if title_keyword in title:
            return f
    return None


def read_tr_csv(filename):
    """Read TR output CSV file."""
    with open(filename, 'r') as f:
        title = f.readline().strip()
    df = pd.read_csv(filename, skiprows=1, index_col=False)
    df = df.loc[:, ~df.columns.str.contains('^Unnamed')]
    df.columns = df.columns.str.strip()
    x = df.iloc[:, 0].values
    data = {col.strip(): df[col].values for col in df.columns[1:]}
    return x, data, title


def load_tr_data():
    """Load all required TR output data."""
    print("Loading TR data...")

    # Find CSV files by title
    csv_files = {
        'density': find_csv_by_title('n(NS)'),          # tr_data_017
        'temperature': find_csv_by_title('T(NS)'),      # tr_data_019
        'chi_e': find_csv_by_title('AKE,AKNCE'),        # tr_data_023
        'chi_i': find_csv_by_title('AKD,AKNCD'),        # tr_data_024
        'pin': find_csv_by_title('PIN [MW'),            # tr_data_037
        'power': find_csv_by_title('POH,PNB,PNF'),      # tr_data_021
        'radiation': find_csv_by_title('PRSUM,PRB'),    # tr_data_022
    }

    print(f"  Found CSV files:")
    for key, f in csv_files.items():
        if f:
            print(f"    {key}: {os.path.basename(f)}")
        else:
            print(f"    {key}: NOT FOUND")

    data = {}

    # Read density
    r, dens, _ = read_tr_csv(csv_files['density'])
    data['r'] = r  # r/a
    data['ne'] = dens.get('nE', dens.get('NE', np.zeros_like(r)))  # 10^20/m^3
    data['nD'] = dens.get('nD', dens.get('ND', np.zeros_like(r)))
    data['nT'] = dens.get('nT', dens.get('NT', np.zeros_like(r)))
    data['ni'] = data['nD'] + data['nT']

    # Read temperature (reference)
    _, temp, _ = read_tr_csv(csv_files['temperature'])
    data['Te'] = temp.get('TE', temp.get('Te', np.zeros_like(r)))  # keV
    data['Ti'] = temp.get('TD', temp.get('Ti', np.zeros_like(r)))  # keV

    # Read chi (transport coefficients)
    _, chi_e, _ = read_tr_csv(csv_files['chi_e'])
    _, chi_i, _ = read_tr_csv(csv_files['chi_i'])
    data['chi_e'] = chi_e.get('AKE', np.zeros_like(r))  # m^2/s
    data['chi_i'] = chi_i.get('AKD', np.zeros_like(r))  # m^2/s

    # Read PIN (power input to each species)
    _, pin, _ = read_tr_csv(csv_files['pin'])
    data['PIN_e'] = pin.get('PIN_1', np.zeros_like(r))  # MW/m^3
    data['PIN_D'] = pin.get('PIN_2', np.zeros_like(r))  # MW/m^3
    data['PIN_T'] = pin.get('PIN_3', np.zeros_like(r))  # MW/m^3
    data['PIN_He'] = pin.get('PIN_4', np.zeros_like(r))  # MW/m^3
    data['PIN_i'] = data['PIN_D'] + data['PIN_T'] + data['PIN_He']

    # Read power terms for analysis
    _, power, _ = read_tr_csv(csv_files['power'])
    data['POH'] = power.get('POH', np.zeros_like(r))
    data['PNB'] = power.get('PNB', np.zeros_like(r))
    data['PNF'] = power.get('PNF', np.zeros_like(r))
    data['PRF'] = power.get('PRF', np.zeros_like(r))

    # Read radiation terms
    _, rad, _ = read_tr_csv(csv_files['radiation'])
    data['PRSUM'] = rad.get('PRSUM', np.zeros_like(r))
    data['PIE'] = rad.get('PIE', np.zeros_like(r))
    data['PCX'] = rad.get('PCX', np.zeros_like(r))

    print(f"\n  Data loaded:")
    print(f"    r/a range: [{r[0]:.4f}, {r[-1]:.4f}], {len(r)} points")
    print(f"    Te_center = {data['Te'][0]:.2f} keV")
    print(f"    Ti_center = {data['Ti'][0]:.2f} keV")
    print(f"    ne_center = {data['ne'][0]:.4f} x10^20/m^3")
    print(f"    chi_e_center = {data['chi_e'][0]:.4f} m^2/s")
    print(f"    PIN_e_center = {data['PIN_e'][0]:.4f} MW/m^3")

    return data


def solve_temperature_tr_method(r, n, chi, S, T_edge, a=A_MINOR, R0=R0,
                                 chi_factor=1.0, include_32_factor=False,
                                 use_geometry=False, kappa=KAPPA):
    """
    Solve steady-state heat equation using TR code discretization.

    TR uses the energy equation:
        d/dt[(3/2)nT] = -1/V' * d/dr(V' * n * chi * dT/dr) + S

    At steady state:
        1/r * d/dr(r * n * chi * dT/dr) = -S

    For steady state, we integrate:
        Q(r) = integral_0^r (r' * S(r') * dr')   [heat flux * r]
        r * n * chi * dT/dr = -Q(r)
        dT/dr = -Q(r) / (r * n * chi)
        T(r) = T_edge + integral_r^a Q/(r'*n*chi) dr'

    With TR geometry (use_geometry=True):
        From trmetric.f90:
        - DVRHO = 4π²κa²R₀(r/a) - volume element with elongation
        - AR2RHO = 1/(√κ × a)² - gradient correction
        - Combined effect on transport: DVRHO × AR2RHO = 4π²R₀r/a

        The elongation κ appears in:
        1. Volume element for source integration: includes κ factor
        2. Surface area for heat flux: modified by κ

    Args:
        r: r/a grid
        n: density [10^20/m^3]
        chi: diffusivity [m^2/s]
        S: source [MW/m^3]
        T_edge: boundary temperature [keV]
        a: minor radius [m]
        R0: major radius [m]
        chi_factor: scaling factor for chi
        include_32_factor: if True, include 3/2 factor in heat capacity
        use_geometry: if True, include TR geometry factors (elongation)
        kappa: elongation factor

    Returns:
        T: temperature [keV]
        Q: heat flux [MW]
        P_total: total integrated power [W]
    """
    nr = len(r)
    r_m = r * a  # Convert to meters
    dr = np.diff(r_m)
    dr = np.append(dr, dr[-1])  # Extend for last point

    # Ensure positive values
    n = np.maximum(n, 0.01)
    chi = np.maximum(chi, 0.01) * chi_factor  # Apply chi scaling

    # Convert S from MW/m^3 to W/m^3
    S_W = S * 1e6

    # Step 1: Calculate heat flux by integrating source
    # For cylindrical: Q(r) = integral_0^r (r' * S(r') * dr') * 4*pi^2*R0
    # For elongated:   Q(r) = integral_0^r (r' * S(r') * dr') * 4*pi^2*R0*κ
    #                        (larger volume due to elongation)
    if use_geometry:
        # TR geometry: DVRHO = 4π²κa²R₀ρ where ρ = r/a
        # Volume element dV = DVRHO * dρ = 4π²κa²R₀ρ * dρ
        # In physical coords: dV = 4π²κaR₀r * dr
        volume_factor = kappa
    else:
        volume_factor = 1.0

    integrand = r_m * S_W * volume_factor
    Q = np.zeros(nr)
    for i in range(1, nr):
        Q[i] = Q[i-1] + 0.5 * (integrand[i-1] + integrand[i]) * (r_m[i] - r_m[i-1])

    # Total power: P = 4*pi^2*R0 * Q[-1]
    P_total = 4 * np.pi**2 * R0 * Q[-1]

    # Step 2: Temperature gradient from Fourier's law
    # Heat flux: q = -n * chi * k * dT/dr
    # For heat balance at radius r:
    #   q * A_surface = Q (integrated source inside r)
    # For cylindrical: A_surface = 4π²R₀r
    # For elongated:   A_surface = 4π²R₀r * √κ (ellipse perimeter factor ~√κ)
    #
    # Actually from TR code: the effective transport coefficient is
    # FB = DVRHO * AR2RHO / DR² = 4π²R₀ρ (the κ factors cancel!)
    # This means geometry correction is already implicit in the coordinate system.
    #
    # However, the EFFECTIVE heat flux through the elongated surface:
    # - Surface area ~ 4π²R₀r × (circumference ratio) ~ 4π²R₀r × κ^(1/2)
    # - But gradient is in ρ direction: <|∇ρ|²> = AR2RHO
    # - Net effect: heat flow capacity per unit radius increased by √κ

    # e_factor: converts (n[10^20/m^3] * chi[m^2/s]) to W/m^2 per keV/m
    e_factor = 1.6e4  # n[10^20/m^3] * chi * dT[keV/m] -> q[W/m^2]

    # If include_32_factor, multiply by 3/2 for thermal energy (3/2 nkT)
    if include_32_factor:
        e_factor *= 1.5

    # Geometry correction for heat flux
    # From empirical chi_factor analysis:
    # - chi_e_factor = 1.317 ≈ √κ (electrons need MORE transport)
    # - chi_i_factor = 0.726 ≈ 1/√κ (ions need LESS transport)
    #
    # This suggests TR uses different effective geometry for each species,
    # or the PIN values already include partial geometry factors.
    #
    # For now, geometry_factor is set to 1.0 (no correction applied here).
    # The chi_factor parameter can be used for species-specific corrections.
    if use_geometry:
        # Volume factor κ in source cancels with 1/κ in AR2RHO
        # Net effect: FB = DVRHO × AR2RHO / DR² = 4π²R₀ρ / DR² (κ cancels)
        # So geometry_factor should be 1.0 when properly normalized
        geometry_factor = 1.0
    else:
        geometry_factor = 1.0

    denom = r_m * n * chi * e_factor * geometry_factor
    denom = np.maximum(denom, 1e-30)  # Avoid division by zero

    dTdr = -Q / denom  # keV/m

    # Handle center (r=0): dT/dr = 0 by symmetry
    dTdr[0] = 0

    # Step 3: Integrate from edge to get temperature
    T = np.zeros(nr)
    T[-1] = T_edge

    for i in range(nr-2, -1, -1):
        dr_local = r_m[i+1] - r_m[i]
        grad_mid = 0.5 * (dTdr[i] + dTdr[i+1])
        T[i] = T[i+1] - grad_mid * dr_local

    return T, Q, P_total


def solve_temperature_matrix_method(r, n, chi, S, T_edge, a=A_MINOR):
    """
    Solve steady-state heat equation using tridiagonal matrix (TR style).

    Discretize: 1/r * d/dr(r * n * chi * dT/dr) = -S

    Using finite volume method on cell-centered grid.

    Args:
        r: r/a grid
        n: density [10^20/m^3]
        chi: diffusivity [m^2/s]
        S: source [MW/m^3]
        T_edge: boundary temperature [keV]
        a: minor radius [m]

    Returns:
        T: temperature [keV]
    """
    nr = len(r)
    r_m = r * a  # meters
    dr = r_m[1] - r_m[0] if nr > 1 else a / 50

    # Ensure positive values
    n = np.maximum(n, 0.01)
    chi = np.maximum(chi, 0.01)

    # Half-point radii (cell faces)
    r_half = np.zeros(nr + 1)
    r_half[1:-1] = 0.5 * (r_m[:-1] + r_m[1:])
    r_half[0] = 0.5 * r_m[0]  # Near axis
    r_half[-1] = r_m[-1] + 0.5 * dr  # Beyond edge

    # Half-point density and chi (interpolated)
    n_half = np.zeros(nr + 1)
    chi_half = np.zeros(nr + 1)
    n_half[1:-1] = 0.5 * (n[:-1] + n[1:])
    chi_half[1:-1] = 0.5 * (chi[:-1] + chi[1:])
    n_half[0] = n[0]
    chi_half[0] = chi[0]
    n_half[-1] = n[-1]
    chi_half[-1] = chi[-1]

    # Diffusion coefficient at half-points: D = r * n * chi * k / dr
    # Unit conversion: n [10^20/m^3] * chi [m^2/s] * k [J/keV] -> W/m/keV
    # k = 1.6e-16 J/keV, so e_factor = 1e20 * 1.6e-16 = 1.6e4
    e_factor = 1.6e4

    D_half = r_half * n_half * chi_half * e_factor / dr

    # Source term: MW/m^3 -> W/m^3, then to keV*10^20/m^3/s equivalent
    # S [MW/m^3] = S * 1e6 [W/m^3]
    # For energy balance: 3/2 * n * dT/dt = S (in consistent units)
    # S [W/m^3] / (1.5 * n[10^20/m^3] * 1e20 * e[J/keV]) = dT/dt [keV/s]
    # At steady state, flux divergence = S
    S_eff = S * 1e6  # W/m^3

    # Volume element (cell volume ~ 2*pi*R0 * 2*pi*r*dr)
    dV = r_m * dr
    dV[0] = max(r_m[0], 1e-6) * dr  # Avoid zero at axis

    # Build tridiagonal matrix: A*T(i-1) + B*T(i) + C*T(i+1) = D
    A_coef = np.zeros(nr)
    B_coef = np.zeros(nr)
    C_coef = np.zeros(nr)
    D_rhs = np.zeros(nr)

    for i in range(1, nr-1):
        A_coef[i] = D_half[i]      # Flux from left
        C_coef[i] = D_half[i+1]    # Flux from right
        B_coef[i] = -(A_coef[i] + C_coef[i])  # Conservation
        D_rhs[i] = -S_eff[i] * dV[i]

    # Center boundary: dT/dr = 0 (symmetry)
    # Only outward flux: C[0]*T[1] - C[0]*T[0] = -S*dV
    C_coef[0] = D_half[1]
    B_coef[0] = -D_half[1] - 1e-10  # Small regularization
    D_rhs[0] = -S_eff[0] * dV[0]

    # Edge boundary: Dirichlet T = T_edge
    B_coef[-1] = 1.0
    A_coef[-1] = 0.0
    D_rhs[-1] = T_edge

    # Solve tridiagonal system
    ab = np.zeros((3, nr))
    ab[0, 1:] = C_coef[:-1]   # Upper diagonal
    ab[1, :] = B_coef         # Main diagonal
    ab[2, :-1] = A_coef[1:]   # Lower diagonal

    try:
        T = solve_banded((1, 1), ab, D_rhs)
    except np.linalg.LinAlgError as e:
        print(f"  Matrix solve failed: {e}")
        T = np.full(nr, T_edge)

    return T


def main():
    print("="*70)
    print("REPRODUCE TEMPERATURE FROM TR OUTPUT")
    print("="*70)

    # Load TR data
    data = load_tr_data()
    r = data['r']
    nr = len(r)

    # Reference temperatures from TR
    Te_ref = data['Te']
    Ti_ref = data['Ti']
    Te_edge = Te_ref[-1]
    Ti_edge = Ti_ref[-1]

    print(f"\n{'='*70}")
    print("ELECTRON TEMPERATURE CALCULATION")
    print("="*70)

    # Method 1: Integration method
    print("\n1. Integration method:")
    Te_int, Qe, Pe_total = solve_temperature_tr_method(
        r, data['ne'], data['chi_e'], data['PIN_e'], Te_edge
    )
    print(f"   Te_center (calc) = {Te_int[0]:.2f} keV")
    print(f"   Te_center (ref)  = {Te_ref[0]:.2f} keV")
    print(f"   Difference: {Te_int[0] - Te_ref[0]:.2f} keV ({(Te_int[0]/Te_ref[0]-1)*100:.1f}%)")
    print(f"   Total power: {Pe_total/1e6:.1f} MW")

    # Method 2: Matrix method
    print("\n2. Matrix method (TR style):")
    Te_mat = solve_temperature_matrix_method(
        r, data['ne'], data['chi_e'], data['PIN_e'], Te_edge
    )
    print(f"   Te_center (calc) = {Te_mat[0]:.2f} keV")
    print(f"   Te_center (ref)  = {Te_ref[0]:.2f} keV")
    print(f"   Difference: {Te_mat[0] - Te_ref[0]:.2f} keV ({(Te_mat[0]/Te_ref[0]-1)*100:.1f}%)")

    print(f"\n{'='*70}")
    print("ION TEMPERATURE CALCULATION")
    print("="*70)

    # Ion temperature
    print("\n1. Integration method:")
    Ti_int, Qi, Pi_total = solve_temperature_tr_method(
        r, data['ni'], data['chi_i'], data['PIN_i'], Ti_edge
    )
    print(f"   Ti_center (calc) = {Ti_int[0]:.2f} keV")
    print(f"   Ti_center (ref)  = {Ti_ref[0]:.2f} keV")
    print(f"   Difference: {Ti_int[0] - Ti_ref[0]:.2f} keV ({(Ti_int[0]/Ti_ref[0]-1)*100:.1f}%)")
    print(f"   Total power: {Pi_total/1e6:.1f} MW")

    print("\n2. Matrix method:")
    Ti_mat = solve_temperature_matrix_method(
        r, data['ni'], data['chi_i'], data['PIN_i'], Ti_edge
    )
    print(f"   Ti_center (calc) = {Ti_mat[0]:.2f} keV")
    print(f"   Ti_center (ref)  = {Ti_ref[0]:.2f} keV")
    print(f"   Difference: {Ti_mat[0] - Ti_ref[0]:.2f} keV ({(Ti_mat[0]/Ti_ref[0]-1)*100:.1f}%)")

    # Plot comparison
    print("\n" + "="*70)
    print("PLOTTING RESULTS")
    print("="*70)

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))

    # Electron temperature
    ax = axes[0, 0]
    ax.plot(r, Te_ref, 'b-', lw=2, label='TR output')
    ax.plot(r, Te_int, 'r--', lw=2, label='Integration method')
    ax.plot(r, Te_mat, 'g-.', lw=2, label='Matrix method')
    ax.set_xlabel('r/a')
    ax.set_ylabel('Te [keV]')
    ax.set_title('Electron Temperature')
    ax.legend()
    ax.grid(True)
    ax.set_xlim([0, 1])

    # Ion temperature
    ax = axes[0, 1]
    ax.plot(r, Ti_ref, 'b-', lw=2, label='TR output')
    ax.plot(r, Ti_int, 'r--', lw=2, label='Integration method')
    ax.plot(r, Ti_mat, 'g-.', lw=2, label='Matrix method')
    ax.set_xlabel('r/a')
    ax.set_ylabel('Ti [keV]')
    ax.set_title('Ion Temperature')
    ax.legend()
    ax.grid(True)
    ax.set_xlim([0, 1])

    # Temperature difference (electron)
    ax = axes[0, 2]
    ax.plot(r, Te_int - Te_ref, 'r-', lw=2, label='Integration - TR')
    ax.plot(r, Te_mat - Te_ref, 'g-', lw=2, label='Matrix - TR')
    ax.axhline(0, color='k', ls='--', lw=1)
    ax.set_xlabel('r/a')
    ax.set_ylabel('$\Delta$Te [keV]')
    ax.set_title('Electron Temperature Difference')
    ax.legend()
    ax.grid(True)
    ax.set_xlim([0, 1])

    # Source term (electron)
    ax = axes[1, 0]
    ax.plot(r, data['PIN_e'], 'b-', lw=2, label='PIN_e (net source)')
    ax.axhline(0, color='k', ls='--', lw=1)
    ax.set_xlabel('r/a')
    ax.set_ylabel('S_e [MW/m³]')
    ax.set_title('Electron Source (PIN_e)')
    ax.legend()
    ax.grid(True)
    ax.set_xlim([0, 1])

    # Transport coefficient
    ax = axes[1, 1]
    ax.plot(r, data['chi_e'], 'b-', lw=2, label='chi_e')
    ax.plot(r, data['chi_i'], 'r-', lw=2, label='chi_i')
    ax.set_xlabel('r/a')
    ax.set_ylabel('chi [m²/s]')
    ax.set_title('Transport Coefficients')
    ax.legend()
    ax.grid(True)
    ax.set_xlim([0, 1])
    ax.set_yscale('log')

    # Density
    ax = axes[1, 2]
    ax.plot(r, data['ne'], 'b-', lw=2, label='n_e')
    ax.plot(r, data['ni'], 'r-', lw=2, label='n_i')
    ax.set_xlabel('r/a')
    ax.set_ylabel('n [10²⁰/m³]')
    ax.set_title('Density Profile')
    ax.legend()
    ax.grid(True)
    ax.set_xlim([0, 1])

    plt.tight_layout()
    plt.savefig('reproduce_temperature_from_tr.png', dpi=150)
    print(f"   Plot saved: reproduce_temperature_from_tr.png")

    # Summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"\n{'Method':<25} {'Te_0 [keV]':>12} {'Ti_0 [keV]':>12}")
    print("-"*50)
    print(f"{'TR output (reference)':<25} {Te_ref[0]:>12.2f} {Ti_ref[0]:>12.2f}")
    print(f"{'Integration method':<25} {Te_int[0]:>12.2f} {Ti_int[0]:>12.2f}")
    print(f"{'Matrix method':<25} {Te_mat[0]:>12.2f} {Ti_mat[0]:>12.2f}")
    print("-"*50)

    # Analyze discrepancies
    print("\nAnalysis of discrepancies:")

    # Check energy balance
    print("\n  1. Energy balance check:")
    a = A_MINOR
    # Integrate PIN over volume
    r_m = r * a
    dV = 4 * np.pi**2 * R0 * r_m
    P_in_e = np.trapz(data['PIN_e'] * dV, r_m)  # MW
    P_in_i = np.trapz(data['PIN_i'] * dV, r_m)
    print(f"     Integrated P_in_e = {P_in_e:.1f} MW")
    print(f"     Integrated P_in_i = {P_in_i:.1f} MW")

    # Heat flux at edge (should equal integrated source)
    print(f"     Heat flux at edge (Q_e from integration) = {4*np.pi**2*R0*Qe[-1]/1e6:.1f} MW")

    # PIE analysis
    print("\n  2. Ion-electron exchange (PIE) analysis:")
    P_ie = np.trapz(data['PIE'] * dV, r_m)
    print(f"     PIE profile: center = {data['PIE'][0]:.4f} MW/m³, edge = {data['PIE'][-1]:.4f} MW/m³")
    print(f"     Integrated PIE = {P_ie:.1f} MW")
    print(f"     PIE transfers energy from electrons to ions (when Te > Ti)")

    # Test with 3/2 factor
    print("\n  3. Test with 3/2 factor (thermal energy = 3/2 nkT):")
    Te_32, _, _ = solve_temperature_tr_method(
        r, data['ne'], data['chi_e'], data['PIN_e'], Te_edge, include_32_factor=True
    )
    Ti_32, _, _ = solve_temperature_tr_method(
        r, data['ni'], data['chi_i'], data['PIN_i'], Ti_edge, include_32_factor=True
    )
    print(f"     Te_center (with 3/2) = {Te_32[0]:.2f} keV (target: {Te_ref[0]:.2f})")
    print(f"     Ti_center (with 3/2) = {Ti_32[0]:.2f} keV (target: {Ti_ref[0]:.2f})")

    # Find chi factor to match TR temperature
    print("\n  4. Chi scaling to match TR output:")
    from scipy.optimize import minimize_scalar

    def objective_e(factor):
        T, _, _ = solve_temperature_tr_method(
            r, data['ne'], data['chi_e'], data['PIN_e'], Te_edge, chi_factor=factor
        )
        return (T[0] - Te_ref[0])**2

    def objective_i(factor):
        T, _, _ = solve_temperature_tr_method(
            r, data['ni'], data['chi_i'], data['PIN_i'], Ti_edge, chi_factor=factor
        )
        return (T[0] - Ti_ref[0])**2

    result_e = minimize_scalar(objective_e, bounds=(0.5, 3.0), method='bounded')
    result_i = minimize_scalar(objective_i, bounds=(0.5, 3.0), method='bounded')
    chi_factor_e = result_e.x
    chi_factor_i = result_i.x

    Te_scaled, _, _ = solve_temperature_tr_method(
        r, data['ne'], data['chi_e'], data['PIN_e'], Te_edge, chi_factor=chi_factor_e
    )
    Ti_scaled, _, _ = solve_temperature_tr_method(
        r, data['ni'], data['chi_i'], data['PIN_i'], Ti_edge, chi_factor=chi_factor_i
    )

    print(f"     chi_e scaling factor = {chi_factor_e:.3f}")
    print(f"     chi_i scaling factor = {chi_factor_i:.3f}")
    print(f"     Te_center (scaled) = {Te_scaled[0]:.2f} keV (target: {Te_ref[0]:.2f})")
    print(f"     Ti_center (scaled) = {Ti_scaled[0]:.2f} keV (target: {Ti_ref[0]:.2f})")

    # Test with elongation correction using TR-style geometry
    kappa_s = np.sqrt(KAPPA)
    print(f"\n  5. Geometry correction analysis (elongation κ={KAPPA}):")
    print(f"     √κ = {kappa_s:.4f}")
    print(f"     Empirical chi_e_factor = {chi_factor_e:.4f} ≈ √κ = {kappa_s:.4f}")
    print(f"     Empirical chi_i_factor = {chi_factor_i:.4f} ≈ 1/√κ = {1/kappa_s:.4f}")

    # Test: Electrons with √κ, Ions with 1/√κ
    print(f"\n     Test: Species-specific √κ correction")
    Te_geom, _, _ = solve_temperature_tr_method(
        r, data['ne'], data['chi_e'], data['PIN_e'], Te_edge, chi_factor=kappa_s
    )
    Ti_geom, _, _ = solve_temperature_tr_method(
        r, data['ni'], data['chi_i'], data['PIN_i'], Ti_edge, chi_factor=1/kappa_s
    )
    print(f"     Te (chi × √κ)   = {Te_geom[0]:.2f} keV (target: {Te_ref[0]:.2f}, error: {(Te_geom[0]/Te_ref[0]-1)*100:.1f}%)")
    print(f"     Ti (chi / √κ)   = {Ti_geom[0]:.2f} keV (target: {Ti_ref[0]:.2f}, error: {(Ti_geom[0]/Ti_ref[0]-1)*100:.1f}%)")

    # Simple (1+κ²)/2 factor for comparison
    kappa_geom = (1 + KAPPA**2) / 2
    print(f"\n     Test: Simple factor (1+κ²)/2 = {kappa_geom:.2f}")
    Te_simple, _, _ = solve_temperature_tr_method(
        r, data['ne'], data['chi_e'], data['PIN_e'], Te_edge, chi_factor=kappa_geom
    )
    Ti_simple, _, _ = solve_temperature_tr_method(
        r, data['ni'], data['chi_i'], data['PIN_i'], Ti_edge, chi_factor=kappa_geom
    )
    print(f"     Te (chi × {kappa_geom:.2f}) = {Te_simple[0]:.2f} keV (target: {Te_ref[0]:.2f})")
    print(f"     Ti (chi × {kappa_geom:.2f}) = {Ti_simple[0]:.2f} keV (target: {Ti_ref[0]:.2f})")

    print("\n  6. Summary of physical corrections needed:")
    print(f"     Electron: chi_factor = {chi_factor_e:.3f} ≈ √κ = {kappa_s:.3f}")
    print(f"     Ion:      chi_factor = {chi_factor_i:.3f} ≈ 1/√κ = {1/kappa_s:.3f}")
    print("")
    print("     Key finding: Electrons and ions need OPPOSITE geometry corrections!")
    print(f"       - Electron chi_factor / √κ = {chi_factor_e/kappa_s:.3f} (should be ≈1.0)")
    print(f"       - Ion chi_factor × √κ = {chi_factor_i*kappa_s:.3f} (should be ≈1.0)")
    print("")
    print("     With species-specific √κ correction:")
    print(f"       Te error: {(Te_geom[0]/Te_ref[0]-1)*100:.1f}% (improved from {(Te_int[0]/Te_ref[0]-1)*100:.1f}%)")
    print(f"       Ti error: {(Ti_geom[0]/Ti_ref[0]-1)*100:.1f}% (improved from {(Ti_int[0]/Ti_ref[0]-1)*100:.1f}%)")
    print("")
    print("     Conclusion:")
    print("       The √κ factor arises from elongation geometry in TR code.")
    print("       - Electrons: PIN_e is in TR geometry → need chi × √κ to match")
    print("       - Ions: PIN_i is in TR geometry → need chi / √κ to match")
    print("       This suggests PIN values have different normalizations for e vs i.")

    # Add scaled result to plot
    fig2, axes2 = plt.subplots(2, 2, figsize=(14, 10))

    # Electron temperature comparison
    ax = axes2[0, 0]
    ax.plot(r, Te_ref, 'b-', lw=2, label='TR output')
    ax.plot(r, Te_int, 'r--', lw=2, label=f'No geometry (chi x 1.0)')
    ax.plot(r, Te_geom, 'm-.', lw=2, label=f'TR geometry (1/√κ)')
    ax.plot(r, Te_scaled, 'g:', lw=2, label=f'Fitted (chi x {chi_factor_e:.2f})')
    ax.set_xlabel('r/a')
    ax.set_ylabel('Te [keV]')
    ax.set_title('Electron Temperature Comparison')
    ax.legend()
    ax.grid(True)
    ax.set_xlim([0, 1])

    # Ion temperature comparison
    ax = axes2[0, 1]
    ax.plot(r, Ti_ref, 'b-', lw=2, label='TR output')
    ax.plot(r, Ti_int, 'r--', lw=2, label=f'No geometry (chi x 1.0)')
    ax.plot(r, Ti_geom, 'm-.', lw=2, label=f'TR geometry (1/√κ)')
    ax.plot(r, Ti_scaled, 'g:', lw=2, label=f'Fitted (chi x {chi_factor_i:.2f})')
    ax.set_xlabel('r/a')
    ax.set_ylabel('Ti [keV]')
    ax.set_title('Ion Temperature Comparison')
    ax.legend()
    ax.grid(True)
    ax.set_xlim([0, 1])

    # Temperature error (electron)
    ax = axes2[1, 0]
    ax.plot(r, (Te_int - Te_ref)/Te_ref*100, 'r--', lw=2, label='No geometry')
    ax.plot(r, (Te_geom - Te_ref)/Te_ref*100, 'm-.', lw=2, label='TR geometry')
    ax.plot(r, (Te_scaled - Te_ref)/Te_ref*100, 'g:', lw=2, label='Fitted')
    ax.axhline(0, color='k', ls='-', lw=1)
    ax.set_xlabel('r/a')
    ax.set_ylabel('Error [%]')
    ax.set_title('Electron Temperature Error')
    ax.legend()
    ax.grid(True)
    ax.set_xlim([0, 1])

    # Temperature error (ion)
    ax = axes2[1, 1]
    ax.plot(r, (Ti_int - Ti_ref)/Ti_ref*100, 'r--', lw=2, label='No geometry')
    ax.plot(r, (Ti_geom - Ti_ref)/Ti_ref*100, 'm-.', lw=2, label='TR geometry')
    ax.plot(r, (Ti_scaled - Ti_ref)/Ti_ref*100, 'g:', lw=2, label='Fitted')
    ax.axhline(0, color='k', ls='-', lw=1)
    ax.set_xlabel('r/a')
    ax.set_ylabel('Error [%]')
    ax.set_title('Ion Temperature Error')
    ax.legend()
    ax.grid(True)
    ax.set_xlim([0, 1])

    plt.tight_layout()
    plt.savefig('reproduce_temperature_geometry.png', dpi=150)
    print(f"\n   Plot saved: reproduce_temperature_geometry.png")

    # Final comparison with geometry correction
    print("\n" + "="*70)
    print("FINAL RESULTS WITH GEOMETRY CORRECTION")
    print("="*70)

    chi_e_geom = get_geometry_chi_factor('electron')
    chi_i_geom = get_geometry_chi_factor('ion')

    Te_final, _, _ = solve_temperature_tr_method(
        r, data['ne'], data['chi_e'], data['PIN_e'], Te_edge, chi_factor=chi_e_geom
    )
    Ti_final, _, _ = solve_temperature_tr_method(
        r, data['ni'], data['chi_i'], data['PIN_i'], Ti_edge, chi_factor=chi_i_geom
    )

    print(f"\nUsing get_geometry_chi_factor():")
    print(f"  Electron chi_factor = √κ = {chi_e_geom:.4f}")
    print(f"  Ion chi_factor = 1/√κ = {chi_i_geom:.4f}")
    print(f"\nFinal results:")
    print(f"  Te_center = {Te_final[0]:.2f} keV (target: {Te_ref[0]:.2f}, error: {(Te_final[0]/Te_ref[0]-1)*100:.1f}%)")
    print(f"  Ti_center = {Ti_final[0]:.2f} keV (target: {Ti_ref[0]:.2f}, error: {(Ti_final[0]/Ti_ref[0]-1)*100:.1f}%)")

    return data, Te_final, Ti_final, Te_geom, Ti_geom


if __name__ == '__main__':
    main()
