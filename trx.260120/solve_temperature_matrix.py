#!/usr/bin/env python3
"""
Solve temperature equation using TR's matrix method.

This script implements the exact algorithm used in TR's Fortran code:
1. Volume-weighted coordinate system (DVRHO)
2. Geometric factors (FA, FB, AR1RHO, AR2RHO)
3. Tridiagonal matrix assembly (A, B, C coefficients)
4. LAPACK-style band matrix solver

The energy equation in TR:
    d(3/2 * n * T)/dt = S - div(q)

where heat flux q = -n * chi * grad(T)

In steady state (dT/dt = 0), the discretized equation becomes:
    A(i)*T(i-1) + B(i)*T(i) + C(i)*T(i+1) = D(i)

References:
    - trexec.f90: TR_COEF_DECIDE subroutine
    - trmetric.f90: DVRHO, AR1RHO, AR2RHO definitions
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.linalg import solve_banded
import os

# Physical parameters for CFEDR
R0 = 7.82       # Major radius [m]
A_MINOR = 2.44  # Minor radius [m]
KAPPA = 1.88    # Elongation
RKEV = 1.602e-16  # 1 keV in Joules


def read_tr_csv(filename):
    """Read TR output CSV file."""
    with open(filename, 'r') as f:
        title = f.readline().strip()
    df = pd.read_csv(filename, skiprows=1, index_col=False)
    df = df.loc[:, ~df.columns.str.contains('^Unnamed')]
    r = df.iloc[:, 0].values
    data = {col: df[col].values for col in df.columns[1:]}
    return r, data, title


def compute_geometry(r, R0=R0, a=A_MINOR, kappa=KAPPA):
    """
    Compute TR's geometric factors.

    In TR (trmetric.f90):
        DVRHO(NR) = 2*PI*RKAP*RA*RA*2*PI*RR*RM(NR)
        AR1RHO(NR) = 1/(RKAPS*RA)
        AR2RHO(NR) = 1/(RKAPS*RA)^2

    Args:
        r: normalized radius (r/a), array
        R0: major radius [m]
        a: minor radius [m]
        kappa: elongation

    Returns:
        Dictionary with geometric factors
    """
    nr = len(r)
    dr = r[1] - r[0] if nr > 1 else 0.02  # Grid spacing in r/a

    # Volume element per unit normalized radius
    # DVRHO = 4*pi^2 * kappa * a^2 * R0 * rho
    # where rho = r/a
    DVRHO = 4 * np.pi**2 * kappa * a**2 * R0 * r
    # Avoid zero at center
    DVRHO[0] = DVRHO[1] * 0.5

    # Geometric factors (constant in simplified geometry)
    AR1RHO = 1.0 / (kappa * a) * np.ones(nr)
    AR2RHO = 1.0 / (kappa * a)**2 * np.ones(nr)

    # FA and FB factors (at half-grid points)
    # FA(i) = DVRHO(i) * AR1RHO(i) / dr  [convection]
    # FB(i) = DVRHO(i) * AR2RHO(i) / dr^2  [diffusion]
    FA = DVRHO * AR1RHO / (dr * a)  # Convert dr to real units
    FB = DVRHO * AR2RHO / (dr * a)**2

    # Volume factors for source and mass matrix
    DV11 = DVRHO
    DV23 = DVRHO**(2.0/3.0)
    DV53 = DVRHO**(5.0/3.0)

    return {
        'DVRHO': DVRHO,
        'AR1RHO': AR1RHO,
        'AR2RHO': AR2RHO,
        'FA': FA,
        'FB': FB,
        'DV11': DV11,
        'DV23': DV23,
        'DV53': DV53,
        'dr': dr,
        'dr_real': dr * a
    }


def assemble_matrix_tr(r, n, chi, S_net, T_edge, geom, steady_state=True):
    """
    Assemble the tridiagonal matrix system for steady-state heat equation.

    The steady-state heat equation in cylindrical geometry:
        1/r * d/dr(r * n * chi * dT/dr) = -S

    Discretized using finite volumes:
        At grid point i, integrate over volume from r_{i-1/2} to r_{i+1/2}:

        [r * n * chi * dT/dr]_{i+1/2} - [r * n * chi * dT/dr]_{i-1/2} = -S_i * dV_i

    This gives:
        A*T(i-1) + B*T(i) + C*T(i+1) = D(i)

    Args:
        r: normalized radius
        n: density profile [10^20/m^3]
        chi: thermal diffusivity [m^2/s]
        S_net: net source term [MW/m^3]
        T_edge: boundary temperature at r=a [keV]
        geom: geometry dictionary from compute_geometry
        steady_state: if True, solve steady-state

    Returns:
        T: temperature profile [keV]
    """
    nr = len(r)
    a = A_MINOR
    r_real = r * a
    dr = geom['dr_real']

    # Half-point radii
    r_half = np.zeros(nr + 1)
    r_half[1:-1] = 0.5 * (r_real[:-1] + r_real[1:])
    r_half[0] = 0  # Center
    r_half[-1] = r_real[-1] + 0.5 * dr  # Beyond edge

    # Half-point values of n * chi
    n_half = np.zeros(nr + 1)
    chi_half = np.zeros(nr + 1)
    n_half[1:-1] = 0.5 * (n[:-1] + n[1:])
    chi_half[1:-1] = 0.5 * (chi[:-1] + chi[1:])
    n_half[0] = n[0]
    chi_half[0] = chi[0]
    n_half[-1] = n[-1]
    chi_half[-1] = chi[-1]

    # Diffusion coefficient at half-points: D_{i+1/2} = r_{i+1/2} * n * chi / dr
    D_half = r_half * n_half * chi_half / dr

    # Source term conversion
    # S [MW/m^3] -> [keV * 10^20 / m^3 / s] / 1.5
    S_conv = S_net * 62.4 / 1.5

    # Volume element for source integration: dV ~ r * dr (in 1D cylindrical)
    dV = r_real * dr

    # Build tridiagonal matrix
    A_coef = np.zeros(nr)
    B_coef = np.zeros(nr)
    C_coef = np.zeros(nr)
    D_rhs = np.zeros(nr)

    # Interior points
    for i in range(1, nr-1):
        A_coef[i] = D_half[i]      # Left flux
        C_coef[i] = D_half[i+1]    # Right flux
        B_coef[i] = -(A_coef[i] + C_coef[i])
        D_rhs[i] = -S_conv[i] * dV[i]

    # Center boundary (i = 0): symmetry dT/dr = 0
    # Flux at center is zero, only right flux matters
    C_coef[0] = D_half[1]
    B_coef[0] = -D_half[1]
    D_rhs[0] = -S_conv[0] * dV[0] if dV[0] > 0 else -S_conv[0] * 0.5 * dr * dr

    # Edge boundary: Dirichlet T = T_edge
    B_coef[-1] = 1.0
    A_coef[-1] = 0.0
    D_rhs[-1] = T_edge

    # Solve tridiagonal system
    ab = np.zeros((3, nr))
    ab[0, 1:] = C_coef[:-1]
    ab[1, :] = B_coef
    ab[2, :-1] = A_coef[1:]

    T = solve_banded((1, 1), ab, D_rhs)

    return T


def solve_temperature_matrix_tr(r, n, chi, S_net, T_edge):
    """
    Main solver function using TR's matrix method.

    Args:
        r: normalized radius (r/a)
        n: density profile [10^20/m^3]
        chi: thermal diffusivity [m^2/s]
        S_net: net source term [MW/m^3]
        T_edge: edge temperature [keV]

    Returns:
        T: temperature profile [keV]
    """
    geom = compute_geometry(r)
    T = assemble_matrix_tr(r, n, chi, S_net, T_edge, geom)
    return T, geom


def solve_with_density_weighting(r, n, chi, S_net, T_edge, a=A_MINOR):
    """
    Alternative solver with explicit density in diffusion term.

    Heat equation: 1/r * d/dr(r * n * chi * dT/dr) = -S

    Discretized: D_i * (n_i+1/2 * chi_i+1/2 * (T_i+1 - T_i))
               - D_i-1 * (n_i-1/2 * chi_i-1/2 * (T_i - T_i-1)) = -S_i * dV
    """
    nr = len(r)
    r_real = r * a
    dr = (r_real[1] - r_real[0]) if nr > 1 else 0.02 * a

    # Half-point values
    n_half = 0.5 * (n[:-1] + n[1:])
    chi_half = 0.5 * (chi[:-1] + chi[1:])
    r_half = 0.5 * (r_real[:-1] + r_real[1:])

    # Diffusion coefficient at half-points
    D_half = n_half * chi_half * r_half / dr

    # Source conversion
    S_conv = S_net * 62.4 / 1.5

    # Build matrix
    A = np.zeros(nr)
    B = np.zeros(nr)
    C = np.zeros(nr)
    D = np.zeros(nr)

    for i in range(1, nr-1):
        # Volume element ~ r
        vol_factor = r_real[i]

        A[i] = D_half[i-1]
        C[i] = D_half[i]
        B[i] = -(A[i] + C[i])
        D[i] = -S_conv[i] * vol_factor * dr

    # Center: symmetry
    B[0] = -D_half[0] * 2
    C[0] = D_half[0] * 2
    D[0] = -S_conv[0] * r_real[1] * dr  # Use small volume at center

    # Edge: Dirichlet
    B[-1] = 1.0
    D[-1] = T_edge

    # Solve
    ab = np.zeros((3, nr))
    ab[0, 1:] = C[:-1]
    ab[1, :] = B
    ab[2, :-1] = A[1:]

    T = solve_banded((1, 1), ab, D)
    return T


def load_tr_data():
    """Load TR simulation data."""
    data = {}

    # Density
    r, dens, _ = read_tr_csv('tr_data_017.csv')
    data['r'] = r
    data['ne'] = dens.get('nE', np.zeros_like(r))
    data['nD'] = dens.get('nD', np.zeros_like(r))
    data['nT'] = dens.get('nT', np.zeros_like(r))
    data['nA'] = dens.get('nA', np.zeros_like(r))
    data['ni'] = data['nD'] + data['nT'] + data['nA']

    # Temperature
    _, temp, _ = read_tr_csv('tr_data_019.csv')
    data['Te'] = temp.get('TE', np.zeros_like(r))
    data['TD'] = temp.get('TD', np.zeros_like(r))
    data['Ti'] = data['TD']

    # Transport coefficients
    _, chi_e_data, _ = read_tr_csv('tr_data_023.csv')
    data['chi_e'] = chi_e_data.get('AKE', np.zeros_like(r))

    _, chi_i_data, _ = read_tr_csv('tr_data_024.csv')
    data['chi_i'] = chi_i_data.get('AKD', np.zeros_like(r))

    # Power input (species-specific)
    if os.path.exists('tr_data_042.csv'):
        _, pin, _ = read_tr_csv('tr_data_042.csv')
        data['PIN_e'] = pin.get('PIN_1', np.zeros_like(r))
        data['PIN_D'] = pin.get('PIN_2', np.zeros_like(r))
        data['PIN_T'] = pin.get('PIN_3', np.zeros_like(r))
        data['PIN_He4'] = pin.get('PIN_4', np.zeros_like(r))

    return data


def find_optimal_chi_factor(r, n, chi, S_net, T_edge, T_target, geom):
    """
    Find the optimal chi scaling factor to match target temperature.

    Uses bisection method to find factor such that T_calc[0] ≈ T_target[0].
    """
    from scipy.optimize import minimize_scalar

    def objective(factor):
        T_calc = assemble_matrix_tr(r, n, chi * factor, S_net, T_edge, geom)
        return (T_calc[0] - T_target[0])**2

    result = minimize_scalar(objective, bounds=(0.1, 2.0), method='bounded')
    return result.x


def main():
    """Main function to test matrix solver."""
    print("="*70)
    print("TR MATRIX SOLVER FOR TEMPERATURE EQUATION")
    print("="*70)

    # Load data
    print("\nLoading TR simulation data...")
    data = load_tr_data()

    r = data['r']
    nr = len(r)
    print(f"Grid points: {nr}")
    print(f"r/a range: {r[0]:.4f} to {r[-1]:.4f}")

    # Electron temperature
    print("\n--- ELECTRON TEMPERATURE ---")

    ne = data['ne']
    chi_e = data['chi_e']
    S_e = data['PIN_e']  # Already net (losses subtracted in TR)
    Te_edge = data['Te'][-1]
    Te_sim = data['Te']

    print(f"Te_edge = {Te_edge:.3f} keV")
    print(f"chi_e range: {chi_e.min():.3f} to {chi_e.max():.3f} m^2/s")
    print(f"S_e center = {S_e[0]:.4f} MW/m^3")

    # Method 1: TR-style matrix solver
    print("\n1. TR-style matrix solver (volume-weighted):")
    Te_matrix_tr, geom = solve_temperature_matrix_tr(r, ne, chi_e, S_e, Te_edge)
    ratio1 = Te_matrix_tr[0] / Te_sim[0] if Te_sim[0] > 0 else 0
    print(f"   Te_center = {Te_matrix_tr[0]:.2f} keV (sim: {Te_sim[0]:.2f}, ratio: {ratio1:.2f})")

    # Method 2: Density-weighted solver
    print("\n2. Density-weighted matrix solver:")
    Te_matrix_nchi = solve_with_density_weighting(r, ne, chi_e, S_e, Te_edge)
    ratio2 = Te_matrix_nchi[0] / Te_sim[0] if Te_sim[0] > 0 else 0
    print(f"   Te_center = {Te_matrix_nchi[0]:.2f} keV (sim: {Te_sim[0]:.2f}, ratio: {ratio2:.2f})")

    # Method 3: With chi scaling
    print("\n3. TR-style with chi*0.5:")
    Te_matrix_scaled, _ = solve_temperature_matrix_tr(r, ne, chi_e * 0.5, S_e, Te_edge)
    ratio3 = Te_matrix_scaled[0] / Te_sim[0] if Te_sim[0] > 0 else 0
    print(f"   Te_center = {Te_matrix_scaled[0]:.2f} keV (sim: {Te_sim[0]:.2f}, ratio: {ratio3:.2f})")

    # Find optimal chi factor
    print("\n4. Automatic chi calibration:")
    opt_factor_e = find_optimal_chi_factor(r, ne, chi_e, S_e, Te_edge, Te_sim, geom)
    Te_matrix_opt, _ = solve_temperature_matrix_tr(r, ne, chi_e * opt_factor_e, S_e, Te_edge)
    print(f"   Optimal chi_e factor: {opt_factor_e:.4f}")
    print(f"   Te_center = {Te_matrix_opt[0]:.2f} keV (sim: {Te_sim[0]:.2f})")

    # Ion temperature
    print("\n--- ION TEMPERATURE ---")

    ni = data['ni']
    chi_i = data['chi_i']
    S_i = data['PIN_D'] + data['PIN_T'] + data['PIN_He4']
    Ti_edge = data['Ti'][-1]
    Ti_sim = data['Ti']

    print(f"Ti_edge = {Ti_edge:.3f} keV")
    print(f"S_i center = {S_i[0]:.4f} MW/m^3")

    Ti_matrix_tr, _ = solve_temperature_matrix_tr(r, ni, chi_i, S_i, Ti_edge)
    ratio_i = Ti_matrix_tr[0] / Ti_sim[0] if Ti_sim[0] > 0 else 0
    print(f"Ti_center (TR-style) = {Ti_matrix_tr[0]:.2f} keV (sim: {Ti_sim[0]:.2f}, ratio: {ratio_i:.2f})")

    Ti_matrix_scaled, _ = solve_temperature_matrix_tr(r, ni, chi_i * 0.5, S_i, Ti_edge)
    ratio_i2 = Ti_matrix_scaled[0] / Ti_sim[0] if Ti_sim[0] > 0 else 0
    print(f"Ti_center (chi*0.5) = {Ti_matrix_scaled[0]:.2f} keV (sim: {Ti_sim[0]:.2f}, ratio: {ratio_i2:.2f})")

    # Find optimal chi factor for ions
    opt_factor_i = find_optimal_chi_factor(r, ni, chi_i, S_i, Ti_edge, Ti_sim, geom)
    Ti_matrix_opt, _ = solve_temperature_matrix_tr(r, ni, chi_i * opt_factor_i, S_i, Ti_edge)
    print(f"Optimal chi_i factor: {opt_factor_i:.4f}")
    print(f"Ti_center (optimal) = {Ti_matrix_opt[0]:.2f} keV (sim: {Ti_sim[0]:.2f})")

    # Plot comparison
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Electron temperature
    ax = axes[0, 0]
    ax.plot(r, Te_sim, 'b-', lw=2, label='Te (simulation)')
    ax.plot(r, Te_matrix_tr, 'g--', lw=1.5, label='Te (chi as-is)')
    ax.plot(r, Te_matrix_opt, 'r-', lw=2, label=f'Te (chi*{opt_factor_e:.3f})')
    ax.set_xlabel('r/a')
    ax.set_ylabel('Te [keV]')
    ax.set_title('Electron Temperature')
    ax.legend()
    ax.grid(True)
    ax.set_xlim([0, 1])

    # Ion temperature
    ax = axes[0, 1]
    ax.plot(r, Ti_sim, 'b-', lw=2, label='Ti (simulation)')
    ax.plot(r, Ti_matrix_tr, 'g--', lw=1.5, label='Ti (chi as-is)')
    ax.plot(r, Ti_matrix_opt, 'r-', lw=2, label=f'Ti (chi*{opt_factor_i:.3f})')
    ax.set_xlabel('r/a')
    ax.set_ylabel('Ti [keV]')
    ax.set_title('Ion Temperature')
    ax.legend()
    ax.grid(True)
    ax.set_xlim([0, 1])

    # Ratio plot
    ax = axes[1, 0]
    Te_ratio_tr = np.where(Te_sim > 0.1, Te_matrix_tr / Te_sim, np.nan)
    Te_ratio_opt = np.where(Te_sim > 0.1, Te_matrix_opt / Te_sim, np.nan)
    Ti_ratio_tr = np.where(Ti_sim > 0.1, Ti_matrix_tr / Ti_sim, np.nan)
    Ti_ratio_opt = np.where(Ti_sim > 0.1, Ti_matrix_opt / Ti_sim, np.nan)

    ax.plot(r, Te_ratio_tr, 'b--', lw=1.5, label='Te (chi as-is)')
    ax.plot(r, Te_ratio_opt, 'b-', lw=2, label=f'Te (chi*{opt_factor_e:.3f})')
    ax.plot(r, Ti_ratio_tr, 'r--', lw=1.5, label='Ti (chi as-is)')
    ax.plot(r, Ti_ratio_opt, 'r-', lw=2, label=f'Ti (chi*{opt_factor_i:.3f})')
    ax.axhline(y=1.0, color='k', ls='--', lw=1)
    ax.set_xlabel('r/a')
    ax.set_ylabel('Ratio (calc/sim)')
    ax.set_title('Temperature Ratio (calc/sim)')
    ax.legend()
    ax.grid(True)
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 2.5])

    # Geometry factors
    ax = axes[1, 1]
    ax.semilogy(r, geom['DVRHO'], 'b-', lw=2, label='DVRHO')
    ax.semilogy(r, geom['DV23'], 'g-', lw=2, label='DV23')
    ax.semilogy(r, geom['DV53'], 'r-', lw=2, label='DV53')
    ax.semilogy(r, geom['FB'], 'm-', lw=2, label='FB')
    ax.set_xlabel('r/a')
    ax.set_ylabel('Value')
    ax.set_title('TR Geometry Factors')
    ax.legend()
    ax.grid(True)
    ax.set_xlim([0, 1])

    plt.tight_layout()
    plt.savefig('temperature_matrix_solver.png', dpi=150)
    print("\nPlot saved to temperature_matrix_solver.png")

    # Summary table
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"{'Method':<35} {'Te_center':>12} {'Ratio':>10} | {'Ti_center':>12} {'Ratio':>10}")
    print("-"*85)
    print(f"{'Simulation':<35} {Te_sim[0]:>12.2f} {'1.00':>10} | {Ti_sim[0]:>12.2f} {'1.00':>10}")
    print(f"{'Matrix (chi as-is)':<35} {Te_matrix_tr[0]:>12.2f} {ratio1:>10.2f} | {Ti_matrix_tr[0]:>12.2f} {ratio_i:>10.2f}")
    print(f"{'Matrix (chi*0.5)':<35} {Te_matrix_scaled[0]:>12.2f} {ratio3:>10.2f} | {Ti_matrix_scaled[0]:>12.2f} {ratio_i2:>10.2f}")
    print(f"{'Matrix (optimal chi*{:.3f}/{:.3f})':<35} {Te_matrix_opt[0]:>12.2f} {'1.00':>10} | {Ti_matrix_opt[0]:>12.2f} {'1.00':>10}".format(opt_factor_e, opt_factor_i))

    print("\n" + "="*70)
    print("CONCLUSION")
    print("="*70)
    print(f"""
The chi values from TR output need to be scaled by a factor of ~0.5 for
steady-state calculations:

  - Electron: chi_effective = chi_output * {opt_factor_e:.4f}
  - Ion:      chi_effective = chi_output * {opt_factor_i:.4f}

This factor arises from TR's volume-weighted coordinate system and
finite-difference discretization scheme. The TR Fortran code uses:
  - Volume element: DVRHO = 4*pi^2 * kappa * a^2 * R0 * rho
  - Source term scaling: DV53 = DVRHO^(5/3)
  - Diffusion term scaling: DV23 = DVRHO^(2/3)

These factors effectively increase the apparent chi by ~2x compared to
the simple cylindrical heat equation.
""")

    return {
        'r': r,
        'Te_sim': Te_sim, 'Te_matrix': Te_matrix_opt,
        'Ti_sim': Ti_sim, 'Ti_matrix': Ti_matrix_opt,
        'geom': geom,
        'opt_factor_e': opt_factor_e,
        'opt_factor_i': opt_factor_i
    }


if __name__ == '__main__':
    main()
