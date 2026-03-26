#!/usr/bin/env python3
"""Strict TR heat-transport comparison based on the Fortran matrix form.

This script does NOT use `verify_heat_flux.py` logic.
It reconstructs the energy-transport balance directly from the TR code form:

For NSVN = 2 in `trexec.f90`:
    DD(NEQ,NEQ  ) = FB * AK * DV23
    DD(NEQ,NEQ-1) = FB * (AD*CC - AK) * T * DV23
    D_rhs_source  = (PIN / (RKEV*1e20)) * DV53

with
    FB   = DVRHO * AR2RHO / DR^2
    DV23 = DVRHO^(2/3)
    DV53 = DVRHO^(5/3)

Using available output files, this script compares:
  1) Source-integrated flux variable from PIN and DV53
  2) Gradient flux variable reconstructed from pressure gradient term only
  3) Gradient flux variable reconstructed from temperature-gradient form

The off-diagonal density-coupling term involving AD is not included because AD
is not currently exported in the standard CSV set for this case.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path


R0 = 8.03
A_MINOR = 2.68
KAPPA = 1.89
RKEV = 1.6022e-16


def read_tr_csv(path: str) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    df = pd.read_csv(path, skiprows=1)
    df.columns = df.columns.str.strip()
    x = df.iloc[:, 0].to_numpy()
    data = {col: df[col].to_numpy() for col in df.columns[1:]}
    return x, data


def find_csv_by_title_fragment(fragment: str) -> str:
    for path in sorted(Path('.').glob('tr_data_*.csv')):
        with path.open() as f:
            title = f.readline().strip()
        if fragment in title:
            return str(path)
    raise FileNotFoundError(f'Could not find CSV title containing: {fragment}')


def get_col(data: dict[str, np.ndarray], names: list[str], x: np.ndarray) -> np.ndarray:
    for name in names:
        if name in data:
            return data[name]
    return np.zeros_like(x)


def cumulative_trapezoid(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    out = np.zeros_like(x, dtype=float)
    for i in range(1, len(x)):
        dx = x[i] - x[i - 1]
        out[i] = out[i - 1] + 0.5 * (y[i - 1] + y[i]) * dx
    return out


def main() -> None:
    # Read TR outputs
    rho, n_data = read_tr_csv(find_csv_by_title_fragment('@n(NS)'))
    _, t_data = read_tr_csv(find_csv_by_title_fragment('@T(NS)'))
    _, ake_data = read_tr_csv(find_csv_by_title_fragment('@AKE,AKNCE,AKDWE'))
    _, akd_data = read_tr_csv(find_csv_by_title_fragment('@AKD,AKNCD,AKDWD'))
    _, ad_data = read_tr_csv(find_csv_by_title_fragment('@AD(NS)'))
    _, pin_data = read_tr_csv(find_csv_by_title_fragment('@PIN [MW/m$+3$=]'))

    # Profiles in TR code units
    nE = get_col(n_data, ['nE', 'NE'], rho)  # 1e20 / m^3
    nD = get_col(n_data, ['nD', 'ND'], rho)
    nT = get_col(n_data, ['nT', 'NT'], rho)
    nA = get_col(n_data, ['nA', 'NA'], rho)

    TE = get_col(t_data, ['TE', 'Te'], rho)  # keV
    TD = get_col(t_data, ['TD', 'Ti'], rho)
    TT = get_col(t_data, ['TT'], rho)
    TA = get_col(t_data, ['TA'], rho)

    AKE = get_col(ake_data, ['AKE'], rho)
    AKD = get_col(akd_data, ['AKD'], rho)
    ADE = get_col(ad_data, ['ADE'], rho)
    ADD = get_col(ad_data, ['ADD'], rho)
    ADT = get_col(ad_data, ['ADT'], rho)
    ADA = get_col(ad_data, ['ADA'], rho)

    PIN_e = get_col(pin_data, ['PIN_1'], rho)             # MW/m^3
    PIN_D = get_col(pin_data, ['PIN_2'], rho)
    PIN_T = get_col(pin_data, ['PIN_3'], rho)
    PIN_A = get_col(pin_data, ['PIN_4'], rho)
    PIN_i = PIN_D + PIN_T + PIN_A

    # Geometry factors used in TR
    rkaps = np.sqrt(KAPPA)
    ar2rho = 1.0 / (rkaps * A_MINOR) ** 2
    dvrho = 4.0 * np.pi**2 * R0 * KAPPA * A_MINOR**2 * rho
    dv23 = dvrho ** (2.0 / 3.0)
    dv53 = dvrho ** (5.0 / 3.0)

    # Pressure variables in code units: n[1e20/m^3] * T[keV]
    p_e = nE * TE
    p_i = nD * TD + nT * TT + nA * TA

    dp_e_drho = np.gradient(p_e, rho)
    dp_i_drho = np.gradient(p_i, rho)

    dT_e_drho = np.gradient(TE, rho)
    dT_i_drho = np.gradient(0.5 * (TD + TT), rho)

    # Strict source-side flux variable from the RHS source term in trexec.f90:1599,1643
    # F_src(ρ) = ∫ [ PIN/(RKEV*1e20) * DV53 ] dρ
    Fsrc_e = cumulative_trapezoid(rho, (PIN_e * 1e6 / (RKEV * 1.0e20)) * dv53)
    Fsrc_i = cumulative_trapezoid(rho, (PIN_i * 1e6 / (RKEV * 1.0e20)) * dv53)

    # Gradient-side flux variable from the diagonal pressure-gradient term only
    # Fgrad_p = DVRHO * AR2RHO * AK * DV23 * d(nT)/dρ
    Fgrad_p_e = dvrho * ar2rho * AKE * dv23 * dp_e_drho
    Fgrad_p_i = dvrho * ar2rho * AKD * dv23 * dp_i_drho

    # Add trexec.f90:1418 off-diagonal density-coupling term:
    # DD(NEQ,NEQ-1)= FB*(AD*CC - AK)*RTV*DV23
    # This combines with the diagonal pressure-gradient term. Using
    # dp/dρ = n dT/dρ + T dn/dρ, the total reconstructed flux becomes
    # proportional to AK*n*dT/dρ + CC*AD*T*dn/dρ.
    CC = 1.5
    dnE_drho = np.gradient(nE, rho)
    dnD_drho = np.gradient(nD, rho)
    dnT_drho = np.gradient(nT, rho)
    dnA_drho = np.gradient(nA, rho)
    Fgrad_strict_e = dvrho * ar2rho * dv23 * (AKE * nE * np.gradient(TE, rho) + CC * ADE * TE * dnE_drho)
    Fgrad_strict_i = dvrho * ar2rho * dv23 * (
        AKD * (nD * np.gradient(TD, rho) + nT * np.gradient(TT, rho) + nA * np.gradient(TA, rho))
        + CC * (ADD * TD * dnD_drho + ADT * TT * dnT_drho + ADA * TA * dnA_drho)
    )

    # Also show the simpler temperature-gradient form for reference
    Fgrad_T_e = dvrho * ar2rho * AKE * dv23 * nE * dT_e_drho
    Fgrad_T_i = dvrho * ar2rho * AKD * dv23 * (nD + nT) * dT_i_drho

    # Sign convention: source-integrated flux should oppose outward gradient term
    # compare with -Fgrad_*
    fig, axes = plt.subplots(2, 2, figsize=(14, 9), sharex=True)

    ax = axes[0, 0]
    ax.plot(rho, Fsrc_e, 'k-', lw=2.2, label='Fsrc_e from PIN*DV53 integral')
    ax.plot(rho, -Fgrad_p_e, 'r--', lw=2.0, label='-Fgrad_p_e (pressure gradient)')
    ax.plot(rho, -Fgrad_strict_e, 'g-.', lw=2.0, label='-Fgrad_strict_e (+AD density term)')
    ax.plot(rho, -Fgrad_T_e, 'b:', lw=2.0, label='-Fgrad_T_e (temperature gradient)')
    ax.set_title('Electron Strict Flux Comparison')
    ax.set_ylabel('Code flux variable')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9)

    ax = axes[0, 1]
    ratio_p_e = np.divide(Fsrc_e, -Fgrad_p_e, out=np.full_like(Fsrc_e, np.nan), where=np.abs(Fgrad_p_e) > 1e-12)
    ratio_s_e = np.divide(Fsrc_e, -Fgrad_strict_e, out=np.full_like(Fsrc_e, np.nan), where=np.abs(Fgrad_strict_e) > 1e-12)
    ratio_t_e = np.divide(Fsrc_e, -Fgrad_T_e, out=np.full_like(Fsrc_e, np.nan), where=np.abs(Fgrad_T_e) > 1e-12)
    ax.plot(rho, ratio_p_e, 'r--', lw=2.0, label='Fsrc / (-Fgrad_p)')
    ax.plot(rho, ratio_s_e, 'g-.', lw=2.0, label='Fsrc / (-Fgrad_strict)')
    ax.plot(rho, ratio_t_e, 'b:', lw=2.0, label='Fsrc / (-Fgrad_T)')
    ax.axhline(1.0, color='k', lw=1)
    ax.set_title('Electron Ratio to Source Integral')
    ax.set_ylabel('Ratio')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9)

    ax = axes[1, 0]
    ax.plot(rho, Fsrc_i, 'k-', lw=2.2, label='Fsrc_i from PIN*DV53 integral')
    ax.plot(rho, -Fgrad_p_i, 'r--', lw=2.0, label='-Fgrad_p_i (pressure gradient)')
    ax.plot(rho, -Fgrad_strict_i, 'g-.', lw=2.0, label='-Fgrad_strict_i (+AD density term)')
    ax.plot(rho, -Fgrad_T_i, 'b:', lw=2.0, label='-Fgrad_T_i (temperature gradient)')
    ax.set_title('Ion Strict Flux Comparison')
    ax.set_xlabel('rho')
    ax.set_ylabel('Code flux variable')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9)

    ax = axes[1, 1]
    ratio_p_i = np.divide(Fsrc_i, -Fgrad_p_i, out=np.full_like(Fsrc_i, np.nan), where=np.abs(Fgrad_p_i) > 1e-12)
    ratio_s_i = np.divide(Fsrc_i, -Fgrad_strict_i, out=np.full_like(Fsrc_i, np.nan), where=np.abs(Fgrad_strict_i) > 1e-12)
    ratio_t_i = np.divide(Fsrc_i, -Fgrad_T_i, out=np.full_like(Fsrc_i, np.nan), where=np.abs(Fgrad_T_i) > 1e-12)
    ax.plot(rho, ratio_p_i, 'r--', lw=2.0, label='Fsrc / (-Fgrad_p)')
    ax.plot(rho, ratio_s_i, 'g-.', lw=2.0, label='Fsrc / (-Fgrad_strict)')
    ax.plot(rho, ratio_t_i, 'b:', lw=2.0, label='Fsrc / (-Fgrad_T)')
    ax.axhline(1.0, color='k', lw=1)
    ax.set_title('Ion Ratio to Source Integral')
    ax.set_xlabel('rho')
    ax.set_ylabel('Ratio')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9)

    fig.tight_layout()
    fig.savefig('heat_flux_strict_comparison.png', dpi=160)

    out = pd.DataFrame({
        'rho': rho,
        'Fsrc_e': Fsrc_e,
        'Fgrad_p_e': Fgrad_p_e,
        'Fgrad_T_e': Fgrad_T_e,
        'Fgrad_strict_e': Fgrad_strict_e,
        'ratio_p_e': ratio_p_e,
        'ratio_s_e': ratio_s_e,
        'ratio_t_e': ratio_t_e,
        'Fsrc_i': Fsrc_i,
        'Fgrad_p_i': Fgrad_p_i,
        'Fgrad_T_i': Fgrad_T_i,
        'Fgrad_strict_i': Fgrad_strict_i,
        'ratio_p_i': ratio_p_i,
        'ratio_s_i': ratio_s_i,
        'ratio_t_i': ratio_t_i,
        'ADE': ADE,
        'ADD': ADD,
        'ADT': ADT,
        'ADA': ADA,
        'dvrho': dvrho,
        'dv23': dv23,
        'p_e': p_e,
        'p_i': p_i,
        'dp_e_drho': dp_e_drho,
        'dp_i_drho': dp_i_drho,
    })
    out.to_csv('heat_flux_strict_comparison.csv', index=False)

    core = (rho >= 0.1) & (rho <= 0.9)
    print('Saved: heat_flux_strict_comparison.png')
    print('Saved: heat_flux_strict_comparison.csv')
    print('Core-region mean ratios:')
    print(f"  Electron pressure-gradient ratio mean = {np.nanmean(ratio_p_e[core]):.4f}")
    print(f"  Electron strict(+AD) ratio mean = {np.nanmean(ratio_s_e[core]):.4f}")
    print(f"  Electron temperature-gradient ratio mean = {np.nanmean(ratio_t_e[core]):.4f}")
    print(f"  Ion pressure-gradient ratio mean = {np.nanmean(ratio_p_i[core]):.4f}")
    print(f"  Ion strict(+AD) ratio mean = {np.nanmean(ratio_s_i[core]):.4f}")
    print(f"  Ion temperature-gradient ratio mean = {np.nanmean(ratio_t_i[core]):.4f}")


if __name__ == '__main__':
    main()
