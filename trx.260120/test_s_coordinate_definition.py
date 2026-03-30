#!/usr/bin/env python3

import os
import glob
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from plot_sumavir2_from_gfile import read_input_profiles_mapping


BASE = "/Users/dengxiaoya/TASK/CFEDR/git/task/trx.260120"
OMFIT_OUT = "/Users/dengxiaoya/CFEDRSW/OMFIT_out"
INPUT_PROFILES = os.path.join(OMFIT_OUT, "input.profiles")
DEFAULT_INPUT = os.path.join(BASE, "tr.CFEDR0114Test1.in")


def find_csv_by_title(title_keyword, search_dir=BASE):
    for f in sorted(glob.glob(os.path.join(search_dir, "tr_data_*.csv"))):
        with open(f, "r") as fh:
            title = fh.readline().strip()
        if title_keyword in title:
            return f, title
    return None, None


def read_tr_profile_csv(filename, prefix, time_col=-1):
    df = pd.read_csv(filename, skiprows=1)
    df.columns = df.columns.str.strip()
    x = df["X"].values
    cols = [c for c in df.columns if c.startswith(prefix + "_")]
    if cols:
        cols = sorted(cols, key=lambda c: int(c.split("_")[-1]))
        y = df[cols[time_col]].values
    elif prefix in df.columns:
        y = df[prefix].values
    else:
        raise ValueError(f"No column for {prefix} in {filename}")
    return x, y


def deriv_like_tr(x, y):
    n = len(x)
    out = np.zeros_like(y)
    dx = x[1] - x[0]
    out[0] = (4.0 * y[1] - 3.0 * y[0] - y[2]) / (2.0 * dx)
    out[-1] = (3.0 * y[-1] - 4.0 * y[-2] + y[-3]) / (2.0 * dx)
    out[1:-1] = (y[2:] - y[:-2]) / (2.0 * dx)
    return out


def analytic_q_profile(rho, q0=1.05, q95=3.5):
    """Same analytic q model used in calculate_cdbm_chi.py fallback."""
    return q0 + (q95 - q0) * rho**2 * (1.0 + 0.5 * rho**2)


def read_ra_rkap(input_path=DEFAULT_INPUT):
    if not os.path.exists(input_path):
        return None, None
    text = open(input_path).read()
    m_ra = re.search(r'\bRA\s*=\s*([0-9.Ee+-]+)', text)
    m_k = re.search(r'\bRKAP\s*=\s*([0-9.Ee+-]+)', text)
    if not m_ra or not m_k:
        return None, None
    return float(m_ra.group(1)), float(m_k.group(1))


def main():
    s_csv_path, s_title = find_csv_by_title('@s  vs r@')
    q_csv_path, q_title = find_csv_by_title('@QP  vs r@')
    if s_csv_path is None or q_csv_path is None:
        raise FileNotFoundError('Required CSV files not found.')

    x_s, s_csv = read_tr_profile_csv(s_csv_path, 's')
    x_q, q_csv = read_tr_profile_csv(q_csv_path, 'QP')
    q_on_x = np.interp(x_s, x_q, q_csv)
    arho_exp = 3.6146142

    # Reconstruction on TR output coordinate X
    dq_dX = deriv_like_tr(x_s, q_on_x)
    s_from_x = x_s / q_on_x * dq_dX

    # Reconstruction on r/a using mapping from input.profiles
    psi_n_map, rhoa_map = read_input_profiles_mapping(INPUT_PROFILES)
    sqrt_psin = np.sqrt(np.clip(psi_n_map, 0.0, None))
    # Treat TR X as sqrt(psit_n) ~ flux coordinate; map to r/a through input.profiles sqrt(psin)
    rhoa_on_x = np.interp(x_s, sqrt_psin, rhoa_map)
    # Since mapped r/a is nonuniform, use coordinate-aware gradient
    dq_drhoa = np.gradient(q_on_x, rhoa_on_x, edge_order=2)
    s_from_rhoa = rhoa_on_x / q_on_x * dq_drhoa

    # Reconstruction using ARHO_EXP * X as effective radius scale
    s_from_arhoexp = arho_exp * x_s / q_on_x * dq_dX

    # More faithful Fortran-like reconstruction:
    #   S(NR) = RHOG(NR)/Q * DQ,  with DQ ~ dq/d(rho_flux)
    # In TRSTGF (default metric), RHOG = RG / RJCB = RG * sqrt(RKAP) * RA.
    ra_in, rkap_in = read_ra_rkap()
    if ra_in is not None and rkap_in is not None:
        rhog_factor = np.sqrt(rkap_in) * ra_in
        dq_dflux_stag = (q_csv[1:] - q_csv[:-1]) / (x_q[1] - x_q[0])
        q_half = 0.5 * (q_csv[:-1] + q_csv[1:])
        rhog_half = rhog_factor * x_q[1:]
        s_fortran_like = rhog_half / q_half * dq_dflux_stag
    else:
        rhog_factor = np.nan
        s_fortran_like = np.full_like(s_csv, np.nan)

    # Analytic-model shear from calculate_cdbm_chi.py fallback
    q_analytic = analytic_q_profile(x_s)
    dqdr_analytic = np.gradient(q_analytic, x_s)
    s_analytic = np.zeros_like(x_s)
    s_analytic[1:] = (x_s[1:] / q_analytic[1:]) * dqdr_analytic[1:]
    s_analytic[0] = s_analytic[1]

    rms_x = np.sqrt(np.mean((s_csv - s_from_x) ** 2))
    rms_rhoa = np.sqrt(np.mean((s_csv - s_from_rhoa) ** 2))
    rms_arhoexp = np.sqrt(np.mean((s_csv - s_from_arhoexp) ** 2))
    rms_analytic = np.sqrt(np.mean((s_csv - s_analytic) ** 2))
    rms_fortran = np.sqrt(np.nanmean((s_csv - s_fortran_like) ** 2))

    out = pd.DataFrame({
        'X_TR': x_s,
        'q': q_on_x,
        's_csv': s_csv,
        's_recon_from_X': s_from_x,
        'r_over_a_from_mapping': rhoa_on_x,
        's_recon_from_r_over_a': s_from_rhoa,
        's_recon_from_ARHO_EXP_X': s_from_arhoexp,
        'q_analytic': q_analytic,
        's_analytic_model': s_analytic,
        's_fortran_like': s_fortran_like,
        'err_X': s_csv - s_from_x,
        'err_rhoa': s_csv - s_from_rhoa,
        'err_arhoexp': s_csv - s_from_arhoexp,
        'err_analytic': s_csv - s_analytic,
        'err_fortran_like': s_csv - s_fortran_like,
    })
    csv_out = os.path.join(BASE, 's_coordinate_test.csv')
    out.to_csv(csv_out, index=False)

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    ax = axes[0, 0]
    ax.plot(x_s, s_csv, 'k-', lw=2.0, label='s from TR CSV')
    ax.plot(x_s, s_from_x, 'r--', lw=2.0, label=r'$s=(X/q)dq/dX$')
    ax.plot(x_s, s_from_rhoa, 'b-.', lw=2.0, label=r'$s=((r/a)/q)dq/d(r/a)$')
    ax.plot(x_s, s_from_arhoexp, color='tab:purple', lw=2.0,
            label=r'$s=(A_{\rho}X/q)\,dq/dX$')
    ax.plot(x_s, s_analytic, 'g:', lw=2.2, label='analytic shear model')
    if np.isfinite(rms_fortran):
        ax.plot(x_s, s_fortran_like, color='tab:orange', lw=2.0, label='Fortran-like S=RHOG/Q*DQ')
    ax.set_xlabel('TR X')
    ax.set_ylabel('s')
    ax.set_title('Magnetic shear reconstruction')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    ax.plot(x_s, s_csv - s_from_x, 'r-', lw=2.0, label='err from X')
    ax.plot(x_s, s_csv - s_from_rhoa, 'b-', lw=2.0, label='err from r/a')
    ax.plot(x_s, s_csv - s_from_arhoexp, color='tab:purple', lw=2.0, label='err from ARHO_EXP*X')
    ax.plot(x_s, s_csv - s_analytic, 'g-', lw=2.0, label='err from analytic model')
    if np.isfinite(rms_fortran):
        ax.plot(x_s, s_csv - s_fortran_like, color='tab:orange', lw=2.0, label='err from Fortran-like')
    ax.axhline(0.0, color='k', lw=1.0, alpha=0.5)
    ax.set_xlabel('TR X')
    ax.set_ylabel('error')
    ax.set_title('Reconstruction error')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    ax = axes[1, 0]
    ax.plot(x_s, q_on_x, 'g-', lw=2.0, label='q from TR CSV')
    ax.plot(x_s, q_analytic, 'm--', lw=2.0, label='analytic q model')
    ax.set_xlabel('TR X')
    ax.set_ylabel('q')
    ax.set_title('Safety factor')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    ax = axes[1, 1]
    ax.plot(x_s, rhoa_on_x, 'm-', lw=2.0, label='mapped r/a')
    ax.plot(x_s, x_s, 'k--', lw=1.0, label='TR X')
    ax.set_xlabel('TR X')
    ax.set_ylabel('coordinate value')
    ax.set_title('TR X to r/a mapping')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    plt.tight_layout()
    fig_out = os.path.join(BASE, 's_coordinate_test.png')
    plt.savefig(fig_out, dpi=170)
    plt.close(fig)

    print(f's CSV : {os.path.basename(s_csv_path)}')
    print(f'q CSV : {os.path.basename(q_csv_path)}')
    print(f'Wrote {fig_out}')
    print(f'Wrote {csv_out}')
    print(f'RMS error using TR X        : {rms_x:.6e}')
    print(f'RMS error using mapped r/a  : {rms_rhoa:.6e}')
    print(f'RMS error using ARHO_EXP*X  : {rms_arhoexp:.6e}')
    print(f'RMS error using analytic q/s: {rms_analytic:.6e}')
    print(f'RMS error using Fortran-like: {rms_fortran:.6e}')
    if np.isfinite(rhog_factor):
        print(f'Fortran-like RHOG factor    : sqrt(RKAP)*RA = {rhog_factor:.6f}')
    if rms_x < rms_rhoa:
        print('Conclusion: s in TR CSV is more consistent with TR output coordinate X.')
    else:
        print('Conclusion: s in TR CSV is more consistent with mapped r/a.')


if __name__ == '__main__':
    main()
