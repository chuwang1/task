#!/usr/bin/env python3

import os
import glob
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


BASE = "/Users/dengxiaoya/TASK/CFEDR/git/task/trx.260120"
RKEV = 1.602176487e-16
RMU0 = 4.0e-7 * np.pi
DEFAULT_INPUT = os.path.join(BASE, "tr.CFEDR0114Test1.in")


def find_csv_by_title(title_keyword, search_dir=BASE):
    for f in sorted(glob.glob(os.path.join(search_dir, "tr_data_*.csv"))):
        with open(f, "r") as fh:
            title = fh.readline().strip()
        if title_keyword in title:
            return f, title
    return None, None


def read_tr_profile_csv(filename, prefix=None, time_col=-1):
    df = pd.read_csv(filename, skiprows=1)
    df.columns = df.columns.str.strip()
    x = df["X"].values
    if prefix is None:
        return x, df
    cols = [c for c in df.columns if c.startswith(prefix + "_")]
    if cols:
        cols = sorted(cols, key=lambda c: int(c.split("_")[-1]))
        y = df[cols[time_col]].values
    elif prefix in df.columns:
        y = df[prefix].values
    else:
        raise ValueError(f"No column for {prefix} in {filename}")
    return x, y


def read_input_params(input_path=DEFAULT_INPUT):
    text = open(input_path).read()
    vals = {}
    for key in ["RA", "RR", "BB"]:
        m = re.search(rf'\b{key}\s*=\s*([0-9.Ee+-]+)', text)
        if not m:
            raise ValueError(f"Missing {key} in {input_path}")
        vals[key] = float(m.group(1))
    return vals


def fit_rhog_scale():
    path, _ = find_csv_by_title('@RHOM,RHOG  vs RHO@')
    if path is None:
        return None
    x, df = read_tr_profile_csv(path, None)
    y = df['RHOG'].values
    mask = x > 0
    a = np.dot(x[mask], y[mask]) / np.dot(x[mask], x[mask])
    return a


def main():
    vals = read_input_params()
    RA, RR, BB = vals["RA"], vals["RR"], vals["BB"]

    q_path, _ = find_csv_by_title('@QP  vs r@')
    dens_path, _ = find_csv_by_title('@n(NS)')
    temp_path, _ = find_csv_by_title('@T(NS)')
    alpha_path, _ = find_csv_by_title('@alpha  vs r@')
    if None in [q_path, dens_path, temp_path, alpha_path]:
        raise FileNotFoundError('Required CSV files not found.')

    x_q, q = read_tr_profile_csv(q_path, 'QP')
    x_alpha, alpha_csv = read_tr_profile_csv(alpha_path, 'alpha')
    _, dens_df = read_tr_profile_csv(dens_path, None)
    _, temp_df = read_tr_profile_csv(temp_path, None)

    # thermal pressure on density/temperature grid (10^20 m^-3 * keV)
    x_p = dens_df['X'].values
    p_th = (dens_df['nE'].values * temp_df['TE'].values +
            dens_df['nD'].values * temp_df['TD'].values +
            dens_df['nT'].values * temp_df['TT'].values +
            dens_df['nA'].values * temp_df['TA'].values)
    q_on_p = np.interp(x_p, x_q, q)

    dp_d_ra = np.gradient(p_th, RA * x_p, edge_order=2) * 1.0e20 * RKEV
    alpha_from_ra = -2.0 * RMU0 * q_on_p**2 * RR / BB**2 * dp_d_ra

    rhog_scale = fit_rhog_scale()
    if rhog_scale is None:
        rhog_scale = 3.6146142
    dp_d_rhog = np.gradient(p_th, rhog_scale * x_p, edge_order=2) * 1.0e20 * RKEV
    alpha_from_rhog = -2.0 * RMU0 * q_on_p**2 * RR / BB**2 * dp_d_rhog

    alpha_ra_on_csv = np.interp(x_alpha, x_p, alpha_from_ra)
    alpha_rhog_on_csv = np.interp(x_alpha, x_p, alpha_from_rhog)

    rms_ra = np.sqrt(np.mean((alpha_csv - alpha_ra_on_csv)**2))
    rms_rhog = np.sqrt(np.mean((alpha_csv - alpha_rhog_on_csv)**2))

    out = pd.DataFrame({
        'X': x_alpha,
        'alpha_csv': alpha_csv,
        'alpha_from_ra': alpha_ra_on_csv,
        'alpha_from_rhog': alpha_rhog_on_csv,
        'err_ra': alpha_csv - alpha_ra_on_csv,
        'err_rhog': alpha_csv - alpha_rhog_on_csv,
    })
    out_csv = os.path.join(BASE, 'alpha_reconstruction_test.csv')
    out.to_csv(out_csv, index=False)

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    ax = axes[0, 0]
    ax.plot(x_alpha, alpha_csv, 'k-', lw=2.0, label='alpha from TR CSV')
    ax.plot(x_alpha, alpha_ra_on_csv, 'r--', lw=2.0, label=r'$\alpha$ from $dp/d(RA\,X)$')
    ax.plot(x_alpha, alpha_rhog_on_csv, 'b-.', lw=2.0, label=r'$\alpha$ from $dp/d(RHOG)$')
    ax.set_xlabel('TR X')
    ax.set_ylabel(r'$\alpha$')
    ax.set_title('Alpha reconstruction')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    ax.plot(x_alpha, alpha_csv - alpha_ra_on_csv, 'r-', lw=2.0, label='err from RA*X')
    ax.plot(x_alpha, alpha_csv - alpha_rhog_on_csv, 'b-', lw=2.0, label='err from RHOG')
    ax.axhline(0.0, color='k', lw=1.0, alpha=0.5)
    ax.set_xlabel('TR X')
    ax.set_ylabel('error')
    ax.set_title('Alpha reconstruction error')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    ax = axes[1, 0]
    ax.plot(x_p, p_th, 'g-', lw=2.0)
    ax.set_xlabel('TR X')
    ax.set_ylabel(r'$p_{th}$ [$10^{20}$ m$^{-3}$ keV]')
    ax.set_title('Thermal pressure used')
    ax.grid(True, alpha=0.3)

    ax = axes[1, 1]
    ax.plot(x_q, q, 'm-', lw=2.0, label='q')
    ax2 = ax.twinx()
    ax2.plot(x_p, p_th, 'c--', lw=1.5, label='p_th')
    ax.set_xlabel('TR X')
    ax.set_ylabel('q', color='m')
    ax2.set_ylabel(r'$p_{th}$', color='c')
    ax.set_title('Inputs for alpha')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    fig_out = os.path.join(BASE, 'alpha_reconstruction_test.png')
    plt.savefig(fig_out, dpi=170)
    plt.close(fig)

    print(f'Wrote {fig_out}')
    print(f'Wrote {out_csv}')
    print(f'Fitted RHOG scale ~ {rhog_scale:.6f}')
    print(f'RMS error using d/d(RA*X)  : {rms_ra:.6e}')
    print(f'RMS error using d/d(RHOG)  : {rms_rhog:.6e}')
    if rms_rhog < rms_ra:
        print('Conclusion: alpha is better reproduced using RHOG-scaled radius.')
    else:
        print('Conclusion: alpha is better reproduced using RA*X radius.')


if __name__ == '__main__':
    main()
