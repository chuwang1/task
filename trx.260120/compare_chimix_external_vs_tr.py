#!/usr/bin/env python3

import os
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


BASE = os.path.dirname(os.path.abspath(__file__))
CHI_FILE = os.path.join(BASE, "chi_div_grad2_drdrho.dat")


def find_csv_by_title(title_keyword):
    matched = []
    for f in sorted(glob.glob(os.path.join(BASE, "tr_data_*.csv"))):
        with open(f, "r") as fh:
            title = fh.readline().strip()
        if title_keyword in title:
            matched.append((f, title))
    if matched:
        return matched[-1]
    return None, None


def read_csv(path):
    df = pd.read_csv(path, skiprows=1)
    df.columns = df.columns.str.strip()
    return df


def main():
    mix_path, _ = find_csv_by_title('CHIMIXW,AKEXT_E,AKEXT_I')
    ake_path, _ = find_csv_by_title('AKE,AKNCE,AKDWE')
    akd_path, _ = find_csv_by_title('AKD,AKNCD,AKDWD')
    if None in [mix_path, ake_path, akd_path]:
        raise FileNotFoundError('Required TR CSV outputs not found.')

    ext = pd.read_csv(CHI_FILE, sep=r'\s+')
    mix = read_csv(mix_path)
    ake = read_csv(ake_path)
    akd = read_csv(akd_path)

    rho = mix['X'].values
    w = mix['CHIMIXW'].values
    akext_e = mix['AKEXT_E'].values
    akext_i = mix['AKEXT_I'].values

    ext_e = np.interp(rho, ext['rho'].values, ext['chi_e'].values)
    ext_i = np.interp(rho, ext['rho'].values, ext['chi_i'].values)
    akdwe = np.interp(rho, ake['X'].values, ake['AKDWE'].values)
    akdwd = np.interp(rho, akd['X'].values, akd['AKDWD'].values)

    edge = rho >= 0.9
    mid = (rho >= 0.88) & (rho <= 0.92)

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    ax = axes[0, 0]
    ax.plot(rho, w, 'k-', lw=2.0)
    ax.set_xlabel(r'$\rho$')
    ax.set_ylabel('CHIMIXW')
    ax.set_title('Mix weight')
    ax.grid(True, alpha=0.3)

    ax = axes[0, 1]
    ax.plot(rho, ext_e, 'r--', lw=2.0, label='external chi_e file')
    ax.plot(rho, akext_e, 'b-', lw=2.0, label='AKEXT_E in TR')
    ax.plot(rho, akdwe, 'g-', lw=1.5, label='TR AKDWE')
    ax.set_xlabel(r'$\rho$')
    ax.set_ylabel(r'$\chi_e$ [m$^2$/s]')
    ax.set_title('Electron chi comparison')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    ax = axes[1, 0]
    ax.plot(rho, ext_i, 'r--', lw=2.0, label='external chi_i file')
    ax.plot(rho, akext_i, 'b-', lw=2.0, label='AKEXT_I in TR')
    ax.plot(rho, akdwd, 'g-', lw=1.5, label='TR AKDWD')
    ax.set_xlabel(r'$\rho$')
    ax.set_ylabel(r'$\chi_i$ [m$^2$/s]')
    ax.set_title('Ion chi comparison')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    ax = axes[1, 1]
    ax.plot(rho[edge], akdwe[edge], 'g-', lw=1.8, label='TR AKDWE (edge)')
    ax.plot(rho[edge], akext_e[edge], 'b--', lw=2.0, label='AKEXT_E (edge)')
    ax.plot(rho[edge], w[edge], 'k-.', lw=1.8, label='CHIMIXW (edge)')
    ax.set_xlabel(r'$\rho$')
    ax.set_title(r'Edge behavior ($\rho \geq 0.9$)')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    plt.tight_layout()
    fig_out = os.path.join(BASE, 'compare_chimix_external_vs_tr.png')
    plt.savefig(fig_out, dpi=170)
    plt.close(fig)

    print(f'Wrote {fig_out}')
    print(f'mix csv  : {os.path.basename(mix_path)}')
    print(f'e chi csv: {os.path.basename(ake_path)}')
    print(f'i chi csv: {os.path.basename(akd_path)}')
    print(f'max CHIMIXW          = {np.max(w):.6f}')
    print(f'max AKEXT_E in TR    = {np.max(akext_e):.6f}')
    print(f'max AKEXT_I in TR    = {np.max(akext_i):.6f}')
    print(f'max ext chi_e (file) = {np.max(ext_e):.6f}')
    print(f'max ext chi_i (file) = {np.max(ext_i):.6f}')
    print(f'edge mean CHIMIXW    = {np.mean(w[edge]) if np.any(edge) else np.nan:.6f}')
    print(f'transition mean CHIMIXW (0.88-0.92)= {np.mean(w[mid]) if np.any(mid) else np.nan:.6f}')

    if np.max(w) < 1e-12 and np.max(akext_e) < 1e-12:
        print('Conclusion: current TR diagnostic CSV shows NO mixing applied yet (weight and external chi are zero).')
    else:
        print('Conclusion: mixing diagnostics are nonzero; compare edge AKDWE/AKEXT_E curves to assess blend.')


if __name__ == '__main__':
    main()
