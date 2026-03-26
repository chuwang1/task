#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


TRX_DIR = "/Users/dengxiaoya/TASK/CFEDR/git/task/trx.260120"


def read_tr_dvrho(path):
    with open(path) as f:
        lines = f.readlines()
    rho = []
    dvrho = []
    for line in lines[2:]:
        a, b = [float(s) for s in line.strip().split(',')]
        rho.append(a)
        dvrho.append(b)
    return np.array(rho), np.array(dvrho)


def main():
    tr_rho, tr_dvrho = read_tr_dvrho(os.path.join(TRX_DIR, 'tr_data_133.csv'))
    gf = pd.read_csv(os.path.join(TRX_DIR, 'gfile_dvrho.csv'))
    allp = pd.read_csv('/Users/dengxiaoya/CFEDRSW/OMFIT_out/outputs_new/all_profiles.csv')

    psita = abs(allp['polflux'].iloc[-1])
    dpsit_drho = 2.0 * psita * tr_rho

    dvpsit_tr = np.full_like(tr_dvrho, np.nan)
    mask = np.abs(dpsit_drho) > 1e-30
    dvpsit_tr[mask] = tr_dvrho[mask] / dpsit_drho[mask]

    gfile_dvrho_on_tr = np.interp(tr_rho, gf['rho'].values, gf['DVRHO_gfile'].values)
    dvpsit_g = np.full_like(tr_dvrho, np.nan)
    dvpsit_g[mask] = gfile_dvrho_on_tr[mask] / dpsit_drho[mask]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True)

    ax1.plot(tr_rho, dpsit_drho, 'k-', lw=2, label=r'$d\psi_t/d\rho = 2\psi_{ta}\rho$')
    ax1.set_ylabel(r'$d\psi_t/d\rho$')
    ax1.set_title(r'Comparison of $d\psi_t/d\rho$ and $dV/d\psi_t$')
    ax1.grid(True, alpha=0.3)
    ax1.legend()

    ax2.plot(tr_rho, dvpsit_tr, 'o-', lw=2, ms=4, label=r'TR inferred $dV/d\psi_t$')
    ax2.plot(tr_rho, dvpsit_g, 's--', lw=1.8, ms=4, label=r'gfile inferred $dV/d\psi_t$')
    ax2.set_xlabel(r'$\rho$')
    ax2.set_ylabel(r'$dV/d\psi_t$')
    ax2.grid(True, alpha=0.3)
    ax2.legend()

    out = os.path.join(TRX_DIR, 'dvpsit_dpsitdrho_comparison.png')
    plt.tight_layout()
    plt.savefig(out, dpi=160)
    plt.close(fig)
    print(f'Wrote {out}')


if __name__ == '__main__':
    main()
