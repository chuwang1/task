#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


BASE = "/Users/dengxiaoya/TASK/CFEDR/git/task/trx.260120"
OMFIT_OUT = "/Users/dengxiaoya/CFEDRSW/OMFIT_out"


def read_tr_eqcalq(path):
    lines = open(path).read().splitlines()
    out = {"nr": [], "rhot": [], "psip": [], "psit": [], "dvdpsit": [], "dvdpsip": [], "qps": []}
    for line in lines[1:]:
        parts = line.split(',')
        vals = [float(x) for x in parts[1].split()]
        out["nr"].append(int(parts[0]))
        out["rhot"].append(vals[0])
        out["psip"].append(vals[1])
        out["psit"].append(vals[2])
        out["dvdpsit"].append(vals[3])
        out["dvdpsip"].append(vals[4])
        out["qps"].append(vals[5])
    return {k: np.array(v) for k, v in out.items()}


def integrate_psit(psip, q):
    psit = np.zeros_like(psip)
    for i in range(1, len(psip)):
        qeff = 2.0 * q[i] * q[i - 1] / (q[i] + q[i - 1])
        psit[i] = psit[i - 1] + qeff * (psip[i] - psip[i - 1])
    return psit


def load_omfit():
    allp = pd.read_csv(os.path.join(OMFIT_OUT, 'outputs_new', 'all_profiles.csv'))
    rho = []
    polflux = []
    with open(os.path.join(OMFIT_OUT, 'input.profiles')) as f:
        for line in f:
            if line.startswith('#') or 'N_ION' in line or '=' in line or not line.strip():
                continue
            parts = line.split()
            if len(parts) == 5:
                try:
                    rho.append(float(parts[0]))
                    polflux.append(float(parts[2]))
                except ValueError:
                    pass
            if len(rho) == 201:
                break
    rho = np.array(rho)
    psip = np.abs(np.array(polflux)) * 2.0 * np.pi
    q = np.abs(np.interp(rho, allp['rho'].values, allp['q'].values))
    psit = integrate_psit(psip, q)
    dvdpsip = np.abs(np.gradient(pd.read_csv(os.path.join(OMFIT_OUT, 'outputs_new', 'input_profiles_extra.csv'))['EXPRO_vol'].values, psip))
    dvdpsit = dvdpsip / q
    return {"rho": rho, "psip": psip, "psit": psit, "qps": q, "dvdpsip": dvdpsip, "dvdpsit": dvdpsit}


def main():
    tr = read_tr_eqcalq(os.path.join(BASE, 'eqcalq_dvdpsit_raw.csv'))
    om = load_omfit()

    x_tr = tr['psip'] / tr['psip'][-1]
    x_om = om['psip'] / om['psip'][-1]
    om_interp = {name: np.interp(x_tr, x_om, om[name]) for name in ['psip', 'psit', 'qps', 'dvdpsip', 'dvdpsit']}

    fig, axes = plt.subplots(5, 2, figsize=(13, 16), sharex='col')
    items = [
        ('psip', 'PSIP'),
        ('psit', 'PSIT'),
        ('qps', 'QPS / q'),
        ('dvdpsip', 'dV/dpsi_p'),
        ('dvdpsit', 'dV/dpsi_t'),
    ]

    for i, (key, title) in enumerate(items):
        ax_l = axes[i, 0]
        ax_r = axes[i, 1]
        ax_l.plot(x_tr, tr[key], 'o-', ms=3, lw=1.8, label=f'TR {title}')
        ax_l.plot(x_tr, om_interp[key], '--', lw=2.0, color='tab:red', label=f'OMFIT {title}')
        ax_l.set_ylabel(title)
        ax_l.set_title(f'{title} on normalized psi_p grid')
        ax_l.grid(True, alpha=0.3)
        ax_l.legend(fontsize=8)

        ratio = tr[key] / np.maximum(om_interp[key], 1e-30)
        ax_r.plot(x_tr, ratio, 'd-', color='tab:blue', lw=1.5, ms=3)
        ax_r.axhline(1.0, color='k', lw=1, alpha=0.5)
        ax_r.set_ylabel('TR / OMFIT')
        ax_r.set_title(f'{title} ratio on normalized psi_p grid')
        ax_r.grid(True, alpha=0.3)

    axes[-1, 0].set_xlabel('psi_p / psi_pa')
    axes[-1, 1].set_xlabel('psi_p / psi_pa')

    out = os.path.join(BASE, 'equilibrium_overview_comparison.png')
    plt.tight_layout()
    plt.savefig(out, dpi=160)
    plt.close(fig)
    print(f'Wrote {out}')


if __name__ == '__main__':
    main()
