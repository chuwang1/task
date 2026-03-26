#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


BASE = "/Users/dengxiaoya/TASK/CFEDR/git/task/trx.260120"


def read_tr_eqcalq(path):
    lines = open(path).read().splitlines()
    out = {"nr": [], "rhot": [], "psip": [], "psit": [], "qps": []}
    for line in lines[1:]:
        parts = line.split(',')
        vals = [float(x) for x in parts[1].split()]
        out["nr"].append(int(parts[0]))
        out["rhot"].append(vals[0])
        out["psip"].append(vals[1])
        out["psit"].append(vals[2])
        out["qps"].append(vals[5])
    return {k: np.array(v) for k, v in out.items()}


def integrate_psit(psip, q):
    psit = np.zeros_like(psip)
    for i in range(1, len(psip)):
        qeff = 2.0 * q[i] * q[i - 1] / (q[i] + q[i - 1])
        psit[i] = psit[i - 1] + qeff * (psip[i] - psip[i - 1])
    return psit


def read_gfile_psipa_q(path):
    with open(path, 'r') as f:
        lines = f.readlines()
    hdr = lines[0].split()
    nw = int(hdr[-2]); nh = int(hdr[-1])
    def parse(line):
        res = []
        for i in range(0, len(line), 16):
            c = line[i:i+16].strip()
            if not c:
                continue
            if '-' in c[1:] and 'E' not in c and 'e' not in c:
                idx = c.rfind('-'); c = c[:idx] + 'E' + c[idx:]
            elif '+' in c[1:] and 'E' not in c and 'e' not in c:
                idx = c.rfind('+'); c = c[:idx] + 'E' + c[idx:]
            try:
                res.append(float(c))
            except ValueError:
                pass
        return res
    l3 = parse(lines[2])
    simag, sibry = l3[2], l3[3]
    data = []
    for line in lines[5:]:
        data.extend(parse(line.rstrip('\n')))
    qpsi = np.abs(np.array(data[4*nw + nw*nh: 4*nw + nw*nh + nw]))
    psi_n = np.linspace(0.0, 1.0, nw)
    psipa = 2.0 * np.pi * (sibry - simag)
    return psi_n, psipa, qpsi


def main():
    tr = read_tr_eqcalq(os.path.join(BASE, "eqcalq_dvdpsit_raw.csv"))

    # gfile-based reference for absolute PSIP scale and q
    psi_n_gf, psipa_gf, qpsi_gf = read_gfile_psipa_q('/Users/dengxiaoya/CFEDRSW/OMFIT_out/g260206.20000_teq_0114')
    psi_n_om = psi_n_gf
    psip_om = psi_n_om * psipa_gf
    q_om = qpsi_gf
    psit_om = integrate_psit(psip_om, q_om)

    x_tr = tr['psip'] / tr['psip'][-1]
    psip_om_on_tr = np.interp(x_tr, psi_n_om, psip_om)
    q_om_on_tr = np.interp(x_tr, psi_n_om, q_om)
    psit_om_on_tr = np.interp(x_tr, psi_n_om, psit_om)

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True)

    axes[0, 0].plot(x_tr, tr['psip'], 'o-', ms=3, lw=2, label='TR PSIP')
    axes[0, 0].plot(x_tr, psip_om_on_tr, '--', lw=2, color='tab:red', label='gfile PSIP')
    axes[0, 0].set_ylabel('PSIP')
    axes[0, 0].set_title('Poloidal Flux')
    axes[0, 0].grid(True, alpha=0.3)
    axes[0, 0].legend(fontsize=8)

    axes[0, 1].plot(x_tr, tr['psit'], 'o-', ms=3, lw=2, label='TR PSIT')
    axes[0, 1].plot(x_tr, psit_om_on_tr, '--', lw=2, color='tab:red', label='gfile PSIT (integrated)')
    axes[0, 1].set_ylabel('PSIT')
    axes[0, 1].set_title('Toroidal Flux')
    axes[0, 1].grid(True, alpha=0.3)
    axes[0, 1].legend(fontsize=8)

    r_psip = tr['psip'] / np.maximum(psip_om_on_tr, 1e-30)
    r_psit = tr['psit'] / np.maximum(psit_om_on_tr, 1e-30)
    axes[1, 0].plot(x_tr, r_psip, 'd-', color='tab:blue', lw=1.5, ms=3)
    axes[1, 0].axhline(1.0, color='k', lw=1, alpha=0.5)
    axes[1, 0].set_xlabel('rho')
    axes[1, 0].set_ylabel('TR / OMFIT')
    axes[1, 0].set_title('PSIP ratio')
    axes[1, 0].grid(True, alpha=0.3)

    axes[1, 1].plot(x_tr, r_psit, 'd-', color='tab:blue', lw=1.5, ms=3)
    axes[1, 1].axhline(1.0, color='k', lw=1, alpha=0.5)
    axes[1, 1].set_xlabel('rho')
    axes[1, 1].set_ylabel('TR / OMFIT')
    axes[1, 1].set_title('PSIT ratio')
    axes[1, 1].grid(True, alpha=0.3)

    out = os.path.join(BASE, 'psit_psip_comparison.png')
    plt.tight_layout()
    plt.savefig(out, dpi=160)
    plt.close(fig)
    print(f'Wrote {out}')
    print('PSIP mean rel err =', float(np.mean(np.abs(r_psip - 1.0))))
    print('PSIT mean rel err =', float(np.mean(np.abs(r_psit - 1.0))))


if __name__ == '__main__':
    main()
