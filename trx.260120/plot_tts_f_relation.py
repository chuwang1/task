#!/usr/bin/env python3

import os
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/Users/dengxiaoya/TASK/CFEDR/git/task/trx.260120"
GFILE = "/Users/dengxiaoya/CFEDRSW/OMFIT_out/g260206.20000_teq_0114"
COMPARE_F_CSV = "/Users/dengxiaoya/TASK/latest/task/eq/in/compare_f_from_gfile_geom_q.csv"


def read_tr_eqcalq(path):
    lines = open(path).read().splitlines()
    rho = []
    psip = []
    tts = []
    for line in lines[1:]:
        parts = line.split(',')
        vals = [float(x) for x in parts[1].split()]
        rho.append(vals[0])
        psip.append(vals[1])
        tts.append(vals[6])
    return np.array(rho), np.array(psip), np.array(tts)


def read_ttps_chain(path):
    lines = open(path).read().splitlines()
    psips = []
    ttps = []
    ttfunc = []
    for line in lines[1:]:
        parts = line.split(',')
        vals = [float(x) for x in parts[1].split()]
        psips.append(vals[0])
        ttps.append(vals[1])
        ttfunc.append(vals[2])
    return np.array(psips), np.array(ttps), np.array(ttfunc)


def read_gfile_fpol(path):
    with open(path, 'r') as f:
        lines = f.readlines()
    header = lines[0].strip().split()
    nw = int(header[-2])

    def parse_floats(line):
        res = []
        for i in range(0, len(line), 16):
            chunk = line[i:i+16].strip()
            if not chunk:
                continue
            if '-' in chunk[1:] and 'E' not in chunk and 'e' not in chunk:
                idx = chunk.rfind('-'); chunk = chunk[:idx] + 'E' + chunk[idx:]
            elif '+' in chunk[1:] and 'E' not in chunk and 'e' not in chunk:
                idx = chunk.rfind('+'); chunk = chunk[:idx] + 'E' + chunk[idx:]
            try:
                res.append(float(chunk))
            except ValueError:
                pass
        return res

    data = []
    for line in lines[5:]:
        data.extend(parse_floats(line.rstrip('\n')))
    fpol = np.array(data[:nw])
    psi_n = np.linspace(0.0, 1.0, nw)
    return psi_n, fpol


def read_f_rebuilt(path):
    psi_n = []
    f_rebuilt = []
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            psi_n.append(float(row['psi_n']))
            f_rebuilt.append(float(row['F_rebuilt']))
    return np.array(psi_n), np.array(f_rebuilt)


def main():
    rho_tr, psip_tr, tts_tr = read_tr_eqcalq(os.path.join(BASE, 'eqcalq_dvdpsit_raw.csv'))
    psips_chain, ttps_raw, ttfunc_on_psips = read_ttps_chain(os.path.join(BASE, 'eqcalq_ttps_chain.csv'))
    psi_n, fpol = read_gfile_fpol(GFILE)
    psi_n_reb, f_rebuilt = read_f_rebuilt(COMPARE_F_CSV)

    # psi-grid comparison
    psip_norm = psips_chain / psips_chain[-1]
    twopi_fpol_on_psi = -2.0 * np.pi * np.interp(psip_norm, psi_n, fpol)
    f_rebuilt_on_psi = np.interp(psip_norm, psi_n_reb, f_rebuilt)
    err_same = np.nanmean((ttps_raw - f_rebuilt_on_psi) ** 2)
    err_flip = np.nanmean((ttps_raw + f_rebuilt_on_psi) ** 2)
    if err_flip < err_same:
        f_rebuilt = -f_rebuilt
        f_rebuilt_on_psi = -f_rebuilt_on_psi

    # rho-grid comparison via PSIP(NR)
    twopi_fpol_on_rho = np.interp(psip_tr, psips_chain, twopi_fpol_on_psi)
    f_rebuilt_on_rho = np.interp(psip_tr, psips_chain, f_rebuilt_on_psi)
    ttfunc_on_rho = np.interp(psip_tr, psips_chain, ttfunc_on_psips)

    fig, axes = plt.subplots(2, 2, figsize=(10, 8), sharex=True)

    axes[0, 0].plot(psip_norm, ttps_raw, 'o-', ms=3, lw=1.8, label='TTPS raw on psi-grid')
    axes[0, 0].plot(psip_norm, ttfunc_on_psips, '--', lw=1.8, color='tab:orange', label='TTFUNC(PSIPS)')
    axes[0, 0].plot(psip_norm, twopi_fpol_on_psi, ':', lw=2.0, color='tab:red', label='-2π * gfile fpol')
    axes[0, 0].plot(psip_norm, f_rebuilt_on_psi, '-.', lw=1.8, color='tab:green', label='F_rebuilt aligned')
    axes[0, 0].set_ylabel('2π R Bt')
    axes[0, 0].set_title('psi-grid comparison')
    axes[0, 0].grid(True, alpha=0.3)
    axes[0, 0].legend(fontsize=8)

    axes[0, 1].plot(rho_tr, tts_tr, 'o-', ms=3, lw=1.8, label='TTS on rho-grid')
    axes[0, 1].plot(rho_tr, ttfunc_on_rho, '--', lw=1.8, color='tab:orange', label='TTFUNC at PSIP(NR)')
    axes[0, 1].plot(rho_tr, twopi_fpol_on_rho, ':', lw=2.0, color='tab:red', label='interp(-2π fpol -> PSIP(NR))')
    axes[0, 1].plot(rho_tr, f_rebuilt_on_rho, '-.', lw=1.8, color='tab:green', label='interp(F_rebuilt -> PSIP(NR))')
    axes[0, 1].set_ylabel('2π R Bt')
    axes[0, 1].set_title('rho-grid comparison')
    axes[0, 1].grid(True, alpha=0.3)
    axes[0, 1].legend(fontsize=8)

    ratio1 = ttfunc_on_psips / np.where(np.abs(ttps_raw) > 1e-30, ttps_raw, np.nan)
    ratio1g = ttps_raw / np.where(np.abs(twopi_fpol_on_psi) > 1e-30, twopi_fpol_on_psi, np.nan)
    ratio1r = ttps_raw / np.where(np.abs(f_rebuilt_on_psi) > 1e-30, f_rebuilt_on_psi, np.nan)
    axes[1, 0].plot(psip_norm, ratio1, 'd-', color='tab:blue', lw=1.5, ms=3, label='TTFUNC(PSIPS) / TTPS raw')
    axes[1, 0].plot(psip_norm, ratio1g, '--', color='tab:red', lw=1.5, label='TTPS raw / (-2π fpol)')
    axes[1, 0].plot(psip_norm, ratio1r, '-.', color='tab:green', lw=1.5, label='TTPS raw / F_rebuilt')
    axes[1, 0].axhline(1.0, color='k', lw=1, alpha=0.5)
    axes[1, 0].set_xlabel('psi_p normalized')
    axes[1, 0].set_ylabel('ratio')
    axes[1, 0].set_title('psi-grid ratios')
    axes[1, 0].grid(True, alpha=0.3)
    axes[1, 0].legend(fontsize=8)

    ratio2 = tts_tr / np.where(np.abs(ttfunc_on_rho) > 1e-30, ttfunc_on_rho, np.nan)
    ratio2g = tts_tr / np.where(np.abs(twopi_fpol_on_rho) > 1e-30, twopi_fpol_on_rho, np.nan)
    ratio2r = tts_tr / np.where(np.abs(f_rebuilt_on_rho) > 1e-30, f_rebuilt_on_rho, np.nan)
    axes[1, 1].plot(rho_tr, ratio2, 'd-', color='tab:blue', lw=1.5, ms=3, label='TTS / TTFUNC(PSIP(NR))')
    axes[1, 1].plot(rho_tr, ratio2g, '--', color='tab:red', lw=1.5, label='TTS / interp(-2π fpol)')
    axes[1, 1].plot(rho_tr, ratio2r, '-.', color='tab:green', lw=1.5, label='TTS / interp(F_rebuilt)')
    axes[1, 1].axhline(1.0, color='k', lw=1, alpha=0.5)
    axes[1, 1].set_xlabel('rho')
    axes[1, 1].set_ylabel('ratio')
    axes[1, 1].set_title('rho-grid ratios')
    axes[1, 1].grid(True, alpha=0.3)
    axes[1, 1].legend(fontsize=8)

    out = os.path.join(BASE, 'tts_fpol_relation.png')
    plt.tight_layout()
    plt.savefig(out, dpi=160)
    plt.close(fig)
    print(f'Wrote {out}')


if __name__ == '__main__':
    main()
