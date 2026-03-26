#!/usr/bin/env python3

import os
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/Users/dengxiaoya/TASK/latest/task/trx"
GFILE = "/Users/dengxiaoya/CFEDRSW/OMFIT_out/g260206.20000_teq_0114"
COMPARE_F_CSV = "/Users/dengxiaoya/TASK/latest/task/eq/in/compare_f_from_gfile_geom_q.csv"


def read_ttps_chain(path):
    rows = open(path).read().splitlines()
    out = {"psips": [], "ttps_raw": [], "ttfunc_on_psips": []}
    for line in rows[1:]:
        parts = line.split(',')
        vals = [float(x) for x in parts[1].split()]
        out["psips"].append(vals[0])
        out["ttps_raw"].append(vals[1])
        out["ttfunc_on_psips"].append(vals[2])
    return {k: np.array(v) for k, v in out.items()}


def read_tts_rhot(path):
    rows = open(path).read().splitlines()
    out = {"rho": [], "psip": [], "tts": []}
    for line in rows[1:]:
        parts = line.split(',')
        vals = [float(x) for x in parts[1].split()]
        out["rho"].append(vals[0])
        out["psip"].append(vals[1])
        out["tts"].append(vals[6])
    return {k: np.array(v) for k, v in out.items()}


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
    psi_n = []; f = []
    with open(path) as fh:
        rdr = csv.DictReader(fh)
        for row in rdr:
            psi_n.append(float(row['psi_n']))
            f.append(float(row['F_rebuilt']))
    return np.array(psi_n), np.array(f)


def main():
    chain = read_ttps_chain(os.path.join(BASE, 'eqcalq_ttps_chain.csv'))
    ttsrho = read_tts_rhot(os.path.join(BASE, 'eqcalq_dvdpsit_raw.csv'))
    psi_n, fpol = read_gfile_fpol(GFILE)
    psi_freb, freb = read_f_rebuilt(COMPARE_F_CSV)

    # sign align rebuilt F to TTPS raw
    ttps_ref = np.interp(psi_freb, chain['psips'] / chain['psips'][-1], chain['ttps_raw'])
    err_same = np.nanmean((ttps_ref - freb) ** 2)
    err_flip = np.nanmean((ttps_ref + freb) ** 2)
    if err_flip < err_same:
        freb = -freb

    fig, axes = plt.subplots(2, 2, figsize=(11, 8))

    # Panel 1: raw TTPS chain on psi_p grid
    psi_pn = chain['psips'] / chain['psips'][-1]
    axes[0, 0].plot(psi_pn, chain['ttps_raw'], lw=2, label='TTPS raw')
    axes[0, 0].plot(psi_pn, chain['ttfunc_on_psips'], '--', lw=1.8, label='TTFUNC(PSIPS)')
    axes[0, 0].plot(psi_n, -2*np.pi*fpol, ':', lw=2, color='tab:red', label='-2π fpol')
    axes[0, 0].plot(psi_freb, freb, '-.', lw=1.8, color='tab:green', label='F_rebuilt aligned')
    axes[0, 0].set_title('Stage 1: raw F data on psi_p grid')
    axes[0, 0].set_ylabel('F-like quantity')
    axes[0, 0].grid(True, alpha=0.3)
    axes[0, 0].legend(fontsize=8)

    # Panel 2: TTPS raw vs TTS on RHOT using PSIP as bridge
    tts_from_ttps = np.interp(ttsrho['psip'], chain['psips'], chain['ttps_raw'])
    axes[0, 1].plot(ttsrho['rho'], ttsrho['tts'], lw=2, label='TTS on RHOT')
    axes[0, 1].plot(ttsrho['rho'], tts_from_ttps, '--', lw=1.8, color='tab:orange', label='interp(TTPS raw -> PSIP(NR))')
    axes[0, 1].set_title('Stage 2: does TTS equal TTPS evaluated at PSIP?')
    axes[0, 1].grid(True, alpha=0.3)
    axes[0, 1].legend(fontsize=8)

    # Panel 3: ratios on psi_p grid
    ratio_tt = chain['ttfunc_on_psips'] / np.maximum(chain['ttps_raw'], 1e-30)
    ratio_fp = chain['ttps_raw'] / np.maximum(np.interp(psi_pn, psi_n, -2*np.pi*fpol), 1e-30)
    axes[1, 0].plot(psi_pn, ratio_tt, lw=1.8, label='TTFUNC(PSIPS)/TTPS raw')
    axes[1, 0].plot(psi_pn, ratio_fp, '--', lw=1.8, label='TTPS raw / (-2π fpol)')
    axes[1, 0].axhline(1.0, color='k', lw=1, alpha=0.5)
    axes[1, 0].set_xlabel('psi_p normalized')
    axes[1, 0].set_ylabel('ratio')
    axes[1, 0].set_title('Consistency on psi_p grid')
    axes[1, 0].grid(True, alpha=0.3)
    axes[1, 0].legend(fontsize=8)

    # Panel 4: ratios on rho grid
    ratio_rho = ttsrho['tts'] / np.maximum(tts_from_ttps, 1e-30)
    axes[1, 1].plot(ttsrho['rho'], ratio_rho, lw=1.8, label='TTS / interp(TTPS raw)')
    axes[1, 1].axhline(1.0, color='k', lw=1, alpha=0.5)
    axes[1, 1].set_xlabel('rho')
    axes[1, 1].set_ylabel('ratio')
    axes[1, 1].set_title('Consistency on rho grid')
    axes[1, 1].grid(True, alpha=0.3)
    axes[1, 1].legend(fontsize=8)

    out = os.path.join(BASE, 'tts_chain_detailed.png')
    plt.tight_layout()
    plt.savefig(out, dpi=160)
    plt.close(fig)
    print(f'Wrote {out}')


if __name__ == '__main__':
    main()
