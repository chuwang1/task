#!/usr/bin/env python3

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/Users/dengxiaoya/TASK/latest/task/trx"
GFILE = "/Users/dengxiaoya/CFEDRSW/OMFIT_out/g260206.20000_teq_0114"


def read_tr_eqcalq(path):
    lines = open(path).read().splitlines()
    out = {"rho": [], "psip": [], "qps": [], "tts": [], "sumavir2": []}
    for line in lines[1:]:
        parts = line.split(',')
        vals = [float(x) for x in parts[1].split()]
        out["rho"].append(vals[0])
        out["psip"].append(vals[1])
        out["qps"].append(vals[5])
        out["tts"].append(vals[6])
        out["sumavir2"].append(vals[7])
    return {k: np.array(v) for k, v in out.items()}


def read_gfile_fpol_q(path):
    with open(path, 'r') as f:
        lines = f.readlines()
    header = lines[0].strip().split()
    nw = int(header[-2]); nh = int(header[-1])

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

    l3 = parse_floats(lines[2])
    simag, sibry = l3[2], l3[3]

    data = []
    for line in lines[5:]:
        data.extend(parse_floats(line.rstrip('\n')))

    off = 0
    fpol = np.array(data[off:off+nw]); off += nw
    pres = np.array(data[off:off+nw]); off += nw
    ffprime = np.array(data[off:off+nw]); off += nw
    pprime = np.array(data[off:off+nw]); off += nw
    off += nw * nh
    qpsi = np.array(data[off:off+nw]); off += nw

    psi_n = np.linspace(0.0, 1.0, nw)
    return psi_n, np.abs(fpol), np.abs(qpsi)


def main():
    tr = read_tr_eqcalq(os.path.join(BASE, 'eqcalq_dvdpsit_raw.csv'))
    psi_n_gf, fpol_gf, q_gf = read_gfile_fpol_q(GFILE)

    x_tr = tr['psip'] / tr['psip'][-1]
    # TTS in TR corresponds approximately to 2*pi*Fpol
    tts_gf = np.interp(x_tr, psi_n_gf, 2.0 * np.pi * fpol_gf)
    q_gf_on_tr = np.interp(x_tr, psi_n_gf, q_gf)
    sumavir2_gf = 4.0 * np.pi**2 * q_gf_on_tr / np.maximum(tts_gf, 1e-30)

    fig, axes = plt.subplots(3, 2, figsize=(12, 11), sharex='col')
    items = [
        ('qps', q_gf_on_tr, 'QPS / q'),
        ('tts', tts_gf, 'TTS ~ 2π Fpol'),
        ('sumavir2', sumavir2_gf, 'SUMAVIR2 = 4π² q / TTS'),
    ]
    for i, (k, ref, title) in enumerate(items):
        axes[i, 0].plot(x_tr, tr[k], 'o-', ms=3, lw=1.8, label=f'TR {title}')
        axes[i, 0].plot(x_tr, ref, '--', lw=2.0, color='tab:red', label=f'gfile/OMFIT {title}')
        axes[i, 0].set_ylabel(title)
        axes[i, 0].grid(True, alpha=0.3)
        axes[i, 0].legend(fontsize=8)

        ratio = tr[k] / np.maximum(ref, 1e-30)
        axes[i, 1].plot(x_tr, ratio, 'd-', color='tab:blue', lw=1.5, ms=3)
        axes[i, 1].axhline(1.0, color='k', lw=1, alpha=0.5)
        axes[i, 1].set_ylabel('TR / ref')
        axes[i, 1].grid(True, alpha=0.3)
        axes[i, 1].set_title(f'{title} ratio')

    axes[-1, 0].set_xlabel('psi_p / psi_pa')
    axes[-1, 1].set_xlabel('psi_p / psi_pa')
    axes[0, 0].set_title('QPS and its component comparison on normalized psi_p grid')

    out = os.path.join(BASE, 'qps_components_comparison.png')
    plt.tight_layout()
    plt.savefig(out, dpi=160)
    plt.close(fig)
    print(f'Wrote {out}')


if __name__ == '__main__':
    main()
