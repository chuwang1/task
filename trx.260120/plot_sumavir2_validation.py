#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/Users/dengxiaoya/TASK/CFEDR/git/task/trx.260120"


def read_eqcalq(path):
    lines = open(path).read().splitlines()
    rho = []
    psip = []
    qps = []
    tts = []
    suma = []
    for line in lines[1:]:
        parts = line.split(',')
        vals = [float(x) for x in parts[1].split()]
        rho.append(vals[0])
        psip.append(vals[1])
        qps.append(vals[5])
        tts.append(vals[6])
        suma.append(vals[7])
    return np.array(rho), np.array(psip), np.array(qps), np.array(tts), np.array(suma)


def main():
    rho, psip, qps, tts, suma = read_eqcalq(os.path.join(BASE, 'eqcalq_dvdpsit_raw.csv'))
    # compare q on the proper normalized psi_p grid, not rho
    x_tr = psip / psip[-1]
    gpath = '/Users/dengxiaoya/CFEDRSW/OMFIT_out/g260206.20000_teq_0114'
    lines = open(gpath).read().splitlines()
    hdr = lines[0].split(); nw = int(hdr[-2]); nh = int(hdr[-1])
    def parse(line):
        res=[]
        for i in range(0,len(line),16):
            c=line[i:i+16].strip()
            if not c: continue
            if '-' in c[1:] and 'E' not in c and 'e' not in c:
                idx=c.rfind('-'); c=c[:idx]+'E'+c[idx:]
            elif '+' in c[1:] and 'E' not in c and 'e' not in c:
                idx=c.rfind('+'); c=c[:idx]+'E'+c[idx:]
            try: res.append(float(c))
            except: pass
        return res
    data=[]
    for line in lines[5:]: data.extend(parse(line.rstrip('\n')))
    qpsi = np.abs(np.array(data[4*nw + nw*nh : 4*nw + nw*nh + nw]))
    psi_n = np.linspace(0.0,1.0,nw)
    q_g = np.interp(x_tr, psi_n, qpsi)

    q_from_tr = suma * tts / (4.0 * np.pi**2)
    suma_from_g = 4.0 * np.pi**2 * q_g / tts

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex='col')

    axes[0, 0].plot(x_tr, qps, 'o-', ms=3, lw=1.8, label='TR QPS')
    axes[0, 0].plot(x_tr, q_from_tr, '--', lw=2.0, color='tab:green', label='SUMAVIR2*TTS/4π²')
    axes[0, 0].plot(x_tr, q_g, ':', lw=2.0, color='tab:red', label='gfile qpsi')
    axes[0, 0].set_ylabel('q')
    axes[0, 0].set_title('QPS validation')
    axes[0, 0].grid(True, alpha=0.3)
    axes[0, 0].legend(fontsize=8)

    axes[0, 1].plot(x_tr, suma, 'o-', ms=3, lw=1.8, label='TR SUMAVIR2')
    axes[0, 1].plot(x_tr, suma_from_g, '--', lw=2.0, color='tab:red', label='4π² q_gfile / TTS')
    axes[0, 1].set_ylabel('SUMAVIR2')
    axes[0, 1].set_title('SUMAVIR2 comparison')
    axes[0, 1].grid(True, alpha=0.3)
    axes[0, 1].legend(fontsize=8)

    ratio_q = qps / np.maximum(q_g, 1e-30)
    ratio_s = suma / np.maximum(suma_from_g, 1e-30)
    axes[1, 0].plot(x_tr, ratio_q, 'd-', color='tab:blue', lw=1.5, ms=3, label='QPS / q_gfile')
    axes[1, 0].plot(x_tr, ratio_s, '--', color='tab:green', lw=1.5, label='SUMAVIR2 / implied_SUMAVIR2')
    axes[1, 0].axhline(1.0, color='k', lw=1, alpha=0.5)
    axes[1, 0].set_xlabel('psi_p / psi_pa')
    axes[1, 0].set_ylabel('ratio')
    axes[1, 0].set_title('Mismatch equivalence')
    axes[1, 0].grid(True, alpha=0.3)
    axes[1, 0].legend(fontsize=8)

    axes[1, 1].plot(x_tr, tts, lw=1.8, color='tab:purple', label='TTS')
    axes[1, 1].set_xlabel('psi_p / psi_pa')
    axes[1, 1].set_ylabel('TTS')
    axes[1, 1].set_title('TTS stays smooth / nearly fixed')
    axes[1, 1].grid(True, alpha=0.3)
    axes[1, 1].legend(fontsize=8)

    out = os.path.join(BASE, 'sumavir2_validation.png')
    plt.tight_layout()
    plt.savefig(out, dpi=160)
    plt.close(fig)
    print(f'Wrote {out}')
    print('max |QPS - SUMAVIR2*TTS/4π²| =', float(np.max(np.abs(qps - q_from_tr))))


if __name__ == '__main__':
    main()
