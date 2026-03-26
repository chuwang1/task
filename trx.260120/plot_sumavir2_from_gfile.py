#!/usr/bin/env python3
"""Compute SUMAVIR2 from gfile contours and compare with TR eqcalq output.

SUMAVIR2 in TR is the contour integral

    SUMAVIR2 = ∮ dl / (Bp * R^2)

with Bp = |grad psi_p| / (2*pi*R).

This script reconstructs the same quantity directly from gfile contours.
"""

from __future__ import annotations

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path


BASE = "/Users/dengxiaoya/TASK/CFEDR/git/task/trx.260120"
OMFIT_BASE = "/Users/dengxiaoya/CFEDRSW/OMFIT_out"
GFILE = os.path.join(OMFIT_BASE, "g260206.20000_teq_0114")
INPUT_PROFILES = os.path.join(OMFIT_BASE, "input.profiles")


def read_gfile_2d(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()

    header = lines[0].strip().split()
    nw = int(header[-2])
    nh = int(header[-1])

    def parse_floats(line):
        res = []
        for i in range(0, len(line), 16):
            chunk = line[i:i + 16].strip()
            if not chunk:
                continue
            if '-' in chunk[1:] and 'E' not in chunk and 'e' not in chunk:
                idx = chunk.rfind('-')
                chunk = chunk[:idx] + 'E' + chunk[idx:]
            elif '+' in chunk[1:] and 'E' not in chunk and 'e' not in chunk:
                idx = chunk.rfind('+')
                chunk = chunk[:idx] + 'E' + chunk[idx:]
            try:
                res.append(float(chunk))
            except ValueError:
                pass
        return res

    l2 = parse_floats(lines[1])
    rdim, zdim, rcent, rleft, zmid = l2[0:5]

    l3 = parse_floats(lines[2])
    rmaxis, zmaxis, simag, sibry, bcent = l3[0:5]

    data = []
    for line in lines[5:]:
        data.extend(parse_floats(line.rstrip('\n')))

    cur = 4 * nw
    psirz_1d = data[cur: cur + nw * nh]
    psirz = np.array(psirz_1d).reshape((nh, nw))

    R = np.linspace(rleft, rleft + rdim, nw)
    Z = np.linspace(zmid - zdim / 2.0, zmid + zdim / 2.0, nh)
    psi_n = (psirz - simag) / (sibry - simag)
    # convert to TR PSIP convention: psi_p in Wb, not Wb/rad
    psi_p = 2.0 * np.pi * (psirz - simag)
    return R, Z, psi_n, psi_p, rmaxis, zmaxis


def read_input_profiles_mapping(filename):
    rho_arr = []
    polflux_arr = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#') or 'N_ION' in line or '=' in line or not line.strip():
                continue
            parts = line.split()
            if len(parts) == 5:
                try:
                    rho = float(parts[0])
                    polflux = float(parts[2])
                    if 0.0 <= rho <= 1.05:
                        rho_arr.append(rho)
                        polflux_arr.append(polflux)
                except ValueError:
                    pass
            if len(rho_arr) == 201:
                break
    rho_arr = np.array(rho_arr)
    polflux_arr = np.array(polflux_arr)
    psi_axis = polflux_arr[0]
    psi_bry = polflux_arr[-1]
    psi_n_arr = (polflux_arr - psi_axis) / (psi_bry - psi_axis)
    return psi_n_arr, rho_arr


def bilinear_interp(xarr, yarr, f2d, x, y):
    if x <= xarr[0]:
        i = 0
    elif x >= xarr[-1]:
        i = len(xarr) - 2
    else:
        i = np.searchsorted(xarr, x) - 1
    if y <= yarr[0]:
        j = 0
    elif y >= yarr[-1]:
        j = len(yarr) - 2
    else:
        j = np.searchsorted(yarr, y) - 1
    x1, x2 = xarr[i], xarr[i + 1]
    y1, y2 = yarr[j], yarr[j + 1]
    tx = 0.0 if x2 == x1 else (x - x1) / (x2 - x1)
    ty = 0.0 if y2 == y1 else (y - y1) / (y2 - y1)
    f11 = f2d[j, i]
    f21 = f2d[j, i + 1]
    f12 = f2d[j + 1, i]
    f22 = f2d[j + 1, i + 1]
    return (1 - tx) * (1 - ty) * f11 + tx * (1 - ty) * f21 + (1 - tx) * ty * f12 + tx * ty * f22


def read_tr_eqcalq(path):
    lines = open(path).read().splitlines()
    out = {"rho": [], "psip": [], "sumavir2": [], "qps": [], "tts": []}
    for line in lines[1:]:
        parts = line.split(',')
        vals = [float(x) for x in parts[1].split()]
        out["rho"].append(vals[0])
        out["psip"].append(vals[1])
        out["qps"].append(vals[5])
        out["tts"].append(vals[6])
        out["sumavir2"].append(vals[7])
    return {k: np.array(v) for k, v in out.items()}


def main():
    R, Z, psi_n_2d, psi_p_2d, rmaxis, zmaxis = read_gfile_2d(GFILE)
    psi_n_map, rho_map = read_input_profiles_mapping(INPUT_PROFILES)
    tr = read_tr_eqcalq(os.path.join(BASE, 'eqcalq_dvdpsit_raw.csv'))

    # gradients of physical psi_p [Wb]
    dpsi_dZ, dpsi_dR = np.gradient(psi_p_2d, Z, R)

    fig_dummy, ax = plt.subplots()
    sumavir2_g = np.zeros_like(rho_map)
    valid = np.ones_like(rho_map, dtype=bool)

    for i, p_target in enumerate(psi_n_map):
        if p_target <= 1e-6:
            sumavir2_g[i] = sumavir2_g[1] if i == 0 and len(rho_map) > 1 else 0.0
            continue
        cs = ax.contour(R, Z, psi_n_2d, levels=[p_target])
        paths = []
        if len(cs.allsegs) > 0 and len(cs.allsegs[0]) > 0:
            for seg in cs.allsegs[0]:
                paths.append(Path(seg))
        if not paths:
            sumavir2_g[i] = np.nan
            valid[i] = False
            continue
        valid_paths = []
        for p in paths:
            v = p.vertices
            if np.min(v[:, 0]) < rmaxis < np.max(v[:, 0]) and np.min(v[:, 1]) < zmaxis < np.max(v[:, 1]):
                valid_paths.append(p)
        if not valid_paths:
            valid_paths = paths
        path = max(valid_paths, key=lambda x: len(x.vertices))
        v = path.vertices
        if not np.allclose(v[0], v[-1]):
            v = np.vstack([v, v[0]])

        s = 0.0
        for k in range(1, len(v)):
            x1, y1 = v[k - 1]
            x2, y2 = v[k]
            dl = np.hypot(x2 - x1, y2 - y1)
            rm = 0.5 * (x1 + x2)
            zm = 0.5 * (y1 + y2)
            dpr = bilinear_interp(R, Z, dpsi_dR, rm, zm)
            dpz = bilinear_interp(R, Z, dpsi_dZ, rm, zm)
            bp = np.sqrt(dpr * dpr + dpz * dpz) / (2.0 * np.pi * rm)
            if bp <= 1e-30 or rm <= 1e-30:
                continue
            s += dl / (bp * rm * rm)
        sumavir2_g[i] = s

    plt.close(fig_dummy)

    # interpolate gfile result to TR normalized psi_p grid
    x_tr = tr['psip'] / tr['psip'][-1]
    x_g = psi_n_map
    sumavir2_g_on_tr = np.interp(x_tr, x_g, sumavir2_g)

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
    axes[0, 0].plot(x_tr, tr['sumavir2'], 'o-', ms=3, lw=1.8, label='TR SUMAVIR2')
    axes[0, 0].plot(x_tr, sumavir2_g_on_tr, '--', lw=2.0, color='tab:red', label='gfile contour SUMAVIR2')
    axes[0, 0].set_ylabel('SUMAVIR2')
    axes[0, 0].set_title('SUMAVIR2 comparison')
    axes[0, 0].grid(True, alpha=0.3)
    axes[0, 0].legend(fontsize=8)

    ratio_s = tr['sumavir2'] / np.maximum(sumavir2_g_on_tr, 1e-30)
    axes[1, 0].plot(x_tr, ratio_s, 'd-', color='tab:blue', lw=1.5, ms=3)
    axes[1, 0].axhline(1.0, color='k', lw=1, alpha=0.5)
    axes[1, 0].set_xlabel('psi_p / psi_pa')
    axes[1, 0].set_ylabel('TR / gfile')
    axes[1, 0].set_title('SUMAVIR2 ratio')
    axes[1, 0].grid(True, alpha=0.3)

    # show how q mismatch follows SUMAVIR2 mismatch when TTS is fixed
    # compare q on proper psi_p grid
    with open(GFILE, 'r') as f:
        lines = f.readlines()
    header = lines[0].strip().split()
    nw = int(header[-2]); nh = int(header[-1])
    def parse_floats(line):
        res = []
        for i in range(0, len(line), 16):
            chunk = line[i:i + 16].strip()
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
    qpsi = np.abs(np.array(data[4 * nw + nw * nh: 4 * nw + nw * nh + nw]))
    psi_n = np.linspace(0.0, 1.0, nw)
    q_g = np.interp(x_tr, psi_n, qpsi)
    q_ratio = tr['qps'] / np.maximum(q_g, 1e-30)
    axes[0, 1].plot(x_tr, q_ratio, 'o-', ms=3, lw=1.8, label='QPS / q_gfile')
    axes[0, 1].plot(x_tr, ratio_s, '--', lw=2.0, color='tab:green', label='SUMAVIR2 ratio')
    axes[0, 1].axhline(1.0, color='k', lw=1, alpha=0.5)
    axes[0, 1].set_ylabel('ratio')
    axes[0, 1].set_title('QPS mismatch vs SUMAVIR2 mismatch')
    axes[0, 1].grid(True, alpha=0.3)
    axes[0, 1].legend(fontsize=8)

    axes[1, 1].plot(x_tr, tr['tts'], lw=1.8, color='tab:purple', label='TTS')
    axes[1, 1].set_xlabel('psi_p / psi_pa')
    axes[1, 1].set_ylabel('TTS')
    axes[1, 1].set_title('TTS stays smooth')
    axes[1, 1].grid(True, alpha=0.3)
    axes[1, 1].legend(fontsize=8)

    out = os.path.join(BASE, 'sumavir2_from_gfile_comparison.png')
    plt.tight_layout()
    plt.savefig(out, dpi=160)
    plt.close(fig)
    print(f'Wrote {out}')


if __name__ == '__main__':
    main()
