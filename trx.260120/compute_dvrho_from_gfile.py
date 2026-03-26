#!/usr/bin/env python3
"""
Compute DVRHO directly from gfile by contour-integrated volume.

Algorithm:
1. Read 2D poloidal flux map from gfile.
2. For each flux surface (using rho/psi_n mapping from input.profiles),
   compute enclosed torus volume by contour integration.
3. Differentiate V(rho) numerically to obtain dV/drho.
4. Compare with TR `tr_data_133.csv` if available.

This is the direct geometric route corresponding to the idea behind
    dV/drho = (dV/dpsi_t) * (dpsi_t/drho)
but implemented through V(rho) from 2D flux surfaces.
"""

from __future__ import annotations

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path


TRX_DIR = "/Users/dengxiaoya/TASK/latest/task/trx"
BASE_DIR = "/Users/dengxiaoya/CFEDRSW/OMFIT_out"
GFILE = os.path.join(BASE_DIR, "g260206.20000_teq_0114")
INPUT_PROFILES = os.path.join(BASE_DIR, "input.profiles")


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
    return R, Z, psi_n, rmaxis, zmaxis


def compute_polygon_volume(R_v, Z_v):
    vol = 0.0
    for i in range(len(R_v) - 1):
        vol += (np.pi / 3.0) * (R_v[i] ** 2 + R_v[i] * R_v[i + 1] + R_v[i + 1] ** 2) * (Z_v[i + 1] - Z_v[i])
    return abs(vol)


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


def gradient(x, y):
    out = np.zeros_like(y)
    out[0] = (y[1] - y[0]) / (x[1] - x[0])
    out[-1] = (y[-1] - y[-2]) / (x[-1] - x[-2])
    out[1:-1] = (y[2:] - y[:-2]) / (x[2:] - x[:-2])
    return out


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


def interp(xq, x, y):
    return np.interp(xq, x, y)


def main():
    print("Reading gfile...")
    R, Z, psi_n_2d, rmaxis, zmaxis = read_gfile_2d(GFILE)
    psi_n_map, rho_map = read_input_profiles_mapping(INPUT_PROFILES)

    print(f"Computing volume for {len(rho_map)} flux surfaces...")
    vol_gfile = np.zeros_like(rho_map)
    fig_dummy, ax = plt.subplots()

    for i, p_target in enumerate(psi_n_map):
        if p_target <= 1e-4:
            vol_gfile[i] = 0.0
            continue

        cs = ax.contour(R, Z, psi_n_2d, levels=[p_target])
        paths = []
        if len(cs.allsegs) > 0 and len(cs.allsegs[0]) > 0:
            for seg in cs.allsegs[0]:
                paths.append(Path(seg))

        if not paths:
            vol_gfile[i] = np.nan
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
        vol_gfile[i] = compute_polygon_volume(v[:, 0], v[:, 1])

    plt.close(fig_dummy)

    # Fill any NaN by interpolation
    mask = np.isfinite(vol_gfile)
    vol_gfile = np.interp(rho_map, rho_map[mask], vol_gfile[mask])
    dvrho_gfile = gradient(rho_map, vol_gfile)

    out_csv = os.path.join(TRX_DIR, 'gfile_dvrho.csv')
    import csv
    with open(out_csv, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['rho', 'V_gfile', 'DVRHO_gfile'])
        for i in range(len(rho_map)):
            w.writerow([rho_map[i], vol_gfile[i], dvrho_gfile[i]])
    print(f"Wrote {out_csv}")

    tr_csv = os.path.join(TRX_DIR, 'tr_data_133.csv')
    tr_rho = tr_dvrho = None
    if os.path.exists(tr_csv):
        tr_rho, tr_dvrho = read_tr_dvrho(tr_csv)
        dvrho_on_tr = interp(tr_rho, rho_map, dvrho_gfile)
        rel = np.abs(dvrho_on_tr - tr_dvrho) / np.maximum(np.abs(tr_dvrho), 1e-30)
        print(f"mean |relative error| = {np.mean(rel):.4f}")
        print(f"max  |relative error| = {np.max(rel):.4f}")

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True,
                                   gridspec_kw={'height_ratios': [3, 1]})
    ax1.plot(rho_map, dvrho_gfile, '-', lw=2, label='gfile contour-integrated DVRHO')
    if tr_rho is not None:
        ax1.plot(tr_rho, tr_dvrho, 'o', ms=4, label='TR tr_data_133 DVRHO')
    ax1.set_ylabel('DVRHO [m^3]')
    ax1.set_title('DVRHO from gfile volume integration')
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=9)

    if tr_rho is not None:
        relpct = 100.0 * (interp(tr_rho, rho_map, dvrho_gfile) - tr_dvrho) / np.maximum(np.abs(tr_dvrho), 1e-30)
        ax2.plot(tr_rho, relpct, 'd-', color='tab:red', lw=1.5, ms=3)
        ax2.axhline(0.0, color='k', lw=1, alpha=0.5)
        ax2.set_ylabel('Rel diff [%]')
        ax2.grid(True, alpha=0.3)
    ax2.set_xlabel('rho')

    out_png = os.path.join(TRX_DIR, 'gfile_dvrho_comparison.png')
    plt.tight_layout()
    plt.savefig(out_png, dpi=160)
    plt.close(fig)
    print(f"Wrote {out_png}")


if __name__ == '__main__':
    main()
