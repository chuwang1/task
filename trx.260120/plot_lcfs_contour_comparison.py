#!/usr/bin/env python3

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path

from plot_sumavir2_from_gfile import read_gfile_2d, read_input_profiles_mapping

BASE = '/Users/dengxiaoya/TASK/CFEDR/git/task/trx.260120'
GFILE = '/Users/dengxiaoya/CFEDRSW/OMFIT_out/g260206.20000_teq_0114'
INPUT_PROFILES = '/Users/dengxiaoya/CFEDRSW/OMFIT_out/input.profiles'


def read_tr_contour(path):
    lines = open(path).read().splitlines()
    R=[]; Z=[]
    for line in lines[2:]:
        parts = line.split(',')
        vals = [float(x) for x in parts[1].split()]
        R.append(vals[0]); Z.append(vals[1])
    return np.array(R), np.array(Z)


def get_gfile_lcfs():
    R, Z, psi_n_2d, _, rmaxis, zmaxis = read_gfile_2d(GFILE)
    psi_n_map, rho_map = read_input_profiles_mapping(INPUT_PROFILES)
    p_target = psi_n_map[-1]
    fig, ax = plt.subplots()
    cs = ax.contour(R, Z, psi_n_2d, levels=[p_target])
    paths=[]
    if len(cs.allsegs)>0 and len(cs.allsegs[0])>0:
        for seg in cs.allsegs[0]:
            paths.append(Path(seg))
    valid=[]
    for p in paths:
        v=p.vertices
        if np.min(v[:,0]) < rmaxis < np.max(v[:,0]) and np.min(v[:,1]) < zmaxis < np.max(v[:,1]):
            valid.append(p)
    if not valid:
        valid=paths
    path=max(valid,key=lambda x: len(x.vertices))
    plt.close(fig)
    v = path.vertices
    if not np.allclose(v[0], v[-1]):
        v = np.vstack([v, v[0]])
    return v[:,0], v[:,1]


def main():
    r_tr, z_tr = read_tr_contour(os.path.join(BASE, 'eqcalq_lcfs_contour.csv'))
    r_g, z_g = get_gfile_lcfs()

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    axes[0].plot(r_tr, z_tr, 'o-', ms=2.5, lw=1.5, label='TR LCFS contour')
    axes[0].plot(r_g, z_g, '--', lw=2.0, color='tab:red', label='gfile LCFS contour')
    axes[0].set_aspect('equal', adjustable='box')
    axes[0].set_xlabel('R [m]')
    axes[0].set_ylabel('Z [m]')
    axes[0].set_title('LCFS contour comparison')
    axes[0].grid(True, alpha=0.3)
    axes[0].legend(fontsize=8)

    axes[1].plot(r_tr, z_tr, 'o-', ms=2.5, lw=1.5, label='TR')
    axes[1].plot(r_g, z_g, '--', lw=2.0, color='tab:red', label='gfile')
    axes[1].set_xlim(min(r_tr.min(), r_g.min())-0.05, max(r_tr.max(), r_g.max())+0.05)
    axes[1].set_ylim(min(z_tr.min(), z_g.min())-0.05, max(z_tr.max(), z_g.max())+0.05)
    axes[1].set_aspect('equal', adjustable='box')
    axes[1].set_xlabel('R [m]')
    axes[1].set_ylabel('Z [m]')
    axes[1].set_title('Same contours, zoomed')
    axes[1].grid(True, alpha=0.3)
    axes[1].legend(fontsize=8)

    out = os.path.join(BASE, 'lcfs_contour_comparison.png')
    plt.tight_layout()
    plt.savefig(out, dpi=170)
    plt.close(fig)
    print(f'Wrote {out}')


if __name__ == '__main__':
    main()
