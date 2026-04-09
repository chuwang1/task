#!/usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path

from plot_sumavir2_from_gfile import read_gfile_2d, read_input_profiles_mapping

BASE = "/Users/dengxiaoya/TASK/CFEDR/git/task/trx.260120"
GFILE = "/Users/dengxiaoya/CFEDRSW/OMFIT_out/g260206.20000_teq_0114"
INPUT = "/Users/dengxiaoya/CFEDRSW/OMFIT_out/input.profiles"
ALLP = "/Users/dengxiaoya/CFEDRSW/OMFIT_out/outputs_new/all_profiles.csv"


def read_tr_rmnrho(path):
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("Title") or line.startswith("X,"):
                continue
            p = [float(x) for x in line.split(",")]
            rows.append((p[0], p[2]))
    arr = np.array(rows)
    return arr[:, 0], arr[:, 1]


def gfile_halfwidth_trace():
    R, Z, psi_n_2d, _, rmaxis, zmaxis = read_gfile_2d(GFILE)
    psi_n_map, rho_map = read_input_profiles_mapping(INPUT)
    fig, ax = plt.subplots()
    rhalf = np.zeros_like(rho_map)
    for i, p_target in enumerate(psi_n_map):
        if p_target <= 1e-6:
            rhalf[i] = 0.0
            continue
        cs = ax.contour(R, Z, psi_n_2d, levels=[p_target])
        paths = []
        if len(cs.allsegs) > 0 and len(cs.allsegs[0]) > 0:
            paths = [Path(seg) for seg in cs.allsegs[0]]
        valid = []
        for p in paths:
            v = p.vertices
            if np.min(v[:, 0]) < rmaxis < np.max(v[:, 0]) and np.min(v[:, 1]) < zmaxis < np.max(v[:, 1]):
                valid.append(p)
        if not valid:
            valid = paths
        path = max(valid, key=lambda x: len(x.vertices))
        v = path.vertices
        rhalf[i] = 0.5 * (np.max(v[:, 0]) - np.min(v[:, 0]))
    plt.close(fig)
    return rho_map, rhalf


def main():
    rho_tr, rmnrho_tr = read_tr_rmnrho(f"{BASE}/tr_data_053.csv")
    rho_g, rhalf_g = gfile_halfwidth_trace()
    allp = pd.read_csv(ALLP)
    rho_om = allp["rho"].values
    rmin_om = allp["rmin"].values

    rhalf_on_tr = np.interp(rho_tr, rho_g, rhalf_g)
    rmin_on_tr = np.interp(rho_tr, rho_om, rmin_om)

    out = pd.DataFrame({
        "rho": rho_tr,
        "TR_RMNRHO": rmnrho_tr,
        "gfile_halfwidth": rhalf_on_tr,
        "OMFIT_rmin": rmin_on_tr,
        "TR_over_gfile_halfwidth": rmnrho_tr / np.maximum(rhalf_on_tr, 1e-30),
        "TR_over_OMFIT_rmin": rmnrho_tr / np.maximum(rmin_on_tr, 1e-30),
        "OMFIT_over_gfile_halfwidth": rmin_on_tr / np.maximum(rhalf_on_tr, 1e-30),
    })
    out_csv = f"{BASE}/rmin_source_trace.csv"
    out.to_csv(out_csv, index=False)

    mask = (rho_tr >= 0.0) & (rho_tr <= 0.25)
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), sharex=True)

    axes[0, 0].plot(rho_tr[mask], rmnrho_tr[mask], "o-", ms=3, lw=1.6, label="TR RMNRHO")
    axes[0, 0].plot(rho_tr[mask], rmin_on_tr[mask], "--", lw=2.0, label="OMFIT rmin")
    axes[0, 0].plot(rho_tr[mask], rhalf_on_tr[mask], ":", lw=2.2, label="gfile contour (Rmax-Rmin)/2")
    axes[0, 0].set_ylabel("r [m]")
    axes[0, 0].set_title("Core minor-radius definitions")
    axes[0, 0].grid(True, alpha=0.3)
    axes[0, 0].legend(fontsize=8)

    axes[0, 1].plot(rho_tr[mask], (rmnrho_tr / np.maximum(rhalf_on_tr, 1e-30))[mask], "o-", ms=3, lw=1.6, label="TR / gfile halfwidth")
    axes[0, 1].plot(rho_tr[mask], (rmin_on_tr / np.maximum(rhalf_on_tr, 1e-30))[mask], "--", lw=2.0, label="OMFIT / gfile halfwidth")
    axes[0, 1].axhline(1.0, color="k", lw=1.0, alpha=0.5)
    axes[0, 1].set_ylabel("ratio")
    axes[0, 1].set_title("Relative to direct gfile half-width")
    axes[0, 1].grid(True, alpha=0.3)
    axes[0, 1].legend(fontsize=8)

    axes[1, 0].plot(rho_tr[mask], (rmnrho_tr - rhalf_on_tr)[mask], "o-", ms=3, lw=1.6, label="TR - gfile halfwidth")
    axes[1, 0].plot(rho_tr[mask], (rmin_on_tr - rhalf_on_tr)[mask], "--", lw=2.0, label="OMFIT - gfile halfwidth")
    axes[1, 0].axhline(0.0, color="k", lw=1.0, alpha=0.5)
    axes[1, 0].set_xlabel(r"$\rho = \sqrt{\psi_{t,n}}$")
    axes[1, 0].set_ylabel("difference [m]")
    axes[1, 0].set_title("Absolute offset from gfile half-width")
    axes[1, 0].grid(True, alpha=0.3)
    axes[1, 0].legend(fontsize=8)

    axes[1, 1].plot(rho_tr[mask], (rmnrho_tr / np.maximum(rmin_on_tr, 1e-30))[mask], "o-", ms=3, lw=1.6)
    axes[1, 1].axhline(1.0, color="k", lw=1.0, alpha=0.5)
    axes[1, 1].set_xlabel(r"$\rho = \sqrt{\psi_{t,n}}$")
    axes[1, 1].set_ylabel("TR / OMFIT")
    axes[1, 1].set_title("Final RMNRHO vs OMFIT rmin ratio")
    axes[1, 1].grid(True, alpha=0.3)

    plt.tight_layout()
    out_png = f"{BASE}/rmin_source_trace.png"
    plt.savefig(out_png, dpi=170)
    plt.close(fig)

    print(f"Wrote {out_png}")
    print(f"Wrote {out_csv}")


if __name__ == "__main__":
    main()
