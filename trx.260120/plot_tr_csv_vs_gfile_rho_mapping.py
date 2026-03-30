#!/usr/bin/env python3

import os
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from plot_sumavir2_from_gfile import read_input_profiles_mapping


BASE = "/Users/dengxiaoya/TASK/CFEDR/git/task/trx.260120"
OMFIT_OUT = "/Users/dengxiaoya/CFEDRSW/OMFIT_out"
INPUT_PROFILES = os.path.join(OMFIT_OUT, "input.profiles")


def find_radial_csv(search_dir):
    for path in sorted(glob.glob(os.path.join(search_dir, "tr_data_*.csv"))):
        with open(path, "r") as fh:
            title = fh.readline().strip()
        if "vs r@" in title:
            return path, title
    return None, None


def read_tr_csv_x(path):
    df = pd.read_csv(path, skiprows=1)
    df.columns = [c.strip() for c in df.columns]
    return df["X"].to_numpy()


def main():
    csv_path, csv_title = find_radial_csv(BASE)
    if csv_path is None:
        raise FileNotFoundError("No radial TR CSV file found (tr_data_*.csv with 'vs r@').")

    x_tr = read_tr_csv_x(csv_path)
    psi_n_map, rhoa_g = read_input_profiles_mapping(INPUT_PROFILES)
    sqrt_psin_g = np.sqrt(np.clip(psi_n_map, 0.0, None))

    # Interpolate gfile/input.profiles mapping onto TR CSV x-grid for direct comparison
    sqrt_psin_on_tr = np.interp(x_tr, rhoa_g, sqrt_psin_g)
    diff = x_tr - sqrt_psin_on_tr

    out_df = pd.DataFrame({
        "x_tr_csv": x_tr,
        "r_over_a_gfile_mapping": np.interp(x_tr, rhoa_g, rhoa_g),
        "sqrt_psi_n_gfile_mapping": sqrt_psin_on_tr,
        "xtr_minus_sqrtpsin": diff,
    })
    csv_out = os.path.join(BASE, "tr_csv_vs_gfile_rho_mapping.csv")
    out_df.to_csv(csv_out, index=False)

    fig, axes = plt.subplots(2, 2, figsize=(11, 8))

    ax = axes[0, 0]
    ax.plot(x_tr, x_tr, "k-", lw=2.0, label="TR CSV X")
    ax.plot(rhoa_g, rhoa_g, "b--", lw=1.6, label="gfile/input.profiles r/a")
    ax.plot(rhoa_g, sqrt_psin_g, "r-", lw=2.0, label=r"gfile/input.profiles $\sqrt{\psi_n}$")
    ax.set_xlabel("reference x")
    ax.set_ylabel("value")
    ax.set_title("TR CSV grid vs gfile/input.profiles mapping")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    ax.plot(x_tr, sqrt_psin_on_tr, "r-", lw=2.0, label=r"$\sqrt{\psi_n}$ on TR grid")
    ax.plot(x_tr, x_tr, "k--", lw=1.2, label="TR CSV X")
    ax.set_xlabel("TR CSV X")
    ax.set_ylabel("value")
    ax.set_title(r"Compare TR CSV X and $\sqrt{\psi_n}$")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    ax = axes[1, 0]
    ax.plot(x_tr, diff, "m-", lw=2.0)
    ax.axhline(0.0, color="k", lw=1.0, alpha=0.5)
    ax.set_xlabel("TR CSV X")
    ax.set_ylabel(r"TR X $-$ $\sqrt{\psi_n}$")
    ax.set_title("Difference on TR grid")
    ax.grid(True, alpha=0.3)

    ax = axes[1, 1]
    dsdx = np.gradient(sqrt_psin_on_tr, x_tr)
    ax.plot(x_tr, dsdx, "g-", lw=2.0)
    ax.axhline(1.0, color="k", lw=1.0, alpha=0.5)
    ax.set_xlabel("TR CSV X")
    ax.set_ylabel(r"$d\sqrt{\psi_n}/dX_{TR}$")
    ax.set_title("Slope on TR grid")
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    fig_out = os.path.join(BASE, "tr_csv_vs_gfile_rho_mapping.png")
    plt.savefig(fig_out, dpi=170)
    plt.close(fig)

    idx = int(np.nanargmax(np.abs(diff)))
    print(f"TR radial CSV: {os.path.basename(csv_path)}")
    print(f"Title: {csv_title}")
    print(f"Wrote {fig_out}")
    print(f"Wrote {csv_out}")
    print(f"max |TR X - sqrt(psi_n)| = {abs(diff[idx]):.6f} at TR X = {x_tr[idx]:.3f}")
    print(f"mean |TR X - sqrt(psi_n)| = {np.mean(np.abs(diff)):.6f}")


if __name__ == "__main__":
    main()
