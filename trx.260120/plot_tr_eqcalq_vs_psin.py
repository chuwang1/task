#!/usr/bin/env python3

import os
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from plot_sumavir2_from_gfile import read_input_profiles_mapping
from plot_equilibrium_chain_big import read_tr_eqcalq


BASE = "/Users/dengxiaoya/TASK/CFEDR/git/task/trx.260120"
OMFIT_OUT = "/Users/dengxiaoya/CFEDRSW/OMFIT_out"
INPUT_PROFILES = os.path.join(OMFIT_OUT, "input.profiles")
EQCALQ = os.path.join(BASE, "eqcalq_dvdpsit_raw.csv")


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
        raise FileNotFoundError("No radial TR CSV file found.")

    x_tr_csv = read_tr_csv_x(csv_path)
    tr = read_tr_eqcalq(EQCALQ)
    psi_n_g, rhoa_g = read_input_profiles_mapping(INPUT_PROFILES)
    sqrt_psin_g = np.sqrt(np.clip(psi_n_g, 0.0, None))

    psip_n_tr = tr["psip"] / tr["psip"][-1]
    psit_n_tr = tr["psit"] / tr["psit"][-1]
    sqrt_psip_tr = np.sqrt(np.clip(psip_n_tr, 0.0, None))
    sqrt_psit_tr = np.sqrt(np.clip(psit_n_tr, 0.0, None))
    rho_tr = tr["rho"]

    # Interpolate everything onto TR CSV X for direct comparison
    sqrt_psip_on_csv = np.interp(x_tr_csv, rho_tr, sqrt_psip_tr)
    sqrt_psit_on_csv = np.interp(x_tr_csv, rho_tr, sqrt_psit_tr)
    sqrt_psin_g_on_csv = np.interp(x_tr_csv, rhoa_g, sqrt_psin_g)

    out = pd.DataFrame({
        "x_tr_csv": x_tr_csv,
        "sqrt_psip_tr_on_csv": sqrt_psip_on_csv,
        "sqrt_psit_tr_on_csv": sqrt_psit_on_csv,
        "sqrt_psin_gfile_on_csv": sqrt_psin_g_on_csv,
        "x_minus_sqrt_psip_tr": x_tr_csv - sqrt_psip_on_csv,
        "x_minus_sqrt_psit_tr": x_tr_csv - sqrt_psit_on_csv,
        "x_minus_sqrt_psin_gfile": x_tr_csv - sqrt_psin_g_on_csv,
    })
    csv_out = os.path.join(BASE, "tr_eqcalq_vs_psin.csv")
    out.to_csv(csv_out, index=False)

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    ax = axes[0, 0]
    ax.plot(x_tr_csv, x_tr_csv, "k-", lw=2.0, label="TR CSV X")
    ax.plot(rho_tr, sqrt_psip_tr, "b--", lw=2.0, label=r"TR $\sqrt{\psi_{p,n}}$")
    ax.plot(rho_tr, sqrt_psit_tr, "g-.", lw=2.0, label=r"TR $\sqrt{\psi_{t,n}}$")
    ax.plot(rhoa_g, sqrt_psin_g, "r-", lw=1.8, label=r"gfile/input $\sqrt{\psi_n}$")
    ax.set_xlabel("radial coordinate")
    ax.set_ylabel("value")
    ax.set_title("TR CSV X vs psi-based coordinates")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    ax.plot(x_tr_csv, x_tr_csv - sqrt_psip_on_csv, "b-", lw=2.0, label=r"$X_{TR}-\sqrt{\psi_{p,n}}$")
    ax.plot(x_tr_csv, x_tr_csv - sqrt_psit_on_csv, "g-", lw=2.0, label=r"$X_{TR}-\sqrt{\psi_{t,n}}$")
    ax.plot(x_tr_csv, x_tr_csv - sqrt_psin_g_on_csv, "r-", lw=2.0, label=r"$X_{TR}-\sqrt{\psi_n}$ gfile")
    ax.axhline(0.0, color="k", lw=1.0, alpha=0.5)
    ax.set_xlabel(r"$X_{TR}$")
    ax.set_ylabel("difference")
    ax.set_title("Difference on TR CSV grid")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    ax = axes[1, 0]
    ax.plot(rho_tr, psip_n_tr, "b--", lw=2.0, label=r"TR $\psi_{p,n}$")
    ax.plot(rho_tr, psit_n_tr, "g-.", lw=2.0, label=r"TR $\psi_{t,n}$")
    ax.plot(rhoa_g, psi_n_g, "r-", lw=1.8, label=r"gfile/input $\psi_n$")
    ax.set_xlabel(r"$\rho$ / $X_{TR}$")
    ax.set_ylabel("normalized flux")
    ax.set_title("Normalized flux profiles")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    ax = axes[1, 1]
    ax.plot(rho_tr, np.gradient(sqrt_psip_tr, rho_tr), "b--", lw=2.0, label=r"$d\sqrt{\psi_{p,n}}/d\rho_{TR}$")
    ax.plot(rho_tr, np.gradient(sqrt_psit_tr, rho_tr), "g-.", lw=2.0, label=r"$d\sqrt{\psi_{t,n}}/d\rho_{TR}$")
    ax.plot(rhoa_g, np.gradient(sqrt_psin_g, rhoa_g), "r-", lw=1.8, label=r"$d\sqrt{\psi_n}/d(r/a)$")
    ax.axhline(1.0, color="k", lw=1.0, alpha=0.5)
    ax.set_xlabel("radial coordinate")
    ax.set_ylabel("slope")
    ax.set_title("Slope comparison")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    plt.tight_layout()
    fig_out = os.path.join(BASE, "tr_eqcalq_vs_psin.png")
    plt.savefig(fig_out, dpi=170)
    plt.close(fig)

    print(f"TR radial CSV: {os.path.basename(csv_path)}")
    print(f"Title: {csv_title}")
    print(f"Wrote {fig_out}")
    print(f"Wrote {csv_out}")
    print(f"max |X_TR - sqrt(psip_n)_TR| = {np.max(np.abs(x_tr_csv - sqrt_psip_on_csv)):.6f}")
    print(f"max |X_TR - sqrt(psit_n)_TR| = {np.max(np.abs(x_tr_csv - sqrt_psit_on_csv)):.6f}")
    print(f"max |X_TR - sqrt(psi_n)_gfile| = {np.max(np.abs(x_tr_csv - sqrt_psin_g_on_csv)):.6f}")


if __name__ == "__main__":
    main()
