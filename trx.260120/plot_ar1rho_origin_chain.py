#!/usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

BASE = "/Users/dengxiaoya/TASK/CFEDR/git/task/trx.260120"
OMFIT = "/Users/dengxiaoya/CFEDRSW/OMFIT_out/outputs_new"


def read_two_column_csv(path):
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("Title") or line.startswith("X,"):
                continue
            parts = [p.strip() for p in line.split(",")]
            rows.append([float(x) for x in parts])
    return np.array(rows)


def main():
    tr_ar = read_two_column_csv(f"{BASE}/tr_data_134.csv")
    tr_rmnr = read_two_column_csv(f"{BASE}/tr_data_053.csv")

    rho_tr = tr_ar[:, 0]
    ar1rho_tr = tr_ar[:, 1]
    rmnrho_tr = tr_rmnr[:, 2]
    jac_tr = np.gradient(rho_tr, rmnrho_tr)
    grad_rmin_implied_tr = ar1rho_tr / jac_tr

    extra = pd.read_csv(f"{OMFIT}/input_profiles_extra.csv")
    allp = pd.read_csv(f"{OMFIT}/all_profiles.csv")
    rho_om = allp["rho"].values
    rmin_om = allp["rmin"].values
    ave_grad_r_om = extra["EXPRO_ave_grad_r"].values
    jac_om = np.gradient(rho_om, rmin_om)
    ar1rho_om = ave_grad_r_om * np.abs(jac_om)

    f_rmin_om = interp1d(rho_om, rmin_om, kind="cubic", fill_value="extrapolate")
    f_jac_om = interp1d(rho_om, jac_om, kind="cubic", fill_value="extrapolate")
    f_gradr_om = interp1d(rho_om, ave_grad_r_om, kind="cubic", fill_value="extrapolate")
    f_ar1_om = interp1d(rho_om, ar1rho_om, kind="cubic", fill_value="extrapolate")

    rmin_om_on_tr = f_rmin_om(rho_tr)
    jac_om_on_tr = f_jac_om(rho_tr)
    gradr_om_on_tr = f_gradr_om(rho_tr)
    ar1_om_on_tr = f_ar1_om(rho_tr)

    mask = (rho_tr > 0.0) & (rho_tr < 0.22)

    fig, axes = plt.subplots(4, 2, figsize=(13, 16), sharex="col")

    axes[0, 0].plot(rho_tr[mask], rmnrho_tr[mask], "o-", ms=3, lw=1.5, label="TR RMNRHO")
    axes[0, 0].plot(rho_tr[mask], rmin_om_on_tr[mask], "--", lw=2.0, label="OMFIT rmin")
    axes[0, 0].set_ylabel("r_min [m]")
    axes[0, 0].set_title("Earliest divergence: r_min mapping")
    axes[0, 0].grid(True, alpha=0.3)
    axes[0, 0].legend(fontsize=9)
    axes[0, 1].plot(rho_tr[mask], rmnrho_tr[mask] / np.maximum(rmin_om_on_tr[mask], 1e-30), "d-", ms=3, lw=1.4)
    axes[0, 1].axhline(1.0, color="k", lw=1.0, alpha=0.5)
    axes[0, 1].set_ylabel("TR / OMFIT")
    axes[0, 1].set_title("r_min ratio")
    axes[0, 1].grid(True, alpha=0.3)

    axes[1, 0].plot(rho_tr[mask], jac_tr[mask], "o-", ms=3, lw=1.5, label=r"TR $d\rho/dr_{min}$")
    axes[1, 0].plot(rho_tr[mask], jac_om_on_tr[mask], "--", lw=2.0, label=r"OMFIT $d\rho/dr_{min}$")
    axes[1, 0].set_ylabel(r"$d\rho/dr_{min}$ [1/m]")
    axes[1, 0].set_title("Jacobian inherited from r_min mapping")
    axes[1, 0].grid(True, alpha=0.3)
    axes[1, 0].legend(fontsize=9)
    axes[1, 1].plot(rho_tr[mask], jac_tr[mask] / np.maximum(jac_om_on_tr[mask], 1e-30), "d-", ms=3, lw=1.4)
    axes[1, 1].axhline(1.0, color="k", lw=1.0, alpha=0.5)
    axes[1, 1].set_ylabel("TR / OMFIT")
    axes[1, 1].set_title("Jacobian ratio")
    axes[1, 1].grid(True, alpha=0.3)

    axes[2, 0].plot(rho_tr[mask], grad_rmin_implied_tr[mask], "o-", ms=3, lw=1.5,
                    label=r"TR implied $\langle|\nabla r_{min}|\rangle = AR1RHO / (d\rho/dr_{min})$")
    axes[2, 0].plot(rho_tr[mask], gradr_om_on_tr[mask], "--", lw=2.0,
                    label=r"OMFIT EXPRO_ave_grad_r")
    axes[2, 0].set_ylabel(r"$\langle|\nabla r_{min}|\rangle$")
    axes[2, 0].set_title("After removing Jacobian, core mismatch is much smaller")
    axes[2, 0].grid(True, alpha=0.3)
    axes[2, 0].legend(fontsize=8)
    axes[2, 1].plot(rho_tr[mask], grad_rmin_implied_tr[mask] / np.maximum(gradr_om_on_tr[mask], 1e-30), "d-", ms=3, lw=1.4)
    axes[2, 1].axhline(1.0, color="k", lw=1.0, alpha=0.5)
    axes[2, 1].set_ylabel("TR / OMFIT")
    axes[2, 1].set_title(r"Implied $\langle|\nabla r_{min}|\rangle$ ratio")
    axes[2, 1].grid(True, alpha=0.3)

    axes[3, 0].plot(rho_tr[mask], ar1rho_tr[mask], "o-", ms=3, lw=1.5, label="TR AR1RHO")
    axes[3, 0].plot(rho_tr[mask], ar1_om_on_tr[mask], "--", lw=2.0, label="OMFIT converted AR1RHO")
    axes[3, 0].set_ylabel("AR1RHO")
    axes[3, 0].set_title("Final comparison")
    axes[3, 0].grid(True, alpha=0.3)
    axes[3, 0].legend(fontsize=9)
    axes[3, 1].plot(rho_tr[mask], ar1rho_tr[mask] / np.maximum(ar1_om_on_tr[mask], 1e-30), "d-", ms=3, lw=1.4)
    axes[3, 1].axhline(1.0, color="k", lw=1.0, alpha=0.5)
    axes[3, 1].set_ylabel("TR / OMFIT")
    axes[3, 1].set_title("AR1RHO ratio")
    axes[3, 1].grid(True, alpha=0.3)

    axes[3, 0].set_xlabel(r"$\rho_{tor}$")
    axes[3, 1].set_xlabel(r"$\rho_{tor}$")

    fig.suptitle("AR1RHO origin chain: locate the first core divergence", fontsize=14, y=0.995)
    plt.tight_layout(rect=(0, 0, 1, 0.985))
    out = f"{BASE}/ar1rho_origin_chain_core.png"
    plt.savefig(out, dpi=170)
    plt.close(fig)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
