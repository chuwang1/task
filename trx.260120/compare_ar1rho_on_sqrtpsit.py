#!/usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

BASE = "/Users/dengxiaoya/TASK/CFEDR/git/task/trx.260120"
OMFIT = "/Users/dengxiaoya/CFEDRSW/OMFIT_out/outputs_new"
ARHO_EXP = 3.6146142


def read_tr_curve(path):
    rows = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("Title") or line.startswith("X,"):
                continue
            a, b = [float(x) for x in line.split(",")]
            rows.append((a, b))
    arr = np.array(rows)
    return arr[:, 0], arr[:, 1]


def main():
    s_tr, ar1rho_tr = read_tr_curve(f"{BASE}/tr_data_134.csv")

    extra = pd.read_csv(f"{OMFIT}/input_profiles_extra.csv")
    allp = pd.read_csv(f"{OMFIT}/all_profiles.csv")

    s_om = allp["rho"].values
    rmin_om = allp["rmin"].values
    ave_grad_r_om = extra["EXPRO_ave_grad_r"].values

    drho_drmin_num = np.gradient(s_om, rmin_om)
    drho_drmin_expro = 1.0 / (extra["EXPRO_drdrho"].values * ARHO_EXP)

    ar1rho_om_num = ave_grad_r_om * np.abs(drho_drmin_num)
    ar1rho_om_expro = ave_grad_r_om * np.abs(drho_drmin_expro)

    f_num = interp1d(s_om, ar1rho_om_num, kind="cubic", fill_value="extrapolate")
    f_expro = interp1d(s_om, ar1rho_om_expro, kind="cubic", fill_value="extrapolate")
    om_num_on_tr = f_num(s_tr)
    om_expro_on_tr = f_expro(s_tr)

    mask = (s_tr > 0.0) & (s_tr < 0.22)
    ratio_num = ar1rho_tr / np.maximum(om_num_on_tr, 1e-30)
    ratio_expro = ar1rho_tr / np.maximum(om_expro_on_tr, 1e-30)

    out = pd.DataFrame({
        "sqrt_psit_n": s_tr,
        "TR_AR1RHO": ar1rho_tr,
        "OMFIT_AR1RHO_num_on_TRs": om_num_on_tr,
        "OMFIT_AR1RHO_expro_on_TRs": om_expro_on_tr,
        "TR_over_OMFIT_num": ratio_num,
        "TR_over_OMFIT_expro": ratio_expro,
    })
    out_csv = f"{BASE}/compare_ar1rho_on_sqrtpsit.csv"
    out.to_csv(out_csv, index=False)

    fig, axes = plt.subplots(2, 2, figsize=(12, 9), sharex=True)

    ax = axes[0, 0]
    ax.plot(s_tr[mask], ar1rho_tr[mask], "o-", ms=3, lw=1.5, label="TR AR1RHO")
    ax.plot(s_tr[mask], om_num_on_tr[mask], "-", lw=2.0, label="OMFIT converted (numerical)")
    ax.plot(s_tr[mask], om_expro_on_tr[mask], "--", lw=2.0, label="OMFIT converted (EXPRO_drdrho)")
    ax.set_ylabel("AR1RHO")
    ax.set_title(r"AR1RHO on common x-axis $\sqrt{\psi_{t,n}}$")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    ax.plot(s_tr[mask], ratio_num[mask], "-", lw=1.8, label="TR / OMFIT (numerical)")
    ax.plot(s_tr[mask], ratio_expro[mask], "--", lw=1.8, label="TR / OMFIT (EXPRO_drdrho)")
    ax.axhline(1.0, color="k", lw=1.0, alpha=0.5)
    ax.set_ylabel("ratio")
    ax.set_title("Ratio on common x-axis")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    ax = axes[1, 0]
    ax.plot(s_tr[mask], ar1rho_tr[mask] - om_num_on_tr[mask], "-", lw=1.8, label="TR - OMFIT (numerical)")
    ax.plot(s_tr[mask], ar1rho_tr[mask] - om_expro_on_tr[mask], "--", lw=1.8, label="TR - OMFIT (EXPRO_drdrho)")
    ax.axhline(0.0, color="k", lw=1.0, alpha=0.5)
    ax.set_xlabel(r"$\sqrt{\psi_{t,n}}$")
    ax.set_ylabel("difference")
    ax.set_title("Absolute difference")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    ax = axes[1, 1]
    ax.plot(s_om[s_om < 0.22], drho_drmin_num[s_om < 0.22], "-", lw=2.0, label=r"$d\rho/dr_{min}$ numerical")
    ax.plot(s_om[s_om < 0.22], drho_drmin_expro[s_om < 0.22], "--", lw=2.0, label=r"$1/(dr_{min}/d\rho)$ from EXPRO")
    ax.set_xlabel(r"$\sqrt{\psi_{t,n}}$")
    ax.set_ylabel(r"$d\rho/dr_{min}$ [1/m]")
    ax.set_title("OMFIT Jacobian consistency")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    plt.tight_layout()
    out_png = f"{BASE}/compare_ar1rho_on_sqrtpsit.png"
    plt.savefig(out_png, dpi=170)
    plt.close(fig)

    print(f"Wrote {out_png}")
    print(f"Wrote {out_csv}")
    core = (s_tr > 0.0) & (s_tr < 0.2)
    print(f"core max |TR/OMFIT_num - 1| = {np.max(np.abs(ratio_num[core] - 1.0)):.6f}")
    print(f"core max |TR/OMFIT_expro - 1| = {np.max(np.abs(ratio_expro[core] - 1.0)):.6f}")


if __name__ == "__main__":
    main()
