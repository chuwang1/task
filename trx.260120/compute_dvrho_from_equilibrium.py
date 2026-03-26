#!/usr/bin/env python3
"""
Compute TR-like DVRHO from equilibrium-derived profile files.

This script extracts the algorithm used in `tr_bpsd_get`:

    rho = sqrt(psit / psita)
    dpsit/drho = 2 * psita * rho
    dV/drho = (dV/dpsit) * (dpsit/drho)

In the current workflow we usually do not have `dV/dpsit` directly in a CSV,
but we do have equilibrium-derived profile outputs from OMFIT/GACODE:

    EXPRO_volp = dV/drmin
    all_profiles.csv gives rho and rmin

So the equivalent chain-rule form is:

    dV/drho = (dV/drmin) * (drmin/drho)

This script computes DVRHO on the equilibrium profile grid and optionally
compares it with TR output `tr_data_133.csv`.
"""

from __future__ import annotations

import csv
import math
import os


TRX_DIR = "/Users/dengxiaoya/TASK/CFEDR/git/task/trx.260120"
OMFIT_DIR = "/Users/dengxiaoya/CFEDRSW/OMFIT_out/outputs_new"


def read_csv_dict(path: str):
    with open(path) as f:
        return list(csv.DictReader(f))


def gradient(x, y):
    n = len(x)
    out = [0.0] * n
    if n < 2:
        return out
    out[0] = (y[1] - y[0]) / (x[1] - x[0])
    for i in range(1, n - 1):
        out[i] = (y[i + 1] - y[i - 1]) / (x[i + 1] - x[i - 1])
    out[-1] = (y[-1] - y[-2]) / (x[-1] - x[-2])
    return out


def interp(xq, x, y):
    if xq <= x[0]:
        return y[0]
    if xq >= x[-1]:
        return y[-1]
    for i in range(len(x) - 1):
        if x[i] <= xq <= x[i + 1]:
            t = (xq - x[i]) / (x[i + 1] - x[i])
            return y[i] * (1.0 - t) + y[i + 1] * t
    return y[-1]


def read_tr_dvrho(path: str):
    with open(path) as f:
        lines = f.readlines()
    rho = []
    dvrho = []
    for line in lines[2:]:
        a, b = [float(s) for s in line.strip().split(",")]
        rho.append(a)
        dvrho.append(b)
    return rho, dvrho


def main():
    all_profiles = read_csv_dict(os.path.join(OMFIT_DIR, "all_profiles.csv"))
    extra = read_csv_dict(os.path.join(OMFIT_DIR, "input_profiles_extra.csv"))

    rho = [float(r["rho"]) for r in all_profiles]
    rmin = [float(r["rmin"]) for r in all_profiles]
    volp = [float(r["EXPRO_volp"]) for r in extra]  # dV/drmin

    # Chain rule equivalent to tr_bpsd_get metric path
    drmin_drho = gradient(rho, rmin)
    dvrho_equ = [volp[i] * drmin_drho[i] for i in range(len(rho))]

    # Cross-check by direct derivative of volume, if available
    dvrho_from_vol = None
    if "EXPRO_vol" in extra[0]:
        vol = [float(r["EXPRO_vol"]) for r in extra]
        dvrho_from_vol = gradient(rho, vol)

    out = os.path.join(TRX_DIR, "equilibrium_dvrho.csv")
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        header = ["rho", "rmin", "drmin_drho", "EXPRO_volp_dVdr", "DVRHO_eq"]
        if dvrho_from_vol is not None:
            header.append("DVRHO_from_EXPRO_vol")
        w.writerow(header)
        for i in range(len(rho)):
            row = [rho[i], rmin[i], drmin_drho[i], volp[i], dvrho_equ[i]]
            if dvrho_from_vol is not None:
                row.append(dvrho_from_vol[i])
            w.writerow(row)

    print(f"Wrote {out}")
    print("Equilibrium DVRHO samples:")
    for i in [0, 20, 40, 80, 120, 160, 200]:
        if i < len(rho):
            msg = (
                f"rho={rho[i]:.3f} rmin={rmin[i]:.4f} drmin/drho={drmin_drho[i]:.4f} "
                f"volp={volp[i]:.4f} DVRHO={dvrho_equ[i]:.4f}"
            )
            if dvrho_from_vol is not None:
                msg += f" d(EXPRO_vol)/drho={dvrho_from_vol[i]:.4f}"
            print(msg)

    # Optional comparison to TR output
    tr_csv = os.path.join(TRX_DIR, "tr_data_133.csv")
    if os.path.exists(tr_csv):
        tr_rho, tr_dvrho = read_tr_dvrho(tr_csv)
        print("\nCompare equilibrium DVRHO to TR tr_data_133.csv:")
        abs_rel = []
        for x, ytr in zip(tr_rho, tr_dvrho):
            yeq = interp(x, rho, dvrho_equ)
            rel = abs(yeq - ytr) / max(abs(ytr), 1e-30)
            abs_rel.append(rel)
        print(f"mean |relative error| = {sum(abs_rel)/len(abs_rel):.4f}")
        print(f"max  |relative error| = {max(abs_rel):.4f}")

        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt

            tr_eq_interp = [interp(x, rho, dvrho_equ) for x in tr_rho]
            rel_pct = [100.0 * (tr_eq_interp[i] - tr_dvrho[i]) / max(abs(tr_dvrho[i]), 1e-30)
                       for i in range(len(tr_rho))]

            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True,
                                           gridspec_kw={"height_ratios": [3, 1]})
            ax1.plot(rho, dvrho_equ, '-', lw=2, label='Equilibrium DVRHO from volp*dr/drho')
            if dvrho_from_vol is not None:
                ax1.plot(rho, dvrho_from_vol, '--', lw=1.8, label='Equilibrium DVRHO from d(EXPRO_vol)/drho')
            ax1.plot(tr_rho, tr_dvrho, 'o', ms=4, label='TR tr_data_133 DVRHO')
            ax1.set_ylabel('DVRHO [m^3]')
            ax1.set_title('DVRHO Comparison: Equilibrium vs TR')
            ax1.grid(True, alpha=0.3)
            ax1.legend(fontsize=9)

            ax2.plot(tr_rho, rel_pct, 'd-', color='tab:red', lw=1.5, ms=3)
            ax2.axhline(0.0, color='k', lw=1, alpha=0.5)
            ax2.set_xlabel('rho')
            ax2.set_ylabel('Rel diff [%]')
            ax2.grid(True, alpha=0.3)

            outplot = os.path.join(TRX_DIR, 'equilibrium_dvrho_comparison.png')
            plt.tight_layout()
            plt.savefig(outplot, dpi=160)
            plt.close(fig)
            print(f"Wrote {outplot}")
        except Exception as e:
            print(f"Plot skipped: {e}")


if __name__ == "__main__":
    main()
