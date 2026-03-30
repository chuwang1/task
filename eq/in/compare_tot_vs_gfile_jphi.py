#!/usr/bin/env python3
"""
Compare OMFIT total current-density profile (tot vs rho)
with jphi profile derived from gfile export.
"""

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def metrics(y_ref, y_cmp):
    d = y_cmp - y_ref
    rmse = float(np.sqrt(np.mean(d**2)))
    mae = float(np.mean(np.abs(d)))
    mx = float(np.max(np.abs(d)))
    rel = rmse / max(float(np.max(np.abs(y_ref))), 1e-12)
    return rmse, mae, mx, rel


def main():
    parser = argparse.ArgumentParser(description="Compare OMFIT tot(rho) with gfile jphi.")
    parser.add_argument(
        "--tot-csv",
        default="/Users/dengxiaoya/CFEDRSW/OMFIT_out/current_profiles_vs_rho.csv",
        help="CSV path containing columns: rho, tot",
    )
    parser.add_argument(
        "--gfile-jphi-csv",
        default="in/g260206.20000_teq_0114_pffneg_R803_current_profiles_from_gfile.csv",
        help="CSV from plot_current_profile_from_gfile.py",
    )
    parser.add_argument(
        "--gfile-col",
        default="jphi_ref_raw_A_per_m2",
        help="Column name in gfile-jphi-csv for jphi profile",
    )
    parser.add_argument(
        "--out-prefix",
        default="in/tot_vs_gfile_jphi",
        help="Output prefix for CSV/PNG",
    )
    parser.add_argument(
        "--sqrt-left-shift",
        type=float,
        default=0.0,
        help="Left shift applied to rho=sqrt(psi_n) curve (in rho units).",
    )
    args = parser.parse_args()

    if not os.path.exists(args.tot_csv):
        raise FileNotFoundError(f"tot csv not found: {args.tot_csv}")
    if not os.path.exists(args.gfile_jphi_csv):
        raise FileNotFoundError(f"gfile jphi csv not found: {args.gfile_jphi_csv}")

    dft = pd.read_csv(args.tot_csv)
    dfg = pd.read_csv(args.gfile_jphi_csv)

    if "rho" not in dft.columns or "tot" not in dft.columns:
        raise ValueError("tot-csv must contain columns: rho, tot")
    if "psi_n" not in dfg.columns or args.gfile_col not in dfg.columns:
        raise ValueError(f"gfile-jphi-csv must contain columns: psi_n, {args.gfile_col}")

    rho_t = dft["rho"].to_numpy(dtype=float)
    tot = dft["tot"].to_numpy(dtype=float)

    psi_n = np.clip(dfg["psi_n"].to_numpy(dtype=float), 0.0, 1.0)
    jg = dfg[args.gfile_col].to_numpy(dtype=float)

    rho_lin = psi_n
    rho_sqrt = np.sqrt(psi_n)

    j_lin = np.interp(rho_t, np.sort(rho_lin), jg[np.argsort(rho_lin)])
    j_sqrt = np.interp(rho_t, np.sort(rho_sqrt), jg[np.argsort(rho_sqrt)])
    # left shift by s: y_shifted(x) = y_original(x + s)
    j_sqrt_shift = np.interp(
        np.clip(rho_t + args.sqrt_left_shift, 0.0, 1.0),
        np.sort(rho_sqrt),
        jg[np.argsort(rho_sqrt)],
    )

    m_lin = metrics(tot, j_lin)
    m_sqrt = metrics(tot, j_sqrt)
    m_sqrt_shift = metrics(tot, j_sqrt_shift)

    print("Comparison metrics (A/m^2):")
    print(
        f"  rho=psi_n      : RMSE={m_lin[0]:.4e}, MAE={m_lin[1]:.4e}, "
        f"MAX={m_lin[2]:.4e}, relRMSE={m_lin[3]:.2%}"
    )
    print(
        f"  rho=sqrt(psi_n): RMSE={m_sqrt[0]:.4e}, MAE={m_sqrt[1]:.4e}, "
        f"MAX={m_sqrt[2]:.4e}, relRMSE={m_sqrt[3]:.2%}"
    )
    print(
        f"  rho=sqrt(psi_n), left shift={args.sqrt_left_shift:.3f}: "
        f"RMSE={m_sqrt_shift[0]:.4e}, MAE={m_sqrt_shift[1]:.4e}, "
        f"MAX={m_sqrt_shift[2]:.4e}, relRMSE={m_sqrt_shift[3]:.2%}"
    )

    out_prefix = os.path.abspath(args.out_prefix)
    os.makedirs(os.path.dirname(out_prefix), exist_ok=True)

    out_csv = f"{out_prefix}.csv"
    out_png = f"{out_prefix}.png"

    out = pd.DataFrame(
        {
            "rho": rho_t,
            "tot_A_per_m2": tot,
            "gfile_jphi_rho_eq_psin_A_per_m2": j_lin,
            "gfile_jphi_rho_eq_sqrt_psin_A_per_m2": j_sqrt,
            "gfile_jphi_rho_eq_sqrt_psin_shifted_A_per_m2": j_sqrt_shift,
            "diff_lin_A_per_m2": j_lin - tot,
            "diff_sqrt_A_per_m2": j_sqrt - tot,
            "diff_sqrt_shift_A_per_m2": j_sqrt_shift - tot,
        }
    )
    out.to_csv(out_csv, index=False)

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8), constrained_layout=True)

    axes[0].plot(rho_t, tot / 1e6, "k-", lw=2, label="OMFIT tot")
    axes[0].plot(rho_t, j_lin / 1e6, "b--", lw=1.8, label="gfile jphi (rho=psi_n)")
    axes[0].plot(rho_t, j_sqrt / 1e6, "r-.", lw=1.8, label="gfile jphi (rho=sqrt(psi_n))")
    axes[0].plot(
        rho_t,
        j_sqrt_shift / 1e6,
        color="tab:green",
        lw=1.8,
        linestyle=":",
        label=f"gfile jphi (sqrt, left shift {args.sqrt_left_shift:.2f})",
    )
    axes[0].set_xlabel("rho")
    axes[0].set_ylabel("Current density (MA/m^2)")
    axes[0].set_title("tot vs gfile jphi")
    axes[0].grid(alpha=0.25)
    axes[0].legend(fontsize=8)

    axes[1].plot(rho_t, (j_lin - tot) / 1e6, "b--", lw=1.8, label="gfile-lin - tot")
    axes[1].plot(rho_t, (j_sqrt - tot) / 1e6, "r-.", lw=1.8, label="gfile-sqrt - tot")
    axes[1].plot(
        rho_t,
        (j_sqrt_shift - tot) / 1e6,
        color="tab:green",
        lw=1.8,
        linestyle=":",
        label="gfile-sqrt-shift - tot",
    )
    axes[1].axhline(0.0, color="k", lw=1)
    axes[1].set_xlabel("rho")
    axes[1].set_ylabel("Difference (MA/m^2)")
    axes[1].set_title("Difference")
    axes[1].grid(alpha=0.25)
    axes[1].legend(fontsize=8)

    fig.suptitle("OMFIT tot vs gfile-derived jphi")
    fig.savefig(out_png, dpi=160)
    plt.close(fig)

    print(f"saved: {out_csv}")
    print(f"saved: {out_png}")


if __name__ == "__main__":
    main()
