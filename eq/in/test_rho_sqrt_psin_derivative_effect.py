#!/usr/bin/env python3
"""
Test derivative-coordinate assumption impact on jphi comparison:
  A) x = psi_N
  B) x = rho, with psi_N = x^2
"""

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from solve_gs_from_gfile import MU0


def safe_grad(y, x):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    idx = np.argsort(x)
    xs = x[idx]
    ys = y[idx]
    xu, iu = np.unique(xs, return_index=True)
    yu = ys[iu]
    if len(xu) < 2:
        return np.zeros_like(y)
    gu = np.gradient(yu, xu)
    gs = np.interp(xs, xu, gu)
    out = np.empty_like(gs)
    out[idx] = gs
    return out


def calc_jphi_from_pf(
    x,
    p,
    f,
    psi_axis,
    psi_bdy,
    r_ref,
    coord_mode="psin",
    p_sign=1,
    f_sign=1,
):
    """
    coord_mode:
      - psin: x == psi_N
      - rho_sqrt_psin: x == rho, psi_N = x^2
    """
    dpsi = float(psi_bdy - psi_axis)
    x = np.asarray(x, dtype=float)
    p = p_sign * np.asarray(p, dtype=float)
    f = f_sign * np.asarray(f, dtype=float)

    dpdx = safe_grad(p, x)
    dfdx = safe_grad(f, x)

    if coord_mode == "psin":
        dpsi_dx = np.full_like(x, dpsi)
    elif coord_mode == "rho_sqrt_psin":
        dpsi_dx = 2.0 * x * dpsi
        if len(dpsi_dx) >= 2:
            dpsi_dx[0] = dpsi_dx[1]
        eps = 1e-10
        dpsi_dx = np.where(np.abs(dpsi_dx) < eps, np.sign(dpsi) * eps, dpsi_dx)
    else:
        raise ValueError(coord_mode)

    pprime = dpdx / dpsi_dx
    ffprime = f * (dfdx / dpsi_dx)
    jphi = r_ref * pprime + ffprime / (MU0 * r_ref)
    return pprime, ffprime, jphi


def metric(y_ref, y_cmp):
    d = y_cmp - y_ref
    rmse = float(np.sqrt(np.mean(d**2)))
    mae = float(np.mean(np.abs(d)))
    mx = float(np.max(np.abs(d)))
    rel = rmse / max(float(np.max(np.abs(y_ref))), 1e-12)
    return rmse, mae, mx, rel


def main():
    parser = argparse.ArgumentParser(description="Test rho=sqrt(psi_N) derivative effect.")
    parser.add_argument(
        "--tot-csv",
        default="/Users/dengxiaoya/CFEDRSW/OMFIT_out/current_profiles_vs_rho.csv",
    )
    parser.add_argument(
        "--gfile-profile-csv",
        default="in/g260206.20000_teq_0114_profiles_from_gfile.csv",
    )
    parser.add_argument("--r-ref", type=float, default=8.03)
    parser.add_argument("--p-sign", type=int, choices=[-1, 1], default=-1)
    parser.add_argument("--f-sign", type=int, choices=[-1, 1], default=-1)
    parser.add_argument("--out-prefix", default="in/jphi_derivative_coord_test")
    args = parser.parse_args()

    if not os.path.exists(args.tot_csv):
        raise FileNotFoundError(args.tot_csv)
    if not os.path.exists(args.gfile_profile_csv):
        raise FileNotFoundError(args.gfile_profile_csv)

    dft = pd.read_csv(args.tot_csv)
    dfg = pd.read_csv(args.gfile_profile_csv)

    rho = dft["rho"].to_numpy(dtype=float)
    tot = dft["tot"].to_numpy(dtype=float)

    x = dfg["psi_n"].to_numpy(dtype=float)  # this x is tested as psi_N or rho
    p = dfg["p_pa"].to_numpy(dtype=float)
    f = dfg["F_tesla_meter"].to_numpy(dtype=float)
    psi_axis = float(dfg["psi_wb"].iloc[0])
    psi_bdy = float(dfg["psi_wb"].iloc[-1])

    _, _, j_psin = calc_jphi_from_pf(
        x, p, f, psi_axis, psi_bdy, args.r_ref, "psin", args.p_sign, args.f_sign
    )
    _, _, j_rho = calc_jphi_from_pf(
        x, p, f, psi_axis, psi_bdy, args.r_ref, "rho_sqrt_psin", args.p_sign, args.f_sign
    )

    j_psin_i = np.interp(rho, x, j_psin)
    j_rho_i = np.interp(rho, x, j_rho)

    m_psin = metric(tot, j_psin_i)
    m_rho = metric(tot, j_rho_i)

    print(
        f"settings: R_ref={args.r_ref}, p_sign={args.p_sign:+d}, f_sign={args.f_sign:+d}"
    )
    print(
        f"x=psi_N      : RMSE={m_psin[0]:.4e}, MAE={m_psin[1]:.4e}, "
        f"MAX={m_psin[2]:.4e}, relRMSE={m_psin[3]:.2%}"
    )
    print(
        f"x=rho(sqrt)  : RMSE={m_rho[0]:.4e}, MAE={m_rho[1]:.4e}, "
        f"MAX={m_rho[2]:.4e}, relRMSE={m_rho[3]:.2%}"
    )

    out_prefix = os.path.abspath(args.out_prefix)
    os.makedirs(os.path.dirname(out_prefix), exist_ok=True)
    out_csv = f"{out_prefix}.csv"
    out_png = f"{out_prefix}.png"

    out = pd.DataFrame(
        {
            "rho": rho,
            "tot_A_per_m2": tot,
            "jphi_assume_x_eq_psin_A_per_m2": j_psin_i,
            "jphi_assume_x_eq_rho_A_per_m2": j_rho_i,
            "diff_psin_minus_tot": j_psin_i - tot,
            "diff_rho_minus_tot": j_rho_i - tot,
        }
    )
    out.to_csv(out_csv, index=False)

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8), constrained_layout=True)
    axes[0].plot(rho, tot / 1e6, "k-", lw=2, label="OMFIT tot")
    axes[0].plot(rho, j_psin_i / 1e6, "b--", lw=1.8, label="jphi (x=psi_N)")
    axes[0].plot(rho, j_rho_i / 1e6, "r-.", lw=1.8, label="jphi (x=rho, psi_N=x^2)")
    axes[0].set_xlabel("rho")
    axes[0].set_ylabel("Current density (MA/m^2)")
    axes[0].set_title("tot vs jphi")
    axes[0].set_xlim(0.2, 1.05)
    axes[0].set_ylim(0., 3)
    axes[0].grid(alpha=0.25)
    axes[0].legend(fontsize=8)

    axes[1].plot(rho, (j_psin_i - tot) / 1e6, "b--", lw=1.8, label="x=psi_N - tot")
    axes[1].plot(rho, (j_rho_i - tot) / 1e6, "r-.", lw=1.8, label="x=rho - tot")
    axes[1].axhline(0.0, color="k", lw=1)
    axes[1].set_xlabel("rho")
    axes[1].set_ylabel("Difference (MA/m^2)")
    axes[1].set_title("Difference")
    axes[1].grid(alpha=0.25)
    axes[1].legend(fontsize=8)

    fig.suptitle("Derivative Coordinate Assumption Test")
    fig.savefig(out_png, dpi=160)
    plt.close(fig)

    print(f"saved: {out_csv}")
    print(f"saved: {out_png}")


if __name__ == "__main__":
    main()
