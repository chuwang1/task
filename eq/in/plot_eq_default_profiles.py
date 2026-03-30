#!/usr/bin/env python3
"""
Plot TASK/EQ default analytic profiles from eqinit.f + eqfunc.f.

This script reproduces the default profile shapes used when MDLEQF=0
and no profile coefficients are overridden in input.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def eqfunc(
    psin: np.ndarray,
    f0: float,
    fs: float,
    f1: float,
    f2: float,
    psiitb: float,
    profr0: float,
    profr1: float,
    profr2: float,
    prof0: float,
    prof1: float,
    prof2: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Python equivalent of Fortran EQFUNC (for psin in [0, 1])."""
    arg0 = np.power(psin, profr0)
    arg1 = np.power(psin, profr1)

    f = fs + (f0 - fs) * np.power(1.0 - arg0, prof0) + f1 * np.power(1.0 - arg1, prof1)
    df = (
        -(f0 - fs) * prof0 * np.power(1.0 - arg0, prof0 - 1.0) * profr0 * np.power(psin, profr0 - 1.0)
        - f1 * prof1 * np.power(1.0 - arg1, prof1 - 1.0) * profr1 * np.power(psin, profr1 - 1.0)
    )

    if psiitb > 0.0:
        mask = psin < psiitb
        if np.any(mask):
            arg2 = np.power(psin[mask] / psiitb, profr2)
            f[mask] += f2 * np.power(1.0 - arg2, prof2)
            df[mask] += (
                -f2
                * prof2
                * np.power(1.0 - arg2, prof2 - 1.0)
                * profr2
                * np.power(psin[mask] / psiitb, profr2 - 1.0)
                / psiitb
            )

    return f, df


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot TASK/EQ default profiles.")
    parser.add_argument("--out-prefix", default="in/eq_default_profiles", help="Output prefix without extension.")
    parser.add_argument("--npts", type=int, default=401, help="Number of psi_N samples.")
    args = parser.parse_args()

    # Defaults from eqinit.f (plus PL_INIT default RHOITB(1)=0).
    rr = 3.0
    bb = 3.0
    q0 = 1.0
    qa = 3.0
    rhoitb = 0.0
    psiitb = rhoitb**2

    pp0, pp1, pp2 = 0.001, 0.0, 0.0
    profp0, profp1, profp2 = 1.5, 1.5, 2.0

    ff0, ff1, ff2 = 1.0, 0.0, 0.0
    proff0, proff1, proff2 = 1.5, 1.5, 1.5

    pj0, pj1, pj2 = 1.0, 0.0, 0.0
    profj0, profj1, profj2 = 1.5, 1.5, 1.5

    profr0, profr1, profr2 = 1.0, 2.0, 2.0

    eps = 1.0e-6
    psin = np.linspace(eps, 1.0, args.npts)

    # P and dP/dpsi_N (MPa space, as in eqinit analytic coefficients).
    p_mpa, dp_mpa_dpsin = eqfunc(
        psin, pp0, 0.0, pp1, pp2, psiitb, profr0, profr1, profr2, profp0, profp1, profp2
    )
    p_pa = p_mpa * 1.0e6
    dp_pa_dpsin = dp_mpa_dpsin * 1.0e6

    # F and dF/dpsi_N.
    fs = 2.0 * np.pi * bb * rr
    f, df_dpsin = eqfunc(
        psin, ff0 + fs, fs, ff1, ff2, psiitb, profr0, profr1, profr2, proff0, proff1, proff2
    )
    ffprime_dpsin = f * df_dpsin

    # HJPSI shape (eqfunc.f::EQJPSI, MDLEQF<5 branch), scaled by 1e6 like Fortran output.
    arg0 = np.power(psin, profr0)
    arg1 = np.power(psin, profr1)
    arg2 = np.power(psin, profr2)
    hjpsi_shape = (
        pj0 * np.power(1.0 - arg0, profj0) * np.power(psin, profr0 - 1.0)
        + pj1 * np.power(1.0 - arg1, profj1) * np.power(psin, profr1 - 1.0)
        + pj2 * np.power(1.0 - arg2, profj2) * np.power(psin, profr2 - 1.0)
    ) * 1.0e6

    # Model q profile used in EQQPSI default analytic branch (if PSITN~PSIN).
    q_model = q0 + (qa - q0) * psin

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    csv_path = out_prefix.with_suffix(".csv")
    with csv_path.open("w", newline="") as fcsv:
        writer = csv.writer(fcsv)
        writer.writerow(
            [
                "psi_n",
                "p_mpa",
                "dp_dpsi_n_mpa",
                "p_pa",
                "dp_dpsi_n_pa",
                "F",
                "dF_dpsi_n",
                "FFprime_dpsi_n",
                "hjpsi_shape_x1e6",
                "q_model_linear",
            ]
        )
        for i in range(len(psin)):
            writer.writerow(
                [
                    f"{psin[i]:.10e}",
                    f"{p_mpa[i]:.10e}",
                    f"{dp_mpa_dpsin[i]:.10e}",
                    f"{p_pa[i]:.10e}",
                    f"{dp_pa_dpsin[i]:.10e}",
                    f"{f[i]:.10e}",
                    f"{df_dpsin[i]:.10e}",
                    f"{ffprime_dpsin[i]:.10e}",
                    f"{hjpsi_shape[i]:.10e}",
                    f"{q_model[i]:.10e}",
                ]
            )

    fig, axs = plt.subplots(2, 3, figsize=(13, 8), constrained_layout=True)
    axs = axs.ravel()

    axs[0].plot(psin, p_mpa, lw=2)
    axs[0].set_title("P(psi_N) [MPa]")
    axs[0].set_xlabel("psi_N")
    axs[0].grid(True, alpha=0.3)

    axs[1].plot(psin, dp_mpa_dpsin, lw=2)
    axs[1].set_title("dP/dpsi_N [MPa]")
    axs[1].set_xlabel("psi_N")
    axs[1].grid(True, alpha=0.3)

    axs[2].plot(psin, f, lw=2)
    axs[2].set_title("F(psi_N) = R*Bt")
    axs[2].set_xlabel("psi_N")
    axs[2].grid(True, alpha=0.3)

    axs[3].plot(psin, df_dpsin, lw=2)
    axs[3].set_title("dF/dpsi_N")
    axs[3].set_xlabel("psi_N")
    axs[3].grid(True, alpha=0.3)

    axs[4].plot(psin, ffprime_dpsin, lw=2)
    axs[4].set_title("F * dF/dpsi_N")
    axs[4].set_xlabel("psi_N")
    axs[4].grid(True, alpha=0.3)

    ax6 = axs[5]
    ax6.plot(psin, hjpsi_shape, lw=2, color="tab:blue", label="HJPSI shape")
    ax6.set_title("HJPSI shape + q model")
    ax6.set_xlabel("psi_N")
    ax6.set_ylabel("HJPSI (x1e6 scale)", color="tab:blue")
    ax6.tick_params(axis="y", labelcolor="tab:blue")
    ax6.grid(True, alpha=0.3)
    ax6r = ax6.twinx()
    ax6r.plot(psin, q_model, lw=2, color="tab:orange", label="q model (linear)")
    ax6r.set_ylabel("q", color="tab:orange")
    ax6r.tick_params(axis="y", labelcolor="tab:orange")

    fig.suptitle("TASK/EQ Default Analytic Profiles (MDLEQF=0)", fontsize=14)
    png_path = out_prefix.with_suffix(".png")
    fig.savefig(png_path, dpi=150)
    plt.close(fig)

    print(f"Saved CSV: {csv_path}")
    print(f"Saved PNG: {png_path}")


if __name__ == "__main__":
    main()
