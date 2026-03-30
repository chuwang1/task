#!/usr/bin/env python3
"""Rebuild F from gfile-derived geometry and q, then compare to gfile F.

The geometry inputs (RSV/VPV/AVIR2/AVRR2) are computed from gfile contours via
`reconstruct_eqcalv_from_gfile.py`, which uses the factor=1 approximation.
"""

from __future__ import annotations

import argparse
import csv
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from reconstruct_eqcalv_from_gfile import compute_profiles, read_gfile

MU0 = 4.0e-7 * np.pi
TWO_PI = 2.0 * np.pi
FOUR_PI2 = 4.0 * np.pi * np.pi
EPS = 1.0e-14


def _signed_clip_den(x: float, eps: float = EPS) -> float:
    if not np.isfinite(x):
        return np.nan
    if abs(x) >= eps:
        return x
    return eps if x >= 0.0 else -eps


def _interp_from_unit_grid(y: np.ndarray, x_new: np.ndarray) -> np.ndarray:
    x_old = np.linspace(0.0, 1.0, y.size)
    return np.interp(x_new, x_old, y)


def rebuild_f_eqipqp(
    psip: np.ndarray,
    q: np.ndarray,
    pprime: np.ndarray,
    avir2: np.ndarray,
    vpv: np.ndarray,
    avrr2: np.ndarray,
    rsv: np.ndarray,
    rr: float,
    bb: float,
) -> np.ndarray:
    n = len(psip)
    f = np.full(n, np.nan, dtype=float)
    if n < 2:
        return f

    f_sign = 1.0 if bb >= 0.0 else -1.0
    f[-1] = TWO_PI * bb * rr

    for i in range(n - 1, 0, -1):
        if not np.isfinite(f[i]):
            continue

        vals = [
            psip[i],
            psip[i - 1],
            q[i],
            q[i - 1],
            pprime[i],
            pprime[i - 1],
            avir2[i],
            avir2[i - 1],
            vpv[i],
            vpv[i - 1],
            avrr2[i],
            avrr2[i - 1],
            rsv[i],
            rsv[i - 1],
        ]
        if not np.isfinite(vals).all():
            continue
        if abs(q[i]) < EPS or abs(q[i - 1]) < EPS:
            continue

        xp = rsv[i] / q[i]
        xm = rsv[i - 1] / q[i - 1]

        den_alp = _signed_clip_den(avir2[i] * vpv[i])
        if i == 1:
            den_alm = _signed_clip_den(avir2[i - 1] * vpv[i])
            alm_num = xp
        else:
            den_alm = _signed_clip_den(avir2[i - 1] * vpv[i - 1])
            alm_num = xm

        den_blp = _signed_clip_den(avir2[i])
        den_blm = _signed_clip_den(avir2[i - 1])
        if not np.isfinite([den_alp, den_alm, den_blp, den_blm]).all():
            continue

        alp = FOUR_PI2 * bb * bb * rr * rr * xp / den_alp
        alm = FOUR_PI2 * bb * bb * rr * rr * alm_num / den_alm
        blp = FOUR_PI2 * MU0 * rr * rr * pprime[i] / den_blp
        blm = FOUR_PI2 * MU0 * rr * rr * pprime[i - 1] / den_blm
        flp = vpv[i] * avrr2[i] * xp
        flm = vpv[i - 1] * avrr2[i - 1] * xm

        yp = 0.5 * f[i] * f[i]
        dpsi = psip[i] - psip[i - 1]
        ym = yp + 0.5 * (blp + blm) * dpsi + 0.5 * (alp + alm) * (flp - flm)
        if not np.isfinite(ym) or ym <= 0.0:
            continue
        f[i - 1] = f_sign * np.sqrt(2.0 * ym)

    return f


def save_csv(path: str, out: dict[str, np.ndarray]) -> None:
    cols = [
        "index",
        "psi_n",
        "psi_wb",
        "q",
        "pprime",
        "RSV",
        "VPV",
        "AVIR2",
        "AVRR2",
        "F_rebuilt",
        "F_gfile",
        "dF",
        "rel_err",
        "contour_ok",
    ]
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        n = len(out["index"])
        for i in range(n):
            w.writerow(
                [
                    int(out["index"][i]),
                    float(out["psi_n"][i]),
                    float(out["psi_wb"][i]),
                    float(out["q"][i]),
                    float(out["pprime"][i]),
                    float(out["RSV"][i]),
                    float(out["VPV"][i]),
                    float(out["AVIR2"][i]),
                    float(out["AVRR2"][i]),
                    float(out["F_rebuilt"][i]) if np.isfinite(out["F_rebuilt"][i]) else np.nan,
                    float(out["F_gfile"][i]) if np.isfinite(out["F_gfile"][i]) else np.nan,
                    float(out["dF"][i]) if np.isfinite(out["dF"][i]) else np.nan,
                    float(out["rel_err"][i]) if np.isfinite(out["rel_err"][i]) else np.nan,
                    bool(out["contour_ok"][i]),
                ]
            )


def save_plot(path: str, out: dict[str, np.ndarray], title_tag: str) -> None:
    psi_n = out["psi_n"]
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
    ax = axes.ravel()

    ax[0].plot(psi_n, out["F_gfile"], lw=1.8, label="F_gfile")
    ax[0].plot(psi_n, out["F_rebuilt"], lw=1.6, ls="--", label="F_rebuilt")
    ax[0].set_title("F profile")
    ax[0].grid(True, alpha=0.3)
    ax[0].legend(loc="best")

    ax[1].plot(psi_n, out["dF"], lw=1.7)
    ax[1].set_title("dF = F_rebuilt - F_gfile")
    ax[1].grid(True, alpha=0.3)

    ax[2].plot(psi_n, 100.0 * out["rel_err"], lw=1.7)
    ax[2].set_title("Relative Error [%]")
    ax[2].grid(True, alpha=0.3)

    ax[3].plot(psi_n, out["q"], lw=1.7, label="q")
    ax[3].plot(psi_n, out["pprime"], lw=1.2, label="pprime")
    bad = np.where(~out["contour_ok"])[0]
    if bad.size > 0:
        ax[3].scatter(psi_n[bad], out["q"][bad], s=16, color="tab:red", alpha=0.6, label="contour failed")
    ax[3].set_title("q and pprime")
    ax[3].grid(True, alpha=0.3)
    ax[3].legend(loc="best")

    ax[2].set_xlabel("psi_n")
    ax[3].set_xlabel("psi_n")
    fig.suptitle(title_tag)
    fig.tight_layout()
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare EQIPQP-style rebuilt F against gfile F.")
    parser.add_argument("--gfile", default="in/g260206.20000_teq_0114", help="Input GEQDSK gfile path")
    parser.add_argument(
        "--out-prefix",
        default="in/compare_f_from_gfile_geom_q",
        help="Output prefix; writes <prefix>.csv and <prefix>.png",
    )
    parser.add_argument("--npsi", type=int, default=None, help="Number of psi_n points (default: gfile qpsi size)")
    args = parser.parse_args()

    gfile_path = os.path.abspath(args.gfile)
    if not os.path.exists(gfile_path):
        raise FileNotFoundError(f"Gfile not found: {gfile_path}")

    out_prefix = args.out_prefix
    out_dir = os.path.dirname(out_prefix) if os.path.dirname(out_prefix) else "."
    os.makedirs(out_dir, exist_ok=True)

    g = read_gfile(gfile_path)
    npsi = int(args.npsi) if args.npsi is not None else int(g.qpsi.size)
    geom = compute_profiles(g, npsi)

    psi_n = geom["psi_n"]
    # GEQDSK pprime is dp/d(psirz) in Pa/(Wb/rad).
    # EQIPQP formulas expect dp/d(Psi_p) in Pa/Wb. So divide by 2*pi.
    pprime = _interp_from_unit_grid(g.pprime, psi_n) / TWO_PI
    f_gfile = _interp_from_unit_grid(TWO_PI * g.fpol, psi_n)
    f_rebuilt = rebuild_f_eqipqp(
        psip=geom["psi_wb"],
        q=geom["q"],
        pprime=pprime,
        avir2=geom["AVIR2"],
        vpv=geom["VPV"],
        avrr2=geom["AVRR2"],
        rsv=geom["RSV"],
        rr=float(g.rcentr),
        bb=float(g.bcentr),
    )

    d_f = f_rebuilt - f_gfile
    rel = np.abs(d_f) / np.maximum(np.abs(f_gfile), 1.0e-12)
    valid = np.isfinite(f_rebuilt) & np.isfinite(f_gfile)
    n_valid = int(np.count_nonzero(valid))
    rmse = float(np.sqrt(np.mean((d_f[valid]) ** 2))) if n_valid > 0 else np.nan
    mae = float(np.mean(np.abs(d_f[valid]))) if n_valid > 0 else np.nan
    mre = float(np.mean(rel[valid])) if n_valid > 0 else np.nan

    out = {
        "index": geom["index"],
        "psi_n": psi_n,
        "psi_wb": geom["psi_wb"],
        "q": geom["q"],
        "pprime": pprime,
        "RSV": geom["RSV"],
        "VPV": geom["VPV"],
        "AVIR2": geom["AVIR2"],
        "AVRR2": geom["AVRR2"],
        "F_rebuilt": f_rebuilt,
        "F_gfile": f_gfile,
        "dF": d_f,
        "rel_err": rel,
        "contour_ok": geom["contour_ok"],
    }

    csv_path = f"{out_prefix}.csv"
    png_path = f"{out_prefix}.png"
    save_csv(csv_path, out)
    save_plot(png_path, out, title_tag="F rebuilt from gfile geometry + q")

    print(f"Gfile: {gfile_path}")
    print(f"N valid: {n_valid}/{len(psi_n)}")
    print(f"RMSE: {rmse:.6e}")
    print(f"MAE: {mae:.6e}")
    print(f"Mean Relative Error: {mre:.6e}")
    print(f"Output CSV: {os.path.abspath(csv_path)}")
    print(f"Output PNG: {os.path.abspath(png_path)}")


if __name__ == "__main__":
    main()
