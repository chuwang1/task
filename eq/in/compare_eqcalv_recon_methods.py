#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import plot_eqipqp_all_variables_big as pb
import reconstruct_eqcalv_from_gfile as rc
from solve_gs_from_gfile import read_gfile


def _resolve_path(path_str: str, in_dir: Path, repo_dir: Path) -> Path:
    p = Path(path_str)
    if p.is_absolute() and p.exists():
        return p
    candidates = [Path.cwd() / p, in_dir / p, repo_dir / p]
    for c in candidates:
        if c.exists():
            return c.resolve()
    return (in_dir / p).resolve()


def _metrics(a: np.ndarray, b: np.ndarray) -> tuple[float, float, float, float, int]:
    m = np.isfinite(a) & np.isfinite(b)
    if np.count_nonzero(m) < 3:
        return np.nan, np.nan, np.nan, np.nan, int(np.count_nonzero(m))
    d = a[m] - b[m]
    rmse = float(np.sqrt(np.mean(d * d)))
    mae = float(np.mean(np.abs(d)))
    max_abs = float(np.max(np.abs(d)))
    rel_rmse = rmse / max(float(np.max(np.abs(a[m]))), 1e-12)
    return rmse, mae, max_abs, rel_rmse, int(np.count_nonzero(m))


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare two gfile-based EQCALV reconstruction methods.")
    parser.add_argument("--prefix", default="eqdata0114", help="eqdata prefix (for parameters/radial csv).")
    parser.add_argument("--gfile", default="g260206.20000_teq_0114", help="Path to gfile.")
    parser.add_argument(
        "--gfile-profile-csv",
        default="g260206.20000_teq_0114_profiles_from_gfile.csv",
        help="Path to gfile profile csv used by plot_eqipqp method.",
    )
    parser.add_argument("--f-sign", type=float, default=1.0)
    parser.add_argument("--pprime-sign", type=float, default=1.0)
    parser.add_argument(
        "--out-prefix",
        default=None,
        help="Output prefix; default: in/<prefix>_eqcalv_method_compare",
    )
    args = parser.parse_args()

    in_dir = Path(__file__).resolve().parent
    repo_dir = in_dir.parent

    params_csv = repo_dir / f"{args.prefix}_parameters.csv"
    radial_csv = repo_dir / f"{args.prefix}_radial_profiles.csv"
    if not params_csv.exists():
        raise FileNotFoundError(f"Missing parameters csv: {params_csv}")
    if not radial_csv.exists():
        raise FileNotFoundError(f"Missing radial csv: {radial_csv}")

    gfile_path = _resolve_path(args.gfile, in_dir, repo_dir)
    profile_csv_path = _resolve_path(args.gfile_profile_csv, in_dir, repo_dir)
    if not gfile_path.exists():
        raise FileNotFoundError(f"Missing gfile: {gfile_path}")
    if not profile_csv_path.exists():
        raise FileNotFoundError(f"Missing gfile profile csv: {profile_csv_path}")

    out_prefix = args.out_prefix
    if out_prefix is None:
        out_prefix = str(in_dir / f"{args.prefix}_eqcalv_method_compare")
    out_prefix_path = Path(out_prefix)
    out_prefix_path.parent.mkdir(parents=True, exist_ok=True)

    params = pd.read_csv(params_csv)
    pmap = dict(zip(params["Parameter"], params["Value"]))
    rr = float(pmap["RR"])
    bb = float(pmap["BB"])
    psip_max = float(pmap["PSIPA"])

    dr = pd.read_csv(radial_csv)
    psi_n = np.clip(dr["PSIPNV"].to_numpy(float), 0.0, 1.0)
    psip = dr["PSIPV"].to_numpy(float)

    g = read_gfile(str(gfile_path))
    g_prof = pb._load_gfile_profiles(
        profile_csv_path,
        g,
        psip_max=psip_max,
        f_sign=args.f_sign,
        pprime_sign=args.pprime_sign,
    )
    if g_prof is None:
        raise RuntimeError("Failed to load gfile profile data for plot_eqipqp method.")

    q_levels = pb._safe_interp(g_prof["psi_n"], g_prof["q"], psi_n)
    f_levels = pb._safe_interp(g_prof["psi_n"], g_prof["f"], psi_n)
    geom_plot = pb._compute_eqcalv_like_from_gfile(
        g,
        psi_n_levels=psi_n,
        psip_levels=psip,
        q_levels=q_levels,
        f_levels=f_levels,
        rr=rr,
        bb=bb,
    )

    g_rc = rc.read_gfile(str(gfile_path))
    geom_rc = rc.compute_profiles(g_rc, npsi=len(psi_n))

    b_rsv = np.interp(psi_n, geom_rc["psi_n"], geom_rc["RSV"])
    b_avir2 = np.interp(psi_n, geom_rc["psi_n"], geom_rc["AVIR2"])
    b_vpv = np.interp(psi_n, geom_rc["psi_n"], geom_rc["VPV"])
    b_avrr2 = np.interp(psi_n, geom_rc["psi_n"], geom_rc["AVRR2"])

    a_rsv = np.asarray(geom_plot["rsv"], dtype=float)
    a_avir2 = np.asarray(geom_plot["avir2"], dtype=float)
    a_vpv = np.asarray(geom_plot["vpv"], dtype=float)
    a_avrr2 = np.asarray(geom_plot["avrr2"], dtype=float)

    rows = [
        ("RSV", a_rsv, b_rsv),
        ("AVIR2", a_avir2, b_avir2),
        ("VPV", a_vpv, b_vpv),
        ("AVRR2", a_avrr2, b_avrr2),
    ]

    fig, axes = plt.subplots(4, 2, figsize=(14, 16), constrained_layout=True)
    for i, (name, a, b) in enumerate(rows):
        rmse, mae, max_abs, rel_rmse, n_ok = _metrics(a, b)
        d = a - b

        ax0 = axes[i, 0]
        pb._plot_line(ax0, psi_n, a, "b-", lw=1.8, label="plot_eqipqp_all_variables_big")
        pb._plot_line(ax0, psi_n, b, "r--", lw=1.8, label="reconstruct_eqcalv_from_gfile")
        ax0.set_title(f"{name}  RMSE={rmse:.3e}  rel={rel_rmse:.3e}  N={n_ok}")
        ax0.set_xlabel("psi_n")
        ax0.set_ylabel(name)
        ax0.grid(alpha=0.3)
        if i == 0:
            ax0.legend(fontsize=8)

        ax1 = axes[i, 1]
        pb._plot_line(ax1, psi_n, d, "k-", lw=1.5)
        ax1.axhline(0.0, color="gray", lw=1.0, ls="--")
        ax1.set_title(f"{name} delta (plot - recon), MAE={mae:.3e}, MAX={max_abs:.3e}")
        ax1.set_xlabel("psi_n")
        ax1.set_ylabel("delta")
        ax1.grid(alpha=0.3)

    fig.suptitle("EQCALV-like geometry reconstruction: method comparison", fontsize=15)
    out_png = out_prefix_path.with_suffix(".png")
    fig.savefig(out_png, dpi=170)

    out_csv = out_prefix_path.with_suffix(".csv")
    df_out = pd.DataFrame(
        {
            "psi_n": psi_n,
            "RSV_plot": a_rsv,
            "RSV_recon": b_rsv,
            "RSV_delta": a_rsv - b_rsv,
            "AVIR2_plot": a_avir2,
            "AVIR2_recon": b_avir2,
            "AVIR2_delta": a_avir2 - b_avir2,
            "VPV_plot": a_vpv,
            "VPV_recon": b_vpv,
            "VPV_delta": a_vpv - b_vpv,
            "AVRR2_plot": a_avrr2,
            "AVRR2_recon": b_avrr2,
            "AVRR2_delta": a_avrr2 - b_avrr2,
        }
    )
    df_out.to_csv(out_csv, index=False)

    print(f"Saved comparison figure: {out_png}")
    print(f"Saved comparison csv: {out_csv}")


if __name__ == "__main__":
    main()
