#!/usr/bin/env python3
"""
GS residual diagnostic on (R,Z) grid.

Residual definition:
  Res(R,Z) = Delta*psi + mu0*R^2*p'(psi) + F(psi)F'(psi)
where ideal consistency gives Res ~= 0 in plasma region.
"""

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from solve_gs_from_gfile import (
    MU0,
    build_lcfs_polygon_mask,
    interpolate_source,
    load_source_profiles,
    read_gfile,
    resolve_profile_csv,
    to_plot_zr,
)


def load_psi_matrix(path, nw, nh):
    """Load psi matrix saved as CSV in [R,Z] or [Z,R] layout."""
    arr = pd.read_csv(path).to_numpy(dtype=float)
    if arr.shape == (nw, nh):
        return arr
    if arr.shape == (nh, nw):
        return arr.T
    raise ValueError(f"Unexpected psi CSV shape {arr.shape}, expected ({nw},{nh}) or ({nh},{nw})")


def gs_operator(psi, r, z):
    """Compute Delta*psi on interior points."""
    dr = float(r[1] - r[0])
    dz = float(z[1] - z[0])
    out = np.full_like(psi, np.nan, dtype=float)

    rc = r[1:-1, None]
    term_rr = (psi[:-2, 1:-1] - 2.0 * psi[1:-1, 1:-1] + psi[2:, 1:-1]) / (dr**2)
    term_r = -(psi[2:, 1:-1] - psi[:-2, 1:-1]) / (2.0 * dr * rc)
    term_zz = (psi[1:-1, :-2] - 2.0 * psi[1:-1, 1:-1] + psi[1:-1, 2:]) / (dz**2)
    out[1:-1, 1:-1] = term_rr + term_r + term_zz
    return out


def apply_source_mask(mask_mode, psi_mask, lcfs_mask):
    use_psi = mask_mode in ("psi", "both")
    use_lcfs = mask_mode in ("lcfs", "both")

    if mask_mode == "none":
        return np.ones_like(psi_mask, dtype=bool)
    m = np.ones_like(psi_mask, dtype=bool)
    if use_psi:
        m &= psi_mask
    if use_lcfs:
        m &= lcfs_mask
    return m


def residual_stats(name, res, interior_mask, lcfs_mask):
    def summarize(tag, m):
        vals = res[m]
        vals = vals[np.isfinite(vals)]
        if vals.size == 0:
            print(f"{name} [{tag}] no valid points")
            return
        rmse = float(np.sqrt(np.mean(vals**2)))
        mae = float(np.mean(np.abs(vals)))
        mx = float(np.max(np.abs(vals)))
        mean = float(np.mean(vals))
        print(
            f"{name} [{tag}] mean={mean:.3e}, rmse={rmse:.3e}, "
            f"mae={mae:.3e}, max|res|={mx:.3e}, n={vals.size}"
        )

    summarize("interior", interior_mask)
    summarize("lcfs_inside", interior_mask & lcfs_mask)


def plot_residual_figure(out_png, g, psi, rhs, res, lcfs_mask):
    r2d, z2d = np.meshgrid(g.R, g.Z)
    interior = np.zeros_like(psi, dtype=bool)
    interior[1:-1, 1:-1] = True

    rhs_in = np.where(interior, rhs, np.nan)
    res_in = np.where(interior, res, np.nan)
    rel = np.where(interior, np.abs(res) / (np.abs(rhs) + 1e-12), np.nan)
    log_rel = np.log10(np.clip(rel, 1e-16, 1e16))

    fig, axes = plt.subplots(2, 2, figsize=(14, 10), constrained_layout=True)

    c0 = axes[0, 0].contourf(r2d, z2d, to_plot_zr(psi, g.R, g.Z), levels=40, cmap="jet")
    fig.colorbar(c0, ax=axes[0, 0])
    axes[0, 0].plot(g.rbbbs, g.zbbbs, "w--", lw=1.2, label="LCFS")
    axes[0, 0].set_title("psi(R,Z)")
    axes[0, 0].set_xlabel("R (m)")
    axes[0, 0].set_ylabel("Z (m)")
    axes[0, 0].set_aspect("equal")
    axes[0, 0].legend(fontsize=8)

    c1 = axes[0, 1].contourf(r2d, z2d, to_plot_zr(rhs_in, g.R, g.Z), levels=40, cmap="viridis")
    fig.colorbar(c1, ax=axes[0, 1])
    axes[0, 1].plot(g.rbbbs, g.zbbbs, "w--", lw=1.2)
    axes[0, 1].set_title("GS RHS = -mu0 R^2 p' - FF'")
    axes[0, 1].set_xlabel("R (m)")
    axes[0, 1].set_ylabel("Z (m)")
    axes[0, 1].set_aspect("equal")

    vmax = np.nanpercentile(np.abs(res_in), 99)
    if not np.isfinite(vmax) or vmax <= 0:
        vmax = 1.0
    c2 = axes[1, 0].contourf(
        r2d,
        z2d,
        to_plot_zr(res_in, g.R, g.Z),
        levels=40,
        cmap="coolwarm",
        vmin=-vmax,
        vmax=vmax,
    )
    fig.colorbar(c2, ax=axes[1, 0])
    axes[1, 0].plot(g.rbbbs, g.zbbbs, "k--", lw=1.0)
    axes[1, 0].set_title("Residual: Delta*psi + mu0 R^2 p' + FF'")
    axes[1, 0].set_xlabel("R (m)")
    axes[1, 0].set_ylabel("Z (m)")
    axes[1, 0].set_aspect("equal")

    c3 = axes[1, 1].contourf(
        r2d, z2d, to_plot_zr(log_rel, g.R, g.Z), levels=40, cmap="magma"
    )
    fig.colorbar(c3, ax=axes[1, 1], label="log10(|res|/(|rhs|+eps))")
    axes[1, 1].plot(g.rbbbs, g.zbbbs, "w--", lw=1.0)
    axes[1, 1].set_title("Relative Residual (log10)")
    axes[1, 1].set_xlabel("R (m)")
    axes[1, 1].set_ylabel("Z (m)")
    axes[1, 1].set_aspect("equal")

    fig.suptitle("Grad-Shafranov Residual Diagnostic")
    fig.savefig(out_png, dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="GS residual diagnostic in R,Z.")
    parser.add_argument("--gfile", default="in/g260206.20000_teq_0114")
    parser.add_argument("--profile-csv", default=None)
    parser.add_argument(
        "--psi-csv",
        default=None,
        help="Optional psi matrix CSV to diagnose (default: use gfile psirz).",
    )
    parser.add_argument(
        "--mask-mode",
        choices=["none", "psi", "lcfs", "both"],
        default="both",
        help="Source mask for p'/FF': none / psi / lcfs / both.",
    )
    parser.add_argument(
        "--source-scale",
        type=float,
        default=1.0,
        help="Scale factor applied to GS source term.",
    )
    parser.add_argument("--out-prefix", default="in/gs_residual")
    args = parser.parse_args()

    gpath = os.path.abspath(args.gfile)
    if not os.path.exists(gpath):
        raise FileNotFoundError(f"gfile not found: {gpath}")

    g = read_gfile(gpath)
    prof_csv = resolve_profile_csv(gpath, args.profile_csv)
    psi_n_src, pprime_src, ffprime_src, src_desc = load_source_profiles(g, prof_csv)

    if args.psi_csv:
        psi = load_psi_matrix(args.psi_csv, g.nw, g.nh)
        psi_name = os.path.abspath(args.psi_csv)
    else:
        psi = g.psirz.copy()
        psi_name = "gfile_psirz"

    lcfs_mask = build_lcfs_polygon_mask(g.R, g.Z, g.rbbbs, g.zbbbs)
    if lcfs_mask is None:
        lcfs_mask = np.zeros_like(psi, dtype=bool)

    pprime, ffprime, psi_inside = interpolate_source(
        psi,
        g.psimag,
        g.psibdy,
        psi_n_src,
        pprime_src,
        ffprime_src,
        lcfs_mask=False,
    )
    src_mask = apply_source_mask(args.mask_mode, psi_inside, lcfs_mask)
    pprime = np.where(src_mask, pprime, 0.0)
    ffprime = np.where(src_mask, ffprime, 0.0)

    lhs = gs_operator(psi, g.R, g.Z)
    rhs = args.source_scale * (-MU0 * (g.R[:, None] ** 2) * pprime - ffprime)
    res = lhs - rhs  # == Delta*psi + mu0 R^2 p' + FF'

    interior_mask = np.zeros_like(psi, dtype=bool)
    interior_mask[1:-1, 1:-1] = True

    print(f"gfile: {gpath}")
    print(f"profile source: {src_desc}")
    print(f"psi source: {psi_name}")
    print(f"mask_mode: {args.mask_mode}, source_scale: {args.source_scale}")
    print(f"lcfs_inside_frac={np.mean(lcfs_mask):.3f}, src_active_frac={np.mean(src_mask):.3f}")

    residual_stats("GS residual", res, interior_mask, lcfs_mask)

    out_prefix = os.path.abspath(args.out_prefix)
    out_csv = f"{out_prefix}_map.csv"
    out_png = f"{out_prefix}_figure.png"
    os.makedirs(os.path.dirname(out_prefix), exist_ok=True)
    pd.DataFrame(res).to_csv(out_csv, index=False)
    plot_residual_figure(out_png, g, psi, rhs, res, lcfs_mask)

    print(f"saved: {out_csv}")
    print(f"saved: {out_png}")


if __name__ == "__main__":
    main()
