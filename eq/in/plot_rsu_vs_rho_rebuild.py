#!/usr/bin/env python3
"""
Overlay original RSU/ZSU boundary and rebuilt rho(theta) boundary.

This reproduces the logic in EQSET_RHOB_FROM_RSU:
1) convert (R, Z) -> (theta, rho) around RR
2) sort by theta and merge near-duplicate theta
3) linearly interpolate rho(theta)
4) rebuild boundary with R = RR + rho*cos(theta), Z = rho*sin(theta)
"""

import argparse
import csv
import math
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def _close_curve(r, z):
    if len(r) == 0:
        return r, z
    if abs(r[0] - r[-1]) < 1e-12 and abs(z[0] - z[-1]) < 1e-12:
        return r, z
    return np.append(r, r[0]), np.append(z, z[0])


def _read_rr_from_parameters_csv(path):
    if not path.exists():
        return None
    with path.open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            if row.get("parameter", "").strip() == "RR":
                return float(row["value"])
    return None


def _build_rho_theta_model(rsu, zsu, rr, eps=1e-8):
    dx = rsu - rr
    dz = zsu
    rho = np.hypot(dx, dz)
    theta = np.arctan2(dz, dx)
    theta = np.where(theta < 0.0, theta + 2.0 * np.pi, theta)

    mask = np.isfinite(theta) & np.isfinite(rho) & (rho > eps)
    theta = theta[mask]
    rho = rho[mask]

    order = np.argsort(theta)
    theta_s = theta[order]
    rho_s = rho[order]

    theta_u = [theta_s[0]]
    rho_u = [rho_s[0]]
    for t, r in zip(theta_s[1:], rho_s[1:]):
        if t - theta_u[-1] > eps:
            theta_u.append(t)
            rho_u.append(r)
        else:
            rho_u[-1] = 0.5 * (rho_u[-1] + r)
    theta_u = np.asarray(theta_u)
    rho_u = np.asarray(rho_u)

    theta_p = np.concatenate(([theta_u[-1] - 2.0 * np.pi], theta_u, [theta_u[0] + 2.0 * np.pi]))
    rho_p = np.concatenate(([rho_u[-1]], rho_u, [rho_u[0]]))

    return theta, rho, theta_u, rho_u, theta_p, rho_p


def _rho_interp(theta, theta_p, rho_p):
    t = np.mod(theta, 2.0 * np.pi)
    return np.interp(t, theta_p, rho_p)


def _nearest_distance(a_xy, b_xy):
    d2 = np.sum((a_xy[:, None, :] - b_xy[None, :, :]) ** 2, axis=2)
    return np.sqrt(np.min(d2, axis=1))


def main():
    parser = argparse.ArgumentParser(description="Compare original RSU boundary vs rebuilt rho(theta) boundary.")
    parser.add_argument("--gfile", default="in/g260206.20000_teq_0114", help="Path to gfile.")
    parser.add_argument(
        "--prefix",
        default="eqdata_modelg5_rebuild_mdleqf9_qmap",
        help="eqdata prefix, used only to auto-load RR from <prefix>_parameters.csv.",
    )
    parser.add_argument("--rr", type=float, default=None, help="Override RR (m).")
    parser.add_argument("--ntgmax", type=int, default=32, help="Poloidal mesh count used for THGM/THGG sampling.")
    parser.add_argument("--out", default="in/rsu_vs_rho_rebuild_overlay.png", help="Output figure path.")
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    eq_dir = script_dir.parent
    gfile_dir = script_dir.parents[3] / "gfile"
    if str(gfile_dir) not in sys.path:
        sys.path.append(str(gfile_dir))

    from geqdsk import geqdsk

    gfile_path = Path(args.gfile)
    if not gfile_path.is_absolute():
        gfile_path = (eq_dir / args.gfile).resolve() if not Path(args.gfile).exists() else Path(args.gfile).resolve()
    g = geqdsk(str(gfile_path))
    rsu = np.asarray(g.rbbbs, dtype=float)
    zsu = np.asarray(g.zbbbs, dtype=float)

    rr = args.rr
    rr_source = "cli --rr"
    if rr is None:
        p_csv = eq_dir / f"{args.prefix}_parameters.csv"
        rr = _read_rr_from_parameters_csv(p_csv)
        rr_source = str(p_csv)
    if rr is None:
        rr = 8.03
        rr_source = "fallback 8.03"

    theta_raw, rho_raw, theta_u, rho_u, theta_p, rho_p = _build_rho_theta_model(rsu, zsu, rr)

    n_dense = max(8 * args.ntgmax, 720)
    th_dense = np.linspace(0.0, 2.0 * np.pi, n_dense, endpoint=False)
    rho_dense = _rho_interp(th_dense, theta_p, rho_p)
    r_dense = rr + rho_dense * np.cos(th_dense)
    z_dense = rho_dense * np.sin(th_dense)

    dtg = 2.0 * np.pi / float(args.ntgmax)
    thgm = dtg * (np.arange(args.ntgmax) + 0.5)
    thgg = dtg * np.arange(args.ntgmax)
    rhgm = _rho_interp(thgm, theta_p, rho_p)
    rhgg = _rho_interp(thgg, theta_p, rho_p)
    rgm = rr + rhgm * np.cos(thgm)
    zgm = rhgm * np.sin(thgm)
    rgg = rr + rhgg * np.cos(thgg)
    zgg = rhgg * np.sin(thgg)

    # Pointwise radial-map residual: raw point vs rho(theta_raw) mapped point.
    rho_hat = _rho_interp(theta_raw, theta_p, rho_p)
    r_hat = rr + rho_hat * np.cos(theta_raw)
    z_hat = rho_hat * np.sin(theta_raw)
    d_map = np.hypot((rr + rho_raw * np.cos(theta_raw)) - r_hat, (rho_raw * np.sin(theta_raw)) - z_hat)

    # Geometric nearest residual: original RSU polyline vs rebuilt dense polyline.
    raw_xy = np.column_stack([rsu, zsu])
    rec_xy = np.column_stack([r_dense, z_dense])
    d_near = _nearest_distance(raw_xy, rec_xy)

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(14, 6), constrained_layout=True)

    r0, z0 = _close_curve(rsu, zsu)
    r1, z1 = _close_curve(r_dense, z_dense)
    ax0.plot(r0, z0, color="tab:blue", lw=1.8, label="Original RSU/ZSU")
    ax0.plot(r1, z1, color="tab:red", lw=1.8, ls="--", label="Rebuilt rho(theta)")
    ax0.plot(rgm, zgm, "o", ms=2.2, color="tab:green", alpha=0.8, label="THGM samples")
    ax0.plot(rgg, zgg, ".", ms=2.0, color="tab:orange", alpha=0.8, label="THGG samples")
    ax0.set_aspect("equal")
    ax0.set_xlabel("R (m)")
    ax0.set_ylabel("Z (m)")
    ax0.set_title("Boundary Overlay in R-Z")
    ax0.grid(alpha=0.25)
    ax0.legend(loc="best", fontsize=9)

    txt = (
        f"RR = {rr:.6f} m ({rr_source})\n"
        f"NTGMAX = {args.ntgmax}\n"
        f"Map residual mean/max = {np.mean(d_map):.3e} / {np.max(d_map):.3e} m\n"
        f"Nearest residual mean/max = {np.mean(d_near):.3e} / {np.max(d_near):.3e} m"
    )
    ax0.text(
        0.02,
        0.02,
        txt,
        transform=ax0.transAxes,
        fontsize=9,
        va="bottom",
        bbox=dict(facecolor="white", edgecolor="0.7", alpha=0.9),
    )

    ax1.plot(theta_raw, rho_raw, ".", ms=2.0, color="0.75", label="Raw (theta, rho) from RSU")
    ax1.plot(theta_u, rho_u, "-", lw=1.4, color="tab:blue", label="Sorted+merged nodes")
    ax1.plot(th_dense, rho_dense, "--", lw=1.8, color="tab:red", label="Interpolated rho(theta)")
    ax1.plot(thgm, rhgm, "o", ms=2.2, color="tab:green", alpha=0.8, label="THGM")
    ax1.plot(thgg, rhgg, ".", ms=2.0, color="tab:orange", alpha=0.8, label="THGG")
    ax1.set_xlim(0.0, 2.0 * np.pi)
    ax1.set_xlabel("theta (rad)")
    ax1.set_ylabel("rho (m)")
    ax1.set_title("Rebuild in theta-rho Space")
    ax1.grid(alpha=0.25)
    ax1.legend(loc="best", fontsize=9)

    out_path = Path(args.out)
    if not out_path.is_absolute():
        out_path = (eq_dir / out_path).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=170)
    print(f"Saved: {out_path}")
    print(
        "Residual stats [m]: "
        f"map mean={np.mean(d_map):.6e}, map max={np.max(d_map):.6e}, "
        f"nearest mean={np.mean(d_near):.6e}, nearest max={np.max(d_near):.6e}"
    )


if __name__ == "__main__":
    main()
