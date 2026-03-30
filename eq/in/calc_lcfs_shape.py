#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np


def lcfs_eq(rr: float, ra: float, rkap: float, rdlt: float, n: int) -> tuple[np.ndarray, np.ndarray]:
    t = np.linspace(0.0, 2.0 * math.pi, n, endpoint=False)
    r = rr + ra * np.cos(t + rdlt * np.sin(t))
    z = ra * rkap * np.sin(t)
    return r, z


def _qx_ex2(rmaj: float, rpla: float, elip: float, trig: float) -> tuple[float, float, float]:
    ra2 = rmaj * rmaj + rpla * rpla
    rx = rmaj - trig * rpla
    zx = elip * rpla
    rx2 = rx * rx
    qx = 1.0 - 2.0 * ra2 * (ra2 - rx2) / (((rx2 - ra2) ** 2) + (2.0 * rmaj * rpla) ** 2)
    ex2 = 4.0 * ra2 * (zx * zx) / (((rx2 - ra2) ** 2) + (2.0 * rmaj * rpla) ** 2)
    return qx, ex2, zx


def _rpair_at_z(ra2: float, rmaj: float, rpla: float, qx: float, ex2: float, zabs: float) -> tuple[float, float]:
    z3 = (zabs * zabs) / max(ex2, 1e-14)
    y = 2.0 * math.sqrt(max(((1.0 - qx) * z3) ** 2 - ra2 * z3 + (rmaj * rpla) ** 2, 0.0))
    rp = math.sqrt(max(ra2 - 2.0 * (1.0 - qx) * z3 + y, 0.0))
    rm = math.sqrt(max(ra2 - 2.0 * (1.0 - qx) * z3 - y, 0.0))
    return rp, rm


def lcfs_equ_like(
    rr: float,
    ra: float,
    rkap: float,
    rdlt: float,
    zpla: float,
    yh: float,
    yd: float,
    elipup: float,
    trigup: float,
    n: int,
) -> tuple[np.ndarray, np.ndarray]:
    # map common params to equ naming
    rmaj = rr
    rpla = ra
    elip = rkap
    trig = rdlt
    if elipup <= 0.0:
        elipup = elip
    # if user leaves trigup NaN-like impossible by argparse; keep explicit value

    qx, ex2, zx = _qx_ex2(rmaj, rpla, elip, trig)
    qxup, ex2up, zxup = _qx_ex2(rmaj, rpla, elipup, trigup)
    ra2 = rmaj * rmaj + rpla * rpla

    # Build marker points following equ/eqinit.f msfx=9 asymmetric branch.
    pts: list[tuple[float, float]] = []
    pts.append((rmaj - rpla, zpla))
    pts.append((rmaj + rpla, zpla))

    # lower edge (yd)
    zd = yd * zx
    r2, r3 = _rpair_at_z(ra2, rmaj, rpla, qx, ex2, abs(zd))
    pts.extend([(r2, -abs(zd) + zpla), (r3, -abs(zd) + zpla)])

    # lower half-height (yh)
    zh = yh * zx
    r4, r5 = _rpair_at_z(ra2, rmaj, rpla, qx, ex2, abs(zh))
    pts.extend([(r4, -abs(zh) + zpla), (r5, -abs(zh) + zpla)])

    # upper edge (yd)
    zdup = yd * zxup
    r6, r7 = _rpair_at_z(ra2, rmaj, rpla, qxup, ex2up, abs(zdup))
    pts.extend([(r6, abs(zdup) + zpla), (r7, abs(zdup) + zpla)])

    # upper half-height (yh)
    zhup = yh * zxup
    r8, r9 = _rpair_at_z(ra2, rmaj, rpla, qxup, ex2up, abs(zhup))
    pts.extend([(r8, abs(zhup) + zpla), (r9, abs(zhup) + zpla)])

    pa = np.asarray(pts, dtype=float)
    cx, cy = rmaj, zpla
    ang = np.arctan2(pa[:, 1] - cy, pa[:, 0] - cx)
    ang = np.where(ang < 0.0, ang + 2.0 * math.pi, ang)
    idx = np.argsort(ang)
    angs = ang[idx]
    rs = pa[idx, 0]
    zs = pa[idx, 1]

    # periodic interpolation on angle
    ang_ext = np.concatenate([angs[-1:] - 2.0 * math.pi, angs, angs[:1] + 2.0 * math.pi])
    r_ext = np.concatenate([rs[-1:], rs, rs[:1]])
    z_ext = np.concatenate([zs[-1:], zs, zs[:1]])

    t = np.linspace(0.0, 2.0 * math.pi, n, endpoint=False)
    r = np.interp(t, ang_ext, r_ext)
    z = np.interp(t, ang_ext, z_ext)
    return r, z


def save_csv(path: Path, r: np.ndarray, z: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["R", "Z"])
        for ri, zi in zip(r, z):
            w.writerow([f"{ri:.12e}", f"{zi:.12e}"])


def save_plot(path: Path, req: np.ndarray, zeq: np.ndarray, requ: np.ndarray, zequ: np.ndarray) -> None:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6.6, 6.2))
    ax.plot(req, zeq, "b-", lw=1.8, label="eq parameterization")
    ax.plot(requ, zequ, "r--", lw=1.8, label="equ-like upper/lower")
    ax.set_xlabel("R [m]")
    ax.set_ylabel("Z [m]")
    ax.set_title("LCFS: eq vs equ-like")
    ax.axis("equal")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160)
    plt.close(fig)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Compute and compare LCFS curves for eq and equ-like parameterizations.")
    p.add_argument("--RR", type=float, required=True, help="Major radius [m] (common)")
    p.add_argument("--RA", type=float, required=True, help="Minor radius [m] (common)")
    p.add_argument("--RKAP", type=float, required=True, help="Elongation (common / lower for equ-like)")
    p.add_argument("--RDLT", type=float, required=True, help="Triangularity (common / lower for equ-like)")
    p.add_argument("--zpla", type=float, default=0.0, help="equ-like axis Z shift")
    p.add_argument("--yh", type=float, default=0.5, help="equ-like half-height marker fraction")
    p.add_argument("--yd", type=float, default=0.995, help="equ-like edge-height marker fraction")
    p.add_argument("--elipup", type=float, default=-1.0, help="equ-like upper elongation; <=0 means use RKAP")
    p.add_argument("--trigup", type=float, default=0.0, help="equ-like upper triangularity")
    p.add_argument("-n", "--num-points", type=int, default=720, help="Number of sampled points")
    return p


def main() -> None:
    args = build_parser().parse_args()

    r_eq, z_eq = lcfs_eq(args.RR, args.RA, args.RKAP, args.RDLT, args.num_points)
    r_equ, z_equ = lcfs_equ_like(
        args.RR,
        args.RA,
        args.RKAP,
        args.RDLT,
        args.zpla,
        args.yh,
        args.yd,
        args.elipup,
        args.trigup,
        args.num_points,
    )

    # fixed output names (requested)
    out_eq = Path("in/lcfs_eq.csv")
    out_equ = Path("in/lcfs_equ.csv")
    out_png = Path("in/lcfs_eq_vs_equ.png")

    save_csv(out_eq, r_eq, z_eq)
    save_csv(out_equ, r_equ, z_equ)
    save_plot(out_png, r_eq, z_eq, r_equ, z_equ)

    print(f"Saved EQ LCFS points : {out_eq} ({len(r_eq)} points)")
    print(f"Saved EQU LCFS points: {out_equ} ({len(r_equ)} points)")
    print(f"Saved comparison plot: {out_png}")


if __name__ == "__main__":
    main()
