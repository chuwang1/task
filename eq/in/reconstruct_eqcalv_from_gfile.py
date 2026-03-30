#!/usr/bin/env python3
"""Reconstruct EQCALV-like geometry quantities from a GEQDSK gfile.

This script computes RSV, VPV, AVIR2, and AVRR2 on a user-selected psi_n
grid using contour integrals on psi(R,Z)=const surfaces extracted from psirz.

Notes
- Uses a factor=1 approximation for VPV/AVRR2 expressions here.
- R0 and B0 are taken as gfile `rcentr` and `bcentr`.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
from dataclasses import dataclass

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.path import Path as MplPath

try:
    from contourpy import contour_generator
except Exception:
    contour_generator = None

TWO_PI = 2.0 * np.pi
EPS = 1.0e-12


@dataclass
class GFile:
    header: str
    nw: int
    nh: int
    rdim: float
    zdim: float
    rcentr: float
    rleft: float
    zmid: float
    rmaxis: float
    zmaxis: float
    psimag: float
    psibdy: float
    bcentr: float
    currentA: float
    fpol: np.ndarray
    pres: np.ndarray
    ffprime: np.ndarray
    pprime: np.ndarray
    psirz: np.ndarray
    qpsi: np.ndarray
    nbbbs: int
    limitr: int
    rbbbs: np.ndarray
    zbbbs: np.ndarray
    R: np.ndarray
    Z: np.ndarray


def _split_nums(text: str) -> list[float]:
    text = re.sub(r"([^Ee])-", r"\1 -", text)
    out = []
    for token in re.split(r"\s+", text.strip()):
        if not token:
            continue
        try:
            out.append(float(token))
        except ValueError:
            continue
    return out


def read_gfile(path: str) -> GFile:
    with open(path, "r", encoding="utf-8") as fh:
        m = re.match(r"^\s*(.*)\s+(\d+)\s+(\d+)\s*$", fh.readline())
        if not m:
            raise ValueError(f"Invalid gfile header: {path}")
        header = m.group(1)
        nw = int(m.group(2))
        nh = int(m.group(3))

        def read5() -> list[float]:
            vals = _split_nums(fh.readline())
            if len(vals) < 5:
                raise ValueError("Failed to parse 5-value gfile header row.")
            return vals[:5]

        row2 = read5()
        row3 = read5()
        row4 = read5()
        _ = read5()
        body = fh.read().replace("\n", " ")

    nums = _split_nums(body)
    ia = 0

    fpol = np.asarray(nums[ia : ia + nw], dtype=float)
    ia += nw
    pres = np.asarray(nums[ia : ia + nw], dtype=float)
    ia += nw
    ffprime = np.asarray(nums[ia : ia + nw], dtype=float)
    ia += nw
    pprime = np.asarray(nums[ia : ia + nw], dtype=float)
    ia += nw

    psirz = np.asarray(nums[ia : ia + nw * nh], dtype=float).reshape(nh, nw).T
    ia += nw * nh
    qpsi = np.asarray(nums[ia : ia + nw], dtype=float)
    ia += nw

    nbbbs = int(nums[ia])
    ia += 1
    limitr = int(nums[ia])
    ia += 1

    rb = []
    zb = []
    for _ in range(nbbbs):
        rb.append(nums[ia])
        zb.append(nums[ia + 1])
        ia += 2

    rdim, zdim, rcentr, rleft, zmid = row2
    rmaxis, zmaxis, psimag, psibdy, bcentr = row3
    
    # GEQDSK psirz is typically poloidal flux divided by 2*pi (Wb/rad).
    # Convert to true Wb for consistency with EQCALV/EQIPQP formulations.
    psimag *= TWO_PI
    psibdy *= TWO_PI
    psirz *= TWO_PI

    currentA = row4[0]
    r = np.linspace(rleft, rleft + rdim, nw)
    z = np.linspace(zmid - 0.5 * zdim, zmid + 0.5 * zdim, nh)

    return GFile(
        header=header,
        nw=nw,
        nh=nh,
        rdim=rdim,
        zdim=zdim,
        rcentr=rcentr,
        rleft=rleft,
        zmid=zmid,
        rmaxis=rmaxis,
        zmaxis=zmaxis,
        psimag=psimag,
        psibdy=psibdy,
        bcentr=bcentr,
        currentA=currentA,
        fpol=fpol,
        pres=pres,
        ffprime=ffprime,
        pprime=pprime,
        psirz=psirz,
        qpsi=qpsi,
        nbbbs=nbbbs,
        limitr=limitr,
        rbbbs=np.asarray(rb, dtype=float),
        zbbbs=np.asarray(zb, dtype=float),
        R=r,
        Z=z,
    )


def _interp_bilinear(field: np.ndarray, r: np.ndarray, z: np.ndarray, rp: np.ndarray, zp: np.ndarray) -> np.ndarray:
    ir = np.searchsorted(r, rp, side="right") - 1
    iz = np.searchsorted(z, zp, side="right") - 1
    ir = np.clip(ir, 0, len(r) - 2)
    iz = np.clip(iz, 0, len(z) - 2)

    r0 = r[ir]
    r1 = r[ir + 1]
    z0 = z[iz]
    z1 = z[iz + 1]
    tr = (rp - r0) / np.maximum(r1 - r0, EPS)
    tz = (zp - z0) / np.maximum(z1 - z0, EPS)

    f00 = field[ir, iz]
    f10 = field[ir + 1, iz]
    f01 = field[ir, iz + 1]
    f11 = field[ir + 1, iz + 1]
    return (1.0 - tr) * (1.0 - tz) * f00 + tr * (1.0 - tz) * f10 + (1.0 - tr) * tz * f01 + tr * tz * f11


def _extract_surface(psirz: np.ndarray, r: np.ndarray, z: np.ndarray, level: float, axis_r: float, axis_z: float) -> np.ndarray | None:
    segs = []
    if contour_generator is not None:
        cg = contour_generator(x=r, y=z, z=psirz.T)
        segs = cg.lines(level)
    else:
        fig, ax = plt.subplots(1, 1)
        try:
            rr2d, zz2d = np.meshgrid(r, z, indexing="xy")
            cs = ax.contour(rr2d, zz2d, psirz.T, levels=[level])
            if cs.allsegs and cs.allsegs[0]:
                segs = cs.allsegs[0]
        finally:
            plt.close(fig)

    if not segs:
        return None

    best = None
    best_score = -1.0
    for seg in segs:
        seg = np.asarray(seg, dtype=float)
        if seg.ndim != 2 or seg.shape[0] < 8:
            continue
        if np.hypot(seg[0, 0] - seg[-1, 0], seg[0, 1] - seg[-1, 1]) > 1e-8:
            seg = np.vstack([seg, seg[0]])
        d = np.diff(seg, axis=0)
        plen = float(np.sum(np.hypot(d[:, 0], d[:, 1])))
        if plen <= 0.0:
            continue
        contains_axis = MplPath(seg).contains_point((axis_r, axis_z))
        score = plen + (1.0e9 if contains_axis else 0.0)
        if score > best_score:
            best_score = score
            best = seg
    return best


def _cumtrapz(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    out = np.zeros_like(x, dtype=float)
    if len(x) < 2:
        return out
    dx = np.diff(x)
    out[1:] = np.cumsum(0.5 * (y[:-1] + y[1:]) * dx)
    return out


def _fill_interior_nans(y: np.ndarray) -> np.ndarray:
    y = np.asarray(y, dtype=float).copy()
    n = y.size
    if n < 3:
        return y
    idx = np.arange(n, dtype=float)
    good = np.isfinite(y)
    if np.count_nonzero(good) < 2:
        return y
    interp = np.interp(idx, idx[good], y[good])
    first = int(idx[good][0])
    last = int(idx[good][-1])
    mask = (~good) & (idx >= first) & (idx <= last)
    y[mask] = interp[mask]
    return y


def compute_profiles(g: GFile, npsi: int) -> dict[str, np.ndarray]:
    npsi = max(3, int(npsi))
    psi_n = np.linspace(0.0, 1.0, npsi)
    psi_wb = g.psimag + psi_n * (g.psibdy - g.psimag)
    q_base_x = np.linspace(0.0, 1.0, g.qpsi.size)
    q = np.interp(psi_n, q_base_x, g.qpsi)

    # Toroidal flux from q(psi): Phi(psi) = integral q dpsi.
    phi = _cumtrapz(q, psi_wb)
    rsv = np.sqrt(np.abs(phi) / (np.pi * max(abs(g.bcentr), EPS)))

    dpsi_dr = np.gradient(g.psirz, g.R, axis=0, edge_order=2)
    dpsi_dz = np.gradient(g.psirz, g.Z, axis=1, edge_order=2)

    s0 = np.full(npsi, np.nan, dtype=float)
    s1 = np.full(npsi, np.nan, dtype=float)
    s2 = np.full(npsi, np.nan, dtype=float)
    contour_ok = np.zeros(npsi, dtype=bool)

    for i in range(1, npsi):
        seg = _extract_surface(g.psirz, g.R, g.Z, psi_wb[i], g.rmaxis, g.zmaxis)
        if seg is None or seg.shape[0] < 8:
            continue

        p0 = seg[:-1]
        p1 = seg[1:]
        rm = 0.5 * (p0[:, 0] + p1[:, 0])
        zm = 0.5 * (p0[:, 1] + p1[:, 1])
        dl = np.hypot(p1[:, 0] - p0[:, 0], p1[:, 1] - p0[:, 1])

        dpr = _interp_bilinear(dpsi_dr, g.R, g.Z, rm, zm)
        dpz = _interp_bilinear(dpsi_dz, g.R, g.Z, rm, zm)
        bp = np.sqrt(dpr * dpr + dpz * dpz) / np.maximum(TWO_PI * rm, EPS)

        valid = np.isfinite(dl) & np.isfinite(rm) & np.isfinite(bp) & (dl > 0.0) & (rm > EPS) & (bp > EPS)
        if np.count_nonzero(valid) < 8:
            continue

        dlv = dl[valid]
        rv = rm[valid]
        bpv = bp[valid]
        s0_i = np.sum(dlv / bpv)
        s1_i = np.sum(dlv / (bpv * rv * rv))
        s2_i = np.sum(bpv * dlv)
        if not (np.isfinite(s0_i) and np.isfinite(s1_i) and np.isfinite(s2_i) and s0_i > EPS):
            continue

        s0[i] = s0_i
        s1[i] = s1_i
        s2[i] = s2_i
        contour_ok[i] = True

    avir2 = (g.rcentr * g.rcentr) * s1 / s0

    vpv = np.full(npsi, np.nan, dtype=float)
    den_q = np.abs(q) > EPS
    valid_vpv = np.isfinite(rsv) & np.isfinite(s0) & den_q
    vpv[valid_vpv] = TWO_PI * rsv[valid_vpv] * g.bcentr * s0[valid_vpv] / q[valid_vpv]

    avrr2 = np.full(npsi, np.nan, dtype=float)
    den = (rsv * rsv) * (g.bcentr * g.bcentr) * s0
    valid_avrr2 = np.isfinite(s2) & np.isfinite(q) & np.isfinite(den) & (np.abs(den) > EPS)
    avrr2[valid_avrr2] = s2[valid_avrr2] * (q[valid_avrr2] ** 2) / den[valid_avrr2]

    # Axis handling: force RSV/VPV to 0, and extrapolate AVIR2/AVRR2 from 1,2.
    rsv[0] = 0.0
    vpv[0] = 0.0
    if npsi >= 3:
        if np.isfinite(avir2[1]) and np.isfinite(avir2[2]):
            avir2[0] = (4.0 * avir2[1] - avir2[2]) / 3.0
        if np.isfinite(avrr2[1]) and np.isfinite(avrr2[2]):
            avrr2[0] = (4.0 * avrr2[1] - avrr2[2]) / 3.0

    # Fill interior failed points from neighboring valid surfaces.
    s0 = _fill_interior_nans(s0)
    s1 = _fill_interior_nans(s1)
    s2 = _fill_interior_nans(s2)
    avir2 = _fill_interior_nans(avir2)
    vpv = _fill_interior_nans(vpv)
    avrr2 = _fill_interior_nans(avrr2)

    return {
        "index": np.arange(npsi, dtype=int),
        "psi_n": psi_n,
        "psi_wb": psi_wb,
        "q": q,
        "RSV": rsv,
        "S0": s0,
        "S1": s1,
        "S2": s2,
        "AVIR2": avir2,
        "VPV": vpv,
        "AVRR2": avrr2,
        "contour_ok": contour_ok,
    }


def save_csv(path: str, data: dict[str, np.ndarray]) -> None:
    cols = ["index", "psi_n", "psi_wb", "q", "RSV", "S0", "S1", "S2", "AVIR2", "VPV", "AVRR2", "contour_ok"]
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        n = len(data["index"])
        for i in range(n):
            w.writerow(
                [
                    int(data["index"][i]),
                    float(data["psi_n"][i]),
                    float(data["psi_wb"][i]),
                    float(data["q"][i]),
                    float(data["RSV"][i]) if np.isfinite(data["RSV"][i]) else np.nan,
                    float(data["S0"][i]) if np.isfinite(data["S0"][i]) else np.nan,
                    float(data["S1"][i]) if np.isfinite(data["S1"][i]) else np.nan,
                    float(data["S2"][i]) if np.isfinite(data["S2"][i]) else np.nan,
                    float(data["AVIR2"][i]) if np.isfinite(data["AVIR2"][i]) else np.nan,
                    float(data["VPV"][i]) if np.isfinite(data["VPV"][i]) else np.nan,
                    float(data["AVRR2"][i]) if np.isfinite(data["AVRR2"][i]) else np.nan,
                    bool(data["contour_ok"][i]),
                ]
            )


def save_plot(path: str, data: dict[str, np.ndarray]) -> None:
    x = data["psi_n"]
    ok = data["contour_ok"]
    fig, ax = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
    axs = ax.ravel()

    axs[0].plot(x, data["RSV"], lw=1.8)
    axs[0].set_title("RSV")
    axs[0].grid(True, alpha=0.3)

    axs[1].plot(x, data["VPV"], lw=1.8)
    axs[1].set_title("VPV")
    axs[1].grid(True, alpha=0.3)

    axs[2].plot(x, data["AVIR2"], lw=1.8, label="AVIR2")
    axs[2].plot(x, data["AVRR2"], lw=1.8, label="AVRR2")
    axs[2].set_title("AVIR2 and AVRR2")
    axs[2].legend(loc="best")
    axs[2].grid(True, alpha=0.3)

    axs[3].plot(x, data["q"], lw=1.8)
    axs[3].set_title("q")
    axs[3].grid(True, alpha=0.3)

    bad = np.where(~ok)[0]
    if bad.size > 0:
        axs[0].scatter(x[bad], data["RSV"][bad], s=14, color="tab:red", alpha=0.6, label="contour failed")
        axs[0].legend(loc="best")

    axs[2].set_xlabel("psi_n")
    axs[3].set_xlabel("psi_n")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Reconstruct EQCALV-like RSV/VPV/AVIR2/AVRR2 from a gfile.")
    parser.add_argument("--gfile", default="in/g260206.20000_teq_0114", help="Input GEQDSK gfile path.")
    parser.add_argument(
        "--out-prefix",
        default="in/reconstruct_eqcalv_from_gfile",
        help="Output prefix; writes <prefix>.csv and <prefix>.png",
    )
    parser.add_argument("--npsi", type=int, default=None, help="Number of psi_n points (default: gfile qpsi size).")
    args = parser.parse_args()

    gfile = os.path.abspath(args.gfile)
    if not os.path.exists(gfile):
        raise FileNotFoundError(f"Gfile not found: {gfile}")

    out_prefix = args.out_prefix
    out_dir = os.path.dirname(out_prefix) if os.path.dirname(out_prefix) else "."
    os.makedirs(out_dir, exist_ok=True)

    g = read_gfile(gfile)
    npsi = args.npsi if args.npsi is not None else int(g.qpsi.size)
    data = compute_profiles(g, npsi)

    csv_path = f"{out_prefix}.csv"
    png_path = f"{out_prefix}.png"
    save_csv(csv_path, data)
    save_plot(png_path, data)

    print(f"Gfile: {gfile}")
    print(f"Output CSV: {os.path.abspath(csv_path)}")
    print(f"Output PNG: {os.path.abspath(png_path)}")


if __name__ == "__main__":
    main()
