#!/usr/bin/env python3
import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from contourpy import contour_generator
from matplotlib.path import Path as MplPath

TWO_PI = 2.0 * np.pi


class GFile:
    pass


def read_gfile(path: str) -> GFile:
    g = GFile()
    with open(path, "r", encoding="utf-8") as f:
        head = re.match(r"^\s*(.*)\s+(\d+)\s+(\d+)\s*$", f.readline())
        if not head:
            raise ValueError(f"Invalid gfile header: {path}")
        g.header = head.group(1)
        g.nw = int(head.group(2))
        g.nh = int(head.group(3))

        def read5() -> list[float]:
            line = re.sub(r"([^Ee])\-", r"\1 -", f.readline())
            vals = [float(x) for x in re.split(r"\s+", line.strip()) if x]
            if len(vals) < 5:
                raise ValueError("Failed to parse 5-value row in gfile.")
            return vals[:5]

        r2 = read5()
        r3 = read5()
        r4 = read5()
        _ = read5()

        g.rdim, g.zdim, g.rcentr, g.rleft, g.zmid = r2
        g.rmaxis, g.zmaxis, g.psimag, g.psibdy, g.bcentr = r3
        g.currentA = r4[0]

        body = f.read().replace("\n", " ")

    body = re.sub(r"(\d)\-", r"\1 -", body)
    nums: list[float] = []
    for tok in re.split(r"\s+", body):
        if not tok:
            continue
        try:
            nums.append(float(tok))
        except ValueError:
            pass

    ia = 0
    nw, nh = g.nw, g.nh
    g.fpol = np.asarray(nums[ia : ia + nw], dtype=float)
    ia += nw
    g.pres = np.asarray(nums[ia : ia + nw], dtype=float)
    ia += nw
    g.ffprime = np.asarray(nums[ia : ia + nw], dtype=float)
    ia += nw
    g.pprime = np.asarray(nums[ia : ia + nw], dtype=float)
    ia += nw

    g.psirz = np.asarray(nums[ia : ia + nw * nh], dtype=float).reshape(nw, nh)
    ia += nw * nh
    g.qpsi = np.asarray(nums[ia : ia + nw], dtype=float)
    ia += nw

    g.nbbbs = int(nums[ia])
    ia += 1
    g.limitr = int(nums[ia])
    ia += 1
    rb, zb = [], []
    for _ in range(g.nbbbs):
        rb.append(nums[ia])
        zb.append(nums[ia + 1])
        ia += 2
    g.rbbbs = np.asarray(rb, dtype=float)
    g.zbbbs = np.asarray(zb, dtype=float)

    g.R = np.linspace(g.rleft, g.rleft + g.rdim, g.nw)
    g.Z = np.linspace(g.zmid - 0.5 * g.zdim, g.zmid + 0.5 * g.zdim, g.nh)
    g.psi_n = np.linspace(0.0, 1.0, g.nw)
    return g


def _float_tokens(line: str) -> list[float]:
    vals = []
    for tok in re.split(r"\s+", line.strip()):
        if not tok:
            continue
        try:
            vals.append(float(tok))
        except ValueError:
            pass
    return vals


def parse_pltout_polylines(path: Path, rmin: float, rmax: float, zmin: float, zmax: float) -> list[np.ndarray]:
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    out: list[np.ndarray] = []

    i = 0
    while i < len(lines):
        m = re.match(r"^\s*([+-]?\d+)\s+([+-]?\d+)\s*$", lines[i])
        if m:
            npts = int(m.group(1))
            if 20 <= npts <= 5000:
                vals: list[float] = []
                j = i + 1
                while j < len(lines) and len(vals) < 2 * npts:
                    vals.extend(_float_tokens(lines[j]))
                    j += 1

                if len(vals) >= 2 * npts:
                    arr = np.asarray(vals[: 2 * npts], dtype=float).reshape(npts, 2)
                    good = np.isfinite(arr[:, 0]) & np.isfinite(arr[:, 1])
                    arr = arr[good]
                    if len(arr) >= 10:
                        inside = (
                            (arr[:, 0] >= rmin)
                            & (arr[:, 0] <= rmax)
                            & (arr[:, 1] >= zmin)
                            & (arr[:, 1] <= zmax)
                        )
                        if np.count_nonzero(inside) >= max(8, int(0.6 * len(arr))):
                            out.append(arr)
                    i = j
                    continue
        i += 1

    return out


def extract_gfile_contours(g, psi_n_levels: np.ndarray) -> list[np.ndarray]:
    psi_total = TWO_PI * (np.asarray(g.psirz, dtype=float) - float(g.psimag))
    dpsi_total = TWO_PI * (float(g.psibdy) - float(g.psimag))
    cg = contour_generator(x=np.asarray(g.R, float), y=np.asarray(g.Z, float), z=psi_total)

    out: list[np.ndarray] = []
    for psin in psi_n_levels:
        level = float(psin) * dpsi_total
        segs = cg.lines(level)
        if not segs:
            continue

        best = None
        best_score = -1.0
        for seg in segs:
            seg = np.asarray(seg, dtype=float)
            if seg.ndim != 2 or seg.shape[0] < 8:
                continue
            d = np.diff(seg, axis=0)
            plen = float(np.sum(np.hypot(d[:, 0], d[:, 1])))
            if plen <= 0.0:
                continue
            contains = MplPath(seg).contains_point((g.rmaxis, g.zmaxis))
            score = plen + (1.0e9 if contains else 0.0)
            if score > best_score:
                best = seg
                best_score = score

        if best is not None:
            out.append(best)

    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="Compare pltout polylines vs gfile contours")
    ap.add_argument("--gfile", default="/Users/dengxiaoya/CFEDRSW/OMFIT_out/CFEDR_v2_0114/EFIT/FILES/g0.30000")
    ap.add_argument("--pltout", default="/Users/dengxiaoya/CFEDRSW/OMFIT_out/CFEDR_v2_0114/EFIT/FILES/pltout.out")
    ap.add_argument("--out", default="in/pltout_vs_gfile_contours_overlay.png")
    args = ap.parse_args()

    g = read_gfile(args.gfile)
    rmin, rmax = float(np.min(g.R)), float(np.max(g.R))
    zmin, zmax = float(np.min(g.Z)), float(np.max(g.Z))

    polys = parse_pltout_polylines(Path(args.pltout), rmin, rmax, zmin, zmax)
    psi_n_levels = np.array([0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.97, 1.0], dtype=float)
    contours = extract_gfile_contours(g, psi_n_levels)

    fig, ax = plt.subplots(1, 1, figsize=(8.6, 8.2), constrained_layout=True)

    for p in polys:
        ax.plot(p[:, 0], p[:, 1], color="0.7", lw=0.7, alpha=0.75)

    cmap = plt.cm.viridis
    n = max(len(contours), 1)
    for i, c in enumerate(contours):
        ax.plot(c[:, 0], c[:, 1], color=cmap(i / (n - 1 if n > 1 else 1)), lw=1.4)

    if hasattr(g, "rbbbs") and hasattr(g, "zbbbs") and len(g.rbbbs) > 5:
        ax.plot(g.rbbbs, g.zbbbs, "k--", lw=1.2, label="LCFS from gfile")

    ax.plot([g.rmaxis], [g.zmaxis], "r*", ms=9, label="Magnetic axis")
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(rmin, rmax)
    ax.set_ylim(zmin, zmax)
    ax.set_xlabel("R [m]")
    ax.set_ylabel("Z [m]")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best", fontsize=8)
    ax.set_title(f"pltout vs gfile contours\npolylines={len(polys)}, gfile-levels={len(contours)}")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=180)
    plt.close(fig)

    print(f"saved: {out}")
    print(f"parsed pltout polylines: {len(polys)}")
    print(f"extracted gfile contours: {len(contours)}")


if __name__ == "__main__":
    main()
