#!/usr/bin/env python3
"""
Export key 1D profiles from a G-file to CSV and plot figures.

Usage example:
  python export_gfile_profiles.py \
      --gfile in/g260206.20000_teq_0114 \
      --out-dir in
"""

import argparse
import os
import re
from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

TWO_PI = 2.0 * np.pi


@dataclass
class GFileData:
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
    qpsi: np.ndarray


def _split_numeric_line(line: str) -> list[float]:
    line = re.sub(r"([^Ee])\-", r"\1 -", line.strip())
    fields = re.split(r"\s+", line)
    vals = []
    for v in fields:
        if not v:
            continue
        try:
            vals.append(float(v))
        except ValueError:
            pass
    return vals


def _read_row_of_5(fh) -> list[float]:
    vals = _split_numeric_line(fh.readline())
    if len(vals) < 5:
        raise ValueError("Failed to parse a required 5-number G-file header row.")
    return vals[:5]


def read_gfile(path: str) -> GFileData:
    with open(path, "r", encoding="utf-8") as fh:
        m = re.match(r"^\s*(.*)\s+(\d+)\s+(\d+)\s*$", fh.readline())
        if not m:
            raise ValueError(f"Invalid G-file header in: {path}")
        nw = int(m.group(2))
        nh = int(m.group(3))

        row2 = _read_row_of_5(fh)
        row3 = _read_row_of_5(fh)
        row4 = _read_row_of_5(fh)
        fh.readline()  # row5, not used by this exporter

        body = fh.read()

    body = body.replace("\n", " ")
    body = re.sub(r"(\d)\-", r"\1 -", body)
    numbers = []
    for token in re.split(r"\s+", body):
        if not token:
            continue
        try:
            numbers.append(float(token))
        except ValueError:
            pass

    ia = 0
    fpol = np.asarray(numbers[ia : ia + nw], dtype=float)
    ia += nw
    pres = np.asarray(numbers[ia : ia + nw], dtype=float)
    ia += nw
    ffprime = np.asarray(numbers[ia : ia + nw], dtype=float)
    ia += nw
    pprime = np.asarray(numbers[ia : ia + nw], dtype=float)
    ia += nw

    # Skip PSIRZ (nw * nh)
    ia += nw * nh
    qpsi = np.asarray(numbers[ia : ia + nw], dtype=float)

    return GFileData(
        nw=nw,
        nh=nh,
        rdim=row2[0],
        zdim=row2[1],
        rcentr=row2[2],
        rleft=row2[3],
        zmid=row2[4],
        rmaxis=row3[0],
        zmaxis=row3[1],
        psimag=row3[2],
        psibdy=row3[3],
        bcentr=row3[4],
        currentA=row4[0],
        fpol=fpol,
        pres=pres,
        ffprime=ffprime,
        pprime=pprime,
        qpsi=qpsi,
    )


def export_profiles(
    g: GFileData, gfile_path: str, out_dir: str, prefix: str
) -> tuple[str, str, str]:
    os.makedirs(out_dir, exist_ok=True)

    psi_n = np.linspace(0.0, 1.0, g.nw)
    psi_wb = g.psimag + psi_n * (g.psibdy - g.psimag)

    profiles = pd.DataFrame(
        {
            "i": np.arange(g.nw, dtype=int),
            "psi_n": psi_n,
            "psi_wb": psi_wb,
            "p_pa": g.pres,
            "pprime_pa_per_wb": g.pprime,
            "F_tesla_meter": g.fpol,
            "F_2pi": TWO_PI * g.fpol,
            "FFprime": g.ffprime,
            "q": g.qpsi,
        }
    )
    profiles_csv = os.path.join(out_dir, f"{prefix}_profiles_from_gfile.csv")
    profiles.to_csv(profiles_csv, index=False)

    summary = pd.DataFrame(
        [
            ("gfile_path", gfile_path, "", "Source gfile path"),
            ("nw", g.nw, "", "R grid count"),
            ("nh", g.nh, "", "Z grid count"),
            ("currentA", g.currentA, "A", "Total plasma current"),
            ("rmaxis", g.rmaxis, "m", "Magnetic-axis R"),
            ("zmaxis", g.zmaxis, "m", "Magnetic-axis Z"),
            ("psimag", g.psimag, "Wb/rad", "Poloidal flux at axis"),
            ("psibdy", g.psibdy, "Wb/rad", "Poloidal flux at LCFS"),
            ("bcentr", g.bcentr, "T", "Vacuum toroidal field at R=rcentr"),
            ("rcentr", g.rcentr, "m", "Reference major radius for bcentr"),
            ("rdim", g.rdim, "m", "R-domain size"),
            ("zdim", g.zdim, "m", "Z-domain size"),
            ("rleft", g.rleft, "m", "R minimum in rectangular mesh"),
            ("zmid", g.zmid, "m", "Z center in rectangular mesh"),
        ],
        columns=["key", "value", "unit", "description"],
    )
    summary_csv = os.path.join(out_dir, f"{prefix}_summary_from_gfile.csv")
    summary.to_csv(summary_csv, index=False)

    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    ax = axes.ravel()

    ax[0].plot(psi_n, g.pres, lw=2)
    ax[0].set_title("p(psi)")
    ax[0].set_xlabel("psi_n")
    ax[0].set_ylabel("p [Pa]")
    ax[0].grid(True, alpha=0.3)

    ax[1].plot(psi_n, g.pprime, lw=2)
    ax[1].set_title("p'(psi)")
    ax[1].set_xlabel("psi_n")
    ax[1].set_ylabel("dp/dpsi")
    ax[1].grid(True, alpha=0.3)

    ax[2].plot(psi_n, g.fpol, lw=2)
    ax[2].set_title("F(psi) = R*Bt")
    ax[2].set_xlabel("psi_n")
    ax[2].set_ylabel("F [T*m]")
    ax[2].grid(True, alpha=0.3)

    ax[3].plot(psi_n, g.ffprime, lw=2)
    ax[3].set_title("FF'(psi)")
    ax[3].set_xlabel("psi_n")
    ax[3].set_ylabel("FF'")
    ax[3].grid(True, alpha=0.3)

    ax[4].plot(psi_n, g.qpsi, lw=2)
    ax[4].set_title("q(psi)")
    ax[4].set_xlabel("psi_n")
    ax[4].set_ylabel("q")
    ax[4].grid(True, alpha=0.3)

    ax[5].axis("off")
    ax[5].text(
        0.03,
        0.95,
        "\n".join(
            [
                f"Ip = {g.currentA:.6e} A",
                f"R_axis = {g.rmaxis:.6f} m",
                f"Z_axis = {g.zmaxis:.6f} m",
                f"psi_axis = {g.psimag:.6e}",
                f"psi_lcfs = {g.psibdy:.6e}",
                f"nw, nh = {g.nw}, {g.nh}",
            ]
        ),
        va="top",
        ha="left",
        fontsize=10,
        family="monospace",
    )

    fig.suptitle(f"G-file Profiles: {prefix}")
    fig.tight_layout()
    fig_path = os.path.join(out_dir, f"{prefix}_profiles_from_gfile.png")
    fig.savefig(fig_path, dpi=160, bbox_inches="tight")
    plt.close(fig)

    return profiles_csv, summary_csv, fig_path


def main():
    parser = argparse.ArgumentParser(
        description="Export pressure/current-related profiles from a G-file."
    )
    parser.add_argument(
        "--gfile",
        default="in/g260206.20000_teq_0114",
        help="Input gfile path",
    )
    parser.add_argument(
        "--out-dir",
        default="in",
        help="Output directory for CSV/PNG",
    )
    parser.add_argument(
        "--prefix",
        default=None,
        help="Output filename prefix; default is basename(gfile)",
    )
    args = parser.parse_args()

    gfile_abs = os.path.abspath(args.gfile)
    if not os.path.exists(gfile_abs):
        raise FileNotFoundError(f"G-file not found: {gfile_abs}")

    prefix = args.prefix if args.prefix else os.path.basename(gfile_abs)
    g = read_gfile(gfile_abs)
    profiles_csv, summary_csv, fig_path = export_profiles(
        g, gfile_abs, args.out_dir, prefix
    )

    print(f"G-file: {gfile_abs}")
    print(f"Saved CSV: {os.path.abspath(profiles_csv)}")
    print(f"Saved CSV: {os.path.abspath(summary_csv)}")
    print(f"Saved PNG: {os.path.abspath(fig_path)}")


if __name__ == "__main__":
    main()
