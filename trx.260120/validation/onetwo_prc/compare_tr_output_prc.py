#!/usr/bin/env python3
"""Compare TR ONETWO PRC output against PyMak recomputation from TR CSV."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import re
import sys
import types
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


REPO_DIR = Path(__file__).resolve().parents[2]
PYMAK_DIR = REPO_DIR / "PyMak"


def load_pymak_onetwo():
    """Load PyMak's ONETWO function without importing PyMak package globals."""

    pkg = types.ModuleType("PyMak")
    pkg.__path__ = [str(PYMAK_DIR)]
    sys.modules.setdefault("PyMak", pkg)

    stub = types.ModuleType("PyMak.gfile_reader")
    stub.build_cytran_input_arrays = lambda *args, **kwargs: None
    stub.map_psi_n_to_rho_tor = lambda *args, **kwargs: None
    sys.modules["PyMak.gfile_reader"] = stub

    spec = importlib.util.spec_from_file_location(
        "PyMak.cyclotron_radiation", PYMAK_DIR / "cyclotron_radiation.py"
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules["PyMak.cyclotron_radiation"] = module
    spec.loader.exec_module(module)
    return module.calculate_cyclotron_radiation_onetwo


def read_gfile_runtime_geometry(path: Path):
    """Return TR's EQDSK-derived RR, RA, and BB convention."""

    with path.open("r", encoding="utf-8") as stream:
        header = stream.readline()
        match = re.match(r"^\s*(.*)\s+(\d+)\s+(\d+)\s*$", header)
        if not match:
            raise ValueError(f"Invalid gfile header: {path}")
        nw = int(match.group(2))
        nh = int(match.group(3))

        def read5():
            line = re.sub(r"([^Ee])-", r"\1 -", stream.readline())
            values = [float(x) for x in re.split(r"\s+", line.strip()) if x]
            if len(values) < 5:
                raise ValueError(f"Failed to read 5-value gfile row from {path}")
            return values[:5]

        r2 = read5()
        r3 = read5()
        _ = read5()
        _ = read5()
        body = stream.read().split("&", 1)[0]

    nums = [float(x) for x in re.split(r"\s+", re.sub(r"(\d)-", r"\1 -", body)) if x]
    cursor = 5 * nw + nw * nh
    nbbbs = int(nums[cursor])
    cursor += 2
    rbbbs = np.asarray([nums[cursor + 2 * i] for i in range(nbbbs)], dtype=float)

    rcentr = float(r2[2])
    bcentr = float(r3[4])
    rr = 0.5 * (float(np.max(rbbbs)) + float(np.min(rbbbs)))
    ra = 0.5 * (float(np.max(rbbbs)) - float(np.min(rbbbs)))
    bb = bcentr * rcentr / rr
    return types.SimpleNamespace(bcentr=bb, rmaxis=rr, rbbbs=rbbbs), rr, ra, bb


def read_tr_csv(path: Path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(row for row in stream if not row.startswith(" Title")))


def csv_title(path: Path) -> str:
    with path.open() as stream:
        first = stream.readline().strip()
    return first.removeprefix("Title:").strip()


def csv_number(path: Path) -> int:
    return int(path.stem.rsplit("_", 1)[1])


def find_direct_profile_group(run_dir: Path):
    """Find the TRGRR1/TRGRR2 group: n(NS), T(NS), then PRC."""

    files = sorted(run_dir.glob("tr_data_*.csv"), key=csv_number)
    numbered = [(csv_number(path), path, csv_title(path)) for path in files]
    density = [(idx, path) for idx, path, title in numbered if "@n(NS)" in title]
    temperature = [(idx, path) for idx, path, title in numbered if "@T(NS)" in title]
    radiation = [(idx, path) for idx, path, title in numbered if "PRSUM,PRB,PRC,PRL" in title]

    for n_idx, n_path in density:
        t_matches = [(idx, path) for idx, path in temperature if n_idx < idx <= n_idx + 4]
        p_matches = [(idx, path) for idx, path in radiation if n_idx < idx <= n_idx + 10]
        if t_matches and p_matches:
            return n_path, t_matches[0][1], p_matches[0][1]

    raise FileNotFoundError(f"Could not locate a direct n/T/PRC CSV group in {run_dir}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--gfile", type=Path, required=True)
    parser.add_argument("--out-prefix", type=Path, required=True)
    parser.add_argument("--wall-reflection", type=float, default=0.8)
    args = parser.parse_args()

    calculate_onetwo = load_pymak_onetwo()
    geometry, rr, ra, bb = read_gfile_runtime_geometry(args.gfile)
    density_csv, temperature_csv, radiation_csv = find_direct_profile_group(args.run_dir)

    density_rows = read_tr_csv(density_csv)
    temperature_rows = read_tr_csv(temperature_csv)
    radiation_rows = read_tr_csv(radiation_csv)

    rho = np.asarray([float(row["X"]) for row in density_rows], dtype=float)
    ne = np.asarray([float(row["nE"]) for row in density_rows], dtype=float)
    te = np.asarray([float(row["TE"]) for row in temperature_rows], dtype=float)
    tr_prc = np.asarray([float(row["PRC"]) for row in radiation_rows], dtype=float)

    result = calculate_onetwo(
        geometry,
        ne,
        te,
        bt_T=bb,
        wall_reflection=args.wall_reflection,
        minor_radius_m=ra,
        major_radius_m=rr,
    )
    pymak_prc = np.asarray(result["power_density_total_MW_m3"], dtype=float)
    abs_err = np.abs(pymak_prc - tr_prc)
    rel_err = abs_err / np.maximum(np.abs(tr_prc), 1.0e-30)
    mask = np.abs(tr_prc) > 1.0e-8

    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    out_csv = args.out_prefix.with_suffix(".csv")
    with out_csv.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "rho",
                "ne_1e20_m3",
                "te_kev",
                "tr_prc_mw_m3",
                "pymak_prc_mw_m3",
                "abs_err_mw_m3",
                "rel_err",
            ]
        )
        for row in zip(rho, ne, te, tr_prc, pymak_prc, abs_err, rel_err):
            writer.writerow([f"{value:.17e}" for value in row])

    out_png = args.out_prefix.with_suffix(".png")
    fig, axes = plt.subplots(2, 1, figsize=(7.5, 7.0), sharex=True)
    axes[0].plot(rho, tr_prc, label="TR output PRC", lw=2.0)
    axes[0].plot(rho, pymak_prc, "--", label="PyMak recompute", lw=1.6)
    axes[0].set_ylabel("PRC [MW/m^3]")
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    axes[1].plot(rho, rel_err, color="tab:red", lw=1.5)
    axes[1].set_xlabel("rho")
    axes[1].set_ylabel("relative error")
    axes[1].grid(True, alpha=0.3)
    fig.suptitle(
        "ONETWO PRC: PyMak recompute vs TR output, 3 steps\n"
        f"RR={rr:.8f} m, RA={ra:.8f} m, |BB|={abs(bb):.8f} T, "
        f"REFRAD={args.wall_reflection:g}"
    )
    fig.tight_layout()
    fig.savefig(out_png, dpi=180)

    print(f"density_csv={density_csv}")
    print(f"temperature_csv={temperature_csv}")
    print(f"radiation_csv={radiation_csv}")
    print(f"geometry_rr_m={rr:.10e}")
    print(f"geometry_ra_m={ra:.10e}")
    print(f"geometry_bb_signed_T={bb:.10e}")
    print(f"geometry_abs_bb_T={abs(bb):.10e}")
    print(f"points={rho.size}")
    print(f"max_abs_err_mw_m3={float(abs_err.max()):.10e}")
    print(f"max_rel_err_masked={float(rel_err[mask].max()):.10e}")
    print(f"mean_rel_err_masked={float(rel_err[mask].mean()):.10e}")
    print(f"core_tr={float(tr_prc[0]):.10e}")
    print(f"core_pymak={float(pymak_prc[0]):.10e}")
    print(f"core_rel_err={float(rel_err[0]):.10e}")
    print(f"out_csv={out_csv}")
    print(f"out_png={out_png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
