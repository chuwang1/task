#!/usr/bin/env python3

from pathlib import Path
import argparse
import pandas as pd
import numpy as np
import sys

BASE_DIR = Path(__file__).resolve().parent
EQ_DIR = BASE_DIR.parent / "eq"
sys.path.insert(0, str(EQ_DIR))

from solve_gs_from_gfile import read_gfile, rho_from_qpsi  # noqa: E402

RKEV = 1.602176634e-16


def read_tr_csv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, skiprows=1)
    df.columns = [str(c).strip() for c in df.columns]
    return df


def build_total_pressure_pa() -> tuple[np.ndarray, np.ndarray]:
    dn = read_tr_csv(BASE_DIR / "tr_data_074.csv")
    tt = read_tr_csv(BASE_DIR / "tr_data_076.csv")
    fn = read_tr_csv(BASE_DIR / "tr_data_075.csv")
    ft = read_tr_csv(BASE_DIR / "tr_data_030.csv")

    rho = dn["X"].to_numpy(float)
    p_total = (
        dn["nE"] * tt["TE"]
        + dn["nD"] * tt["TD"]
        + dn["nT"] * tt["TT"]
        + dn["nA"] * tt["TA"]
        + fn["NF"] * ft["TF"]
    ).to_numpy(float) * 1e20 * RKEV
    return rho, p_total


def build_total_pressure_from_snapshot(snapshot_csv: Path) -> tuple[np.ndarray, np.ndarray]:
    df = pd.read_csv(snapshot_csv)
    rho = df["rho"].to_numpy(float)
    p_total = (
        df["nE"] * df["TE"]
        + df["nD"] * df["TD"]
        + df["nT"] * df["TT"]
        + df["nA"] * df["TA"]
        + df["NF"] * df["TF"]
    ).to_numpy(float) * 1e20 * RKEV
    return rho, p_total


def build_ttrho_f() -> tuple[np.ndarray, np.ndarray]:
    ftr = read_tr_csv(BASE_DIR / "tr_data_168.csv")
    rho = ftr["X"].to_numpy(float)
    f_rho = -ftr["TTRHO"].to_numpy(float)
    return rho, f_rho


def build_ttrho_f_from_snapshot(snapshot_csv: Path) -> tuple[np.ndarray, np.ndarray]:
    df = pd.read_csv(snapshot_csv)
    rho = df["rho"].to_numpy(float)
    f_rho = -df["TTRHO"].to_numpy(float)
    return rho, f_rho


def main() -> None:
    parser = argparse.ArgumentParser(description="Build GS profile CSV for offline tests")
    parser.add_argument("--gfile", required=True, help="Reference gfile path")
    parser.add_argument(
        "--pressure-mode",
        choices=["gfile", "tr_total"],
        default="tr_total",
        help="Pressure source",
    )
    parser.add_argument(
        "--f-mode",
        choices=["gfile", "ttrho"],
        default="gfile",
        help="F source",
    )
    parser.add_argument("--snapshot-csv", help="Runtime snapshot CSV from TR")
    parser.add_argument("--output", required=True, help="Output CSV path")
    parser.add_argument("--debug", action="store_true", help="Print debug output")
    args = parser.parse_args()

    g = read_gfile(args.gfile)
    psi_wb = g.psimag + g.psi_n * (g.psibdy - g.psimag)

    if args.pressure_mode == "gfile":
        p_pa = np.asarray(g.pres, dtype=float)
    else:
        if args.snapshot_csv:
            rho_p, p_total = build_total_pressure_from_snapshot(Path(args.snapshot_csv))
        else:
            rho_p, p_total = build_total_pressure_pa()
        rho_g = rho_from_qpsi(g.qpsi, g.psimag, g.psibdy)
        p_pa = np.interp(rho_g, rho_p, p_total)

    if args.f_mode == "gfile":
        f_tesla_meter = np.asarray(g.fpol, dtype=float)
        ffprime = np.asarray(g.ffprime, dtype=float)
    else:
        if args.snapshot_csv:
            rho_f, f_rho = build_ttrho_f_from_snapshot(Path(args.snapshot_csv))
        else:
            rho_f, f_rho = build_ttrho_f()
        rho_g = rho_from_qpsi(g.qpsi, g.psimag, g.psibdy)
        f_tesla_meter = np.interp(rho_g, rho_f, f_rho)
        ffprime = f_tesla_meter * np.gradient(f_tesla_meter, psi_wb)

    out = pd.DataFrame(
        {
            "psi_n": g.psi_n,
            "psi_wb": psi_wb,
            "p_pa": p_pa,
            "F_tesla_meter": f_tesla_meter,
            "FFprime": ffprime,
        }
    )
    out_path = Path(args.output)
    out.to_csv(out_path, index=False)
    if args.debug:
        print(out_path)


if __name__ == "__main__":
    main()
