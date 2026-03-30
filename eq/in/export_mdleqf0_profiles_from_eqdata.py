#!/usr/bin/env python3
"""
Export EQ (MDLEQF=0-style) profile inputs from eqdata binary.

Outputs:
1) profile CSV on psi grid: P, J_tor (reconstructed), T, Vphi
2) scalar CSV: Ip
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Export P + J_tor + T + Vphi + Ip from eqdata."
    )
    parser.add_argument(
        "--eqdata",
        default="eqdata0114",
        help="Path to eqdata binary file (default: eqdata0114).",
    )
    parser.add_argument(
        "--out-dir",
        default="in",
        help="Output directory for CSV files (default: in).",
    )
    parser.add_argument(
        "--prefix",
        default=None,
        help="Output filename prefix (default: basename(eqdata)).",
    )
    args = parser.parse_args()

    eqdata_path = Path(args.eqdata)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = args.prefix if args.prefix else eqdata_path.name

    root_dir = Path(__file__).resolve().parents[1]
    if str(root_dir) not in sys.path:
        sys.path.insert(0, str(root_dir))
    from read_eqdata_to_csv import read_eqdata  # pylint: disable=import-error

    data = read_eqdata(str(eqdata_path))

    # Profiles on psi-poloidal grid from EQCALP output.
    psip = np.asarray(data["PSIPS"], dtype=float)
    p_pa = np.asarray(data["PPPS"], dtype=float)
    tt = np.asarray(data["TTPS"], dtype=float)
    t_kev = np.asarray(data["TEPS"], dtype=float)
    omega = np.asarray(data["OMPS"], dtype=float)

    psipa = float(data["PSIPA"])
    rr = float(data["RR"])
    raxis = float(data["RAXIS"])
    rip_ma = float(data["RIP"])
    mu0 = 4.0 * np.pi * 1.0e-7

    # Normalized poloidal flux psi_N in [0,1].
    if abs(psipa) > 0.0:
        psi_n = np.clip(psip / psipa, 0.0, 1.0)
    else:
        psi_n = np.zeros_like(psip)

    # Derivatives on the native psi grid.
    dp_dpsi = np.gradient(p_pa, psip, edge_order=2)
    dtt_dpsi = np.gradient(tt, psip, edge_order=2)

    # Same definition as EQCALQ AVEJTR:
    # AVEJTR = -RR*DPPS - TTS*DTTS/(4*pi^2*mu0*RR)
    jtor_ave = -rr * dp_dpsi - tt * dtt_dpsi / (4.0 * np.pi**2 * mu0 * rr)

    # EQ stores OMPS as angular velocity; Vphi = Omega * R.
    vphi = omega * raxis

    profiles_csv = out_dir / f"{prefix}_mdleqf0_input_profiles.csv"
    with profiles_csv.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "index",
                "psi_pol",
                "psi_n",
                "P_pa",
                "dP_dpsi",
                "TT_2piBR",
                "dTT_dpsi",
                "Jtor_ave_A_per_m2",
                "T_keV",
                "Omega_rad_per_s",
                "Vphi_m_per_s",
            ]
        )
        for i in range(psip.size):
            writer.writerow(
                [
                    i,
                    f"{psip[i]:.12e}",
                    f"{psi_n[i]:.12e}",
                    f"{p_pa[i]:.12e}",
                    f"{dp_dpsi[i]:.12e}",
                    f"{tt[i]:.12e}",
                    f"{dtt_dpsi[i]:.12e}",
                    f"{jtor_ave[i]:.12e}",
                    f"{t_kev[i]:.12e}",
                    f"{omega[i]:.12e}",
                    f"{vphi[i]:.12e}",
                ]
            )

    p_csv = out_dir / f"{prefix}_P_profile.csv"
    with p_csv.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["index", "psi_pol", "psi_n", "P_pa", "dP_dpsi"])
        for i in range(psip.size):
            writer.writerow(
                [i, f"{psip[i]:.12e}", f"{psi_n[i]:.12e}", f"{p_pa[i]:.12e}", f"{dp_dpsi[i]:.12e}"]
            )

    j_csv = out_dir / f"{prefix}_Jtor_profile.csv"
    with j_csv.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            ["index", "psi_pol", "psi_n", "Jtor_ave_A_per_m2", "TT_2piBR", "dTT_dpsi"]
        )
        for i in range(psip.size):
            writer.writerow(
                [
                    i,
                    f"{psip[i]:.12e}",
                    f"{psi_n[i]:.12e}",
                    f"{jtor_ave[i]:.12e}",
                    f"{tt[i]:.12e}",
                    f"{dtt_dpsi[i]:.12e}",
                ]
            )

    t_csv = out_dir / f"{prefix}_T_profile.csv"
    with t_csv.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["index", "psi_pol", "psi_n", "T_keV"])
        for i in range(psip.size):
            writer.writerow([i, f"{psip[i]:.12e}", f"{psi_n[i]:.12e}", f"{t_kev[i]:.12e}"])

    vphi_csv = out_dir / f"{prefix}_Vphi_profile.csv"
    with vphi_csv.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["index", "psi_pol", "psi_n", "Omega_rad_per_s", "Vphi_m_per_s"])
        for i in range(psip.size):
            writer.writerow(
                [i, f"{psip[i]:.12e}", f"{psi_n[i]:.12e}", f"{omega[i]:.12e}", f"{vphi[i]:.12e}"]
            )

    ip_csv = out_dir / f"{prefix}_mdleqf0_ip.csv"
    with ip_csv.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Ip_MA", "Ip_A"])
        writer.writerow([f"{rip_ma:.12e}", f"{rip_ma * 1.0e6:.12e}"])

    print(f"Saved: {profiles_csv}")
    print(f"Saved: {p_csv}")
    print(f"Saved: {j_csv}")
    print(f"Saved: {t_csv}")
    print(f"Saved: {vphi_csv}")
    print(f"Saved: {ip_csv}")


if __name__ == "__main__":
    main()
