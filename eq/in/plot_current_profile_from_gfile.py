#!/usr/bin/env python3
"""
Plot current-related profiles from a gfile.

Outputs:
  - <prefix>_current_profiles_from_gfile.csv
  - <prefix>_current_midplane_from_gfile.csv
  - <prefix>_current_profiles_from_gfile.png
"""

import argparse
import os
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from solve_gs_from_gfile import MU0, build_lcfs_polygon_mask, read_gfile


def gs_operator(psi, r, z):
    """Compute Delta*psi on interior points for [R,Z]-ordered matrix."""
    dr = float(r[1] - r[0])
    dz = float(z[1] - z[0])
    out = np.full_like(psi, np.nan, dtype=float)

    rc = r[1:-1, None]
    term_rr = (psi[:-2, 1:-1] - 2.0 * psi[1:-1, 1:-1] + psi[2:, 1:-1]) / (dr**2)
    term_r = -(psi[2:, 1:-1] - psi[:-2, 1:-1]) / (2.0 * dr * rc)
    term_zz = (psi[1:-1, :-2] - 2.0 * psi[1:-1, 1:-1] + psi[1:-1, 2:]) / (dz**2)
    out[1:-1, 1:-1] = term_rr + term_r + term_zz
    return out


def rho_from_qpsi(qpsi, psimag, psibdy):
    """
    Compute rho = sqrt(normalized toroidal flux) from q(psi) on uniform psi_N grid.

    Method is consistent with compare_q_geqdsk_namelist.py:
      phi_tor = integral(q * dpsi_poloidal), then rho = sqrt(|phi_tor / phi_tor_edge|).
    """
    q = np.asarray(qpsi, dtype=float)
    n = q.size
    if n == 0:
        return np.asarray([], dtype=float)
    if n == 1:
        return np.asarray([0.0], dtype=float)

    dpsi_pol = float(psibdy - psimag) / float(n - 1)
    phi_tor = np.zeros(n, dtype=float)
    for i in range(1, n):
        phi_tor[i] = phi_tor[i - 1] + 0.5 * (q[i - 1] + q[i]) * dpsi_pol

    edge = float(phi_tor[-1])
    if abs(edge) < 1e-30:
        return np.zeros(n, dtype=float)
    phi_norm = phi_tor / edge
    return np.sqrt(np.abs(phi_norm))


def compute_avebb2_from_gfile(g, psi_n):
    """
    Estimate <B^2>(psi) from gfile on each psi_n bin inside LCFS.

    B is reconstructed using:
      BR = -(1/R) dpsi/dZ, BZ = (1/R) dpsi/dR, BT = F(psi)/R
    where psi is the raw gfile poloidal flux and F=fpol.
    """
    psi = np.asarray(g.psirz, dtype=float)
    rr = np.asarray(g.R, dtype=float)[:, None]
    zz = np.asarray(g.Z, dtype=float)
    npsi = len(psi_n)
    if npsi < 2:
        return np.full(npsi, np.nan, dtype=float)

    dpsi = float(g.psibdy - g.psimag)
    if abs(dpsi) < 1e-30:
        return np.full(npsi, np.nan, dtype=float)

    r_safe = np.maximum(rr, 1e-12)
    dpsi_dr = np.gradient(psi, g.R, axis=0, edge_order=2)
    dpsi_dz = np.gradient(psi, zz, axis=1, edge_order=2)
    br = -dpsi_dz / r_safe
    bz = dpsi_dr / r_safe

    psi_n_2d = np.clip((psi - g.psimag) / dpsi, 0.0, 1.0)
    f_2d = np.interp(psi_n_2d.ravel(), psi_n, g.fpol).reshape(psi_n_2d.shape)
    bt = f_2d / r_safe
    b2 = br * br + bz * bz + bt * bt

    inside = build_lcfs_polygon_mask(g.R, g.Z, g.rbbbs, g.zbbbs)
    if inside is None:
        inside = np.isfinite(b2)
    valid = inside & np.isfinite(b2) & np.isfinite(psi_n_2d)
    if not np.any(valid):
        return np.full(npsi, np.nan, dtype=float)

    dpsi_n = 1.0 / float(npsi - 1)
    edges = np.linspace(-0.5 * dpsi_n, 1.0 + 0.5 * dpsi_n, npsi + 1)
    bins = np.digitize(psi_n_2d[valid], edges) - 1
    bins = np.clip(bins, 0, npsi - 1)
    sums = np.bincount(bins, weights=b2[valid], minlength=npsi).astype(float)
    cnts = np.bincount(bins, minlength=npsi).astype(float)

    avebb2 = np.full(npsi, np.nan, dtype=float)
    good = cnts > 0.0
    avebb2[good] = sums[good] / cnts[good]

    good_idx = np.where(np.isfinite(avebb2))[0]
    if len(good_idx) == 0:
        return avebb2
    if len(good_idx) == 1:
        avebb2[:] = avebb2[good_idx[0]]
        return avebb2
    all_idx = np.arange(npsi)
    avebb2 = np.interp(all_idx, good_idx, avebb2[good_idx])
    return avebb2


def parse_top_level_array(filepath, varname):
    """Parse a top-level array variable from namelist-style text."""
    with open(filepath, "r", encoding="utf-8") as f:
        text = f.read()

    pattern = re.compile(rf"^\s*{re.escape(varname)}\s*=\s*", re.MULTILINE | re.IGNORECASE)
    m = pattern.search(text)
    if m is None:
        return None

    values = []
    rest = text[m.end() :]
    for line in rest.splitlines():
        s = line.strip()
        if not s:
            continue
        if s.startswith("&") or re.match(r"^[A-Za-z_]\w*\s*=", s):
            break
        for token in s.replace(",", " ").split():
            if "*" in token:
                c_str, v_str = token.split("*", 1)
                try:
                    count = int(c_str)
                    val = float(v_str)
                    values.extend([val] * count)
                except ValueError:
                    continue
            else:
                try:
                    values.append(float(token))
                except ValueError:
                    continue

    if not values:
        return None
    return np.asarray(values, dtype=float)


def load_namelist_tot_profile(namelist_path):
    """Load rho and total current-density profile from namelist."""
    if namelist_path is None:
        return None, None, None
    p = os.path.abspath(namelist_path)
    if not os.path.exists(p):
        return None, None, None

    rho = parse_top_level_array(p, "rho")
    if rho is None:
        return None, None, None

    tot = parse_top_level_array(p, "current_density")
    if tot is None:
        # Fallback name used by some outputs.
        tot = parse_top_level_array(p, "curden")
    if tot is None:
        return None, None, None

    n = min(len(rho), len(tot))
    if n <= 1:
        return None, None, None
    return rho[:n], tot[:n], p


def main():
    parser = argparse.ArgumentParser(description="Plot current profiles from gfile.")
    parser.add_argument("--gfile", default="in/g260206.20000_teq_0114")
    parser.add_argument("--out-dir", default="in")
    parser.add_argument("--prefix", default=None)
    parser.add_argument(
        "--r-ref",
        type=float,
        default=8.03,
        help="Reference R (m) for jphi(psi)=R*p'+FF'/(mu0*R). Default: magnetic-axis R.",
    )
    parser.add_argument(
        "--align-ip",
        action="store_true",
        help="Align jphi sign/scale so integrated Ip matches target.",
    )
    parser.add_argument(
        "--target-ip-ma",
        type=float,
        default=None,
        help="Target total plasma current in MA (default: gfile currentA).",
    )
    parser.add_argument(
        "--ffprime-sign",
        type=int,
        choices=[-1, 1],
        default=-1,
        help="Use FF' with sign factor (+1 original, -1 flipped).",
    )
    parser.add_argument(
        "--pprime-sign",
        type=int,
        choices=[-1, 1],
        default=-1,
        help="Use p' with sign factor (+1 original, -1 flipped).",
    )
    parser.add_argument(
        "--x-coordinate",
        choices=["psi_n", "rho_tor"],
        default="rho_tor",
        help="Profile x-axis: psi_n or rho=sqrt(normalized toroidal flux). Default: rho_tor.",
    )
    parser.add_argument(
        "--namelist",
        default="/Users/dengxiaoya/CFEDRSW/OMFIT_out/profiles_CFEDR.namelist",
        help="Namelist path for overlaying tot current_density vs rho.",
    )
    args = parser.parse_args()

    gpath = os.path.abspath(args.gfile)
    if not os.path.exists(gpath):
        raise FileNotFoundError(f"gfile not found: {gpath}")

    g = read_gfile(gpath)
    os.makedirs(args.out_dir, exist_ok=True)
    prefix = args.prefix if args.prefix else os.path.basename(gpath)

    psi_n = np.linspace(0.0, 1.0, g.nw)
    psi_wb = g.psimag + psi_n * (g.psibdy - g.psimag)
    rho_tor = rho_from_qpsi(g.qpsi, g.psimag, g.psibdy)

    pprime_used = args.pprime_sign * g.pprime
    ffprime_used = args.ffprime_sign * g.ffprime
    r_ref = float(args.r_ref) if args.r_ref is not None else float(g.rmaxis)
    b0_ref = abs(float(g.bcentr) * float(g.rcentr) / max(r_ref, 1e-12))
    jphi_ref = r_ref * pprime_used + ffprime_used / (MU0 * r_ref)
    # GS RHS evaluated at reference R.
    rhs_ref = -MU0 * (r_ref**2) * pprime_used - ffprime_used

    # Reconstruct AVEJTR/AVEJPR using eqcalq.f formulas from gfile profiles.
    # TTS = 2*pi*F, DTT = d(TTS)/dpsi = (4*pi^2*FF')/TTS.
    tts = 2.0 * np.pi * np.asarray(g.fpol, dtype=float)
    ttdtt_raw = 4.0 * (np.pi**2) * np.asarray(g.ffprime, dtype=float)
    ttdtt_used = 4.0 * (np.pi**2) * ffprime_used
    dtt_raw = np.divide(ttdtt_raw, tts, out=np.zeros_like(tts), where=np.abs(tts) > 1e-14)
    dtt_used = np.divide(ttdtt_used, tts, out=np.zeros_like(tts), where=np.abs(tts) > 1e-14)
    avebb2_fs = compute_avebb2_from_gfile(g, psi_n)

    avejtr_raw = -r_ref * g.pprime - g.ffprime / (MU0 * r_ref)
    avejtr_used = -r_ref * pprime_used - ffprime_used / (MU0 * r_ref)
    denom_pr = 2.0 * np.pi * max(b0_ref, 1e-12)
    avejpr_raw = (-tts * g.pprime - dtt_raw * avebb2_fs / MU0) / denom_pr
    avejpr_used = (-tts * pprime_used - dtt_used * avebb2_fs / MU0) / denom_pr

    # 2D toroidal current density from psi via GS operator.
    dstar = gs_operator(g.psirz, g.R, g.Z)
    rr = g.R[:, None]
    jphi_2d_raw_a = -(1.0 / (MU0 * rr)) * dstar
    jphi_2d_raw_b = +(1.0 / (MU0 * rr)) * dstar
    lcfs_mask = build_lcfs_polygon_mask(g.R, g.Z, g.rbbbs, g.zbbbs)
    if lcfs_mask is None:
        lcfs_mask = np.isfinite(jphi_2d_raw_a)
    dr = float(g.R[1] - g.R[0])
    dz = float(g.Z[1] - g.Z[0])

    ip_target_a = g.currentA if args.target_ip_ma is None else args.target_ip_ma * 1e6
    ip_a = float(np.nansum(jphi_2d_raw_a[lcfs_mask]) * dr * dz)
    ip_b = float(np.nansum(jphi_2d_raw_b[lcfs_mask]) * dr * dz)

    # Pick sign convention closer to target.
    if abs(ip_a - ip_target_a) <= abs(ip_b - ip_target_a):
        jphi_2d = jphi_2d_raw_a
        ip_raw = ip_a
        sign_label = "-Delta*psi/(mu0 R)"
    else:
        jphi_2d = jphi_2d_raw_b
        ip_raw = ip_b
        sign_label = "+Delta*psi/(mu0 R)"

    scale = 1.0
    if args.align_ip and abs(ip_raw) > 1e-12:
        scale = ip_target_a / ip_raw
    jphi_2d = jphi_2d * scale
    jphi_ref_aligned = jphi_ref * scale
    avejtr_used_aligned = avejtr_used * scale
    avejpr_used_aligned = avejpr_used * scale

    jmid = int(np.argmin(np.abs(g.Z - g.zmaxis)))
    jphi_mid = jphi_2d[:, jmid]
    psi_mid = g.psirz[:, jmid]
    psi_mid_n = (psi_mid - g.psimag) / (g.psibdy - g.psimag + 1e-30)
    rho_mid = np.interp(np.clip(psi_mid_n, 0.0, 1.0), psi_n, rho_tor)

    if args.x_coordinate == "rho_tor":
        x_prof = rho_tor
        x_label = r"$\rho=\sqrt{\Phi_{tor,norm}}$"
    else:
        x_prof = psi_n
        x_label = "psi_n"

    rho_nml, jtot_nml, namelist_used = load_namelist_tot_profile(args.namelist)

    df_prof = pd.DataFrame(
        {
            "psi_n": psi_n,
            "rho_tor": rho_tor,
            "psi_wb": psi_wb,
            "qpsi": g.qpsi,
            "F_fpol_Tm": g.fpol,
            "TTS_2piF": tts,
            "pprime": g.pprime,
            "pprime_used": pprime_used,
            "FFprime": g.ffprime,
            "FFprime_used": ffprime_used,
            "DTT_raw": dtt_raw,
            "DTT_used": dtt_used,
            "AVEBB2_fs_T2": avebb2_fs,
            "jphi_ref_raw_A_per_m2": jphi_ref,
            "jphi_ref_aligned_A_per_m2": jphi_ref_aligned,
            "AVEJTR_raw_A_per_m2": avejtr_raw,
            "AVEJTR_used_A_per_m2": avejtr_used,
            "AVEJTR_used_aligned_A_per_m2": avejtr_used_aligned,
            "AVEJPR_raw_A_per_m2": avejpr_raw,
            "AVEJPR_used_A_per_m2": avejpr_used,
            "AVEJPR_used_aligned_A_per_m2": avejpr_used_aligned,
            "rhs_ref": rhs_ref,
        }
    )
    df_mid = pd.DataFrame(
        {
            "R_m": g.R,
            "Z_mid_m": np.full_like(g.R, g.Z[jmid], dtype=float),
            "psi_mid": psi_mid,
            "psi_mid_n": psi_mid_n,
            "rho_mid": rho_mid,
            "jphi_mid_aligned_A_per_m2": jphi_mid,
        }
    )

    csv1 = os.path.join(args.out_dir, f"{prefix}_current_profiles_from_gfile.csv")
    csv2 = os.path.join(args.out_dir, f"{prefix}_current_midplane_from_gfile.csv")
    df_prof.to_csv(csv1, index=False)
    df_mid.to_csv(csv2, index=False)

    fig, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)

    axes[0, 0].plot(x_prof, pprime_used, lw=2)
    axes[0, 0].set_title(f"p'(psi) used (sign={args.pprime_sign:+d})")
    axes[0, 0].set_xlabel(x_label)
    axes[0, 0].grid(alpha=0.25)

    axes[0, 1].plot(x_prof, ffprime_used, lw=2)
    axes[0, 1].set_title(f"FF'(psi) used (sign={args.ffprime_sign:+d})")
    axes[0, 1].set_xlabel(x_label)
    axes[0, 1].grid(alpha=0.25)

    axes[1, 0].plot(x_prof, jphi_ref, lw=1.2, ls="--", label="jphi_ref raw")
    axes[1, 0].plot(x_prof, jphi_ref_aligned, lw=1.6, label="jphi_ref aligned")
    axes[1, 0].plot(x_prof, avejtr_used_aligned, lw=2.0, label="AVEJTR (eq formula)")
    axes[1, 0].plot(x_prof, avejpr_used_aligned, lw=2.0, label="AVEJPR (eq formula)")

    if rho_nml is not None and jtot_nml is not None:
        # Use namelist rho as x-axis on the current panel; interpolate gfile curves to this x.
        idx = np.argsort(rho_tor)
        rho_sorted = rho_tor[idx]
        j_raw_sorted = jphi_ref[idx]
        j_aligned_sorted = jphi_ref_aligned[idx]
        avejtr_sorted = avejtr_used_aligned[idx]
        avejpr_sorted = avejpr_used_aligned[idx]
        j_raw_on_nml = np.interp(rho_nml, rho_sorted, j_raw_sorted)
        j_aligned_on_nml = np.interp(rho_nml, rho_sorted, j_aligned_sorted)
        avejtr_on_nml = np.interp(rho_nml, rho_sorted, avejtr_sorted)
        avejpr_on_nml = np.interp(rho_nml, rho_sorted, avejpr_sorted)

        axes[1, 0].cla()
        axes[1, 0].plot(rho_nml, j_raw_on_nml, lw=1.2, ls="--", label="jphi_ref raw")
        axes[1, 0].plot(rho_nml, j_aligned_on_nml, lw=1.6, label="jphi_ref aligned")
        axes[1, 0].plot(rho_nml, avejtr_on_nml, lw=2.0, label="AVEJTR (eq formula)")
        axes[1, 0].plot(rho_nml, avejpr_on_nml, lw=2.0, label="AVEJPR (eq formula)")
        axes[1, 0].plot(rho_nml, jtot_nml, "k--", lw=2, label="namelist tot current_density")
        axes[1, 0].set_xlabel("rho (namelist)")
    else:
        axes[1, 0].set_xlabel(x_label)

    axes[1, 0].set_title("Current profiles: jphi_ref / AVEJTR / AVEJPR")
    axes[1, 0].set_ylabel("A/m^2")
    axes[1, 0].grid(alpha=0.25)
    axes[1, 0].legend(fontsize=8, title=f"R_ref={r_ref:.3f} m, B0_ref={b0_ref:.3f} T")

    axes[1, 1].plot(g.R, jphi_mid, lw=2)
    axes[1, 1].set_title(f"jphi(R) at Z={g.Z[jmid]:.3f} m (near axis)")
    axes[1, 1].set_xlabel("R (m)")
    axes[1, 1].set_ylabel("A/m^2")
    axes[1, 1].grid(alpha=0.25)

    fig.suptitle(f"Current Profiles from gfile: {prefix}")
    fig_path = os.path.join(args.out_dir, f"{prefix}_current_profiles_from_gfile.png")
    fig.savefig(fig_path, dpi=160)
    plt.close(fig)

    print(f"gfile: {gpath}")
    print(f"R_ref for jphi(psi): {r_ref:.6f} m")
    print(f"B0_ref for AVEJPR: {b0_ref:.6f} T")
    print(f"profile x-axis = {args.x_coordinate}")
    if namelist_used is not None:
        print(f"namelist overlay = {namelist_used}")
    else:
        print("namelist overlay = disabled/not found")
    print(f"p' sign factor = {args.pprime_sign:+d}")
    print(f"FF' sign factor = {args.ffprime_sign:+d}")
    print(f"target Ip = {ip_target_a:.6e} A")
    print(f"selected sign = {sign_label}")
    print(f"raw Ip (integrated) = {ip_raw:.6e} A")
    print(f"scale = {scale:.6e}")
    print(f"aligned Ip = {ip_raw*scale:.6e} A")
    print(f"saved: {os.path.abspath(csv1)}")
    print(f"saved: {os.path.abspath(csv2)}")
    print(f"saved: {os.path.abspath(fig_path)}")


if __name__ == "__main__":
    main()
