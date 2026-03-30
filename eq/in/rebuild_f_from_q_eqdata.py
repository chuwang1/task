#!/usr/bin/env python3
"""
Rebuild F(psi) from q(psi) using exported EQ metrics.

Two methods:
1) direct:
     q = F * (dV/dpsi_p) * <1/R^2> / (4*pi^2)
     F = 4*pi^2 * q / ((dV/dpsi_p) * <1/R^2>)
2) eqipqp:
     mirror Fortran EQIPQP (eqcalv.f) on the saved EQ grid.
"""

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

FOUR_PI2 = 4.0 * np.pi * np.pi
MU0 = 4.0e-7 * np.pi


def _pick_col(df, candidates, tag):
    for c in candidates:
        if c in df.columns:
            return c
    raise KeyError(f"Cannot find column for {tag}. candidates={candidates}")


def _safe_interp(x_new, x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    idx = np.argsort(x)
    x = x[idx]
    y = y[idx]
    xu, iu = np.unique(x, return_index=True)
    yu = y[iu]
    return np.interp(x_new, xu, yu)


def _safe_grad(y, x):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    idx = np.argsort(x)
    xs = x[idx]
    ys = y[idx]
    xu, iu = np.unique(xs, return_index=True)
    yu = ys[iu]
    if len(xu) < 2:
        return np.full_like(yu, np.nan, dtype=float), xu
    return np.gradient(yu, xu), xu


def _load_single_xy(path, x_cands, y_cands, tag):
    df = pd.read_csv(path)
    cx = _pick_col(df, x_cands, f"{tag}:x")
    cy = _pick_col(df, y_cands, f"{tag}:y")
    return df[cx].to_numpy(dtype=float), df[cy].to_numpy(dtype=float)


def _load_geom_compare(path, psip):
    """Load paired geometry sets (gfile/input) and interpolate to psip grid.

    Expected columns:
      - axis: one of [PSIP, PSIPV] or [psi_n]
      - gfile set: RSV_gf, AVIR2_gf, VPV_gf, AVRR2_gf
      - input set: RSV_ip, AVIR2_ip, VPV_ip, AVRR2_ip
    """
    dg = pd.read_csv(path)

    if "PSIP" in dg.columns:
        x = dg["PSIP"].to_numpy(dtype=float)
    elif "PSIPV" in dg.columns:
        x = dg["PSIPV"].to_numpy(dtype=float)
    elif "psi_n" in dg.columns:
        psip_max = float(np.nanmax(psip))
        x = dg["psi_n"].to_numpy(dtype=float) * psip_max
    else:
        raise KeyError("geom-compare csv requires one axis column: PSIP/PSIPV/psi_n")

    out = {}
    for key in [
        "RSV_gf",
        "AVIR2_gf",
        "VPV_gf",
        "AVRR2_gf",
        "RSV_ip",
        "AVIR2_ip",
        "VPV_ip",
        "AVRR2_ip",
    ]:
        if key not in dg.columns:
            raise KeyError(f"Missing column in geom-compare csv: {key}")
        out[key] = _safe_interp(psip, x, dg[key].to_numpy(dtype=float))
    return out


def _load_gfile_qp(path, psip):
    dg = pd.read_csv(path)
    if "psi_n" in dg.columns:
        x = dg["psi_n"].to_numpy(dtype=float) * float(np.nanmax(psip))
    else:
        x = dg[_pick_col(dg, ["PSIP", "PSIPV"], "gfile q/p axis")].to_numpy(dtype=float)

    q = dg[_pick_col(dg, ["q", "QPSI", "QPV", "QPS"], "gfile q")].to_numpy(dtype=float)

    if "pprime_pa_per_wb" in dg.columns:
        pprime_wb = dg["pprime_pa_per_wb"].to_numpy(dtype=float)
        if "psi_wb" in dg.columns:
            psi_wb = dg["psi_wb"].to_numpy(dtype=float)
            dpsi_wb_dx = np.gradient(psi_wb, x, edge_order=1)
            pprime = pprime_wb * dpsi_wb_dx
        else:
            pprime = pprime_wb
    elif "pprime" in dg.columns:
        pprime = dg["pprime"].to_numpy(dtype=float)
    else:
        pcol = _pick_col(dg, ["p_pa", "P", "p"], "gfile p for gradient")
        p = dg[pcol].to_numpy(dtype=float)
        pprime = np.gradient(p, x, edge_order=1)

    return _safe_interp(psip, x, q), _safe_interp(psip, x, pprime)


def _rebuild_f_eqipqp(psip, q, pprime, avir2, vpv, avrr2, rsv, rr, bb):
    """
    Python mirror of Fortran EQIPQP (eqcalv.f):
      F(N)=2*pi*BB*RR
      integrate inward for F(i-1) from P', q and geometry metrics.
    """
    n = len(psip)
    f = np.full(n, np.nan, dtype=float)
    if n < 2:
        return f

    f[-1] = 2.0 * np.pi * bb * rr

    for i in range(n - 1, 0, -1):
        psip_p = psip[i]
        psip_m = psip[i - 1]
        q_p = q[i]
        q_m = q[i - 1]
        dp_p = pprime[i]
        dp_m = pprime[i - 1]

        if not np.isfinite(
            [q_p, q_m, dp_p, dp_m, avir2[i], avir2[i - 1], vpv[i], vpv[i - 1], avrr2[i], avrr2[i - 1]]
        ).all():
            continue
        if abs(q_p) < 1e-14 or abs(q_m) < 1e-14:
            continue

        xp = rsv[i] / q_p
        xm = rsv[i - 1] / q_m

        alp = FOUR_PI2 * bb * bb * rr * rr * xp / max(avir2[i] * vpv[i], 1e-14)
        blp = FOUR_PI2 * MU0 * rr * rr * dp_p / max(avir2[i], 1e-14)
        flp = vpv[i] * avrr2[i] * xp

        if i == 1:
            alm = FOUR_PI2 * bb * bb * rr * rr * xp / max(avir2[i - 1] * vpv[i], 1e-14)
        else:
            alm = FOUR_PI2 * bb * bb * rr * rr * xm / max(avir2[i - 1] * vpv[i - 1], 1e-14)
        blm = FOUR_PI2 * MU0 * rr * rr * dp_m / max(avir2[i - 1], 1e-14)
        flm = vpv[i - 1] * avrr2[i - 1] * xm

        yp = 0.5 * f[i] * f[i]
        ym = yp + 0.5 * (blp + blm) * (psip_p - psip_m) + 0.5 * (alp + alm) * (flp - flm)
        if not np.isfinite(ym):
            continue
        if ym <= 0.0:
            f[i - 1] = np.nan
        else:
            f[i - 1] = np.sqrt(2.0 * ym)

    return f


def main():
    parser = argparse.ArgumentParser(description="Rebuild F(psi) from q(psi) using eq flux-surface metrics.")
    parser.add_argument("--prefix", default="eqdata_modelg5_rebuild_mdleqf9_qmap")
    parser.add_argument("--radial-csv", default=None, help="Defaults to <prefix>_radial_profiles.csv")
    parser.add_argument("--profile1d-csv", default=None, help="Defaults to <prefix>_1D_profiles.csv")
    parser.add_argument("--dvdpsip-csv", default="eqgs1d_21_DVDPSIP.csv")
    parser.add_argument("--aveir2-csv", default="eqgs1d_10_AVEIR2.csv")
    parser.add_argument("--averr2-csv", default="eqgs1d_09_AVERR2.csv", help="Legacy input (not used by eqipqp).")
    parser.add_argument("--vpv-csv", default="eqgs1d_05_VPS.csv", help="Legacy input (not used by eqipqp).")
    parser.add_argument("--avegvr2-csv", default="eqgs1d_14_AVEGVR2.csv")
    parser.add_argument("--eqipqp-input-csv", default="eqgs1d_24_EQIPQP_INPUTS.csv")
    parser.add_argument(
        "--use-eqipqp-input",
        action="store_true",
        default=True,
        help="Use eqgs1d_24_EQIPQP_INPUTS.csv for DPPSI/QPSI/geometry (default: True).",
    )
    parser.add_argument(
        "--no-eqipqp-input",
        action="store_true",
        help="Force fallback path: reconstruct inputs from eqdata exports (uses QPV, not QPSI).",
    )
    parser.add_argument(
        "--method",
        choices=["direct", "eqipqp"],
        default="eqipqp",
        help="direct: algebraic F from q,dV/dpsi,<1/R^2>; eqipqp: Fortran-like backward integration.",
    )
    parser.add_argument("--out-prefix", default=None, help="Defaults to in/<prefix>_f_from_q")
    parser.add_argument(
        "--geom-compare-csv",
        default=None,
        help=(
            "Optional csv with paired geometry sets (RSV/AVIR2/VPV/AVRR2 for gfile+input) "
            "used to rebuild and compare F under identical q,pprime recursion."
        ),
    )
    parser.add_argument(
        "--use-gfile-qp",
        action="store_true",
        help="Use q and p' from gfile profile csv for eqipqp recursion.",
    )
    parser.add_argument(
        "--gfile-qp-csv",
        default="in/g260206.20000_teq_0114_profiles_from_gfile.csv",
        help="Gfile profile csv containing psi_n,q,pprime data.",
    )
    args = parser.parse_args()

    if args.no_eqipqp_input:
        args.use_eqipqp_input = False

    radial_csv = args.radial_csv or f"{args.prefix}_radial_profiles.csv"
    profile1d_csv = args.profile1d_csv or f"{args.prefix}_1D_profiles.csv"
    params_csv = f"{args.prefix}_parameters.csv"
    out_prefix = args.out_prefix or os.path.join("in", f"{args.prefix}_f_from_q")
    out_dir = os.path.dirname(os.path.abspath(out_prefix))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    dr = pd.read_csv(radial_csv)
    ddv = pd.read_csv(args.dvdpsip_csv)
    dir2 = pd.read_csv(args.aveir2_csv)
    dpar = pd.read_csv(params_csv)

    c_psip = _pick_col(dr, ["PSIPV", "PSIP", "PSIPS"], "psi_p")
    c_q = _pick_col(dr, ["QPV", "QPS", "q"], "q")
    c_ttv = None
    for c in ["TTV", "TTS", "F"]:
        if c in dr.columns:
            c_ttv = c
            break

    c_psip_dv = _pick_col(ddv, ["PSIP", "PSIPV"], "PSIP for dV/dpsi")
    c_dv = _pick_col(ddv, ["DVDPSIP_m3perWb", "DVDPSIP"], "dV/dpsi")
    c_psip_ir2 = _pick_col(dir2, ["PSIP", "PSIPV"], "PSIP for <1/R^2>")
    c_ir2 = _pick_col(dir2, ["AVEIR2_per_m2", "AVEIR2"], "<1/R^2>")

    psip = dr[c_psip].to_numpy(dtype=float)
    q = dr[c_q].to_numpy(dtype=float)
    dv = _safe_interp(psip, ddv[c_psip_dv].to_numpy(dtype=float), ddv[c_dv].to_numpy(dtype=float))
    ir2 = _safe_interp(psip, dir2[c_psip_ir2].to_numpy(dtype=float), dir2[c_ir2].to_numpy(dtype=float))

    denom = dv * ir2
    good = np.isfinite(psip) & np.isfinite(q) & np.isfinite(denom) & (np.abs(denom) > 1e-14)
    f_direct = np.full_like(q, np.nan, dtype=float)
    f_direct[good] = FOUR_PI2 * q[good] / denom[good]  # F_direct always uses q(psi) from radial CSV (not q_use)

    fipv_e = None  # Fortran EQIPQP result, set only if available
    f_rec_gfgeom = None
    f_rec_ipgeom = None
    q_gfile_out = np.full_like(q, np.nan, dtype=float)  # populated if --use-gfile-qp
    if args.method == "direct":
        f_rec = f_direct.copy()
    else:
        pmap = {str(k): float(v) for k, v in zip(dpar["Parameter"], dpar["Value"])}
        rr = pmap["RR"]
        bb = pmap["BB"]
        used_exact_file = False
        fipv_e = None
        q_use = q.copy()
        pprime_i = np.full_like(psip, np.nan, dtype=float)
        rsv = np.full_like(psip, np.nan, dtype=float)
        avir2_calv = np.full_like(psip, np.nan, dtype=float)
        vpv_calv = np.full_like(psip, np.nan, dtype=float)
        avrr2_calv = np.full_like(psip, np.nan, dtype=float)

        if args.use_eqipqp_input and os.path.exists(args.eqipqp_input_csv):
            de = pd.read_csv(args.eqipqp_input_csv)
            c_psip_e = _pick_col(de, ["PSIPV", "PSIP"], "EQIPQP inputs: PSIPV")
            c_dp_e = _pick_col(de, ["DPPSI_PaPerWb", "DPPSI"], "EQIPQP inputs: DPPSI")
            c_q_e = _pick_col(de, ["QPSI", "QPV", "QPS"], "EQIPQP inputs: Q")
            c_rsv_e = _pick_col(de, ["RSV_m", "RSV"], "EQIPQP inputs: RSV")
            c_avir2_e = _pick_col(de, ["AVIR2"], "EQIPQP inputs: AVIR2")
            c_vpv_e = _pick_col(de, ["VPV_m3", "VPV"], "EQIPQP inputs: VPV")
            c_avrr2_e = _pick_col(de, ["AVRR2_m2", "AVRR2"], "EQIPQP inputs: AVRR2")

            psip_e = de[c_psip_e].to_numpy(dtype=float)
            dp_e = de[c_dp_e].to_numpy(dtype=float)
            q_e = de[c_q_e].to_numpy(dtype=float)
            rsv_e = de[c_rsv_e].to_numpy(dtype=float)
            avir2_e = de[c_avir2_e].to_numpy(dtype=float)
            vpv_e = de[c_vpv_e].to_numpy(dtype=float)
            avrr2_e = de[c_avrr2_e].to_numpy(dtype=float)

            exact_ok = (
                np.isfinite(dp_e).all()
                and np.isfinite(q_e).all()
                and np.isfinite(rsv_e).all()
                and np.isfinite(avir2_e).all()
                and np.isfinite(vpv_e).all()
                and np.isfinite(avrr2_e).all()
                and np.count_nonzero(np.abs(vpv_e) > 1e-14) > max(3, len(vpv_e) // 4)
                and np.count_nonzero(np.abs(avir2_e) > 1e-14) > max(3, len(avir2_e) // 4)
            )
            # Also try to read FIPV from Fortran EQIPQP for reference
            c_fipv_e = None
            for c in ["FIPV", "FIPV_Tm"]:
                if c in de.columns:
                    c_fipv_e = c
                    break
            fipv_e = None
            if c_fipv_e is not None:
                _fipv_raw = de[c_fipv_e].to_numpy(dtype=float)
                if np.isfinite(_fipv_raw).any() and np.count_nonzero(np.abs(_fipv_raw) > 1e-10) > 3:
                    fipv_e = _safe_interp(psip, psip_e, _fipv_raw)

            if exact_ok:
                pprime_i = _safe_interp(psip, psip_e, dp_e)
                q_use = _safe_interp(psip, psip_e, q_e)
                rsv = _safe_interp(psip, psip_e, rsv_e)
                avir2_calv = _safe_interp(psip, psip_e, avir2_e)
                vpv_calv = _safe_interp(psip, psip_e, vpv_e)
                avrr2_calv = _safe_interp(psip, psip_e, avrr2_e)
                used_exact_file = True
            else:
                print(
                    f"[eqipqp] exact file looks uninitialized (likely zero RSV/VPV): {args.eqipqp_input_csv}; fallback."
                )
        elif args.use_eqipqp_input:
            print(f"[eqipqp] exact file not found: {args.eqipqp_input_csv}; fallback.")

        if not used_exact_file:
            # Fallback path: reconstruct missing EQCALV metrics from EQCALQ exports.
            p1d = pd.read_csv(profile1d_csv)
            c_psips = _pick_col(p1d, ["PSIPS", "PSIP", "PSIPV"], "PSIPS in 1D profile")
            c_ppps = _pick_col(p1d, ["PPPS", "p", "P"], "PPPS in 1D profile")
            psips = p1d[c_psips].to_numpy(dtype=float)
            p = p1d[c_ppps].to_numpy(dtype=float)
            psitv = dr["PSITV"].to_numpy(dtype=float) if "PSITV" in dr.columns else None
            if psitv is None:
                raise ValueError("eqipqp requires PSITV in radial CSV to align with Fortran mapping.")

            psipa = float(np.nanmax(psip))
            psita = float(np.nanmax(psitv))
            psipn_r = np.clip(psip / max(psipa, 1e-14), 0.0, 1.0)
            psitn_r = psitv / max(psita, 1e-14)

            psipn_p = np.clip(psips / max(psipa, 1e-14), 0.0, 1.0)
            psit_p = _safe_interp(psipn_p, psipn_r, psitv)
            psitn_p = psit_p / max(psita, 1e-14)
            p_mpa = p * 1.0e-6
            dp_dpsitn_u, psitn_u = _safe_grad(p_mpa, psitn_p)
            dp_dpsitn = _safe_interp(psitn_r, psitn_u, dp_dpsitn_u)
            pprime_i = (q / max(psita, 1e-14)) * dp_dpsitn * 1.0e6
            q_use = q.copy()

            x_dv, y_dv = _load_single_xy(
                args.dvdpsip_csv, ["PSIP", "PSIPV"], ["DVDPSIP_m3perWb", "DVDPSIP"], "DVDPSIP"
            )
            x_aveir2, y_aveir2 = _load_single_xy(
                args.aveir2_csv, ["PSIP", "PSIPV"], ["AVEIR2_per_m2", "AVEIR2"], "AVEIR2"
            )
            x_avegvr2, y_avegvr2 = _load_single_xy(
                args.avegvr2_csv, ["PSIP", "PSIPV"], ["AVEGVR2_m2", "AVEGVR2"], "AVEGVR2"
            )
            dvdpsip_i = _safe_interp(psip, x_dv, y_dv)
            aveir2_i = _safe_interp(psip, x_aveir2, y_aveir2)
            avegvr2_i = _safe_interp(psip, x_avegvr2, y_avegvr2)

            rsv = np.sqrt(np.abs(psitv / max(bb * np.pi, 1e-14)))
            avir2_calv = aveir2_i * rr * rr

            vpv_calv = np.full_like(psip, np.nan, dtype=float)
            ok_q = np.abs(q_use) > 1e-14
            vpv_calv[ok_q] = 2.0 * np.pi * rsv[ok_q] * bb * dvdpsip_i[ok_q] / q_use[ok_q]

            sumavgr2 = np.full_like(psip, np.nan, dtype=float)
            ok_dv = np.abs(dvdpsip_i) > 1e-14
            sumavgr2[ok_dv] = avegvr2_i[ok_dv] / (dvdpsip_i[ok_dv] * FOUR_PI2)

            avrr2_calv = np.full_like(psip, np.nan, dtype=float)
            den = (rsv * rsv) * (bb * bb) * dvdpsip_i
            ok = np.abs(den) > 1e-14
            avrr2_calv[ok] = sumavgr2[ok] * q_use[ok] * q_use[ok] / den[ok]

        f_rec = _rebuild_f_eqipqp(psip, q_use, pprime_i, avir2_calv, vpv_calv, avrr2_calv, rsv, rr, bb)
        if args.use_gfile_qp:
            q_gf, pprime_gf = _load_gfile_qp(args.gfile_qp_csv, psip)
            q_use = q_gf
            q_gfile_out = q_gf.copy()
            pprime_i = pprime_gf
            f_rec = _rebuild_f_eqipqp(psip, q_use, pprime_i, avir2_calv, vpv_calv, avrr2_calv, rsv, rr, bb)
            print(f"[eqipqp] using gfile q/p' from: {args.gfile_qp_csv}")
        if args.geom_compare_csv:
            geo = _load_geom_compare(args.geom_compare_csv, psip)
            f_rec_gfgeom = _rebuild_f_eqipqp(
                psip,
                q_use,
                pprime_i,
                geo["AVIR2_gf"],
                geo["VPV_gf"],
                geo["AVRR2_gf"],
                geo["RSV_gf"],
                rr,
                bb,
            )
            f_rec_ipgeom = _rebuild_f_eqipqp(
                psip,
                q_use,
                pprime_i,
                geo["AVIR2_ip"],
                geo["VPV_ip"],
                geo["AVRR2_ip"],
                geo["RSV_ip"],
                rr,
                bb,
            )

            v_cmp = np.isfinite(f_rec_gfgeom) & np.isfinite(f_rec_ipgeom)
            if np.count_nonzero(v_cmp) > 0:
                rmse_cmp = float(np.sqrt(np.mean((f_rec_gfgeom[v_cmp] - f_rec_ipgeom[v_cmp]) ** 2)))
                print(f"[eqipqp] F compare (gfile-geom vs input-geom): RMSE={rmse_cmp:.6e}, N={int(np.count_nonzero(v_cmp))}")

        if used_exact_file:
            print(f"[eqipqp] using exact Fortran inputs from: {args.eqipqp_input_csv}")
            if fipv_e is not None:
                d_fipv = f_rec - fipv_e
                v_fipv = np.isfinite(f_rec) & np.isfinite(fipv_e)
                n_fipv = int(np.count_nonzero(v_fipv))
                if n_fipv > 0:
                    rmse_fipv = float(np.sqrt(np.mean((d_fipv[v_fipv]) ** 2)))
                    print(f"[eqipqp] F vs Fortran FIPV: RMSE={rmse_fipv:.6e}, N={n_fipv}")

    if c_ttv is not None:
        f_ref = dr[c_ttv].to_numpy(dtype=float)
        d = f_rec - f_ref
        valid = np.isfinite(f_rec) & np.isfinite(f_ref)
        n_valid = int(np.count_nonzero(valid))
        if n_valid == 0:
            print(f"[{args.method}] No valid points for error metrics (all NaN).")
        else:
            rmse = float(np.sqrt(np.mean((d[valid]) ** 2)))
            mae = float(np.mean(np.abs(d[valid])))
            print(f"[{args.method}] F compare vs {c_ttv}: RMSE={rmse:.6e}, MAE={mae:.6e}, N={n_valid}/{len(d)}")
            if n_valid < len(d):
                print(f"[{args.method}] Warning: {len(d)-n_valid} points are NaN/invalid in comparison.")
        d_direct = f_direct - f_ref
        valid_d = np.isfinite(f_direct) & np.isfinite(f_ref)
        n_valid_d = int(np.count_nonzero(valid_d))
        if n_valid_d > 0:
            rmse_d = float(np.sqrt(np.mean((d_direct[valid_d]) ** 2)))
            mae_d = float(np.mean(np.abs(d_direct[valid_d])))
            print(f"[direct] F compare vs {c_ttv}: RMSE={rmse_d:.6e}, MAE={mae_d:.6e}, N={n_valid_d}/{len(d_direct)}")

        if f_rec_gfgeom is not None:
            vg = np.isfinite(f_rec_gfgeom) & np.isfinite(f_ref)
            if np.count_nonzero(vg) > 0:
                rmse_g = float(np.sqrt(np.mean((f_rec_gfgeom[vg] - f_ref[vg]) ** 2)))
                print(f"[eqipqp] F(gfile-geom) vs {c_ttv}: RMSE={rmse_g:.6e}, N={int(np.count_nonzero(vg))}")
        if f_rec_ipgeom is not None:
            vi = np.isfinite(f_rec_ipgeom) & np.isfinite(f_ref)
            if np.count_nonzero(vi) > 0:
                rmse_i = float(np.sqrt(np.mean((f_rec_ipgeom[vi] - f_ref[vi]) ** 2)))
                print(f"[eqipqp] F(input-geom) vs {c_ttv}: RMSE={rmse_i:.6e}, N={int(np.count_nonzero(vi))}")
    else:
        f_ref = np.full_like(f_rec, np.nan)
        print("No reference F column found in radial CSV; skip error metrics.")

    q_source_str = "gfile" if args.use_gfile_qp else "radial"
    out_df = pd.DataFrame(
        {
            "PSIP": psip,
            "q": q,
            "q_gfile": q_gfile_out,
            "q_source": q_source_str,
            "dvdpsip": dv,
            "aveir2": ir2,
            "F_direct": f_direct,
            "F_rebuilt": f_rec,
            "F_rebuilt_gfile_geom": f_rec_gfgeom if f_rec_gfgeom is not None else np.full_like(f_rec, np.nan),
            "F_rebuilt_input_geom": f_rec_ipgeom if f_rec_ipgeom is not None else np.full_like(f_rec, np.nan),
            "F_ref": f_ref,
            "dF_direct": f_direct - f_ref,
            "dF": f_rec - f_ref,
            "method": args.method,
        }
    )
    out_csv = f"{out_prefix}.csv"
    out_df.to_csv(out_csv, index=False)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)

    # Check if fipv_e is available (Fortran EQIPQP result)
    has_fipv = args.method == "eqipqp" and fipv_e is not None and np.isfinite(fipv_e).any()

    axes[0].plot(psip, f_rec, "r--", lw=2.0, label=f"F rebuilt from q ({args.method})")
    if f_rec_gfgeom is not None:
        axes[0].plot(psip, f_rec_gfgeom, color="tab:blue", ls="-.", lw=1.6, label="F rebuilt (gfile geom)")
    if f_rec_ipgeom is not None:
        axes[0].plot(psip, f_rec_ipgeom, color="tab:green", ls=":", lw=1.8, label="F rebuilt (input geom)")
    if np.isfinite(f_ref).any():
        axes[0].plot(psip, f_ref, "k-", lw=1.6, label=f"F ref ({c_ttv})")
    if has_fipv:
        axes[0].plot(psip, fipv_e, color="tab:purple", ls=":", lw=1.6, label="F Fortran EQIPQP")
    axes[0].set_xlabel("psi_p")
    axes[0].set_ylabel("F(psi)")
    axes[0].set_title("F(psi) Rebuild from q")
    axes[0].set_yscale('log')
    axes[0].grid(alpha=0.3)
    axes[0].legend()

    axes[1].plot(psip, f_rec - f_ref, color="tab:blue", lw=1.8, label=f"{args.method} - ref")
    if f_rec_gfgeom is not None:
        axes[1].plot(psip, f_rec_gfgeom - f_ref, color="tab:purple", ls="-.", lw=1.4, label="gfile geom - ref")
    if f_rec_ipgeom is not None:
        axes[1].plot(psip, f_rec_ipgeom - f_ref, color="tab:green", ls=":", lw=1.6, label="input geom - ref")
    if has_fipv:
        axes[1].plot(psip, f_rec - fipv_e, color="tab:purple", ls=":", lw=1.6, label=f"{args.method} - Fortran EQIPQP")
    axes[1].axhline(0.0, color="k", ls="--", lw=1.0)
    axes[1].set_xlabel("psi_p")
    axes[1].set_ylabel("dF")
    axes[1].set_title("Rebuild Error (F_rebuilt - F_ref)")
    axes[1].grid(alpha=0.3)
    axes[1].legend()

    out_png = f"{out_prefix}.png"
    fig.savefig(out_png, dpi=160)
    plt.close(fig)

    print(f"saved: {out_csv}")
    print(f"saved: {out_png}")


if __name__ == "__main__":
    main()
