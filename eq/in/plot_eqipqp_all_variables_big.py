#!/usr/bin/env python3
"""Plot EQIPQP-related variables into one large overview figure."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.path import Path as MplPath

from solve_gs_from_gfile import read_gfile

try:
    from contourpy import contour_generator
except Exception:  # pragma: no cover
    contour_generator = None

TWO_PI = 2.0 * np.pi
MU0 = 4.0e-7 * np.pi
EPS = 1.0e-14


def _col(df: pd.DataFrame, candidates: list[str], tag: str) -> str:
    for c in candidates:
        if c in df.columns:
            return c
    raise KeyError(f"Missing column for {tag}. candidates={candidates}")


def _safe_interp(x: np.ndarray, y: np.ndarray, x_new: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x_new = np.asarray(x_new, dtype=float)
    good = np.isfinite(x) & np.isfinite(y)
    if np.count_nonzero(good) < 2:
        return np.full_like(x_new, np.nan, dtype=float)
    xg = x[good]
    yg = y[good]
    idx = np.argsort(xg)
    xg = xg[idx]
    yg = yg[idx]
    xu, inv = np.unique(xg, return_inverse=True)
    if len(xu) != len(xg):
        acc = np.bincount(inv, weights=yg)
        cnt = np.bincount(inv)
        yg = acc / np.maximum(cnt, 1.0)
        xg = xu
    if len(xg) < 2:
        return np.full_like(x_new, np.nan, dtype=float)
    return np.interp(x_new, xg, yg)


def _plot_line(ax, x, y, *args, positive=False, **kwargs):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    if positive:
        m &= y > 0.0
    if np.count_nonzero(m) < 2:
        return
    ax.plot(x[m], y[m], *args, **kwargs)


def _load_external_geom(path: Path, psi_n_target: np.ndarray) -> dict[str, np.ndarray] | None:
    if not path.exists():
        print(f"warning: external geometry csv not found: {path}")
        return None
    try:
        dg = pd.read_csv(path)
    except Exception as exc:
        print(f"warning: failed to read external geometry csv {path}: {exc}")
        return None

    psi_col = _col(dg, ["psi_n", "PSIPNV"], "external geometry psi_n")
    psi = dg[psi_col].to_numpy(float)

    def _pick(cands: list[str], tag: str) -> np.ndarray:
        c = _col(dg, cands, tag)
        return _safe_interp(psi, dg[c].to_numpy(float), psi_n_target)

    return {
        "rsv": _pick(["RSV_ip", "RSV_omfit", "RSV_new"], "external RSV"),
        "avir2": _pick(["AVIR2_ip", "AVIR2_omfit", "AVIR2_new"], "external AVIR2"),
        "vpv": _pick(["VPV_ip", "VPV_omfit", "VPV_new"], "external VPV"),
        "avrr2": _pick(["AVRR2_ip", "AVRR2_omfit_proxy", "AVRR2_new"], "external AVRR2"),
    }


def _load_gfile_profiles(
    path: Path, g, psip_max: float, f_sign: float, pprime_sign: float
) -> dict[str, np.ndarray] | None:
    dg = None
    if path.exists():
        try:
            dg = pd.read_csv(path)
        except Exception as exc:  # pragma: no cover
            print(f"warning: failed to read gfile profile csv {path}: {exc}")
            dg = None
    else:
        print(f"warning: gfile profile csv not found, fallback to raw gfile arrays: {path}")

    if dg is not None:
        psi_col = _col(dg, ["psi_n"], "gfile psi_n")
        q_col = _col(dg, ["q"], "gfile q")
        p_col = _col(dg, ["p_pa"], "gfile p")
        if "F_2pi" in dg.columns:
            f = dg["F_2pi"].to_numpy(float)
        elif "F_tesla_meter" in dg.columns:
            f = TWO_PI * dg["F_tesla_meter"].to_numpy(float)
        else:
            raise KeyError("Missing gfile F column: require F_2pi or F_tesla_meter.")
        psi_n = dg[psi_col].to_numpy(float)
        q = dg[q_col].to_numpy(float)
        p = dg[p_col].to_numpy(float)
        x_tmp = psi_n * psip_max
        if "pprime_pa_per_wb" in dg.columns:
            pprime_wb = dg["pprime_pa_per_wb"].to_numpy(float)
            if "psi_wb" in dg.columns:
                psi_wb = dg["psi_wb"].to_numpy(float)
                dpsi_wb_dx = np.gradient(psi_wb, x_tmp, edge_order=1)
                pprime = pprime_wb * dpsi_wb_dx
            else:
                pprime = pprime_wb
        else:
            # Finite-difference fallback on current x convention.
            pprime = np.gradient(p, x_tmp, edge_order=1)
        ffprime = dg["FFprime"].to_numpy(float) if "FFprime" in dg.columns else np.full_like(pprime, np.nan)
    else:
        psi_n = np.linspace(0.0, 1.0, g.nw, dtype=float)
        q = np.asarray(g.qpsi, dtype=float)
        p = np.asarray(g.pres, dtype=float)
        pprime = np.asarray(g.pprime, dtype=float)
        ffprime = np.asarray(g.ffprime, dtype=float)
        f = TWO_PI * np.asarray(g.fpol, dtype=float)

    f = f_sign * f
    pprime = pprime_sign * pprime
    ffprime = f_sign * ffprime
    valid = np.isfinite(psi_n) & np.isfinite(q) & np.isfinite(p) & np.isfinite(f)
    if np.count_nonzero(valid) < 3:
        print(f"warning: too few valid points in gfile profiles: {path}")
        return None

    psi_n = psi_n[valid]
    q = q[valid]
    p = p[valid]
    f = f[valid]
    pprime = pprime[valid]
    ffprime = ffprime[valid]

    idx = np.argsort(psi_n)
    psi_n = psi_n[idx]
    q = q[idx]
    p = p[idx]
    f = f[idx]
    pprime = pprime[idx]
    ffprime = ffprime[idx]

    psi_u, inv = np.unique(psi_n, return_inverse=True)
    if len(psi_u) != len(psi_n):
        def _avg(a):
            s = np.bincount(inv, weights=a)
            c = np.bincount(inv)
            return s / np.maximum(c, 1.0)

        q = _avg(q)
        p = _avg(p)
        f = _avg(f)
        pprime = _avg(pprime)
        ffprime = _avg(ffprime)
        psi_n = psi_u

    x_psip = psi_n * psip_max
    return {
        "psi_n": psi_n,
        "x": x_psip,
        "q": q,
        "p": p,
        "f": f,
        "pprime": pprime,
        "ffprime": ffprime,
    }


def _integrate_psit_from_q(psip: np.ndarray, q: np.ndarray) -> np.ndarray:
    psip = np.asarray(psip, dtype=float)
    q = np.asarray(q, dtype=float)
    n = len(psip)
    out = np.full(n, np.nan, dtype=float)
    if n == 0:
        return out
    out[0] = 0.0
    for i in range(1, n):
        if not (np.isfinite(psip[i]) and np.isfinite(psip[i - 1]) and np.isfinite(q[i]) and np.isfinite(q[i - 1])):
            continue
        dq = q[i] + q[i - 1]
        if abs(dq) < EPS:
            qh = 0.5 * (q[i] + q[i - 1])
        else:
            qh = 2.0 * q[i] * q[i - 1] / dq
        out[i] = out[i - 1] + qh * (psip[i] - psip[i - 1])
    return out


def _interp_bilinear(field: np.ndarray, r: np.ndarray, z: np.ndarray, rp: np.ndarray, zp: np.ndarray) -> np.ndarray:
    field = np.asarray(field, dtype=float)
    r = np.asarray(r, dtype=float)
    z = np.asarray(z, dtype=float)
    rp = np.asarray(rp, dtype=float)
    zp = np.asarray(zp, dtype=float)
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


def _extract_contour_segment(psi_total: np.ndarray, r: np.ndarray, z: np.ndarray, level: float, axis_r: float, axis_z: float):
    segs = []
    if contour_generator is not None:
        cg = contour_generator(x=r, y=z, z=psi_total.T)
        segs = cg.lines(level)
    else:  # pragma: no cover
        fig, ax = plt.subplots(1, 1)
        try:
            rr2d, zz2d = np.meshgrid(r, z, indexing="xy")
            cs = ax.contour(rr2d, zz2d, psi_total.T, levels=[level])
            if cs.allsegs and cs.allsegs[0]:
                segs = cs.allsegs[0]
        finally:
            plt.close(fig)

    if not segs:
        return None
    best = None
    best_len = -1.0
    for seg in segs:
        seg = np.asarray(seg, dtype=float)
        if seg.ndim != 2 or seg.shape[0] < 8:
            continue
        d = np.diff(seg, axis=0)
        plen = float(np.sum(np.hypot(d[:, 0], d[:, 1])))
        if plen <= 0.0:
            continue
        contains = MplPath(seg).contains_point((axis_r, axis_z))
        score = plen + (1.0e9 if contains else 0.0)
        if score > best_len:
            best_len = score
            best = seg
    return best


def _compute_eqcalv_like_from_gfile(
    g,
    psi_n_levels: np.ndarray,
    psip_levels: np.ndarray,
    q_levels: np.ndarray,
    f_levels: np.ndarray,
    rr: float,
    bb: float,
) -> dict[str, np.ndarray]:
    psi_n_levels = np.asarray(psi_n_levels, dtype=float)
    psip_levels = np.asarray(psip_levels, dtype=float)
    q_levels = np.asarray(q_levels, dtype=float)
    f_levels = np.asarray(f_levels, dtype=float)
    n = len(psi_n_levels)
    out = {
        "psitv": np.full(n, np.nan, dtype=float),
        "rsv": np.full(n, np.nan, dtype=float),
        "avir2": np.full(n, np.nan, dtype=float),
        "vpv": np.full(n, np.nan, dtype=float),
        "avrr2": np.full(n, np.nan, dtype=float),
        "dvdpsip": np.full(n, np.nan, dtype=float),
        "aveir2": np.full(n, np.nan, dtype=float),
    }
    if n < 2 or not np.isfinite(rr) or not np.isfinite(bb) or abs(bb) < EPS:
        return out

    psi_total = TWO_PI * (np.asarray(g.psirz, dtype=float) - float(g.psimag))
    dpsi_dr = np.gradient(psi_total, np.asarray(g.R, dtype=float), axis=0, edge_order=2)
    dpsi_dz = np.gradient(psi_total, np.asarray(g.Z, dtype=float), axis=1, edge_order=2)
    dpsi_total = TWO_PI * (float(g.psibdy) - float(g.psimag))
    if abs(dpsi_total) < EPS:
        return out

    psitv = _integrate_psit_from_q(psip_levels, q_levels)
    out["psitv"] = psitv
    out["rsv"] = np.sqrt(np.abs(psitv / (np.pi * bb)))

    for i in range(n):
        if not (np.isfinite(psi_n_levels[i]) and np.isfinite(out["rsv"][i]) and np.isfinite(q_levels[i]) and np.isfinite(f_levels[i])):
            continue
        level = psi_n_levels[i] * dpsi_total
        seg = _extract_contour_segment(psi_total, g.R, g.Z, level, g.rmaxis, g.zmaxis)
        if seg is None:
            continue
        if np.hypot(seg[0, 0] - seg[-1, 0], seg[0, 1] - seg[-1, 1]) > 1e-10:
            seg = np.vstack([seg, seg[0]])

        p0 = seg[:-1]
        p1 = seg[1:]
        rm = 0.5 * (p0[:, 0] + p1[:, 0])
        zm = 0.5 * (p0[:, 1] + p1[:, 1])
        dl = np.hypot(p1[:, 0] - p0[:, 0], p1[:, 1] - p0[:, 1])
        dpr = _interp_bilinear(dpsi_dr, g.R, g.Z, rm, zm)
        dpz = _interp_bilinear(dpsi_dz, g.R, g.Z, rm, zm)

        bpl = np.sqrt(dpr * dpr + dpz * dpz) / np.maximum(TWO_PI * rm, EPS)
        btl = f_levels[i] / np.maximum(TWO_PI * rm, EPS)
        b2 = bpl * bpl + btl * btl
        valid = np.isfinite(dl) & np.isfinite(rm) & np.isfinite(bpl) & np.isfinite(b2) & (dl > 0.0) & (bpl > EPS) & (rm > EPS)
        if np.count_nonzero(valid) < 8:
            continue

        dlv = dl[valid]
        rv = rm[valid]
        bplv = bpl[valid]
        b2v = b2[valid]

        sumv = np.sum(dlv / bplv)
        sumavir2 = np.sum(dlv / (bplv * rv * rv))
        sumavrr2 = np.sum(dlv * bplv)
        if not np.isfinite(sumv) or sumv <= EPS:
            continue

        out["dvdpsip"][i] = sumv
        out["aveir2"][i] = sumavir2 / sumv
        out["avir2"][i] = out["aveir2"][i] * rr * rr

        ql = q_levels[i]
        rsv = out["rsv"][i]
        if abs(ql) > EPS and rsv > EPS:
            out["vpv"][i] = TWO_PI * rsv * bb * sumv / ql
            den = rsv * rsv * bb * bb * sumv
            if abs(den) > EPS:
                out["avrr2"][i] = sumavrr2 * ql * ql / den

    # Axis-end treatment consistent with EQCALV behavior.
    if n >= 3:
        if np.isfinite(out["avir2"][1]) and np.isfinite(out["avir2"][2]):
            out["avir2"][0] = (4.0 * out["avir2"][1] - out["avir2"][2]) / 3.0
        if np.isfinite(out["avrr2"][1]) and np.isfinite(out["avrr2"][2]):
            out["avrr2"][0] = (4.0 * out["avrr2"][1] - out["avrr2"][2]) / 3.0
        out["vpv"][0] = 0.0
        out["rsv"][0] = 0.0
    return out


def _rebuild_f_eqipqp(
    psip: np.ndarray,
    q: np.ndarray,
    pprime: np.ndarray,
    avir2: np.ndarray,
    vpv: np.ndarray,
    avrr2: np.ndarray,
    rsv: np.ndarray,
    rr: float,
    bb: float,
) -> np.ndarray:
    n = len(psip)
    f = np.full(n, np.nan, dtype=float)
    if n < 2 or not np.isfinite(rr) or not np.isfinite(bb):
        return f

    f[-1] = TWO_PI * bb * rr
    four_pi2 = 4.0 * np.pi * np.pi

    for i in range(n - 1, 0, -1):
        if not np.isfinite(f[i]):
            continue
        psip_p = psip[i]
        psip_m = psip[i - 1]
        q_p = q[i]
        q_m = q[i - 1]
        dp_p = pprime[i]
        dp_m = pprime[i - 1]

        vals = np.array([q_p, q_m, dp_p, dp_m, avir2[i], avir2[i - 1], vpv[i], vpv[i - 1], avrr2[i], avrr2[i - 1]])
        if not np.isfinite(vals).all():
            continue
        if abs(q_p) < EPS or abs(q_m) < EPS:
            continue

        xp = rsv[i] / q_p
        xm = rsv[i - 1] / q_m

        alp = four_pi2 * bb * bb * rr * rr * xp / max(avir2[i] * vpv[i], EPS)
        blp = four_pi2 * MU0 * rr * rr * dp_p / max(avir2[i], EPS)
        flp = vpv[i] * avrr2[i] * xp

        if i == 1:
            alm = four_pi2 * bb * bb * rr * rr * xp / max(avir2[i - 1] * vpv[i], EPS)
        else:
            alm = four_pi2 * bb * bb * rr * rr * xm / max(avir2[i - 1] * vpv[i - 1], EPS)
        blm = four_pi2 * MU0 * rr * rr * dp_m / max(avir2[i - 1], EPS)
        flm = vpv[i - 1] * avrr2[i - 1] * xm

        yp = 0.5 * f[i] * f[i]
        ym = yp + 0.5 * (blp + blm) * (psip_p - psip_m) + 0.5 * (alp + alm) * (flp - flm)
        if not np.isfinite(ym) or ym <= 0.0:
            continue
        f[i - 1] = np.sqrt(2.0 * ym)

    return f


def _calc_abg(
    q: np.ndarray,
    pprime: np.ndarray,
    rsv: np.ndarray,
    avir2: np.ndarray,
    vpv: np.ndarray,
    avrr2: np.ndarray,
    rr: float,
    bb: float,
) -> dict[str, np.ndarray]:
    q = np.asarray(q, dtype=float)
    pprime = np.asarray(pprime, dtype=float)
    rsv = np.asarray(rsv, dtype=float)
    avir2 = np.asarray(avir2, dtype=float)
    vpv = np.asarray(vpv, dtype=float)
    avrr2 = np.asarray(avrr2, dtype=float)

    x = np.full_like(q, np.nan, dtype=float)
    a = np.full_like(q, np.nan, dtype=float)
    b = np.full_like(q, np.nan, dtype=float)
    g = np.full_like(q, np.nan, dtype=float)
    if not np.isfinite(rr) or not np.isfinite(bb):
        return {"x": x, "a": a, "b": b, "g": g}

    good_x = np.isfinite(rsv) & np.isfinite(q) & (np.abs(q) > EPS)
    x[good_x] = rsv[good_x] / q[good_x]

    good_a = (
        np.isfinite(x)
        & np.isfinite(avir2)
        & np.isfinite(vpv)
        & (np.abs(avir2 * vpv) > EPS)
    )
    a[good_a] = 4.0 * np.pi * np.pi * bb * bb * rr * rr * x[good_a] / (avir2[good_a] * vpv[good_a])

    good_b = np.isfinite(pprime) & np.isfinite(avir2) & (np.abs(avir2) > EPS)
    b[good_b] = 4.0 * np.pi * np.pi * MU0 * rr * rr * pprime[good_b] / avir2[good_b]

    good_g = np.isfinite(vpv) & np.isfinite(avrr2) & np.isfinite(x)
    g[good_g] = vpv[good_g] * avrr2[good_g] * x[good_g]
    return {"x": x, "a": a, "b": b, "g": g}


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot all EQIPQP variables into one big figure.")
    parser.add_argument("--prefix", default="eqdata_modelg5_rebuild_mdleqf9_qmap")
    parser.add_argument("--eqipqp-input-csv", default="eqgs1d_24_EQIPQP_INPUTS.csv")
    parser.add_argument("--rebuild-csv", default=None, help="Default: in/<prefix>_f_from_q_eqipqp.csv")
    parser.add_argument("--gfile", default="in/g260206.20000_teq_0114", help="Path to gfile.")
    parser.add_argument(
        "--gfile-profile-csv",
        default="in/g260206.20000_teq_0114_profiles_from_gfile.csv",
        help="Gfile profiles csv. If present, overlays q/p/F on comparable panels.",
    )
    parser.add_argument(
        "--gfile-f-sign",
        type=float,
        default=-1.0,
        help="Multiply gfile F by this factor for sign convention alignment.",
    )
    parser.add_argument(
        "--gfile-pprime-sign",
        type=float,
        default=1.0,
        help="Multiply gfile p' by this factor for sign convention alignment.",
    )
    parser.add_argument("--out", default=None, help="Default: in/<prefix>_eqipqp_all_variables_big.png")
    parser.add_argument(
        "--new-geom-csv",
        default=None,
        help="Optional external geometry csv to overlay RSV/AVIR2/VPV/AVRR2 on geometry panels.",
    )
    parser.add_argument(
        "--out-abg",
        default=None,
        help="Optional output path for A/B/G and deltaF diagnostics figure."
             " Default: in/<prefix>_gfile_rebuild_abg_deltaf.png",
    )
    parser.add_argument("--dpi", type=int, default=180)
    args = parser.parse_args()

    prefix = args.prefix
    rebuild_csv = Path(args.rebuild_csv or f"in/{prefix}_f_from_q_eqipqp.csv")
    eqipqp_csv = Path(args.eqipqp_input_csv)
    radial_csv = Path(f"{prefix}_radial_profiles.csv")
    profile1d_csv = Path(f"{prefix}_1D_profiles.csv")
    params_csv = Path(f"{prefix}_parameters.csv")
    out_png = Path(args.out or f"in/{prefix}_eqipqp_all_variables_big.png")
    out_abg_png = Path(args.out_abg or f"in/{prefix}_gfile_rebuild_abg_deltaf.png")
    out_png.parent.mkdir(parents=True, exist_ok=True)
    out_abg_png.parent.mkdir(parents=True, exist_ok=True)

    deq = pd.read_csv(eqipqp_csv)
    drb = pd.read_csv(rebuild_csv)
    drad = pd.read_csv(radial_csv)
    d1d = pd.read_csv(profile1d_csv)
    dpar = pd.read_csv(params_csv)

    x_eq = deq[_col(deq, ["PSIPV", "PSIP"], "PSIPV")].to_numpy(float)
    x_rb = drb[_col(drb, ["PSIP", "PSIPV"], "PSIP in rebuild")].to_numpy(float)
    x_rd = drad[_col(drad, ["PSIPV", "PSIP"], "PSIPV radial")].to_numpy(float)
    x_1d = d1d[_col(d1d, ["PSIPS", "PSIP", "PSIPV"], "PSIPS")].to_numpy(float)
    psip_max = float(np.nanmax(x_rd))
    g = read_gfile(str(args.gfile))
    gfile_data = _load_gfile_profiles(
        Path(args.gfile_profile_csv),
        g,
        psip_max=psip_max,
        f_sign=args.gfile_f_sign,
        pprime_sign=args.gfile_pprime_sign,
    )

    rr = float(dpar.loc[dpar["Parameter"] == "RR", "Value"].iloc[0]) if "RR" in dpar["Parameter"].values else np.nan
    bb = float(dpar.loc[dpar["Parameter"] == "BB", "Value"].iloc[0]) if "BB" in dpar["Parameter"].values else np.nan
    f_edge = 2.0 * np.pi * rr * bb if np.isfinite(rr) and np.isfinite(bb) else np.nan

    q_eq = deq[_col(deq, ["QPSI", "QPV", "QPS"], "Q for A/B/G")].to_numpy(float)
    dp_eq = deq[_col(deq, ["DPPSI", "DPPSI_PaPerWb"], "DPPSI for A/B/G")].to_numpy(float)
    rsv_eq = deq[_col(deq, ["RSV", "RSV_m"], "RSV for A/B/G")].to_numpy(float)
    avir2_eq = deq[_col(deq, ["AVIR2"], "AVIR2 for A/B/G")].to_numpy(float)
    vpv_eq = deq[_col(deq, ["VPV", "VPV_m3"], "VPV for A/B/G")].to_numpy(float)
    avrr2_eq = deq[_col(deq, ["AVRR2", "AVRR2_m2"], "AVRR2 for A/B/G")].to_numpy(float)
    abg_eq = _calc_abg(q_eq, dp_eq, rsv_eq, avir2_eq, vpv_eq, avrr2_eq, rr, bb)

    g_eq = None
    g_rb = None
    g_rd = None
    g_1d = None
    g_geom = None
    abg_g = None
    ext_geom = None
    deltaf_g = None
    g_f_rebuilt_eq = None
    g_f_rebuilt_rb = None
    g_f_rebuilt_rmse = np.nan
    fdfdpsin_rb = None
    fdfdpsin_gf = None
    if gfile_data is not None:
        psi_n_eq = np.clip(x_eq / max(psip_max, EPS), 0.0, 1.0)
        psi_n_rb = np.clip(x_rb / max(psip_max, EPS), 0.0, 1.0)
        psi_n_rd = np.clip(x_rd / max(psip_max, EPS), 0.0, 1.0)
        psi_n_1d = np.clip(x_1d / max(psip_max, EPS), 0.0, 1.0)

        g_eq = {
            "x": x_eq,
            "psi_n": psi_n_eq,
            "q": _safe_interp(gfile_data["psi_n"], gfile_data["q"], psi_n_eq),
            "p": _safe_interp(gfile_data["psi_n"], gfile_data["p"], psi_n_eq),
            "f": _safe_interp(gfile_data["psi_n"], gfile_data["f"], psi_n_eq),
            "pprime": _safe_interp(gfile_data["psi_n"], gfile_data["pprime"], psi_n_eq),
            "ffprime": _safe_interp(gfile_data["psi_n"], gfile_data["ffprime"], psi_n_eq),
        }
        g_rb = {
            "x": x_rb,
            "f": _safe_interp(gfile_data["x"], gfile_data["f"], x_rb),
        }
        g_rd = {
            "x": x_rd,
            "f": _safe_interp(gfile_data["x"], gfile_data["f"], x_rd),
            "q": _safe_interp(gfile_data["x"], gfile_data["q"], x_rd),
            "p": _safe_interp(gfile_data["x"], gfile_data["p"], x_rd),
        }
        g_1d = {
            "x": x_1d,
            "f": _safe_interp(gfile_data["x"], gfile_data["f"], x_1d),
            "q": _safe_interp(gfile_data["x"], gfile_data["q"], x_1d),
            "p": _safe_interp(gfile_data["x"], gfile_data["p"], x_1d),
        }
        g_geom = _compute_eqcalv_like_from_gfile(
            g,
            psi_n_levels=g_eq["psi_n"],
            psip_levels=g_eq["x"],
            q_levels=g_eq["q"],
            f_levels=g_eq["f"],
            rr=rr,
            bb=bb,
        )
        # q_eq = _safe_interp(gfile_data["x"], gfile_data["q"], x_eq)
        g_f_rebuilt_eq = _rebuild_f_eqipqp(
            x_eq,
            q_eq,
            g_eq["pprime"],
            g_geom["avir2"],
            g_geom["vpv"],
            g_geom["avrr2"],
            g_geom["rsv"],
            rr,
            bb,
        )
        g_f_rebuilt_eq = _rebuild_f_eqipqp(
            x_eq,
            q_eq,
            g_eq["pprime"],
            deq[_col(deq, ["AVIR2"], "AVIR2")],
            deq[_col(deq, ["VPV", "VPV_m3"], "VPV")],
            deq[_col(deq, ["AVRR2", "AVRR2_m2"], "AVRR2")],
            rsv_eq,
            rr,
            bb,
        )

        g_f_rebuilt_rb = _safe_interp(x_eq, g_f_rebuilt_eq, x_rb)
        v_gf = np.isfinite(g_f_rebuilt_eq) & np.isfinite(g_eq["f"])
        if np.count_nonzero(v_gf) > 0:
            g_f_rebuilt_rmse = float(np.sqrt(np.mean((g_f_rebuilt_eq[v_gf] - g_eq["f"][v_gf]) ** 2)))
        abg_g = _calc_abg(g_eq["q"], g_eq["pprime"], g_geom["rsv"], g_geom["avir2"], g_geom["vpv"], g_geom["avrr2"], rr, bb)
        deltaf_g = g_f_rebuilt_eq - g_eq["f"]

        # Compute F*dF/dpsi_n for comparison panel
        f_rb_on_eq = _safe_interp(x_rb, drb["F_rebuilt"].to_numpy(float), x_eq)
        psi_n_eq_arr = np.clip(x_eq / max(psip_max, EPS), 0.0, 1.0)
        valid_rb = np.isfinite(f_rb_on_eq) & np.isfinite(psi_n_eq_arr)
        if np.count_nonzero(valid_rb) >= 4:
            dfdpsin_rb = np.full_like(f_rb_on_eq, np.nan)
            dfdpsin_rb[valid_rb] = np.gradient(f_rb_on_eq[valid_rb], psi_n_eq_arr[valid_rb], edge_order=1)
            fdfdpsin_rb = f_rb_on_eq * dfdpsin_rb
        # F_gfile * dF_gfile/dpsi_n on eq grid
        f_gf_on_eq = g_eq["f"]
        valid_gf = np.isfinite(f_gf_on_eq) & np.isfinite(psi_n_eq_arr)
        if np.count_nonzero(valid_gf) >= 4:
            dfdpsin_gf = np.full_like(f_gf_on_eq, np.nan)
            dfdpsin_gf[valid_gf] = np.gradient(f_gf_on_eq[valid_gf], psi_n_eq_arr[valid_gf], edge_order=1)
            fdfdpsin_gf = f_gf_on_eq * dfdpsin_gf

    if args.new_geom_csv is not None:
        psi_n_eq = np.clip(x_eq / max(psip_max, EPS), 0.0, 1.0)
        ext_geom = _load_external_geom(Path(args.new_geom_csv), psi_n_eq)

    fig, ax = plt.subplots(4, 4, figsize=(24, 16), constrained_layout=True)
    axs = ax.ravel()

    # 1) F comparison
    # _plot_line(axs[0], x_rb, drb["F_ref"], "k-", lw=1.8, label="F_ref (TTV)", positive=True)
    _plot_line(axs[0], x_rb, drb["F_rebuilt"], "r--", lw=1.8, label="F_rebuilt", positive=True)
    if gfile_data is not None:
        _plot_line(axs[0], gfile_data["x"], gfile_data["f"], color="tab:blue", ls=":", lw=1.8, label="F_gfile", positive=True)
    # if g_f_rebuilt_rb is not None:
    # if 1:
    _plot_line(
            axs[0],
            x_rb,
            g_f_rebuilt_rb,
            color="tab:orange",
            ls="-.",
            lw=1.7,
            label="F_rebuilt_from_gfile",
            positive=True,
        )
    axs[0].set_title("F(psi) Comparison")
    axs[0].set_xlabel("PSIP")
    axs[0].set_ylabel("F")
    axs[0].set_yscale("log")
    axs[0].grid(alpha=0.3)
    axs[0].legend(fontsize=8)

    # 2) Error comparison
    axs[1].axhline(0.0, color="k", ls="--", lw=1.0)
    _plot_line(axs[1], x_rb, drb["dF"], color="tab:blue", lw=1.7, label="dF = rebuilt-ref")
    if g_rb is not None:
        dfg = g_rb["f"] - drb["F_ref"].to_numpy(float)
        _plot_line(axs[1], x_rb, dfg, color="tab:purple", ls=":", lw=1.5, label="dF_gfile = gfile-ref")
    axs[1].set_title("F Error")
    axs[1].set_xlabel("PSIP")
    axs[1].set_ylabel("dF")
    axs[1].grid(alpha=0.3)
    axs[1].legend(fontsize=8)

    # 3) q profiles
    axs[2].plot(x_eq, deq[_col(deq, ["QPSI", "QPV", "QPS"], "QPSI")], "m-", lw=1.6, label="QPSI(eqipqp)")
    axs[2].plot(x_rd, drad[_col(drad, ["QPV", "QPS"], "QPV")], "c--", lw=1.6, label="QPV(radial)")
    if gfile_data is not None:
        _plot_line(axs[2], gfile_data["x"], gfile_data["q"], color="k", ls=":", lw=1.7, label="q(gfile)")
    axs[2].set_title("q(psi)")
    axs[2].set_xlabel("PSIP")
    axs[2].set_ylabel("q")
    axs[2].grid(alpha=0.3)
    axs[2].legend(fontsize=8)

    # 4) DPPSI
    _plot_line(axs[3], x_eq, deq[_col(deq, ["DPPSI", "DPPSI_PaPerWb"], "DPPSI")], color="tab:red", lw=1.6, label="DPPSI(eq)")
    if g_eq is not None:
        _plot_line(axs[3], x_eq, g_eq["pprime"], color="k", ls=":", lw=1.5, label="DPPSI(gfile)")
    axs[3].set_title("DPPSI")
    axs[3].set_xlabel("PSIP")
    axs[3].set_ylabel("Pa/Wb")
    axs[3].grid(alpha=0.3)
    axs[3].legend(fontsize=8)

    # 5) RSV
    _plot_line(axs[4], x_eq, deq[_col(deq, ["RSV", "RSV_m"], "RSV")], color="tab:orange", lw=1.6, label="RSV(eq)")
    if g_geom is not None:
        _plot_line(axs[4], x_eq, g_geom["rsv"], color="k", ls=":", lw=1.5, label="RSV(gfile)")
    if ext_geom is not None:
        _plot_line(axs[4], x_eq, ext_geom["rsv"], color="tab:green", ls="--", lw=1.5, label="RSV(new geom)")
    axs[4].set_title("RSV")
    axs[4].set_xlabel("PSIP")
    axs[4].set_ylabel("m")
    axs[4].grid(alpha=0.3)
    axs[4].legend(fontsize=8)

    # 6) AVIR2
    _plot_line(axs[5], x_eq, deq[_col(deq, ["AVIR2"], "AVIR2")], color="tab:purple", lw=1.6, label="AVIR2(eq)")
    if g_geom is not None:
        _plot_line(axs[5], x_eq, g_geom["avir2"], color="k", ls=":", lw=1.5, label="AVIR2(gfile)")
    if ext_geom is not None:
        _plot_line(axs[5], x_eq, ext_geom["avir2"], color="tab:green", ls="--", lw=1.5, label="AVIR2(new geom)")
    axs[5].set_title("AVIR2")
    axs[5].set_xlabel("PSIP")
    axs[5].set_ylabel("arb.")
    axs[5].grid(alpha=0.3)
    axs[5].legend(fontsize=8)

    # 7) VPV
    _plot_line(axs[6], x_eq, deq[_col(deq, ["VPV", "VPV_m3"], "VPV")], color="tab:brown", lw=1.6, label="VPV(eq)")
    if g_geom is not None:
        _plot_line(axs[6], x_eq, g_geom["vpv"], color="k", ls=":", lw=1.5, label="VPV(gfile)")
    if ext_geom is not None:
        _plot_line(axs[6], x_eq, ext_geom["vpv"], color="tab:green", ls="--", lw=1.5, label="VPV(new geom)")
    axs[6].set_title("VPV")
    axs[6].set_xlabel("PSIP")
    axs[6].set_ylabel("m^3")
    axs[6].grid(alpha=0.3)
    axs[6].legend(fontsize=8)

    # 8) AVRR2
    _plot_line(axs[7], x_eq, deq[_col(deq, ["AVRR2", "AVRR2_m2"], "AVRR2")], color="tab:gray", lw=1.6, label="AVRR2(eq)")
    if g_geom is not None:
        _plot_line(axs[7], x_eq, g_geom["avrr2"], color="k", ls=":", lw=1.5, label="AVRR2(gfile)")
    if ext_geom is not None:
        _plot_line(axs[7], x_eq, ext_geom["avrr2"], color="tab:green", ls="--", lw=1.5, label="AVRR2(new geom)")
    axs[7].set_title("AVRR2")
    axs[7].set_xlabel("PSIP")
    axs[7].set_ylabel("m^2")
    axs[7].grid(alpha=0.3)
    axs[7].legend(fontsize=8)

    # 9) FIPV from eqipqp inputs vs radial TTV
    _plot_line(axs[8], x_eq, deq[_col(deq, ["FIPV"], "FIPV")], "r-", lw=1.6, label="FIPV(eqipqp csv)")
    _plot_line(axs[8], x_rd, drad[_col(drad, ["TTV", "TTS"], "TTV")], "k--", lw=1.6, label="TTV(radial)")
    if g_eq is not None:
        _plot_line(axs[8], x_eq, g_eq["f"], color="tab:blue", ls=":", lw=1.5, label="F(gfile)")
    axs[8].set_title("FIPV vs TTV")
    axs[8].set_xlabel("PSIP")
    axs[8].set_ylabel("F")
    axs[8].grid(alpha=0.3)
    axs[8].legend(fontsize=8)

    # 10) PSITV
    _plot_line(axs[9], x_rd, drad[_col(drad, ["PSITV", "TTS"], "PSITV")], color="tab:blue", lw=1.6, label="PSITV(eq)")
    if g_geom is not None:
        g_psitv_rd = _safe_interp(x_eq, g_geom["psitv"], x_rd)
        _plot_line(axs[9], x_rd, g_psitv_rd, color="k", ls=":", lw=1.5, label="PSITV(gfile)")
    axs[9].set_title("PSITV")
    axs[9].set_xlabel("PSIP")
    axs[9].set_ylabel("Wb")
    axs[9].grid(alpha=0.3)
    axs[9].legend(fontsize=8)

    # 11) PSIPNV
    _plot_line(axs[10], x_rd, drad[_col(drad, ["PSIPNV"], "PSIPNV")], color="tab:olive", lw=1.6, label="PSIPNV(eq)")
    if g_rd is not None:
        _plot_line(axs[10], x_rd, np.clip(x_rd / max(psip_max, EPS), 0.0, 1.0), color="k", ls=":", lw=1.5, label="PSIPNV(gfile)")
    axs[10].set_title("PSIPNV")
    axs[10].set_xlabel("PSIP")
    axs[10].set_ylabel("normalized")
    axs[10].grid(alpha=0.3)
    axs[10].legend(fontsize=8)

    # 12) p(psi)
    _plot_line(axs[11], x_1d, d1d[_col(d1d, ["PPPS"], "PPPS")], color="tab:cyan", lw=1.6, label="PPPS(eq)")
    if gfile_data is not None:
        _plot_line(axs[11], gfile_data["x"], gfile_data["p"], color="k", ls=":", lw=1.7, label="p(gfile)")
    axs[11].set_title("PPPS(psi)")
    axs[11].set_xlabel("PSIPS")
    axs[11].set_ylabel("Pa")
    axs[11].grid(alpha=0.3)
    axs[11].legend(fontsize=8)

    # 13) TTPS(psi)
    _plot_line(axs[12], x_1d, d1d[_col(d1d, ["TTPS"], "TTPS")], color="tab:pink", lw=1.6, label="TTPS(eq)")
    if g_1d is not None:
        _plot_line(axs[12], x_1d, g_1d["f"], color="k", ls=":", lw=1.5, label="TTPS(gfile)")
    axs[12].set_title("TTPS(psi)")
    axs[12].set_xlabel("PSIPS")
    axs[12].set_ylabel("F")
    axs[12].grid(alpha=0.3)
    axs[12].legend(fontsize=8)

    # 14) F*dF/dpsi_n comparison (replaces Direct Inputs panel)
    axs[13].axhline(0.0, color="k", ls="--", lw=1.0)
    if fdfdpsin_rb is not None:
        _plot_line(axs[13], x_eq, fdfdpsin_rb, color="tab:red", lw=1.7, label="F_rebuilt*dF/dpsi_n")
    if fdfdpsin_gf is not None:
        _plot_line(axs[13], x_eq, fdfdpsin_gf, color="tab:blue", ls=":", lw=1.7, label="F_gfile*dF/dpsi_n")
    axs[13].set_title("F*dF/dpsi_n")
    axs[13].set_xlabel("PSIP")
    axs[13].set_ylabel("F*dF/dpsi_n")
    axs[13].grid(alpha=0.3)
    axs[13].legend(fontsize=8)

    # 15) constants panel
    axs[14].axis("off")
    info = [
        f"prefix: {prefix}",
        f"RR = {rr:.6f}" if np.isfinite(rr) else "RR = NaN",
        f"BB = {bb:.6f}" if np.isfinite(bb) else "BB = NaN",
        f"2*pi*BB*RR = {f_edge:.6f}" if np.isfinite(f_edge) else "2*pi*BB*RR = NaN",
        f"eqipqp N = {len(deq)}",
        f"rebuild N = {len(drb)}",
        f"gfile csv = {args.gfile_profile_csv}" if gfile_data is not None else "gfile csv = (not loaded)",
        f"gfile path = {args.gfile}",
    ]
    if gfile_data is not None:
        info.append(f"gfile F sign = {args.gfile_f_sign:.1f}")
        info.append(f"gfile p' sign = {args.gfile_pprime_sign:.1f}")
    if args.new_geom_csv is not None:
        info.append(f"new geom csv = {args.new_geom_csv}")
    if g_geom is not None:
        n_geom = int(np.count_nonzero(np.isfinite(g_geom["rsv"])))
        info.append(f"gfile geom valid = {n_geom}/{len(g_geom['rsv'])}")
        if np.isfinite(g_f_rebuilt_rmse):
            info.append(f"RMSE(Freb_gf,F_gf) = {g_f_rebuilt_rmse:.3e}")
    else:
        info.append("gfile geom = failed")
    axs[14].text(0.02, 0.98, "\n".join(info), va="top", ha="left", fontsize=11, family="monospace")

    # 16) reserved for method string / quick stats
    axs[15].axis("off")
    method = drb["method"].iloc[0] if "method" in drb.columns and len(drb) > 0 else "unknown"
    valid = np.isfinite(drb["F_rebuilt"].to_numpy(float)) & np.isfinite(drb["F_ref"].to_numpy(float))
    if np.count_nonzero(valid) > 0:
        d = drb["F_rebuilt"].to_numpy(float) - drb["F_ref"].to_numpy(float)
        rmse = np.sqrt(np.mean((d[valid]) ** 2))
        mae = np.mean(np.abs(d[valid]))
        txt = f"method = {method}\nRMSE = {rmse:.6e}\nMAE  = {mae:.6e}\nNvalid = {np.count_nonzero(valid)}/{len(valid)}"
    else:
        txt = f"method = {method}\nNo valid F_rebuilt/F_ref points"
    axs[15].text(0.02, 0.98, txt, va="top", ha="left", fontsize=11, family="monospace")

    fig.suptitle("EQIPQP Variables Overview", fontsize=18, y=1.01)
    fig.savefig(out_png, dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"saved: {out_png}")

    if abg_g is not None and deltaf_g is not None:
        fig2, ax2 = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)
        aa = ax2.ravel()

        _plot_line(aa[0], x_eq, abg_eq["a"], color="tab:blue", lw=1.7, label="A(eq snapshot)")
        _plot_line(aa[0], x_eq, abg_g["a"], color="tab:orange", ls="--", lw=1.7, label="A(gfile rebuild)")
        aa[0].set_title("A Coefficient")
        aa[0].set_xlabel("PSIP")
        aa[0].set_ylabel("A")
        aa[0].grid(alpha=0.3)
        aa[0].legend(fontsize=9)

        _plot_line(aa[1], x_eq, abg_eq["b"], color="tab:blue", lw=1.7, label="B(eq snapshot)")
        _plot_line(aa[1], x_eq, abg_g["b"], color="tab:orange", ls="--", lw=1.7, label="B(gfile rebuild)")
        aa[1].set_title("B Coefficient")
        aa[1].set_xlabel("PSIP")
        aa[1].set_ylabel("B")
        aa[1].grid(alpha=0.3)
        aa[1].legend(fontsize=9)

        _plot_line(aa[2], x_eq, abg_eq["g"], color="tab:blue", lw=1.7, label="G(eq snapshot)")
        _plot_line(aa[2], x_eq, abg_g["g"], color="tab:orange", ls="--", lw=1.7, label="G(gfile rebuild)")
        aa[2].set_title("G Combination")
        aa[2].set_xlabel("PSIP")
        aa[2].set_ylabel("G")
        aa[2].grid(alpha=0.3)
        aa[2].legend(fontsize=9)

        aa[3].axhline(0.0, color="k", ls="--", lw=1.0)
        _plot_line(aa[3], x_eq, deltaf_g, color="tab:red", lw=1.8, label="deltaF = F_rebuilt_from_gfile - F_gfile")
        if "dF" in drb.columns:
            d_eq_on_eq = _safe_interp(x_rb, drb["dF"].to_numpy(float), x_eq)
            _plot_line(aa[3], x_eq, d_eq_on_eq, color="tab:purple", ls=":", lw=1.6, label="deltaF(eq rebuilt-ref)")
        aa[3].set_title("deltaF Diagnostics")
        aa[3].set_xlabel("PSIP")
        aa[3].set_ylabel("dF")
        aa[3].grid(alpha=0.3)
        aa[3].legend(fontsize=9)

        vdf = np.isfinite(deltaf_g)
        if np.count_nonzero(vdf) > 0:
            rmse_df = float(np.sqrt(np.mean((deltaf_g[vdf]) ** 2)))
            fig2.suptitle(f"Gfile Rebuild Diagnostics (RMSE deltaF={rmse_df:.3e})", fontsize=14)
        else:
            fig2.suptitle("Gfile Rebuild Diagnostics", fontsize=14)
        fig2.savefig(out_abg_png, dpi=args.dpi, bbox_inches="tight")
        plt.close(fig2)
        print(f"saved: {out_abg_png}")


if __name__ == "__main__":
    main()
