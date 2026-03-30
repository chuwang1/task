#!/usr/bin/env python3
"""Compare gfile-derived EQCALV geometry against input.profiles proxies,
and rebuild F(psi) via EQIPQP recursion using both geometry sources."""

import argparse
import csv
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from contourpy import contour_generator
from matplotlib.path import Path as MplPath

TWO_PI = 2.0 * np.pi
MU0 = 4.0e-7 * np.pi
EPS = 1.0e-14


# ---------------------------------------------------------------------------
# gfile reader
# ---------------------------------------------------------------------------
class GFile:
    pass


def read_gfile(path: Path) -> GFile:
    g = GFile()
    with path.open("r", encoding="utf-8") as f:
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
                raise ValueError("Failed to parse 5-value gfile row")
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
    nums = [float(x) for x in re.split(r"\s+", body) if x]

    ia = 0
    g.fpol = np.asarray(nums[ia : ia + g.nw], dtype=float);  ia += g.nw
    g.pres = np.asarray(nums[ia : ia + g.nw], dtype=float);  ia += g.nw
    g.ffprime = np.asarray(nums[ia : ia + g.nw], dtype=float); ia += g.nw
    g.pprime = np.asarray(nums[ia : ia + g.nw], dtype=float);  ia += g.nw
    g.psirz = np.asarray(nums[ia : ia + g.nw * g.nh], dtype=float).reshape(g.nw, g.nh)
    ia += g.nw * g.nh
    g.qpsi = np.asarray(nums[ia : ia + g.nw], dtype=float)

    g.R = np.linspace(g.rleft, g.rleft + g.rdim, g.nw)
    g.Z = np.linspace(g.zmid - 0.5 * g.zdim, g.zmid + 0.5 * g.zdim, g.nh)
    g.psi_n = np.linspace(0.0, 1.0, g.nw)
    return g


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def rho_from_qpsi(qpsi, psimag, psibdy):
    q = np.asarray(qpsi, dtype=float)
    n = len(q)
    out = np.zeros(n, dtype=float)
    if n < 2:
        return out
    dpsi = float(psibdy - psimag) / float(n - 1)
    for i in range(1, n):
        out[i] = out[i - 1] + 0.5 * (q[i - 1] + q[i]) * dpsi
    edge = out[-1]
    if abs(edge) < EPS:
        return np.zeros(n, dtype=float)
    return np.sqrt(np.abs(out / edge))


def read_csv_numeric(path: Path) -> dict[str, np.ndarray]:
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        cols = {k: [] for k in reader.fieldnames or []}
        for row in reader:
            for k, v in row.items():
                try:
                    cols[k].append(float(v))
                except Exception:
                    cols[k].append(np.nan)
    return {k: np.asarray(v, dtype=float) for k, v in cols.items()}


def interp_safe(x, y, x_new):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    good = np.isfinite(x) & np.isfinite(y)
    if np.count_nonzero(good) < 2:
        return np.full_like(x_new, np.nan, dtype=float)
    xg, yg = x[good], y[good]
    idx = np.argsort(xg); xg, yg = xg[idx], yg[idx]
    xu, inv = np.unique(xg, return_inverse=True)
    if len(xu) != len(xg):
        yg = np.bincount(inv, weights=yg) / np.maximum(np.bincount(inv), 1.0)
        xg = xu
    return np.interp(x_new, xg, yg)


def _interp_bilinear(field, r, z, rp, zp):
    ir = np.clip(np.searchsorted(r, rp, side="right") - 1, 0, len(r) - 2)
    iz = np.clip(np.searchsorted(z, zp, side="right") - 1, 0, len(z) - 2)
    r0, r1 = r[ir], r[ir + 1]
    z0, z1 = z[iz], z[iz + 1]
    tr = (rp - r0) / np.maximum(r1 - r0, EPS)
    tz = (zp - z0) / np.maximum(z1 - z0, EPS)
    return ((1 - tr) * (1 - tz) * field[ir, iz] + tr * (1 - tz) * field[ir + 1, iz]
            + (1 - tr) * tz * field[ir, iz + 1] + tr * tz * field[ir + 1, iz + 1])


def _integrate_psit_from_q(psip, q):
    n = len(psip)
    out = np.full(n, np.nan, dtype=float)
    if n == 0:
        return out
    out[0] = 0.0
    for i in range(1, n):
        if not (np.isfinite(psip[i]) and np.isfinite(psip[i - 1])
                and np.isfinite(q[i]) and np.isfinite(q[i - 1])):
            continue
        dq = q[i] + q[i - 1]
        qh = 0.5 * (q[i] + q[i - 1]) if abs(dq) < EPS else 2.0 * q[i] * q[i - 1] / dq
        out[i] = out[i - 1] + qh * (psip[i] - psip[i - 1])
    return out


def _extract_contour_segment(psi_total, r, z, level, axis_r, axis_z):
    cg = contour_generator(x=r, y=z, z=psi_total.T)
    segs = cg.lines(level)
    if not segs:
        return None
    best, best_score = None, -1.0
    for seg in segs:
        seg = np.asarray(seg, dtype=float)
        if seg.ndim != 2 or seg.shape[0] < 8:
            continue
        d = np.diff(seg, axis=0)
        plen = float(np.sum(np.hypot(d[:, 0], d[:, 1])))
        if plen <= 0.0:
            continue
        contains = MplPath(seg).contains_point((axis_r, axis_z))
        score = plen + (1e9 if contains else 0.0)
        if score > best_score:
            best_score, best = score, seg
    return best


# ---------------------------------------------------------------------------
# gfile contour-based geometry (EQCALV-like)
# ---------------------------------------------------------------------------
def compute_eqcalv_like_from_gfile(g, rr, bb):
    psi_n = np.asarray(g.psi_n, dtype=float)
    psip = psi_n * (TWO_PI * (g.psibdy - g.psimag))
    q = np.asarray(g.qpsi, dtype=float)
    f = TWO_PI * np.asarray(g.fpol, dtype=float)

    n = len(psi_n)
    out = {k: np.full(n, np.nan, dtype=float)
           for k in ("rsv", "avir2", "vpv", "avrr2")}
    if n < 2 or not np.isfinite(rr) or not np.isfinite(bb) or abs(bb) < EPS:
        return out

    psi_total = TWO_PI * (np.asarray(g.psirz, dtype=float) - float(g.psimag))
    dpsi_dr = np.gradient(psi_total, g.R, axis=0, edge_order=2)
    dpsi_dz = np.gradient(psi_total, g.Z, axis=1, edge_order=2)
    dpsi_total = TWO_PI * (float(g.psibdy) - float(g.psimag))
    if abs(dpsi_total) < EPS:
        return out

    psit = _integrate_psit_from_q(psip, q)
    out["rsv"] = np.sqrt(np.abs(psit / (np.pi * bb)))

    for i in range(n):
        if not (np.isfinite(psi_n[i]) and np.isfinite(out["rsv"][i])
                and np.isfinite(q[i]) and np.isfinite(f[i])):
            continue
        level = psi_n[i] * dpsi_total
        seg = _extract_contour_segment(psi_total, g.R, g.Z, level, g.rmaxis, g.zmaxis)
        if seg is None:
            continue
        if np.hypot(seg[0, 0] - seg[-1, 0], seg[0, 1] - seg[-1, 1]) > 1e-10:
            seg = np.vstack([seg, seg[0]])

        p0, p1 = seg[:-1], seg[1:]
        rm = 0.5 * (p0[:, 0] + p1[:, 0])
        dl = np.hypot(p1[:, 0] - p0[:, 0], p1[:, 1] - p0[:, 1])
        dpr = _interp_bilinear(dpsi_dr, g.R, g.Z, rm,
                               0.5 * (p0[:, 1] + p1[:, 1]))
        dpz = _interp_bilinear(dpsi_dz, g.R, g.Z, rm,
                               0.5 * (p0[:, 1] + p1[:, 1]))
        bpl = np.sqrt(dpr**2 + dpz**2) / np.maximum(TWO_PI * rm, EPS)

        ok = (np.isfinite(dl) & np.isfinite(rm) & np.isfinite(bpl)
              & (dl > 0) & (bpl > EPS) & (rm > EPS))
        if np.count_nonzero(ok) < 8:
            continue
        dlv, rv, bplv = dl[ok], rm[ok], bpl[ok]

        sumv = np.sum(dlv / bplv)
        sumavir2 = np.sum(dlv / (bplv * rv * rv))
        sumavrr2 = np.sum(dlv * bplv)
        if not np.isfinite(sumv) or sumv <= EPS:
            continue

        out["avir2"][i] = (sumavir2 / sumv) * rr * rr
        rsv_i = out["rsv"][i]
        ql = q[i]
        if abs(ql) > EPS and rsv_i > EPS:
            out["vpv"][i] = TWO_PI * rsv_i * bb * sumv / ql
            den = rsv_i**2 * bb**2 * sumv
            if abs(den) > EPS:
                out["avrr2"][i] = sumavrr2 * ql**2 / den

    if n >= 3:
        for k in ("avir2", "avrr2"):
            if np.isfinite(out[k][1]) and np.isfinite(out[k][2]):
                out[k][0] = (4.0 * out[k][1] - out[k][2]) / 3.0
        out["vpv"][0] = 0.0
        out["rsv"][0] = 0.0
    return out


# ---------------------------------------------------------------------------
# input.profiles proxy geometry
# ---------------------------------------------------------------------------
def compute_input_profiles_geometry(ag, ex, rr, bb):
    """Derive EQCALV-equivalent RSV/AVIR2/VPV/AVRR2 from input.profiles data.

    Better proxies than the first-pass raw dump:
      RSV  = rmin                           (geometric minor radius)
      AVIR2 = (R0/Rmaj)^2                   (zeroth-order <1/R^2>*R0^2)
      VPV  = 2*pi*rmin*B0/|q| * (dV/dpsip)  where dV/dpsip = EXPRO_volp*EXPRO_drdrho / (2*pi*|dpolflux/drho|)
      AVRR2 = q^2 * Bp0^2 / (rmin^2 * B0^2) (constant-Bp contour approx)
    """
    rho = ag["rho"]
    rmin = ag["rmin"]
    rmaj = ag["rmaj"]
    q_ip = ag["q"]
    polflux = ag["polflux"]  # Wb/rad, measured from axis

    volp = ex["EXPRO_volp"]       # dV/dr  (m^2)
    drdrho = ex["EXPRO_drdrho"]   # dr/drho (dimensionless)
    bp0 = ex["EXPRO_bp0"]        # B_pol at theta=0 (T)

    n = min(len(rho), len(volp), len(drdrho), len(bp0))
    rho = rho[:n]; rmin = rmin[:n]; rmaj = rmaj[:n]
    q_ip = q_ip[:n]; polflux = polflux[:n]
    volp = volp[:n]; drdrho = drdrho[:n]; bp0 = bp0[:n]

    # --- RSV = rmin ---
    rsv_ip = rmin.copy()

    # --- AVIR2 = (R0/Rmaj)^2 ---
    avir2_ip = (rr / np.maximum(np.abs(rmaj), EPS)) ** 2

    # --- VPV from dV/dpsip ---
    # psip_code = 2*pi * polflux (Wb)
    # dpsip/drho = 2*pi * d(polflux)/d(rho)
    dpolflux_drho = np.gradient(polflux, rho, edge_order=2)
    dpsip_drho = TWO_PI * dpolflux_drho   # Wb  per  rho-unit

    # dV/dpsip = (dV/dr * dr/drho) / (dpsip/drho)
    dV_dpsip = volp * drdrho / np.where(np.abs(dpsip_drho) > EPS, dpsip_drho, np.nan)

    # VPV = 2*pi * RSV * B0 * (dV/dpsip) / q
    vpv_ip = TWO_PI * rmin * bb * dV_dpsip / np.where(np.abs(q_ip) > EPS, q_ip, np.nan)
    vpv_ip = np.abs(vpv_ip)  # sign convention: VPV > 0
    vpv_ip[0] = 0.0

    # --- AVRR2 approx: constant-Bp around contour ---
    # AVRR2 ~ q^2 * Bp0^2 / (rmin^2 * B0^2)
    avrr2_ip = q_ip**2 * bp0**2 / np.where(
        (rmin**2 * bb**2) > EPS, rmin**2 * bb**2, np.nan)
    # axis extrapolation
    if n >= 3 and np.isfinite(avrr2_ip[1]) and np.isfinite(avrr2_ip[2]):
        avrr2_ip[0] = (4.0 * avrr2_ip[1] - avrr2_ip[2]) / 3.0

    return {
        "rsv": rsv_ip,
        "avir2": avir2_ip,
        "vpv": vpv_ip,
        "avrr2": avrr2_ip,
        "rho": rho,
    }


# ---------------------------------------------------------------------------
# EQIPQP F rebuild
# ---------------------------------------------------------------------------
def rebuild_f_eqipqp(psip, q, pprime, avir2, vpv, avrr2, rsv, rr, bb):
    """Backward-sweep EQIPQP: reconstruct F from edge to axis."""
    n = len(psip)
    f = np.full(n, np.nan, dtype=float)
    if n < 2:
        return f
    f[-1] = TWO_PI * bb * rr
    four_pi2 = 4.0 * np.pi ** 2

    for i in range(n - 1, 0, -1):
        if not np.isfinite(f[i]):
            continue
        pp, pm = psip[i], psip[i - 1]
        qp, qm = q[i], q[i - 1]
        dp_p, dp_m = pprime[i], pprime[i - 1]
        vals = np.array([qp, qm, dp_p, dp_m,
                         avir2[i], avir2[i - 1],
                         vpv[i], vpv[i - 1],
                         avrr2[i], avrr2[i - 1]])
        if not np.isfinite(vals).all():
            continue
        if abs(qp) < EPS or abs(qm) < EPS:
            continue

        xp = rsv[i] / qp
        xm = rsv[i - 1] / qm

        alp = four_pi2 * bb**2 * rr**2 * xp / max(avir2[i] * vpv[i], EPS)
        blp = four_pi2 * MU0 * rr**2 * dp_p / max(avir2[i], EPS)
        flp = vpv[i] * avrr2[i] * xp

        if i == 1:
            alm = four_pi2 * bb**2 * rr**2 * xp / max(avir2[i - 1] * vpv[i], EPS)
        else:
            alm = four_pi2 * bb**2 * rr**2 * xm / max(avir2[i - 1] * vpv[i - 1], EPS)
        blm = four_pi2 * MU0 * rr**2 * dp_m / max(avir2[i - 1], EPS)
        flm = vpv[i - 1] * avrr2[i - 1] * xm

        yp = 0.5 * f[i] ** 2
        ym = yp + 0.5 * (blp + blm) * (pp - pm) + 0.5 * (alp + alm) * (flp - flm)
        if not np.isfinite(ym) or ym <= 0.0:
            continue
        f[i - 1] = np.sqrt(2.0 * ym)
    return f


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gfile", default="in/g260206.20000_teq_0114")
    ap.add_argument("--all-profiles",
                    default="/Users/dengxiaoya/CFEDRSW/OMFIT_out/outputs_new/all_profiles.csv")
    ap.add_argument("--extra-profiles",
                    default="/Users/dengxiaoya/CFEDRSW/OMFIT_out/outputs_new/input_profiles_extra.csv")
    ap.add_argument("--out-png",
                    default="in/gfile_vs_inputprofiles_geom_and_F_compare.png")
    ap.add_argument("--out-csv",
                    default="in/gfile_vs_inputprofiles_geom_and_F_compare.csv")
    args = ap.parse_args()

    # ---- read gfile ----
    g = read_gfile(Path(args.gfile))
    rr = float(g.rcentr)           # R0  (m)
    bb = float(abs(g.bcentr))      # |B0|  (T)
    psi_n_g = g.psi_n              # 0..1, nw points
    psip_g = psi_n_g * (TWO_PI * (g.psibdy - g.psimag))  # code-unit psip
    q_g = np.asarray(g.qpsi, dtype=float)
    # pprime from gfile: dp/d(psi_Wb/rad).  Convert to dp/dpsip_code:
    #   dpsip_code = 2*pi * dpsi_Wb  =>  dp/dpsip = pprime_gfile / (2*pi)
    #   BUT in EQIPQP the psip grid *already* has the 2*pi factor baked in
    #   (psip = 2*pi * psi_Wb/rad), so dp/dpsip_code = pprime_gfile.
    #   The gfile stores pprime in Pa / (Wb/rad).
    #   Fortran DPPSI = dp/d(PSIPV) where PSIPV = 2*pi*(psi-psi0).
    #   pprime_Wb/rad * dpsi_Wb/rad = pprime_code * dpsip_code  =>
    #   pprime_code = pprime_Wb/rad / (2*pi).
    pprime_code = np.asarray(g.pprime, dtype=float) / TWO_PI
    f_gfile = TWO_PI * np.asarray(g.fpol, dtype=float)  # F = 2*pi*R*Bt

    # sign convention: use absolute q,F for recursion
    sign_q = np.sign(q_g[len(q_g) // 2])
    sign_f = np.sign(f_gfile[len(f_gfile) // 2])
    sign_pp = np.sign(pprime_code[len(pprime_code) // 4])
    q_g_abs = sign_q * q_g
    f_gfile_abs = sign_f * f_gfile
    pprime_abs = sign_pp * pprime_code

    rho_g = rho_from_qpsi(g.qpsi, g.psimag, g.psibdy)

    # ---- gfile contour geometry ----
    print("computing gfile contour geometry ...")
    geom_g = compute_eqcalv_like_from_gfile(g, rr=rr, bb=bb)

    # ---- read input.profiles ----
    ag = read_csv_numeric(Path(args.all_profiles))
    ex = read_csv_numeric(Path(args.extra_profiles))

    # ---- input.profiles proxy geometry (on rho grid) ----
    ip_geom = compute_input_profiles_geometry(ag, ex, rr, bb)
    rho_ip = ip_geom["rho"]

    # ---- map everything onto gfile psi_n grid (257 pts) ----
    # input.profiles geometry -> interpolate rho_ip -> rho_g -> psi_n_g
    rsv_ip_g  = interp_safe(rho_ip, ip_geom["rsv"],   rho_g)
    avir2_ip_g = interp_safe(rho_ip, ip_geom["avir2"], rho_g)
    vpv_ip_g  = interp_safe(rho_ip, ip_geom["vpv"],   rho_g)
    avrr2_ip_g = interp_safe(rho_ip, ip_geom["avrr2"], rho_g)

    # gfile geometry (already on psi_n_g grid)
    rsv_gf   = geom_g["rsv"]
    avir2_gf = geom_g["avir2"]
    vpv_gf   = geom_g["vpv"]
    avrr2_gf = geom_g["avrr2"]

    # ---- EQIPQP rebuild: gfile geometry ----
    print("rebuilding F with gfile geometry ...")
    f_rb_gf = rebuild_f_eqipqp(
        psip_g, q_g_abs, pprime_abs,
        avir2_gf, vpv_gf, avrr2_gf, rsv_gf, rr, bb)

    # ---- EQIPQP rebuild: input.profiles geometry ----
    print("rebuilding F with input.profiles geometry ...")
    f_rb_ip = rebuild_f_eqipqp(
        psip_g, q_g_abs, pprime_abs,
        avir2_ip_g, vpv_ip_g, avrr2_ip_g, rsv_ip_g, rr, bb)

    # ---- metrics ----
    def rmse(a, b):
        ok = np.isfinite(a) & np.isfinite(b)
        return float(np.sqrt(np.mean((a[ok] - b[ok]) ** 2))) if ok.any() else np.nan

    rmse_gf = rmse(f_rb_gf, f_gfile_abs)
    rmse_ip = rmse(f_rb_ip, f_gfile_abs)
    rmse_gf_ip = rmse(f_rb_gf, f_rb_ip)
    print(f"RMSE  F_rebuilt(gfile_geom)  vs F_gfile       = {rmse_gf:.6f}")
    print(f"RMSE  F_rebuilt(input_geom)  vs F_gfile       = {rmse_ip:.6f}")
    print(f"RMSE  F_rebuilt(gfile_geom)  vs F_rebuilt(ip)  = {rmse_gf_ip:.6f}")

    # ---- also map gfile geometry to rho_ip for comparison ----
    g_rsv_onip   = interp_safe(rho_g, rsv_gf,   rho_ip)
    g_avir2_onip = interp_safe(rho_g, avir2_gf, rho_ip)
    g_vpv_onip   = interp_safe(rho_g, vpv_gf,   rho_ip)
    g_avrr2_onip = interp_safe(rho_g, avrr2_gf, rho_ip)

    # ---- CSV output ----
    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as fout:
        w = csv.writer(fout)
        w.writerow(["psi_n", "rho_tor",
                     "F_gfile", "F_rebuilt_gfile_geom", "F_rebuilt_input_geom",
                     "RSV_gf", "RSV_ip", "AVIR2_gf", "AVIR2_ip",
                     "VPV_gf", "VPV_ip", "AVRR2_gf", "AVRR2_ip"])
        for i in range(len(psi_n_g)):
            w.writerow([f"{psi_n_g[i]:.10e}", f"{rho_g[i]:.10e}",
                         f"{f_gfile_abs[i]:.10e}",
                         f"{f_rb_gf[i]:.10e}", f"{f_rb_ip[i]:.10e}",
                         f"{rsv_gf[i]:.10e}", f"{rsv_ip_g[i]:.10e}",
                         f"{avir2_gf[i]:.10e}", f"{avir2_ip_g[i]:.10e}",
                         f"{vpv_gf[i]:.10e}", f"{vpv_ip_g[i]:.10e}",
                         f"{avrr2_gf[i]:.10e}", f"{avrr2_ip_g[i]:.10e}"])

    # ---- PLOT  (3 rows x 3 cols) ----
    fig, axes = plt.subplots(3, 3, figsize=(18, 14), constrained_layout=True)

    # --- row 0: geometry comparison (on rho_ip grid) ---
    ax = axes[0]
    ax[0].plot(rho_ip, g_rsv_onip, "k--", lw=1.8, label="RSV(gfile)")
    ax[0].plot(rho_ip, ip_geom["rsv"], color="tab:blue", lw=1.4, label="rmin(input)")
    ax[0].set_title("RSV"); ax[0].legend(fontsize=8)

    ax[1].plot(rho_ip, g_avir2_onip, "k--", lw=1.8, label="AVIR2(gfile)")
    ax[1].plot(rho_ip, ip_geom["avir2"], color="tab:orange", lw=1.4, label=r"$(R_0/R_{maj})^2$")
    ax[1].set_title("AVIR2"); ax[1].legend(fontsize=8)

    ax[2].plot(rho_ip, g_vpv_onip, "k--", lw=1.8, label="VPV(gfile)")
    ax[2].plot(rho_ip, ip_geom["vpv"], color="tab:green", lw=1.4, label="VPV(input)")
    ax[2].set_title("VPV"); ax[2].legend(fontsize=8)

    # --- row 1: AVRR2 + F comparison ---
    ax = axes[1]
    ax[0].plot(rho_ip, g_avrr2_onip, "k--", lw=1.8, label="AVRR2(gfile)")
    ax[0].plot(rho_ip, ip_geom["avrr2"], color="tab:red", lw=1.4,
               label=r"$q^2 B_{p0}^2 / (r_{min}^2 B_0^2)$")
    ax[0].set_title("AVRR2"); ax[0].legend(fontsize=8)

    ax[1].plot(psi_n_g, f_gfile_abs, "k-", lw=2, label="F(gfile)")
    ax[1].plot(psi_n_g, f_rb_gf, "b--", lw=1.6, label="F_rebuilt(gfile geom)")
    ax[1].plot(psi_n_g, f_rb_ip, "r:", lw=1.6, label="F_rebuilt(input geom)")
    ax[1].set_title(f"F comparison (vs psi_n)\n"
                    f"RMSE gf={rmse_gf:.3f}  ip={rmse_ip:.3f}")
    ax[1].legend(fontsize=8)

    ax[2].plot(rho_g, f_gfile_abs, "k-", lw=2, label="F(gfile)")
    ax[2].plot(rho_g, f_rb_gf, "b--", lw=1.6, label="F_rebuilt(gfile geom)")
    ax[2].plot(rho_g, f_rb_ip, "r:", lw=1.6, label="F_rebuilt(input geom)")
    ax[2].set_title("F comparison (vs rho)")
    ax[2].legend(fontsize=8)

    # --- row 2: residuals ---
    ax = axes[2]
    res_gf = f_rb_gf - f_gfile_abs
    res_ip = f_rb_ip - f_gfile_abs
    res_diff = f_rb_ip - f_rb_gf

    ok_gf = np.isfinite(res_gf)
    ok_ip = np.isfinite(res_ip)
    ok_d  = np.isfinite(res_diff)

    ax[0].plot(psi_n_g[ok_gf], res_gf[ok_gf], "b-", lw=1.4, label="gfile geom - F_gfile")
    ax[0].plot(psi_n_g[ok_ip], res_ip[ok_ip], "r-", lw=1.4, label="input geom - F_gfile")
    ax[0].axhline(0, color="gray", lw=0.5)
    ax[0].set_title("F residual vs psi_n"); ax[0].legend(fontsize=8)

    ax[1].plot(rho_g[ok_gf], res_gf[ok_gf], "b-", lw=1.4, label="gfile geom - F_gfile")
    ax[1].plot(rho_g[ok_ip], res_ip[ok_ip], "r-", lw=1.4, label="input geom - F_gfile")
    ax[1].axhline(0, color="gray", lw=0.5)
    ax[1].set_title("F residual vs rho"); ax[1].legend(fontsize=8)

    ax[2].plot(psi_n_g[ok_d], res_diff[ok_d], "m-", lw=1.4,
               label="F(input geom) - F(gfile geom)")
    ax[2].axhline(0, color="gray", lw=0.5)
    ax[2].set_title("F rebuild difference"); ax[2].legend(fontsize=8)

    for row in axes:
        for a in row:
            a.grid(True, alpha=0.3)
            a.set_xlabel("rho_tor" if "rho" in a.get_title().lower() else "psi_n")

    fig.suptitle("Gfile vs input.profiles geometry → F rebuild comparison", fontsize=13)
    out_png = Path(args.out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=170)
    plt.close(fig)

    print(f"saved figure: {out_png}")
    print(f"saved table : {out_csv}")


if __name__ == "__main__":
    main()
