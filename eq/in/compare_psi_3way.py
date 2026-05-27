#!/usr/bin/env python3
"""
Compare equilibrium metrics between a G-file and TASK/EQ outputs.

Primary comparison targets:
1. q(psi)
2. p(psi)
3. F(psi)
4. LCFS geometry
5. Magnetic axis

Supplementary plots of PSI(R,Z) are also generated.

Usage:
  python compare_psi_3way.py --prefix eqdata0114Test
  python compare_psi_3way.py --prefix eqdata0114Test --no-gfile
"""

import argparse
import glob
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

TWO_PI = 2.0 * np.pi
# Empirical sign alignment: gfile F(psi) uses opposite sign vs eqdata TTPS.
GFILE_F_SIGN = -1.0


def polygon_area(r, z):
    """Signed polygon area by shoelace formula."""
    if len(r) < 3:
        return 0.0
    return 0.5 * np.sum(r * np.roll(z, -1) - z * np.roll(r, -1))


def normalize_01(x):
    x = np.asarray(x, dtype=float)
    xmin = np.nanmin(x)
    xmax = np.nanmax(x)
    if not np.isfinite(xmin) or not np.isfinite(xmax) or abs(xmax - xmin) < 1e-12:
        return np.zeros_like(x)
    return (x - xmin) / (xmax - xmin)


def find_first_existing(paths):
    for p in paths:
        if os.path.exists(p):
            return p
    return None


def safe_gradient(y, x):
    """Numerically stable dy/dx for mostly monotonic x."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 2 or len(y) < 2:
        return np.zeros_like(y)

    order = np.argsort(x)
    xs = x[order]
    ys = y[order]
    # Collapse duplicate x points to avoid divide-by-zero in gradient.
    xu, idx = np.unique(xs, return_index=True)
    yu = ys[idx]
    if len(xu) < 2:
        return np.zeros_like(y)

    gu = np.gradient(yu, xu)
    gs = np.interp(xs, xu, gu)
    out = np.empty_like(gs)
    out[order] = gs
    return out


def bilinear_interp(x, y, f, xp, yp):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    f = np.asarray(f, dtype=float)
    if xp < x[0] or xp > x[-1] or yp < y[0] or yp > y[-1]:
        return np.nan
    ix = np.searchsorted(x, xp) - 1
    iy = np.searchsorted(y, yp) - 1
    ix = int(np.clip(ix, 0, len(x) - 2))
    iy = int(np.clip(iy, 0, len(y) - 2))
    x1, x2 = x[ix], x[ix + 1]
    y1, y2 = y[iy], y[iy + 1]
    q11 = f[iy, ix]
    q21 = f[iy, ix + 1]
    q12 = f[iy + 1, ix]
    q22 = f[iy + 1, ix + 1]
    tx = 0.0 if abs(x2 - x1) < 1e-14 else (xp - x1) / (x2 - x1)
    ty = 0.0 if abs(y2 - y1) < 1e-14 else (yp - y1) / (y2 - y1)
    return (
        (1 - tx) * (1 - ty) * q11
        + tx * (1 - ty) * q21
        + (1 - tx) * ty * q12
        + tx * ty * q22
    )


def integrate_psit(psip, q):
    psip = np.asarray(psip, dtype=float)
    q = np.asarray(q, dtype=float)
    psit = np.zeros_like(psip)
    for i in range(1, len(psip)):
        qeff = 2.0 * q[i] * q[i - 1] / max(q[i] + q[i - 1], 1e-30)
        psit[i] = psit[i - 1] + qeff * (psip[i] - psip[i - 1])
    return psit


def reconstruct_q_from_psit_psip(psit, psip):
    psit = np.asarray(psit, dtype=float)
    psip = np.asarray(psip, dtype=float)
    q = np.full_like(psit, np.nan)
    if len(psit) < 2:
        return q
    for i in range(len(psit)):
        if i == 0:
            dpsit = psit[1] - psit[0]
            dpsip = psip[1] - psip[0]
        elif i == len(psit) - 1:
            dpsit = psit[-1] - psit[-2]
            dpsip = psip[-1] - psip[-2]
        else:
            dpsit = psit[i + 1] - psit[i - 1]
            dpsip = psip[i + 1] - psip[i - 1]
        if abs(dpsip) > 1e-30:
            q[i] = dpsit / dpsip
    return q


def read_eqcalq_raw(path):
    lines = open(path).read().splitlines()
    out = {
        "rho": [], "psip": [], "psit": [], "dvdpsit": [], "dvdpsip": [],
        "qps": [], "tts": [], "sumavir2": [], "inside": []
    }
    for line in lines:
        if not line.strip() or line.startswith("#") or line.startswith("nr,"):
            continue
        parts = [p.strip() for p in line.split(",")]
        vals = [float(x) for x in parts[1].split()]
        out["rho"].append(vals[0])
        out["psip"].append(vals[1])
        out["psit"].append(vals[2])
        out["dvdpsit"].append(vals[3])
        out["dvdpsip"].append(vals[4])
        out["qps"].append(vals[5])
        out["tts"].append(vals[6])
        out["sumavir2"].append(vals[7])
        inside = int(parts[2]) if len(parts) > 2 else 1
        out["inside"].append(inside)
    return {k: np.asarray(v, dtype=float) for k, v in out.items()}


def compute_gfile_geometry_chain(g, psi_n_samples, axis_shift=None):
    if g is None:
        return None
    axis_shift = axis_shift or (0.0, 0.0)
    dR, dZ = axis_shift
    r = np.asarray(g["R"], dtype=float) + dR
    z = np.asarray(g["Z"], dtype=float) + dZ
    psi_rz = np.asarray(g["PSIRZ"], dtype=float)
    psipa = float(-np.nanmin(psi_rz))
    if not np.isfinite(psipa) or psipa <= 0:
        return None

    # [Z,R] layout in compare script
    dpsi_dz, dpsi_dr = np.gradient(psi_rz, z, r)

    sumavir2 = np.full_like(psi_n_samples, np.nan, dtype=float)
    dvdpsip = np.full_like(psi_n_samples, np.nan, dtype=float)
    for i, psi_n in enumerate(psi_n_samples):
        if psi_n <= 1e-6 or psi_n >= 1.0:
            continue
        level = psipa * (psi_n - 1.0)
        fig_dummy, ax = plt.subplots()
        cs = ax.contour(*np.meshgrid(r, z), psi_rz, levels=[level])
        segs = cs.allsegs[0] if cs.allsegs else []
        plt.close(fig_dummy)
        if not segs:
            continue
        # Prefer contour enclosing shifted magnetic axis.
        rax = float(g["axis"][0] + dR)
        zax = float(g["axis"][1] + dZ)
        valid = []
        for seg in segs:
            if np.min(seg[:, 0]) < rax < np.max(seg[:, 0]) and np.min(seg[:, 1]) < zax < np.max(seg[:, 1]):
                valid.append(seg)
        seg = max(valid if valid else segs, key=lambda s: len(s))
        if not np.allclose(seg[0], seg[-1]):
            seg = np.vstack([seg, seg[0]])
        s1 = 0.0
        s2 = 0.0
        for k in range(1, len(seg)):
            x1, y1 = seg[k - 1]
            x2, y2 = seg[k]
            dl = float(np.hypot(x2 - x1, y2 - y1))
            rm = 0.5 * (x1 + x2)
            zm = 0.5 * (y1 + y2)
            dpr = bilinear_interp(r, z, dpsi_dr, rm, zm)
            dpz = bilinear_interp(r, z, dpsi_dz, rm, zm)
            if not np.isfinite(dpr) or not np.isfinite(dpz) or rm <= 1e-30:
                continue
            bp = np.sqrt(dpr * dpr + dpz * dpz) / (2.0 * np.pi * rm)
            if bp <= 1e-30:
                continue
            s1 += dl / (bp * rm * rm)
            s2 += dl / bp
        sumavir2[i] = s1
        dvdpsip[i] = s2
    return {"psi_n": psi_n_samples, "sumavir2": sumavir2, "dvdpsip": dvdpsip, "psipa": psipa}


def interp_field_to_grid(x_src, y_src, f_src, x_tgt, y_tgt):
    x_src = np.asarray(x_src, dtype=float)
    y_src = np.asarray(y_src, dtype=float)
    f_src = np.asarray(f_src, dtype=float)
    x_tgt = np.asarray(x_tgt, dtype=float)
    y_tgt = np.asarray(y_tgt, dtype=float)
    out = np.full((len(y_tgt), len(x_tgt)), np.nan, dtype=float)
    for j, yy in enumerate(y_tgt):
        for i, xx in enumerate(x_tgt):
            out[j, i] = bilinear_interp(x_src, y_src, f_src, xx, yy)
    return out


def export_surface_point_psi_differences(out_csv, eq, g, axis_shift, sigma_levels, ntheta=72):
    g_axis = np.asarray(g["axis"], dtype=float) + np.asarray(axis_shift, dtype=float)
    lcfs_r = np.asarray(g["lcfs_r"], dtype=float) + float(axis_shift[0])
    lcfs_z = np.asarray(g["lcfs_z"], dtype=float) + float(axis_shift[1])

    # Resample shifted gfile LCFS uniformly by curve index.
    idx = np.linspace(0, len(lcfs_r) - 1, int(ntheta), dtype=int)
    rows = []
    for ip, ii in enumerate(idx):
        rb = lcfs_r[ii]
        zb = lcfs_z[ii]
        for sg in sigma_levels:
            r = g_axis[0] + sg * (rb - g_axis[0])
            z = g_axis[1] + sg * (zb - g_axis[1])
            psi_g = bilinear_interp(g["R"] + axis_shift[0], g["Z"] + axis_shift[1], g["PSIRZ"], r, z)
            psi_e = bilinear_interp(eq["R"], eq["Z"], eq["PSIRZ"], r, z)
            rows.append(
                {
                    "point_id": ip + 1,
                    "sigma_geom": float(sg),
                    "R": float(r),
                    "Z": float(z),
                    "psi_gfile_shifted": float(psi_g) if np.isfinite(psi_g) else np.nan,
                    "psi_eqdata": float(psi_e) if np.isfinite(psi_e) else np.nan,
                    "delta_psi": float(psi_e - psi_g) if np.isfinite(psi_g) and np.isfinite(psi_e) else np.nan,
                }
            )
    pd.DataFrame(rows).to_csv(out_csv, index=False)
    return out_csv


def reconstruct_rho_theta_boundary(rsu, zsu, axis_ref, center_ref, ntgmax):
    rsu = np.asarray(rsu, dtype=float)
    zsu = np.asarray(zsu, dtype=float)
    rax, zax = axis_ref
    rr0, zz0 = center_ref
    dx = rsu - rax
    dz = zsu - zax
    rho = np.hypot(dx, dz)
    theta = np.arctan2(dz, dx)
    theta = np.where(theta < 0.0, theta + 2.0 * np.pi, theta)
    m = rho > 1e-8
    theta = theta[m]
    rho = rho[m]
    if theta.size < 3:
        return None, None, None, None
    order = np.argsort(theta)
    theta = theta[order]
    rho = rho[order]
    # merge near-duplicate theta values like EQSET_RHOB_FROM_RSU
    th_u = [theta[0]]
    rh_u = [rho[0]]
    for t, r in zip(theta[1:], rho[1:]):
        if t - th_u[-1] > 1e-8:
            th_u.append(t)
            rh_u.append(r)
        else:
            rh_u[-1] = 0.5 * (rh_u[-1] + r)
    th_u = np.asarray(th_u)
    rh_u = np.asarray(rh_u)
    thp = np.concatenate(([th_u[-1] - 2.0 * np.pi], th_u, [th_u[0] + 2.0 * np.pi]))
    rhp = np.concatenate(([rh_u[-1]], rh_u, [rh_u[0]]))
    thg = np.linspace(0.0, 2.0 * np.pi, int(ntgmax) + 1)
    thg[-1] = 2.0 * np.pi
    rhg = np.interp(thg, thp, rhp)
    r_bnd = rr0 + rhg * np.cos(thg)
    z_bnd = zz0 + rhg * np.sin(thg)
    return thg, rhg, r_bnd, z_bnd


def load_gfile(current_dir, eq_dir):
    from geqdsk import geqdsk

    gfile_path = find_first_existing(
        [
            os.path.join(current_dir, "in/g260206.20000_teq_0114"),
            os.path.join(eq_dir, "in/g260206.20000_teq_0114"),
        ]
    )
    if gfile_path is None:
        return None

    g = geqdsk(gfile_path)
    r = np.linspace(g.rleft, g.rleft + g.rdim, g.nw)
    z = np.linspace(g.zmid - 0.5 * g.zdim, g.zmid + 0.5 * g.zdim, g.nh)

    # TASK/EQ internal PSIRZ convention from gfile values.
    psi_rz = TWO_PI * (g.psirz - g.psibdy)

    return {
        "path": gfile_path,
        "obj": g,
        "R": r,
        "Z": z,
        "PSIRZ": psi_rz,
        "axis": (float(g.rmaxis), float(g.zmaxis)),
        "lcfs_r": np.asarray(g.rbbbs, dtype=float),
        "lcfs_z": np.asarray(g.zbbbs, dtype=float),
        "psi_n_profile": np.linspace(0.0, 1.0, g.nw),
        "q": np.asarray(g.qpsi, dtype=float),
        "p": np.asarray(g.pres, dtype=float),
        # gfile fpol is F=R*Bt, while TASK/EQ TTPS is 2*pi*F.
        # Apply sign alignment for direct comparison with eqdata.
        "F": GFILE_F_SIGN * TWO_PI * np.asarray(g.fpol, dtype=float),
        "pprime_raw": np.asarray(getattr(g, "pprime", np.full(g.nw, np.nan)), dtype=float),
        "ffprime_raw": np.asarray(getattr(g, "ffprime", np.full(g.nw, np.nan)), dtype=float),
    }


def resolve_gfile_profiles_csv(current_dir, eq_dir, csv_arg):
    candidates = []
    if csv_arg:
        if os.path.isabs(csv_arg):
            candidates.append(csv_arg)
        else:
            candidates.append(os.path.join(current_dir, csv_arg))
            candidates.append(os.path.join(eq_dir, csv_arg))

    candidates.extend(
        [
            os.path.join(current_dir, "g260206.20000_teq_0114_profiles_from_gfile.csv"),
            os.path.join(eq_dir, "in", "g260206.20000_teq_0114_profiles_from_gfile.csv"),
            os.path.join(eq_dir, "g260206.20000_teq_0114_profiles_from_gfile.csv"),
        ]
    )

    for p in candidates:
        if p and os.path.exists(p):
            return p

    for d in [current_dir, os.path.join(eq_dir, "in"), eq_dir]:
        for p in sorted(glob.glob(os.path.join(d, "*_profiles_from_gfile.csv"))):
            if os.path.exists(p):
                return p
    return None


def load_gfile_profiles_csv(path):
    df = pd.read_csv(path)
    required = ["psi_n", "p_pa", "F_2pi", "q"]
    for c in required:
        if c not in df.columns:
            raise ValueError(f"Missing column '{c}' in {path}")

    out = {
        "path_profile_csv": path,
        "psi_n_profile": df["psi_n"].to_numpy(dtype=float),
        "q": df["q"].to_numpy(dtype=float),
        "p": df["p_pa"].to_numpy(dtype=float),
        "F": GFILE_F_SIGN * df["F_2pi"].to_numpy(dtype=float),
    }

    if "pprime_pa_per_wb" in df.columns:
        out["pprime_raw"] = df["pprime_pa_per_wb"].to_numpy(dtype=float)
    if "FFprime" in df.columns:
        out["ffprime_raw"] = df["FFprime"].to_numpy(dtype=float)
    return out


def resolve_eqdata_prefix_paths(eq_dir, current_dir, prefix_arg):
    prefix_candidates = [prefix_arg] if prefix_arg else ["eqdata0114", "eqdata"]
    for prefix in prefix_candidates:
        if not prefix:
            continue

        parent = {
            "psirz": os.path.join(eq_dir, f"{prefix}_PSIRZ.csv"),
            "rg": os.path.join(eq_dir, f"{prefix}_RG_grid.csv"),
            "zg": os.path.join(eq_dir, f"{prefix}_ZG_grid.csv"),
            "prof1d": os.path.join(eq_dir, f"{prefix}_1D_profiles.csv"),
            "radial": os.path.join(eq_dir, f"{prefix}_radial_profiles.csv"),
            "params": os.path.join(eq_dir, f"{prefix}_parameters.csv"),
        }
        if all(os.path.exists(v) for v in parent.values()):
            return prefix, parent

        local = {
            "psirz": os.path.join(current_dir, f"{prefix}_PSIRZ.csv"),
            "rg": os.path.join(current_dir, f"{prefix}_RG_grid.csv"),
            "zg": os.path.join(current_dir, f"{prefix}_ZG_grid.csv"),
            "prof1d": os.path.join(current_dir, f"{prefix}_1D_profiles.csv"),
            "radial": os.path.join(current_dir, f"{prefix}_radial_profiles.csv"),
            "params": os.path.join(current_dir, f"{prefix}_parameters.csv"),
        }
        if all(os.path.exists(v) for v in local.values()):
            return prefix, local

    return None, None


def pick_first_column(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None


def load_eq_q_profile_from_csv(path):
    df = pd.read_csv(path)

    # Prefer QPS from eqgs1d_04_QPS.csv for MODELG=5 comparisons.
    if "QPS" in df.columns and "PSIP" in df.columns:
        psi_n = normalize_01(df["PSIP"].to_numpy(dtype=float))
        q = df["QPS"].to_numpy(dtype=float)
        return psi_n, q, f"{path}:QPS"

    q_col = pick_first_column(df, ["QPSI", "QPV", "QPS"])
    if q_col is None:
        return None

    if "PSIPNV" in df.columns:
        psi_n = df["PSIPNV"].to_numpy(dtype=float)
    elif "PSIPV" in df.columns:
        psi_n = normalize_01(df["PSIPV"].to_numpy(dtype=float))
    elif "PSIP" in df.columns:
        psi_n = normalize_01(df["PSIP"].to_numpy(dtype=float))
    else:
        return None

    q = df[q_col].to_numpy(dtype=float)
    return psi_n, q, f"{path}:{q_col}"


def resolve_eq_q_profile(current_dir, eq_dir, csv_arg, prefix=None):
    candidates = []
    if csv_arg:
        if os.path.isabs(csv_arg):
            candidates.append(csv_arg)
        else:
            candidates.append(os.path.join(current_dir, csv_arg))
            candidates.append(os.path.join(eq_dir, csv_arg))
    else:
        candidates.extend(
            [
                os.path.join(current_dir, f"{prefix}_eqgs1d_04_QPS.csv") if prefix else "",
                os.path.join(eq_dir, f"eq/{prefix}_eqgs1d_04_QPS.csv") if prefix else "",
                os.path.join(current_dir, f"{prefix}_eqgs1d_24_EQIPQP_INPUTS.csv") if prefix else "",
                os.path.join(eq_dir, f"eq/{prefix}_eqgs1d_24_EQIPQP_INPUTS.csv") if prefix else "",
                os.path.join(current_dir, "eqgs1d_04_QPS.csv"),
                os.path.join(eq_dir, "eqgs1d_04_QPS.csv"),
                os.path.join(current_dir, "eqgs1d_24_EQIPQP_INPUTS.csv"),
                os.path.join(eq_dir, "eqgs1d_24_EQIPQP_INPUTS.csv"),
            ]
        )

    for p in candidates:
        if not p or not os.path.exists(p):
            continue
        try:
            out = load_eq_q_profile_from_csv(p)
            if out is None:
                continue
            psi_n, q, source = out
            q_span = float(np.nanmax(q) - np.nanmin(q))
            if np.isfinite(q_span) and q_span > 1e-10:
                return {"psi_n": psi_n, "q": q, "source": source}
        except Exception as e:
            print(f"Warning: failed to load eq q profile from {p}: {e}")

    return None


def resolve_eqcalq_raw(current_dir, eq_dir, prefix, csv_arg=None):
    candidates = []
    if csv_arg:
        if os.path.isabs(csv_arg):
            candidates.append(csv_arg)
        else:
            candidates.append(os.path.join(current_dir, csv_arg))
            candidates.append(os.path.join(eq_dir, csv_arg))
    if prefix:
        candidates.extend(
            [
                os.path.join(current_dir, f"{prefix}_eqcalq_dvdpsit_raw.csv"),
                os.path.join(eq_dir, f"{prefix}_eqcalq_dvdpsit_raw.csv"),
            ]
        )
    candidates.extend(
        [
            os.path.join(current_dir, "eqcalq_dvdpsit_raw.csv"),
            os.path.join(eq_dir, "eq", "eqcalq_dvdpsit_raw.csv"),
        ]
    )
    return find_first_existing(candidates)


def load_eqdata(paths, eq_q_override=None):
    r = pd.read_csv(paths["rg"])["RG"].to_numpy(dtype=float)
    z = pd.read_csv(paths["zg"])["ZG"].to_numpy(dtype=float)

    # CSV matrix is stored as [R, Z]. contourf expects [Z, R].
    psi_matrix = pd.read_csv(paths["psirz"], comment="#").to_numpy(dtype=float).T

    df1d = pd.read_csv(paths["prof1d"])
    dfr = pd.read_csv(paths["radial"])
    params = pd.read_csv(paths["params"])
    pmap = {k: v for k, v in zip(params["Parameter"], params["Value"])}

    psi_p = df1d["PSIPS"].to_numpy(dtype=float)
    psi_n_p = normalize_01(psi_p)
    psi_n_q = dfr["PSIPNV"].to_numpy(dtype=float)
    p = df1d["PPPS"].to_numpy(dtype=float)
    f = df1d["TTPS"].to_numpy(dtype=float)
    dp_dpsin = safe_gradient(p, psi_n_p)
    ffprime_n = f * safe_gradient(f, psi_n_p)
    q = dfr["QPV"].to_numpy(dtype=float)
    q_source = f"{os.path.basename(paths['radial'])}:QPV"

    if eq_q_override is not None:
        q = np.asarray(eq_q_override["q"], dtype=float)
        psi_n_q = np.asarray(eq_q_override["psi_n"], dtype=float)
        q_source = str(eq_q_override["source"])

    return {
        "R": r,
        "Z": z,
        "PSIRZ": psi_matrix,
        "psi_p": psi_p,
        "psi_n_p": psi_n_p,
        "psi_n_q": psi_n_q,
        "p": p,
        "F": f,
        "dp_dpsin": dp_dpsin,
        "ffprime_n": ffprime_n,
        "q": q,
        "q_source": q_source,
        "axis": (float(pmap["RAXIS"]), float(pmap["ZAXIS"])),
        "params": pmap,
    }


def save_eq_psirz_csv(eq, out_path):
    r2d, z2d = np.meshgrid(eq["R"], eq["Z"])
    df = pd.DataFrame(
        {
            "R": r2d.ravel(),
            "Z": z2d.ravel(),
            "PSI": np.asarray(eq["PSIRZ"], dtype=float).ravel(),
        }
    )
    df.to_csv(out_path, index=False)
    print(f"Saved eq PSI(R,Z) CSV: {out_path}")


def extract_lcfs_from_psirz(r, z, psi_rz, level=0.0):
    """Extract largest contour line for a given level."""
    r2d, z2d = np.meshgrid(r, z)
    fig, ax = plt.subplots()
    cs = ax.contour(r2d, z2d, psi_rz, levels=[level])
    segs = cs.allsegs[0] if cs.allsegs else []
    plt.close(fig)

    if not segs:
        return None, None

    best = max(segs, key=lambda s: abs(polygon_area(s[:, 0], s[:, 1])))
    return best[:, 0], best[:, 1]


def extract_all_contours_from_psirz(r, z, psi_rz, level=0.0):
    r2d, z2d = np.meshgrid(r, z)
    fig, ax = plt.subplots()
    cs = ax.contour(r2d, z2d, psi_rz, levels=[level])
    segs = cs.allsegs[0] if cs.allsegs else []
    plt.close(fig)
    out = []
    for seg in segs:
        seg = np.asarray(seg, dtype=float)
        closed = bool(len(seg) >= 3 and np.allclose(seg[0], seg[-1]))
        out.append({"r": seg[:, 0], "z": seg[:, 1], "closed": closed, "n": len(seg)})
    return out


def boundary_shape_metrics(r, z):
    r = np.asarray(r, dtype=float)
    z = np.asarray(z, dtype=float)
    rmin, rmax = np.min(r), np.max(r)
    zmin, zmax = np.min(z), np.max(z)
    rgeo = 0.5 * (rmax + rmin)
    a = 0.5 * (rmax - rmin)
    if a <= 0:
        return None

    i_top = int(np.argmax(z))
    i_bot = int(np.argmin(z))
    delta_u = (rgeo - r[i_top]) / a
    delta_l = (rgeo - r[i_bot]) / a

    return {
        "Rmin": rmin,
        "Rmax": rmax,
        "Zmin": zmin,
        "Zmax": zmax,
        "Rgeo": rgeo,
        "a": a,
        "kappa": (zmax - zmin) / (2.0 * a),
        "delta_u": delta_u,
        "delta_l": delta_l,
        "delta": 0.5 * (delta_u + delta_l),
        "area": abs(polygon_area(r, z)),
    }


def nearest_distances(pa, pb):
    # pa: (N,2), pb: (M,2)
    d2 = np.sum((pa[:, None, :] - pb[None, :, :]) ** 2, axis=2)
    return np.sqrt(np.min(d2, axis=1))


def lcfs_distance_metrics(r1, z1, r2, z2):
    p1 = np.column_stack([r1, z1])
    p2 = np.column_stack([r2, z2])
    d12 = nearest_distances(p1, p2)
    d21 = nearest_distances(p2, p1)
    return {
        "mean_g_to_e": float(np.mean(d12)),
        "max_g_to_e": float(np.max(d12)),
        "mean_e_to_g": float(np.mean(d21)),
        "max_e_to_g": float(np.max(d21)),
        "chamfer_mean": float(0.5 * (np.mean(d12) + np.mean(d21))),
        "chamfer_max": float(max(np.max(d12), np.max(d21))),
    }


def compare_profile(psi_ref, y_ref, psi_cmp, y_cmp, psi_max=None):
    psi_ref = np.asarray(psi_ref, dtype=float)
    y_ref = np.asarray(y_ref, dtype=float)
    psi_cmp = np.asarray(psi_cmp, dtype=float)
    y_cmp = np.asarray(y_cmp, dtype=float)

    idx_ref = np.argsort(psi_ref)
    idx_cmp = np.argsort(psi_cmp)
    psi_ref, y_ref = psi_ref[idx_ref], y_ref[idx_ref]
    psi_cmp, y_cmp = psi_cmp[idx_cmp], y_cmp[idx_cmp]

    lo = max(float(np.min(psi_ref)), float(np.min(psi_cmp)))
    hi = min(float(np.max(psi_ref)), float(np.max(psi_cmp)))
    if psi_max is not None:
        hi = min(hi, float(psi_max))
    if hi <= lo:
        return None

    psi_common = np.linspace(lo, hi, 300)
    ref_i = np.interp(psi_common, psi_ref, y_ref)
    cmp_i = np.interp(psi_common, psi_cmp, y_cmp)
    diff = cmp_i - ref_i

    rmse = float(np.sqrt(np.mean(diff ** 2)))
    mae = float(np.mean(np.abs(diff)))
    max_abs = float(np.max(np.abs(diff)))
    rel_rmse = rmse / max(float(np.max(np.abs(ref_i))), 1e-12)

    return {
        "psi": psi_common,
        "ref": ref_i,
        "cmp": cmp_i,
        "rmse": rmse,
        "mae": mae,
        "max_abs": max_abs,
        "rel_rmse": rel_rmse,
    }


def print_profile_metrics(name, m):
    if m is None:
        print(f"{name}: overlap not found")
        return
    print(
        f"{name}: RMSE={m['rmse']:.4e}, MAE={m['mae']:.4e}, "
        f"MAX={m['max_abs']:.4e}, rel_RMSE={m['rel_rmse']:.4e}"
    )


def plot_profile_pair(ax, title, m, y_label):
    if m is None:
        ax.text(0.5, 0.5, "no overlap", ha="center", va="center")
        ax.set_title(title)
        ax.grid(alpha=0.25)
        return

    ax.plot(m["psi"], m["ref"], "b-", lw=2, label="gfile")
    ax.plot(m["psi"], m["cmp"], "r--", lw=2, label="eqdata")
    ax.set_title(title)
    ax.set_xlabel("normalized psi")
    ax.set_ylabel(y_label)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)


def make_overview_figure(eq, g, metrics, lcfs_eq_r, lcfs_eq_z, lcfs_m, axis_d, axis_dist):
    fig, axes = plt.subplots(2, 3, figsize=(18, 10), constrained_layout=True)

    plot_profile_pair(axes[0, 0], "q(psi)", metrics["q"], "q")
    plot_profile_pair(axes[0, 1], "p(psi)", metrics["p"], "Pa")
    plot_profile_pair(axes[0, 2], "F(psi)=2pi*R*Bt", metrics["F"], "T*m")
    plot_profile_pair(axes[1, 0], "dp/dpsi_n", metrics["dp"], "Pa")
    plot_profile_pair(axes[1, 1], "F*dF/dpsi_n", metrics["ffn"], "arb")

    ax = axes[1, 2]
    if g is None:
        ax.text(0.5, 0.5, "G-file missing", ha="center", va="center")
        ax.set_title("LCFS / Axis")
        ax.grid(alpha=0.25)
    else:
        g_axis = np.asarray(g["axis"])
        e_axis = np.asarray(eq["axis"])
        dR_s = e_axis[0] - g_axis[0]
        dZ_s = e_axis[1] - g_axis[1]
        ax.plot(g["lcfs_r"] + dR_s, g["lcfs_z"] + dZ_s, "b-", lw=1.8,
                label="gfile LCFS (shifted)")
        if lcfs_eq_r is not None:
            ax.plot(lcfs_eq_r, lcfs_eq_z, "r--", lw=1.8, label="eqdata LCFS")
        ax.plot([e_axis[0]], [e_axis[1]], "ro", ms=5, label="eqdata axis")
        ax.plot([g_axis[0] + dR_s], [g_axis[1] + dZ_s], "b+", ms=9, mew=2,
                label="gfile axis (shifted)")
        shift_str = f"dR={dR_s:+.3f}, dZ={dZ_s:+.3f}"
        ax.set_title(f"LCFS (gfile shifted {shift_str})")
        ax.set_xlabel("R (m)")
        ax.set_ylabel("Z (m)")
        ax.set_aspect("equal")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)

        txt_lines = []
        if lcfs_m is not None:
            txt_lines.append(f"LCFS mean={lcfs_m['chamfer_mean']:.3e} m")
            txt_lines.append(f"LCFS max={lcfs_m['chamfer_max']:.3e} m")
            txt_lines.append(f"(after axis-shift)")
        if txt_lines:
            ax.text(
                0.02,
                0.02,
                "\n".join(txt_lines),
                transform=ax.transAxes,
                va="bottom",
                ha="left",
                fontsize=8,
                bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none"},
            )

    fig.suptitle("G-file CSV vs EQDATA: Profile + Geometry Overview", fontsize=15)
    out = "comparison_overview_big.png"
    fig.savefig(out, dpi=160)
    print(f"Saved: {out}")


def main():
    parser = argparse.ArgumentParser(
        description="Compare q(psi), p(psi), F(psi), LCFS and axis between gfile and eqdata."
    )
    parser.add_argument(
        "--prefix",
        type=str,
        default=None,
        help="eqdata CSV prefix, e.g. 'eqdata0114Test' for eqdata0114Test_*.csv",
    )
    parser.add_argument(
        "--no-gfile",
        action="store_true",
        help="Disable loading/comparing external gfile data.",
    )
    parser.add_argument(
        "--gfile-profile-csv",
        type=str,
        default=None,
        help="Path to *_profiles_from_gfile.csv; if omitted, auto-detect.",
    )
    parser.add_argument(
        "--eq-psi-csv",
        type=str,
        default=None,
        help="Output CSV path for eq PSI(R,Z); default: <prefix>_eq_psirz_rz.csv",
    )
    parser.add_argument(
        "--eq-q-csv",
        type=str,
        default=None,
        help=(
            "Optional q-profile csv for EQ side. "
            "Default auto-detect priority: eqgs1d_04_QPS.csv, then eqgs1d_24_EQIPQP_INPUTS.csv"
        ),
    )
    parser.add_argument(
        "--q-psi-max",
        type=float,
        default=0.98,
        help="Upper psi_n limit used for q-profile comparison (default: 0.98).",
    )
    parser.add_argument(
        "--eqcalq-raw-csv",
        type=str,
        default=None,
        help="Optional eqcalq raw chain CSV; default auto-detects prefix-specific then generic file.",
    )
    parser.add_argument(
        "--eq-psi-scale",
        type=float,
        default=1.907,
        help="Global multiplicative scale applied to eq PSI-chain quantities for test overlay.",
    )
    args = parser.parse_args()

    current_dir = os.getcwd()
    eq_dir = os.path.abspath(os.path.join(current_dir, ".."))
    if eq_dir not in sys.path:
        sys.path.append(eq_dir)

    gfile_dir = os.path.abspath(os.path.join(current_dir, "./"))
    if gfile_dir not in sys.path:
        sys.path.append(gfile_dir)

    print(f"Current Dir: {current_dir}")
    print(f"EQ Dir: {eq_dir}")
    print(f"GFile Dir: {gfile_dir}")

    # Load gfile geometry data
    g = None
    if args.no_gfile:
        print("Skip G-file comparison (--no-gfile).")
    else:
        try:
            g = load_gfile(current_dir, eq_dir)
            if g is None:
                print("G-file not found")
            else:
                print(f"Loaded G-file: {g['path']}")
        except Exception as e:
            print(f"Error loading G-file: {e}")
            g = None

    # Load gfile profile CSV (preferred for profile comparison).
    g_profile = None
    if not args.no_gfile:
        profile_csv = resolve_gfile_profiles_csv(current_dir, eq_dir, args.gfile_profile_csv)
        if profile_csv is not None:
            try:
                g_profile = load_gfile_profiles_csv(profile_csv)
                print(f"Loaded G-file profile CSV: {profile_csv}")
                if g is not None:
                    # Keep geometry from raw gfile, but override profiles by CSV.
                    g.update(g_profile)
            except Exception as e:
                print(f"Error loading gfile profile CSV: {e}")

    # Fallback: use raw gfile profiles if CSV not available.
    if g_profile is None and g is not None:
        g_profile = {
            "psi_n_profile": g["psi_n_profile"],
            "q": g["q"],
            "p": g["p"],
            "F": g["F"],
            "pprime_raw": g.get("pprime_raw"),
            "ffprime_raw": g.get("ffprime_raw"),
        }

    if g_profile is not None:
        g_profile["dp_dpsin"] = safe_gradient(g_profile["p"], g_profile["psi_n_profile"])
        g_profile["ffprime_n"] = g_profile["F"] * safe_gradient(
            g_profile["F"], g_profile["psi_n_profile"]
        )

    # Load eqdata data first so prefix-scoped supplementary files can be found.
    prefix, paths = resolve_eqdata_prefix_paths(eq_dir, current_dir, args.prefix)
    if paths is None:
        print("eqdata CSV files not found")
        return

    # Load eqgs2d PSI (supplementary)
    psi_gs2d = None
    r_gs2d = None
    z_gs2d = None
    try:
        from read_eqgs2d_csv import load_psirz_grid

        csv_path = find_first_existing(
            [
                os.path.join(eq_dir, f"{prefix}_eqgs2d_01_PSIRZ_grid.csv") if prefix else "",
                os.path.join(current_dir, f"{prefix}_eqgs2d_01_PSIRZ_grid.csv") if prefix else "",
                os.path.join(eq_dir, "eqgs2d_01_PSIRZ_grid.csv"),
                os.path.join(current_dir, "eqgs2d_01_PSIRZ_grid.csv"),
            ]
        )
        if csv_path:
            r_gs2d, z_gs2d, psi_gs2d = load_psirz_grid(csv_path)
            print(f"Loaded eqgs2d CSV: {csv_path}")
    except Exception as e:
        print(f"Error loading eqgs2d CSV: {e}")

    eq_q_override = resolve_eq_q_profile(current_dir, eq_dir, args.eq_q_csv, prefix=prefix)
    eq = load_eqdata(paths, eq_q_override=eq_q_override)
    print(f"Loaded eqdata prefix: {prefix},{paths}")
    print(f"EQ q source: {eq['q_source']}")

    eq_psi_csv = args.eq_psi_csv or os.path.join(current_dir, f"{prefix}_eq_psirz_rz.csv")
    save_eq_psirz_csv(eq, eq_psi_csv)

    metrics = {"q": None, "p": None, "F": None, "dp": None, "ffn": None}
    # ------------------------------
    # 1) Profile comparisons
    # ------------------------------
    if g_profile is not None:
        q_m = compare_profile(
            g_profile["psi_n_profile"],
            g_profile["q"],
            eq["psi_n_q"],
            eq["q"],
            psi_max=args.q_psi_max,
        )
        p_m = compare_profile(g_profile["psi_n_profile"], g_profile["p"], eq["psi_n_p"], eq["p"])
        f_m = compare_profile(g_profile["psi_n_profile"], g_profile["F"], eq["psi_n_p"], eq["F"])
        dp_m = compare_profile(
            g_profile["psi_n_profile"], g_profile["dp_dpsin"], eq["psi_n_p"], eq["dp_dpsin"]
        )
        ffn_m = compare_profile(
            g_profile["psi_n_profile"], g_profile["ffprime_n"], eq["psi_n_p"], eq["ffprime_n"]
        )
        metrics.update({"q": q_m, "p": p_m, "F": f_m, "dp": dp_m, "ffn": ffn_m})

        print("\n=== Profile Metrics (gfile ref vs eqdata) ===")
        print_profile_metrics("q(psi)", q_m)
        print_profile_metrics("p(psi)", p_m)
        print_profile_metrics("F(psi)", f_m)
        print_profile_metrics("dp/dpsi_n", dp_m)
        print_profile_metrics("F*dF/dpsi_n", ffn_m)

        fig_pf, axes_pf = plt.subplots(1, 5, figsize=(24, 4.8), constrained_layout=True)
        plot_items = [
            ("q(psi)", q_m, "q"),
            ("p(psi) [Pa]", p_m, "Pa"),
            ("F(psi)=2pi*R*Bt", f_m, "T*m"),
            ("dp/dpsi_n", dp_m, "Pa"),
            ("F*dF/dpsi_n", ffn_m, "arb"),
        ]
        for ax, (title, m, ylab) in zip(axes_pf, plot_items):
            plot_profile_pair(ax, title, m, ylab)

        profile_plot = "profile_q_p_f_dp_ff_comparison.png"
        fig_pf.savefig(profile_plot, dpi=150)
        print(f"Saved: {profile_plot}")
        # Backward-compatible filename used by earlier workflow.
        fig_pf.savefig("profile_q_p_f_comparison.png", dpi=150)

    # ------------------------------
    # 2) LCFS and magnetic axis
    # ------------------------------
    lcfs_eq_r, lcfs_eq_z = extract_lcfs_from_psirz(eq["R"], eq["Z"], eq["PSIRZ"], level=0.0)

    lcfs_m = None
    axis_d = None
    axis_dist = None
    if g is not None and lcfs_eq_r is not None and "lcfs_r" in g and "axis" in g:
        g_axis = np.asarray(g["axis"])
        e_axis = np.asarray(eq["axis"])
        axis_d = e_axis - g_axis
        axis_dist = float(np.hypot(axis_d[0], axis_d[1]))

        # Shift gfile LCFS to align axes for shape comparison
        dR_s = axis_d[0]
        dZ_s = axis_d[1]
        g_lcfs_r_shifted = g["lcfs_r"] + dR_s
        g_lcfs_z_shifted = g["lcfs_z"] + dZ_s

        # Compute LCFS distance metrics on shifted coordinates
        lcfs_m = lcfs_distance_metrics(
            g_lcfs_r_shifted, g_lcfs_z_shifted, lcfs_eq_r, lcfs_eq_z
        )
        g_shape = boundary_shape_metrics(g_lcfs_r_shifted, g_lcfs_z_shifted)
        e_shape = boundary_shape_metrics(lcfs_eq_r, lcfs_eq_z)

        print("\n=== LCFS / Axis Metrics ===")
        print(
            f"Axis: dR={axis_d[0]:.4e} m, dZ={axis_d[1]:.4e} m, "
            f"distance={axis_dist:.4e} m"
        )
        print(
            f"LCFS Chamfer (after axis-shift): "
            f"mean={lcfs_m['chamfer_mean']:.4e} m, "
            f"max={lcfs_m['chamfer_max']:.4e} m"
        )
        if g_shape and e_shape:
            print(
                "LCFS shape delta (after axis-shift): "
                f"da={e_shape['a']-g_shape['a']:.4e} m, "
                f"dkappa={e_shape['kappa']-g_shape['kappa']:.4e}, "
                f"ddelta={e_shape['delta']-g_shape['delta']:.4e}, "
                f"dArea={e_shape['area']-g_shape['area']:.4e} m^2"
            )

        shift_str = f"dR={dR_s:+.3f}, dZ={dZ_s:+.3f} m"
        fig_lcfs, ax_lcfs = plt.subplots(figsize=(7, 7))
        ax_lcfs.plot(g_lcfs_r_shifted, g_lcfs_z_shifted, "b-", lw=2,
                     label="gfile LCFS (shifted)")
        ax_lcfs.plot(lcfs_eq_r, lcfs_eq_z, "r--", lw=2, label="eqdata LCFS")
        ax_lcfs.plot([e_axis[0]], [e_axis[1]], "ro", ms=6, label="eqdata axis")
        ax_lcfs.plot([g_axis[0] + dR_s], [g_axis[1] + dZ_s], "b+", ms=10, mew=2,
                     label="gfile axis (shifted)")
        ax_lcfs.set_xlabel("R (m)")
        ax_lcfs.set_ylabel("Z (m)")
        ax_lcfs.set_aspect("equal")
        ax_lcfs.set_title(f"LCFS Comparison (gfile shifted {shift_str})")
        ax_lcfs.grid(alpha=0.25)
        ax_lcfs.legend()
        lcfs_plot = "lcfs_axis_comparison.png"
        fig_lcfs.savefig(lcfs_plot, dpi=150)
        print(f"Saved: {lcfs_plot}")

        # Diagnostic: compare original shifted RSU/ZSU to all psi=0 contours from EQTORZ/PSIRZ.
        eq_zero_contours = extract_all_contours_from_psirz(eq["R"], eq["Z"], eq["PSIRZ"], level=0.0)
        ntg_overlay = int(float(eq["params"].get("NTGMAX", 64)))
        nsg_overlay = int(float(eq["params"].get("NSGMAX", 64)))
        rr_eq = float(eq["params"].get("RR", e_axis[0]))
        thg, rhg, r_bnd_rr, z_bnd_rr = reconstruct_rho_theta_boundary(
            g_lcfs_r_shifted, g_lcfs_z_shifted,
            axis_ref=(e_axis[0], e_axis[1]),
            center_ref=(rr_eq, 0.0),
            ntgmax=ntg_overlay,
        )
        _, _, r_bnd_axis, z_bnd_axis = reconstruct_rho_theta_boundary(
            g_lcfs_r_shifted, g_lcfs_z_shifted,
            axis_ref=(e_axis[0], e_axis[1]),
            center_ref=(e_axis[0], e_axis[1]),
            ntgmax=ntg_overlay,
        )
        fig_map, ax_map = plt.subplots(figsize=(7.5, 7.5))
        ax_map.plot(g_lcfs_r_shifted, g_lcfs_z_shifted, "b-", lw=2.2, label="gfile RSU/ZSU (shifted)")
        th_overlay = np.linspace(0.0, 2.0 * np.pi, ntg_overlay + 1)
        if thg is not None and rhg is not None:
            rh_overlay = np.interp(th_overlay, thg, rhg)
            # sigma=const curves
            for s in np.linspace(0.1, 1.0, 10):
                ax_map.plot(
                    rr_eq + s * rh_overlay * np.cos(th_overlay),
                    0.0 + s * rh_overlay * np.sin(th_overlay),
                    color="tab:orange",
                    lw=0.7 if s < 1.0 else 1.2,
                    alpha=0.5,
                )
            # theta=const rays using actual NTGMAX count
            for th, rho_l in zip(th_overlay[:-1], rh_overlay[:-1]):
                ss = np.linspace(0.0, 1.0, nsg_overlay)
                ax_map.plot(
                    rr_eq + ss * rho_l * np.cos(th),
                    0.0 + ss * rho_l * np.sin(th),
                    color="0.4",
                    lw=0.35,
                    alpha=0.25,
                )
        if eq_zero_contours:
            for i, seg in enumerate(eq_zero_contours):
                label = None
                if i == 0:
                    label = "eq psi=0 contour(s)"
                ls = "--" if seg["closed"] else ":"
                ax_map.plot(seg["r"], seg["z"], color="tab:red", lw=1.6, ls=ls, label=label)
        ax_map.plot([e_axis[0]], [e_axis[1]], "ro", ms=6, label="eq axis")
        ax_map.plot([g_axis[0] + dR_s], [g_axis[1] + dZ_s], "b+", ms=10, mew=2, label="gfile axis (shifted)")
        ax_map.set_xlabel("R (m)")
        ax_map.set_ylabel("Z (m)")
        ax_map.set_aspect("equal")
        ax_map.grid(alpha=0.25)
        ax_map.set_title("RSU/ZSU vs EQTORZ implied psi=0 contour")
        n_open = sum(0 if seg["closed"] else 1 for seg in eq_zero_contours)
        n_tot = len(eq_zero_contours)
        ax_map.text(
            0.02,
            0.98,
            f"eq psi=0 segments: {n_tot}\nopen segments: {n_open}\nNSGMAX={nsg_overlay}, NTGMAX={ntg_overlay}",
            transform=ax_map.transAxes,
            va="top",
            ha="left",
            fontsize=8,
            bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "none"},
        )
        ax_map.legend(fontsize=8)
        lcfs_map_plot = "lcfs_rsu_vs_eqtorz_map.png"
        fig_map.savefig(lcfs_map_plot, dpi=150)
        print(f"Saved: {lcfs_map_plot}")

        # Zoomed view around the midplane for easier inspection.
        fig_map_zoom, ax_map_zoom = plt.subplots(figsize=(8, 4.5))
        ax_map_zoom.plot(g_lcfs_r_shifted, g_lcfs_z_shifted, "b-", lw=2.2, label="gfile RSU/ZSU (shifted)")
        if thg is not None and rhg is not None:
            rh_overlay = np.interp(th_overlay, thg, rhg)
            for s in np.linspace(0.1, 1.0, 10):
                ax_map_zoom.plot(
                    rr_eq + s * rh_overlay * np.cos(th_overlay),
                    0.0 + s * rh_overlay * np.sin(th_overlay),
                    color="tab:orange",
                    lw=0.7 if s < 1.0 else 1.1,
                    alpha=0.45,
                )
            for th, rho_l in zip(th_overlay[:-1], rh_overlay[:-1]):
                ss = np.linspace(0.0, 1.0, nsg_overlay)
                ax_map_zoom.plot(
                    rr_eq + ss * rho_l * np.cos(th),
                    0.0 + ss * rho_l * np.sin(th),
                    color="0.4",
                    lw=0.3,
                    alpha=0.2,
                )
        if eq_zero_contours:
            for i, seg in enumerate(eq_zero_contours):
                ax_map_zoom.plot(seg["r"], seg["z"], color="tab:red", ls="-.", lw=1.4,
                                 label="eq psi=0 contour" if i == 0 else None)
        ax_map_zoom.plot([e_axis[0]], [e_axis[1]], "ro", ms=5, label="eq axis")
        ax_map_zoom.set_xlim(rr_eq - 3.0, rr_eq + 3.0)
        ax_map_zoom.set_ylim(-0.5, 0.5)
        ax_map_zoom.set_xlabel("R (m)")
        ax_map_zoom.set_ylabel("Z (m)")
        ax_map_zoom.set_title("LCFS / psi=0 / sigma-theta grid zoom (-0.5<Z<0.5)")
        ax_map_zoom.grid(alpha=0.25)
        ax_map_zoom.legend(fontsize=8, loc="upper right")
        lcfs_map_zoom_plot = "lcfs_rsu_vs_eqtorz_map_zoom_midplane.png"
        fig_map_zoom.savefig(lcfs_map_zoom_plot, dpi=170)
        print(f"Saved: {lcfs_map_zoom_plot}")

        # Full-resolution local zoom near the outboard midplane boundary.
        fig_map_full, ax_map_full = plt.subplots(figsize=(8, 4.5))
        ax_map_full.plot(g_lcfs_r_shifted, g_lcfs_z_shifted, "b-", lw=2.0, label="gfile RSU/ZSU (shifted)")
        if thg is not None and rhg is not None:
            # actual sigma count
            for s in np.linspace(1.0 / max(nsg_overlay, 1), 1.0, nsg_overlay):
                ax_map_full.plot(
                    rr_eq + s * rh_overlay * np.cos(th_overlay),
                    0.0 + s * rh_overlay * np.sin(th_overlay),
                    color="tab:orange",
                    lw=0.35 if s < 1.0 else 0.9,
                    alpha=0.35,
                )
            # actual theta count
            for th, rho_l in zip(th_overlay[:-1], rh_overlay[:-1]):
                ss = np.linspace(0.0, 1.0, nsg_overlay)
                ax_map_full.plot(
                    rr_eq + ss * rho_l * np.cos(th),
                    0.0 + ss * rho_l * np.sin(th),
                    color="0.25",
                    lw=0.22,
                    alpha=0.18,
                )
        if eq_zero_contours:
            for i, seg in enumerate(eq_zero_contours):
                ax_map_full.plot(seg["r"], seg["z"], color="tab:red", ls="-.", lw=1.4,
                                 label="eq psi=0 contour" if i == 0 else None)
        ax_map_full.plot([e_axis[0]], [e_axis[1]], "ro", ms=4, label="eq axis")
        # ax_map_full.set_xlim(10.0, max(np.nanmax(g_lcfs_r_shifted), np.nanmax(lcfs_eq_r) if lcfs_eq_r is not None else 10.5) + 0.2)
        ax_map_full.set_ylim(-6, -4)
        ax_map_full.set_xlabel("R (m)")
        ax_map_full.set_ylabel("Z (m)")
        ax_map_full.set_title(f"Outboard midplane full grid (NSGMAX={nsg_overlay}, NTGMAX={ntg_overlay})")
        ax_map_full.grid(alpha=0.25)
        ax_map_full.legend(fontsize=8, loc="upper right")
        lcfs_map_full_plot = "lcfs_rsu_vs_eqtorz_map_zoom_outboard_fullgrid.png"
        fig_map_full.savefig(lcfs_map_full_plot, dpi=200)
        print(f"Saved: {lcfs_map_full_plot}")

        # Reconstruct RHOG(theta) boundary from shifted RSU/ZSU and show how
        # using different centers changes the mapped boundary.
        fig_rhog, ax_rhog = plt.subplots(figsize=(8, 8))
        ax_rhog.plot(g_lcfs_r_shifted, g_lcfs_z_shifted, "b-", lw=2.2,
                     label="gfile RSU/ZSU (shifted)")
        if r_bnd_rr is not None:
            ax_rhog.plot(r_bnd_rr, z_bnd_rr, color="tab:orange", ls="--", lw=2.0,
                         label="RHOG(theta) mapped about (RR,0)")
        if r_bnd_axis is not None:
            ax_rhog.plot(r_bnd_axis, z_bnd_axis, color="tab:green", ls=":", lw=2.0,
                         label="RHOG(theta) mapped about axis")
        if eq_zero_contours:
            for i, seg in enumerate(eq_zero_contours):
                label = "eq psi=0 contour" if i == 0 else None
                ax_rhog.plot(seg["r"], seg["z"], color="tab:red", ls="-.", lw=1.4, label=label)
        # Diagnostic spokes from axis to LCFS points to compare with radial artifacts
        # seen in |grad psi| / Bp maps.
        n_spoke = 24
        idx_spoke = np.linspace(0, len(g_lcfs_r_shifted) - 1, n_spoke, dtype=int)
        for ii in idx_spoke:
            ax_rhog.plot(
                [e_axis[0], g_lcfs_r_shifted[ii]],
                [e_axis[1], g_lcfs_z_shifted[ii]],
                color="0.6",
                lw=0.8,
                alpha=0.55,
            )
        ax_rhog.plot([e_axis[0]], [e_axis[1]], "ro", ms=6, label="eq axis")
        ax_rhog.set_xlabel("R (m)")
        ax_rhog.set_ylabel("Z (m)")
        ax_rhog.set_aspect("equal")
        ax_rhog.grid(alpha=0.25)
        ax_rhog.set_title("RSU/ZSU -> RHOG(theta) -> mapped boundary + spokes")
        ax_rhog.legend(fontsize=8)
        rhog_plot = "lcfs_rhog_mapping_compare.png"
        fig_rhog.savefig(rhog_plot, dpi=150)
        print(f"Saved: {rhog_plot}")

        # Mapping-formula explanation figure.
        fig_mapexp, axes_mapexp = plt.subplots(1, 2, figsize=(14, 7), constrained_layout=True)
        ax0, ax1 = axes_mapexp
        ax0.plot(g_lcfs_r_shifted, g_lcfs_z_shifted, "b-", lw=2.2, label="shifted gfile LCFS")
        if r_bnd_rr is not None:
            ax0.plot(r_bnd_rr, z_bnd_rr, color="tab:orange", ls="--", lw=2.0,
                     label="sigma=1 from RHOG(theta)")
            for frac, col in zip([0.25, 0.5, 0.75], ["0.6", "0.45", "0.3"]):
                ax0.plot(
                    rr_eq + frac * rhg * np.cos(thg),
                    frac * rhg * np.sin(thg),
                    color=col,
                    lw=1.0,
                    label=f"sigma={frac:.2f}" if frac == 0.25 else None,
                )
        for ii in idx_spoke:
            ax0.plot(
                [rr_eq, g_lcfs_r_shifted[ii]],
                [0.0, g_lcfs_z_shifted[ii]],
                color="0.55",
                lw=0.8,
                alpha=0.6,
            )
        ax0.plot([rr_eq], [0.0], marker="o", color="tab:orange", ms=6, label="(RR,0)")
        ax0.plot([e_axis[0]], [e_axis[1]], marker="*", color="tab:red", ms=10, label="eq axis")
        ax0.set_title("Geometric mapping used by TASK/EQ")
        ax0.set_xlabel("R (m)")
        ax0.set_ylabel("Z (m)")
        ax0.set_aspect("equal")
        ax0.grid(alpha=0.25)
        ax0.legend(fontsize=8)
        ax0.text(
            0.02,
            0.98,
            "theta = atan2(Z, R-RR)\n"
            "sigma = sqrt((R-RR)^2 + Z^2) / RHOL(theta)",
            transform=ax0.transAxes,
            va="top",
            ha="left",
            fontsize=9,
            bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "none"},
        )

        ax1.plot(g_lcfs_r_shifted, g_lcfs_z_shifted, "b-", lw=2.2, label="shifted gfile LCFS")
        if r_bnd_rr is not None:
            ax1.plot(r_bnd_rr, z_bnd_rr, color="tab:orange", ls="--", lw=2.0,
                     label="sigma=1 boundary")
        if eq_zero_contours:
            for i, seg in enumerate(eq_zero_contours):
                ax1.plot(seg["r"], seg["z"], color="tab:red", ls="-.", lw=1.4,
                         label="psi=0 contour" if i == 0 else None)
        for ii in idx_spoke:
            ax1.plot(
                [rr_eq, g_lcfs_r_shifted[ii]],
                [0.0, g_lcfs_z_shifted[ii]],
                color="0.55",
                lw=0.8,
                alpha=0.5,
            )
        ax1.plot([rr_eq], [0.0], marker="o", color="tab:orange", ms=6, label="(RR,0)")
        ax1.plot([e_axis[0]], [e_axis[1]], marker="*", color="tab:red", ms=10, label="eq axis")
        ax1.set_title("Where the mapping geometry and psi=0 differ")
        ax1.set_xlabel("R (m)")
        ax1.set_ylabel("Z (m)")
        ax1.set_aspect("equal")
        ax1.grid(alpha=0.25)
        ax1.legend(fontsize=8)

        mapexp_plot = "sigma_theta_mapping_explain.png"
        fig_mapexp.savefig(mapexp_plot, dpi=160)
        print(f"Saved: {mapexp_plot}")

        # Explicit sigma/theta grid visualizations.
        sigma_levels = np.linspace(0.0, 1.0, 9)
        theta_samples = np.linspace(0.0, 2.0 * np.pi, 25)[:-1]

        fig_grid_rz, ax_grid_rz = plt.subplots(figsize=(8, 8))
        ax_grid_rz.plot(g_lcfs_r_shifted, g_lcfs_z_shifted, "k-", lw=2.0, label="LCFS")
        if thg is not None and rhg is not None:
            # sigma = const curves
            for s in sigma_levels[1:]:
                ax_grid_rz.plot(
                    rr_eq + s * rhg * np.cos(thg),
                    0.0 + s * rhg * np.sin(thg),
                    color="tab:blue",
                    lw=0.9 if s < 1.0 else 1.6,
                    alpha=0.8,
                )
            # theta = const rays
            rho_theta = np.interp(theta_samples, thg, rhg)
            for th, rho_l in zip(theta_samples, rho_theta):
                ss = np.linspace(0.0, 1.0, 80)
                ax_grid_rz.plot(
                    rr_eq + ss * rho_l * np.cos(th),
                    0.0 + ss * rho_l * np.sin(th),
                    color="tab:orange",
                    lw=0.8,
                    alpha=0.55,
                )
        ax_grid_rz.plot([rr_eq], [0.0], marker="o", color="tab:red", ms=6, label="mapping center (RR,0)")
        ax_grid_rz.set_title("Internal sigma/theta geometry grid in (R,Z)")
        ax_grid_rz.set_xlabel("R (m)")
        ax_grid_rz.set_ylabel("Z (m)")
        ax_grid_rz.set_aspect("equal")
        ax_grid_rz.grid(alpha=0.25)
        ax_grid_rz.legend(fontsize=8)
        grid_rz_plot = "sigma_theta_grid_rz.png"
        fig_grid_rz.savefig(grid_rz_plot, dpi=160)
        print(f"Saved: {grid_rz_plot}")

        fig_grid_st, ax_grid_st = plt.subplots(figsize=(8, 5.5))
        for s in sigma_levels:
            ax_grid_st.plot([0.0, 2.0 * np.pi], [s, s], color="tab:blue", lw=0.9)
        for th in np.linspace(0.0, 2.0 * np.pi, 25):
            ax_grid_st.plot([th, th], [0.0, 1.0], color="tab:orange", lw=0.8)
        ax_grid_st.set_xlim(0.0, 2.0 * np.pi)
        ax_grid_st.set_ylim(0.0, 1.0)
        ax_grid_st.set_xlabel(r"theta")
        ax_grid_st.set_ylabel(r"sigma")
        ax_grid_st.set_title("Internal computational grid in (sigma, theta)")
        ax_grid_st.grid(alpha=0.25)
        ax_grid_st.text(
            0.02,
            0.98,
            "Horizontal lines: sigma = const\nVertical lines: theta = const",
            transform=ax_grid_st.transAxes,
            va="top",
            ha="left",
            fontsize=9,
            bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "none"},
        )
        grid_st_plot = "sigma_theta_grid_plane.png"
        fig_grid_st.savefig(grid_st_plot, dpi=160)
        print(f"Saved: {grid_st_plot}")

    elif lcfs_eq_r is None:
        print("Could not extract eqdata LCFS contour at PSI=0.")

    g_geo = g if (g is not None and "lcfs_r" in g and "axis" in g) else None
    if g_profile is not None:
        make_overview_figure(eq, g_geo, metrics, lcfs_eq_r, lcfs_eq_z, lcfs_m, axis_d, axis_dist)

    # ------------------------------
    # 2b) Equilibrium chain diagnostics for q
    # ------------------------------
    if g is not None and g_profile is not None:
        eqcalq_raw_path = resolve_eqcalq_raw(
            current_dir, eq_dir, prefix, csv_arg=args.eqcalq_raw_csv
        )
        if eqcalq_raw_path is not None:
            try:
                tr = read_eqcalq_raw(eqcalq_raw_path)
                # Use eq radial profile psi_p grid as the common abscissa.
                # This matches the selected prefix instead of the potentially stale
                # generic eqcalq raw sampling.
                x_eq = np.asarray(eq["psi_n_q"], dtype=float)
                psip_eq = np.asarray(pd.read_csv(paths["radial"])["PSIPV"].to_numpy(dtype=float), dtype=float)
                psit_eq = np.asarray(pd.read_csv(paths["radial"])["PSITV"].to_numpy(dtype=float), dtype=float)
                q_eq = np.asarray(pd.read_csv(paths["radial"])["QPV"].to_numpy(dtype=float), dtype=float)
                tts_eq = np.asarray(pd.read_csv(paths["radial"])["TTV"].to_numpy(dtype=float), dtype=float)

                psipa_g = float(-np.nanmin(g["PSIRZ"]))
                psi_n_g = np.asarray(g_profile["psi_n_profile"], dtype=float)
                q_g = np.asarray(g_profile["q"], dtype=float)
                tts_g = np.asarray(g_profile["F"], dtype=float)
                psip_g = psi_n_g * psipa_g
                psit_g = integrate_psit(psip_g, q_g)

                axis_shift = (float(axis_d[0]), float(axis_d[1])) if axis_d is not None else (0.0, 0.0)
                g_geom = compute_gfile_geometry_chain(g, x_eq, axis_shift=axis_shift)
                if g_geom is not None:
                    sumavir2_g = g_geom["sumavir2"]
                    dvdpsip_g = g_geom["dvdpsip"]
                else:
                    sumavir2_g = np.full_like(x_eq, np.nan)
                    dvdpsip_g = np.full_like(x_eq, np.nan)

                refs = {
                    "PSIP": np.interp(x_eq, psi_n_g, psip_g),
                    "PSIT": np.interp(x_eq, psi_n_g, psit_g),
                    "q": np.interp(x_eq, psi_n_g, q_g),
                    "TTS": np.interp(x_eq, psi_n_g, tts_g),
                    "SUMAVIR2": sumavir2_g,
                    "dV/dpsi_p": dvdpsip_g,
                }
                eqvals = {
                    "PSIP": psip_eq,
                    "PSIT": psit_eq,
                    "q": q_eq,
                    "TTS": tts_eq,
                    "SUMAVIR2": np.interp(x_eq, tr["psip"] / max(float(tr["psip"][-1]), 1e-30), tr["sumavir2"]),
                    "dV/dpsi_p": np.interp(x_eq, tr["psip"] / max(float(tr["psip"][-1]), 1e-30), tr["dvdpsip"]),
                }
                # Keep the two different eq q definitions explicit:
                #   QPV: from EQCALV radial_profiles
                #   QPS/raw qps: from EQCALQ raw chain
                q_qps_raw = np.interp(x_eq, tr["psip"] / max(float(tr["psip"][-1]), 1e-30), tr["qps"])
                # Reconstruct q from the *same* raw EQCALQ chain.
                q_from_geom_ref = refs["TTS"] * refs["SUMAVIR2"] / (4.0 * np.pi ** 2)
                q_from_geom_eq = np.interp(x_eq, tr["psip"] / max(float(tr["psip"][-1]), 1e-30), tr["tts"]) * eqvals["SUMAVIR2"] / (4.0 * np.pi ** 2)
                q_from_psit_psip_eq = reconstruct_q_from_psit_psip(eqvals["PSIT"], eqvals["PSIP"])
                s_eq = float(args.eq_psi_scale)
                eqvals_scaled = {
                    "PSIP": s_eq * eqvals["PSIP"],
                    "PSIT": s_eq * eqvals["PSIT"],
                    "q": eqvals["q"],
                    "TTS": eqvals["TTS"],
                    "SUMAVIR2": eqvals["SUMAVIR2"] / s_eq,
                    "dV/dpsi_p": eqvals["dV/dpsi_p"] / s_eq,
                }
                eqvals_psip_only = {
                    "PSIP": s_eq * eqvals["PSIP"],
                    "PSIT": eqvals["PSIT"],
                    "q": eqvals["q"],
                    "TTS": eqvals["TTS"],
                    "SUMAVIR2": eqvals["SUMAVIR2"],
                    "dV/dpsi_p": eqvals["dV/dpsi_p"],
                }

                print(f"\n=== q-Chain Diagnostics (common psi_p grid, raw source: {eqcalq_raw_path}) ===")
                for key in ["PSIP", "PSIT", "q", "TTS", "SUMAVIR2", "dV/dpsi_p"]:
                    ref = np.asarray(refs[key], dtype=float)
                    val = np.asarray(eqvals[key], dtype=float)
                    m = np.isfinite(ref) & np.isfinite(val)
                    if not np.any(m):
                        print(f"{key}: no valid overlap")
                        continue
                    diff = val[m] - ref[m]
                    rmse = float(np.sqrt(np.mean(diff ** 2)))
                    mae = float(np.mean(np.abs(diff)))
                    print(f"{key}: RMSE={rmse:.4e}, MAE={mae:.4e}")
                m = np.isfinite(refs["q"]) & np.isfinite(q_from_geom_ref)
                if np.any(m):
                    print(
                        f"q_ref vs TTS*SUMAVIR2/(4pi^2): RMSE="
                        f"{float(np.sqrt(np.mean((refs['q'][m]-q_from_geom_ref[m])**2))):.4e}"
                    )
                m = np.isfinite(q_qps_raw) & np.isfinite(q_from_geom_eq)
                if np.any(m):
                    print(
                        f"q_raw(QPS) vs TTS*SUMAVIR2/(4pi^2): RMSE="
                        f"{float(np.sqrt(np.mean((q_qps_raw[m]-q_from_geom_eq[m])**2))):.4e}"
                    )
                m = np.isfinite(eqvals['q']) & np.isfinite(q_qps_raw)
                if np.any(m):
                    print(
                        f"QPV vs QPS(raw): RMSE="
                        f"{float(np.sqrt(np.mean((eqvals['q'][m]-q_qps_raw[m])**2))):.4e}"
                    )
                m = np.isfinite(eqvals['q']) & np.isfinite(q_from_psit_psip_eq)
                if np.any(m):
                    print(
                        f"QPV vs dPSIT/dPSIP: RMSE="
                        f"{float(np.sqrt(np.mean((eqvals['q'][m]-q_from_psit_psip_eq[m])**2))):.4e}"
                    )
                print(f"--- With eq PSI globally scaled by {s_eq:.3f} ---")
                for key in ["PSIP", "PSIT", "q", "TTS", "SUMAVIR2", "dV/dpsi_p"]:
                    ref = np.asarray(refs[key], dtype=float)
                    val = np.asarray(eqvals_scaled[key], dtype=float)
                    m = np.isfinite(ref) & np.isfinite(val)
                    if not np.any(m):
                        continue
                    diff = val[m] - ref[m]
                    rmse = float(np.sqrt(np.mean(diff ** 2)))
                    mae = float(np.mean(np.abs(diff)))
                    print(f"{key} (scaled): RMSE={rmse:.4e}, MAE={mae:.4e}")
                print(f"--- With only eq PSIP scaled by {s_eq:.3f} ---")
                for key in ["PSIP", "PSIT", "q", "TTS", "SUMAVIR2", "dV/dpsi_p"]:
                    ref = np.asarray(refs[key], dtype=float)
                    val = np.asarray(eqvals_psip_only[key], dtype=float)
                    m = np.isfinite(ref) & np.isfinite(val)
                    if not np.any(m):
                        continue
                    diff = val[m] - ref[m]
                    rmse = float(np.sqrt(np.mean(diff ** 2)))
                    mae = float(np.mean(np.abs(diff)))
                    print(f"{key} (PSIP-only): RMSE={rmse:.4e}, MAE={mae:.4e}")

                fig_chain, axes_chain = plt.subplots(3, 2, figsize=(14, 14), constrained_layout=True)
                items = ["PSIP", "PSIT", "q", "TTS", "SUMAVIR2", "dV/dpsi_p"]
                for ax, key in zip(axes_chain.ravel(), items):
                    ax.plot(x_eq, refs[key], "b-", lw=2.0, label="gfile shifted-geom ref")
                    ax.plot(x_eq, eqvals[key], "r--", lw=1.8, label="eqdata / eqcalq")
                    if key == "q":
                        ax.cla()
                        ax.plot(x_eq, refs["q"], "b-", lw=2.0, label="gfile q (shifted geom ref)")
                        ax.plot(x_eq, eqvals["q"], "r--", lw=1.8, label="eq QPV (EQCALV)")
                        ax.plot(
                            x_eq,
                            q_qps_raw,
                            color="tab:brown",
                            ls="--",
                            lw=1.8,
                            alpha=0.9,
                            label="eq QPS raw (EQCALQ)",
                        )
                        ax.plot(
                            x_eq,
                            q_from_geom_eq,
                            color="tab:orange",
                            ls=":",
                            lw=2.2,
                            alpha=0.9,
                            label="eq q reconstructed from TTS*SUMAVIR2",
                        )
                        ax.plot(
                            x_eq,
                            q_from_psit_psip_eq,
                            color="tab:green",
                            ls="-.",
                            lw=1.8,
                            alpha=0.95,
                            label="eq q reconstructed from dPSIT/dPSIP",
                        )
                    ax.set_title(key)
                    ax.set_xlabel("normalized psi_p")
                    ax.grid(alpha=0.25)
                    ax.legend(fontsize=8)
                fig_chain.suptitle("Equilibrium Chain: PSIP, PSIT, q, TTS, SUMAVIR2, dV/dpsi_p", fontsize=14)
                chain_plot = "equilibrium_chain_compare.png"
                fig_chain.savefig(chain_plot, dpi=160)
                print(f"Saved: {chain_plot}")

                # Dedicated QPV chain figure: QPV <-> dPSIT/dPSIP on radial_profiles grid.
                fig_qpv, axes_qpv = plt.subplots(2, 2, figsize=(12, 10), constrained_layout=True)
                axes_qpv[0, 0].plot(x_eq, refs["PSIP"], "b-", lw=2.0, label="gfile PSIP")
                axes_qpv[0, 0].plot(x_eq, eqvals["PSIP"], "r--", lw=1.8, label="eq PSIP")
                axes_qpv[0, 0].set_title("PSIP (QPV chain grid)")
                axes_qpv[0, 0].grid(alpha=0.25)
                axes_qpv[0, 0].legend(fontsize=8)

                axes_qpv[0, 1].plot(x_eq, refs["PSIT"], "b-", lw=2.0, label="gfile PSIT")
                axes_qpv[0, 1].plot(x_eq, eqvals["PSIT"], "r--", lw=1.8, label="eq PSIT")
                axes_qpv[0, 1].set_title("PSIT (QPV chain grid)")
                axes_qpv[0, 1].grid(alpha=0.25)
                axes_qpv[0, 1].legend(fontsize=8)

                axes_qpv[1, 0].plot(x_eq, refs["q"], "b-", lw=2.0, label="gfile q")
                axes_qpv[1, 0].plot(x_eq, eqvals["q"], "r--", lw=1.8, label="eq QPV")
                axes_qpv[1, 0].plot(x_eq, q_from_psit_psip_eq, color="tab:green", ls="-.", lw=1.8,
                                    label="dPSIT/dPSIP")
                axes_qpv[1, 0].set_title("QPV Chain: q")
                axes_qpv[1, 0].grid(alpha=0.25)
                axes_qpv[1, 0].legend(fontsize=8)

                axes_qpv[1, 1].plot(x_eq, refs["TTS"], "b-", lw=2.0, label="gfile TTS")
                axes_qpv[1, 1].plot(x_eq, eqvals["TTS"], "r--", lw=1.8, label="eq TTV")
                axes_qpv[1, 1].set_title("TTS/TTV (context for QPV chain)")
                axes_qpv[1, 1].grid(alpha=0.25)
                axes_qpv[1, 1].legend(fontsize=8)

                for ax in axes_qpv.ravel():
                    ax.set_xlabel("normalized psi_p")
                fig_qpv.suptitle("QPV Chain: PSIP -> PSIT -> QPV", fontsize=14)
                qpv_plot = "equilibrium_chain_qpv.png"
                fig_qpv.savefig(qpv_plot, dpi=160)
                print(f"Saved: {qpv_plot}")

                # Dedicated QPS chain figure: QPS(raw) <-> TTS*SUMAVIR2/(4pi^2) on EQCALQ raw grid.
                inside_raw = tr["inside"] > 0.5
                x_raw = tr["psip"][inside_raw]
                tts_raw = tr["tts"][inside_raw]
                sumavir2_raw = tr["sumavir2"][inside_raw]
                qps_raw = tr["qps"][inside_raw]
                dvdpsip_raw = tr["dvdpsip"][inside_raw]
                qps_geom_raw = tts_raw * sumavir2_raw / (4.0 * np.pi ** 2)
                x_raw_n = x_raw / max(float(np.nanmax(x_raw)), 1e-30)
                tts_g_raw = np.interp(x_raw_n, psi_n_g, tts_g)
                q_g_raw = np.interp(x_raw_n, psi_n_g, q_g)

                # Build explicit TTPS/UTTPS-chain reference on raw PSIP grid.
                psips_eq = np.asarray(eq["psi_p"], dtype=float)
                ttps_eq = np.asarray(eq["F"], dtype=float)
                tts_from_ttps_on_raw = np.interp(x_raw, psips_eq, ttps_eq)
                # Gfile support/abscissa for comparison on raw-PSIP panels.
                psips_g = psip_g
                tts_from_g_on_raw = np.interp(x_raw, psips_g, tts_g, left=np.nan, right=np.nan)

                fig_qps, axes_qps = plt.subplots(3, 2, figsize=(12, 14), constrained_layout=True)
                axes_qps[0, 0].plot(x_raw, q_g_raw, "b-", lw=2.0, label="gfile q")
                axes_qps[0, 0].plot(x_raw, qps_raw, "r--", lw=1.8, label="eq QPS raw")
                axes_qps[0, 0].plot(x_raw, qps_geom_raw, color="tab:orange", ls=":", lw=2.2,
                                    label="raw TTS*SUMAVIR2/(4pi^2)")
                axes_qps[0, 0].set_title("QPS Chain: q")
                axes_qps[0, 0].grid(alpha=0.25)
                axes_qps[0, 0].legend(fontsize=8)

                axes_qps[0, 1].plot(x_raw, tts_g_raw, "b-", lw=2.0, label="gfile TTS")
                axes_qps[0, 1].plot(x_raw, tts_raw, "r--", lw=1.8, label="eq raw TTS")
                axes_qps[0, 1].plot(x_raw, tts_from_ttps_on_raw, color="tab:purple", ls="-", lw=1.6,
                                    label="TTFUNC via TTPS/UTTPS on raw PSIP")
                # Diagnostic only: test whether a raw-PSIP x-scale factor of 2 helps alignment.
                tts_raw_x2 = np.interp(x_raw, 2.0 * tr["psip"], tr["tts"], left=np.nan, right=np.nan)
                axes_qps[0, 1].plot(x_raw, tts_raw_x2, color="tab:green", ls="-.", lw=1.8,
                                    label="eq raw TTS (x-axis x2 test)")
                axes_qps[0, 1].set_title("TTS (QPS raw chain)")
                axes_qps[0, 1].grid(alpha=0.25)
                axes_qps[0, 1].legend(fontsize=8)

                axes_qps[1, 0].plot(x_raw, np.interp(x_raw_n, x_eq, refs["SUMAVIR2"]), "b-", lw=2.0,
                                    label="gfile shifted SUMAVIR2")
                axes_qps[1, 0].plot(x_raw, sumavir2_raw, "r--", lw=1.8, label="eq raw SUMAVIR2")
                axes_qps[1, 0].set_title("SUMAVIR2 (QPS raw chain)")
                axes_qps[1, 0].grid(alpha=0.25)
                axes_qps[1, 0].legend(fontsize=8)

                axes_qps[1, 1].plot(x_raw, np.interp(x_raw_n, x_eq, refs["dV/dpsi_p"]), "b-", lw=2.0,
                                    label="gfile shifted dV/dpsi_p")
                axes_qps[1, 1].plot(x_raw, dvdpsip_raw, "r--", lw=1.8, label="eq raw dV/dpsi_p")
                axes_qps[1, 1].set_title("dV/dpsi_p (QPS raw chain)")
                axes_qps[1, 1].grid(alpha=0.25)
                axes_qps[1, 1].legend(fontsize=8)

                # Explicit mapping diagnostics: raw PSIP and TTPS backbone.
                axes_qps[2, 0].plot(x_raw, x_raw, "r--", lw=1.8, label="raw PSIP (EQCALQ, inside LCFS)")
                axes_qps[2, 0].plot(psips_eq, psips_eq, color="tab:purple", lw=1.6,
                                    label="PSIPS support (TTPS spline)")
                axes_qps[2, 0].plot(psips_g, psips_g, color="tab:blue", lw=1.4,
                                    label="gfile PSIP support")
                axes_qps[2, 0].set_title("PSIP(raw) and TTPS support grid")
                axes_qps[2, 0].grid(alpha=0.25)
                axes_qps[2, 0].legend(fontsize=8)

                axes_qps[2, 1].plot(psips_eq, ttps_eq, color="tab:purple", lw=1.8,
                                    label="TTPS on PSIPS support")
                axes_qps[2, 1].plot(psips_g, tts_g, color="tab:blue", lw=1.6,
                                    label="gfile TTS on gfile PSIP")
                axes_qps[2, 1].plot(x_raw, tts_from_ttps_on_raw, "k--", lw=1.5,
                                    label="TTFUNC(PSIP raw)")
                axes_qps[2, 1].plot(x_raw, tts_from_g_on_raw, color="tab:cyan", ls="-.", lw=1.5,
                                    label="gfile TTS sampled on raw PSIP")
                axes_qps[2, 1].plot(x_raw, tts_raw, "r:", lw=2.0,
                                    label="stored TTS raw")
                axes_qps[2, 1].set_title("PSIP(raw) -> TTFUNC -> TTPS/UTTPS -> TTS(raw)")
                axes_qps[2, 1].grid(alpha=0.25)
                axes_qps[2, 1].legend(fontsize=8)

                for ax in axes_qps.ravel():
                    ax.set_xlabel("raw PSIP")
                fig_qps.suptitle("QPS Raw Chain: TTS + SUMAVIR2 -> QPS", fontsize=14)
                qps_plot = "equilibrium_chain_qps.png"
                fig_qps.savefig(qps_plot, dpi=160)
                print(f"Saved: {qps_plot}")
            except Exception as e:
                print(f"Warning: failed to build equilibrium chain diagnostics: {e}")

    # ------------------------------
    # 3) Supplementary PSI(R,Z) plots
    # ------------------------------
    datasets = []
    if g is not None:
        datasets.append(g["PSIRZ"])
    if psi_gs2d is not None:
        datasets.append(psi_gs2d)
    datasets.append(eq["PSIRZ"])

    vmin = min(np.percentile(d, 1) for d in datasets)
    vmax = max(np.percentile(d, 99) for d in datasets)
    if vmax == vmin:
        vmax += 1.0
    levels = np.linspace(vmin, vmax, 40)

    fig, axes = plt.subplots(1, 3, figsize=(20, 6), constrained_layout=True)

    ax = axes[0]
    if g is not None:
        cf = ax.contourf(np.meshgrid(g["R"], g["Z"])[0], np.meshgrid(g["R"], g["Z"])[1], g["PSIRZ"], levels=levels, cmap="jet")
        fig.colorbar(cf, ax=ax)
        ax.plot(g["lcfs_r"], g["lcfs_z"], "w--", lw=1)
        ax.set_title("G-file PSI(R,Z)")
    else:
        ax.text(0.5, 0.5, "G-file missing", ha="center", va="center")
    ax.set_xlabel("R (m)")
    ax.set_ylabel("Z (m)")
    ax.set_aspect("equal")

    ax = axes[1]
    if psi_gs2d is not None:
        r2, z2 = np.meshgrid(r_gs2d, z_gs2d)
        cf = ax.contourf(r2, z2, psi_gs2d, levels=levels, cmap="jet")
        fig.colorbar(cf, ax=ax)
        ax.set_title("EQGS2D PSI(R,Z)")
    else:
        ax.text(0.5, 0.5, "EQGS2D missing", ha="center", va="center")
    ax.set_xlabel("R (m)")
    ax.set_aspect("equal")

    ax = axes[2]
    r3, z3 = np.meshgrid(eq["R"], eq["Z"])
    cf = ax.contourf(r3, z3, eq["PSIRZ"], levels=levels, cmap="jet")
    fig.colorbar(cf, ax=ax)
    ax.set_title("EQDATA PSI(R,Z)")
    ax.set_xlabel("R (m)")
    ax.set_aspect("equal")

    fig.suptitle("PSI(R,Z) Supplementary Comparison", fontsize=15)
    psi_plot = "psi_comparison_3way.png"
    fig.savefig(psi_plot, dpi=150)
    print(f"Saved: {psi_plot}")

    # Keep the legacy overlay figure for quick visual comparison.
    # Shift gfile grid so that gfile magnetic axis aligns with eqdata axis.
    if g is not None:
        g_axis = np.asarray(g["axis"], dtype=float)
        e_axis = np.asarray(eq["axis"], dtype=float)
        dR_shift = e_axis[0] - g_axis[0]
        dZ_shift = e_axis[1] - g_axis[1]
        R_g_shifted = g["R"] + dR_shift
        Z_g_shifted = g["Z"] + dZ_shift

        fig_ov, ax_ov = plt.subplots(figsize=(8, 10))
        r_g2d, z_g2d = np.meshgrid(R_g_shifted, Z_g_shifted)
        r_e2d, z_e2d = np.meshgrid(eq["R"], eq["Z"])
        # Use a sparse labeled subset so psi values are readable on the overlay.
        label_levels = np.unique(
            np.concatenate(
                [
                    np.linspace(vmin, vmax, 8),
                    np.array([0.0, float(np.nanmin(g["PSIRZ"])), float(np.nanmin(eq["PSIRZ"]))]),
                ]
            )
        )
        label_levels = label_levels[(label_levels >= vmin) & (label_levels <= vmax)]

        cs_g = ax_ov.contour(r_g2d, z_g2d, g["PSIRZ"], levels=levels, colors="blue", linewidths=1.2)
        cs_e = ax_ov.contour(
            r_e2d,
            z_e2d,
            eq["PSIRZ"],
            levels=levels,
            colors="red",
            linestyles="dashed",
            linewidths=1.2,
        )
        cs_g_lbl = ax_ov.contour(r_g2d, z_g2d, g["PSIRZ"], levels=label_levels, colors="blue", linewidths=0.0)
        cs_e_lbl = ax_ov.contour(r_e2d, z_e2d, eq["PSIRZ"], levels=label_levels, colors="red", linewidths=0.0)
        ax_ov.clabel(cs_g_lbl, fmt=lambda x: f"{x:.0f}", fontsize=7, colors="blue")
        ax_ov.clabel(cs_e_lbl, fmt=lambda x: f"{x:.0f}", fontsize=7, colors="red")
        shift_str = f"dR={dR_shift:+.3f}, dZ={dZ_shift:+.3f} m"
        psipa_g = float(-np.nanmin(g["PSIRZ"]))
        psipa_e = float(-np.nanmin(eq["PSIRZ"]))
        ax_ov.set_title(
            f"Overlay: G-file (Blue, shifted {shift_str}) vs Eqdata (Red)"
        )
        ax_ov.set_xlabel("R (m)")
        ax_ov.set_ylabel("Z (m)")
        ax_ov.set_aspect("equal")
        ax_ov.grid(alpha=0.25)
        # Mark both axes (after shift they should overlap)
        ax_ov.plot(
            [e_axis[0]], [e_axis[1]], "r*", ms=12, zorder=10,
            label=f"eqdata axis ({e_axis[0]:.3f},{e_axis[1]:.3f})"
        )
        ax_ov.plot(
            [g_axis[0] + dR_shift], [g_axis[1] + dZ_shift], "b+", ms=12, mew=2, zorder=10,
            label=f"gfile axis (shifted)"
        )
        from matplotlib.lines import Line2D

        ax_ov.legend(
            handles=[
                Line2D([0], [0], color="blue", lw=1.5, label="G-file (axis-shifted)"),
                Line2D([0], [0], color="red", lw=1.5, linestyle="--", label="Eqdata"),
                Line2D([0], [0], marker="*", color="red", lw=0, ms=10, label="eqdata axis"),
            ],
            loc="upper right",
        )
        ax_ov.text(
            0.02,
            0.98,
            f"psi=0 is LCFS\nPSIPA_gfile={psipa_g:.2f}\nPSIPA_eq={psipa_e:.2f}",
            transform=ax_ov.transAxes,
            va="top",
            ha="left",
            fontsize=8,
            bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "none"},
        )

        # Show psi values evaluated at each magnetic axis.
        psi_g_axis = float(np.nanmin(g["PSIRZ"]))
        psi_e_axis = float(np.nanmin(eq["PSIRZ"]))
        ax_ov.text(
            0.02,
            0.84,
            (
                f"gfile axis psi={psi_g_axis:.2f}\n"
                f"eq axis psi={psi_e_axis:.2f}"
            ),
            transform=ax_ov.transAxes,
            va="top",
            ha="left",
            fontsize=8,
            color="black",
            bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "none"},
        )
        overlay_plot = "gfile_eqdata_overlay.png"
        fig_ov.savefig(overlay_plot, dpi=150)
        print(f"Saved: {overlay_plot}  (gfile shifted: {shift_str})")

        psi_diff_csv = export_surface_point_psi_differences(
            f"{prefix}_psi_surface_point_differences.csv",
            eq,
            g,
            axis_shift=(dR_shift, dZ_shift),
            sigma_levels=np.linspace(0.1, 1.0, 10),
            ntheta=72,
        )
        print(f"Saved: {psi_diff_csv}")

        # Direct delta-psi map on the eq grid after shifting gfile to the eq axis.
        psi_g_on_eq = interp_field_to_grid(
            R_g_shifted, Z_g_shifted, np.asarray(g["PSIRZ"], dtype=float), eq["R"], eq["Z"]
        )
        delta_psi = np.asarray(eq["PSIRZ"], dtype=float) - psi_g_on_eq
        fig_dpsi, ax_dpsi = plt.subplots(figsize=(8, 10), constrained_layout=True)
        vv = np.nanmax(np.abs(delta_psi))
        if not np.isfinite(vv) or vv <= 0:
            vv = 1.0
        cs = ax_dpsi.contourf(np.meshgrid(eq["R"], eq["Z"])[0], np.meshgrid(eq["R"], eq["Z"])[1],
                              delta_psi, levels=40, cmap="coolwarm", vmin=-vv, vmax=vv)
        fig_dpsi.colorbar(cs, ax=ax_dpsi, label="psi_eq - psi_gfile_shifted")
        ax_dpsi.plot(g["lcfs_r"] + dR_shift, g["lcfs_z"] + dZ_shift, "k--", lw=1.2, label="shifted gfile LCFS")
        if lcfs_eq_r is not None:
            ax_dpsi.plot(lcfs_eq_r, lcfs_eq_z, "g-", lw=1.0, label="eq psi=0 contour")
        ax_dpsi.plot([e_axis[0]], [e_axis[1]], "wo", ms=4)
        ax_dpsi.set_title("Step 3: delta psi(R,Z) on common geometry")
        ax_dpsi.set_xlabel("R (m)")
        ax_dpsi.set_ylabel("Z (m)")
        ax_dpsi.set_aspect("equal")
        ax_dpsi.legend(fontsize=8)
        delta_plot = "delta_psi_map.png"
        fig_dpsi.savefig(delta_plot, dpi=160)
        print(f"Saved: {delta_plot}")

        # Local flux-gradient / Bp diagnostic.
        dpsi_dz_g, dpsi_dr_g = np.gradient(np.asarray(g["PSIRZ"], dtype=float), g["Z"], g["R"])
        gradpsi_g = np.sqrt(dpsi_dr_g ** 2 + dpsi_dz_g ** 2)
        bp_g = gradpsi_g / np.maximum(2.0 * np.pi * g["R"][None, :], 1e-12)

        dpsi_dz_e, dpsi_dr_e = np.gradient(np.asarray(eq["PSIRZ"], dtype=float), eq["Z"], eq["R"])
        gradpsi_e = np.sqrt(dpsi_dr_e ** 2 + dpsi_dz_e ** 2)
        bp_e = gradpsi_e / np.maximum(2.0 * np.pi * eq["R"][None, :], 1e-12)

        gradpsi_g_on_eq = interp_field_to_grid(
            R_g_shifted, Z_g_shifted, gradpsi_g, eq["R"], eq["Z"]
        )
        bp_g_on_eq = interp_field_to_grid(
            R_g_shifted, Z_g_shifted, bp_g, eq["R"], eq["Z"]
        )

        fig_bp, axes_bp = plt.subplots(2, 3, figsize=(18, 10), constrained_layout=True)
        r_e2d, z_e2d = np.meshgrid(eq["R"], eq["Z"])
        n_spoke = 24
        idx_spoke = np.linspace(0, len(g["lcfs_r"]) - 1, n_spoke, dtype=int)

        for ax, arr, title in [
            (axes_bp[0, 0], gradpsi_g_on_eq, "|grad psi| gfile shifted -> eq grid"),
            (axes_bp[0, 1], gradpsi_e, "|grad psi| eqdata"),
            (axes_bp[1, 0], bp_g_on_eq, "Bp gfile shifted -> eq grid"),
            (axes_bp[1, 1], bp_e, "Bp eqdata"),
        ]:
            vv = np.nanmax(np.abs(arr))
            if not np.isfinite(vv) or vv <= 0:
                vv = 1.0
            cs = ax.contourf(r_e2d, z_e2d, arr, levels=40, cmap="viridis")
            fig_bp.colorbar(cs, ax=ax)
            ax.plot(g["lcfs_r"] + dR_shift, g["lcfs_z"] + dZ_shift, "w--", lw=1.0)
            if lcfs_eq_r is not None:
                ax.plot(lcfs_eq_r, lcfs_eq_z, "r-", lw=1.0)
            for ii in idx_spoke:
                ax.plot(
                    [e_axis[0], g["lcfs_r"][ii] + dR_shift],
                    [e_axis[1], g["lcfs_z"][ii] + dZ_shift],
                    color="w",
                    lw=0.7,
                    alpha=0.45,
                )
            ax.set_title(title)
            ax.set_xlabel("R (m)")
            ax.set_ylabel("Z (m)")
            ax.set_aspect("equal")

        grad_ratio = gradpsi_e / np.where(np.abs(gradpsi_g_on_eq) > 1e-12, gradpsi_g_on_eq, np.nan)
        bp_ratio = bp_e / np.where(np.abs(bp_g_on_eq) > 1e-12, bp_g_on_eq, np.nan)
        for ax, arr, title in [
            (axes_bp[0, 2], grad_ratio, "|grad psi| ratio eq/gfile"),
            (axes_bp[1, 2], bp_ratio, "Bp ratio eq/gfile"),
        ]:
            cs = ax.contourf(r_e2d, z_e2d, arr, levels=40, cmap="coolwarm")
            fig_bp.colorbar(cs, ax=ax)
            ax.plot(g["lcfs_r"] + dR_shift, g["lcfs_z"] + dZ_shift, "k--", lw=1.0)
            if lcfs_eq_r is not None:
                ax.plot(lcfs_eq_r, lcfs_eq_z, "g-", lw=1.0)
            for ii in idx_spoke:
                ax.plot(
                    [e_axis[0], g["lcfs_r"][ii] + dR_shift],
                    [e_axis[1], g["lcfs_z"][ii] + dZ_shift],
                    color="k",
                    lw=0.7,
                    alpha=0.35,
                )
            ax.set_title(title)
            ax.set_xlabel("R (m)")
            ax.set_ylabel("Z (m)")
            ax.set_aspect("equal")

        bp_plot = "psi_gradient_bp_compare.png"
        fig_bp.savefig(bp_plot, dpi=160)
        print(f"Saved: {bp_plot}")


if __name__ == "__main__":
    main()
