#!/usr/bin/env python3
"""
Solve fixed-boundary Grad-Shafranov equation from gfile/profile data.

Equation:
  d2psi/dR2 - (1/R) dpsi/dR + d2psi/dZ2 = -mu0 * R^2 * p'(psi) - F(psi)F'(psi)

Supports two source modes:
  - "pprime_ffprime" (default/legacy): directly use p'(psi) and FF'(psi)
  - "mdleqf6": use P(psi) and F(psi) profiles with TJ scaling to match Ip
    (equivalent to TASK/EQ MDLEQF=6: given spline P,F + Ip constraint)

This script uses:
  - Dirichlet boundary on either rectangular box or LCFS polygon
  - Picard outer iteration for nonlinear source update
  - Weighted Jacobi inner iteration for elliptic solve (0 < omega <= 1)
"""

import argparse
import os
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.path import Path as MplPath

MU0 = 4.0e-7 * np.pi


class GFile:
    pass


def rho_from_qpsi(qpsi, psimag, psibdy):
    """
    Compute rho_tor = sqrt(normalized toroidal flux) from q(psi) on uniform psi_n grid.
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


def build_lcfs_polygon_mask(r, z, r_lcfs, z_lcfs):
    """Return boolean mask [R,Z] for points inside LCFS polygon."""
    if len(r_lcfs) < 3:
        return None
    rr, zz = np.meshgrid(r, z, indexing="ij")
    pts = np.column_stack([rr.ravel(), zz.ravel()])
    poly = MplPath(np.column_stack([r_lcfs, z_lcfs]))
    inside = poly.contains_points(pts)
    return inside.reshape(rr.shape)


def largest_connected_component(mask):
    """Keep only the largest 4-connected True component in a [R,Z] boolean mask."""
    m = np.asarray(mask, dtype=bool)
    nr, nz = m.shape
    visited = np.zeros_like(m, dtype=bool)
    best = None
    best_count = 0

    for i0 in range(nr):
        for j0 in range(nz):
            if (not m[i0, j0]) or visited[i0, j0]:
                continue
            stack = [(i0, j0)]
            visited[i0, j0] = True
            comp = []
            while stack:
                i, j = stack.pop()
                comp.append((i, j))
                if i > 0 and m[i - 1, j] and not visited[i - 1, j]:
                    visited[i - 1, j] = True
                    stack.append((i - 1, j))
                if i + 1 < nr and m[i + 1, j] and not visited[i + 1, j]:
                    visited[i + 1, j] = True
                    stack.append((i + 1, j))
                if j > 0 and m[i, j - 1] and not visited[i, j - 1]:
                    visited[i, j - 1] = True
                    stack.append((i, j - 1))
                if j + 1 < nz and m[i, j + 1] and not visited[i, j + 1]:
                    visited[i, j + 1] = True
                    stack.append((i, j + 1))
            if len(comp) > best_count:
                best_count = len(comp)
                best = comp

    out = np.zeros_like(m, dtype=bool)
    if best is not None:
        ii, jj = zip(*best)
        out[np.asarray(ii), np.asarray(jj)] = True
    return out


def fill_single_cell_holes(mask):
    """Fill isolated 1-cell holes inside True region on [R,Z] grid."""
    m = np.asarray(mask, dtype=bool).copy()
    if m.shape[0] < 3 or m.shape[1] < 3:
        return m
    hole = (~m).copy()
    four_true = np.zeros_like(m, dtype=bool)
    four_true[1:-1, 1:-1] = (
        m[:-2, 1:-1] & m[2:, 1:-1] & m[1:-1, :-2] & m[1:-1, 2:]
    )
    m[hole & four_true] = True
    return m


def build_lcfs_update_masks(inside):
    """
    From inside-domain mask [R,Z], build:
      - interior_update: points with all 4 neighbors also inside
      - fixed_dirichlet: all other points (outside + LCFS boundary band)
    """
    # Array layout is [R,Z] -> axis 0: R (west/east), axis 1: Z (south/north).
    north = np.zeros_like(inside, dtype=bool)
    south = np.zeros_like(inside, dtype=bool)
    west = np.zeros_like(inside, dtype=bool)
    east = np.zeros_like(inside, dtype=bool)

    # Z-direction neighbors (north/south): shift along axis 1.
    north[:, :-1] = inside[:, 1:]
    south[:, 1:] = inside[:, :-1]

    # R-direction neighbors (west/east): shift along axis 0.
    west[1:, :] = inside[:-1, :]
    east[:-1, :] = inside[1:, :]
    all4_inside = north & south & west & east
    interior_update = inside & all4_inside
    # Dirichlet ring on plasma-side LCFS band only.
    lcfs_ring = inside & (~all4_inside)
    outside = ~inside
    return interior_update, lcfs_ring, outside


def segment_intersection_fraction(p0, p1, q0, q1, eps=1e-12):
    """Return t in [0,1] where segment p(t)=p0+t*(p1-p0) intersects q-segment; NaN if none."""
    x1, y1 = p0
    x2, y2 = p1
    x3, y3 = q0
    x4, y4 = q1
    den = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
    if abs(den) < eps:
        return np.nan
    t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / den
    u = ((x1 - x3) * (y1 - y2) - (y1 - y3) * (x1 - x2)) / den
    if -eps <= t <= 1.0 + eps and -eps <= u <= 1.0 + eps:
        return float(min(max(t, 0.0), 1.0))
    return np.nan


def polyline_first_hit_fraction(p0, p1, r_lcfs, z_lcfs):
    """Find first intersection fraction t from p0 to p1 against LCFS polyline."""
    ts = []
    n = len(r_lcfs)
    if n < 2:
        return np.nan
    for k in range(n):
        q0 = (r_lcfs[k], z_lcfs[k])
        q1 = (r_lcfs[(k + 1) % n], z_lcfs[(k + 1) % n])
        t = segment_intersection_fraction(p0, p1, q0, q1)
        if np.isfinite(t):
            ts.append(t)
    if not ts:
        return np.nan
    tmin = min(ts)
    return float(min(max(tmin, 1e-3), 1.0))


def precompute_lcfs_edge_fractions(r, z, inside, r_lcfs, z_lcfs, lambda_min=0.2):
    """
    For inside cells near LCFS, precompute intersection fractions to LCFS along
    west/east/south/north center-to-center segments.
    """
    nr = len(r)
    nz = len(z)
    lam_w = np.full((nr, nz), np.nan, dtype=float)
    lam_e = np.full((nr, nz), np.nan, dtype=float)
    lam_s = np.full((nr, nz), np.nan, dtype=float)
    lam_n = np.full((nr, nz), np.nan, dtype=float)

    def hit(i0, j0, i1, j1):
        p0 = (r[i0], z[j0])
        p1 = (r[i1], z[j1])
        t = polyline_first_hit_fraction(p0, p1, r_lcfs, z_lcfs)
        if np.isfinite(t):
            return float(min(max(t, lambda_min), 1.0))
        return 0.5

    for i in range(1, nr - 1):
        for j in range(1, nz - 1):
            if not inside[i, j]:
                continue
            if not inside[i - 1, j]:
                lam_w[i, j] = hit(i, j, i - 1, j)
            if not inside[i + 1, j]:
                lam_e[i, j] = hit(i, j, i + 1, j)
            if not inside[i, j - 1]:
                lam_s[i, j] = hit(i, j, i, j - 1)
            if not inside[i, j + 1]:
                lam_n[i, j] = hit(i, j, i, j + 1)
    return lam_w, lam_e, lam_s, lam_n


def to_plot_zr(matrix, r, z):
    """
    Return matrix in [Z,R] layout expected by contour/contourf with
    meshgrid(r, z) (default indexing="xy").

    Internal solver arrays use [R,Z] layout, so transpose them for plotting.
    """
    arr = np.asarray(matrix, dtype=float)
    nz = len(z)
    nr = len(r)
    if arr.shape == (nz, nr):
        return arr
    if arr.shape == (nr, nz):
        return arr
    raise ValueError(
        f"Unexpected matrix shape {arr.shape}, expected ({nz},{nr}) or ({nr},{nz})"
    )


def read_gfile(path):
    g = GFile()
    with open(path, "r", encoding="utf-8") as f:
        head = re.match(r"^\s*(.*)\s+(\d+)\s+(\d+)\s*$", f.readline())
        if not head:
            raise ValueError(f"Invalid gfile header: {path}")
        g.header = head.group(1)
        g.nw = int(head.group(2))
        g.nh = int(head.group(3))

        def read5():
            line = re.sub(r"([^Ee])\-", r"\1 -", f.readline())
            vals = [float(x) for x in re.split(r"\s+", line.strip()) if x]
            if len(vals) < 5:
                raise ValueError("Failed to parse 5-value row in gfile.")
            return vals[:5]

        r2 = read5()
        r3 = read5()
        r4 = read5()
        _ = read5()  # row5 not needed directly

        g.rdim, g.zdim, g.rcentr, g.rleft, g.zmid = r2
        g.rmaxis, g.zmaxis, g.psimag, g.psibdy, g.bcentr = r3
        g.currentA = r4[0]

        body = f.read().replace("\n", " ")
    body = re.sub(r"(\d)\-", r"\1 -", body)
    nums = [float(x) for x in re.split(r"\s+", body) if x]

    nw, nh = g.nw, g.nh
    ia = 0
    g.fpol = np.asarray(nums[ia : ia + nw], dtype=float)
    ia += nw
    g.pres = np.asarray(nums[ia : ia + nw], dtype=float)
    ia += nw
    g.ffprime = np.asarray(nums[ia : ia + nw], dtype=float)
    ia += nw
    g.pprime = np.asarray(nums[ia : ia + nw], dtype=float)
    ia += nw

    # GEQDSK psirz block is ordered with Z running outer and R running inner.
    # Convert to internal [R,Z] layout used throughout this solver.
    g.psirz = np.asarray(nums[ia : ia + nw * nh], dtype=float).reshape(nh, nw).T
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


def resolve_profile_csv(gfile_path, profile_csv):
    if profile_csv:
        if os.path.isabs(profile_csv):
            return profile_csv
        return os.path.abspath(profile_csv)
    b = os.path.basename(gfile_path)
    cands = [
        os.path.join(os.path.dirname(gfile_path), f"{b}_profiles_from_gfile.csv"),
        os.path.join(os.getcwd(), f"{b}_profiles_from_gfile.csv"),
    ]
    for c in cands:
        if os.path.exists(c):
            return c
    return None


def sanitize_axis_and_profiles(x, pprime, ffprime, axis_name):
    x = np.asarray(x, dtype=float)
    pprime = np.asarray(pprime, dtype=float)
    ffprime = np.asarray(ffprime, dtype=float)
    good = np.isfinite(x) & np.isfinite(pprime) & np.isfinite(ffprime)
    x = x[good]
    pprime = pprime[good]
    ffprime = ffprime[good]
    if x.size < 2:
        raise ValueError(f"Insufficient valid source points on axis '{axis_name}'.")
    idx = np.argsort(x)
    x = x[idx]
    pprime = pprime[idx]
    ffprime = ffprime[idx]

    # Merge duplicate x by arithmetic mean to make np.interp safe.
    xu, inv = np.unique(x, return_inverse=True)
    if xu.size != x.size:
        p_acc = np.bincount(inv, weights=pprime)
        f_acc = np.bincount(inv, weights=ffprime)
        c_acc = np.bincount(inv)
        pprime = p_acc / np.maximum(c_acc, 1.0)
        ffprime = f_acc / np.maximum(c_acc, 1.0)
        x = xu
    if x.size < 2:
        raise ValueError(f"Degenerate source axis '{axis_name}' after duplicate merge.")
    if not np.all(np.diff(x) > 0.0):
        raise ValueError(f"Source axis '{axis_name}' is not strictly increasing.")
    return x, pprime, ffprime


def load_source_profiles(g, profile_csv_path, profile_x="auto"):
    if profile_csv_path and os.path.exists(profile_csv_path):
        df = pd.read_csv(profile_csv_path)
        if profile_x == "auto":
            if "psi_n" in df.columns:
                axis_name = "psi_n"
            elif "rho_tor" in df.columns:
                axis_name = "rho_tor"
            else:
                raise ValueError(
                    f"Missing profile axis in {profile_csv_path}: require psi_n or rho_tor."
                )
        else:
            axis_name = profile_x
            if axis_name not in df.columns:
                # Compatibility: derive rho_tor axis from psi_n when requested.
                if axis_name == "rho_tor" and "psi_n" in df.columns:
                    rho_of_psin = rho_from_qpsi(g.qpsi, g.psimag, g.psibdy)
                    psin_ref = np.linspace(0.0, 1.0, len(rho_of_psin))
                    src_axis = np.interp(
                        df["psi_n"].to_numpy(dtype=float), psin_ref, rho_of_psin
                    )
                else:
                    raise ValueError(f"Missing '{axis_name}' in {profile_csv_path}")
            else:
                src_axis = df[axis_name].to_numpy(dtype=float)
        if profile_x == "auto":
            src_axis = df[axis_name].to_numpy(dtype=float)

        if "pprime_pa_per_wb" in df.columns:
            pprime = df["pprime_pa_per_wb"].to_numpy(dtype=float)
        else:
            p = df["p_pa"].to_numpy(dtype=float)
            psi_wb = df["psi_wb"].to_numpy(dtype=float)
            pprime = np.gradient(p, psi_wb)

        if "FFprime" in df.columns:
            ffprime = df["FFprime"].to_numpy(dtype=float)
        else:
            f = df["F_tesla_meter"].to_numpy(dtype=float)
            psi_wb = df["psi_wb"].to_numpy(dtype=float)
            ffprime = f * np.gradient(f, psi_wb)

        src_axis, pprime, ffprime = sanitize_axis_and_profiles(
            src_axis, pprime, ffprime, axis_name
        )
        return src_axis, pprime, ffprime, profile_csv_path, axis_name

    if profile_x == "rho_tor":
        src_axis = rho_from_qpsi(g.qpsi, g.psimag, g.psibdy)
        axis_name = "rho_tor"
    else:
        src_axis = g.psi_n
        axis_name = "psi_n"
    src_axis, pprime, ffprime = sanitize_axis_and_profiles(
        src_axis, g.pprime, g.ffprime, axis_name
    )
    return src_axis, pprime, ffprime, "from gfile arrays", axis_name


def load_initial_psi_from_csv(path, r_target, z_target):
    """
    Load initial psi on solver grid [R,Z].

    Supported CSV formats:
      1) flat table with columns: R, Z, PSI
      2) raw matrix saved by this solver (shape [R,Z] or [Z,R])
    """
    df = pd.read_csv(path)

    # Preferred flat format: R,Z,PSI
    need = {"R", "Z", "PSI"}
    if need.issubset(df.columns):
        flat = df[["R", "Z", "PSI"]].astype(float)
        flat = flat[np.isfinite(flat["R"]) & np.isfinite(flat["Z"]) & np.isfinite(flat["PSI"])]
        if flat.empty:
            raise ValueError(f"No valid R,Z,PSI rows in {path}")

        r_src = np.unique(flat["R"].to_numpy(dtype=float))
        z_src = np.unique(flat["Z"].to_numpy(dtype=float))
        if r_src.size < 2 or z_src.size < 2:
            raise ValueError(f"R,Z,PSI CSV has insufficient grid points: {path}")

        piv = flat.pivot_table(index="Z", columns="R", values="PSI", aggfunc="mean")
        piv = piv.reindex(index=z_src, columns=r_src)
        psi_zr = piv.to_numpy(dtype=float)
        if not np.isfinite(psi_zr).all():
            raise ValueError(f"R,Z,PSI CSV does not form a full rectangular grid: {path}")

        same_r = (len(r_src) == len(r_target)) and np.allclose(r_src, r_target, rtol=0.0, atol=1e-10)
        same_z = (len(z_src) == len(z_target)) and np.allclose(z_src, z_target, rtol=0.0, atol=1e-10)
        if same_r and same_z:
            return psi_zr.T, "direct-grid"

        # Bilinear interpolation on rectilinear (R,Z) grid.
        tmp = np.empty((len(z_src), len(r_target)), dtype=float)
        for j in range(len(z_src)):
            tmp[j, :] = np.interp(r_target, r_src, psi_zr[j, :])
        psi_interp_zr = np.empty((len(z_target), len(r_target)), dtype=float)
        for i in range(len(r_target)):
            psi_interp_zr[:, i] = np.interp(z_target, z_src, tmp[:, i])
        return psi_interp_zr.T, "interpolated-grid"

    # Compatibility: plain matrix CSV saved from earlier runs.
    mat = df.to_numpy(dtype=float)
    if mat.shape == (len(r_target), len(z_target)):
        return mat, "matrix-rz"
    if mat.shape == (len(z_target), len(r_target)):
        return mat.T, "matrix-zr"
    raise ValueError(
        f"Unsupported initial psi CSV format or shape for {path}; "
        "expect R,Z,PSI columns or matrix matching target grid."
    )


def interpolate_source(
    psi,
    psi_axis,
    psi_bdy,
    x_src,
    pprime_src,
    ffprime_src,
    source_x_kind="psi_n",
    rho_of_psin=None,
    lcfs_mask=True,
):
    dpsi = psi_bdy - psi_axis
    if abs(dpsi) < 1e-14:
        raise ValueError("psi_bdy == psi_axis, cannot normalize psi.")
    psi_n_raw = (psi - psi_axis) / dpsi
    inside = (psi_n_raw >= 0.0) & (psi_n_raw <= 1.0)
    psi_n = np.clip(psi_n_raw, 0.0, 1.0)
    if source_x_kind == "psi_n":
        x_eval = psi_n
    elif source_x_kind == "rho_tor":
        if rho_of_psin is None:
            raise ValueError("rho_of_psin is required when source_x_kind='rho_tor'.")
        psin_ref = np.linspace(0.0, 1.0, len(rho_of_psin))
        x_eval = np.interp(psi_n, psin_ref, rho_of_psin)
    else:
        raise ValueError(f"Unsupported source_x_kind: {source_x_kind}")
    pprime = np.interp(x_eval, x_src, pprime_src)
    ffprime = np.interp(x_eval, x_src, ffprime_src)
    if lcfs_mask:
        pprime = np.where(inside, pprime, 0.0)
        ffprime = np.where(inside, ffprime, 0.0)
    return pprime, ffprime, inside


def load_pf_profiles(g):
    """
    Load P(psi_n) and F(psi_n) profiles from gfile arrays.

    Returns (psi_n, P_pa, F_tm, dPdpsi_wb, dFdpsi_wb):
      - psi_n: normalized poloidal flux [0,1], shape (nw,)
      - P_pa: pressure in Pa, shape (nw,)
      - F_tm: poloidal current function F = R*Bphi in T*m, shape (nw,)
      - dPdpsi_wb: dP/dpsi in Pa/Wb (derivative w.r.t. physical psi)
      - dFdpsi_wb: dF/dpsi in T*m/Wb
    """
    psi_n = g.psi_n.copy()  # uniform [0,1], size nw
    P_pa = g.pres.copy()    # pressure in Pa
    F_tm = g.fpol.copy()    # F = R*Bphi in T*m

    # Compute derivatives w.r.t. physical psi (Wb)
    dpsi_span = g.psibdy - g.psimag  # psi_bdy - psi_axis in Wb
    # psi_phys = psimag + psi_n * dpsi_span
    # d(anything)/dpsi_phys = d(anything)/dpsi_n / dpsi_span
    dPdpsi_n = np.gradient(P_pa, psi_n)
    dFdpsi_n = np.gradient(F_tm, psi_n)
    if abs(dpsi_span) > 1e-30:
        dPdpsi_wb = dPdpsi_n / dpsi_span
        dFdpsi_wb = dFdpsi_n / dpsi_span
    else:
        dPdpsi_wb = np.zeros_like(psi_n)
        dFdpsi_wb = np.zeros_like(psi_n)

    return psi_n, P_pa, F_tm, dPdpsi_wb, dFdpsi_wb


def interpolate_pf_profiles(psi, psi_axis, psi_bdy, psi_n_src, P_src, F_src,
                            dPdpsi_src, dFdpsi_src, lcfs_mask_psi=True,
                            lcfs_geom_mask=None):
    """
    Interpolate P(psi), F(psi), dP/dpsi, dF/dpsi onto 2D grid from psi field.

    Returns (P_2d, F_2d, dPdpsi_2d, dFdpsi_2d, inside_mask).
    All derivatives are w.r.t. physical psi (Wb).
    """
    dpsi = psi_bdy - psi_axis
    if abs(dpsi) < 1e-14:
        raise ValueError("psi_bdy == psi_axis in interpolate_pf_profiles")
    psi_n_raw = (psi - psi_axis) / dpsi
    inside = (psi_n_raw >= 0.0) & (psi_n_raw <= 1.0)
    psi_n = np.clip(psi_n_raw, 0.0, 1.0)

    P_2d = np.interp(psi_n, psi_n_src, P_src)
    F_2d = np.interp(psi_n, psi_n_src, F_src)
    dPdpsi_2d = np.interp(psi_n, psi_n_src, dPdpsi_src)
    dFdpsi_2d = np.interp(psi_n, psi_n_src, dFdpsi_src)

    if lcfs_mask_psi:
        P_2d = np.where(inside, P_2d, 0.0)
        F_2d = np.where(inside, F_2d, F_src[-1])  # F_vac at boundary
        dPdpsi_2d = np.where(inside, dPdpsi_2d, 0.0)
        dFdpsi_2d = np.where(inside, dFdpsi_2d, 0.0)

    if lcfs_geom_mask is not None:
        P_2d = np.where(lcfs_geom_mask, P_2d, 0.0)
        F_2d = np.where(lcfs_geom_mask, F_2d, F_src[-1])
        dPdpsi_2d = np.where(lcfs_geom_mask, dPdpsi_2d, 0.0)
        dFdpsi_2d = np.where(lcfs_geom_mask, dFdpsi_2d, 0.0)
        inside = inside & lcfs_geom_mask

    return P_2d, F_2d, dPdpsi_2d, dFdpsi_2d, inside


def compute_tj_and_rhs_mdleqf6(R, psi, psi_axis, psi_bdy,
                                psi_n_src, P_src, F_src,
                                dPdpsi_src, dFdpsi_src,
                                Ip_target, F_vac,
                                lcfs_mask_psi=True,
                                lcfs_geom_mask=None,
                                dr=None, dz=None):
    """
    Compute TJ scaling factor and GS right-hand side using MDLEQF=6 algorithm.

    The GS equation source (toroidal current * mu0 * R) is decomposed as:
      mu0*R*j_phi = mu0*R^2*dP/dpsi + F*dF/dpsi / R   (... wait, in the GS eq it's)
      Delta*psi = -mu0*R^2*dP/dpsi - TT*dTT/dpsi
      where TT = F_vac + TJ*(F - F_vac)

    The source decomposes into:
      HJP  = mu0*R^2*dP/dpsi                          (pressure, no TJ)
      HJT1 = F_vac * dF/dpsi * TJ                     (linear in TJ)
      HJT2 = (F - F_vac) * dF/dpsi * TJ^2             (quadratic in TJ)

    Total: rhs = -(HJP + TJ*HJT1_coeff + TJ^2*HJT2_coeff)
    where HJT1_coeff and HJT2_coeff are integrated to get FJT1, FJT2.

    Ip constraint: FJP + TJ*FJT1 + TJ^2*FJT2 = Ip_target
    => quadratic equation for TJ.

    Returns (rhs, TJ, Ip_computed, inside_mask, F_2d, TT_2d).
    """
    P_2d, F_2d, dPdpsi_2d, dFdpsi_2d, inside = interpolate_pf_profiles(
        psi, psi_axis, psi_bdy, psi_n_src, P_src, F_src,
        dPdpsi_src, dFdpsi_src,
        lcfs_mask_psi=lcfs_mask_psi,
        lcfs_geom_mask=lcfs_geom_mask,
    )

    R2d = R[:, None] if R.ndim == 1 else R

    # --- Decompose toroidal current density ---
    # j_phi = R*dP/dpsi + F*dF/dpsi / (mu0*R)
    # But we need to match Ip = integral(j_phi * dA) over the cross-section.
    #
    # With TT = F_vac + TJ*(F - F_vac):
    #   TT * dTT/dpsi = [F_vac + TJ*(F-F_vac)] * TJ * dF/dpsi
    #                 = TJ * F_vac * dF/dpsi + TJ^2 * (F-F_vac) * dF/dpsi
    #
    # GS equation: Delta*psi = -mu0*R^2*dP/dpsi - TT*dTT/dpsi
    #   = -mu0*R^2*dP/dpsi - TJ*F_vac*dF/dpsi - TJ^2*(F-F_vac)*dF/dpsi
    #
    # j_phi * mu0 * R = mu0*R^2*dP/dpsi + TT*dTT/dpsi / ???
    # Actually from GS: -mu0*R*j_phi = Delta*psi (with the sign convention)
    # So: j_phi = -Delta*psi / (mu0*R)
    #          = [mu0*R^2*dP/dpsi + TT*dTT/dpsi] / (mu0*R)
    #          = R*dP/dpsi + TT*dTT/dpsi / (mu0*R)
    #
    # Ip = integral(j_phi dA) = integral over grid cells:
    #    = sum_ij [R_i * dP/dpsi_ij + TT_ij * dTT/dpsi_ij / (mu0 * R_i)] * dr * dz
    #
    # Decompose:
    #   j_phi_P   = R * dP/dpsi                                        (no TJ)
    #   j_phi_T1  = F_vac * dF/dpsi / (mu0 * R)          * TJ         (linear)
    #   j_phi_T2  = (F - F_vac) * dF/dpsi / (mu0 * R)    * TJ^2       (quadratic)

    # Per-cell contributions to Ip:
    hjp_cell = R2d * dPdpsi_2d                              # [A/m^2] pressure
    hjt1_cell = F_vac * dFdpsi_2d / (MU0 * R2d)            # [A/m^2] * TJ
    hjt2_cell = (F_2d - F_vac) * dFdpsi_2d / (MU0 * R2d)   # [A/m^2] * TJ^2

    # Zero outside plasma
    hjp_cell = np.where(inside, hjp_cell, 0.0)
    hjt1_cell = np.where(inside, hjt1_cell, 0.0)
    hjt2_cell = np.where(inside, hjt2_cell, 0.0)

    # Integrate for Ip constraint
    dA = dr * dz
    FJP = float(np.sum(hjp_cell) * dA)
    FJT1 = float(np.sum(hjt1_cell) * dA)
    FJT2 = float(np.sum(hjt2_cell) * dA)

    # Solve quadratic: FJT2*TJ^2 + FJT1*TJ + (FJP - Ip) = 0
    a_coeff = FJT2
    b_coeff = FJT1
    c_coeff = FJP - Ip_target

    TJ = 1.0  # default
    discriminant = b_coeff ** 2 - 4.0 * a_coeff * c_coeff
    if abs(a_coeff) < 1e-30:
        # Linear case
        if abs(b_coeff) > 1e-30:
            TJ = -c_coeff / b_coeff
        else:
            TJ = 1.0
            print(f"  warning: degenerate TJ equation, FJT1={FJT1:.3e}, FJT2={FJT2:.3e}")
    elif discriminant < 0:
        # No real solution -- use TJ=1 as fallback
        TJ = 1.0
        print(f"  warning: TJ discriminant < 0 ({discriminant:.3e}), using TJ=1.0")
    else:
        sqrt_d = np.sqrt(discriminant)
        if b_coeff > 0:
            TJ = (-b_coeff + sqrt_d) / (2.0 * a_coeff)
        else:
            TJ = (-b_coeff - sqrt_d) / (2.0 * a_coeff)

    # Compute actual TT = F_vac + TJ*(F - F_vac) and the full RHS
    TT_2d = F_vac + TJ * (F_2d - F_vac)
    # TT * dTT/dpsi = TT * TJ * dF/dpsi (since dTT/dpsi = TJ * dF/dpsi)
    # But more precisely following Fortran decomposition:
    #   total HJT = HJP + TJ*HJT1 + TJ^2*HJT2  (this is j_phi per cell)
    jphi_2d = hjp_cell + TJ * hjt1_cell + TJ ** 2 * hjt2_cell

    # GS RHS: Delta*psi = -mu0 * R * j_phi
    rhs = -MU0 * R2d * jphi_2d

    Ip_computed = float(np.sum(jphi_2d) * dA)

    return rhs, TJ, Ip_computed, inside, F_2d, TT_2d


def locate_axis_from_psi(psi, r, z, mask=None, dpsi_sign=1.0):
    arr = np.asarray(psi, dtype=float)
    if mask is not None:
        dom = np.where(mask, arr, np.nan)
    else:
        dom = arr
    if not np.any(np.isfinite(dom)):
        raise ValueError("No finite psi values for axis localization.")
    if dpsi_sign >= 0.0:
        flat_idx = int(np.nanargmin(dom))
    else:
        flat_idx = int(np.nanargmax(dom))
    i_ax, j_ax = np.unravel_index(flat_idx, dom.shape)
    return float(r[i_ax]), float(z[j_ax]), float(arr[i_ax, j_ax])


def solve_gs_fixed_boundary(
    g,
    x_src,
    pprime_src,
    ffprime_src,
    source_x_kind="psi_n",
    rho_of_psin=None,
    max_outer=20,
    max_inner=400,
    omega=0.6,
    tol_inner=1e-6,
    freeze_source=False,
    mask_mode="both",
    source_scale=1.0,
    boundary_mode="rectangle",
    ghost_lambda_min=0.2,
    init_psi=None,
    recompute_axis=False,
    # --- MDLEQF=6 parameters ---
    mdleqf6=False,
    ip_target=None,       # Ip in Amperes (if None, use g.currentA)
    pf_profiles=None,     # tuple (psi_n, P, F, dPdpsi, dFdpsi) from load_pf_profiles
):
    R = g.R
    Z = g.Z
    dr = float(R[1] - R[0])
    dz = float(Z[1] - Z[0])

    if init_psi is None:
        psi = g.psirz.copy()  # default initialize from gfile
    else:
        psi = np.asarray(init_psi, dtype=float).copy()
        if psi.shape != g.psirz.shape:
            raise ValueError(
                f"init_psi shape {psi.shape} does not match gfile grid {g.psirz.shape}"
            )
    psi_bc = g.psirz.copy()

    Ri = R[1:-1, None]
    # Delta*psi = d2psi/dR2 - (1/R)dpsi/dR + d2psi/dZ2
    # Coefficients on (i-1) and (i+1) must keep the minus sign on first derivative term.
    A = 1.0 / dr**2 + 1.0 / (2.0 * Ri * dr)  # west  (i-1)
    B = 1.0 / dr**2 - 1.0 / (2.0 * Ri * dr)  # east  (i+1)
    C = 1.0 / dz**2
    D = A + B + 2.0 * C
    psi_span = abs(g.psibdy - g.psimag)
    psi_lo = min(g.psibdy, g.psimag) - 10.0 * psi_span
    psi_hi = max(g.psibdy, g.psimag) + 10.0 * psi_span

    use_psi_mask = mask_mode in ("psi", "both")
    use_lcfs_mask = mask_mode in ("lcfs", "both")
    src_geom_mask = None
    if use_lcfs_mask:
        src_geom_mask = build_lcfs_polygon_mask(R, Z, g.rbbbs, g.zbbbs)
        if src_geom_mask is None:
            print("warning: LCFS polygon unavailable, disable source lcfs mask.")
            use_lcfs_mask = False
        else:
            src_geom_mask = fill_single_cell_holes(src_geom_mask)
            src_geom_mask = largest_connected_component(src_geom_mask)

    if boundary_mode == "rsu_polygon":
        boundary_mode = "lcfs_curve"

    if boundary_mode in ("lcfs", "lcfs_curve"):
        dom_mask = build_lcfs_polygon_mask(R, Z, g.rbbbs, g.zbbbs)
        if dom_mask is None:
            raise ValueError("LCFS boundary mode requested but LCFS polygon is unavailable.")
        dom_mask = fill_single_cell_holes(dom_mask)
        dom_mask = largest_connected_component(dom_mask)
        update_mask, fixed_mask, outside_mask = build_lcfs_update_masks(dom_mask)
        psi_bc_lcfs = np.full_like(psi, g.psibdy)
        if boundary_mode == "lcfs":
            # Ring-Dirichlet mode.
            psi[fixed_mask] = psi_bc_lcfs[fixed_mask]
            boundary_idx = None
            lam_w = lam_e = lam_s = lam_n = None
        else:
            # Curve-Dirichlet mode with ghost-cell interpolation.
            boundary_idx = np.argwhere(fixed_mask)
            lam_w, lam_e, lam_s, lam_n = precompute_lcfs_edge_fractions(
                R, Z, dom_mask, g.rbbbs, g.zbbbs, lambda_min=ghost_lambda_min
            )
        print(
            f"LCFS/RSU Dirichlet domain: inside_frac={np.mean(dom_mask):.3f}, "
            f"update_frac={np.mean(update_mask):.3f}, ring_frac={np.mean(fixed_mask):.3f}"
        )
    else:
        dom_mask = None
        update_mask = None
        fixed_mask = None
        outside_mask = None
        psi_bc_lcfs = None
        boundary_idx = None
        lam_w = lam_e = lam_s = lam_n = None

    compare_mask = dom_mask if dom_mask is not None else np.ones_like(psi, dtype=bool)
    dpsi_sign = 1.0 if float(g.psibdy - g.psimag) >= 0.0 else -1.0

    axis_r = float(g.rmaxis)
    axis_z = float(g.zmaxis)
    axis_psi = float(g.psimag)
    if recompute_axis:
        axis_r, axis_z, axis_psi = locate_axis_from_psi(
            psi, R, Z, mask=compare_mask, dpsi_sign=dpsi_sign
        )
        print(
            "axis mode: dynamic (from current psi), "
            f"initial axis=({axis_r:.5f},{axis_z:.5f}), psi_axis={axis_psi:.6e}"
        )

    # --- MDLEQF=6 setup ---
    if mdleqf6:
        if pf_profiles is None:
            raise ValueError("mdleqf6 mode requires pf_profiles from load_pf_profiles()")
        psi_n_pf, P_pf, F_pf, dPdpsi_pf, dFdpsi_pf = pf_profiles
        F_vac = float(F_pf[-1])  # F at boundary = R0*B0
        if ip_target is None:
            ip_target = float(g.currentA)

        # Rescale init psi to gfile psi range.
        # The init CSV may come from TASK/EQ which uses a 2*pi normalization
        # or other scaling. Without rescaling, psi_n = (psi-psi_axis)/dpsi
        # maps most points outside [0,1], yielding tiny inside_frac and
        # forcing TJ to extreme values to match Ip.
        if init_psi is not None:
            dpsi_gf = g.psibdy - g.psimag
            psi_n_init = (psi - g.psimag) / dpsi_gf if abs(dpsi_gf) > 1e-14 else psi * 0.0
            init_inside = (psi_n_init >= 0.0) & (psi_n_init <= 1.0)
            init_frac = float(np.mean(init_inside))
            if init_frac < 0.10:
                # Rescale: map init psi extremes to gfile psi extremes
                # Use the LCFS domain to find the psi range that matters
                if dom_mask is not None:
                    psi_in = psi[dom_mask]
                else:
                    psi_in = psi.ravel()
                psi_min_init = float(np.min(psi_in))
                psi_max_init = float(np.max(psi_in))
                psi_min_gf = float(np.min(g.psirz[dom_mask])) if dom_mask is not None else float(g.psimag)
                psi_max_gf = float(np.max(g.psirz[dom_mask])) if dom_mask is not None else float(g.psibdy)
                denom = psi_max_init - psi_min_init
                if abs(denom) > 1e-14:
                    psi[:] = psi_min_gf + (psi - psi_min_init) / denom * (psi_max_gf - psi_min_gf)
                    psi_n_after = (psi - g.psimag) / dpsi_gf
                    frac_after = float(np.mean((psi_n_after >= 0) & (psi_n_after <= 1)))
                    print(
                        f"MDLEQF=6: rescaled init psi to gfile range "
                        f"(inside_frac {init_frac:.3f} -> {frac_after:.3f})"
                    )

        print(f"MDLEQF=6 mode: Ip_target={ip_target:.6e} A, F_vac={F_vac:.6f} T*m")

    hist = []
    if freeze_source:
        pprime0, ffprime0, inside0 = interpolate_source(
            psi,
            axis_psi,
            g.psibdy,
            x_src,
            pprime_src,
            ffprime_src,
            source_x_kind=source_x_kind,
            rho_of_psin=rho_of_psin,
            lcfs_mask=use_psi_mask,
        )
        if use_lcfs_mask:
            pprime0 = np.where(src_geom_mask, pprime0, 0.0)
            ffprime0 = np.where(src_geom_mask, ffprime0, 0.0)
            inside0 = inside0 & src_geom_mask
        rhs0 = source_scale * (-MU0 * (R[:, None] ** 2) * pprime0 - ffprime0)
        inside_frac0 = float(np.mean(inside0))
    else:
        rhs0 = None
        inside_frac0 = float("nan")

    TJ_current = 1.0  # track TJ across iterations (for mdleqf6)
    Ip_current = float(g.currentA) if not mdleqf6 else ip_target
    old_full_outer = None  # psi snapshot at start of each outer iteration

    for k in range(max_outer):
        old_full_outer = psi.copy()
        if recompute_axis:
            axis_r, axis_z, axis_psi = locate_axis_from_psi(
                psi, R, Z, mask=compare_mask, dpsi_sign=dpsi_sign
            )

        if mdleqf6:
            # MDLEQF=6: compute TJ and RHS from P(psi), F(psi) + Ip constraint
            rhs, TJ_current, Ip_current, inside, _F2d, _TT2d = \
                compute_tj_and_rhs_mdleqf6(
                    R, psi, axis_psi, g.psibdy,
                    psi_n_pf, P_pf, F_pf, dPdpsi_pf, dFdpsi_pf,
                    ip_target, F_vac,
                    lcfs_mask_psi=use_psi_mask,
                    lcfs_geom_mask=src_geom_mask if use_lcfs_mask else None,
                    dr=dr, dz=dz,
                )
            rhs = source_scale * rhs
            inside_frac = float(np.mean(inside))
        elif rhs0 is None:
            pprime, ffprime, inside = interpolate_source(
                psi,
                axis_psi,
                g.psibdy,
                x_src,
                pprime_src,
                ffprime_src,
                source_x_kind=source_x_kind,
                rho_of_psin=rho_of_psin,
                lcfs_mask=use_psi_mask,
            )
            if use_lcfs_mask:
                pprime = np.where(src_geom_mask, pprime, 0.0)
                ffprime = np.where(src_geom_mask, ffprime, 0.0)
                inside = inside & src_geom_mask
            rhs = source_scale * (-MU0 * (R[:, None] ** 2) * pprime - ffprime)
            inside_frac = float(np.mean(inside))
        else:
            rhs = rhs0
            inside_frac = inside_frac0

        max_change = 0.0
        for _ in range(max_inner):
            old_full = psi.copy()
            old = old_full[1:-1, 1:-1].copy()
            cand = (
                A * old_full[:-2, 1:-1]
                + B * old_full[2:, 1:-1]
                + C * old_full[1:-1, :-2]
                + C * old_full[1:-1, 2:]
                - rhs[1:-1, 1:-1]
            ) / D
            mixed = (1.0 - omega) * old + omega * cand
            if boundary_mode == "lcfs":
                um = update_mask[1:-1, 1:-1]
                cur = psi[1:-1, 1:-1]
                psi[1:-1, 1:-1] = np.where(um, mixed, cur)
                psi[fixed_mask] = psi_bc_lcfs[fixed_mask]
                diff = np.abs(psi[1:-1, 1:-1] - old)
                step_change = float(np.max(diff[um])) if np.any(um) else 0.0
            elif boundary_mode == "lcfs_curve":
                um = update_mask[1:-1, 1:-1]
                cur = psi[1:-1, 1:-1]
                psi[1:-1, 1:-1] = np.where(um, mixed, cur)

                # Boundary-band update using ghost-cell values from LCFS Dirichlet.
                max_bnd = 0.0
                for i, j in boundary_idx:
                    if i <= 0 or i >= len(R) - 1 or j <= 0 or j >= len(Z) - 1:
                        continue
                    p0 = old_full[i, j]
                    ri = R[i]
                    ai = 1.0 / dr**2 + 1.0 / (2.0 * ri * dr)  # west
                    bi = 1.0 / dr**2 - 1.0 / (2.0 * ri * dr)  # east
                    di = ai + bi + 2.0 * C

                    # West/east values (ghost if neighbor outside LCFS).
                    if dom_mask[i - 1, j]:
                        wv = old_full[i - 1, j]
                    else:
                        lw = lam_w[i, j] if np.isfinite(lam_w[i, j]) else 0.5
                        wv = (g.psibdy - (1.0 - lw) * p0) / lw
                    if dom_mask[i + 1, j]:
                        ev = old_full[i + 1, j]
                    else:
                        le = lam_e[i, j] if np.isfinite(lam_e[i, j]) else 0.5
                        ev = (g.psibdy - (1.0 - le) * p0) / le

                    # South/north values (ghost if neighbor outside LCFS).
                    if dom_mask[i, j - 1]:
                        sv = old_full[i, j - 1]
                    else:
                        ls = lam_s[i, j] if np.isfinite(lam_s[i, j]) else 0.5
                        sv = (g.psibdy - (1.0 - ls) * p0) / ls
                    if dom_mask[i, j + 1]:
                        nv = old_full[i, j + 1]
                    else:
                        ln = lam_n[i, j] if np.isfinite(lam_n[i, j]) else 0.5
                        nv = (g.psibdy - (1.0 - ln) * p0) / ln

                    # Stability guard for cut-cell ghost values.
                    if not np.isfinite(wv):
                        wv = g.psibdy
                    if not np.isfinite(ev):
                        ev = g.psibdy
                    if not np.isfinite(sv):
                        sv = g.psibdy
                    if not np.isfinite(nv):
                        nv = g.psibdy
                    wv = float(np.clip(wv, psi_lo, psi_hi))
                    ev = float(np.clip(ev, psi_lo, psi_hi))
                    sv = float(np.clip(sv, psi_lo, psi_hi))
                    nv = float(np.clip(nv, psi_lo, psi_hi))

                    cand_ij = (ai * wv + bi * ev + C * sv + C * nv - rhs[i, j]) / di
                    new_ij = (1.0 - omega) * p0 + omega * cand_ij
                    if not np.isfinite(new_ij):
                        new_ij = p0
                    new_ij = float(np.clip(new_ij, psi_lo, psi_hi))
                    psi[i, j] = new_ij
                    max_bnd = max(max_bnd, abs(new_ij - p0))

                # Keep outside frozen.
                psi[outside_mask] = old_full[outside_mask]
                # Keep LCFS ring pinned to Dirichlet value as a stability guard.
                # Ghost-cell values are still used in boundary-neighbor stencil updates.
                psi[fixed_mask] = psi_bc_lcfs[fixed_mask]

                diff = np.abs(psi[1:-1, 1:-1] - old)
                core_change = float(np.max(diff[um])) if np.any(um) else 0.0
                step_change = max(core_change, max_bnd)
            else:
                psi[1:-1, 1:-1] = mixed
                # Keep Dirichlet boundaries fixed from gfile rectangle boundary.
                psi[0, :] = psi_bc[0, :]
                psi[-1, :] = psi_bc[-1, :]
                psi[:, 0] = psi_bc[:, 0]
                psi[:, -1] = psi_bc[:, -1]
                step_change = float(np.max(np.abs(psi[1:-1, 1:-1] - old)))

            max_change = max(max_change, step_change)
            if step_change < tol_inner:
                break

        err = psi - g.psirz
        e = err[compare_mask]
        rmse = float(np.sqrt(np.mean(e**2)))
        max_abs = float(np.max(np.abs(e)))
        if fixed_mask is not None and np.any(fixed_mask):
            ring_delta = psi[fixed_mask] - g.psibdy
            lcfs_bc_max = float(np.max(np.abs(ring_delta)))
            lcfs_bc_rms = float(np.sqrt(np.mean(ring_delta**2)))
        else:
            lcfs_bc_max = float("nan")
            lcfs_bc_rms = float("nan")
        outer_delta_now = float(np.max(np.abs(
            psi[compare_mask] - old_full_outer[compare_mask]
        )))
        hist.append((k + 1, max_change, rmse, max_abs, TJ_current, inside_frac, outer_delta_now))
        tj_str = f", TJ={TJ_current:.6f}, Ip={Ip_current:.3e}" if mdleqf6 else ""
        print(
            f"[outer {k+1:03d}] inner_chg={max_change:.3e}, outer_chg={outer_delta_now:.3e}, "
            f"rmse={rmse:.3e}, maxerr={max_abs:.3e}, "
            f"in_frac={inside_frac:.3f}, "
            f"axis=({axis_r:.4f},{axis_z:.4f})"
            f"{tj_str}"
        )

        # Stop when nonlinear outer loop is stable enough.
        # Use the total psi change across this entire outer iteration
        # (not just the last inner step) to avoid premature exit when
        # max_inner is large enough for the inner loop to fully converge.
        if outer_delta_now < max(tol_inner * 5.0, 1e-8):
            break

    return psi, np.asarray(hist, dtype=float), compare_mask


def save_outputs(
    out_prefix,
    g,
    psi_sol,
    history,
    profile_source_desc,
    compare_mask=None,
    boundary_mode="rectangle",
    jphi_map=None,
    jphi_gfile_map=None,
    pf_profiles=None,
):
    out_prefix = os.path.abspath(out_prefix)
    out_dir = os.path.dirname(out_prefix)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    # Save solved psi matrix in [R,Z] index order.
    df_psi = pd.DataFrame(psi_sol)
    psi_csv = f"{out_prefix}_psirz_solved.csv"
    df_psi.to_csv(psi_csv, index=False)

    # Save solver history.
    hist_cols = ["outer_iter", "inner_max_change", "rmse_vs_gfile", "max_abs_vs_gfile",
                 "TJ", "inside_frac", "outer_delta"]
    df_h = pd.DataFrame(history, columns=hist_cols[:history.shape[1]] if history.ndim == 2 else hist_cols)
    hist_csv = f"{out_prefix}_solver_history.csv"
    df_h.to_csv(hist_csv, index=False)

    # Plot summary.
    err = psi_sol - g.psirz
    if compare_mask is None:
        compare_mask = np.ones_like(err, dtype=bool)
    err_plot = np.where(compare_mask, err, np.nan)
    # Use ij-indexed mesh to match internal array layout [R,Z] exactly.
    r2d, z2d = np.meshgrid(g.R, g.Z, indexing="ij")

    # Locate solved magnetic axis from psi extremum inside compare domain.
    psi_dom = np.where(compare_mask, psi_sol, np.nan)
    if float(g.psibdy - g.psimag) >= 0.0:
        flat_idx = int(np.nanargmin(psi_dom))
    else:
        flat_idx = int(np.nanargmax(psi_dom))
    i_ax, j_ax = np.unravel_index(flat_idx, psi_dom.shape)
    r_axis_sol = float(g.R[i_ax])
    z_axis_sol = float(g.Z[j_ax])
    vmax = np.nanmax(np.abs(err_plot))
    if not np.isfinite(vmax) or vmax <= 0:
        vmax = 1.0

    fig, axes = plt.subplots(2, 4, figsize=(24, 10), constrained_layout=True)

    cs0 = axes[0, 0].contourf(r2d, z2d, g.psirz, levels=40, cmap="jet")
    fig.colorbar(cs0, ax=axes[0, 0])
    axes[0, 0].plot(g.rbbbs, g.zbbbs, "w--", lw=1.4, label="LCFS (gfile)")
    axes[0, 0].set_title("Gfile psi(R,Z)")
    axes[0, 0].set_xlabel("R (m)")
    axes[0, 0].set_ylabel("Z (m)")
    axes[0, 0].set_aspect("equal")
    axes[0, 0].legend(loc="upper right", fontsize=8)

    cs2 = axes[0, 2].contourf(
        r2d,
        z2d,
        err_plot,
        levels=40,
        cmap="coolwarm",
        vmin=-vmax,
        vmax=vmax,
    )
    fig.colorbar(cs2, ax=axes[0, 2])
    axes[0, 2].plot(g.rbbbs, g.zbbbs, "k--", lw=1.2, alpha=0.9, label="LCFS (gfile)")
    axes[0, 2].set_title("psi_solved - psi_gfile")
    axes[0, 2].set_xlabel("R (m)")
    axes[0, 2].set_ylabel("Z (m)")
    axes[0, 2].set_aspect("equal")
    axes[0, 2].legend(loc="upper right", fontsize=8)

    # Summary text panel (moved to top-right).
    ax_txt = axes[0, 3]
    ax_txt.axis("off")
    if history.size > 0:
        last = history[-1]
        txt = (
            f"Final Iteration: {int(last[0])}\n"
            f"Inner Max Change: {last[1]:.3e}\n"
            f"RMSE vs gfile: {last[2]:.3e}\n"
            f"MAX|err| vs gfile: {last[3]:.3e}\n"
            f"psi_axis: {g.psimag:.6e}\n"
            f"psi_bdy:  {g.psibdy:.6e}\n"
            f"axis(gfile): ({g.rmaxis:.4f}, {g.zmaxis:.4f})\n"
            f"axis(solved): ({r_axis_sol:.4f}, {z_axis_sol:.4f})"
        )
    else:
        txt = "No history data."
    ax_txt.text(0.02, 0.98, txt, va="top", ha="left", fontsize=10, family="monospace")
    ax_txt.set_title("Summary Stats")

    # Keep solved psi panel but move to bottom-left.
    cs1 = axes[1, 0].contourf(r2d, z2d, psi_sol, levels=40, cmap="jet")
    fig.colorbar(cs1, ax=axes[1, 0])
    axes[1, 0].plot(g.rbbbs, g.zbbbs, "w--", lw=1.4, label="LCFS (gfile)")
    axes[1, 0].set_title("Solved psi(R,Z)")
    axes[1, 0].set_xlabel("R (m)")
    axes[1, 0].set_ylabel("Z (m)")
    axes[1, 0].set_aspect("equal")
    axes[1, 0].legend(loc="upper right", fontsize=8)

    # Overlay contour comparison is now embedded in summary.
    ax_ov = axes[0, 1]
    psi_min = min(float(g.psimag), float(g.psibdy))
    psi_max = max(float(g.psimag), float(g.psibdy))
    # Use interior levels to avoid singular contour near axis/boundary.
    levels = np.linspace(psi_min, psi_max, 10)[1:-1]
    if len(levels) >= 2:
        ax_ov.contour(
            r2d,
            z2d,
            g.psirz,
            levels=levels,
            colors="tab:blue",
            linewidths=1.4,
            linestyles="-",
        )
        ax_ov.contour(
            r2d,
            z2d,
            psi_sol,
            levels=levels,
            colors="tab:red",
            linewidths=1.2,
            linestyles="--",
        )

    # Explicit LCFS-level contour from solved psi to visualize effective boundary.
    # IMPORTANT: restrict to compare_mask/domain to avoid spurious psi=psibdy contours
    # in outside/frozen regions (especially for lcfs_curve mode).
    psi_plot = np.where(compare_mask, psi_sol, np.nan)
    pmin = float(np.nanmin(psi_plot))
    pmax = float(np.nanmax(psi_plot))
    has_solved_lcfs = pmin <= float(g.psibdy) <= pmax
    if has_solved_lcfs and boundary_mode in ("lcfs", "lcfs_curve", "rsu_polygon"):
        ax_ov.contour(
            r2d,
            z2d,
            psi_plot,
            levels=[float(g.psibdy)],
            colors="tab:green",
            linewidths=1.5,
            linestyles="-.",
        )

    ax_ov.plot(g.rbbbs, g.zbbbs, "k-", lw=1.2, alpha=0.9, label="LCFS (gfile)")
    ax_ov.plot([g.rmaxis], [g.zmaxis], marker="*", color="tab:blue", ms=10, label="Axis (gfile)")
    ax_ov.plot(
        [r_axis_sol],
        [z_axis_sol],
        marker="*",
        color="tab:orange",
        ms=10,
        label="Axis (solved)",
    )
    ax_ov.set_title("Flux-Surface Contour Overlay (gfile vs solved)")
    ax_ov.set_xlabel("R (m)")
    ax_ov.set_ylabel("Z (m)")
    ax_ov.set_aspect("equal")
    from matplotlib.lines import Line2D

    legend_handles = [
        Line2D([0], [0], color="tab:blue", lw=1.4, ls="-", label="gfile psi contours"),
        Line2D([0], [0], color="tab:red", lw=1.2, ls="--", label="solved psi contours"),
        Line2D([0], [0], color="k", lw=1.2, ls="-", label="LCFS (gfile)"),
        Line2D([0], [0], marker="*", color="tab:blue", lw=0, ms=10, label="Axis (gfile)"),
        Line2D([0], [0], marker="*", color="tab:orange", lw=0, ms=10, label="Axis (solved)"),
    ]
    if has_solved_lcfs and boundary_mode in ("lcfs", "lcfs_curve", "rsu_polygon"):
        legend_handles.append(
            Line2D([0], [0], color="tab:green", lw=1.5, ls="-.", label="LCFS from solved psi=psibdy")
        )
    ax_ov.legend(handles=legend_handles, loc="best", fontsize=9)

    # Convergence panel at bottom-middle-left.
    axes[1, 1].plot(history[:, 0], history[:, 1], "o-", label="inner max change")
    axes[1, 1].plot(history[:, 0], history[:, 2], "s-", label="RMSE vs gfile")
    axes[1, 1].plot(history[:, 0], history[:, 3], "^-", label="MAX|err| vs gfile")
    axes[1, 1].set_yscale("log")
    axes[1, 1].set_xlabel("outer iteration")
    axes[1, 1].set_title("Convergence History")
    axes[1, 1].grid(alpha=0.3)
    axes[1, 1].legend()

    # 2D jphi map panel for visualizing in/outside source behavior.
    ax_j = axes[1, 2]
    if jphi_map is not None:
        j_plot = np.where(compare_mask, jphi_map, np.nan)
        jv = np.nanmax(np.abs(j_plot))
        if not np.isfinite(jv) or jv <= 0:
            jv = 1.0
        csj = ax_j.contourf(
            r2d,
            z2d,
            j_plot / 1e6,
            levels=40,
            cmap="coolwarm",
            vmin=-jv / 1e6,
            vmax=jv / 1e6,
        )
        fig.colorbar(csj, ax=ax_j, label="jphi (MA/m^2)")
        ax_j.plot(g.rbbbs, g.zbbbs, "k-", lw=1.2, alpha=0.9, label="LCFS (gfile)")
        ax_j.set_title("jphi(R,Z) inside compare domain")
        ax_j.set_xlabel("R (m)")
        ax_j.set_ylabel("Z (m)")
        ax_j.set_aspect("equal")
        ax_j.legend(loc="upper right", fontsize=8)
    else:
        ax_j.axis("off")
        ax_j.set_title("jphi(R,Z) not available")

    # 2D jphi map from gfile psi (same profile convention/mask as solve input).
    ax_jg = axes[1, 3]
    if jphi_gfile_map is not None:
        jg_plot = np.where(compare_mask, jphi_gfile_map, np.nan)
        jgv = np.nanmax(np.abs(jg_plot))
        if not np.isfinite(jgv) or jgv <= 0:
            jgv = 1.0
        csjg = ax_jg.contourf(
            r2d,
            z2d,
            jg_plot / 1e6,
            levels=40,
            cmap="coolwarm",
            vmin=-jgv / 1e6,
            vmax=jgv / 1e6,
        )
        fig.colorbar(csjg, ax=ax_jg, label="jphi (MA/m^2)")
        ax_jg.plot(g.rbbbs, g.zbbbs, "k-", lw=1.2, alpha=0.9, label="LCFS (gfile)")
        ax_jg.set_title("jphi_gfile(R,Z) inside compare domain")
        ax_jg.set_xlabel("R (m)")
        ax_jg.set_ylabel("Z (m)")
        ax_jg.set_aspect("equal")
        ax_jg.legend(loc="upper right", fontsize=8)
    else:
        ax_jg.axis("off")
        ax_jg.set_title("jphi_gfile(R,Z) not available")

    fig.suptitle(f"GS Solve From Gfile Profiles ({profile_source_desc})")
    fig_png = f"{out_prefix}_summary.png"
    fig.savefig(fig_png, dpi=160)
    plt.close(fig)

    return psi_csv, hist_csv, fig_png


def save_tj_comparison(out_prefix, g, history, pf_profiles=None):
    """
    Save a dedicated MDLEQF=6 diagnostic figure:
      - TJ convergence history (should approach 1 for self-consistent gfile)
      - TT(psi_n) = F_vac + TJ*(F - F_vac)  vs  F_gfile(psi_n)
    """
    out_prefix = os.path.abspath(out_prefix)
    out_dir = os.path.dirname(out_prefix)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    hist = np.asarray(history, dtype=float)
    has_tj = hist.ndim == 2 and hist.shape[1] >= 5

    if not has_tj:
        return None

    iters = hist[:, 0]
    tj_hist = hist[:, 4]
    inside_hist = hist[:, 5] if hist.shape[1] >= 6 else None
    outer_delta_hist = hist[:, 6] if hist.shape[1] >= 7 else None

    fig, axes = plt.subplots(2, 2, figsize=(14, 10), constrained_layout=True)

    # --- Panel (0,0): TJ convergence ---
    ax = axes[0, 0]
    ax.plot(iters, tj_hist, "o-", color="tab:blue", lw=2, ms=4, label="TJ (solver)")
    ax.axhline(1.0, color="tab:red", ls="--", lw=1.5, label="TJ=1 (gfile reference)")
    ax.set_xlabel("outer iteration")
    ax.set_ylabel("TJ")
    ax.set_title("TJ Convergence (MDLEQF=6)")
    ax.legend()
    ax.grid(alpha=0.3)

    # --- Panel (0,1): TT(psi_n) vs F_gfile(psi_n) ---
    ax = axes[0, 1]
    psi_n = g.psi_n
    F_gfile = g.fpol
    F_vac = float(F_gfile[-1])
    final_tj = float(tj_hist[-1])
    if pf_profiles is not None:
        psi_n_pf, P_pf, F_pf, _, _ = pf_profiles
        TT_final = F_vac + final_tj * (F_pf - F_vac)
        TT_tj1 = F_pf.copy()  # TJ=1 means TT = F_input exactly
        ax.plot(psi_n_pf, F_gfile[:len(psi_n_pf)] if len(F_gfile) == len(psi_n_pf) else
                np.interp(psi_n_pf, psi_n, F_gfile),
                "k-", lw=2.0, label="F gfile (reference)")
        ax.plot(psi_n_pf, TT_final, "r--", lw=2.0,
                label=f"TT = F_vac + TJ*dF (TJ={final_tj:.2f})")
        ax.plot(psi_n_pf, TT_tj1, "b:", lw=1.5, label="TT at TJ=1 (=F input)")
    else:
        TT_final_simple = F_vac + final_tj * (F_gfile - F_vac)
        ax.plot(psi_n, F_gfile, "k-", lw=2.0, label="F gfile (reference)")
        ax.plot(psi_n, TT_final_simple, "r--", lw=2.0,
                label=f"TT = F_vac + TJ*dF (TJ={final_tj:.2f})")
    ax.axhline(F_vac, color="gray", ls=":", lw=1, alpha=0.6, label=f"F_vac={F_vac:.3f}")
    ax.set_xlabel("psi_n")
    ax.set_ylabel("F or TT (T*m)")
    ax.set_title("Toroidal Field Function: TT(solved) vs F(gfile)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # --- Panel (1,0): inside_frac history ---
    ax = axes[1, 0]
    if inside_hist is not None:
        ax.plot(iters, inside_hist, "s-", color="tab:green", lw=2, ms=4)
        ax.set_xlabel("outer iteration")
        ax.set_ylabel("inside fraction")
        ax.set_title("Plasma Inside Fraction")
        ax.grid(alpha=0.3)
        ax.set_ylim(bottom=0)
    else:
        ax.axis("off")
        ax.set_title("inside_frac not recorded")

    # --- Panel (1,1): outer_delta and convergence metrics ---
    ax = axes[1, 1]
    ax.plot(iters, hist[:, 2], "s-", color="tab:orange", lw=1.5, ms=3, label="RMSE vs gfile")
    ax.plot(iters, hist[:, 3], "^-", color="tab:red", lw=1.5, ms=3, label="MAX|err| vs gfile")
    if outer_delta_hist is not None:
        ax.plot(iters, outer_delta_hist, "o-", color="tab:blue", lw=1.5, ms=3, label="outer delta")
    ax.set_yscale("log")
    ax.set_xlabel("outer iteration")
    ax.set_title("Convergence Metrics")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    fig.suptitle(f"MDLEQF=6 TJ Diagnostic (final TJ={final_tj:.4f}, gfile TJ=1)")
    fig_png = f"{out_prefix}_tj_comparison.png"
    fig.savefig(fig_png, dpi=160)
    plt.close(fig)
    return fig_png


def save_source_current_profiles(
    out_prefix, g, x_src, pprime_src, ffprime_src, source_x_kind="psi_n", r0=8.03, rho_of_psin=None
):
    """
    Save a dedicated diagnostic figure for source terms and toroidal current profile.
    Profiles are shown using a representative major radius R0 = magnetic-axis R.
    """
    out_prefix = os.path.abspath(out_prefix)
    out_dir = os.path.dirname(out_prefix)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    x = np.asarray(x_src, dtype=float)
    pprime = np.asarray(pprime_src, dtype=float)
    ffprime = np.asarray(ffprime_src, dtype=float)
    if x.size < 2:
        raise ValueError("Not enough source points to plot source/current profiles.")

    r0 = float(r0)
    if not np.isfinite(r0) or r0 <= 0.0:
        raise ValueError("r0 must be finite and > 0 for source/current profile plotting.")

    # GS source decomposition at representative major radius R0:
    # Delta*psi = -mu0*R^2*p' - FF'
    src_p = -MU0 * (r0**2) * pprime
    src_f = -ffprime
    src_tot = src_p + src_f

    # Toroidal current density from GS relation:
    # -mu0*R*jphi = -mu0*R^2*p' - FF' => jphi = R*p' + FF'/(mu0*R)
    jphi = r0 * pprime + ffprime / (MU0 * r0)

    xlab = "psi_n" if source_x_kind == "psi_n" else "rho_tor"

    fig, axes = plt.subplots(2, 3, figsize=(20, 10), constrained_layout=True)

    axes[0, 0].plot(x, src_p, lw=2.0, label="-mu0 R0^2 p'")
    axes[0, 0].plot(x, src_f, lw=2.0, label="-FF'")
    axes[0, 0].plot(x, src_tot, "k--", lw=2.0, label="total source")
    axes[0, 0].set_title(f"GS Source Profiles (R0={r0:.3f} m)")
    axes[0, 0].set_xlabel(xlab)
    axes[0, 0].set_ylabel("source term")
    axes[0, 0].grid(alpha=0.3)
    axes[0, 0].legend()

    axes[0, 1].plot(x, jphi / 1e6, color="tab:red", lw=2.0)
    axes[0, 1].set_title(f"Toroidal Current Density jphi (R0={r0:.3f} m)")
    axes[0, 1].set_xlabel(xlab)
    axes[0, 1].set_ylabel("jphi (MA/m^2)")
    axes[0, 1].grid(alpha=0.3)

    # Reconstruct p and F from p' and FF' assuming boundary values:
    # F_bdy = R0 * B0 (g.bcentr * g.rcentr), p_bdy = 0
    from scipy.integrate import cumulative_trapezoid
    
    # In GS equation, pprime is dp/dpsi and ffprime is F dF/dpsi.
    # We must integrate over true psi (Wb), not psi_n or rho.
    dpsi_span = g.psibdy - g.psimag
    
    if source_x_kind == "psi_n":
        psi_x = g.psimag + x * dpsi_span
        psi_n_x = x
    else:
        psin_ref = np.linspace(0.0, 1.0, len(g.qpsi))
        from scipy.interpolate import interp1d
        rho_to_psin = interp1d(rho_of_psin, psin_ref, fill_value="extrapolate")
        psi_n_x = rho_to_psin(x)
        psi_x = g.psimag + psi_n_x * dpsi_span

    p_profile = cumulative_trapezoid(pprime, psi_x, initial=0.0)
    p_profile = p_profile - p_profile[-1]

    F2_profile = cumulative_trapezoid(2.0 * ffprime, psi_x, initial=0.0)
    F2_bdy = (g.bcentr * g.rcentr) ** 2
    F2_profile = F2_profile - F2_profile[-1] + F2_bdy
    F_profile = np.sqrt(np.maximum(F2_profile, 0.0))

    # Match the sign of F_profile with gfile's fpol (since F = R*Bphi can be negative)
    # Most tokamaks use negative F if Bphi is negative.
    if g.fpol[0] < 0:
        F_profile = -F_profile

    # For gfile arrays, we need to map them to the same x-axis for plotting
    # g.pres and g.fpol are defined on uniform psi_n grid (g.psi_n)
    from scipy.interpolate import interp1d
    pres_interp = interp1d(g.psi_n, g.pres, bounds_error=False, fill_value="extrapolate")
    fpol_interp = interp1d(g.psi_n, g.fpol, bounds_error=False, fill_value="extrapolate")
    
    gfile_p_at_x = pres_interp(psi_n_x)
    gfile_F_at_x = fpol_interp(psi_n_x)

    axes[1, 0].plot(x, p_profile / 1e6, color="tab:blue", lw=2.0, label="Reconstructed ($\int p' d\psi$)")
    axes[1, 0].plot(x, gfile_p_at_x / 1e6, color="tab:orange", lw=2.0, ls="--", label="gfile PRES")
    axes[1, 0].set_title("Pressure p")
    axes[1, 0].set_xlabel(xlab)
    axes[1, 0].set_ylabel("p (MPa)")
    axes[1, 0].grid(alpha=0.3)
    axes[1, 0].legend()

    axes[1, 1].plot(x, F_profile, color="tab:green", lw=2.0, label="Reconstructed ($\int FF' d\psi$)")
    axes[1, 1].plot(x, gfile_F_at_x, color="tab:red", lw=2.0, ls="--", label="gfile FPOL")
    axes[1, 1].set_title("Poloidal Current Function F")
    axes[1, 1].set_xlabel(xlab)
    axes[1, 1].set_ylabel("F = R B_phi (T*m)")
    axes[1, 1].grid(alpha=0.3)
    axes[1, 1].legend()

    # Third column: Explicitly plotted vs psi_n
    axes[0, 2].plot(psi_n_x, p_profile / 1e6, color="tab:blue", lw=2.0, label="Reconstructed")
    axes[0, 2].plot(g.psi_n, g.pres / 1e6, color="tab:orange", lw=2.0, ls="--", label="gfile PRES")
    axes[0, 2].set_title("Pressure p (vs psi_n)")
    axes[0, 2].set_xlabel("psi_n")
    axes[0, 2].set_ylabel("p (MPa)")
    axes[0, 2].grid(alpha=0.3)
    axes[0, 2].legend()

    axes[1, 2].plot(psi_n_x, F_profile, color="tab:green", lw=2.0, label="Reconstructed")
    axes[1, 2].plot(g.psi_n, g.fpol, color="tab:red", lw=2.0, ls="--", label="gfile FPOL")
    axes[1, 2].set_title("Poloidal Current F (vs psi_n)")
    axes[1, 2].set_xlabel("psi_n")
    axes[1, 2].set_ylabel("F = R B_phi (T*m)")
    axes[1, 2].grid(alpha=0.3)
    axes[1, 2].legend()

    fig.suptitle("Input Source and Current Profiles")
    fig_png = f"{out_prefix}_source_current_profiles.png"
    fig.savefig(fig_png, dpi=160)
    plt.close(fig)
    return fig_png


def main():
    parser = argparse.ArgumentParser(description="Solve fixed-boundary GS from gfile profiles.")
    parser.add_argument(
        "--gfile",
        default="in/g260206.20000_teq_0114",
        help="Input gfile path.",
    )
    parser.add_argument(
        "--profile-csv",
        default=None,
        help="Optional *_profiles_from_gfile.csv path. Auto-detect by default.",
    )
    parser.add_argument(
        "--profile-x",
        choices=["auto", "psi_n", "rho_tor"],
        default="auto",
        help="Source-profile axis convention: auto-detect / psi_n / rho_tor.",
    )
    parser.add_argument(
        "--out-prefix",
        default="in/gs_from_gfile",
        help="Output prefix for CSV/PNG.",
    )
    parser.add_argument("--max-outer", type=int, default=20)
    parser.add_argument("--max-inner", type=int, default=400)
    parser.add_argument("--omega", type=float, default=0.6)
    parser.add_argument("--tol-inner", type=float, default=1e-6)
    parser.add_argument(
        "--source-scale",
        type=float,
        default=1.0,
        help="Multiply GS source term by this factor (for unit/normalization diagnosis).",
    )
    parser.add_argument(
        "--freeze-source",
        action="store_true",
        help="Use source term from initial gfile psi only (linear fixed-source solve).",
    )
    parser.add_argument(
        "--no-lcfs-mask",
        action="store_true",
        help="Deprecated shortcut: same as --mask-mode none.",
    )
    parser.add_argument(
        "--mask-mode",
        choices=["none", "psi", "lcfs", "both"],
        default="both",
        help="Source mask mode: none / psi-normalized / LCFS polygon / both.",
    )
    parser.add_argument(
        "--boundary-mode",
        choices=["rectangle", "lcfs", "lcfs_curve", "rsu_polygon"],
        default="rectangle",
        help=(
            "Dirichlet boundary mode: rectangle / LCFS ring / LCFS curve ghost-cell / "
            "RSU polygon (alias of lcfs_curve)."
        ),
    )
    parser.add_argument(
        "--ghost-lambda-min",
        type=float,
        default=0.2,
        help="Minimum boundary intersection fraction for lcfs_curve ghost-cell.",
    )
    parser.add_argument(
        "--r0-profile",
        type=float,
        default=8.03,
        help="Major radius R0 (m) used in source/current profile diagnostics.",
    )
    parser.add_argument(
        "--init-psi-csv",
        default=None,
        help=(
            "Optional initial psi CSV path. Supports flat R,Z,PSI table "
            "(e.g. in/eqdata0114_eq_psirz_rz.csv) or solver matrix CSV."
        ),
    )
    parser.add_argument(
        "--recompute-axis",
        action="store_true",
        help="Recompute psi-axis from current psi each outer iteration.",
    )
    parser.add_argument(
        "--mdleqf6",
        action="store_true",
        help=(
            "Use MDLEQF=6 algorithm: given P(psi) and F(psi) profiles from gfile, "
            "with TJ scaling factor to match Ip constraint."
        ),
    )
    parser.add_argument(
        "--ip-target",
        type=float,
        default=None,
        help=(
            "Target plasma current Ip in Amperes for MDLEQF=6 mode. "
            "If not specified, uses gfile currentA."
        ),
    )
    args = parser.parse_args()

    if not (0.0 < args.omega <= 1.0):
        raise ValueError("For weighted Jacobi, omega must satisfy 0 < omega <= 1.")
    if not (0.0 < args.ghost_lambda_min <= 1.0):
        raise ValueError("--ghost-lambda-min must satisfy 0 < value <= 1.")

    gfile_path = os.path.abspath(args.gfile)
    if not os.path.exists(gfile_path):
        raise FileNotFoundError(f"gfile not found: {gfile_path}")

    g = read_gfile(gfile_path)
    profile_csv = resolve_profile_csv(gfile_path, args.profile_csv)
    x_src, pprime_src, ffprime_src, src_desc, source_x_kind = load_source_profiles(
        g, profile_csv, profile_x=args.profile_x
    )
    rho_of_psin = rho_from_qpsi(g.qpsi, g.psimag, g.psibdy) if source_x_kind == "rho_tor" else None

    # Load P(psi), F(psi) profiles for MDLEQF=6 mode
    pf_profiles = None
    if args.mdleqf6:
        pf_profiles = load_pf_profiles(g)
        src_desc += " [MDLEQF=6: P,F + Ip]"

    init_psi = None
    init_desc = "gfile psirz"
    if args.init_psi_csv:
        init_path = os.path.abspath(args.init_psi_csv)
        if not os.path.exists(init_path):
            raise FileNotFoundError(f"initial psi CSV not found: {init_path}")
        init_psi, init_mode = load_initial_psi_from_csv(init_path, g.R, g.Z)
        init_desc = f"{init_path} ({init_mode})"

    print(f"gfile: {gfile_path}")
    print(f"profile source: {src_desc}")
    print(f"profile axis convention: {source_x_kind}")
    print(f"initial psi: {init_desc}")
    mask_mode = "none" if args.no_lcfs_mask else args.mask_mode
    print(f"mask_mode: {mask_mode}")
    print(f"boundary_mode: {args.boundary_mode}")
    print(f"ghost_lambda_min: {args.ghost_lambda_min}")
    print(f"source_scale: {args.source_scale}")
    print(f"r0_profile: {args.r0_profile}")
    print(f"recompute_axis: {args.recompute_axis}")
    print(f"mdleqf6: {args.mdleqf6}")
    if args.mdleqf6:
        ip_disp = args.ip_target if args.ip_target is not None else g.currentA
        print(f"ip_target: {ip_disp:.6e} A")
    print(
        f"grid: nw={g.nw}, nh={g.nh}, psi_axis={g.psimag:.6e}, "
        f"psi_bdy={g.psibdy:.6e}, Ip={g.currentA:.6e} A"
    )

    psi_sol, hist, compare_mask = solve_gs_fixed_boundary(
        g,
        x_src,
        pprime_src,
        ffprime_src,
        source_x_kind=source_x_kind,
        rho_of_psin=rho_of_psin,
        max_outer=args.max_outer,
        max_inner=args.max_inner,
        omega=args.omega,
        tol_inner=args.tol_inner,
        freeze_source=args.freeze_source,
        mask_mode=mask_mode,
        source_scale=args.source_scale,
        boundary_mode=args.boundary_mode,
        ghost_lambda_min=args.ghost_lambda_min,
        init_psi=init_psi,
        recompute_axis=args.recompute_axis,
        mdleqf6=args.mdleqf6,
        ip_target=args.ip_target,
        pf_profiles=pf_profiles,
    )

    # Build final 2D jphi(R,Z) map for summary visualization.
    use_psi_mask = mask_mode in ("psi", "both")
    use_lcfs_mask = mask_mode in ("lcfs", "both")
    pprime_fin, ffprime_fin, inside_fin = interpolate_source(
        psi_sol,
        g.psimag,
        g.psibdy,
        x_src,
        pprime_src,
        ffprime_src,
        source_x_kind=source_x_kind,
        rho_of_psin=rho_of_psin,
        lcfs_mask=use_psi_mask,
    )
    if use_lcfs_mask:
        src_geom_mask = build_lcfs_polygon_mask(g.R, g.Z, g.rbbbs, g.zbbbs)
        if src_geom_mask is not None:
            src_geom_mask = fill_single_cell_holes(src_geom_mask)
            src_geom_mask = largest_connected_component(src_geom_mask)
            pprime_fin = np.where(src_geom_mask, pprime_fin, 0.0)
            ffprime_fin = np.where(src_geom_mask, ffprime_fin, 0.0)
            inside_fin = inside_fin & src_geom_mask
    jphi_map = g.R[:, None] * pprime_fin + ffprime_fin / (MU0 * g.R[:, None])
    jphi_map = np.where(inside_fin, jphi_map, 0.0)

    # Build 2D jphi(R,Z) map from gfile psi for direct visual comparison.
    pprime_gf, ffprime_gf, inside_gf = interpolate_source(
        g.psirz,
        g.psimag,
        g.psibdy,
        x_src,
        pprime_src,
        ffprime_src,
        source_x_kind=source_x_kind,
        rho_of_psin=rho_of_psin,
        lcfs_mask=use_psi_mask,
    )
    if use_lcfs_mask:
        src_geom_mask = build_lcfs_polygon_mask(g.R, g.Z, g.rbbbs, g.zbbbs)
        if src_geom_mask is not None:
            src_geom_mask = fill_single_cell_holes(src_geom_mask)
            src_geom_mask = largest_connected_component(src_geom_mask)
            pprime_gf = np.where(src_geom_mask, pprime_gf, 0.0)
            ffprime_gf = np.where(src_geom_mask, ffprime_gf, 0.0)
            inside_gf = inside_gf & src_geom_mask
    jphi_gfile_map = g.R[:, None] * pprime_gf + ffprime_gf / (MU0 * g.R[:, None])
    jphi_gfile_map = np.where(inside_gf, jphi_gfile_map, 0.0)

    psi_csv, hist_csv, fig_png = save_outputs(
        args.out_prefix,
        g,
        psi_sol,
        hist,
        src_desc,
        compare_mask=compare_mask,
        boundary_mode=args.boundary_mode,
        jphi_map=jphi_map,
        jphi_gfile_map=jphi_gfile_map,
        pf_profiles=pf_profiles,
    )
    profile_fig_png = save_source_current_profiles(
        args.out_prefix,
        g,
        x_src,
        pprime_src,
        ffprime_src,
        source_x_kind=source_x_kind,
        r0=args.r0_profile,
        rho_of_psin=rho_of_psin,
    )
    print(f"saved: {psi_csv}")
    print(f"saved: {hist_csv}")
    print(f"saved: {fig_png}")
    print(f"saved: {profile_fig_png}")

    # MDLEQF=6 TJ diagnostic figure
    if args.mdleqf6:
        tj_fig = save_tj_comparison(args.out_prefix, g, hist, pf_profiles=pf_profiles)
        if tj_fig:
            print(f"saved: {tj_fig}")


if __name__ == "__main__":
    main()
