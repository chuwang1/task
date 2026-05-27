#!/usr/bin/env python3
"""
Solve fixed-boundary Grad-Shafranov equation from gfile/profile data.

Equation:
  d2psi/dR2 - (1/R) dpsi/dR + d2psi/dZ2 = -mu0 * R^2 * p'(psi) - F(psi)F'(psi)

This script uses:
  - Dirichlet boundary on either rectangular box or LCFS polygon
  - Picard outer iteration for nonlinear source update
  - Weighted Jacobi inner iteration for elliptic solve (0 < omega <= 1)
"""

import argparse
import os
import re
import struct

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.path import Path as MplPath

MU0 = 4.0e-7 * np.pi
DEBUG = False


def dprint(*args, **kwargs):
    if DEBUG:
        print(*args, **kwargs)


class GFile:
    pass


def read_fortran_record(f):
    marker = f.read(4)
    if len(marker) == 0:
        return None
    if len(marker) < 4:
        raise IOError(f"Incomplete record marker: got {len(marker)} bytes")
    record_len = struct.unpack('<i', marker)[0]
    data = f.read(record_len)
    if len(data) < record_len:
        raise IOError(f"Incomplete record data: expected {record_len}, got {len(data)}")
    trailing = f.read(4)
    if len(trailing) < 4:
        raise IOError("Incomplete trailing record marker")
    trailing_len = struct.unpack('<i', trailing)[0]
    if trailing_len != record_len:
        raise IOError(f"Record length mismatch: {record_len} vs {trailing_len}")
    return data


def read_eqdata_binary(filename):
    data = {}
    with open(filename, 'rb') as f:
        rec = read_fortran_record(f)
        vals = np.frombuffer(rec, dtype=np.float64)
        data['RR'], data['BB'], data['RIP'] = vals[:3]

        rec = read_fortran_record(f)
        vals = np.frombuffer(rec, dtype=np.int32)
        nrgmax, nzgmax = int(vals[0]), int(vals[1])
        data['NRGMAX'], data['NZGMAX'] = nrgmax, nzgmax

        data['RG'] = np.frombuffer(read_fortran_record(f), dtype=np.float64)
        data['ZG'] = np.frombuffer(read_fortran_record(f), dtype=np.float64)
        psirz_flat = np.frombuffer(read_fortran_record(f), dtype=np.float64)
        data['PSIRZ'] = psirz_flat.reshape((nrgmax, nzgmax), order='F')

        npsmax = int(np.frombuffer(read_fortran_record(f), dtype=np.int32)[0])
        data['NPSMAX'] = npsmax
        for name in ['PSIPS', 'PPPS', 'TTPS', 'DPPPS', 'DTTPS', 'TTDTTPS', 'QQPS', 'TEPS', 'OMPS']:
            data[name] = np.frombuffer(read_fortran_record(f), dtype=np.float64)

        vals = np.frombuffer(read_fortran_record(f), dtype=np.int32)
        nsgmax, ntgmax, nugmax, nrmax, nthmax, nsumax, nrvmax, ntvmax = [int(v) for v in vals]
        data['NSGMAX'], data['NTGMAX'], data['NSUMAX'], data['NRVMAX'] = nsgmax, ntgmax, nsumax, nrvmax

        for name, shape in [('PSI', (nsgmax, ntgmax)), ('DELPSI', (nsgmax, ntgmax)), ('HJT', (nsgmax, ntgmax))]:
            arr = np.frombuffer(read_fortran_record(f), dtype=np.float64)
            data[name] = arr.reshape(shape, order='F')

        vals = np.frombuffer(read_fortran_record(f), dtype=np.float64)
        data['RAXIS'], data['ZAXIS'], data['PSITA'], data['PSIPA'], data['PSI0'] = vals[:5]

        for name in ['PSIPNV', 'PSIPV', 'PSITV', 'QPV', 'TTV']:
            data[name] = np.frombuffer(read_fortran_record(f), dtype=np.float64)

        vals = np.frombuffer(read_fortran_record(f), dtype=np.float64)
        data['RA'], data['RKAP'], data['RDLT'], data['RB'], data['FRBIN'] = vals[:5]

        # skip 5 parameter blocks
        for _ in range(5):
            read_fortran_record(f)

        try:
            hjtrz_flat = np.frombuffer(read_fortran_record(f), dtype=np.float64)
            data['HJTRZ'] = hjtrz_flat.reshape((nrgmax, nzgmax), order='F')
        except Exception:
            data['HJTRZ'] = np.zeros((nrgmax, nzgmax), dtype=float)

        optional_names = ['RSV', 'VPV', 'AVIR2', 'AVRR2', 'FIPV', 'QIPV', 'DPIPV']
        for name in optional_names:
            try:
                rec = read_fortran_record(f)
                if rec is None:
                    data[name] = np.full(nrvmax, np.nan)
                else:
                    arr = np.frombuffer(rec, dtype=np.float64)
                    data[name] = arr if arr.size == nrvmax else np.full(nrvmax, np.nan)
            except Exception:
                data[name] = np.full(nrvmax, np.nan)

        try:
            rec = read_fortran_record(f)
            if rec is not None:
                nsum_saved = int(np.frombuffer(rec, dtype=np.int32)[0])
                data['NSUMAX_SAVED'] = nsum_saved
                if nsum_saved > 0:
                    data['RSU'] = np.frombuffer(read_fortran_record(f), dtype=np.float64)
                    data['ZSU'] = np.frombuffer(read_fortran_record(f), dtype=np.float64)
                    data['RSW'] = np.frombuffer(read_fortran_record(f), dtype=np.float64)
                    data['ZSW'] = np.frombuffer(read_fortran_record(f), dtype=np.float64)
                else:
                    data['RSU'] = np.asarray([], dtype=float)
                    data['ZSU'] = np.asarray([], dtype=float)
                    data['RSW'] = np.asarray([], dtype=float)
                    data['ZSW'] = np.asarray([], dtype=float)
        except Exception:
            data['NSUMAX_SAVED'] = 0
            data['RSU'] = np.asarray([], dtype=float)
            data['ZSU'] = np.asarray([], dtype=float)
            data['RSW'] = np.asarray([], dtype=float)
            data['ZSW'] = np.asarray([], dtype=float)
    return data


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

    If the input matrix is stored as [R,Z], transpose it for plotting;
    if it is already [Z,R], leave it unchanged.
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


def downsample_gfile(g, nw_new, nh_new):
    from scipy.interpolate import RectBivariateSpline
    
    # Store original grid
    r_old = g.R
    z_old = g.Z
    
    # Create new grid (maintaining physical bounds)
    r_new = np.linspace(g.rleft, g.rleft + g.rdim, nw_new)
    z_new = np.linspace(g.zmid - 0.5 * g.zdim, g.zmid + 0.5 * g.zdim, nh_new)
    
    # Interpolate psirz (which is stored as [R, Z])
    # RectBivariateSpline takes x, y, and z(x,y)
    spline = RectBivariateSpline(r_old, z_old, g.psirz)
    psirz_new = spline(r_new, z_new)
    
    # Interpolate 1D profiles (which are on uniform psi_n grid)
    psi_n_old = g.psi_n
    psi_n_new = np.linspace(0.0, 1.0, nw_new)
    
    # For 1D profiles we can just use np.interp
    fpol_new = np.interp(psi_n_new, psi_n_old, g.fpol)
    pres_new = np.interp(psi_n_new, psi_n_old, g.pres)
    ffprime_new = np.interp(psi_n_new, psi_n_old, g.ffprime)
    pprime_new = np.interp(psi_n_new, psi_n_old, g.pprime)
    qpsi_new = np.interp(psi_n_new, psi_n_old, g.qpsi)
    
    # Update g object
    g.nw = nw_new
    g.nh = nh_new
    g.R = r_new
    g.Z = z_new
    g.psirz = psirz_new
    g.fpol = fpol_new
    g.pres = pres_new
    g.ffprime = ffprime_new
    g.pprime = pprime_new
    g.qpsi = qpsi_new
    g.psi_n = psi_n_new
    
    return g

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


def load_eq_modeg5_profiles(eq_prefix):
    path = f"{eq_prefix}_modeg5_profiles.csv"
    rows = []
    with open(path, "r", encoding="utf-8") as f:
        next(f)
        for line in f:
            left, right = line.split(",", 1)
            vals = [float(x) for x in right.split()]
            rows.append([int(left)] + vals)
    arr = np.asarray(rows, dtype=float)
    psips = arr[:, 1]
    ppps = arr[:, 2]
    ttps = arr[:, 3]
    dppps = arr[:, 4]
    dttps = arr[:, 5]
    ttdttps = arr[:, 6]
    qqps = arr[:, 7]
    return {
        "path": path,
        "psips": psips,
        "ppps": ppps,
        "ttps": ttps,
        "dppps": dppps,
        "dttps": dttps,
        "ttdttps": ttdttps,
        "qqps": qqps,
    }


def derive_gfile_geometry_from_grid(r, z, rcentr=None):
    r = np.asarray(r, dtype=float)
    z = np.asarray(z, dtype=float)
    out = {
        "rdim": float(r[-1] - r[0]),
        "zdim": float(z[-1] - z[0]),
        "rleft": float(r[0]),
        "zmid": 0.5 * float(z[0] + z[-1]),
    }
    if rcentr is not None:
        out["rcentr"] = float(rcentr)
    return out


def derive_gfile_axis_flux_fields(psi, r, z, psi_bdy, bcentr=None, mask=None, dpsi_sign=1.0):
    rmaxis, zmaxis, psimag = locate_axis_from_psi(
        psi, r, z, mask=mask, dpsi_sign=dpsi_sign
    )
    out = {
        "rmaxis": float(rmaxis),
        "zmaxis": float(zmaxis),
        "psimag": float(psimag),
        "psibdy": float(psi_bdy),
    }
    if bcentr is not None:
        out["bcentr"] = float(bcentr)
    return out


def geqdsk_chunks(values, per_line=5):
    vals = [float(v) for v in values]
    for i in range(0, len(vals), per_line):
        chunk = vals[i : i + per_line]
        yield "".join(f"{v:16.9E}" for v in chunk)


def write_gfile(path, g_template, psi_sol, axis_info, lcfs_r=None, lcfs_z=None):
    r = np.asarray(g_template.R, dtype=float)
    z = np.asarray(g_template.Z, dtype=float)
    nw = len(r)
    nh = len(z)
    geom = derive_gfile_geometry_from_grid(r, z, rcentr=getattr(g_template, "rcentr", None))

    rbbbs = np.asarray(lcfs_r if lcfs_r is not None else g_template.rbbbs, dtype=float)
    zbbbs = np.asarray(lcfs_z if lcfs_z is not None else g_template.zbbbs, dtype=float)
    nbbbs = int(len(rbbbs))

    limitr = int(getattr(g_template, "limitr", 0))
    if limitr > 0 and hasattr(g_template, "rlim") and hasattr(g_template, "zlim"):
        rlim = np.asarray(g_template.rlim, dtype=float)
        zlim = np.asarray(g_template.zlim, dtype=float)
    else:
        limitr = 1
        rlim = np.asarray([r[0]], dtype=float)
        zlim = np.asarray([z[0]], dtype=float)

    with open(path, "w", encoding="utf-8") as f:
        header = getattr(g_template, "header", "PYTHON SOLVED GFILE")[:48]
        f.write(f"{header:<48}{3:4d}{nw:4d}{nh:4d}\n")
        f.write(
            "".join(
                f"{v:16.9E}"
                for v in [geom["rdim"], geom["zdim"], geom.get("rcentr", 0.0), geom["rleft"], geom["zmid"]]
            )
            + "\n"
        )
        f.write(
            "".join(
                f"{v:16.9E}"
                for v in [axis_info["rmaxis"], axis_info["zmaxis"], axis_info["psimag"], axis_info["psibdy"], axis_info.get("bcentr", getattr(g_template, "bcentr", 0.0))]
            )
            + "\n"
        )
        f.write(
            "".join(
                f"{v:16.9E}"
                for v in [float(g_template.currentA), axis_info["psimag"], 0.0, axis_info["rmaxis"], 0.0]
            )
            + "\n"
        )
        f.write(
            "".join(
                f"{v:16.9E}"
                for v in [axis_info["zmaxis"], 0.0, axis_info["psibdy"], 0.0, 0.0]
            )
            + "\n"
        )

        for arr in [g_template.fpol, g_template.pres, g_template.ffprime, g_template.pprime]:
            for line in geqdsk_chunks(arr):
                f.write(line + "\n")

        psirz_for_write = np.asarray(psi_sol, dtype=float).T.reshape(-1)
        for line in geqdsk_chunks(psirz_for_write):
            f.write(line + "\n")

        for line in geqdsk_chunks(g_template.qpsi):
            f.write(line + "\n")

        f.write(f"{nbbbs:5d}{limitr:5d}\n")
        bbbs = np.column_stack([rbbbs, zbbbs]).reshape(-1)
        for line in geqdsk_chunks(bbbs):
            f.write(line + "\n")
        lim = np.column_stack([rlim, zlim]).reshape(-1)
        for line in geqdsk_chunks(lim):
            f.write(line + "\n")


def load_eq_parameters(eq_prefix):
    candidates = [
        f"{eq_prefix}_parameters.csv",
        f"{eq_prefix}_eqgs2d_08_parameters_basic.csv",
        f"{eq_prefix}_eqgs2d_08_parameters.csv",
    ]
    path = None
    for cand in candidates:
        if os.path.exists(cand):
            path = cand
            break
    if path is None:
        raise FileNotFoundError(f"No EQ parameter CSV found for prefix {eq_prefix}")
    params = {}
    df = pd.read_csv(path)
    for _, row in df.iterrows():
        params[str(row.iloc[0]).strip()] = float(row.iloc[1])
    params["_path"] = path
    return params


def load_eq_grid_vector(path, value_name):
    df = pd.read_csv(path)
    return df[value_name].to_numpy(dtype=float)


def load_eq_psirz_matrix(path, r_target=None, z_target=None):
    df = pd.read_csv(path)
    arr = df.to_numpy(dtype=float)
    if r_target is not None and z_target is not None:
        nr = len(r_target)
        nz = len(z_target)
        if arr.shape == (nr, nz):
            return arr
        if arr.shape == (nz, nr):
            return arr.T
        raise ValueError(
            f"Unexpected EQ PSIRZ shape {arr.shape}, expected ({nr},{nz}) or ({nz},{nr})"
        )
    return arr


def pick_existing_path(candidates, what):
    for cand in candidates:
        if os.path.exists(cand):
            return cand
    raise FileNotFoundError(f"No {what} file found in candidates: {candidates}")


def try_load_eq_separatrix(eq_prefix):
    candidates = [
        f"{eq_prefix}_eqgs2d_07_separatrix.csv",
    ]
    for path in candidates:
        if os.path.exists(path):
            df = pd.read_csv(path)
            r = df.iloc[:, 1].to_numpy(dtype=float)
            z = df.iloc[:, 2].to_numpy(dtype=float)
            return path, r, z
    raise FileNotFoundError(
        f"Missing LCFS file for --from-eq prefix {eq_prefix}: expected one of {candidates}"
    )


def configure_from_eq_exports(g, eq_prefix):
    prof = load_eq_modeg5_profiles(eq_prefix)
    params = load_eq_parameters(eq_prefix)
    rg_path = pick_existing_path(
        [f"{eq_prefix}_RG_grid.csv"],
        "EQ RG grid",
    )
    zg_path = pick_existing_path(
        [f"{eq_prefix}_ZG_grid.csv"],
        "EQ ZG grid",
    )
    rg = load_eq_grid_vector(rg_path, "RG")
    zg = load_eq_grid_vector(zg_path, "ZG")
    psirz_path = pick_existing_path(
        [f"{eq_prefix}_PSIRZ.csv"],
        "EQ PSIRZ",
    )
    psirz_eq = load_eq_psirz_matrix(psirz_path, rg, zg)

    # In --from-eq mode, all geometry/state must come from the EQ export prefix.
    # Convert EQ full-Wb convention back to the raw gfile-like psi convention used
    # internally by this Python solver via the exported PSI0/PSIPA only.
    psi0_eq = float(params["PSI0"])
    psipa_eq = float(params["PSIPA"])
    psibdy_raw = 0.0
    psimag_raw = psi0_eq / (2.0 * np.pi)
    psirz = psirz_eq / (2.0 * np.pi) + psibdy_raw

    # Keep Python internal convention aligned with raw gfile arrays.
    psipa = float(params["PSIPA"])
    psi_n = prof["psips"] / psipa
    x_src = np.asarray(psi_n, dtype=float)
    pprime_src = np.asarray(prof["dppps"], dtype=float)
    ffprime_src = np.asarray(prof["ttdttps"] / (4.0 * np.pi * np.pi), dtype=float)
    x_src, pprime_src, ffprime_src = sanitize_axis_and_profiles(
        x_src, pprime_src, ffprime_src, "psi_n"
    )

    g.R = np.asarray(rg, dtype=float)
    g.Z = np.asarray(zg, dtype=float)
    g.nw = len(g.R)
    g.nh = len(g.Z)
    g.psirz = np.asarray(psirz, dtype=float)
    g.psimag = psimag_raw
    g.psibdy = psibdy_raw
    g.rmaxis = float(params["RAXIS"])
    g.zmaxis = float(params["ZAXIS"])
    g.psi_n = np.linspace(0.0, 1.0, g.nw)
    g.pres = np.interp(g.psi_n, x_src, prof["ppps"])
    g.fpol = np.interp(g.psi_n, x_src, -prof["ttps"] / (2.0 * np.pi))
    g.pprime = np.interp(g.psi_n, x_src, pprime_src)
    g.ffprime = np.interp(g.psi_n, x_src, ffprime_src)
    g.qpsi = np.interp(g.psi_n, x_src, prof["qqps"])

    sep_path, rbbbs, zbbbs = try_load_eq_separatrix(eq_prefix)
    g.rbbbs = rbbbs
    g.zbbbs = zbbbs
    g.nbbbs = len(rbbbs)

    src_desc = (
        f"{prof['path']} [converted from EQ MODELG=5 conventions; "
        f"psi restored to raw gfile convention]"
    )
    return g, x_src, pprime_src, ffprime_src, src_desc, "psi_n", prof, params, sep_path


def configure_from_eqdata(g, eqdata_path):
    data = read_eqdata_binary(eqdata_path)

    psi0_eq = float(data['PSI0'])
    psipa_eq = float(data['PSIPA'])
    psibdy_raw = 0.0
    psimag_raw = psi0_eq / (2.0 * np.pi)
    psirz = np.asarray(data['PSIRZ'], dtype=float) / (2.0 * np.pi) + psibdy_raw

    psips = np.asarray(data['PSIPS'], dtype=float)
    ppps = np.asarray(data['PPPS'], dtype=float)

    ttps = np.asarray(data['TTPS'], dtype=float)
    x_src = np.asarray(psips / psipa_eq, dtype=float)
    pprime_src = np.asarray(data['DPPPS'], dtype=float)
    ffprime_src = np.asarray(data['TTDTTPS'] / (4.0 * np.pi * np.pi), dtype=float)
    # print(psips,ttps,pprime_src,ffprime_src)
    x_src, pprime_src, ffprime_src = sanitize_axis_and_profiles(
        x_src, pprime_src, ffprime_src, 'psi_n'
    )

    g.R = np.asarray(data['RG'], dtype=float)
    g.Z = np.asarray(data['ZG'], dtype=float)
    g.nw = len(g.R)
    g.nh = len(g.Z)
    g.psirz = psirz
    g.psimag = psimag_raw
    g.psibdy = psibdy_raw
    g.rmaxis = float(data['RAXIS'])
    g.zmaxis = float(data['ZAXIS'])
    g.psi_n = np.linspace(0.0, 1.0, g.nw)
    g.pres = np.interp(g.psi_n, x_src, ppps)
    g.fpol = np.interp(g.psi_n, x_src, -ttps / (2.0 * np.pi))
    g.pprime = np.interp(g.psi_n, x_src, pprime_src)
    g.ffprime = np.interp(g.psi_n, x_src, ffprime_src)
    g.qpsi = np.interp(g.psi_n, x_src, np.asarray(data['QQPS'], dtype=float))
    g.currentA = float(data['RIP'])
    g.rcentr = float(data['RR'])
    g.bcentr = float(data['BB'])

    if int(data.get('NSUMAX_SAVED', 0)) <= 0 or len(data['RSU']) == 0:
        raise FileNotFoundError(f"eqdata {eqdata_path} does not contain saved separatrix arrays")
    g.rbbbs = np.asarray(data['RSU'], dtype=float)
    g.zbbbs = np.asarray(data['ZSU'], dtype=float)
    g.nbbbs = len(g.rbbbs)
    g.limitr = 1
    g.rlim = np.asarray([g.R[0]], dtype=float)
    g.zlim = np.asarray([g.Z[0]], dtype=float)

    src_desc = f"{eqdata_path} [direct eqdata binary]"
    return g, x_src, pprime_src, ffprime_src, src_desc, 'psi_n', data, data


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

    if 0 < i_ax < len(r) - 1 and 0 < j_ax < len(z) - 1:
        pts = []
        vals = []
        ok = True
        for di in (-1, 0, 1):
            for dj in (-1, 0, 1):
                ii = i_ax + di
                jj = j_ax + dj
                if mask is not None and not bool(mask[ii, jj]):
                    ok = False
                    break
                pts.append((float(r[ii] - r[i_ax]), float(z[jj] - z[j_ax])))
                vals.append(float(arr[ii, jj]))
            if not ok:
                break
        if ok:
            a = np.asarray(
                [[1.0, dR, dZ, dR * dR, dR * dZ, dZ * dZ] for dR, dZ in pts],
                dtype=float,
            )
            b = np.asarray(vals, dtype=float)
            coeff, *_ = np.linalg.lstsq(a, b, rcond=None)
            c0, c1, c2, c3, c4, c5 = coeff
            h = np.array([[2.0 * c3, c4], [c4, 2.0 * c5]], dtype=float)
            g = np.array([c1, c2], dtype=float)
            try:
                shift = -np.linalg.solve(h, g)
                dR, dZ = float(shift[0]), float(shift[1])
                dr_lim = abs(float(r[1] - r[0]))
                dz_lim = abs(float(z[1] - z[0]))
                if abs(dR) <= dr_lim and abs(dZ) <= dz_lim:
                    psi_ref = c0 + c1 * dR + c2 * dZ + c3 * dR * dR + c4 * dR * dZ + c5 * dZ * dZ
                    return float(r[i_ax] + dR), float(z[j_ax] + dZ), float(psi_ref)
            except np.linalg.LinAlgError:
                pass

    return float(r[i_ax]), float(z[j_ax]), float(arr[i_ax, j_ax])


def solve_gs_fixed_boundary(
    g,
    x_src,
    pprime_src,
    ffprime_src,
    source_x_kind="psi_n",
    rho_of_psin=None,
    max_outer=60,
    max_inner=400,
    omega=0.6,
    tol_inner=1e-6,
    freeze_source=False,
    mask_mode="both",
    source_scale=1.0,
    boundary_mode="rectangle",
    ghost_lambda_min=0.2,
    init_psi=None,
    recompute_axis=True,
    save_interval=0,
    out_prefix=None,
    g_for_plot=None,
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
            dprint("warning: LCFS polygon unavailable, disable source lcfs mask.")
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
        psi_bc_lcfs = g.psirz.copy()
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
        dprint(
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
        dprint(
            "axis mode: dynamic (from current psi), "
            f"initial axis=({axis_r:.5f},{axis_z:.5f}), psi_axis={axis_psi:.6e}"
        )

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

    for k in range(max_outer):
        psi_prev_outer = psi.copy()
        if recompute_axis:
            axis_r, axis_z, axis_psi = locate_axis_from_psi(
                psi, R, Z, mask=compare_mask, dpsi_sign=dpsi_sign
            )
        if rhs0 is None:
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
        outer_delta_now = float(np.max(np.abs(psi[compare_mask] - psi_prev_outer[compare_mask])))
        hist.append((k + 1, max_change, outer_delta_now, rmse, max_abs))
        dprint(
            f"[outer {k+1:02d}] inner_max_change={max_change:.3e}, "
            f"outer_delta={outer_delta_now:.3e}, "
            f"rmse_vs_gfile={rmse:.3e}, max_abs_vs_gfile={max_abs:.3e}, "
            f"inside_frac={inside_frac:.3f}, "
            f"lcfs_bc_max={lcfs_bc_max:.3e}, lcfs_bc_rms={lcfs_bc_rms:.3e}, "
            f"axis=({axis_r:.4f},{axis_z:.4f}), psi_axis={axis_psi:.6e}"
        )

        if outer_delta_now < 1.0e-5:
            break

        if save_interval > 0 and (k + 1) % save_interval == 0 and out_prefix and g_for_plot:
            save_outputs(
                f"{out_prefix}_step{k + 1:03d}",
                g_for_plot,
                psi,
                np.asarray(hist, dtype=float),
                "intermediate",
                compare_mask=compare_mask,
                boundary_mode=boundary_mode,
                is_intermediate=True,
            )

        # Stop when nonlinear outer loop is stable enough.
        # if max_change < max(tol_inner * 5.0, 1e-8):
        #     break

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
    is_intermediate=False,
    p_rebuilt=None,
    f_rebuilt=None,
    p_gfile=None,
    f_gfile=None,
    q_rebuilt=None,
    q_gfile=None,
    psi_n_eval=None,
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
    df_h = pd.DataFrame(
        history,
        columns=["outer_iter", "inner_max_change", "outer_delta", "rmse_vs_gfile", "max_abs_vs_gfile"],
    )
    hist_csv = f"{out_prefix}_solver_history.csv"
    df_h.to_csv(hist_csv, index=False)

    # Plot summary.
    err = psi_sol - g.psirz
    if compare_mask is None:
        compare_mask = np.ones_like(err, dtype=bool)
    err_plot = np.where(compare_mask, err, np.nan)
    r2d, z2d = np.meshgrid(g.R, g.Z, indexing="ij")

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

    if is_intermediate:
        fig, axes_2d = plt.subplots(2, 3, figsize=(18, 10), constrained_layout=True)
        ax_ov = axes_2d[0, 0]
        ax_sol = axes_2d[0, 1]
        ax_err = axes_2d[0, 2]
        ax_p = axes_2d[1, 0]
        ax_f = axes_2d[1, 1]
        ax_txt_int = axes_2d[1, 2]
        
        if psi_n_eval is not None:
            if p_rebuilt is not None:
                ax_p.plot(psi_n_eval, p_rebuilt/1e6, label="Target Input p", color="tab:blue", lw=2)
            if p_gfile is not None:
                ax_p.plot(psi_n_eval, p_gfile/1e6, label="Gfile p", color="tab:orange", ls="--", lw=2)
            ax_p.set_title("Pressure Profile")
            ax_p.set_xlabel("psi_n")
            ax_p.set_ylabel("p (MPa)")
            ax_p.grid(alpha=0.3)
            ax_p.legend()
            
            if f_rebuilt is not None:
                ax_f.plot(psi_n_eval, f_rebuilt, label="Rebuilt F", color="tab:green", lw=2)
            if f_gfile is not None:
                ax_f.plot(psi_n_eval, f_gfile, label="Gfile F", color="tab:red", ls="--", lw=2)
            ax_f.set_title("F Profile")
            ax_f.set_xlabel("psi_n")
            ax_f.set_ylabel("F (T*m)")
            ax_f.grid(alpha=0.3)
            ax_f.legend()
            
            if q_rebuilt is not None:
                ax_txt_int.plot(psi_n_eval, q_rebuilt, label="Rebuilt q", color="tab:purple", lw=2)
            if q_gfile is not None:
                ax_txt_int.plot(psi_n_eval, q_gfile, label="Gfile q", color="black", ls="--", lw=2)
            ax_txt_int.set_title("Safety Factor q")
            ax_txt_int.set_xlabel("psi_n")
            ax_txt_int.set_ylabel("q")
            ax_txt_int.grid(alpha=0.3)
            ax_txt_int.legend()
        else:
            ax_p.axis("off")
            ax_f.axis("off")
            ax_txt_int.axis("off")
    else:
        fig, axes_2d = plt.subplots(2, 4, figsize=(24, 10), constrained_layout=True)
        ax_gfile = axes_2d[0, 0]
        ax_ov = axes_2d[0, 1]
        ax_err = axes_2d[0, 2]
        ax_txt = axes_2d[0, 3]
        ax_sol = axes_2d[1, 0]
        ax_conv = axes_2d[1, 1]
        ax_j = axes_2d[1, 2]
        ax_jg = axes_2d[1, 3]

        cs0 = ax_gfile.contourf(r2d, z2d, g.psirz, levels=40, cmap="jet")
        fig.colorbar(cs0, ax=ax_gfile)
        ax_gfile.plot(g.rbbbs, g.zbbbs, "w--", lw=1.4, label="LCFS (gfile)")
        ax_gfile.set_title("Gfile psi(R,Z)")
        ax_gfile.set_xlabel("R (m)")
        ax_gfile.set_ylabel("Z (m)")
        ax_gfile.set_aspect("equal")
        ax_gfile.legend(loc="upper right", fontsize=8)

        # Summary text panel (moved to top-right).
        ax_txt.axis("off")
        if history.size > 0:
            last = history[-1]
            txt = (
                f"Final Iteration: {int(last[0])}\n"
                f"Inner Max Change: {last[1]:.3e}\n"
                f"Outer Delta: {last[2]:.3e}\n"
                f"RMSE vs gfile: {last[3]:.3e}\n"
                f"MAX|err| vs gfile: {last[4]:.3e}\n"
                f"psi_axis: {g.psimag:.6e}\n"
                f"psi_bdy:  {g.psibdy:.6e}\n"
                f"axis(gfile): ({g.rmaxis:.4f}, {g.zmaxis:.4f})\n"
                f"axis(solved): ({r_axis_sol:.4f}, {z_axis_sol:.4f})"
            )
        else:
            txt = "No history data."
        ax_txt.text(0.02, 0.98, txt, va="top", ha="left", fontsize=10, family="monospace")
        ax_txt.set_title("Summary Stats")

        # Convergence panel at bottom-middle-left.
        ax_conv.plot(history[:, 0], history[:, 1], "o-", label="inner max change")
        ax_conv.plot(history[:, 0], history[:, 2], "s-", label="outer delta")
        ax_conv.plot(history[:, 0], history[:, 3], "d-", label="RMSE vs gfile")
        ax_conv.plot(history[:, 0], history[:, 4], "^-", label="MAX|err| vs gfile")
        ax_conv.set_yscale("log")
        ax_conv.set_xlabel("outer iteration")
        ax_conv.set_title("Convergence History")
        ax_conv.grid(alpha=0.3)
        ax_conv.legend()

        # 2D jphi map panel for visualizing in/outside source behavior.
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


    cs2 = ax_err.contourf(
        r2d,
        z2d,
        err_plot,
        levels=40,
        cmap="coolwarm",
        vmin=-vmax,
        vmax=vmax,
    )
    fig.colorbar(cs2, ax=ax_err)
    ax_err.plot(g.rbbbs, g.zbbbs, "k--", lw=1.2, alpha=0.9, label="LCFS (gfile)")
    ax_err.set_title("psi_solved - psi_gfile")
    ax_err.set_xlabel("R (m)")
    ax_err.set_ylabel("Z (m)")
    ax_err.set_aspect("equal")
    ax_err.legend(loc="upper right", fontsize=8)

    # Solved psi panel. For intermediate, mask with compare_mask.
    if is_intermediate:
        psi_sol_plot = np.where(compare_mask, psi_sol, np.nan)
    else:
        psi_sol_plot = psi_sol
    cs1 = ax_sol.contourf(r2d, z2d, psi_sol_plot, levels=40, cmap="jet")
    fig.colorbar(cs1, ax=ax_sol)
    ax_sol.plot(g.rbbbs, g.zbbbs, "w--", lw=1.4, label="LCFS (gfile)")
    ax_sol.set_title("Solved psi(R,Z)")
    ax_sol.set_xlabel("R (m)")
    ax_sol.set_ylabel("Z (m)")
    ax_sol.set_aspect("equal")
    ax_sol.legend(loc="upper right", fontsize=8)

    # Overlay contour comparison is now embedded in summary.
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

    fig.suptitle(f"GS Solve From Gfile Profiles ({profile_source_desc})")
    fig_png = f"{out_prefix}_summary.png"
    fig.savefig(fig_png, dpi=160)
    plt.close(fig)

    return psi_csv, hist_csv, fig_png

def save_source_current_profiles(
    out_prefix, g, x_src, pprime_src, ffprime_src, source_x_kind="psi_n", r0=8.03, rho_of_psin=None,
    profile_csv_path=None,
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

    # If the profile CSV contains raw p or F inputs, overlay them so we can
    # distinguish "input profile" from "reconstructed from derivatives".
    if profile_csv_path and os.path.exists(profile_csv_path):
        try:
            df_in = pd.read_csv(profile_csv_path)
            if "psi_n" in df_in.columns:
                xin = df_in["psi_n"].to_numpy(dtype=float)
                if "p_pa" in df_in.columns:
                    pin = df_in["p_pa"].to_numpy(dtype=float)
                    axes[1, 0].plot(xin, pin / 1e6, color="tab:purple", lw=1.5, ls=":", label="Input p_pa")
                    axes[0, 2].plot(xin, pin / 1e6, color="tab:purple", lw=1.5, ls=":", label="Input p_pa")
                    axes[1, 0].legend()
                    axes[0, 2].legend()
                if "F_tesla_meter" in df_in.columns:
                    fin = df_in["F_tesla_meter"].to_numpy(dtype=float)
                    axes[1, 1].plot(xin, fin, color="tab:brown", lw=1.5, ls=":", label="Input F")
                    axes[1, 2].plot(xin, fin, color="tab:brown", lw=1.5, ls=":", label="Input F")
                    axes[1, 1].legend()
                    axes[1, 2].legend()
        except Exception:
            pass

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
        "--from-eq",
        default=None,
        help=(
            "Optional EQ export prefix. If set, read profiles and initial PSIRZ from "
            "<prefix>_modeg5_profiles.csv, <prefix>_RG_grid.csv, <prefix>_ZG_grid.csv, "
            "and <prefix>_PSIRZ.csv, with automatic convention conversion."
        ),
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
    parser.add_argument("--max-outer", type=int, default=60)
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
    parser.add_argument("--debug", action="store_true", help="Print debug/progress output")
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
        default="lcfs",
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
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Recompute psi-axis from current psi each outer iteration (default: enabled).",
    )
    parser.add_argument(
        "--save-interval",
        type=int,
        default=0,
        help="Save intermediate results every N outer iterations (0 to disable).",
    )
    parser.add_argument(
        "--downsample",
        type=int,
        nargs=2,
        metavar=("NW", "NH"),
        default=None,
        help="Downsample the gfile grid to (NW, NH) to speed up calculation.",
    )
    parser.add_argument(
        "--write-gfile",
        default=None,
        help="Optional output GEQDSK file path written from the current Python state.",
    )
    parser.add_argument(
        "--from-eqdata",
        default=None,
        help="Optional eqdata binary file. If set, read geometry/state directly from eqdata.",
    )
    args = parser.parse_args()

    global DEBUG
    DEBUG = bool(args.debug)

    if not (0.0 < args.omega <= 1.0):
        raise ValueError("For weighted Jacobi, omega must satisfy 0 < omega <= 1.")
    if not (0.0 < args.ghost_lambda_min <= 1.0):
        raise ValueError("--ghost-lambda-min must satisfy 0 < value <= 1.")

    gfile_path = os.path.abspath(args.gfile)
    if not os.path.exists(gfile_path):
        raise FileNotFoundError(f"gfile not found: {gfile_path}")

    if args.from_eq and args.from_eqdata:
        raise ValueError("Use only one of --from-eq or --from-eqdata.")

    g = read_gfile(gfile_path)
    if args.from_eqdata:
        eqdata_path = os.path.abspath(args.from_eqdata)
        g, x_src, pprime_src, ffprime_src, src_desc, source_x_kind, _eq_data, _eq_data2 = \
            configure_from_eqdata(g, eqdata_path)
        rho_of_psin = None
        profile_csv = eqdata_path
    elif args.from_eq:
        eq_prefix = os.path.abspath(args.from_eq)
        g, x_src, pprime_src, ffprime_src, src_desc, source_x_kind, _eq_prof, _eq_parm, sep_path = \
            configure_from_eq_exports(g, eq_prefix)
        rho_of_psin = None
        profile_csv = _eq_prof["path"]
        if sep_path is not None:
            src_desc += f", LCFS={sep_path}"
    else:
        profile_csv = resolve_profile_csv(gfile_path, args.profile_csv)
        x_src, pprime_src, ffprime_src, src_desc, source_x_kind = load_source_profiles(
            g, profile_csv, profile_x=args.profile_x
        )
        rho_of_psin = rho_from_qpsi(g.qpsi, g.psimag, g.psibdy) if source_x_kind == "rho_tor" else None

    if args.downsample:
        nw_new, nh_new = args.downsample
        dprint(f"Downsampling grid from {g.nw}x{g.nh} to {nw_new}x{nh_new}...")
        g = downsample_gfile(g, nw_new, nh_new)
        if args.from_eq or args.from_eqdata:
            raise ValueError("--downsample is not supported together with --from-eq/--from-eqdata yet.")

    init_psi = None
    init_desc = "gfile psirz"
    if args.init_psi_csv:
        init_path = os.path.abspath(args.init_psi_csv)
        if not os.path.exists(init_path):
            raise FileNotFoundError(f"initial psi CSV not found: {init_path}")
        init_psi, init_mode = load_initial_psi_from_csv(init_path, g.R, g.Z)
        init_desc = f"{init_path} ({init_mode})"
    elif args.from_eq:
        init_desc = f"{os.path.abspath(args.from_eq)}_PSIRZ.csv (from EQ export)"
    elif args.from_eqdata:
        init_desc = f"{os.path.abspath(args.from_eqdata)} [from eqdata binary]"

    dprint(f"gfile: {gfile_path}")
    if args.from_eq:
        dprint(f"from_eq: {os.path.abspath(args.from_eq)}")
    if args.from_eqdata:
        dprint(f"from_eqdata: {os.path.abspath(args.from_eqdata)}")
    dprint(f"profile source: {src_desc}")
    dprint(f"profile axis convention: {source_x_kind}")
    dprint(f"initial psi: {init_desc}")
    mask_mode = "none" if args.no_lcfs_mask else args.mask_mode
    dprint(f"mask_mode: {mask_mode}")
    dprint(f"boundary_mode: {args.boundary_mode}")
    dprint(f"ghost_lambda_min: {args.ghost_lambda_min}")
    dprint(f"source_scale: {args.source_scale}")
    dprint(f"r0_profile: {args.r0_profile}")
    dprint(f"recompute_axis: {args.recompute_axis}")
    dprint(
        f"grid: nw={g.nw}, nh={g.nh}, psi_axis={g.psimag:.6e}, "
        f"psi_bdy={g.psibdy:.6e}, Ip={g.currentA:.6e} A"
    )
    ggeom = derive_gfile_geometry_from_grid(g.R, g.Z, rcentr=getattr(g, 'rcentr', None))
    dprint(
        "gfile export geometry prep: "
        f"rdim={ggeom['rdim']:.8f}, zdim={ggeom['zdim']:.8f}, "
        f"rleft={ggeom['rleft']:.8f}, zmid={ggeom['zmid']:.8f}"
    )
    if 'rcentr' in ggeom:
        dprint(f"gfile export geometry prep: rcentr={ggeom['rcentr']:.8f}")
    gaxis = derive_gfile_axis_flux_fields(
        g.psirz,
        g.R,
        g.Z,
        g.psibdy,
        bcentr=getattr(g, 'bcentr', None),
        mask=np.ones_like(g.psirz, dtype=bool),
        dpsi_sign=1.0 if float(g.psibdy - g.psimag) >= 0.0 else -1.0,
    )
    dprint(
        "gfile export axis/flux prep: "
        f"rmaxis={gaxis['rmaxis']:.8f}, zmaxis={gaxis['zmaxis']:.8f}, "
        f"psimag={gaxis['psimag']:.8e}, psibdy={gaxis['psibdy']:.8e}"
    )
    if 'bcentr' in gaxis:
        dprint(f"gfile export axis/flux prep: bcentr={gaxis['bcentr']:.8f}")

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
        save_interval=args.save_interval,
        out_prefix=args.out_prefix,
        g_for_plot=g,
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
    )
    if args.write_gfile:
        axis_info_final = derive_gfile_axis_flux_fields(
            psi_sol,
            g.R,
            g.Z,
            g.psibdy,
            bcentr=getattr(g, 'bcentr', None),
            mask=compare_mask,
            dpsi_sign=1.0 if float(g.psibdy - g.psimag) >= 0.0 else -1.0,
        )
        write_gfile(args.write_gfile, g, psi_sol, axis_info_final)
        dprint(f"saved: {os.path.abspath(args.write_gfile)}")
    profile_fig_png = save_source_current_profiles(
        args.out_prefix,
        g,
        x_src,
        pprime_src,
        ffprime_src,
        source_x_kind=source_x_kind,
        r0=args.r0_profile,
        rho_of_psin=rho_of_psin,
        profile_csv_path=profile_csv,
    )
    dprint(f"saved: {psi_csv}")
    dprint(f"saved: {hist_csv}")
    dprint(f"saved: {fig_png}")
    dprint(f"saved: {profile_fig_png}")


if __name__ == "__main__":
    main()
