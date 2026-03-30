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


def resolve_eq_q_profile(current_dir, eq_dir, csv_arg):
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
        ax.plot(g["lcfs_r"], g["lcfs_z"], "b-", lw=1.8, label="gfile LCFS")
        if lcfs_eq_r is not None:
            ax.plot(lcfs_eq_r, lcfs_eq_z, "r--", lw=1.8, label="eqdata LCFS")
        g_axis = np.asarray(g["axis"])
        e_axis = np.asarray(eq["axis"])
        ax.plot([g_axis[0]], [g_axis[1]], "b*", ms=10, label="gfile axis")
        ax.plot([e_axis[0]], [e_axis[1]], "ro", ms=5, label="eqdata axis")
        ax.set_title("LCFS + Magnetic Axis")
        ax.set_xlabel("R (m)")
        ax.set_ylabel("Z (m)")
        ax.set_aspect("equal")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)

        txt_lines = []
        if axis_d is not None:
            txt_lines.append(f"dR={axis_d[0]:.3e} m")
            txt_lines.append(f"dZ={axis_d[1]:.3e} m")
            txt_lines.append(f"|dAxis|={axis_dist:.3e} m")
        if lcfs_m is not None:
            txt_lines.append(f"LCFS mean={lcfs_m['chamfer_mean']:.3e} m")
            txt_lines.append(f"LCFS max={lcfs_m['chamfer_max']:.3e} m")
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

    # Load eqgs2d PSI (supplementary)
    psi_gs2d = None
    r_gs2d = None
    z_gs2d = None
    try:
        from read_eqgs2d_csv import load_psirz_grid

        csv_path = find_first_existing(
            [
                os.path.join(eq_dir, "eqgs2d_01_PSIRZ_grid.csv"),
                os.path.join(current_dir, "eqgs2d_01_PSIRZ_grid.csv"),
            ]
        )
        if csv_path:
            r_gs2d, z_gs2d, psi_gs2d = load_psirz_grid(csv_path)
            print(f"Loaded eqgs2d CSV: {csv_path}")
    except Exception as e:
        print(f"Error loading eqgs2d CSV: {e}")

    # Load eqdata data
    prefix, paths = resolve_eqdata_prefix_paths(eq_dir, current_dir, args.prefix)
    if paths is None:
        print("eqdata CSV files not found")
        return

    eq_q_override = resolve_eq_q_profile(current_dir, eq_dir, args.eq_q_csv)
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
        lcfs_m = lcfs_distance_metrics(g["lcfs_r"], g["lcfs_z"], lcfs_eq_r, lcfs_eq_z)
        g_shape = boundary_shape_metrics(g["lcfs_r"], g["lcfs_z"])
        e_shape = boundary_shape_metrics(lcfs_eq_r, lcfs_eq_z)

        g_axis = np.asarray(g["axis"])
        e_axis = np.asarray(eq["axis"])
        axis_d = e_axis - g_axis
        axis_dist = float(np.hypot(axis_d[0], axis_d[1]))

        print("\n=== LCFS / Axis Metrics ===")
        print(
            f"Axis: dR={axis_d[0]:.4e} m, dZ={axis_d[1]:.4e} m, "
            f"distance={axis_dist:.4e} m"
        )
        print(
            f"LCFS Chamfer: mean={lcfs_m['chamfer_mean']:.4e} m, "
            f"max={lcfs_m['chamfer_max']:.4e} m"
        )
        if g_shape and e_shape:
            print(
                "LCFS shape delta: "
                f"da={e_shape['a']-g_shape['a']:.4e} m, "
                f"dkappa={e_shape['kappa']-g_shape['kappa']:.4e}, "
                f"ddelta={e_shape['delta']-g_shape['delta']:.4e}, "
                f"dArea={e_shape['area']-g_shape['area']:.4e} m^2"
            )

        fig_lcfs, ax_lcfs = plt.subplots(figsize=(7, 7))
        ax_lcfs.plot(g["lcfs_r"], g["lcfs_z"], "b-", lw=2, label="gfile LCFS")
        ax_lcfs.plot(lcfs_eq_r, lcfs_eq_z, "r--", lw=2, label="eqdata LCFS")
        ax_lcfs.plot([g_axis[0]], [g_axis[1]], "b*", ms=12, label="gfile axis")
        ax_lcfs.plot([e_axis[0]], [e_axis[1]], "ro", ms=6, label="eqdata axis")
        ax_lcfs.set_xlabel("R (m)")
        ax_lcfs.set_ylabel("Z (m)")
        ax_lcfs.set_aspect("equal")
        ax_lcfs.set_title("LCFS and Magnetic Axis Comparison")
        ax_lcfs.grid(alpha=0.25)
        ax_lcfs.legend()
        lcfs_plot = "lcfs_axis_comparison.png"
        fig_lcfs.savefig(lcfs_plot, dpi=150)
        print(f"Saved: {lcfs_plot}")

    elif lcfs_eq_r is None:
        print("Could not extract eqdata LCFS contour at PSI=0.")

    g_geo = g if (g is not None and "lcfs_r" in g and "axis" in g) else None
    if g_profile is not None:
        make_overview_figure(eq, g_geo, metrics, lcfs_eq_r, lcfs_eq_z, lcfs_m, axis_d, axis_dist)

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
    if g is not None:
        fig_ov, ax_ov = plt.subplots(figsize=(8, 10))
        r_g2d, z_g2d = np.meshgrid(g["R"], g["Z"])
        r_e2d, z_e2d = np.meshgrid(eq["R"], eq["Z"])
        ax_ov.contour(r_g2d, z_g2d, g["PSIRZ"], levels=levels, colors="blue", linewidths=1.2)
        ax_ov.contour(
            r_e2d,
            z_e2d,
            eq["PSIRZ"],
            levels=levels,
            colors="red",
            linestyles="dashed",
            linewidths=1.2,
        )
        ax_ov.set_title("Overlay: G-file (Blue/Solid) vs Eqdata (Red/Dashed)")
        ax_ov.set_xlabel("R (m)")
        ax_ov.set_ylabel("Z (m)")
        ax_ov.set_aspect("equal")
        ax_ov.grid(alpha=0.25)
        from matplotlib.lines import Line2D

        ax_ov.legend(
            handles=[
                Line2D([0], [0], color="blue", lw=1.5, label="G-file"),
                Line2D([0], [0], color="red", lw=1.5, linestyle="--", label="Eqdata"),
            ],
            loc="upper right",
        )
        overlay_plot = "gfile_eqdata_overlay.png"
        fig_ov.savefig(overlay_plot, dpi=150)
        print(f"Saved: {overlay_plot}")


if __name__ == "__main__":
    main()
