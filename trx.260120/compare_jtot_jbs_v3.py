#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Compare JTOT and JBS from TR CSV output with recalculation (Version 3).

This version uses exact geometry quantities from TR code:
- TTRHOG: <|grad_rho|^2/R^2> on grid
- ABVRHOG: <|grad_V|^2/R^2> on grid
- RDP: dPsi/drho
- RDPVRHOG: dPsi/dV on grid
- EPSRHO: epsilon = r_minor / R_major (inverse aspect ratio)
- RJCB: Jacobian

The Sauter model formula (from TRAJBSSAUTER in trcalc.f90):
  JBS = -PBSCD * TTRHOG * Pe * RKEV * 1e20
        * (L31*(dPe/Pe + dPi/Pe) + L32*dTe/Te + L34*alpha*(1-RPe)/RPe*dTi/Ti)
        / RDP / BB

Version 3 updates:
- choose a self-consistent CSV snapshot bundle (nearest-file matching)
- use TR-like 3-point derivative at the outer boundary
- print explicit warning that missing RW/PADD terms can keep JBS imperfect
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
import os
import re
import netCDF4 as nc

# ============================================================
# Physical constants (matching TR code trcomm.f90)
# ============================================================
RKEV = 1.602176634e-16   # J/keV (= 1.6e-19 * 1e3)
AME = 9.1093837e-31      # electron mass [kg]
AMP = 1.6726219e-27      # proton mass [kg]
RMU0 = 4.0 * np.pi * 1e-7  # vacuum permeability
PI = np.pi
EE = 1.602176634e-19     # electron charge [C]

# ============================================================
# Machine parameters (from tr.CFEDR0114.in)
# ============================================================
RR = 8.03       # major radius [m]
RA = 2.72       # minor radius [m]
BB = 6.0        # toroidal field [T]
RKAP = 1.89     # elongation
PBSCD = 1.0     # bootstrap current drive factor

# OMFIT current profile inputs
OMFIT_BASE_DIR = "/Users/dengxiaoya/CFEDRSW/OMFIT_out"
OMFIT_NAMELIST_PATH = os.path.join(OMFIT_BASE_DIR, "profiles_CFEDR.namelist")
OMFIT_STATEFILE_PATH = os.path.join(OMFIT_BASE_DIR, "statefile_3.000000E+01.nc")

# ============================================================
# Helper functions
# ============================================================

def find_csv_by_title(pattern, data_dir='.'):
    """Find first CSV file whose first line contains the given pattern."""
    files = sorted(glob.glob(os.path.join(data_dir, 'tr_data_*.csv')))
    for f in files:
        with open(f, 'r') as fh:
            title = fh.readline().strip()
            if pattern in title:
                return f
    return None


def find_csv_candidates_by_title(pattern, data_dir='.'):
    """Find all CSV files whose first line contains the given pattern."""
    matches = []
    files = sorted(glob.glob(os.path.join(data_dir, "tr_data_*.csv")))
    for f in files:
        with open(f, "r") as fh:
            title = fh.readline().strip()
        if pattern in title:
            matches.append(f)
    return matches


def csv_index(path):
    """Extract integer index from tr_data_XXX.csv."""
    name = os.path.basename(path)
    m = re.match(r"tr_data_(\d+)\.csv$", name)
    if m is None:
        return None
    return int(m.group(1))


def nearest_csv_distance(path, anchor_idx):
    idx = csv_index(path)
    if idx is None:
        return 10**9
    return abs(idx - anchor_idx)


def choose_snapshot_files(data_dir='.'):
    """
    Pick a consistent snapshot bundle by anchoring on current-density CSV
    and selecting nearest-index files for all required patterns.
    """
    patterns = {
        "current": "JTOT,JOH,JNB,JRF,JBS",
        "density": "n(NS)",
        "temperature": "T(NS)",
        "q": "QP  vs r",
        "zeff": "ZEFF  vs r",
        "fast_energy": "WB,WF [MJ]  vs r",
        "geom1": "TTRHOG,ABVRHOG",
        "geom2": "RDP,RDPVRHOG",
        "geom3": "EPSRHO,RJCB",
        "geom4": "RMJRHO,RMNRHO",
    }
    required = ("current", "density", "temperature", "q", "zeff")

    candidates = {
        key: find_csv_candidates_by_title(pattern, data_dir)
        for key, pattern in patterns.items()
    }
    if not candidates["current"]:
        raise RuntimeError("Cannot find current density CSV by title pattern.")

    best_bundle = None
    best_score = None
    for current_file in candidates["current"]:
        idx0 = csv_index(current_file)
        if idx0 is None:
            continue
        bundle = {"current": current_file}
        missing_required = 0
        score = 0
        idx_list = [idx0]

        for key in patterns:
            if key == "current":
                continue
            cands = candidates[key]
            if not cands:
                bundle[key] = None
                if key in required:
                    missing_required += 1
                continue

            nearest = min(cands, key=lambda p: nearest_csv_distance(p, idx0))
            idxn = csv_index(nearest)
            bundle[key] = nearest
            if idxn is None:
                if key in required:
                    missing_required += 1
            else:
                idx_list.append(idxn)
                score += abs(idxn - idx0)

        spread = max(idx_list) - min(idx_list) if idx_list else 10**9
        tie_break = (missing_required, score, spread, idx0)
        if best_score is None or tie_break < best_score:
            best_score = tie_break
            best_bundle = bundle

    if best_bundle is None:
        raise RuntimeError("Failed to choose a valid CSV snapshot bundle.")
    return best_bundle

def read_tr_csv(filepath):
    """Read a TR CSV file. Returns (r, columns_dict, title)."""
    with open(filepath, 'r') as f:
        title_line = f.readline().strip()
        header_line = f.readline().strip()

    headers = [h.strip() for h in header_line.split(',')]
    data = pd.read_csv(filepath, skiprows=2, header=None, names=headers)
    r = data.iloc[:, 0].values
    cols = {}
    for i, name in enumerate(headers):
        cols[name] = data.iloc[:, i].values
    return r, cols, title_line


def deriv3p(y0, y1, y2, x0, x1, x2):
    """
    3-point nonuniform derivative at x0 from quadratic interpolation.
    Equivalent to TR-style DERIV3P usage at the edge.
    """
    den0 = (x0 - x1) * (x0 - x2)
    den1 = (x1 - x0) * (x1 - x2)
    den2 = (x2 - x0) * (x2 - x1)
    eps = 1e-30
    if abs(den0) < eps or abs(den1) < eps or abs(den2) < eps:
        return (y0 - y1) / (x0 - x1 + eps)
    term0 = y0 * (2.0 * x0 - x1 - x2) / den0
    term1 = y1 * (x0 - x2) / den1
    term2 = y2 * (x0 - x1) / den2
    return term0 + term1 + term2

def trapped_fraction_kim(eps):
    """Trapped particle fraction: Kim et al., PoF B 3 2050 (1991) eq(C18)."""
    return 1.46 * np.sqrt(eps) - 0.46 * eps**1.5


def parse_top_level_array(filepath, varname):
    """Parse a top-level namelist array variable."""
    if not os.path.exists(filepath):
        return None
    with open(filepath, "r", encoding="utf-8") as f:
        text = f.read()

    pattern = re.compile(rf"^\s*{re.escape(varname)}\s*=\s*", re.MULTILINE | re.IGNORECASE)
    match = pattern.search(text)
    if match is None:
        return None

    values = []
    rest = text[match.end():]
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
                    values.extend([float(v_str)] * int(c_str))
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


def first_namelist_match(path, candidates):
    for name in candidates:
        arr = parse_top_level_array(path, name)
        if arr is not None:
            return arr, name
    return None, None


def first_state_var(ds, candidates):
    for name in candidates:
        if name in ds.variables:
            return np.asarray(ds.variables[name][:], dtype=float), name
    return None, None


def interpolate_to_rho(rho_src, val_src, rho_dst):
    idx = np.argsort(rho_src)
    return np.interp(rho_dst, rho_src[idx], val_src[idx])


def load_omfit_current_profiles():
    """Load OMFIT current profiles from namelist/statefile."""
    nm_rho, _ = first_namelist_match(OMFIT_NAMELIST_PATH, ["rho"])

    nm_map = {
        "tot": ["current_density"],
        "oh": ["ohmic_current_density_onetwo", "ohmic_current_density"],
        "bs": ["bootstrap_current_density_onetwo", "bootstrap_current_density"],
        "ec": [
            "ec_current_density_onetwo",
            "eccd_current_density_onetwo",
            "ec_current_density",
            "eccd_current_density",
        ],
        "ic": [
            "ic_current_density_onetwo",
            "iccd_current_density_onetwo",
            "ic_current_density",
            "iccd_current_density",
        ],
        "rf": ["rfcd_current_density_onetwo", "rf_current_density_onetwo"],
    }
    sf_map = {
        "tot": ["curden"],
        "oh": ["curohm"],
        "bs": ["curboot"],
        "ec": ["curec", "cureccd", "curech", "eccd", "j_eccd"],
        "ic": ["curic", "curiccd", "curich", "iccd", "j_iccd"],
        "rf": ["currf"],
    }

    data = {}
    sources = {}
    for key, candidates in nm_map.items():
        arr, used = first_namelist_match(OMFIT_NAMELIST_PATH, candidates)
        if arr is not None:
            data[key] = arr
            sources[key] = f"namelist:{used}"

    sf_rho = None
    if os.path.exists(OMFIT_STATEFILE_PATH):
        try:
            with nc.Dataset(OMFIT_STATEFILE_PATH, "r") as ds:
                sf_rho_raw = np.asarray(ds.variables["rho_grid"][:], dtype=float)
                if np.nanmax(np.abs(sf_rho_raw)) > 1.0:
                    sf_rho = sf_rho_raw / np.nanmax(np.abs(sf_rho_raw))
                else:
                    sf_rho = sf_rho_raw

                for key, candidates in sf_map.items():
                    if key in data:
                        continue
                    arr, used = first_state_var(ds, candidates)
                    if arr is not None:
                        data[key] = arr
                        sources[key] = f"statefile:{used}"
        except Exception as exc:
            print(f"   WARNING: Failed to read OMFIT statefile: {exc}")

    if nm_rho is not None:
        rho = nm_rho
    elif sf_rho is not None:
        rho = sf_rho
    else:
        return None, None

    for key, arr in list(data.items()):
        if len(arr) == len(rho):
            continue
        if sources[key].startswith("statefile:") and sf_rho is not None:
            data[key] = interpolate_to_rho(sf_rho, arr, rho)
        elif sources[key].startswith("namelist:") and nm_rho is not None and len(nm_rho) == len(arr):
            data[key] = interpolate_to_rho(nm_rho, arr, rho)
        else:
            del data[key]
            del sources[key]

    if not data:
        return None, None

    for key in list(data.keys()):
        data[key] = data[key] / 1e6
    return rho, data

def F31_func(x, Z):
    """Sauter L31 coefficient."""
    return ((1.0 + 1.4/(Z+1.0))*x - 1.9/(Z+1.0)*x**2
            + 0.3/(Z+1.0)*x**3 + 0.2/(Z+1.0)*x**4)

def F32EE_func(x, Z):
    """Sauter L32 electron-electron part."""
    return ((0.05 + 0.62*Z)/(Z*(1.0 + 0.44*Z))*(x - x**4)
            + 1.0/(1.0 + 0.22*Z)*(x**2 - x**4 - 1.2*(x**3 - x**4))
            + 1.2/(1.0 + 0.5*Z)*x**4)

def F32EI_func(x, Z):
    """Sauter L32 electron-ion part."""
    return (-(0.56 + 1.93*Z)/(Z*(1.0 + 0.44*Z))*(x - x**4)
            + 4.95/(1.0 + 2.48*Z)*(x**2 - x**4 - 0.55*(x**3 - x**4))
            - 1.2/(1.0 + 0.5*Z)*x**4)


def main():
    data_dir = '.'

    print("=" * 70)
    print("JTOT and JBS COMPARISON V3: Snapshot-Consistent + TR Edge Derivative")
    print("=" * 70)

    # ============================================================
    # 1. Read CSV data
    # ============================================================
    print("\n1. Loading CSV data...")

    bundle = choose_snapshot_files(data_dir)
    print("   Selected snapshot bundle (nearest-index matching):")
    for key in ("current", "density", "temperature", "q", "zeff", "fast_energy", "geom1", "geom2", "geom3", "geom4"):
        fpath = bundle.get(key)
        if fpath is None:
            print(f"     {key:11s}: MISSING")
        else:
            print(f"     {key:11s}: {os.path.basename(fpath)}")

    # Current densities
    f_j = bundle["current"]
    if f_j is None:
        print("ERROR: Cannot find current density CSV")
        return
    r_j, j_cols, _ = read_tr_csv(f_j)
    print(f"   Current densities loaded: {os.path.basename(f_j)}")

    # Densities
    f_n = bundle["density"]
    if f_n is None:
        print("ERROR: Cannot find density CSV")
        return
    r_n, n_cols, _ = read_tr_csv(f_n)
    print(f"   Densities loaded: {os.path.basename(f_n)}")

    # Temperatures
    f_t = bundle["temperature"]
    if f_t is None:
        print("ERROR: Cannot find temperature CSV")
        return
    r_t, t_cols, _ = read_tr_csv(f_t)
    print(f"   Temperatures loaded: {os.path.basename(f_t)}")

    # Safety factor
    f_q = bundle["q"]
    if f_q is None:
        print("ERROR: Cannot find safety factor CSV")
        return
    r_q, q_cols, _ = read_tr_csv(f_q)
    print(f"   Safety factor loaded: {os.path.basename(f_q)}")

    # Zeff
    f_z = bundle["zeff"]
    if f_z is None:
        print("ERROR: Cannot find Zeff CSV")
        return
    r_z, z_cols, _ = read_tr_csv(f_z)
    print(f"   Zeff loaded: {os.path.basename(f_z)}")

    # Fast ion pressure proxy (WB/WF -> RW)
    f_fast = bundle["fast_energy"]
    has_fast = f_fast is not None
    if has_fast:
        r_f, fast_cols, _ = read_tr_csv(f_fast)
        print(f"   Fast energy loaded: {os.path.basename(f_fast)}")
    else:
        r_f, fast_cols = None, None
        print("   Fast energy CSV not found: RW contribution will be set to 0")

    # Geometry quantities (new in V2)
    f_geom1 = bundle["geom1"]
    f_geom2 = bundle["geom2"]
    f_geom3 = bundle["geom3"]
    f_geom4 = bundle["geom4"]

    has_geometry = all([f_geom1, f_geom2, f_geom3, f_geom4])

    if has_geometry:
        r_g1, g1_cols, _ = read_tr_csv(f_geom1)
        r_g2, g2_cols, _ = read_tr_csv(f_geom2)
        r_g3, g3_cols, _ = read_tr_csv(f_geom3)
        r_g4, g4_cols, _ = read_tr_csv(f_geom4)
        print(f"   Geometry (TTRHOG,ABVRHOG): {os.path.basename(f_geom1)}")
        print(f"   Geometry (RDP,RDPVRHOG): {os.path.basename(f_geom2)}")
        print(f"   Geometry (EPSRHO,RJCB): {os.path.basename(f_geom3)}")
        print(f"   Geometry (RMJRHO,RMNRHO): {os.path.basename(f_geom4)}")
    else:
        print("   WARNING: Geometry CSV files not found!")
        print("   Please run 'G R E' in TR menu to generate geometry output.")
        print("   Falling back to approximate geometry...")

    # ============================================================
    # 2. Extract data arrays
    # ============================================================
    print("\n2. Extracting data arrays...")

    j_keys = list(j_cols.keys())
    n_keys = list(n_cols.keys())
    t_keys = list(t_cols.keys())
    q_keys = list(q_cols.keys())
    z_keys = list(z_cols.keys())

    # Current densities [MA/m^2]
    r = j_cols[j_keys[0]]
    JTOT_csv = j_cols[j_keys[1]]
    JOH_csv  = j_cols[j_keys[2]]
    JNB_csv  = j_cols[j_keys[3]]
    JRF_csv  = j_cols[j_keys[4]]
    JBS_csv  = j_cols[j_keys[5]]

    NRMAX = len(r)
    DR = 1.0 / NRMAX
    RHOG = r.copy()
    RHOM = np.clip(r + 0.5 * DR, 0.0, 1.0)

    # Densities [10^20/m^3]
    ne = n_cols[n_keys[1]]
    nD = n_cols[n_keys[2]]
    nT = n_cols[n_keys[3]]
    nA = n_cols[n_keys[4]]

    # Temperatures [keV]
    Te = t_cols[t_keys[1]]
    TD = t_cols[t_keys[2]]
    TT = t_cols[t_keys[3]]
    TA = t_cols[t_keys[4]]

    # Safety factor (interpolate if needed)
    r_qp = q_cols[q_keys[0]]
    QP_val = q_cols[q_keys[1]]
    if len(QP_val) != NRMAX:
        QP = np.interp(r, r_qp, QP_val)
    else:
        QP = QP_val

    # Zeff
    Zeff = z_cols[z_keys[1]]

    # Fast ion energy density RW [same unit as RN*RT, i.e. 1e20 m^-3 * keV]
    if has_fast:
        f_keys = list(fast_cols.keys())
        rw_sum = np.zeros(NRMAX)
        # CSV columns after X are WB/WF for each fast species in MJ.
        for k in f_keys[1:]:
            vals = fast_cols[k]
            if len(vals) != NRMAX:
                vals = np.interp(r, fast_cols[f_keys[0]], vals)
            rw_sum += vals / (1.5 * RKEV * 1.0e14)
        print(f"   Using RW reconstructed from WB/WF (max={rw_sum.max():.3f})")
    else:
        rw_sum = np.zeros(NRMAX)
        print("   RW set to zero (no WB/WF radial CSV)")

    # Geometry quantities
    if has_geometry:
        g1_keys = list(g1_cols.keys())
        g2_keys = list(g2_cols.keys())
        g3_keys = list(g3_cols.keys())
        g4_keys = list(g4_cols.keys())

        TTRHOG = g1_cols[g1_keys[1]]
        ABVRHOG = g1_cols[g1_keys[2]]
        RDP = g2_cols[g2_keys[1]]
        RDPVRHOG = g2_cols[g2_keys[2]]
        EPSRHO = g3_cols[g3_keys[1]]
        RJCB = g3_cols[g3_keys[2]]
        RMJRHO = g4_cols[g4_keys[1]]
        RMNRHO = g4_cols[g4_keys[2]]

        eps = EPSRHO
        print(f"   Using exact EPSRHO from TR")
    else:
        # Fallback to approximate geometry
        eps = r * RA / RR
        TTRHOG = np.ones(NRMAX) / RR**2  # Approximate
        RDP = np.zeros(NRMAX)  # Will estimate from Bp
        print(f"   Using approximate geometry (eps = r*a/R)")

    eps = np.maximum(eps, 1e-10)

    print(f"   Grid: NRMAX={NRMAX}, DR={DR:.4f}")
    print(f"   ne center={ne[0]:.4f}, Te center={Te[0]:.4f} keV")
    print(f"   eps center={eps[0]:.4f}, eps edge={eps[-1]:.4f}")

    # ============================================================
    # 3. Verify JTOT = JOH + JNB + JRF + JBS
    # ============================================================
    print("\n3. Verifying JTOT = JOH + JNB + JRF + JBS...")

    JTOT_sum = JOH_csv + JNB_csv + JRF_csv + JBS_csv
    diff_jtot = JTOT_csv - JTOT_sum
    max_diff = np.max(np.abs(diff_jtot))
    rel_diff = np.max(np.abs(diff_jtot / (np.abs(JTOT_csv) + 1e-20)))

    print(f"   Max absolute difference: {max_diff:.6e} MA/m^2")
    print(f"   Max relative difference: {rel_diff:.6e}")
    if rel_diff < 1e-3:
        print("   ✓ JTOT = JOH + JNB + JRF + JBS is consistent")

    # ============================================================
    # 4. Recalculate JBS using Sauter model with exact geometry
    # ============================================================
    print("\n4. Recalculating JBS using Sauter model...")
    print("   NOTE: Pi includes thermal ions + reconstructed RW from WB/WF.")
    print("         PADD is still unavailable in CSV, so residual mismatch may remain.")

    # Pressures [10^20/m^3 * keV]
    Pe = ne * Te
    Pi_th = nD * TD + nT * TT + nA * TA
    Pi = Pi_th + rw_sum
    Ti = Pi / (nD + nT + nA + 1e-20)
    ANI_arr = nD + nT + 2*nA  # Charge-weighted ion density

    # JBS calculation using Sauter formula - ALL calculations inside loop to match TR exactly
    # TR code uses average values between adjacent grid points for ALL quantities
    JBS_calc = np.zeros(NRMAX)

    # Arrays to store coefficients for diagnostics
    ft = np.zeros(NRMAX)
    nue_star = np.zeros(NRMAX)
    nui_star = np.zeros(NRMAX)
    RL31 = np.zeros(NRMAX)
    RL32 = np.zeros(NRMAX)
    RL34 = np.zeros(NRMAX)
    salfa = np.zeros(NRMAX)

    for i in range(NRMAX - 1):
        if not has_geometry or abs(RDP[i]) < 1e-20:
            continue

        # EPS at grid point (not averaged) - matches TR: EPS=EPSRHO(NR)
        EPS_i = eps[i]
        EPSS = np.sqrt(EPS_i)**3

        # Trapped fraction at grid point - matches TR: FT=FTPF(MDLTPF,EPS)
        ft[i] = 1.46 * np.sqrt(EPS_i) - 0.46 * EPS_i**1.5

        # Zeff averaged - matches TR: ZEFFL=0.5D0*(ZEFF(NR)+ZEFF(NR+1))
        Zeff_avg = 0.5 * (Zeff[i+1] + Zeff[i])

        # Averaged quantities for densities and temperatures
        ne_avg = 0.5 * (ne[i+1] + ne[i])
        Te_avg = 0.5 * (Te[i+1] + Te[i])
        Ti_avg = 0.5 * (Ti[i+1] + Ti[i])
        ANI_avg = 0.5 * (ANI_arr[i+1] + ANI_arr[i])

        Pe_avg = 0.5 * (Pe[i+1] + Pe[i])
        Pi_avg = 0.5 * (Pi[i+1] + Pi[i])
        RPe_avg = Pe_avg / (Pe_avg + Pi_avg + 1e-20)

        # Gradients: (f[i+1] - f[i]) / DR (forward difference per normalized rho)
        DPe_i = (Pe[i+1] - Pe[i]) / DR
        DPi_i = (Pi[i+1] - Pi[i]) / DR
        DTe_i = (Te[i+1] - Te[i]) / DR
        DTi_i = (Ti[i+1] - Ti[i]) / DR

        # Safety factor at grid point
        QL = abs(QP[i])

        # Electron collisionality - using averaged ne and Te
        ln_Lambda_e = 31.3 - np.log(np.sqrt(ne_avg * 1e20) / abs(Te_avg * 1e3 + 1e-10))
        nue_star[i] = (6.921e-18 * QL * RR * ne_avg * 1e20 * Zeff_avg * ln_Lambda_e
                      / (abs(Te_avg * 1e3)**2 * EPSS + 1e-30))

        # Ion collisionality - using averaged ANI and Ti
        PZ_D = 1.0
        ln_Lambda_ii = 30.0 - np.log(PZ_D**3 * np.sqrt(ANI_avg * 1e20) / (abs(Ti_avg * 1e3 + 1e-10)**1.5))
        nui_star[i] = (4.90e-18 * QL * RR * ANI_avg * 1e20 * PZ_D**4 * ln_Lambda_ii
                      / (abs(Ti_avg * 1e3)**2 * EPSS + 1e-30))

        # Effective trapped fractions for Sauter coefficients
        F31_teff = ft[i] / (1.0 + (1.0 - 0.1*ft[i])*np.sqrt(nue_star[i]) + 0.5*(1.0 - ft[i])*nue_star[i]/Zeff_avg)
        F32ee_teff = ft[i] / (1.0 + 0.26*(1.0-ft[i])*np.sqrt(nue_star[i]) + 0.18*(1.0-0.37*ft[i])*nue_star[i]/np.sqrt(Zeff_avg))
        F32ei_teff = ft[i] / (1.0 + (1.0+0.6*ft[i])*np.sqrt(nue_star[i]) + 0.85*(1.0-0.37*ft[i])*nue_star[i]*(1.0+Zeff_avg))
        F34_teff = ft[i] / (1.0 + (1.0 - 0.1*ft[i])*np.sqrt(nue_star[i]) + 0.5*(1.0 - 0.5*ft[i])*nue_star[i]/Zeff_avg)

        # Sauter alpha coefficient
        salfa0 = -1.17*(1.0 - ft[i])/(1.0 - 0.22*ft[i] - 0.19*ft[i]**2 + 1e-20)
        salfa[i] = ((salfa0 + 0.25*(1.0 - ft[i]**2)*np.sqrt(nui_star[i])) / (1.0 + 0.5*np.sqrt(nui_star[i]))
                   + 0.315*nui_star[i]**2 * ft[i]**6) / (1.0 + 0.15*nui_star[i]**2 * ft[i]**6)

        # Neoclassical coefficients
        RL31[i] = F31_func(F31_teff, Zeff_avg)
        RL32[i] = F32EE_func(F32ee_teff, Zeff_avg) + F32EI_func(F32ei_teff, Zeff_avg)
        RL34[i] = F31_func(F34_teff, Zeff_avg)

        # Gradient terms
        term1 = RL31[i] * (DPe_i / (Pe_avg + 1e-20) + DPi_i / (Pe_avg + 1e-20))
        term2 = RL32[i] * DTe_i / (Te_avg + 1e-20)
        term3 = RL34[i] * salfa[i] * (1.0 - RPe_avg) / (RPe_avg + 1e-20) * DTi_i / (Ti_avg + 1e-20)

        # Pe in SI: Pe_avg * 1e20 * RKEV [Pa]
        Pe_SI = Pe_avg * 1e20 * RKEV

        # Sauter formula: JBS = -PBSCD * TTRHOG * Pe * (terms) / RDP / BB
        JBS_calc[i] = -PBSCD * TTRHOG[i] * Pe_SI * (term1 + term2 + term3) / RDP[i] / BB

    # Boundary point (i = NRMAX-1) with TR-like 3-point derivative
    i = NRMAX - 1
    if has_geometry and abs(RDP[i]) > 1e-20:
        EPS_i = eps[i]
        EPSS = np.sqrt(EPS_i)**3
        ft[i] = 1.46 * np.sqrt(EPS_i) - 0.46 * EPS_i**1.5

        # Extrapolate Zeff: ZEFFL=2.D0*ZEFF(NR-1)-ZEFF(NR-2)
        Zeff_avg = 2.0 * Zeff[i-1] - Zeff[i-2]

        # Approximate boundary values (CSV has no PNSS/PTS boundary arrays)
        ne_edge = 2.0 * ne[i] - ne[i-1]
        Te_edge = 2.0 * Te[i] - Te[i-1]
        Ti_edge = 2.0 * Ti[i] - Ti[i-1]
        Pe_edge = ne_edge * Te_edge
        Pi_edge = 2.0 * Pi[i] - Pi[i-1]
        ANI_edge = 2.0 * ANI_arr[i] - ANI_arr[i-1]

        Pe_avg = Pe_edge
        Pi_avg = Pi_edge
        Te_avg = Te_edge
        Ti_avg = Ti_edge
        ne_avg = ne_edge
        ANI_avg = ANI_edge
        RPe_avg = Pe_avg / (Pe_avg + Pi_avg + 1e-20)

        # TR-like nonuniform 3-point derivative at outer grid
        x0 = RHOG[i]
        x1 = RHOM[i]
        x2 = RHOM[i - 1]
        DPe_i = deriv3p(Pe_edge, Pe[i], Pe[i - 1], x0, x1, x2)
        DPi_i = deriv3p(Pi_edge, Pi[i], Pi[i - 1], x0, x1, x2)
        DTe_i = deriv3p(Te_edge, Te[i], Te[i - 1], x0, x1, x2)
        DTi_i = deriv3p(Ti_edge, Ti[i], Ti[i - 1], x0, x1, x2)

        QL = abs(QP[i])
        ln_Lambda_e = 31.3 - np.log(np.sqrt(ne_avg * 1e20) / abs(Te_avg * 1e3 + 1e-10))
        nue_star[i] = (6.921e-18 * QL * RR * ne_avg * 1e20 * Zeff_avg * ln_Lambda_e
                      / (abs(Te_avg * 1e3)**2 * EPSS + 1e-30))

        PZ_D = 1.0
        ln_Lambda_ii = 30.0 - np.log(PZ_D**3 * np.sqrt(ANI_avg * 1e20) / (abs(Ti_avg * 1e3 + 1e-10)**1.5))
        nui_star[i] = (4.90e-18 * QL * RR * ANI_avg * 1e20 * PZ_D**4 * ln_Lambda_ii
                      / (abs(Ti_avg * 1e3)**2 * EPSS + 1e-30))

        F31_teff = ft[i] / (1.0 + (1.0 - 0.1*ft[i])*np.sqrt(nue_star[i]) + 0.5*(1.0 - ft[i])*nue_star[i]/Zeff_avg)
        F32ee_teff = ft[i] / (1.0 + 0.26*(1.0-ft[i])*np.sqrt(nue_star[i]) + 0.18*(1.0-0.37*ft[i])*nue_star[i]/np.sqrt(Zeff_avg))
        F32ei_teff = ft[i] / (1.0 + (1.0+0.6*ft[i])*np.sqrt(nue_star[i]) + 0.85*(1.0-0.37*ft[i])*nue_star[i]*(1.0+Zeff_avg))
        F34_teff = ft[i] / (1.0 + (1.0 - 0.1*ft[i])*np.sqrt(nue_star[i]) + 0.5*(1.0 - 0.5*ft[i])*nue_star[i]/Zeff_avg)

        salfa0 = -1.17*(1.0 - ft[i])/(1.0 - 0.22*ft[i] - 0.19*ft[i]**2 + 1e-20)
        salfa[i] = ((salfa0 + 0.25*(1.0 - ft[i]**2)*np.sqrt(nui_star[i])) / (1.0 + 0.5*np.sqrt(nui_star[i]))
                   + 0.315*nui_star[i]**2 * ft[i]**6) / (1.0 + 0.15*nui_star[i]**2 * ft[i]**6)

        RL31[i] = F31_func(F31_teff, Zeff_avg)
        RL32[i] = F32EE_func(F32ee_teff, Zeff_avg) + F32EI_func(F32ei_teff, Zeff_avg)
        RL34[i] = F31_func(F34_teff, Zeff_avg)

        term1 = RL31[i] * (DPe_i / (Pe_avg + 1e-20) + DPi_i / (Pe_avg + 1e-20))
        term2 = RL32[i] * DTe_i / (Te_avg + 1e-20)
        term3 = RL34[i] * salfa[i] * (1.0 - RPe_avg) / (RPe_avg + 1e-20) * DTi_i / (Ti_avg + 1e-20)

        Pe_SI = Pe_avg * 1e20 * RKEV
        JBS_calc[i] = -PBSCD * TTRHOG[i] * Pe_SI * (term1 + term2 + term3) / RDP[i] / BB

    # Convert to MA/m^2
    JBS_calc_MA = JBS_calc * 1e-6

    # Average adjacent points (as in TR code)
    JBS_avg = np.zeros(NRMAX)
    JBS_avg[0] = 0.5 * JBS_calc_MA[0]
    for i in range(1, NRMAX):
        JBS_avg[i] = 0.5 * (JBS_calc_MA[i] + JBS_calc_MA[i-1])

    # ============================================================
    # 5. Compare results
    # ============================================================
    print("\n5. Comparison results:")
    print(f"\n   {'rho':>6} {'JBS_csv':>12} {'JBS_calc':>12} {'Ratio':>10} {'Diff%':>10}")
    print(f"   {'-'*54}")
    for idx in [0, NRMAX//10, NRMAX//4, NRMAX//2, 3*NRMAX//4, NRMAX-1]:
        ratio = JBS_avg[idx] / (JBS_csv[idx] + 1e-30) if abs(JBS_csv[idx]) > 1e-10 else float('nan')
        diff_pct = (JBS_avg[idx] - JBS_csv[idx]) / (abs(JBS_csv[idx]) + 1e-30) * 100 if abs(JBS_csv[idx]) > 1e-10 else float('nan')
        print(f"   {r[idx]:6.3f} {JBS_csv[idx]:12.4f} {JBS_avg[idx]:12.4f} {ratio:10.4f} {diff_pct:10.1f}%")

    # Integrated currents
    # Use simple trapezoidal integration: I = 2*pi * integral(J * r * dr)
    I_BS_csv = 2 * PI * np.sum(JBS_csv * r * RA * DR * RA) * 1e6
    I_BS_calc = 2 * PI * np.sum(JBS_avg * r * RA * DR * RA) * 1e6
    I_TOT_csv = 2 * PI * np.sum(JTOT_csv * r * RA * DR * RA) * 1e6

    print(f"\n   Integrated currents:")
    print(f"   I_TOT (CSV)  = {I_TOT_csv/1e6:.3f} MA")
    print(f"   I_BS  (CSV)  = {I_BS_csv/1e6:.3f} MA")
    print(f"   I_BS  (calc) = {I_BS_calc/1e6:.3f} MA")
    print(f"   Ratio        = {I_BS_calc/I_BS_csv:.3f}")

    BS_fraction = I_BS_csv / (I_TOT_csv + 1e-20) * 100
    print(f"\n   Bootstrap fraction: {BS_fraction:.1f}%")

    # ============================================================
    # 6. Diagnostics
    # ============================================================
    print(f"\n6. Diagnostic: Key quantities at rho=0.5:")
    idx_half = NRMAX // 2
    print(f"   eps = {eps[idx_half]:.4f}")
    print(f"   ft  = {ft[idx_half]:.4f}")
    print(f"   nu_e* = {nue_star[idx_half]:.4f}")
    print(f"   QP  = {QP[idx_half]:.4f}")
    print(f"   Zeff = {Zeff[idx_half]:.4f}")
    print(f"   L31 = {RL31[idx_half]:.4f}")
    print(f"   L32 = {RL32[idx_half]:.4f}")
    print(f"   L34 = {RL34[idx_half]:.4f}")
    print(f"   alpha = {salfa[idx_half]:.4f}")
    if has_geometry:
        print(f"   TTRHOG = {TTRHOG[idx_half]:.6f}")
        print(f"   RDP = {RDP[idx_half]:.6f}")
        print(f"   ABVRHOG = {ABVRHOG[idx_half]:.6f}")

    # Detailed diagnostic at rho=0.25 (index 12)
    print(f"\n   Detailed diagnostic at rho=0.25 (index 12):")
    i = 12
    if i < NRMAX - 1:
        Pe_avg = 0.5 * (Pe[i+1] + Pe[i])
        Pi_avg = 0.5 * (Pi[i+1] + Pi[i])
        Te_avg = 0.5 * (Te[i+1] + Te[i])
        Ti_avg = 0.5 * (Ti[i+1] + Ti[i])
        RPe_avg = Pe_avg / (Pe_avg + Pi_avg + 1e-20)
        DPe_i = (Pe[i+1] - Pe[i]) / DR
        DPi_i = (Pi[i+1] - Pi[i]) / DR
        DTe_i = (Te[i+1] - Te[i]) / DR
        DTi_i = (Ti[i+1] - Ti[i]) / DR

        print(f"   Pe_avg = {Pe_avg:.4f}, Pi_avg = {Pi_avg:.4f}")
        print(f"   Te_avg = {Te_avg:.4f}, Ti_avg = {Ti_avg:.4f}")
        print(f"   RPe = {RPe_avg:.4f}")
        print(f"   DPe/drho = {DPe_i:.4f}, DPi/drho = {DPi_i:.4f}")
        print(f"   DTe/drho = {DTe_i:.4f}, DTi/drho = {DTi_i:.4f}")
        print(f"   DPe/Pe = {DPe_i/(Pe_avg+1e-20):.4f}, DPi/Pe = {DPi_i/(Pe_avg+1e-20):.4f}")
        print(f"   DTe/Te = {DTe_i/(Te_avg+1e-20):.4f}, DTi/Ti = {DTi_i/(Ti_avg+1e-20):.4f}")

        term1 = RL31[i] * (DPe_i / (Pe_avg + 1e-20) + DPi_i / (Pe_avg + 1e-20))
        term2 = RL32[i] * DTe_i / (Te_avg + 1e-20)
        term3 = RL34[i] * salfa[i] * (1.0 - RPe_avg) / (RPe_avg + 1e-20) * DTi_i / (Ti_avg + 1e-20)
        print(f"   Term1 (L31*dp/p) = {term1:.6f}")
        print(f"   Term2 (L32*dTe/Te) = {term2:.6f}")
        print(f"   Term3 (L34*alpha*...) = {term3:.6f}")
        print(f"   Sum of terms = {term1+term2+term3:.6f}")
        print(f"   JBS_calc[{i}] = {JBS_calc[i]*1e-6:.4f} MA/m^2")
        print(f"   JBS_avg[{i}] = {JBS_avg[i]:.4f} MA/m^2")
        print(f"   JBS_csv[{i}] = {JBS_csv[i]:.4f} MA/m^2")

    # ============================================================
    # 7. Plot results
    # ============================================================
    print("\n7. Generating plots...")

    omfit_rho, omfit_currents = load_omfit_current_profiles()
    if omfit_currents is None:
        print("   OMFIT current profiles not found: figure 1 keeps CSV-only curves")
    else:
        print("   OMFIT current profiles loaded for figure 1 overlay")

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.suptitle('JTOT and JBS Analysis V3: Snapshot-Consistent + TR Edge Derivative', fontsize=14, fontweight='bold')

    # Plot 1: Current components
    ax = axes[0, 0]
    ax.plot(r, JTOT_csv, 'k-', lw=2, label='JTOT')
    ax.plot(r, JOH_csv, 'b-', lw=1.5, label='JOH')
    ax.plot(r, JNB_csv, 'g-', lw=1.5, label='JNB')
    ax.plot(r, JRF_csv, 'c-', lw=1.5, label='JRF')
    ax.plot(r, JBS_csv, 'r-', lw=1.5, label='JBS')
    if omfit_currents is not None and omfit_rho is not None:
        omfit_styles = {
            'tot': ('k--', 'OMFIT tot'),
            'oh': ('b--', 'OMFIT oh'),
            'bs': ('r--', 'OMFIT bs'),
            'ec': ('m--', 'OMFIT ec'),
            'ic': ('y--', 'OMFIT ic'),
            'rf': ('0.4', 'OMFIT rf'),
        }
        for key, (style, label) in omfit_styles.items():
            if key in omfit_currents:
                if key == 'rf':
                    ax.plot(omfit_rho, omfit_currents[key], '--', color=style, lw=1.2, label=label)
                else:
                    ax.plot(omfit_rho, omfit_currents[key], style, lw=1.2, label=label)
    ax.axhline(y=0, color='gray', ls='--', alpha=0.5)
    ax.set_xlabel('rho')
    ax.set_ylabel('J [MA/m²]')
    ax.set_title('Current Density Components (CSV + OMFIT)')
    ax.legend(loc='lower right', fontsize=9)
    ax.grid(True, alpha=0.3)

    # Plot 2: JBS comparison
    ax = axes[0, 1]
    ax.plot(r, JBS_csv, 'r-', lw=2, label='JBS (CSV/TR)')
    ax.plot(r, JBS_avg, 'b--', lw=1.5, label='JBS (recalculated)')
    ax.axhline(y=0, color='gray', ls='--', alpha=0.5)
    ax.set_xlabel('rho')
    ax.set_ylabel('J [MA/m²]')
    ax.set_title('Bootstrap Current: CSV vs Recalculation')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)

    # Plot 3: Ratio
    ax = axes[0, 2]
    ratio_arr = JBS_avg / (JBS_csv + 1e-30)
    ratio_arr[np.abs(JBS_csv) < 1e-10] = np.nan
    ax.plot(r, ratio_arr, 'g-', lw=1.5)
    ax.axhline(y=1.0, color='r', ls='--', alpha=0.7, label='Perfect match')
    ax.set_xlabel('rho')
    ax.set_ylabel('JBS_calc / JBS_csv')
    ax.set_title('Recalculation Ratio')
    ax.set_ylim([0, 2])
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)

    # Plot 4: Sauter coefficients
    ax = axes[1, 0]
    ax.plot(r, RL31, 'b-', lw=1.5, label='L31')
    ax.plot(r, RL32, 'r-', lw=1.5, label='L32')
    ax.plot(r, RL34, 'g-', lw=1.5, label='L34')
    ax.plot(r, salfa, 'm--', lw=1.5, label='α')
    ax.set_xlabel('rho')
    ax.set_ylabel('Coefficient')
    ax.set_title('Sauter Neoclassical Coefficients')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)

    # Plot 5: Geometry quantities
    ax = axes[1, 1]
    if has_geometry:
        ax2 = ax.twinx()
        l1, = ax.plot(r, TTRHOG, 'b-', lw=1.5, label='TTRHOG')
        l2, = ax.plot(r, EPSRHO, 'r-', lw=1.5, label='EPSRHO')
        l3, = ax2.plot(r, RDP, 'g--', lw=1.5, label='RDP')
        ax.set_xlabel('rho')
        ax.set_ylabel('TTRHOG, EPSRHO')
        ax2.set_ylabel('RDP [Wb]')
        ax.set_title('Geometry Quantities from TR')
        lines = [l1, l2, l3]
        ax.legend(lines, [l.get_label() for l in lines], loc='best')
    else:
        ax.plot(r, eps, 'r-', lw=1.5, label='eps (approx)')
        ax.set_xlabel('rho')
        ax.set_ylabel('epsilon')
        ax.set_title('Approximate Geometry')
        ax.legend(loc='best')
    ax.grid(True, alpha=0.3)

    # Plot 6: Collisionality and trapped fraction
    ax = axes[1, 2]
    ax2 = ax.twinx()
    l1, = ax.plot(r, nue_star, 'b-', lw=1.5, label='ν_e*')
    l2, = ax.plot(r, nui_star, 'r-', lw=1.5, label='ν_i*')
    l3, = ax2.plot(r, ft, 'g--', lw=1.5, label='f_t')
    ax.set_xlabel('rho')
    ax.set_ylabel('Collisionality (ν*)')
    ax2.set_ylabel('Trapped fraction')
    ax.set_title('Collisionality and Trapped Fraction')
    ax.set_yscale('log')
    lines = [l1, l2, l3]
    ax.legend(lines, [l.get_label() for l in lines], loc='best')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('jtot_jbs_comparison_v3.png', dpi=150, bbox_inches='tight')
    print("   Plot saved to jtot_jbs_comparison_v3.png")

    # Save data
    df = pd.DataFrame({
        'rho': r,
        'JTOT_csv': JTOT_csv,
        'JBS_csv': JBS_csv,
        'JBS_recalc': JBS_avg,
        'eps': eps,
        'ft': ft,
        'L31': RL31,
        'L32': RL32,
        'L34': RL34,
    })
    if has_geometry:
        df['TTRHOG'] = TTRHOG
        df['RDP'] = RDP
        df['EPSRHO'] = EPSRHO
    df.to_csv('jtot_jbs_comparison_v3.csv', index=False)
    print("   Data saved to jtot_jbs_comparison_v3.csv")

    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)

    if not has_geometry:
        print("""
NOTE: Geometry data not found. To generate:
  1. Run TR code
  2. In graphics menu: G R E
  3. This generates CSV with TTRHOG, RDP, EPSRHO, etc.
  4. Re-run this script
""")


if __name__ == '__main__':
    main()
