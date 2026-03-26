#!/usr/bin/env python3
"""
Verify energy balance in TASK/TR at the final time step.

Uses the integral form: at steady state, integrating the energy equation
from 0 to rho, the cumulative source should equal the total heat flux
(diffusive + convective) through the surface at rho.

This script auto-discovers tr_data_*.csv files by their Title line,
so it is robust to changes in CSV numbering.
"""

import os
import re
import math

BASE = os.path.dirname(os.path.abspath(__file__))
RKEV = 1.602176634e-16  # keV -> J
MW_FACTOR = 1e20 * RKEV / 1e6  # [10^20/m^3 * keV * m^2/s] -> MW


# =====================================================================
#  CSV discovery and reading
# =====================================================================

def build_title_index():
    """Scan all tr_data_*.csv files. Return {number: title_string}."""
    index = {}
    for f in os.listdir(BASE):
        m = re.match(r"tr_data_(\d+)\.csv", f)
        if not m:
            continue
        num = int(m.group(1))
        path = os.path.join(BASE, f)
        with open(path) as fh:
            title = fh.readline().strip()
        index[num] = title
    return index


def find_csv(index, pattern, which="last"):
    """Find csv file number whose title contains `pattern`.
    which='last' returns the highest-numbered match (most recent snapshot).
    """
    matches = sorted(num for num, t in index.items() if pattern in t)
    if not matches:
        raise FileNotFoundError(f"No CSV with title containing '{pattern}'")
    return matches[-1] if which == "last" else matches[0]


def try_find_csv(index, pattern, which="last"):
    matches = sorted(num for num, t in index.items() if pattern in t)
    if not matches:
        return None
    return matches[-1] if which == "last" else matches[0]


def read_csv(num):
    """Read tr_data_NNN.csv. Returns dict of column arrays."""
    path = os.path.join(BASE, f"tr_data_{num:03d}.csv")
    with open(path) as f:
        lines = f.readlines()
    header = [h.strip() for h in lines[1].strip().split(",")]
    data = {h: [] for h in header}
    for line in lines[2:]:
        parts = line.strip().split(",")
        if len(parts) < len(header):
            continue
        for i, h in enumerate(header):
            data[h].append(float(parts[i]))
    return data


def avg_to_grid(arr):
    """Average half-mesh array to grid points; extrapolate at boundary."""
    n = len(arr)
    g = [0.5 * (arr[j] + arr[j + 1]) for j in range(n - 1)]
    g.append(1.5 * arr[-1] - 0.5 * arr[-2])
    return g


def integrate_source_trapezoid(x, vprime, source):
    """Return Q(x) where Q = (1/V') * integral_0^x source*V' dx."""
    n = len(x)
    q = [0.0] * n
    for i in range(n - 1):
        dx = x[i + 1] - x[i]
        dvs = 0.5 * (vprime[i] * source[i] + vprime[i + 1] * source[i + 1]) * dx
        vpq = vprime[i] * q[i] + dvs
        q[i + 1] = vpq / vprime[i + 1] if abs(vprime[i + 1]) > 1e-30 else 0.0
    return q


# =====================================================================
#  Main
# =====================================================================

def main():
    idx = build_title_index()

    # --- Discover files by title (take LAST occurrence = final snapshot) ---
    # Use specific patterns to avoid matching multi-time-slice versions
    f_pin  = find_csv(idx, "@PIN [MW/m$+3$=]")      # single-snapshot PIN (4 species cols)
    f_src  = find_csv(idx, "POH,PNB,PNF,PRSUM")
    f_loss = find_csv(idx, "PRSUM,PRB,PRC,PRL,PCX,PIE,QEI")
    f_T    = find_csv(idx, "T(NS)")
    f_n    = find_csv(idx, "n(NS)")
    f_chie = find_csv(idx, "AKE,AKNCE,AKDWE")
    f_chii = find_csv(idx, "AKD,AKNCD,AKDWD")
    f_dv   = find_csv(idx, "DVRHO vs RHO")
    f_ar1  = find_csv(idx, "AR1RHO vs RHO")
    f_ar2  = find_csv(idx, "AR2RHO vs RHO")
    f_av   = find_csv(idx, "@AV(NS) [m/s]")
    f_avk  = find_csv(idx, "@AVK(NS) [m/s]")
    f_qx   = try_find_csv(idx, "@QX(NS) [MW]")
    f_qtr  = try_find_csv(idx, "@QTR(NS) [MW]")

    print("--- CSV File Mapping (auto-discovered, last snapshot) ---")
    for label, num in [("PIN", f_pin), ("Sources", f_src), ("Losses/QEI", f_loss),
                       ("T(NS)", f_T), ("n(NS)", f_n),
                       ("chi_e", f_chie), ("chi_i", f_chii),
                       ("DVRHO", f_dv), ("AR1RHO", f_ar1), ("AR2RHO", f_ar2),
                       ("AV", f_av), ("AVK", f_avk)]:
        print(f"  {label:12s} -> tr_data_{num:03d}.csv")
    if f_qx is not None:
        print(f"  {'QX':12s} -> tr_data_{f_qx:03d}.csv")
    if f_qtr is not None:
        print(f"  {'QTR':12s} -> tr_data_{f_qtr:03d}.csv")
    print()

    # --- Read data ---
    d_pin  = read_csv(f_pin)
    d_src  = read_csv(f_src)
    d_loss = read_csv(f_loss)
    d_T    = read_csv(f_T)
    d_n    = read_csv(f_n)
    d_chie = read_csv(f_chie)
    d_chii = read_csv(f_chii)
    d_dv   = read_csv(f_dv)
    d_ar1  = read_csv(f_ar1)
    d_ar2  = read_csv(f_ar2)
    d_ad   = read_csv(find_csv(idx, "@AD(NS) [m$+2$=/s]"))
    d_av   = read_csv(f_av)
    d_avk  = read_csv(f_avk)
    d_qx   = read_csv(f_qx) if f_qx is not None else None
    d_qtr  = read_csv(f_qtr) if f_qtr is not None else None

    # Determine nrmax from d_T
    nrmax_len = len(d_T["X"]) if d_T is not None else 50

    # Adjust QX and QTR to match nrmax boundaries
    if d_qx is not None:
        for k in d_qx.keys():
            if len(d_qx[k]) == nrmax_len + 1:
                d_qx[k] = list(d_qx[k][1:])
    if d_qtr is not None:
        for k in d_qtr.keys():
            if len(d_qtr[k]) == nrmax_len + 1:
                d_qtr[k] = list(d_qtr[k][1:])

    # --- Grid info ---
    rho_m = d_T["X"]
    rho_g = d_chie["X"]
    nrmax = len(rho_m)
    dr = rho_m[1] - rho_m[0]

    dvrho  = d_dv["DVRHO"]
    ar1rho = d_ar1["AR1RHO"]
    ar2rho = d_ar2["AR2RHO"]
    dvrhog  = avg_to_grid(dvrho)
    ar1rhog = avg_to_grid(ar1rho)
    ar2rhog = avg_to_grid(ar2rho)

    # QEI
    qei_profile = d_loss["QEI"]
    Q_qei = sum(qei_profile[i] * dvrho[i] * dr for i in range(nrmax))

    # --- Identify AV/AVK column names ---
    ad_cols  = [c for c in d_ad  if c != "X"]
    av_cols  = [c for c in d_av  if c != "X"]
    avk_cols = [c for c in d_avk if c != "X"]

    # --- Species definitions ---
    # PIN_1=e, PIN_2=D, PIN_3=T, PIN_4=He4
    species_list = [
        ("e",   "PIN_1", "TE", "nE", d_chie["AKE"], ad_cols[0], av_cols[0], avk_cols[0]),
        ("D",   "PIN_2", "TD", "nD", d_chii["AKD"], ad_cols[1], av_cols[1], avk_cols[1]),
        ("T",   "PIN_3", "TT", "nT", d_chii["AKD"], ad_cols[2], av_cols[2], avk_cols[2]),
        ("He4", "PIN_4", "TA", "nA", d_chii["AKD"], ad_cols[3], av_cols[3], avk_cols[3]),
    ]

    T_dict = {k: d_T[k] for k in ["TE", "TD", "TT", "TA"]}
    n_dict = {k: d_n[k] for k in ["nE", "nD", "nT", "nA"]}

    CC = 1.5  # = 3/2, matches the code's CC constant

    results = []

    for name, pin_key, T_key, n_key, chi_g, ad_key, av_key, avk_key in species_list:
        pin = list(d_pin[pin_key])
        
        # Target Se_net = PIN_e + QEI, Si_net = PIN_i - QEI  (distribute QEI among ions by density)
        if name == "e":
            pin = [p + qei_profile[j] for j, p in enumerate(pin)]
        else:
            pin_adj = []
            for j in range(nrmax):
                local_ni_tot = d_n["nD"][j] + d_n["nT"][j] + d_n["nA"][j]
                frac = d_n[n_key][j] / local_ni_tot if local_ni_tot > 1e-20 else 0.0
                pin_adj.append(pin[j] - frac * qei_profile[j])
            pin = pin_adj
        T   = T_dict[T_key]
        n   = n_dict[n_key]
        ad_arr  = d_ad[ad_key]
        av_arr  = d_av[av_key]    # on grid-mesh (NRMAX+1 points: 0, dr, 2*dr, ...)
        avk_arr = d_avk[avk_key]  # on grid-mesh

        # chi on grid: skip axis value, use [1:] for grid points j=1..NRMAX
        chi_at_grid = list(chi_g[1:])
        chi_at_grid.append(2 * chi_g[-1] - chi_g[-2])

        # AD, AV, AVK on grid: same structure, skip axis value [0]
        ad_at_grid = list(ad_arr[1:])
        if len(ad_at_grid) < nrmax:
            ad_at_grid.append(2 * ad_arr[-1] - ad_arr[-2])
        av_at_grid = list(av_arr[1:])
        if len(av_at_grid) < nrmax:
            av_at_grid.append(2 * av_arr[-1] - av_arr[-2])
        avk_at_grid = list(avk_arr[1:])
        if len(avk_at_grid) < nrmax:
            avk_at_grid.append(2 * avk_arr[-1] - avk_arr[-2])

        # n, T at grid points (average adjacent half-mesh)
        n_g = avg_to_grid(n)
        T_g = avg_to_grid(T)

        # --- Cumulative source ---
        cum_src = []
        s = 0.0
        for j in range(nrmax):
            s += pin[j] * dvrho[j] * dr
            cum_src.append(s)

        # --- Fluxes at grid points ---
        # Diffusive:  F_diff = V' * <|grad rho|^2> * n * chi * (-dT/drho) * MW_FACTOR
        # Cross:      F_cross = - V' * <|grad rho|^2> * (3/2) * D * T * (dn/drho) * MW_FACTOR
        # Convective: F_conv = V' * <|grad rho|> * (avk + CC*av) * n*T * MW_FACTOR
        flux_diff = []
        flux_cross = []
        flux_conv = []
        flux_total = []
        for j in range(nrmax):
            if j < nrmax - 1:
                dTdr = (T[j + 1] - T[j]) / dr
                dndr = (n[j + 1] - n[j]) / dr
            else:
                dTdr = (T[j] - T[j - 1]) / dr
                dndr = (n[j] - n[j - 1]) / dr

            fd = dvrhog[j] * ar2rhog[j] * n_g[j] *chi_at_grid[j] * (-dTdr) * MW_FACTOR
            fx = -dvrhog[j] * ar2rhog[j] * CC * ad_at_grid[j] * T_g[j] * dndr * MW_FACTOR
            fc = dvrhog[j] * ar1rhog[j] * (avk_at_grid[j] + CC * av_at_grid[j]) \
                 * n_g[j] * T_g[j] * MW_FACTOR
            flux_diff.append(fd)
            flux_cross.append(fx)
            flux_conv.append(fc)
            flux_total.append(fd + fx + fc)

        results.append({
            "name": name,
            "q_src":  cum_src[-1],
            "q_diff": flux_diff[-1],
            "q_cross": flux_cross[-1],
            "q_conv": flux_conv[-1],
            "q_total": flux_total[-1],
            "cum_src": cum_src,
            "flux_diff": flux_diff,
            "flux_cross": flux_cross,
            "flux_conv": flux_conv,
            "flux_total": flux_total,
        })

    if d_qx is not None and d_qtr is not None:
        qx_cols = [c for c in d_qx if c != "X"]
        qtr_cols = [c for c in d_qtr if c != "X"]
        for i, r in enumerate(results):
            r["qx_csv"] = d_qx[qx_cols[i]][-1]
            r["qtr_csv"] = d_qtr[qtr_cols[i]][-1]

    # --- Integrated source components ---
    poh = d_src["POH"]
    pnb = d_src["PNB"]
    pnf = d_src["PNF"]
    neg_prs = d_src["PRSUM"]
    Q_poh = sum(poh[i] * dvrho[i] * dr for i in range(nrmax))
    Q_pnb = sum(pnb[i] * dvrho[i] * dr for i in range(nrmax))
    Q_pnf = sum(pnf[i] * dvrho[i] * dr for i in range(nrmax))
    Q_rad = sum((neg_prs[i]) * dvrho[i] * dr for i in range(nrmax))

    # --- Optional OMFIT integrated source flux ---
    omfit_rho = None
    omfit_qe_power_mw = None
    omfit_qi_power_mw = None
    omfit_qtot_power_mw = None
    omfit_qe_int = None
    omfit_qi_int = None
    omfit_qtot_int = None
    omfit_qe_grad = None
    omfit_qi_grad = None
    omfit_drmin_drho = None
    omfit_ne = None
    omfit_ni = None
    omfit_chie = None
    omfit_chii = None
    omfit_dte_dr = None
    omfit_dti_dr = None
    try:
        import pandas as pd
        import numpy as np

        omfit_base = "/Users/dengxiaoya/CFEDRSW/OMFIT_out"
        all_prof = pd.read_csv(os.path.join(omfit_base, "outputs_new", "all_profiles.csv"))
        vol_df = pd.read_csv(os.path.join(omfit_base, "outputs_new", "input_profiles_extra.csv"))
        src_path = os.path.join(omfit_base, "rho_rmin_Si_Se_qe_qi_sum.csv")
        if os.path.exists(src_path):
            src_df = pd.read_csv(src_path)
            omfit_rho = all_prof["rho"].values
            omfit_r = all_prof["rmin"].values
            omfit_vprime = vol_df["EXPRO_volp"].values  # dV/drmin [m^2]
            omfit_drmin_drho = np.gradient(omfit_r, omfit_rho)
            Se = src_df["qe_sum"].values               # [W/m^3]
            Si = src_df["qi_sum"].values               # [W/m^3]

            plasma_df = pd.read_csv(os.path.join(omfit_base, "outputs_new", "plasma.csv"))
            chi_df = pd.read_csv(os.path.join(omfit_base, "chi_transport_coefficients.csv"))

            e_charge = 1.602176634e-19
            ne = plasma_df["ne"].values * 1e19
            te_ev = plasma_df["Te"].values * 1000.0
            if "ni_2" in all_prof.columns and "ni_1" in all_prof.columns:
                ni = (all_prof["ni_2"].values + all_prof["ni_1"].values) * 1e19
            else:
                ni = plasma_df["ni"].values * 1e19 if "ni" in plasma_df.columns else ne
            if "Ti_2" in all_prof.columns:
                ti_ev = all_prof["Ti_2"].values * 1000.0
            elif "Ti_1" in all_prof.columns:
                ti_ev = all_prof["Ti_1"].values * 1000.0
            else:
                ti_ev = te_ev
            chi_e_omfit = chi_df["chi_e_m2s"].values
            chi_i_omfit = chi_df["chi_i_m2s"].values
            dte_dr = np.gradient(te_ev, omfit_r)
            dti_dr = np.gradient(ti_ev, omfit_r)
            omfit_qe_grad = -e_charge * ne * chi_e_omfit * dte_dr
            omfit_qi_grad = -e_charge * ni * chi_i_omfit * dti_dr
            omfit_ne = ne
            omfit_ni = ni
            omfit_chie = chi_e_omfit
            omfit_chii = chi_i_omfit
            omfit_dte_dr = dte_dr
            omfit_dti_dr = dti_dr

            # Qe_int/Qi_int are heat-flux-density-like quantities [W/m^2]
            omfit_qe_int = integrate_source_trapezoid(omfit_r, omfit_vprime, Se)
            omfit_qi_int = integrate_source_trapezoid(omfit_r, omfit_vprime, Si)
            omfit_qtot_int = [omfit_qe_int[i] + omfit_qi_int[i] for i in range(len(omfit_qe_int))]

            # Convert to total power through the surface: (dV/dr) * Q  [W]
            omfit_qe_power_mw = [omfit_vprime[i] * omfit_qe_int[i] / 1e6 for i in range(len(omfit_r))]
            omfit_qi_power_mw = [omfit_vprime[i] * omfit_qi_int[i] / 1e6 for i in range(len(omfit_r))]
            omfit_qtot_power_mw = [omfit_qe_power_mw[i] + omfit_qi_power_mw[i] for i in range(len(omfit_r))]
            # --- Generate chi_div_grad2_drdrho.dat ---
            # AKE_new = chi_OMFIT / (AR2 * dr_min/drho)
            # so that TASK/TR flux = V'*AR2*n*AKE_new*(-dT/drho) matches OMFIT flux
            rho_plot_tr = [rho_m[j] + dr / 2 if j < nrmax - 1 else 1.0 for j in range(nrmax)]
            ar2_on_omfit = np.interp(omfit_rho, rho_plot_tr, ar2rhog)
            # denom = ar2_on_omfit * omfit_drmin_drho* omfit_drmin_drho
            # # Avoid division by zero near axis
            # denom_safe = np.where(np.abs(denom) > 1e-30, denom, 1e-30)
            # chi_e_for_tr = chi_e_omfit / denom_safe
            # chi_i_for_tr = chi_i_omfit / denom_safe
            # # At axis (rho=0), AR2->0 causes singularity; extrapolate from neighbors
            # if abs(ar2_on_omfit[0]) < 1e-10:
            #     chi_e_for_tr[0] = chi_e_for_tr[1]
            #     chi_i_for_tr[0] = chi_i_for_tr[1]

            # dat_path = os.path.join(BASE, "chi_div_grad2_drdrho.dat")
            # with open(dat_path, "w") as f:
            #     f.write("rho ne ni chi_e chi_i\n")
            #     for k in range(len(omfit_rho)):
            #         f.write(f"{omfit_rho[k]:.8e} {ne[k]:.8e} {ni[k]:.8e} "
            #                 f"{chi_e_for_tr[k]:.8e} {chi_i_for_tr[k]:.8e}\n")
            # print(f"chi_div_grad2_drdrho.dat written: {dat_path}")
            # print(f"  AKE_new = chi_OMFIT / (AR2 * dr_min/drho)")
            # print(f"  so that V'*AR2*n*AKE_new*(-dT/drho) = OMFIT flux")
            # print()

    except Exception as e:
        print(f"OMFIT source integration skipped: {e}")

    # ===== Print Report =====
    print("=" * 82)
    print("  TASK/TR Energy Balance Verification (with Convection)")
    print("=" * 82)
    print(f"  Grid: NRMAX={nrmax}, DR={dr:.4f}, rho=[{rho_m[0]:.2f}..{rho_m[-1]:.2f}]")
    print()

    print("--- Volume-Integrated Source Components [MW] ---")
    print(f"  POH   = {Q_poh:10.3f}")
    print(f"  PNB   = {Q_pnb:10.3f}")
    print(f"  PNF   = {Q_pnf:10.3f}")
    print(f"  PRSUM = {Q_rad:10.3f}  (radiation loss)")
    print(f"  QEI   = {Q_qei:10.3f}  (negative = e heats i)")
    print()

    sS = sum(r["q_src"]   for r in results)
    sD = sum(r["q_diff"]  for r in results)
    sC = sum(r["q_conv"]  for r in results)
    sT = sum(r["q_total"] for r in results)

    print("--- Per-Species Energy Balance at Boundary [MW] ---")
    hdr = f"  {'Spec':<5s} {'Q_src':>10s} {'Q_diff':>10s} {'Q_conv':>10s} {'Q_total':>10s} {'Resid':>10s} {'%':>7s}"
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))
    for r in results:
        resid = r["q_src"] - r["q_total"]
        pct = resid / r["q_src"] * 100 if abs(r["q_src"]) > 0.01 else float("nan")
        print(f"  {r['name']:<5s} {r['q_src']:10.3f} {r['q_diff']:10.3f} "
              f"{r['q_conv']:10.3f} {r['q_total']:10.3f} {resid:10.3f} {pct:7.1f}")
    print("  " + "-" * (len(hdr) - 2))
    resid_t = sS - sT
    pct_t = resid_t / sS * 100 if abs(sS) > 0.01 else float("nan")
    print(f"  {'Total':<5s} {sS:10.3f} {sD:10.3f} {sC:10.3f} {sT:10.3f} {resid_t:10.3f} {pct_t:7.1f}")
    print()
    print("  Q_src   = integral{ PIN_s * V' drho }")
    print("  Q_diff  = V'*<|grad rho|^2>*n*chi*(-dT/drho)  at boundary")
    print("  Q_conv  = V'*<|grad rho|>*(AVK + 3/2*AV)*n*T  at boundary")
    print("  Q_total = Q_diff + Q_conv")
    print("  Resid   = Q_src - Q_total  (= time deriv + numerical + other)")
    if d_qx is not None and d_qtr is not None:
        print()
        print("--- Direct CSV Flux Output Check [MW] ---")
        print(f"  {'Spec':<5s} {'QX(csv)':>10s} {'QTR(csv)':>10s}")
        print(f"  {'-'*5} {'-'*10} {'-'*10}")
        for r in results:
            print(f"  {r['name']:<5s} {r['qx_csv']:10.3f} {r['qtr_csv']:10.3f}")
    print()

    # --- Radial profile ---
    print("--- Radial Profile: Cumulative Source vs Total Flux (all species) ---")
    hdr2 = f"  {'rho':>6s} {'Q_src':>10s} {'Q_recon':>10s} {'QTR_csv':>10s} {'Resid(rec)':>12s} {'Resid(csv)':>12s}"
    print(hdr2)
    print(f"  {'-'*6} {'-'*10} {'-'*10} {'-'*10} {'-'*12} {'-'*12}")
    for j in list(range(0, nrmax, 5)) + [nrmax - 1]:
        rg = rho_m[j] + dr / 2 if j < nrmax - 1 else 1.0
        qs = sum(r["cum_src"][j]    for r in results)
        qt = sum(r["flux_total"][j] for r in results)
        qcsv = 0.0
        if d_qtr is not None:
            qtr_cols = [c for c in d_qtr if c != "X"]
            qcsv = sum(d_qtr[c][j] for c in qtr_cols)
        print(f"  {rg:6.2f} {qs:10.3f} {qt:10.3f} {qcsv:10.3f} {qs-qt:12.3f} {qs-qcsv:12.3f}")
    print()

    # --- CSV output ---
    outpath = os.path.join(BASE, "energy_balance_check.csv")
    with open(outpath, "w") as f:
        f.write("rho_grid,"
                "cum_src_e,cum_src_D,cum_src_T,cum_src_A,cum_src_total,"
                "flux_diff_e,flux_diff_D,flux_diff_T,flux_diff_A,flux_diff_total,"
                "flux_conv_e,flux_conv_D,flux_conv_T,flux_conv_A,flux_conv_total,"
                "flux_total_e,flux_total_D,flux_total_T,flux_total_A,flux_total_total,"
                "resid_e,resid_D,resid_T,resid_A,resid_total\n")
        for j in range(nrmax):
            rg = rho_m[j] + dr / 2 if j < nrmax - 1 else 1.0
            cs = [r["cum_src"][j]    for r in results]
            fd = [r["flux_diff"][j]  for r in results]
            fc = [r["flux_conv"][j]  for r in results]
            ft = [r["flux_total"][j] for r in results]
            rs = [cs[k] - ft[k] for k in range(4)]
            vals = [f"{rg:.4f}"]
            for arr in [cs, fd, fc, ft, rs]:
                for v in arr:
                    vals.append(f"{v:.6e}")
                vals.append(f"{sum(arr):.6e}")
            f.write(",".join(vals) + "\n")
    print(f"Detailed CSV: {outpath}")

    # --- Compact comparison table for source vs total transport ---
    table_path = os.path.join(BASE, "qsource_vs_qtr_table.csv")
    with open(table_path, "w") as f:
        f.write("rho_grid,Q_source_total_MW,QTR_reconstructed_total_MW,QTR_csv_total_MW,Residual_reconstructed_MW,Residual_csv_MW,Relative_error_reconstructed_percent,Relative_error_csv_percent\n")
        for j in range(nrmax):
            rg = rho_m[j] + dr / 2 if j < nrmax - 1 else 1.0
            qs = sum(r["cum_src"][j] for r in results)
            qt = sum(r["flux_total"][j] for r in results)
            qcsv = 0.0
            if d_qtr is not None:
                qtr_cols = [c for c in d_qtr if c != "X"]
                qcsv = sum(d_qtr[c][j] for c in qtr_cols)
            resid = qs - qt
            resid_csv = qs - qcsv
            rel = resid / qs * 100 if abs(qs) > 1e-12 else 0.0
            rel_csv = resid_csv / qs * 100 if abs(qs) > 1e-12 else 0.0
            f.write(f"{rg:.4f},{qs:.6e},{qt:.6e},{qcsv:.6e},{resid:.6e},{resid_csv:.6e},{rel:.6e},{rel_csv:.6e}\n")
    print(f"Table CSV: {table_path}")

    # --- Plot source vs total transport ---
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        rho_plot = []
        qsrc_plot = []
        qtr_recon_plot = []
        qtr_csv_plot = []
        resid_plot = []
        resid_csv_plot = []
        for j in range(nrmax):
            rg = rho_m[j] + dr / 2 if j < nrmax - 1 else 1.0
            qs = sum(r["cum_src"][j] for r in results)
            qt = sum(r["flux_total"][j] for r in results)
            qcsv = 0.0
            if d_qtr is not None:
                qtr_cols = [c for c in d_qtr if c != "X"]
                qcsv = sum(d_qtr[c][j] for c in qtr_cols)
            rho_plot.append(rg)
            qsrc_plot.append(qs)
            qtr_recon_plot.append(qt)
            qtr_csv_plot.append(qcsv)
            resid_plot.append(qs - qt)
            resid_csv_plot.append(qs - qcsv)

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True,
                                       gridspec_kw={"height_ratios": [3, 1]})

        ax1.plot(rho_plot, qsrc_plot, "o-", lw=2, ms=4, label="Q_source(rho)")
        ax1.plot(rho_plot, qtr_recon_plot, "s-", lw=2, ms=4, label="QTR reconstructed")
        if d_qtr is not None:
            ax1.plot(rho_plot, qtr_csv_plot, "^-", lw=2, ms=4, label="QTR from CSV")
        if omfit_rho is not None and omfit_qe_power_mw is not None:
            ax1.plot(omfit_rho, omfit_qe_power_mw, "--", lw=1.8, color="tab:red",
                     label="OMFIT Qe,int -> power")
        if omfit_rho is not None and omfit_qtot_power_mw is not None:
            ax1.plot(omfit_rho, omfit_qtot_power_mw, ":", lw=2.0, color="tab:red",
                     label="OMFIT (Qe,int+Qi,int) -> power")
        ax1.set_ylabel("Power [MW]")
        ax1.set_title("Q_source(rho) vs QTR_total(rho)")
        ax1.grid(True, alpha=0.3)
        ax1.legend()

        ax2.plot(rho_plot, resid_plot, "d-", color="tab:red", lw=1.8, ms=3, label="Residual (reconstructed)")
        if d_qtr is not None:
            ax2.plot(rho_plot, resid_csv_plot, "x-", color="tab:green", lw=1.5, ms=4, label="Residual (CSV)")
        ax2.axhline(0.0, color="k", lw=1, alpha=0.5)
        ax2.set_xlabel("rho")
        ax2.set_ylabel("Residual [MW]")
        ax2.grid(True, alpha=0.3)
        ax2.legend()

        plot_path = os.path.join(BASE, "qsource_vs_qtr_plot.png")
        plt.tight_layout()
        plt.savefig(plot_path, dpi=160)
        plt.close(fig)
        print(f"Plot PNG: {plot_path}")

        # --- V'-normalized comparison (using dV/drmin) ---
        if omfit_rho is not None and omfit_drmin_drho is not None:
            import numpy as np

            drmin_drho_tr = np.interp(rho_plot, omfit_rho, omfit_drmin_drho)
            dv_drmin_tr = []
            qsrc_over_vp = []
            qtr_recon_over_vp = []
            qtr_csv_over_vp = []
            for j in range(nrmax):
                dvdr = dvrhog[j] / drmin_drho_tr[j] if abs(drmin_drho_tr[j]) > 1e-30 else float("nan")
                dv_drmin_tr.append(dvdr)
                qsrc_over_vp.append(qsrc_plot[j] / dvdr if abs(dvdr) > 1e-30 else float("nan"))
                qtr_recon_over_vp.append(qtr_recon_plot[j] / dvdr if abs(dvdr) > 1e-30 else float("nan"))
                qtr_csv_over_vp.append(qtr_csv_plot[j] / dvdr if abs(dvdr) > 1e-30 else float("nan"))

            fig_vp, ax_vp = plt.subplots(1, 1, figsize=(8, 5.5))
            ax_vp.plot(rho_plot, qsrc_over_vp, "o-", lw=2, ms=4, label="TR Q_source / (dV/drmin)")
            ax_vp.plot(rho_plot, qtr_recon_over_vp, "s-", lw=2, ms=4, label="TR QTR reconstructed / (dV/drmin)")
            if d_qtr is not None:
                ax_vp.plot(rho_plot, qtr_csv_over_vp, "^-", lw=2, ms=4, label="TR QTR CSV / (dV/drmin)")
            if omfit_rho is not None and omfit_qe_int is not None:
                ax_vp.plot(omfit_rho, np.array(omfit_qe_int) / 1e6, "--", lw=1.8, color="tab:red",
                           label="OMFIT Qe,int")
            if omfit_rho is not None and omfit_qtot_int is not None:
                ax_vp.plot(omfit_rho, np.array(omfit_qtot_int) / 1e6, ":", lw=2.0, color="tab:red",
                           label="OMFIT (Qe,int+Qi,int)")
            ax_vp.set_xlabel("rho")
            ax_vp.set_ylabel("Power / (dV/drmin) [MW/m$^2$]")
            ax_vp.set_title("V'-Normalized Comparison using dV/drmin")
            ax_vp.grid(True, alpha=0.3)
            ax_vp.legend()
            plot_path_vp = os.path.join(BASE, "qsource_vs_qtr_div_vprime_plot.png")
            plt.tight_layout()
            plt.savefig(plot_path_vp, dpi=160)
            plt.close(fig_vp)
            print(f"Vprime-normalized Plot PNG: {plot_path_vp}")

            # --- Species-separated V'-normalized comparison ---
            qe_src_over_vp = []
            qe_recon_over_vp = []
            qe_csv_over_vp = []
            qe_diff_over_vp = []
            qe_cross_over_vp = []
            qe_conv_over_vp = []
            qi_src_over_vp = []
            qi_recon_over_vp = []
            qi_csv_over_vp = []
            qi_diff_over_vp = []
            qi_cross_over_vp = []
            qi_conv_over_vp = []
            for j in range(nrmax):
                dvdr = dv_drmin_tr[j]
                qe_src_over_vp.append(results[0]["cum_src"][j] / dvdr if abs(dvdr) > 1e-30 else float("nan"))
                qe_recon_over_vp.append(results[0]["flux_total"][j] / dvdr if abs(dvdr) > 1e-30 else float("nan"))
                qe_diff_over_vp.append(results[0]["flux_diff"][j] / dvdr if abs(dvdr) > 1e-30 else float("nan"))
                qe_cross_over_vp.append(results[0]["flux_cross"][j] / dvdr if abs(dvdr) > 1e-30 else float("nan"))
                qe_conv_over_vp.append(results[0]["flux_conv"][j] / dvdr if abs(dvdr) > 1e-30 else float("nan"))

                qi_src = results[1]["cum_src"][j] + results[2]["cum_src"][j] + results[3]["cum_src"][j]
                qi_tot = results[1]["flux_total"][j] + results[2]["flux_total"][j] + results[3]["flux_total"][j]
                qi_diff = results[1]["flux_diff"][j] + results[2]["flux_diff"][j] + results[3]["flux_diff"][j]
                qi_cross = results[1]["flux_cross"][j] + results[2]["flux_cross"][j] + results[3]["flux_cross"][j]
                qi_conv = results[1]["flux_conv"][j] + results[2]["flux_conv"][j] + results[3]["flux_conv"][j]
                qi_src_over_vp.append(qi_src / dvdr if abs(dvdr) > 1e-30 else float("nan"))
                qi_recon_over_vp.append(qi_tot / dvdr if abs(dvdr) > 1e-30 else float("nan"))
                qi_diff_over_vp.append(qi_diff / dvdr if abs(dvdr) > 1e-30 else float("nan"))
                qi_cross_over_vp.append(qi_cross / dvdr if abs(dvdr) > 1e-30 else float("nan"))
                qi_conv_over_vp.append(qi_conv / dvdr if abs(dvdr) > 1e-30 else float("nan"))

                if d_qtr is not None:
                    qtr_cols = [c for c in d_qtr if c != "X"]
                    if len(qtr_cols) >= 4:
                        qcsv_e = d_qtr[qtr_cols[0]][j]
                        qcsv_i = d_qtr[qtr_cols[1]][j] + d_qtr[qtr_cols[2]][j] + d_qtr[qtr_cols[3]][j]
                        qe_csv_over_vp.append(qcsv_e / dvdr if abs(dvdr) > 1e-30 else float("nan"))
                        qi_csv_over_vp.append(qcsv_i / dvdr if abs(dvdr) > 1e-30 else float("nan"))
                    else:
                        qe_csv_over_vp.append(float("nan"))
                        qi_csv_over_vp.append(float("nan"))
                else:
                    qe_csv_over_vp.append(float("nan"))
                    qi_csv_over_vp.append(float("nan"))

            # --- Prepare Chi plot data ---
            import numpy as np
            chi_e_g = list(d_chie["AKE"][1:]) + [2 * d_chie["AKE"][-1] - d_chie["AKE"][-2]]
            aknce_g = list(d_chie["AKNCE"][1:]) + [2 * d_chie["AKNCE"][-1] - d_chie["AKNCE"][-2]]
            akdwe_g = list(d_chie["AKDWE"][1:]) + [2 * d_chie["AKDWE"][-1] - d_chie["AKDWE"][-2]]
            
            chi_i_g = list(d_chii["AKD"][1:]) + [2 * d_chii["AKD"][-1] - d_chii["AKD"][-2]]
            akncd_g = list(d_chii["AKNCD"][1:]) + [2 * d_chii["AKNCD"][-1] - d_chii["AKNCD"][-2]]
            akdwd_g = list(d_chii["AKDWD"][1:]) + [2 * d_chii["AKDWD"][-1] - d_chii["AKDWD"][-2]]
            
            chi_src_e_eff = []
            chi_src_i_eff = []
            chi_csv_e_eff = []
            chi_csv_i_eff = []
            te_tr_loc = d_T["TE"]
            td_tr_loc = d_T["TD"]
            tt_tr_loc = d_T["TT"]
            ta_tr_loc = d_T["TA"]
            ne_tr_g_loc = avg_to_grid(d_n["nE"])
            ni_tr_g_loc = avg_to_grid([d_n["nD"][j] + d_n["nT"][j] + d_n["nA"][j] for j in range(nrmax)])
            
            for j in range(nrmax):
                if j < nrmax - 1:
                    dte_drho = (te_tr_loc[j + 1] - te_tr_loc[j]) / dr
                    dti_drho_num = (avg_to_grid(d_n["nD"])[j] * (td_tr_loc[j+1]-td_tr_loc[j])/dr
                                  + avg_to_grid(d_n["nT"])[j] * (tt_tr_loc[j+1]-tt_tr_loc[j])/dr
                                  + avg_to_grid(d_n["nA"])[j] * (ta_tr_loc[j+1]-ta_tr_loc[j])/dr)
                else:
                    dte_drho = (te_tr_loc[j] - te_tr_loc[j - 1]) / dr
                    dti_drho_num = (avg_to_grid(d_n["nD"])[j] * (td_tr_loc[j]-td_tr_loc[j-1])/dr
                                  + avg_to_grid(d_n["nT"])[j] * (tt_tr_loc[j]-tt_tr_loc[j-1])/dr
                                  + avg_to_grid(d_n["nA"])[j] * (ta_tr_loc[j]-ta_tr_loc[j-1])/dr)
                dti_drho = dti_drho_num / ni_tr_g_loc[j] if abs(ni_tr_g_loc[j]) > 1e-30 else 0.0
                
                denom_chi_e = ne_tr_g_loc[j] * (-dte_drho) * ar2rhog[j] * drmin_drho_tr[j] * MW_FACTOR
                chi_e_val = qe_src_over_vp[j] / denom_chi_e if abs(denom_chi_e) > 1e-30 else float('nan')
                chi_src_e_eff.append(chi_e_val)
                if len(qe_csv_over_vp) > j:
                    chi_csv_e_eff.append(qe_csv_over_vp[j] / denom_chi_e if abs(denom_chi_e) > 1e-30 else float('nan'))
                
                denom_chi_i = ni_tr_g_loc[j] * (-dti_drho) * ar2rhog[j] * drmin_drho_tr[j] * MW_FACTOR
                chi_i_val = qi_src_over_vp[j] / denom_chi_i if abs(denom_chi_i) > 1e-30 else float('nan')
                chi_src_i_eff.append(chi_i_val)
                if len(qi_csv_over_vp) > j:
                    chi_csv_i_eff.append(qi_csv_over_vp[j] / denom_chi_i if abs(denom_chi_i) > 1e-30 else float('nan'))

            fig_vp_sp, axes = plt.subplots(2, 2, figsize=(16, 10), sharex=True)
            ax_vpe = axes[0, 0]
            ax_vpi = axes[1, 0]
            ax_che = axes[0, 1]
            ax_chi = axes[1, 1]

            ax_vpe.plot(rho_plot, qe_src_over_vp, "o-", lw=2, ms=4, color="black", label="TR Qsrc_e / (dV/drmin)")
            if d_qtr is not None:
                ax_vpe.plot(rho_plot, qe_csv_over_vp, "^-", lw=2, ms=5, color="magenta", label="TR QTR_e CSV / (dV/drmin)")
            ax_vpe.plot(rho_plot, qe_recon_over_vp, "s-", lw=2, ms=4, label="TR QTR_e recon total / (dV/drmin)")
            ax_vpe.plot(rho_plot, qe_diff_over_vp, "-.", lw=1.8, color="tab:blue", label="TR Qdiff_e recon / (dV/drmin)")
            ax_vpe.plot(rho_plot, qe_cross_over_vp, ":", lw=2.0, color="tab:green", label="TR Qcross_e / (dV/drmin)")
            ax_vpe.plot(rho_plot, qe_conv_over_vp, "--", lw=1.6, color="tab:gray", label="TR Qconv_e / (dV/drmin)")
            if omfit_rho is not None and omfit_qe_int is not None:
                ax_vpe.plot(omfit_rho, np.array(omfit_qe_int) / 1e6, "--", lw=1.8, color="tab:red",
                            label="OMFIT Qe,int")
            if omfit_rho is not None and omfit_qe_grad is not None:
                ax_vpe.plot(omfit_rho, np.array(omfit_qe_grad) / 1e6, "-.", lw=1.8, color="firebrick",
                            label="OMFIT Qe,eff=-enchi dT/dr")
            ax_vpe.set_ylabel("Power / (dV/drmin) [MW/m$^2$]")
            ax_vpe.set_title("Electron V'-Normalized Comparison")
            ax_vpe.grid(True, alpha=0.3)
            ax_vpe.legend()

            ax_vpi.plot(rho_plot, qi_src_over_vp, "o-", lw=2, ms=4, color="black", label="TR Qsrc_i / (dV/drmin)")
            if d_qtr is not None:
                ax_vpi.plot(rho_plot, qi_csv_over_vp, "^-", lw=2, ms=5, color="magenta", label="TR QTR_i CSV / (dV/drmin)")
            ax_vpi.plot(rho_plot, qi_recon_over_vp, "s-", lw=2, ms=4, label="TR QTR_i recon total / (dV/drmin)")
            ax_vpi.plot(rho_plot, qi_diff_over_vp, "-.", lw=1.8, color="tab:blue", label="TR Qdiff_i recon / (dV/drmin)")
            ax_vpi.plot(rho_plot, qi_cross_over_vp, ":", lw=2.0, color="tab:green", label="TR Qcross_i / (dV/drmin)")
            ax_vpi.plot(rho_plot, qi_conv_over_vp, "--", lw=1.6, color="tab:gray", label="TR Qconv_i / (dV/drmin)")
            if omfit_rho is not None and omfit_qi_int is not None:
                ax_vpi.plot(omfit_rho, np.array(omfit_qi_int) / 1e6, "--", lw=1.8, color="tab:red",
                            label="OMFIT Qi,int")
            if omfit_rho is not None and omfit_qi_grad is not None:
                ax_vpi.plot(omfit_rho, np.array(omfit_qi_grad) / 1e6, "-.", lw=1.8, color="firebrick",
                            label="OMFIT Qi,eff=-enchi dT/dr")
            ax_vpi.set_xlabel("rho")
            ax_vpi.set_ylabel("Power / (dV/drmin) [MW/m$^2$]")
            ax_vpi.set_title("Ion V'-Normalized Comparison (D+T+He4)")
            ax_vpi.grid(True, alpha=0.3)
            ax_vpi.legend()

            ax_che.plot(rho_plot, chi_e_g, "s-", lw=2, ms=4, label="TR AKE (Total $\chi_e$)")
            ax_che.plot(rho_plot, aknce_g, ":", lw=1.8, label="TR AKNCE (Neoclassical $\chi_e$)")
            ax_che.plot(rho_plot, akdwe_g, "--", lw=1.8, label="TR AKDWE (Anomalous $\chi_e$)")
            ax_che.plot(rho_plot, chi_src_e_eff, "o-", lw=2, ms=4, color="black", label=r"Derived $\chi_e$ from Q_src")
            if d_qtr is not None:
                ax_che.plot(rho_plot, chi_csv_e_eff, "^-", lw=2, ms=5, color="magenta", label=r"Derived $\chi_e$ from QTR_csv")
            ax_che.set_ylabel("Transport Coefficient [m$^2$/s]")
            ax_che.set_title("Electron Transport Coefficients")
            ax_che.grid(True, alpha=0.3)
            ymax_che = max(1.0, np.nanmax([x for x in chi_e_g if not np.isnan(x) and x < 100]) * 1.5)
            ax_che.set_ylim(0, ymax_che)
            ax_che.legend()
            
            ax_chi.plot(rho_plot, chi_i_g, "s-", lw=2, ms=4, label="TR AKD (Total $\chi_i$)")
            ax_chi.plot(rho_plot, akncd_g, ":", lw=1.8, label="TR AKNCD (Neoclassical $\chi_i$)")
            ax_chi.plot(rho_plot, akdwd_g, "--", lw=1.8, label="TR AKDWD (Anomalous $\chi_i$)")
            ax_chi.plot(rho_plot, chi_src_i_eff, "o-", lw=2, ms=4, color="black", label=r"Derived $\chi_i$ from Q_src")
            if d_qtr is not None:
                ax_chi.plot(rho_plot, chi_csv_i_eff, "^-", lw=2, ms=5, color="magenta", label=r"Derived $\chi_i$ from QTR_csv")
            ax_chi.set_xlabel("rho")
            ax_chi.set_ylabel("Transport Coefficient [m$^2$/s]")
            ax_chi.set_title("Ion Transport Coefficients")
            ax_chi.grid(True, alpha=0.3)
            ymax_chi = max(1.0, np.nanmax([x for x in chi_i_g if not np.isnan(x) and x < 100]) * 1.5)
            ax_chi.set_ylim(0, ymax_chi)
            ax_chi.legend()

            plot_path_vp_sp = os.path.join(BASE, "qsource_vs_qtr_div_vprime_species_plot.png")
            plt.tight_layout()
            plt.savefig(plot_path_vp_sp, dpi=160)
            plt.close(fig_vp_sp)
            print(f"Vprime-normalized Species Plot PNG: {plot_path_vp_sp}")

            # --- Detailed Qeff vs TR Qdiff comparison ---
            e_charge = 1.602176634e-19
            ne_tr_g = avg_to_grid(d_n["nE"])
            ni_tr_g = avg_to_grid([d_n["nD"][j] + d_n["nT"][j] + d_n["nA"][j] for j in range(nrmax)])
            chi_e_tr_g = list(d_chie["AKE"][1:]) + [2 * d_chie["AKE"][-1] - d_chie["AKE"][-2]]
            chi_i_tr_g = list(d_chii["AKD"][1:]) + [2 * d_chii["AKD"][-1] - d_chii["AKD"][-2]]

            te_tr = d_T["TE"]
            td_tr = d_T["TD"]
            tt_tr = d_T["TT"]
            ta_tr = d_T["TA"]
            dte_drmin_tr = []
            dte_drho_tr = []
            dti_eff_drmin_tr = []
            dti_eff_drho_tr = []
            qdiff_e_over_vp = []
            qdiff_i_over_vp = []
            geom_e = []
            geom_i = []
            for j in range(nrmax):
                if j < nrmax - 1:
                    dte_drho = (te_tr[j + 1] - te_tr[j]) / dr
                    dtd_drho = (td_tr[j + 1] - td_tr[j]) / dr
                    dtt_drho = (tt_tr[j + 1] - tt_tr[j]) / dr
                    dta_drho = (ta_tr[j + 1] - ta_tr[j]) / dr
                else:
                    dte_drho = (te_tr[j] - te_tr[j - 1]) / dr
                    dtd_drho = (td_tr[j] - td_tr[j - 1]) / dr
                    dtt_drho = (tt_tr[j] - tt_tr[j - 1]) / dr
                    dta_drho = (ta_tr[j] - ta_tr[j - 1]) / dr
                drdrho = drmin_drho_tr[j]
                dte_dr_m = dte_drho / drdrho if abs(drdrho) > 1e-30 else float("nan")
                dti_num = (avg_to_grid(d_n["nD"])[j] * dtd_drho + avg_to_grid(d_n["nT"])[j] * dtt_drho + avg_to_grid(d_n["nA"])[j] * dta_drho)
                dti_den = ni_tr_g[j]
                dti_eff_drho = dti_num / dti_den if abs(dti_den) > 1e-30 else float("nan")
                dti_dr_m = dti_eff_drho / drdrho if abs(drdrho) > 1e-30 else float("nan")
                dte_drmin_tr.append(dte_dr_m)
                dte_drho_tr.append(dte_drho)
                dti_eff_drmin_tr.append(dti_dr_m)
                dti_eff_drho_tr.append(dti_eff_drho)

                qde = results[0]["flux_diff"][j] / dv_drmin_tr[j] if abs(dv_drmin_tr[j]) > 1e-30 else float("nan")
                qdi = (results[1]["flux_diff"][j] + results[2]["flux_diff"][j] + results[3]["flux_diff"][j]) / dv_drmin_tr[j] if abs(dv_drmin_tr[j]) > 1e-30 else float("nan")
                qdiff_e_over_vp.append(qde)
                qdiff_i_over_vp.append(qdi)

                base_e = -e_charge * ne_tr_g[j] * chi_e_tr_g[j] * dte_dr_m
                base_i = -e_charge * ni_tr_g[j] * chi_i_tr_g[j] * dti_dr_m
                geom_e.append((qde * 1e6) / base_e if abs(base_e) > 1e-30 else float("nan"))
                geom_i.append((qdi * 1e6) / base_i if abs(base_i) > 1e-30 else float("nan"))

            # chi_eff: pairs with dT/drho  =>  chi * AR2 * (dr_min/drho)
            chi_e_tr_eff = [chi_e_tr_g[j] * ar2rhog[j] * drmin_drho_tr[j] for j in range(nrmax)]
            chi_i_tr_eff = [chi_i_tr_g[j] * ar2rhog[j] * drmin_drho_tr[j] for j in range(nrmax)]
            # chi_r: pairs with dT/dr_min  =>  chi * AR2 * (dr_min/drho)^2
            # so that q = n * chi_r * (-dT/dr_min) is directly comparable to OMFIT
            chi_e_tr_r = [chi_e_tr_g[j] * ar2rhog[j] * drmin_drho_tr[j]**2 for j in range(nrmax)]
            chi_i_tr_r = [chi_i_tr_g[j] * ar2rhog[j] * drmin_drho_tr[j]**2 for j in range(nrmax)]

            # Flux density from chi_eff in V'-normalized form:
            #   q = n * chi_eff * (-dT/drho) * MW_FACTOR  [MW/m^2]
            # where chi_eff = chi * <|grad rho|^2> * dr_min/drho
            # This equals Qdiff/(dV/dr_min) when cross & conv terms are zero.
            # Note: dT/drho (NOT dT/dr_min) must be used here because
            #   Qdiff/(dV/dr_min) = <|grad rho|^2> * (dr_min/drho) * n * chi * (-dT/drho) * MW
            # and chi_eff already absorbs <|grad rho|^2> * (dr_min/drho), so dT stays in drho.
            qe_from_chi_eff = []
            qi_from_chi_eff = []
            for j in range(nrmax):
                if j < nrmax - 1:
                    dte_drho_j = (te_tr[j + 1] - te_tr[j]) / dr
                    dti_drho_num = (avg_to_grid(d_n["nD"])[j] * (td_tr[j+1]-td_tr[j])/dr
                                  + avg_to_grid(d_n["nT"])[j] * (tt_tr[j+1]-tt_tr[j])/dr
                                  + avg_to_grid(d_n["nA"])[j] * (ta_tr[j+1]-ta_tr[j])/dr)
                else:
                    dte_drho_j = (te_tr[j] - te_tr[j - 1]) / dr
                    dti_drho_num = (avg_to_grid(d_n["nD"])[j] * (td_tr[j]-td_tr[j-1])/dr
                                  + avg_to_grid(d_n["nT"])[j] * (tt_tr[j]-tt_tr[j-1])/dr
                                  + avg_to_grid(d_n["nA"])[j] * (ta_tr[j]-ta_tr[j-1])/dr)
                dti_drho_j = dti_drho_num / ni_tr_g[j] if abs(ni_tr_g[j]) > 1e-30 else 0.0
                qe_from_chi_eff.append(ne_tr_g[j] * chi_e_tr_eff[j] * (-dte_drho_j) * MW_FACTOR)
                qi_from_chi_eff.append(ni_tr_g[j] * chi_i_tr_eff[j] * (-dti_drho_j) * MW_FACTOR)

            fig_det, axes = plt.subplots(2, 5, figsize=(18, 8), sharex='col')
            # Electron row
            axes[0, 0].plot(rho_plot, qdiff_e_over_vp, "s-", lw=2, ms=3, label="TR Qdiff_e/(dV/drmin)")
            axes[0, 0].plot(rho_plot, qe_from_chi_eff, "d-", lw=1.8, ms=3, color="tab:green",
                            label=r"TR $n_e\chi_{e,eff}(-dT_e/d\rho)$")
            axes[0, 0].plot(omfit_rho, np.array(omfit_qe_grad) / 1e6, "--", lw=1.8, color="tab:red", label="OMFIT Qe,eff")
            axes[0, 0].set_title("Electron Flux")
            axes[0, 0].set_ylabel("MW/m$^2$")
            axes[0, 0].grid(True, alpha=0.3)
            axes[0, 0].legend(fontsize=8)

            axes[0, 1].plot(rho_plot, np.array(ne_tr_g) * 10.0, "s-", lw=2, ms=3, label="TR n_e")
            axes[0, 1].plot(omfit_rho, np.array(omfit_ne) / 1e19, "--", lw=1.8, color="tab:red", label="OMFIT n_e")
            axes[0, 1].set_title("Electron Density")
            axes[0, 1].set_ylabel(r"$10^{19} m^{-3}$")
            axes[0, 1].grid(True, alpha=0.3)
            axes[0, 1].legend(fontsize=8)

            omfit_chie_rho = omfit_chie / omfit_drmin_drho  # chi_OMFIT / (dr_min/drho)
            axes[0, 2].plot(rho_plot, chi_e_tr_eff, "s-", lw=2, ms=3,
                            label=r"TR $\chi\langle|\nabla\rho|^2\rangle\,dr/d\rho$")
            axes[0, 2].plot(omfit_rho, omfit_chie_rho, "--", lw=1.8, color="tab:red",
                            label=r"OMFIT $\chi_e\,/\,(dr/d\rho)$")
            axes[0, 2].set_title(r"Electron $\chi$ (in $d\rho$ coord)")
            axes[0, 2].set_ylabel("m$^2$/s $\cdot$ m$^{-1}$")
            axes[0, 2].grid(True, alpha=0.3)
            axes[0, 2].legend(fontsize=7)

            omfit_dte_drho = np.array(omfit_dte_dr) / 1e3 * omfit_drmin_drho  # eV/m -> keV/m -> keV/rho
            axes[0, 3].plot(rho_plot, np.array(dte_drho_tr), "d-", lw=2, ms=3, color="tab:green",
                            label=r"TR $dT_e/d\rho$ [keV]")
            axes[0, 3].plot(omfit_rho, omfit_dte_drho, ":", lw=2, color="tab:orange",
                            label=r"OMFIT $dT_e/d\rho$ [keV]")
            axes[0, 3].plot(rho_plot, np.array(dte_drmin_tr), "s-", lw=2, ms=3, label=r"TR $dT_e/dr_{min}$ [keV/m]")
            axes[0, 3].plot(omfit_rho, np.array(omfit_dte_dr) / 1e3, "--", lw=1.8, color="tab:red",
                            label=r"OMFIT $dT_e/dr_{min}$ [keV/m]")
            axes[0, 3].set_title("Electron Temp Gradient")
            axes[0, 3].set_ylabel("keV / [rho or m]")
            axes[0, 3].grid(True, alpha=0.3)
            axes[0, 3].legend(fontsize=7)

            axes[0, 4].plot(rho_plot, drmin_drho_tr, "s-", lw=2, ms=3, label=r"$dr_{min}/d\rho$ (TR interp)")
            axes[0, 4].plot(omfit_rho, omfit_drmin_drho, "--", lw=1.8, color="tab:red",
                            label=r"$dr_{min}/d\rho$ (OMFIT)")
            axes[0, 4].set_title(r"$dr_{min}/d\rho$ [m]")
            axes[0, 4].set_ylabel("m")
            axes[0, 4].grid(True, alpha=0.3)
            axes[0, 4].legend(fontsize=8)

            # Ion row
            axes[1, 0].plot(rho_plot, qdiff_i_over_vp, "s-", lw=2, ms=3, label="TR Qdiff_i/(dV/drmin)")
            axes[1, 0].plot(rho_plot, qi_from_chi_eff, "d-", lw=1.8, ms=3, color="tab:green",
                            label=r"TR $n_i\chi_{i,eff}(-dT_i/d\rho)$")
            axes[1, 0].plot(omfit_rho, np.array(omfit_qi_grad) / 1e6, "--", lw=1.8, color="tab:red", label="OMFIT Qi,eff")
            axes[1, 0].set_title("Ion Flux")
            axes[1, 0].set_ylabel("MW/m$^2$")
            axes[1, 0].set_xlabel("rho")
            axes[1, 0].grid(True, alpha=0.3)
            axes[1, 0].legend(fontsize=8)

            axes[1, 1].plot(rho_plot, np.array(ni_tr_g) * 10.0, "s-", lw=2, ms=3, label="TR n_i")
            axes[1, 1].plot(omfit_rho, np.array(omfit_ni) / 1e19, "--", lw=1.8, color="tab:red", label="OMFIT n_i")
            axes[1, 1].set_title("Ion Density")
            axes[1, 1].set_xlabel("rho")
            axes[1, 1].set_ylabel(r"$10^{19} m^{-3}$")
            axes[1, 1].grid(True, alpha=0.3)
            axes[1, 1].legend(fontsize=8)

            omfit_chii_rho = omfit_chii / omfit_drmin_drho
            axes[1, 2].plot(rho_plot, chi_i_tr_eff, "s-", lw=2, ms=3,
                            label=r"TR $\chi\langle|\nabla\rho|^2\rangle\,dr/d\rho$")
            axes[1, 2].plot(omfit_rho, omfit_chii_rho, "--", lw=1.8, color="tab:red",
                            label=r"OMFIT $\chi_i\,/\,(dr/d\rho)$")
            axes[1, 2].set_title(r"Ion $\chi$ (in $d\rho$ coord)")
            axes[1, 2].set_xlabel("rho")
            axes[1, 2].set_ylabel("m$^2$/s $\cdot$ m$^{-1}$")
            axes[1, 2].grid(True, alpha=0.3)
            axes[1, 2].legend(fontsize=7)

            omfit_dti_drho = np.array(omfit_dti_dr) / 1e3 * omfit_drmin_drho
            axes[1, 3].plot(rho_plot, np.array(dti_eff_drho_tr), "d-", lw=2, ms=3, color="tab:green",
                            label=r"TR $dT_i/d\rho$ [keV]")
            axes[1, 3].plot(omfit_rho, omfit_dti_drho, ":", lw=2, color="tab:orange",
                            label=r"OMFIT $dT_i/d\rho$ [keV]")
            axes[1, 3].plot(rho_plot, np.array(dti_eff_drmin_tr), "s-", lw=2, ms=3,
                            label=r"TR $dT_i/dr_{min}$ [keV/m]")
            axes[1, 3].plot(omfit_rho, np.array(omfit_dti_dr) / 1e3, "--", lw=1.8, color="tab:red",
                            label=r"OMFIT $dT_i/dr_{min}$ [keV/m]")
            axes[1, 3].set_title("Ion Temp Gradient")
            axes[1, 3].set_xlabel("rho")
            axes[1, 3].set_ylabel("keV / [rho or m]")
            axes[1, 3].grid(True, alpha=0.3)
            axes[1, 3].legend(fontsize=7)

            axes[1, 4].plot(rho_plot, ar2rhog, "s-", lw=2, ms=3, label=r"$\langle|\nabla\rho|^2\rangle$")
            ar2_times_drdrho = [ar2rhog[j] * drmin_drho_tr[j] for j in range(nrmax)]
            axes[1, 4].plot(rho_plot, ar2_times_drdrho, "d-", lw=2, ms=3, color="tab:green",
                            label=r"$\langle|\nabla\rho|^2\rangle \cdot dr_{min}/d\rho$")
            axes[1, 4].set_title("Geometry Factors")
            axes[1, 4].set_xlabel("rho")
            axes[1, 4].grid(True, alpha=0.3)
            axes[1, 4].legend(fontsize=7)

            plot_path_det = os.path.join(BASE, "qeff_vs_qdiff_detailed_plot.png")
            plt.tight_layout()
            plt.savefig(plot_path_det, dpi=160)
            plt.close(fig_det)
            print(f"Detailed Qeff/Qdiff Plot PNG: {plot_path_det}")

            # --- Ratio / compensation plot ---
            eps = 1e-30
            omfit_chie_tr = np.interp(rho_plot, omfit_rho, omfit_chie)
            omfit_chii_tr = np.interp(rho_plot, omfit_rho, omfit_chii)
            omfit_dte_tr = np.interp(rho_plot, omfit_rho, np.array(omfit_dte_dr) / 1e3)
            omfit_dti_tr = np.interp(rho_plot, omfit_rho, np.array(omfit_dti_dr) / 1e3)
            omfit_ne_tr = np.interp(rho_plot, omfit_rho, omfit_ne)
            omfit_ni_tr = np.interp(rho_plot, omfit_rho, omfit_ni)

            chi_ratio_e = np.array(chi_e_tr_eff) / (omfit_chie_tr + eps)
            chi_ratio_i = np.array(chi_i_tr_eff) / (omfit_chii_tr + eps)
            grad_ratio_e = np.abs(np.array(dte_drmin_tr)) / (np.abs(omfit_dte_tr) + eps)
            grad_ratio_i = np.abs(np.array(dti_eff_drmin_tr)) / (np.abs(omfit_dti_tr) + eps)
            dens_ratio_e = np.array(ne_tr_g) / (omfit_ne_tr + eps)
            dens_ratio_i = np.array(ni_tr_g) / (omfit_ni_tr + eps)
            prod_ratio_e = chi_ratio_e * grad_ratio_e * dens_ratio_e
            prod_ratio_i = chi_ratio_i * grad_ratio_i * dens_ratio_i

            fig_rat, (ax_re, ax_ri) = plt.subplots(2, 1, figsize=(8, 9), sharex=True)
            ax_re.plot(rho_plot, chi_ratio_e, "-", lw=2, label=r"$\chi_{TR,eff}/\chi_{OMFIT}$")
            ax_re.plot(rho_plot, grad_ratio_e, "--", lw=2, label=r"$|dT/dr|_{TR}/|dT/dr|_{OMFIT}$")
            ax_re.plot(rho_plot, dens_ratio_e, ":", lw=2, label=r"$n_{TR}/n_{OMFIT}$")
            ax_re.plot(rho_plot, prod_ratio_e, "-.", lw=2, label="product ratio")
            ax_re.axhline(1.0, color="k", lw=1, alpha=0.5)
            ax_re.set_ylabel("Ratio")
            ax_re.set_title("Electron Compensation: why flux can still match")
            ax_re.grid(True, alpha=0.3)
            ax_re.legend(fontsize=8)

            ax_ri.plot(rho_plot, chi_ratio_i, "-", lw=2, label=r"$\chi_{TR,eff}/\chi_{OMFIT}$")
            ax_ri.plot(rho_plot, grad_ratio_i, "--", lw=2, label=r"$|dT/dr|_{TR}/|dT/dr|_{OMFIT}$")
            ax_ri.plot(rho_plot, dens_ratio_i, ":", lw=2, label=r"$n_{TR}/n_{OMFIT}$")
            ax_ri.plot(rho_plot, prod_ratio_i, "-.", lw=2, label="product ratio")
            ax_ri.axhline(1.0, color="k", lw=1, alpha=0.5)
            ax_ri.set_xlabel("rho")
            ax_ri.set_ylabel("Ratio")
            ax_ri.set_title("Ion Compensation: why flux can still match")
            ax_ri.grid(True, alpha=0.3)
            ax_ri.legend(fontsize=8)

            plot_path_rat = os.path.join(BASE, "qeff_vs_qdiff_compensation_plot.png")
            plt.tight_layout()
            plt.savefig(plot_path_rat, dpi=160)
            plt.close(fig_rat)
            print(f"Compensation Plot PNG: {plot_path_rat}")

        # --- Electron / ion separated plots ---
        rho_plot = []
        qe_src = []
        qe_recon = []
        qe_csv = []
        qi_src = []
        qi_recon = []
        qi_csv = []
        for j in range(nrmax):
            rg = rho_m[j] + dr / 2 if j < nrmax - 1 else 1.0
            rho_plot.append(rg)

            # electron = species 0
            qe_src.append(results[0]["cum_src"][j])
            qe_recon.append(results[0]["flux_total"][j])

            # ion total = D + T + He4
            qi_src.append(results[1]["cum_src"][j] + results[2]["cum_src"][j] + results[3]["cum_src"][j])
            qi_recon.append(results[1]["flux_total"][j] + results[2]["flux_total"][j] + results[3]["flux_total"][j])

            qcsv_e = 0.0
            qcsv_i = 0.0
            if d_qtr is not None:
                qtr_cols = [c for c in d_qtr if c != "X"]
                if len(qtr_cols) >= 4:
                    qcsv_e = d_qtr[qtr_cols[0]][j]
                    qcsv_i = d_qtr[qtr_cols[1]][j] + d_qtr[qtr_cols[2]][j] + d_qtr[qtr_cols[3]][j]
            qe_csv.append(qcsv_e)
            qi_csv.append(qcsv_i)

        fig2, (ax_e, ax_i) = plt.subplots(2, 1, figsize=(8, 10), sharex=True)

        ax_e.plot(rho_plot, qe_src, "o-", lw=2, ms=4, label="TR Q_source,e")
        ax_e.plot(rho_plot, qe_recon, "s-", lw=2, ms=4, label="TR QTR_e reconstructed")
        if d_qtr is not None:
            ax_e.plot(rho_plot, qe_csv, "^-", lw=2, ms=4, label="TR QTR_e from CSV")
        if omfit_rho is not None and omfit_qe_power_mw is not None:
            ax_e.plot(omfit_rho, omfit_qe_power_mw, "--", lw=1.8, color="tab:red",
                      label="OMFIT Qe,int -> power")
        ax_e.set_ylabel("Power [MW]")
        ax_e.set_title("Electron Source vs Transport")
        ax_e.grid(True, alpha=0.3)
        ax_e.legend()

        ax_i.plot(rho_plot, qi_src, "o-", lw=2, ms=4, label="TR Q_source,i (D+T+He4)")
        ax_i.plot(rho_plot, qi_recon, "s-", lw=2, ms=4, label="TR QTR_i reconstructed")
        if d_qtr is not None:
            ax_i.plot(rho_plot, qi_csv, "^-", lw=2, ms=4, label="TR QTR_i from CSV")
        if omfit_rho is not None and omfit_qi_power_mw is not None:
            ax_i.plot(omfit_rho, omfit_qi_power_mw, "--", lw=1.8, color="tab:red",
                      label="OMFIT Qi,int -> power")
        ax_i.set_xlabel("rho")
        ax_i.set_ylabel("Power [MW]")
        ax_i.set_title("Ion Source vs Transport (D+T+He4)")
        ax_i.grid(True, alpha=0.3)
        ax_i.legend()

        plot_path2 = os.path.join(BASE, "qsource_vs_qtr_species_plot.png")
        plt.tight_layout()
        plt.savefig(plot_path2, dpi=160)
        plt.close(fig2)
        print(f"Species Plot PNG: {plot_path2}")
    except Exception as e:
        print(f"Plot generation skipped: {e}")


if __name__ == "__main__":
    main()
