#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Compare JTOT and JBS from TR CSV output with recalculation (Version 2).

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

Author: Auto-generated for JTOT/JBS analysis (v2 with geometry)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
import os

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

# ============================================================
# Helper functions
# ============================================================

def find_csv_by_title(pattern, data_dir='.'):
    """Find CSV file whose first line contains the given pattern."""
    files = sorted(glob.glob(os.path.join(data_dir, 'tr_data_*.csv')))
    for f in files:
        with open(f, 'r') as fh:
            title = fh.readline().strip()
            if pattern in title:
                return f
    return None

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

def trapped_fraction_kim(eps):
    """Trapped particle fraction: Kim et al., PoF B 3 2050 (1991) eq(C18)."""
    return 1.46 * np.sqrt(eps) - 0.46 * eps**1.5

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
    print("JTOT and JBS COMPARISON V2: Using Exact Geometry from TR")
    print("=" * 70)

    # ============================================================
    # 1. Read CSV data
    # ============================================================
    print("\n1. Loading CSV data...")

    # Current densities
    f_j = find_csv_by_title('JTOT,JOH,JNB,JRF,JBS', data_dir)
    if f_j is None:
        print("ERROR: Cannot find current density CSV")
        return
    r_j, j_cols, _ = read_tr_csv(f_j)
    print(f"   Current densities: {os.path.basename(f_j)}")

    # Densities
    f_n = find_csv_by_title('n(NS)', data_dir)
    r_n, n_cols, _ = read_tr_csv(f_n)
    print(f"   Densities: {os.path.basename(f_n)}")

    # Temperatures
    f_t = find_csv_by_title('T(NS)', data_dir)
    r_t, t_cols, _ = read_tr_csv(f_t)
    print(f"   Temperatures: {os.path.basename(f_t)}")

    # Safety factor
    f_q = find_csv_by_title('QP  vs r', data_dir)
    r_q, q_cols, _ = read_tr_csv(f_q)
    print(f"   Safety factor: {os.path.basename(f_q)}")

    # Zeff
    f_z = find_csv_by_title('ZEFF  vs r', data_dir)
    r_z, z_cols, _ = read_tr_csv(f_z)
    print(f"   Zeff: {os.path.basename(f_z)}")

    # Geometry quantities (new in V2)
    f_geom1 = find_csv_by_title('TTRHOG,ABVRHOG', data_dir)
    f_geom2 = find_csv_by_title('RDP,RDPVRHOG', data_dir)
    f_geom3 = find_csv_by_title('EPSRHO,RJCB', data_dir)
    f_geom4 = find_csv_by_title('RMJRHO,RMNRHO', data_dir)

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

    # Pressures [10^20/m^3 * keV]
    Pe = ne * Te
    Pi = nD * TD + nT * TT + nA * TA
    Ti = (nD * TD + nT * TT + nA * TA) / (nD + nT + nA + 1e-20)
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

    # Boundary point (i = NRMAX-1) - TR uses 3-point extrapolation, we use backward difference
    i = NRMAX - 1
    if has_geometry and abs(RDP[i]) > 1e-20:
        EPS_i = eps[i]
        EPSS = np.sqrt(EPS_i)**3
        ft[i] = 1.46 * np.sqrt(EPS_i) - 0.46 * EPS_i**1.5

        # Extrapolate Zeff: ZEFFL=2.D0*ZEFF(NR-1)-ZEFF(NR-2)
        Zeff_avg = 2.0 * Zeff[i-1] - Zeff[i-2]

        Pe_avg = Pe[i]
        Pi_avg = Pi[i]
        Te_avg = Te[i]
        Ti_avg = Ti[i]
        ne_avg = ne[i]
        ANI_avg = ANI_arr[i]
        RPe_avg = Pe_avg / (Pe_avg + Pi_avg + 1e-20)

        # Backward difference for gradients
        DPe_i = (Pe[i] - Pe[i-1]) / DR
        DPi_i = (Pi[i] - Pi[i-1]) / DR
        DTe_i = (Te[i] - Te[i-1]) / DR
        DTi_i = (Ti[i] - Ti[i-1]) / DR

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
    print(f"\n   {'r/a':>6} {'JBS_csv':>12} {'JBS_calc':>12} {'Ratio':>10} {'Diff%':>10}")
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
    print(f"\n6. Diagnostic: Key quantities at r/a=0.5:")
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

    # Detailed diagnostic at r/a=0.25 (index 12)
    print(f"\n   Detailed diagnostic at r/a=0.25 (index 12):")
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

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.suptitle('JTOT and JBS Analysis V2: Using Exact Geometry', fontsize=14, fontweight='bold')

    # Plot 1: Current components
    ax = axes[0, 0]
    ax.plot(r, JTOT_csv, 'k-', lw=2, label='JTOT')
    ax.plot(r, JOH_csv, 'b-', lw=1.5, label='JOH')
    ax.plot(r, JNB_csv, 'g-', lw=1.5, label='JNB')
    ax.plot(r, JRF_csv, 'c-', lw=1.5, label='JRF')
    ax.plot(r, JBS_csv, 'r-', lw=1.5, label='JBS')
    ax.axhline(y=0, color='gray', ls='--', alpha=0.5)
    ax.set_xlabel('r/a')
    ax.set_ylabel('J [MA/m²]')
    ax.set_title('Current Density Components (from CSV)')
    ax.legend(loc='best', fontsize=9)
    ax.grid(True, alpha=0.3)

    # Plot 2: JBS comparison
    ax = axes[0, 1]
    ax.plot(r, JBS_csv, 'r-', lw=2, label='JBS (CSV/TR)')
    ax.plot(r, JBS_avg, 'b--', lw=1.5, label='JBS (recalculated)')
    ax.axhline(y=0, color='gray', ls='--', alpha=0.5)
    ax.set_xlabel('r/a')
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
    ax.set_xlabel('r/a')
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
    ax.set_xlabel('r/a')
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
        ax.set_xlabel('r/a')
        ax.set_ylabel('TTRHOG, EPSRHO')
        ax2.set_ylabel('RDP [Wb]')
        ax.set_title('Geometry Quantities from TR')
        lines = [l1, l2, l3]
        ax.legend(lines, [l.get_label() for l in lines], loc='best')
    else:
        ax.plot(r, eps, 'r-', lw=1.5, label='eps (approx)')
        ax.set_xlabel('r/a')
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
    ax.set_xlabel('r/a')
    ax.set_ylabel('Collisionality (ν*)')
    ax2.set_ylabel('Trapped fraction')
    ax.set_title('Collisionality and Trapped Fraction')
    ax.set_yscale('log')
    lines = [l1, l2, l3]
    ax.legend(lines, [l.get_label() for l in lines], loc='best')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('jtot_jbs_comparison_v2.png', dpi=150, bbox_inches='tight')
    print("   Plot saved to jtot_jbs_comparison_v2.png")

    # Save data
    df = pd.DataFrame({
        'r_a': r,
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
    df.to_csv('jtot_jbs_comparison_v2.csv', index=False)
    print("   Data saved to jtot_jbs_comparison_v2.csv")

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
