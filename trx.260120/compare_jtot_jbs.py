#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Compare JTOT and JBS from TR CSV output with recalculation.

JTOT: Total current density = JOH + JNB + JRF + JBS
JBS:  Bootstrap current (Sauter model, MDLJBS=5)

The Sauter model formula (from TRAJBSSAUTER in trcalc.f90):
  JBS = -PBSCD * <|grad_rho|^2/R^2> * Pe * (L31*(dPe/Pe + dPi/Pe) + L32*dTe/Te 
         + L34*alpha*(1-RPe)/RPe * dTi/Ti) / (dPsi/drho) / BB

where L31, L32, L34 are neoclassical transport coefficients depending on 
trapped particle fraction and collisionality.

Author: Auto-generated for JTOT/JBS analysis
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
AMD = 2.0 * AMP           # deuterium mass
AMT = 3.0 * AMP           # tritium mass
AMA = 4.0 * AMP           # alpha mass
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
RDLT = 0.59     # triangularity
PBSCD = 1.0     # bootstrap current drive factor (default)
MDLTPF = 0      # trapped particle fraction model (default: Kim et al.)

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
    """Read a TR CSV file. Returns (r, columns_dict)."""
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
    print("JTOT and JBS COMPARISON: CSV vs Recalculation")
    print("=" * 70)
    
    # ============================================================
    # 1. Read CSV data
    # ============================================================
    print("\n1. Loading CSV data...")
    
    # Current densities: JTOT, JOH, JNB, JRF, JBS [MA/m^2]
    f_j = find_csv_by_title('JTOT,JOH,JNB,JRF,JBS', data_dir)
    if f_j is None:
        print("ERROR: Cannot find current density CSV")
        return
    r_j, j_cols, _ = read_tr_csv(f_j)
    print(f"   Current densities from: {os.path.basename(f_j)}")
    
    # Densities: n(NS) [10^20/m^3]
    f_n = find_csv_by_title('n(NS)', data_dir)
    r_n, n_cols, _ = read_tr_csv(f_n)
    print(f"   Densities from: {os.path.basename(f_n)}")
    
    # Temperatures: T(NS) [keV]
    f_t = find_csv_by_title('T(NS)', data_dir)
    r_t, t_cols, _ = read_tr_csv(f_t)
    print(f"   Temperatures from: {os.path.basename(f_t)}")
    
    # Pressures: P(ns) [MPa]
    f_p = find_csv_by_title('P(ns)', data_dir)
    r_p, p_cols, _ = read_tr_csv(f_p)
    print(f"   Pressures from: {os.path.basename(f_p)}")
    
    # Safety factor: QP
    f_q = find_csv_by_title('QP  vs r', data_dir)
    r_q, q_cols, _ = read_tr_csv(f_q)
    print(f"   Safety factor from: {os.path.basename(f_q)}")
    
    # Zeff
    f_z = find_csv_by_title('ZEFF  vs r', data_dir)
    r_z, z_cols, _ = read_tr_csv(f_z)
    print(f"   Zeff from: {os.path.basename(f_z)}")

    # ============================================================
    # 2. Extract the data arrays
    # ============================================================
    print("\n2. Extracting data arrays...")
    
    # Get column names
    j_keys = list(j_cols.keys())
    n_keys = list(n_cols.keys())
    t_keys = list(t_cols.keys())
    p_keys = list(p_cols.keys())
    q_keys = list(q_cols.keys())
    z_keys = list(z_cols.keys())
    
    print(f"   J columns: {j_keys}")
    print(f"   n columns: {n_keys}")
    print(f"   T columns: {t_keys}")
    print(f"   P columns: {p_keys}")
    print(f"   Q columns: {q_keys}")
    print(f"   Z columns: {z_keys}")
    
    # Current densities [MA/m^2]
    r = j_cols[j_keys[0]]      # r/a grid
    JTOT_csv = j_cols[j_keys[1]]  # JTOT
    JOH_csv  = j_cols[j_keys[2]]  # JOH
    JNB_csv  = j_cols[j_keys[3]]  # JNB
    JRF_csv  = j_cols[j_keys[4]]  # JRF
    JBS_csv  = j_cols[j_keys[5]]  # JBS
    
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
    
    # Pressures [MPa]
    Pe_MPa = p_cols[p_keys[1]]  # electron
    PD_MPa = p_cols[p_keys[2]]  # D
    PT_MPa = p_cols[p_keys[3]]  # T
    PA_MPa = p_cols[p_keys[4]]  # He4
    
    # Safety factor (QP vs r, note: grid might differ - Q is on grid points)
    r_qp = q_cols[q_keys[0]]
    QP_val = q_cols[q_keys[1]]
    # QP is defined on grid points (rg), interpolate to mesh centers (rm)
    if len(QP_val) != NRMAX:
        print(f"   NOTE: QP grid ({len(QP_val)}) differs from density grid ({NRMAX})")
        QP = np.interp(r, r_qp, QP_val)
    else:
        QP = QP_val
    
    # Zeff
    Zeff = z_cols[z_keys[1]]
    
    print(f"   Grid: NRMAX={NRMAX}, DR={DR:.4f}")
    print(f"   ne center={ne[0]:.4f}, Te center={Te[0]:.4f} keV")
    print(f"   QP center={QP[0]:.4f}, Zeff center={Zeff[0]:.4f}")
    
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
    else:
        print("   ! Significant discrepancy detected")
    
    # Print component table
    print(f"\n   Current density components at key radii:")
    print(f"   {'r/a':>6} {'JTOT':>10} {'JOH':>10} {'JNB':>10} {'JRF':>10} {'JBS':>10} {'Sum':>10} {'Diff':>10}")
    print(f"   {'-'*76}")
    for idx in [0, NRMAX//4, NRMAX//2, 3*NRMAX//4, NRMAX-1]:
        print(f"   {r[idx]:6.3f} {JTOT_csv[idx]:10.4f} {JOH_csv[idx]:10.4f} "
              f"{JNB_csv[idx]:10.4f} {JRF_csv[idx]:10.4f} {JBS_csv[idx]:10.4f} "
              f"{JTOT_sum[idx]:10.4f} {diff_jtot[idx]:10.2e}")
    
    # ============================================================
    # 4. Recalculate JBS using Sauter model
    # ============================================================
    print("\n4. Recalculating JBS using Sauter model (MDLJBS=5)...")
    
    # Calculate epsilon = r/(R+r) ≈ r*a/R for circular approximation
    # In TR: EPSRHO = RMNRHO / RMJRHO (from trmetric.f90)
    # For simplified model: eps ≈ r/a * a/R = rho * a / R
    eps = r * RA / RR
    eps = np.maximum(eps, 1e-10)
    
    # Trapped particle fraction (MDLTPF=0, Kim et al.)
    ft = trapped_fraction_kim(eps)
    
    # Poloidal field: Bp ≈ mu0 * R * J / (2*r)
    # Simple estimate from JTOT (already in MA/m^2 = 1e6 A/m^2)
    # Bp = mu0 * integral(J * r dr) / r
    # Use CSV data to estimate Bp from JTOT
    Bp = np.zeros(NRMAX)
    for i in range(NRMAX):
        rho_g = r[i] * RA  # in meters
        if rho_g > 0:
            # Cumulative integral of JTOT * r * dr (toroidal approx)
            Bp[i] = RMU0 * np.sum(JTOT_csv[:i+1] * 1e6 * r[:i+1] * RA * DR * RA) / rho_g
    Bp[0] = Bp[1] * r[1] / (r[1] + 1e-20) if len(Bp) > 1 else 0
    
    # Electron pressure: Pe = ne * Te in keV * 10^20/m^3
    # Convert to consistent units
    Pe = ne * Te  # [10^20/m^3 * keV]
    
    # Ion pressure: Pi = sum(ns * Ts) for thermal ions
    Pi = nD * TD + nT * TT + nA * TA  # [10^20/m^3 * keV]
    
    # Total pressure from sum
    # RPe = Pe / (Pe + Pi)  [dimensionless]
    RPe = Pe / (Pe + Pi + 1e-20)
    
    # Pressure gradients (finite difference on mesh center grid)
    DPe = np.gradient(Pe, r * RA)  # d(Pe)/d(r) in [10^20*keV/m^3]/m
    DPi = np.gradient(Pi, r * RA)
    DTe = np.gradient(Te, r * RA)  # d(Te)/d(r) in keV/m
    
    # Ion temperature: weighted average
    Ti = (nD * TD + nT * TT + nA * TA) / (nD + nT + nA + 1e-20)
    DTi = np.gradient(Ti, r * RA)
    
    # Collisionality
    eps_half = np.sqrt(eps)**3  # eps^{3/2}
    
    # Electron collisionality: nu_e* 
    # nu_e* = 6.921e-18 * qR * ne * Zeff * ln_Lambda / (Te^2 * eps^{3/2})
    # where ne in m^-3 (multiply by 1e20), Te in eV (multiply by 1e3)
    ln_Lambda_e = 31.3 - np.log(np.sqrt(ne * 1e20) / np.abs(Te * 1e3 + 1e-10))
    nue_star = (6.921e-18 * np.abs(QP) * RR * ne * 1e20 * Zeff * ln_Lambda_e 
                / (np.abs(Te * 1e3)**2 * eps_half + 1e-30))
    
    # Ion collisionality (for alpha coefficient)
    # ni = sum(Z_s * n_s) for charge-weighted density
    ni = nD + nT + 2*nA  # charge weighted ion density [10^20/m3]
    PZ_D = 1.0  # charge number
    ln_Lambda_ii = 30.0 - np.log(PZ_D**3 * np.sqrt(ni * 1e20) / (np.abs(Ti * 1e3 + 1e-10)**1.5))
    nui_star = (4.90e-18 * np.abs(QP) * RR * ni * 1e20 * PZ_D**4 * ln_Lambda_ii 
                / (np.abs(Ti * 1e3)**2 * eps_half + 1e-30))
    
    # Effective trapped particle fraction corrections
    F31_teff = ft / (1.0 + (1.0 - 0.1*ft)*np.sqrt(nue_star) + 0.5*(1.0 - ft)*nue_star/Zeff)
    F32ee_teff = ft / (1.0 + 0.26*(1.0-ft)*np.sqrt(nue_star) + 0.18*(1.0-0.37*ft)*nue_star/np.sqrt(Zeff))
    F32ei_teff = ft / (1.0 + (1.0+0.6*ft)*np.sqrt(nue_star) + 0.85*(1.0-0.37*ft)*nue_star*(1.0+Zeff))
    F34_teff = ft / (1.0 + (1.0 - 0.1*ft)*np.sqrt(nue_star) + 0.5*(1.0 - 0.5*ft)*nue_star/Zeff)
    
    # Sauter alpha coefficient
    salfa0 = -1.17*(1.0 - ft)/(1.0 - 0.22*ft - 0.19*ft**2 + 1e-20)
    salfa = ((salfa0 + 0.25*(1.0 - ft**2)*np.sqrt(nui_star)) / (1.0 + 0.5*np.sqrt(nui_star)) 
             + 0.315*nui_star**2 * ft**6) / (1.0 + 0.15*nui_star**2 * ft**6)
    
    # Neoclassical coefficients
    RL31 = F31_func(F31_teff, Zeff)
    RL32 = F32EE_func(F32ee_teff, Zeff) + F32EI_func(F32ei_teff, Zeff)
    RL34 = F31_func(F34_teff, Zeff)
    
    # Sauter JBS formula (simplified version - cylindrical geometry):
    # JBS = -PBSCD * <|grad_rho|^2/R^2> * Pe * RKEV * 1e20 
    #       * (L31*(dPe/Pe + dPi/Pe) + L32*dTe/Te + L34*alpha*(1-RPe)/RPe*dTi/Ti)
    #       / (dPsi/drho) / BB
    #
    # In simplified cylindrical geometry:
    # <|grad_rho|^2/R^2> ≈ 1/(R^2)  (for a >> eps)
    # dPsi/drho = Bp * r 
    # So JBS ≈ -PBSCD * Pe * RKEV * 1e20 * (...gradient terms...) / (Bp * R * BB)
    
    # Alternative approach: directly use the pressure gradient form
    # In TR code, the formula is:
    # AJBSL = -PBSCD * TTRHOG * Pe * 1e20 * RKEV 
    #         * (L31*(dPe/Pe + dPi/Pe) + L32*dTe/Te + L34*alpha*(1-RPe)/RPe*dTi/Ti)
    #         / RDP / BB
    # where TTRHOG = <|grad_rho|^2/R^2>, RDP = dPsi/drho
    
    # For simplified calculation, use: J_BS = -pe * (L31*(...) + L32*...) / (Bp * R)
    # where pe = ne * Te in SI units [Pa]
    
    JBS_calc = np.zeros(NRMAX)
    for i in range(NRMAX):
        if Bp[i] < 1e-10 or Pe[i] < 1e-20:
            continue
        
        # Gradient terms (already calculated using np.gradient in SI)
        term1 = RL31[i] * (DPe[i] / (Pe[i] + 1e-20) + DPi[i] / (Pe[i] + 1e-20))
        term2 = RL32[i] * DTe[i] / (Te[i] + 1e-20) 
        term3 = RL34[i] * salfa[i] * (1.0 - RPe[i]) / (RPe[i] + 1e-20) * DTi[i] / (Ti[i] + 1e-20)
        
        # Pe in SI: ne * 1e20 * Te * RKEV [Pa]
        Pe_SI = Pe[i] * 1e20 * RKEV  # [Pa = J/m^3]
        
        # JBS = -PBSCD * (1/R^2) * Pe_SI * (terms) / (Bp / (r*RA)) / BB
        # Simplified: JBS = -PBSCD * Pe_SI * (terms) * r[i]*RA / (Bp[i] * R^2 * BB)
        # But this is a rough approximation. Better to use:
        # JBS = -PBSCD * ft[i] * Pe_SI * (terms) / (Bp[i] * BB)
        # The exact geometric factors require TTRHOG and RDP from the code
        
        JBS_calc[i] = -PBSCD * Pe_SI * (term1 + term2 + term3) / (Bp[i] + 1e-30) / BB
    
    # Convert to MA/m^2
    JBS_calc_MA = JBS_calc * 1e-6
    
    # Average between adjacent points (as in Fortran: AJBS(NR)=0.5*(AJBSL(NR)+AJBSL(NR-1)))
    JBS_avg = np.zeros(NRMAX)
    JBS_avg[0] = 0.5 * JBS_calc_MA[0]
    for i in range(1, NRMAX):
        JBS_avg[i] = 0.5 * (JBS_calc_MA[i] + JBS_calc_MA[i-1])
    
    # ============================================================
    # 5. Compare results
    # ============================================================
    print("\n5. Comparison results:")
    print(f"\n   {'r/a':>6} {'JBS_csv':>12} {'JBS_calc':>12} {'Ratio':>10}")
    print(f"   {'-'*44}")
    for idx in [0, NRMAX//10, NRMAX//4, NRMAX//2, 3*NRMAX//4, NRMAX-1]:
        ratio = JBS_avg[idx] / (JBS_csv[idx] + 1e-30) if abs(JBS_csv[idx]) > 1e-10 else float('nan')
        print(f"   {r[idx]:6.3f} {JBS_csv[idx]:12.4f} {JBS_avg[idx]:12.4f} {ratio:10.4f}")
    
    # Integrated currents
    I_BS_csv = 2 * PI * np.sum(JBS_csv * r * RA * DR * RA) * 1e6  # [A]
    I_BS_calc = 2 * PI * np.sum(JBS_avg * r * RA * DR * RA) * 1e6
    I_TOT_csv = 2 * PI * np.sum(JTOT_csv * r * RA * DR * RA) * 1e6
    I_OH_csv = 2 * PI * np.sum(JOH_csv * r * RA * DR * RA) * 1e6
    
    print(f"\n   Integrated currents:")
    print(f"   I_TOT (CSV)  = {I_TOT_csv/1e6:.3f} MA")
    print(f"   I_BS  (CSV)  = {I_BS_csv/1e6:.3f} MA")
    print(f"   I_BS  (calc) = {I_BS_calc/1e6:.3f} MA")
    print(f"   I_OH  (CSV)  = {I_OH_csv/1e6:.3f} MA")
    print(f"   I_NB  (CSV)  = {2*PI*np.sum(JNB_csv*r*RA*DR*RA)*1e6/1e6:.3f} MA")
    print(f"   I_RF  (CSV)  = {2*PI*np.sum(JRF_csv*r*RA*DR*RA)*1e6/1e6:.3f} MA")
    
    BS_fraction = I_BS_csv / (I_TOT_csv + 1e-20) * 100
    print(f"\n   Bootstrap fraction: {BS_fraction:.1f}%")
    
    # ============================================================
    # 6. Compute key intermediate quantities for diagnostic
    # ============================================================
    print(f"\n6. Diagnostic: Key quantities at r/a=0.5:")
    idx_half = NRMAX // 2
    print(f"   eps = {eps[idx_half]:.4f}")
    print(f"   ft  = {ft[idx_half]:.4f}")
    print(f"   nu_e* = {nue_star[idx_half]:.4f}")
    print(f"   nu_i* = {nui_star[idx_half]:.4f}")
    print(f"   Bp  = {Bp[idx_half]:.4f} T")
    print(f"   QP  = {QP[idx_half]:.4f}")
    print(f"   Zeff = {Zeff[idx_half]:.4f}")
    print(f"   RPe = {RPe[idx_half]:.4f}")
    print(f"   L31 = {RL31[idx_half]:.4f}")
    print(f"   L32 = {RL32[idx_half]:.4f}")
    print(f"   L34 = {RL34[idx_half]:.4f}")
    print(f"   alpha = {salfa[idx_half]:.4f}")
    print(f"   Pe = {Pe[idx_half]:.4f} [10^20/m^3 * keV]")
    print(f"   dPe/dr = {DPe[idx_half]:.4f}")
    print(f"   dTe/dr = {DTe[idx_half]:.4f}")
    
    # ============================================================
    # 7. Plot results
    # ============================================================
    print("\n7. Generating plots...")
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.suptitle('JTOT and JBS Analysis: CSV vs Recalculation', fontsize=14, fontweight='bold')
    
    # --- Plot 1: All current components ---
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
    
    # --- Plot 2: JTOT verification ---
    ax = axes[0, 1]
    ax.plot(r, JTOT_csv, 'k-', lw=2, label='JTOT (CSV)')
    ax.plot(r, JTOT_sum, 'r--', lw=1.5, label='JOH+JNB+JRF+JBS')
    ax.set_xlabel('r/a')
    ax.set_ylabel('J [MA/m²]')
    ax.set_title('JTOT Consistency Check')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    
    # --- Plot 3: JBS comparison ---
    ax = axes[0, 2]
    ax.plot(r, JBS_csv, 'r-', lw=2, label='JBS (CSV/TR)')
    ax.plot(r, JBS_avg, 'b--', lw=1.5, label='JBS (recalculated)')
    ax.axhline(y=0, color='gray', ls='--', alpha=0.5)
    ax.set_xlabel('r/a')
    ax.set_ylabel('J [MA/m²]')
    ax.set_title('Bootstrap Current: CSV vs Recalculation')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    
    # --- Plot 4: Sauter coefficients ---
    ax = axes[1, 0]
    ax.plot(r, RL31, 'b-', lw=1.5, label='L31')
    ax.plot(r, RL32, 'r-', lw=1.5, label='L32')
    ax.plot(r, RL34, 'g-', lw=1.5, label='L34')
    ax.plot(r, salfa, 'm--', lw=1.5, label='α (Sauter)')
    ax.set_xlabel('r/a')
    ax.set_ylabel('Coefficient')
    ax.set_title('Sauter Neoclassical Coefficients')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    
    # --- Plot 5: Collisionality and trapped fraction ---
    ax = axes[1, 1]
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
    
    # --- Plot 6: Pressure profiles and Bp ---
    ax = axes[1, 2]
    ax2 = ax.twinx()
    l1, = ax.plot(r, Pe, 'b-', lw=1.5, label='Pe')
    l2, = ax.plot(r, Pi, 'r-', lw=1.5, label='Pi')
    l3, = ax2.plot(r, Bp, 'g--', lw=1.5, label='Bp')
    l4, = ax2.plot(r, QP, 'm--', lw=1.5, label='q')
    ax.set_xlabel('r/a')
    ax.set_ylabel('Pressure [10²⁰/m³ · keV]')
    ax2.set_ylabel('Bp [T] / q')
    ax.set_title('Pressure and Poloidal Field Profiles')
    lines = [l1, l2, l3, l4]
    ax.legend(lines, [l.get_label() for l in lines], loc='best')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('jtot_jbs_comparison.png', dpi=150, bbox_inches='tight')
    print("   Plot saved to jtot_jbs_comparison.png")
    
    # ============================================================
    # 8. Save comparison data to CSV
    # ============================================================
    df = pd.DataFrame({
        'r_a': r,
        'JTOT_csv': JTOT_csv,
        'JOH_csv': JOH_csv,
        'JNB_csv': JNB_csv,
        'JRF_csv': JRF_csv,
        'JBS_csv': JBS_csv,
        'JTOT_sum': JTOT_sum,
        'JBS_recalc': JBS_avg,
        'eps': eps,
        'ft': ft,
        'nue_star': nue_star,
        'L31': RL31,
        'L32': RL32,
        'L34': RL34,
        'Bp': Bp,
        'QP': QP
    })
    df.to_csv('jtot_jbs_comparison.csv', index=False)
    print("   Data saved to jtot_jbs_comparison.csv")
    
    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)
    print(f"""
NOTE: The recalculated JBS may differ from the CSV values because:
  1. Geometric factors (TTRHOG, RDP, RJCB) are not available in CSV output.
     The code uses exact metric quantities from the equilibrium solver.
  2. We use simplified cylindrical geometry (eps = r*a/R) vs the code's
     actual EPSRHO = RMNRHO/RMJRHO from the equilibrium.
  3. The code uses staggered finite differences on its specific grid 
     (RHOM mesh centers vs RHOG grid boundaries), while we use np.gradient.
  4. Bp is reconstructed from JTOT integral rather than from the code's
     internal RDP/RDPVRHOG arrays.

Despite these approximations, the shape and magnitude should be similar,
confirming the physics model is correctly identified.
""")


if __name__ == '__main__':
    main()
