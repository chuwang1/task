#!/usr/bin/env python3
"""
CDBM (Current Diffusivity Ballooning Mode) Transport Model
Simple Python implementation to calculate chi profile shape.

Includes a closer reproduction of the current Fortran MDLKAI=134 setup.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import pandas as pd

# Physical constants
AEE = 1.602176487e-19   # elementary charge [C]
AME = 9.10938215e-31    # electron mass [kg]
AMP = 1.672621637e-27   # proton mass [kg]
VC = 2.99792458e8       # speed of light [m/s]
RMU0 = 4.0e-7 * np.pi   # permeability [H/m]
EPS0 = 1.0 / (VC**2 * RMU0)  # permittivity [F/m]
RKEV = 1.0e3 * AEE      # keV to J

def trcofs(shear, alpha, curv):
    """
    Form factor in CDBM model

    Parameters:
    -----------
    shear : float or array
        Magnetic shear s = (r/q)(dq/dr)
    alpha : float or array
        Normalized pressure gradient
    curv : float or array
        Magnetic curvature

    Returns:
    --------
    fs : float or array
        Form factor
    """
    shear_min = 0.01

    # Handle arrays
    shear = np.atleast_1d(shear)
    alpha = np.atleast_1d(alpha)
    curv = np.atleast_1d(curv)

    fs = np.zeros_like(shear)

    for i in range(len(shear)):
        if alpha[i] >= 0:
            sa = shear[i] - alpha[i]

            if sa >= 0:
                fs1 = (1.0 + 9.0 * np.sqrt(2.0) * sa**2.5) / \
                      (np.sqrt(2.0) * (1.0 - 2.0*sa + 3.0*sa*sa + 2.0*sa**3))
            else:
                denom = 2.0 * (1.0 - 2.0*sa) * (1.0 - 2.0*sa + 3.0*sa*sa)
                if denom > 0:
                    fs1 = 1.0 / np.sqrt(denom)
                else:
                    fs1 = 0.0

            if curv[i] > 0:
                if abs(shear[i]) > shear_min:
                    fs2 = np.sqrt(curv[i])**3 / (shear[i]**2)
                else:
                    fs2 = np.sqrt(curv[i])**3 / shear_min**2
            else:
                fs2 = 0.0
        else:
            # alpha < 0 case
            sa = alpha[i] - shear[i]
            if sa >= 0:
                fs1 = (1.0 + 9.0 * np.sqrt(2.0) * sa**2.5) / \
                      (np.sqrt(2.0) * (1.0 - 2.0*sa + 3.0*sa*sa + 2.0*sa**3))
            else:
                denom = 2.0 * (1.0 - 2.0*sa) * (1.0 - 2.0*sa + 3.0*sa*sa)
                if denom > 0:
                    fs1 = 1.0 / np.sqrt(denom)
                else:
                    fs1 = 0.0

            if curv[i] < 0:
                if abs(shear[i]) > shear_min:
                    fs2 = np.sqrt(-curv[i])**3 / (shear[i]**2)
                else:
                    fs2 = np.sqrt(-curv[i])**3 / shear_min**2
            else:
                fs2 = 0.0

        fs[i] = max(fs1, fs2)

    return fs


def fexb(x, shear, alpha):
    """Fortran `FEXB` from `trcoef.f90`."""
    x = np.asarray(x)
    shear = np.asarray(shear)
    alpha = np.asarray(alpha)

    alpha_abs = np.where(np.abs(alpha) < 1.0e-3, 1.0e-3, np.abs(alpha))
    beta = (0.5 * alpha_abs ** (-0.602) *
            (13.018 - 22.28915 * shear + 17.018 * shear**2) /
            (1.0 - 0.277584 * shear + 1.42913 * shear**2))

    a = -10.0 / 3.0 * alpha + 16.0 / 3.0
    gamma = np.empty_like(shear, dtype=float)
    mask_neg = shear < 0.0
    denom_neg = np.maximum(1.0 - shear - 2.0 * shear**2 - 3.0 * shear**3, 1.0e-12)
    gamma[mask_neg] = 1.0 / (1.1 * np.sqrt(denom_neg[mask_neg])) + 0.75
    gamma[~mask_neg] = ((1.0 - 0.5 * shear[~mask_neg]) /
                        (1.1 - 2.0 * shear[~mask_neg] + a[~mask_neg] * shear[~mask_neg]**2 +
                         4.0 * shear[~mask_neg]**3) + 0.75)

    xg = np.empty_like(x, dtype=float)
    mask_low = gamma <= -20.0
    mask_high = gamma >= 20.0
    mask_mid = ~(mask_low | mask_high)
    xg[mask_low] = 0.0
    xg[mask_high] = 1.0e10
    xg[mask_mid] = x[mask_mid] ** gamma[mask_mid]
    arg = -beta * xg
    return np.where(arg <= -20.0, 0.0,
           np.where(arg >= 20.0, 1.0e10, np.exp(arg)))


def smooth_boxcar_5(arr):
    """Match the 5-point smoothing used for `S_HM` in `trcoef.f90`."""
    arr = np.asarray(arr, dtype=float)
    out = arr.copy()
    if len(arr) < 3:
        return out
    out[0] = arr[0]
    out[1] = (arr[0] + arr[1] + arr[2]) / 3.0
    if len(arr) > 4:
        for i in range(2, len(arr) - 2):
            out[i] = (arr[i-2] + arr[i-1] + arr[i] + arr[i+1] + arr[i+2]) / 5.0
    out[-2] = (arr[-3] + arr[-2] + arr[-1]) / 3.0
    out[-1] = arr[-1]
    return out


def calculate_cdbm_chi(BB, RR, rs, qp, shear, ne, dpdr, rhoni,
                       dvexbdr=None, calf=1.0, ckap=1.0, cexb=1.0, model=2,
                       shear_factor=1.0, shear_override=None,
                       shear_min=None, fs_floor=None):
    """
    Calculate CDBM thermal diffusivity

    Parameters:
    -----------
    BB : float
        Magnetic field strength [T]
    RR : float
        Major radius [m]
    rs : array
        Minor radius at each grid point [m]
    qp : array
        Safety factor profile
    shear : array
        Magnetic shear profile s = (r/q)(dq/dr)
    ne : array
        Electron density [m^-3]
    dpdr : array
        Pressure gradient [Pa/m]
    rhoni : array
        Ion mass density [kg/m^3]
    dvexbdr : array, optional
        ExB drift velocity gradient [1/s]
    calf : float
        Factor for s-alpha effects
    ckap : float
        Factor for magnetic curvature effects
    cexb : float
        Factor for ExB drift effects
    model : int
        Model ID (0-5)
    shear_factor : float
        Factor to scale magnetic shear (0=no shear, 1=full shear, default=1.0)
        - Set to 0 to completely eliminate shear effect
        - Set to 0.5 to reduce shear effect by half
        - Set to 1.0 for standard CDBM behavior
    shear_override : float or array, optional
        If provided, use this value(s) instead of the input shear profile
        - Useful for testing specific shear values
    shear_min : float, optional
        Minimum shear value to prevent s from going too negative
        - Set to 0.5 to eliminate jump (Method 1)
    fs_floor : float, optional
        Minimum form factor value to prevent fs from dropping too low
        - Set to 0.3 to eliminate jump (Method 4)

    Returns:
    --------
    chi_cdbm : array
        Thermal diffusivity [m^2/s]
    """
    ckcdbm = 20.0  # Fixed numerical factor

    # Apply shear modification
    if shear_override is not None:
        # Use override value(s) instead of input shear
        if np.isscalar(shear_override):
            shear_eff = np.full_like(shear, shear_override)
        else:
            shear_eff = np.asarray(shear_override)
    else:
        # Apply scaling factor to input shear
        shear_eff = shear * shear_factor

    # Apply minimum shear threshold (Method 1 to eliminate jump)
    if shear_min is not None:
        shear_eff = np.maximum(shear_eff, shear_min)

    # Alfven velocity
    va = np.sqrt(BB**2 / (RMU0 * rhoni))

    # Square of plasma frequency
    wpe2 = ne * AEE**2 / (AME * EPS0)

    # Square of collisionless skin depth
    delta2 = VC**2 / wpe2

    # Normalized pressure gradient (Shafranov shift factor)
    alpha = -2.0 * RMU0 * qp**2 * RR / BB**2 * dpdr

    # Magnetic curvature
    curv = -(rs / RR) * (1.0 - 1.0 / (qp**2))

    # Rotational shear (use effective shear with modification)
    shearl = np.sqrt(shear_eff**2 + 0.1**2)

    if dvexbdr is None:
        dvexbdr = np.zeros_like(rs)
    wexb = -qp * RR / (shearl * va) * dvexbdr

    # Model-dependent factors
    # fk: elongation factor
    if model % 2 == 0:
        fk = 1.0
    else:
        # Would need elongation input for CDBM05
        fk = 1.0

    # fe: ExB shear reduction factor
    model_type = (model // 2) % 3
    if model_type == 0:
        fe = 1.0
    elif model_type == 1:
        # Weak ExB shear (model=2,3)
        fe = 1.0 / (1.0 + cexb * wexb**2)
    else:
        # Strong ExB shear (model=4,5)
        fe = np.ones_like(wexb)  # Simplified

    # Form factor (use effective shear with modification)
    fs = trcofs(shear_eff, calf * alpha, ckap * curv)

    # Apply form factor floor (Method 4 to eliminate jump)
    if fs_floor is not None:
        fs = np.maximum(fs, fs_floor)

    # CDBM chi
    chi_cdbm = ckcdbm * fs * fk * fe * np.abs(alpha)**1.5 * delta2 * va / (qp * RR)

    return chi_cdbm, alpha, fs, fe


def calculate_mdlkai134_chi(BB, RR, RA, rs, qp, shear, ne, dpdr, rhoni,
                            alpha_input=None, dvexbdr=None, calf=1.0, cweb=1.0,
                            ck0=12.0, smooth_shear=True):
    """
    Best-effort reproduction of current Fortran `MDLKAI=134` setup.

    Matches the logic around `CASE(130:139)` in `trcoef.f90`:
    - `MODEL = MDLKAI - 130 = 4`
    - optional 5-point smoothing for `S_HM`
    - `SHEARL = MAX(S_HM, -0.5)` when smoothing is enabled
    - `SL = S^2 + 0.1^2`
    - `cexb = CWEB * FEXB(abs(WE1), S, ALPHA)`

    Notes:
    - This script still uses the simplified Python CDBM kernel for `MODEL=4`.
    - If `dvexbdr` is unavailable, zero ExB shear is assumed.
    """
    if dvexbdr is None:
        dvexbdr = np.zeros_like(rs)

    if alpha_input is None:
        alpha = -2.0 * RMU0 * qp**2 * RR / BB**2 * dpdr
    else:
        alpha = np.asarray(alpha_input, dtype=float)
    va = np.sqrt(BB**2 / (RMU0 * rhoni))
    wpe2 = ne * AEE**2 / (AME * EPS0)
    delta2 = VC**2 / wpe2
    curv = -(rs / RR) * (1.0 - 1.0 / (qp**2))

    shear_sm = smooth_boxcar_5(shear) if smooth_shear else np.asarray(shear, dtype=float)
    shearl = np.maximum(shear_sm, -0.5) if smooth_shear else np.asarray(shear, dtype=float)

    sl = shear**2 + 0.1**2
    dve = dvexbdr * RA
    we1 = -qp * RR / (sl * va) * dve
    cexb = cweb * fexb(np.abs(we1), shear, alpha)

    fs = trcofs(shearl, calf * alpha, curv)
    fe = np.ones_like(fs)
    chi_cdbm = 12.0 * fs * fe * np.abs(alpha)**1.5 * delta2 * va / (qp * RR)
    chi_e = (ck0 / 12.0) * chi_cdbm

    return {
        'chi_e': chi_e,
        'alpha': alpha,
        'fs': fs,
        'fe': fe,
        'cexb': cexb,
        'we1': we1,
        'shear_used': shearl,
        'shear_smoothed': shear_sm,
        'curv': curv,
    }


def read_chi_corediv(filename='chi.dat'):
    """Read chi.dat file"""
    data = []
    with open(filename, 'r') as f:
        _ = f.readline()  # Skip header
        for line in f:
            parts = line.split()
            if len(parts) >= 5:
                data.append([float(x) for x in parts[:5]])
    data = np.array(data)
    return {
        'r': data[:, 0],      # r in meters
        'ne': data[:, 1],     # electron density [m^-3]
        'ni': data[:, 2],     # ion density [m^-3]
        'chi_e': data[:, 3],  # electron chi [m^2/s]
        'chi_i': data[:, 4],  # ion chi [m^2/s]
    }


def read_omfit_profiles(filename='omfit_profiles_for_tr.csv'):
    """Read OMFIT profiles CSV file"""
    df = pd.read_csv(filename)
    return df


def read_tr_profile_csv(filename, prefix, time_col=-1):
    """
    Read TR profile CSV file (e.g., tr_data_053.csv for QP)

    TR CSV format:
    - Row 1: Title
    - Row 2: Column names (X, prefix_1, prefix_2, ...)
    - Row 3+: Data where X is r/a, prefix_i is value at time step i

    Parameters:
    -----------
    filename : str
        CSV filename
    prefix : str
        Variable prefix (e.g., 'QP', 's', 'alpha', 'NE', 'TE', 'AKD')
    time_col : int
        Which time column to use (-1 = last)

    Returns:
    --------
    rho : array
        Radial coordinate r/a
    profile : array
        Profile values at specified time
    """
    df = pd.read_csv(filename, skiprows=1)  # Skip title row
    df.columns = df.columns.str.strip()

    # X column is r/a
    rho = df['X'].values

    # Get columns with the prefix: try numbered columns first (time-series),
    # then fall back to exact column name (snapshot)
    cols = [c for c in df.columns if c.startswith(prefix + '_')]
    if len(cols) > 0:
        # Time-series format: prefix_1, prefix_2, ...
        cols = sorted(cols, key=lambda x: int(x.split('_')[-1]))
        profile = df[cols[time_col]].values
    elif prefix in df.columns:
        # Snapshot format: single column named exactly 'prefix'
        profile = df[prefix].values
    else:
        raise ValueError(f"No columns found for '{prefix}' in {filename} "
                         f"(columns: {list(df.columns)})")

    return rho, profile


def find_csv_by_title(title_keyword, search_dir='.'):
    """Find a TR CSV file by searching for a keyword in its title line."""
    import glob, os
    for f in sorted(glob.glob(os.path.join(search_dir, 'tr_data_*.csv'))):
        with open(f, 'r') as fh:
            title = fh.readline().strip()
        if title_keyword in title:
            return f
    return None


def main():
    # Machine parameters for CFEDR
    RR = 8.04  # Major radius [m]
    RA = 2.44  # Minor radius [m]
    BB = 5.3   # Toroidal field [T]

    # Try to read profiles from TR CSV files
    try:
        from scipy.interpolate import interp1d

        # Auto-find CSV files by title content
        csv_qp    = find_csv_by_title('@QP  vs r@')
        csv_s     = find_csv_by_title('@s  vs r@')
        csv_alpha = find_csv_by_title('@alpha  vs r@')
        csv_ne    = find_csv_by_title('@NE [')
        csv_te    = find_csv_by_title('@TE [')
        csv_td    = find_csv_by_title('@TD [')

        missing = []
        for name, path in [('QP', csv_qp), ('s', csv_s), ('alpha', csv_alpha),
                           ('NE', csv_ne), ('TE', csv_te), ('TD', csv_td)]:
            if path is None:
                missing.append(name)
        if missing:
            raise FileNotFoundError(f"Could not find CSV files for: {missing}")

        print(f"Auto-detected CSV files:")
        print(f"  QP    -> {csv_qp}")
        print(f"  s     -> {csv_s}")
        print(f"  alpha -> {csv_alpha}")
        print(f"  NE    -> {csv_ne}")
        print(f"  TE    -> {csv_te}")
        print(f"  TD    -> {csv_td}")

        # Read q, shear, alpha from TR output (use QP grid as reference)
        rho, q = read_tr_profile_csv(csv_qp, 'QP')
        rho_s, shear = read_tr_profile_csv(csv_s, 's')
        rho_a, alpha_tr = read_tr_profile_csv(csv_alpha, 'alpha')
        rho_ne, ne = read_tr_profile_csv(csv_ne, 'NE')
        rho_Te, Te = read_tr_profile_csv(csv_te, 'TE')
        rho_Ti, Ti = read_tr_profile_csv(csv_td, 'TD')

        # Interpolate all to common grid (use rho from QP)
        f_s = interp1d(rho_s, shear, bounds_error=False, fill_value='extrapolate')
        f_alpha = interp1d(rho_a, alpha_tr, bounds_error=False, fill_value='extrapolate')
        f_ne = interp1d(rho_ne, ne, bounds_error=False, fill_value='extrapolate')
        f_Te = interp1d(rho_Te, Te, bounds_error=False, fill_value='extrapolate')
        f_Ti = interp1d(rho_Ti, Ti, bounds_error=False, fill_value='extrapolate')

        shear = f_s(rho)
        alpha_tr = f_alpha(rho)
        ne = f_ne(rho) * 1e20  # Convert to m^-3
        Te = f_Te(rho)
        Ti = f_Ti(rho)

        # Update rs based on actual rho
        rs = rho * RA

        # Try to read chi from TR
        try:
            csv_akd = find_csv_by_title('AKD,AKNCD,AKDWD')
            csv_akdwe = find_csv_by_title('AKE,AKNCE,AKDWE')
            if csv_akd is None:
                raise FileNotFoundError("AKD CSV not found")
            print(f"  AKD   -> {csv_akd}")
            rho_chi, chi_tr_raw = read_tr_profile_csv(csv_akd, 'AKD')
            f_chi = interp1d(rho_chi, chi_tr_raw, bounds_error=False, fill_value='extrapolate')
            chi_tr = f_chi(rho)
            if csv_akdwe is not None:
                print(f"  AKDWE -> {csv_akdwe}")
                rho_akdwe, akdwe_raw = read_tr_profile_csv(csv_akdwe, 'AKDWE')
                f_akdwe = interp1d(rho_akdwe, akdwe_raw, bounds_error=False, fill_value='extrapolate')
                akdwe_tr = f_akdwe(rho)
            else:
                akdwe_tr = None
        except Exception as e:
            print(f"  AKD   -> not found ({e})")
            chi_tr = None
            akdwe_tr = None

        print("Read profiles from TR CSV files (last time step):")
        print(f"  rho points: {len(rho)}")
        print(f"  q range: {q.min():.2f} - {q.max():.2f}")
        print(f"  s range: {shear.min():.3f} - {shear.max():.3f}")
        print(f"  alpha range: {alpha_tr.min():.3f} - {alpha_tr.max():.3f}")
        print(f"  ne range: {ne.min():.2e} - {ne.max():.2e} m^-3")
        print(f"  Te range: {Te.min():.2f} - {Te.max():.2f} keV")

    except Exception as e:
        print(f"Could not read TR CSV files: {e}")
        print("Using analytical profiles.")

        # Define grid
        NRMAX = 41
        rho = np.linspace(0, 1, NRMAX)
        rs = rho * RA

        # Analytical profiles
        q0 = 1.05
        q95 = 3.5
        q = q0 + (q95 - q0) * rho**2 * (1 + 0.5 * rho**2)

        ne0 = 1.8e20
        ne = ne0 * (1 - 0.8 * rho**2)

        Te0 = 25.0
        Te = Te0 * (1 - rho**2)**1.5
        Ti = Te * 0.9

        # Calculate shear from q
        dqdr = np.gradient(q, rs)
        shear = np.zeros_like(rho)
        shear[1:] = (rs[1:] / q[1:]) * dqdr[1:]
        shear[0] = shear[1]

        alpha_tr = None
        chi_tr = None
        akdwe_tr = None

    # Read OMFIT chi for comparison
    chi_omfit_e = None
    chi_omfit_i = None
    try:
        data = read_chi_corediv('chi_div_grad2_drdrho.dat')
        # Interpolate to TR grid
        # print(data)
        from scipy.interpolate import interp1d
        rho_omfit = data['r']
       
        f_chi_e = interp1d(rho_omfit, data['chi_e'], bounds_error=False, fill_value='extrapolate')
        f_chi_i = interp1d(rho_omfit, data['chi_i'], bounds_error=False, fill_value='extrapolate')
        chi_omfit_e = f_chi_e(rho)
        chi_omfit_i = f_chi_i(rho)
        print(f"Read OMFIT chi from chi.dat")
    except Exception as e:
        print(f"Could not read chi.dat: {e}")

    # Pressure and gradient (for CDBM calculation)
    P = ne * (Te + Ti) * RKEV  # Total pressure [Pa]
    dpdr = np.gradient(P, rs)  # Pressure gradient [Pa/m]
    dpdr[0] = 0  # Fix center

    # Ion mass density (assume 50% D, 50% T)
    ni = ne * 0.9  # Allow for impurities
    rhoni = ni * (2.0 + 3.0) / 2.0 * AMP  # Average ion mass

    # Calculate CDBM chi with standard shear (legacy simplified model)
    chi_cdbm, alpha, fs, fe = calculate_cdbm_chi(
        BB, RR, rs, q, shear, ne, dpdr, rhoni, model=2, shear_factor=1.0
    )

    # Closer reproduction of the current Fortran MDLKAI=134 branch
    mdlkai134 = calculate_mdlkai134_chi(
        BB, RR, RA, rs, q, shear, ne, dpdr, rhoni,
        alpha_input=alpha_tr, dvexbdr=None, calf=1.0, cweb=1.0, ck0=12.0, smooth_shear=True
    )

    # Calculate CDBM chi with reduced shear (shear_factor=0.5)
    chi_cdbm_half, _, fs_half, _ = calculate_cdbm_chi(
        BB, RR, rs, q, shear, ne, dpdr, rhoni, model=2, shear_factor=0.5
    )

    # Calculate CDBM chi with no shear (shear_factor=0.0)
    chi_cdbm_noshear, _, fs_noshear, _ = calculate_cdbm_chi(
        BB, RR, rs, q, shear, ne, dpdr, rhoni, model=2, shear_factor=0.0
    )

    # Method 1: Eliminate jump with shear_min=0.5
    chi_cdbm_smin, _, fs_smin, _ = calculate_cdbm_chi(
        BB, RR, rs, q, shear, ne, dpdr, rhoni, model=2, shear_min=1
    )

    # Method 4: Eliminate jump with fs_floor=0.3
    chi_cdbm_fsfloor, _, fs_fsfloor, _ = calculate_cdbm_chi(
        BB, RR, rs, q, shear, ne, dpdr, rhoni, model=2, fs_floor=0.3
    )

    # Apply CK0/12 factor
    # TR default: CK0=12, CK1=12, so effective factor = 12/12 = 1.0
    CK0 = 12.0  # TR default value
    chi_e = (CK0 / 12.0) * chi_cdbm  # = 1.0 * chi_cdbm when CK0=12
    chi_e_134 = mdlkai134['chi_e']
    chi_e_half = (CK0 / 12.0) * chi_cdbm_half
    chi_e_noshear = (CK0 / 12.0) * chi_cdbm_noshear
    chi_e_smin = (CK0 / 12.0) * chi_cdbm_smin
    chi_e_fsfloor = (CK0 / 12.0) * chi_cdbm_fsfloor

    # Limit extreme values
    chi_e = np.clip(chi_e, 0, 100)
    chi_e_134 = np.clip(chi_e_134, 0, 100)
    chi_e_half = np.clip(chi_e_half, 0, 100)
    chi_e_noshear = np.clip(chi_e_noshear, 0, 100)
    chi_e_smin = np.clip(chi_e_smin, 0, 100)
    chi_e_fsfloor = np.clip(chi_e_fsfloor, 0, 100)

    # Create figure
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))

    # Plot 1: Chi profile with jump elimination methods
    ax1 = axes[0, 0]
    ax1.plot(rho, chi_e, 'b-', linewidth=2, label=r'$\chi$ (Legacy model=2)')
    ax1.plot(rho, chi_e_134, 'r-', linewidth=2, label=r'$\chi$ (Fortran-like MDLKAI=134)')
    if akdwe_tr is not None:
        ax1.plot(rho, akdwe_tr, 'g--', linewidth=2, label=r'$AKDWE$ (TR)')
    # ax1.plot(rho, chi_e_smin, 'r--', linewidth=2, label=r'$\chi$ (s_min=0.5)')
    # ax1.plot(rho, chi_e_fsfloor, 'g--', linewidth=2, label=r'$\chi$ (fs_floor=0.3)')
    if chi_omfit_e is not None:
        ax1.plot(rho, chi_omfit_e, 'k:', linewidth=2, label=r'$\chi_e$ (OMFIT)')
    ax1.set_xlabel('rho')
    ax1.set_ylabel(r'$\chi$ [m$^2$/s]')
    ax1.set_title('Chi: Jump Elimination Methods')
    ax1.legend(fontsize=8)
    ax1.grid(True)
    ax1.set_yscale('log')
    ax1.set_xlim(0, 1)

    # Plot 2: Alpha (normalized pressure gradient)
    ax2 = axes[0, 1]
    ax2.plot(rho, alpha, 'b-', linewidth=2, label=r'$\alpha$ (CDBM calc)')
    if alpha_tr is not None:
        ax2.plot(rho, alpha_tr, 'r--', linewidth=2, label=r'$\alpha$ (TR output)')
    ax2.set_xlabel('rho')
    ax2.set_ylabel(r'$\alpha$')
    ax2.set_title('Normalized Pressure Gradient')
    ax2.legend()
    ax2.grid(True)
    ax2.set_xlim(0, 1)

    # Plot 3: Form factor fs with different shear factors
    ax3 = axes[0, 2]
    ax3.plot(rho, fs, 'b-', linewidth=2, label='legacy fs')
    ax3.plot(rho, mdlkai134['fs'], 'r-', linewidth=2, label='MDLKAI=134 fs')
    ax3.plot(rho, fs_half, 'c-', linewidth=2, label='s×0.5')
    ax3.plot(rho, fs_noshear, 'm-', linewidth=2, label='s×0.0')
    ax3.set_xlabel('rho')
    ax3.set_ylabel(r'$f_s$')
    ax3.set_title('CDBM Form Factor')
    ax3.legend(fontsize=8)
    ax3.grid(True)
    ax3.set_xlim(0, 1)

    # Plot 4: Safety factor
    ax4 = axes[1, 0]
    ax4.plot(rho, q, 'b-', linewidth=2)
    ax4.set_xlabel('rho')
    ax4.set_ylabel('q')
    ax4.set_title('Safety Factor Profile')
    ax4.grid(True)
    ax4.set_xlim(0, 1)

    # Plot 5: Magnetic shear
    ax5 = axes[1, 1]
    ax5.plot(rho, shear, 'm-', linewidth=2, label='raw s')
    ax5.plot(rho, mdlkai134['shear_smoothed'], 'k--', linewidth=2, label='S_HM (5-pt)')
    ax5.set_xlabel('rho')
    ax5.set_ylabel('s = (r/q)(dq/dr)')
    ax5.set_title('Magnetic Shear Profile')
    ax5.grid(True)
    ax5.set_xlim(0, 1)
    ax5.legend(fontsize=8)

    # Plot 6: Density and pressure
    ax6 = axes[1, 2]
    ax6.plot(rho, ne/1e20, 'b-', linewidth=2, label=r'$n_e$ [$10^{20}$ m$^{-3}$]')
    ax6.set_xlabel('rho')
    ax6.set_ylabel(r'$n_e$ [$10^{20}$ m$^{-3}$]')
    ax6.set_title('Density Profile')
    ax6.legend()
    ax6.grid(True)
    ax6.set_xlim(0, 1)

    plt.tight_layout()
    plt.savefig('cdbm_chi_profile.png', dpi=150)
    print("Plot saved to cdbm_chi_profile.png")

    print("\nCDBM Chi - Jump Elimination Methods:")
    print("-" * 80)
    print(f"{'rho':<8} {'Original':<12} {'s_min=0.5':<12} {'fs_floor=0.3':<12} {'Jump eliminated?'}")
    print("-" * 80)
    for r_val in [0.25, 0.30, 0.32, 0.35, 0.40]:
        idx = np.argmin(np.abs(rho - r_val))
        orig = chi_e[idx]
        smin = chi_e_smin[idx]
        fsfl = chi_e_fsfloor[idx]
        # Check if jump is eliminated (chi should not drop below 0.2)
        elim = "Yes" if (smin > 0.2 and fsfl > 0.2) else "Partial"
        print(f"{rho[idx]:<8.3f} {orig:<12.3f} {smin:<12.3f} {fsfl:<12.3f} {elim}")
    print("-" * 80)
    print("\nRecommendation: Use shear_min=0.5 or fs_floor=0.3 to eliminate chi jump")
    print("\nFortran-like MDLKAI=134 summary:")
    print(f"  chi range: {chi_e_134.min():.3f} - {chi_e_134.max():.3f} m^2/s")
    print(f"  fs range : {mdlkai134['fs'].min():.3e} - {mdlkai134['fs'].max():.3e}")
    print(f"  cexb range: {mdlkai134['cexb'].min():.3e} - {mdlkai134['cexb'].max():.3e}")
    if akdwe_tr is not None:
        edge = (rho >= 0.8) & (rho <= 0.9)
        rms = np.sqrt(np.mean((chi_e_134[edge] - akdwe_tr[edge])**2))
        print(f"  edge RMS vs TR AKDWE (0.8<=rho<=0.9): {rms:.3f} m^2/s")

    return rho, chi_e, chi_e_134, chi_e_smin, chi_e_fsfloor, alpha, fs


if __name__ == "__main__":
    rho, chi_e, chi_e_134, chi_e_smin, chi_e_fsfloor, alpha, fs = main()
