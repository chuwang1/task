
#!/usr/bin/env python3
import re
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# -------------------------
# 1) CSV (EQGS2D export)
# -------------------------
def load_psirz_csv(filename='eqgs2d_01_PSIRZ_grid.csv'):
    """
    Read eqgs2d PSIRZ grid CSV written as:
      index: Z (float)
      columns: R (float)
      values: PSIRZ(Z,R)
    Return:
      R (nw,), Z (nh,), psi (nh,nw)
    """
    df = pd.read_csv(filename, index_col=0)
    R = df.columns.astype(float).to_numpy()
    Z = df.index.astype(float).to_numpy()
    psi = df.to_numpy()  # (nh, nw)
    return R, Z, psi


def load_separatrix_csv(filename='eqgs2d_07_separatrix.csv'):
    if not os.path.exists(filename):
        return None
    return pd.read_csv(filename)


def load_parameters_csv(filename='eqgs2d_08_parameters.csv'):
    if not os.path.exists(filename):
        return {}
    df = pd.read_csv(filename, skipinitialspace=True)
    params = {}
    for _, row in df.iterrows():
        name = row['parameter'].strip() if isinstance(row['parameter'], str) else row['parameter']
        params[name] = row
    return params


# -------------------------
# 2) GEQDSK (gfile_efit)
# -------------------------
_header1 = re.compile(r'^\s*(.*)\s+(\d+)\s+(\d+)\s*$', re.I)
_header2 = re.compile(r'^\s*(\S+)\s*(\S+)\s*(\S+)\s*(\S+)\s*(\S+)')


def _split_fortran_numbers(s: str) -> str:
    # separate negative numbers unless exponent notation
    s = re.sub(r'\n', ' ', s)
    s = re.sub(r'([^Ee])-', r'\1 -', s)
    s = re.sub(r'(\d)\-', r'\1 -', s)
    return s


def load_geqdsk_psirz(filename):
    """
    Minimal GEQDSK reader for R,Z grid and psirz.
    Return:
      R (nw,), Z (nh,), psi (nh,nw), boundary (rbbbs,zbbbs), mag axis (rmaxis,zmaxis),
      psimag, psibdy
    """
    with open(filename, 'r', encoding='utf-8', errors='ignore') as f:
        line = f.readline()
        m = _header1.match(line)
        if not m:
            raise ValueError(f"GEQDSK header parse failed: {line!r}")
        header = m.group(1).strip()
        nw = int(m.group(2))
        nh = int(m.group(3))

        # line 2
        line = _split_fortran_numbers(f.readline())
        m = _header2.match(line)
        if not m:
            raise ValueError(f"GEQDSK line2 parse failed: {line!r}")
        rdim, zdim, rcentr, rleft, zmid = map(float, m.groups())

        # line 3
        line = _split_fortran_numbers(f.readline())
        m = _header2.match(line)
        if not m:
            raise ValueError(f"GEQDSK line3 parse failed: {line!r}")
        rmaxis, zmaxis, psimag, psibdy, bcentr = map(float, m.groups())

        # line 4 (only currentA first number used, rest ignored)
        line = _split_fortran_numbers(f.readline())
        m = _header2.match(line)
        if not m:
            # some files may have fewer fields; do a fallback
            parts = line.split()
            currentA = float(parts[0]) if parts else np.nan
        else:
            currentA = float(m.group(1))

        # remaining content -> numbers
        rest = f.read()
    rest = _split_fortran_numbers(rest)
    fields = re.split(r'\s+', rest.strip())
    numbers = []
    for t in fields:
        try:
            numbers.append(float(t))
        except Exception:
            pass
    numbers = np.asarray(numbers, dtype=float)

    ia = 0
    fpol = numbers[ia:ia+nw]; ia += nw
    pres = numbers[ia:ia+nw]; ia += nw
    ffprime = numbers[ia:ia+nw]; ia += nw
    pprime = numbers[ia:ia+nw]; ia += nw

    # psirz stored as nw*nh values; we want (nh, nw) to match meshgrid(R,Z)
    psirz_flat = numbers[ia:ia+nw*nh]; ia += nw*nh
    # common convention: data is written with R as fastest index (i=1..nw) inside Z loop,
    # but implementations vary. The safest is try reshape to (nh,nw) first.
    psi = psirz_flat.reshape((nh, nw), order='C')

    qpsi = numbers[ia:ia+nw]; ia += nw
    nbbbs = int(numbers[ia]); ia += 1
    limitr = int(numbers[ia]); ia += 1

    rbbbs = numbers[ia:ia+2*nbbbs:2]
    zbbbs = numbers[ia+1:ia+2*nbbbs:2]
    ia += 2*nbbbs

    # limiter geometry not used here, but consume if present
    # rlim = numbers[ia:ia+2*limitr:2]
    # zlim = numbers[ia+1:ia+2*limitr:2]
    # ia += 2*limitr

    # build grids
    R = np.linspace(rleft, rleft + rdim, nw)
    Z = np.linspace(zmid - 0.5*zdim, zmid + 0.5*zdim, nh)

    meta = {
        "header": header,
        "nw": nw, "nh": nh,
        "rdim": rdim, "zdim": zdim,
        "rleft": rleft, "zmid": zmid,
        "rmaxis": rmaxis, "zmaxis": zmaxis,
        "psimag": psimag, "psibdy": psibdy,
        "bcentr": bcentr,
        "currentA": currentA
    }
    return R, Z, psi, (rbbbs, zbbbs), (rmaxis, zmaxis), meta


# -------------------------
# 3) Comparison plot
# -------------------------
def compare_psirz_contours(
    csv_psirz='eqgs2d_01_PSIRZ_grid.csv',
    gfile='gfile_efit',
    levels=20,
    use_common_levels=True,
    show_boundary=False,
    show_separatrix_csv=True,
    savefig=None,
):
    # CSV
    R_csv, Z_csv, psi_csv = load_psirz_csv(csv_psirz)

    # GEQDSK
    R_g, Z_g, psi_g, (rbbbs, zbbbs), (rmaxis, zmaxis), meta = load_geqdsk_psirz(gfile)

    # choose contour levels
    if use_common_levels:
        vmin = max(np.nanmin(psi_csv), np.nanmin(psi_g))
        vmax = min(np.nanmax(psi_csv), np.nanmax(psi_g))
        if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin >= vmax:
            # fallback
            vmin = np.nanmin(psi_csv)
            vmax = np.nanmax(psi_csv)
        lev = np.linspace(vmin, vmax, levels)
    else:
        lev = levels  # let matplotlib decide if int, or pass array yourself

    fig, ax = plt.subplots(figsize=(12, 10))

    # plot CSV as solid
    R2D_csv, Z2D_csv = np.meshgrid(R_csv, Z_csv)
    c1 = ax.contour(R2D_csv, Z2D_csv, psi_csv, levels=lev, linewidths=1.0, linestyles='-', colors='C0')
    # plot GEQDSK as dashed
    R2D_g, Z2D_g = np.meshgrid(R_g, Z_g)
    c2 = ax.contour(R2D_g, Z2D_g, psi_g, levels=lev, linewidths=1.0, linestyles='--', colors='C3')

    # labels (optional; too many labels can clutter)
    # ax.clabel(c1, inline=True, fontsize=8, fmt="%.3g")
    # ax.clabel(c2, inline=True, fontsize=8, fmt="%.3g")

    # boundary from GEQDSK
    if show_boundary and rbbbs.size > 0:
        ax.plot(rbbbs, zbbbs, 'k-', linewidth=2.0, label='GEQDSK boundary')

    # separatrix from CSV export (if present)
    if show_separatrix_csv:
        sep = load_separatrix_csv()
        if sep is not None:
            # keep your original column names
            if all(col in sep.columns for col in ['RSU_m', 'ZSU_m', 'RSW_m', 'ZSW_m']):
                ax.plot(sep['RSU_m'], sep['ZSU_m'], 'C2-', linewidth=2, label='CSV separatrix (upper)')
                ax.plot(sep['RSW_m'], sep['ZSW_m'], 'C2--', linewidth=2, label='CSV separatrix (lower)')

    # magnetic axes: CSV params + GEQDSK
    params = load_parameters_csv()
    if 'RAXIS' in params and 'ZAXIS' in params:
        ax.plot(float(params['RAXIS']['value']), float(params['ZAXIS']['value']),
                marker='*', markersize=14, linestyle='None', color='C2', label='CSV axis')

    ax.plot(rmaxis, zmaxis, marker='*', markersize=14, linestyle='None', color='C3', label='GEQDSK axis')

    # legend proxies for contour line styles
    ax.plot([], [], linestyle='-', color='C0', label='PSI from CSV (solid)')
    ax.plot([], [], linestyle='--', color='C3', label='PSI from GEQDSK (dashed)')

    ax.set_xlabel('R (m)')
    ax.set_ylabel('Z (m)')
    ax.set_title('RZ Psi comparison: CSV vs GEQDSK')
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='best')

    if savefig:
        fig.savefig(savefig, dpi=150, bbox_inches='tight')
        print(f"Saved: {savefig}")

    return fig, ax


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description="Compare RZ psi contours from EQGS2D CSV and GEQDSK gfile")
    p.add_argument("--csv", default="eqgs2d_01_PSIRZ_grid.csv", help="EQGS2D PSIRZ grid csv")
    p.add_argument("--gfile", required=True, help="GEQDSK file (gfile_efit)")
    p.add_argument("--levels", type=int, default=25)
    p.add_argument("--no-common-levels", action="store_true", help="do not force common contour levels")
    p.add_argument("--save", default=None, help="output png filename")
    args = p.parse_args()

    compare_psirz_contours(
        csv_psirz=args.csv,
        gfile=args.gfile,
        levels=args.levels,
        use_common_levels=not args.no_common_levels,
        savefig=args.save
    )
    plt.show()