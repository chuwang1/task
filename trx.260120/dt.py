import numpy as np
import pandas as pd
import netCDF4 as nc4
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

E_ALPHA = 3.5e6 * 1.602176634e-19  # J

# ---------------------------
# TR table
# ---------------------------
T_tab = np.array([1, 2, 5, 10, 20, 50, 100, 200, 500, 1000], dtype=float)
sv_tab_cm3 = np.array([5.5e-21, 2.6e-19, 1.3e-17, 1.1e-16, 4.2e-16,
                       8.7e-16, 8.5e-16, 6.3e-16, 3.7e-16, 2.7e-16], dtype=float)
sv_tab = sv_tab_cm3 * 1e-6  # m^3/s

T_dense = np.linspace(20.0, 50.0, 500)
sv_dense = 10 ** np.interp(np.log10(T_dense), np.log10(T_tab), np.log10(sv_tab))

# ---------------------------
# OMFIT inferred <sigma v>
# ---------------------------
ds = nc4.Dataset('/Users/dengxiaoya/CFEDRSW/OMFIT_out/statefile_3.000000E+01.nc')
Ti_om = ds.variables['Ti'][:]
stfuse = ds.variables['stfuse'][:]
enion = ds.variables['enion'][:]
rho_om = ds.variables['rho_grid'][:]
rho_om = rho_om / rho_om.max()
ds.close()

nD_om = enion[0, :]
nT_om = enion[1, :]
mask_om = (nD_om > 0) & (nT_om > 0) & (stfuse > 0)
sv_om = stfuse[mask_om] / (nD_om[mask_om] * nT_om[mask_om])
Ti_om = Ti_om[mask_om]
rho_om_v = rho_om[mask_om]

zoom_om = (Ti_om >= 20.0) & (Ti_om <= 50.0)
Ti_om_z = Ti_om[zoom_om]
sv_om_z = sv_om[zoom_om]
rho_om_z = rho_om_v[zoom_om]
ord_om = np.argsort(Ti_om_z)

# ---------------------------
# TR inferred <sigma v> from outputs
# ---------------------------
tr = '/Users/dengxiaoya/TASK/latest/task/trx'
import glob, os

def find_csv(pattern):
    csvs = glob.glob(f"{tr}/tr_data_*.csv")
    idx = {}
    for c in csvs:
        try:
            num = int(c.split("_")[-1].replace(".csv",""))
            with open(c) as f: idx[num] = f.readline().strip()
        except:
            pass
    matches = sorted(n for n,t in idx.items() if pattern in t)
    if not matches: raise FileNotFoundError(f"Missing '{pattern}'")
    return f"{tr}/tr_data_{matches[-1]:03d}.csv"

src_path = find_csv("POH,PNB,PNF")
n_path   = find_csv("n(NS)")
T_path   = find_csv("T(NS)")

# PNF profile on its own grid (MW/m^3)
df_src = pd.read_csv(src_path, skiprows=1)
df_src.columns = df_src.columns.str.strip()
r_tr = df_src['X'].values
pnf_tr = df_src['PNF'].values * 1e6  # W/m^3

# nD, nT, TD, TT from snapshot files
df_n = pd.read_csv(n_path, skiprows=1)
df_n.columns = df_n.columns.str.strip()
df_T = pd.read_csv(T_path, skiprows=1)
df_T.columns = df_T.columns.str.strip()

# Interpolate to r_tr (though they are usually all on the same grid)
r_n, nd, nt = df_n['X'].values, df_n['nD'].values, df_n['nT'].values
r_T, td, tt = df_T['X'].values, df_T['TD'].values, df_T['TT'].values

nd_i = np.interp(r_tr, r_n, nd) * 1e20
nt_i = np.interp(r_tr, r_n, nt) * 1e20
td_i = np.interp(r_tr, r_T, td)
tt_i = np.interp(r_tr, r_T, tt)

# Use TR's effective temperature definition in SIGMAM: TI=(3*TD+2*TT)/5
Ti_eff_tr = (3.0 * np.abs(td_i) + 2.0 * np.abs(tt_i)) / 5.0
sv_tr_inf = pnf_tr / (nd_i * nt_i * E_ALPHA)

mask_tr = (nd_i > 0) & (nt_i > 0) & (pnf_tr > 0)
Ti_eff_tr = Ti_eff_tr[mask_tr]
sv_tr_inf = sv_tr_inf[mask_tr]
rho_tr = r_tr[mask_tr]

zoom_tr = (Ti_eff_tr >= 20.0) & (Ti_eff_tr <= 50.0)
Ti_tr_z = Ti_eff_tr[zoom_tr]
sv_tr_z = sv_tr_inf[zoom_tr]
rho_tr_z = rho_tr[zoom_tr]
ord_tr = np.argsort(Ti_tr_z)

# ---------------------------
# Plot
# ---------------------------
fig, ax = plt.subplots(figsize=(9.2, 5.6))

ax.semilogy(T_dense, sv_dense, 'b-', lw=2.2, label='TR built-in table (interp)')
ax.plot(T_tab[(T_tab>=20)&(T_tab<=50)], sv_tab[(T_tab>=20)&(T_tab<=50)],
        'ro', ms=6, label='TR tabulated points')

ax.semilogy(Ti_om_z[ord_om], sv_om_z[ord_om], color='green', lw=1.8, alpha=0.9,
            label='OMFIT inferred (sorted by Ti)')
ax.scatter(Ti_om_z, sv_om_z, c=rho_om_z, cmap='viridis', s=24, alpha=0.75,
           edgecolors='none', label='OMFIT samples')

ax.semilogy(Ti_tr_z[ord_tr], sv_tr_z[ord_tr], color='magenta', lw=1.8, alpha=0.9,
            label='TR inferred from PNF,nD,nT (sorted by T_eff)')
ax.scatter(Ti_tr_z, sv_tr_z, marker='x', color='magenta', s=28, alpha=0.8,
           label='TR inferred samples')

ax.set_xlim(20, 50)
all_y = np.concatenate([sv_dense, sv_om_z if len(sv_om_z)>0 else sv_dense, sv_tr_z if len(sv_tr_z)>0 else sv_dense])
ax.set_ylim(all_y.min()*0.9, all_y.max()*1.12)

ax.set_xlabel('Ion Temperature T [keV]')
ax.set_ylabel(r'DT Reactivity $\langle\sigma v\rangle$ [m$^3$/s]')
ax.set_title('Zoom: DT Reactivity (20-50 keV) with OMFIT and TR inferred')
ax.grid(True, which='both', ls='--', alpha=0.35)
ax.legend(fontsize=9, loc='lower right')

out = '/Users/dengxiaoya/TASK/latest/task/trx/dt_reactivity_zoom_20_50keV_with_omfit_tr_inferred.png'
fig.tight_layout()
fig.savefig(out, dpi=170)
print('Saved:', out)