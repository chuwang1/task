#!/usr/bin/env python3
"""Export OMFIT qrad profile to CSV for TR to read"""
import numpy as np
import netCDF4 as nc4
import pandas as pd

omfit_dir = '/Users/dengxiaoya/CFEDRSW/OMFIT_out'
statefile_path = f'{omfit_dir}/statefile_3.000000E+01.nc'
ds = nc4.Dataset(statefile_path)

rho_grid = ds.variables['rho_grid'][:]
rho = rho_grid / rho_grid.max()
qrad = ds.variables['qrad'][:]  # W/m^3
ds.close()

# Convert to MW/m^3
qrad_mw = qrad * 1e-6

# Also get r/a from all_profiles
df_profiles = pd.read_csv(f'{omfit_dir}/outputs_new/all_profiles.csv')
rho_profiles = df_profiles['rho'].values if 'rho' in df_profiles.columns else np.linspace(0, 1, len(df_profiles))
rmin = df_profiles['rmin'].values
r_a = rmin / 2.68

# Interpolate qrad to all_profiles grid
qrad_interp = np.interp(rho_profiles, rho, qrad_mw)

# Save as CSV: r/a, qrad [MW/m^3]
out = pd.DataFrame({'r_a': r_a, 'qrad_MW_m3': qrad_interp})
out.to_csv('omfit_qrad_for_tr.csv', index=False, header=True)

print(f"Exported {len(out)} points")
print(f"  r/a range: {r_a.min():.3f} - {r_a.max():.3f}")
print(f"  qrad range: {qrad_interp.min():.4f} - {qrad_interp.max():.4f} MW/m^3")
print(f"  Integrated: ~{np.trapz(qrad_interp, r_a)*2*np.pi:.1f} (rough)")
