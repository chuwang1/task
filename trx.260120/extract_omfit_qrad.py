#!/usr/bin/env python3
"""
从 OMFIT statefile 提取 r/a 和 qrad，保存到 CSV
"""

import numpy as np
import pandas as pd
import netCDF4 as nc4

# 路径配置
omfit_dir = '/Users/dengxiaoya/CFEDRSW/OMFIT_out'
tr_dir = '/Users/dengxiaoya/TASK/latest/task/trx'

# 小半径 a (m)
a_minor = 2.68

# 读取 OMFIT statefile
statefile_path = f'{omfit_dir}/statefile_3.000000E+01.nc'
ds = nc4.Dataset(statefile_path)

rho_grid = ds.variables['rho_grid'][:]
rho_omfit = rho_grid / rho_grid.max()

# qrad (单位: W/m^3)
qrad_wm3 = ds.variables['qrad'][:]
ds.close()

# 单位转换: W/m^3 -> MW/m^3
qrad_mwm3 = qrad_wm3 * 1e-6

# 读取 rmin 并转换为 r/a
df_profiles = pd.read_csv(f'{omfit_dir}/outputs_new/all_profiles.csv')
rmin_omfit = df_profiles['rmin'].values  # [m]
r_a_omfit = rmin_omfit / a_minor         # r/a = rmin / 2.68

# 获取 rho_profiles 用于插值
rho_profiles = df_profiles['rho'].values if 'rho' in df_profiles.columns else np.linspace(0, 1, len(rmin_omfit))

# 将 qrad 插值到 all_profiles 的网格上
qrad_interp = np.interp(rho_profiles, rho_omfit, qrad_mwm3)

# 保存到 CSV
df_output = pd.DataFrame({
    'r_a': r_a_omfit,
    'qrad_MW_m3': qrad_interp
})

output_path = f'{tr_dir}/omfit_qrad_extracted.csv'
df_output.to_csv(output_path, index=False)

print(f"Data saved to: {output_path}")
print(f"\nr/a range: {r_a_omfit.min():.4f} - {r_a_omfit.max():.4f}")
print(f"qrad range: {qrad_interp.min():.6f} - {qrad_interp.max():.6f} MW/m^3")
print(f"\nFirst 5 rows:")
print(df_output.head())
print(f"\nLast 5 rows:")
print(df_output.tail())
