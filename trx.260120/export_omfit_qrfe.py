#!/usr/bin/env python3
"""
从 OMFIT statefile 提取 qrfe (RF对电子加热)，保存到 CSV
"""

import numpy as np
import pandas as pd
import netCDF4 as nc4

# 路径配置
omfit_dir = '/Users/dengxiaoya/CFEDRSW/OMFIT_out'
tr_dir = '/Users/dengxiaoya/TASK/CFEDR/git/task/trx.260120'

# 小半径 a (m)
a_minor = 2.68

# 读取 OMFIT statefile
statefile_path = f'{omfit_dir}/statefile_3.000000E+01.nc'
ds = nc4.Dataset(statefile_path)

rho_grid = ds.variables['rho_grid'][:]
rho_omfit = rho_grid / rho_grid.max()

# qrfe (单位: W/m^3)
qrfe_wm3 = ds.variables['qrfe'][:]
ds.close()

# 单位转换: W/m^3 -> MW/m^3
qrfe_mwm3 = qrfe_wm3 * 1e-6

# 读取 rmin 并转换为 r/a
df_profiles = pd.read_csv(f'{omfit_dir}/outputs_new/all_profiles.csv')
print
rmin_omfit = df_profiles['rmin'].values  # [m]
r_a_omfit = rmin_omfit / a_minor         # r/a = rmin / 2.68
r_a_omfit = df_profiles['rho'].values        

# 获取 rho_profiles 用于插值
rho_profiles = df_profiles['rho'].values if 'rho' in df_profiles.columns else np.linspace(0, 1, len(rmin_omfit))

# 将 qrfe 插值到 all_profiles 的网格上
qrfe_interp = np.interp(rho_profiles, rho_omfit, qrfe_mwm3)

# 移除 r/a > 1.0 的点（避免边界尖峰问题）
mask = r_a_omfit <= 1.0
r_a_clean = r_a_omfit[mask]
qrfe_clean = qrfe_interp[mask]

# 确保有 r/a = 1.0 的点
if r_a_clean.max() < 1.0:
    r_a_clean = np.append(r_a_clean, 1.0)
    qrfe_clean = np.append(qrfe_clean, qrfe_clean[-1])

# 保存到 CSV
df_output = pd.DataFrame({
    'r_a': r_a_clean,
    'qrfe_MW_m3': qrfe_clean
})

output_path = f'{tr_dir}/omfit_qrfe_for_tr.csv'
df_output.to_csv(output_path, index=False)

print(f"Data saved to: {output_path}")
print(f"\nr/a range: {r_a_clean.min():.4f} - {r_a_clean.max():.4f}")
print(f"qrfe range: {qrfe_clean.min():.6f} - {qrfe_clean.max():.6f} MW/m^3")
print(f"\nFirst 5 rows:")
print(df_output.head())
print(f"\nLast 5 rows:")
print(df_output.tail())

# 计算总功率 (体积积分)
# 简单估算: P_total = sum(qrfe * dV) ≈ sum(qrfe * 4*pi*R*r*dr * kappa)
# 这里只做简单输出
print(f"\nNote: qrfe is RF heating power density to electrons [MW/m³]")
