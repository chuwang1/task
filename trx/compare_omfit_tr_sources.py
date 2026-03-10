#!/usr/bin/env python3
"""
对比 OMFIT statefile 热源与 TR 代码输出的对应变量

OMFIT -> TR 对应关系:
  电子:
    qfuse  (聚变对电子加热)    -> PNFCL_e (tr_data_028.csv)
    qdelt  (离子-电子能量交换)  -> PIE     (tr_data_022.csv)
    qrad   (辐射损失)          -> PRSUM   (tr_data_022.csv)
    qrfe   (RF加热电子)        -> PRF_e   (tr_data_021.csv)

  离子:
    qfusi  (聚变对离子加热)    -> PNFCL_i (tr_data_028.csv)
    qcx    (电荷交换损失)      -> PCX     (tr_data_022.csv)
"""

import numpy as np
import pandas as pd
import netCDF4 as nc4
import matplotlib
matplotlib.use('Agg')  # 非交互式后端，不弹出窗口
import matplotlib.pyplot as plt

# ============================================================================
# 路径配置
# ============================================================================
omfit_dir = '/Users/dengxiaoya/CFEDRSW/OMFIT_out'
tr_dir = '/Users/dengxiaoya/TASK/latest/task/trx'

# 小半径 a (m)
a_minor = 2.68

# ============================================================================
# 1. 读取 TR 数据
# ============================================================================
import glob
import os

def find_csv_by_title(title_keyword, search_dir=tr_dir):
    """通过标题行内容查找 TR CSV 文件"""
    for f in sorted(glob.glob(os.path.join(search_dir, 'tr_data_*.csv'))):
        with open(f, 'r') as fh:
            title = fh.readline().strip()
        if title_keyword in title:
            return os.path.basename(f)
    return None

def read_tr_csv(filename):
    """读取TR输出的CSV文件，跳过标题行"""
    df = pd.read_csv(f'{tr_dir}/{filename}', skiprows=1)
    # 清理列名中的空格
    df.columns = df.columns.str.strip()
    return df

# 通过标题内容自动查找CSV文件（不依赖固定编号）
csv_021 = find_csv_by_title('POH,PNB,PNF,-PRSUM,PRF')  # power vs r
csv_022 = find_csv_by_title('PRSUM,PRB,PRC,PRL,PCX,PIE')  # radiation/exchange vs r
csv_028 = find_csv_by_title('PNFIN,PNFCL')  # fusion power vs r

print(f"CSV file mapping:")
print(f"  Power sources  -> {csv_021}")
print(f"  Radiation/exch -> {csv_022}")
print(f"  Fusion detail  -> {csv_028}")

missing = []
if csv_021 is None: missing.append('POH/PNB/PNF/PRF')
if csv_022 is None: missing.append('PRSUM/PRB/PCX/PIE')
if csv_028 is None: missing.append('PNFIN/PNFCL')
if missing:
    raise FileNotFoundError(f"Could not find CSV files for: {missing}")

# 读取各个TR CSV文件
df_021 = read_tr_csv(csv_021)
df_022 = read_tr_csv(csv_022)
df_028 = read_tr_csv(csv_028)

# TR 径向网格对齐说明：
# tr_data_*.csv 的 X 为 RM=(nr-0.5)*dr；
# 但外部 qrfe 在 TR 中是按 RG=nr*dr 读入并插值（见 tr_prep_prffixed）。
# 为避免与外部 qrfe 对比时出现半网格偏移，统一使用 RG 对齐。
r_a_tr_rm = df_022['X'].values
dr_tr = float(np.median(np.diff(r_a_tr_rm))) if len(r_a_tr_rm) > 1 else 0.02
r_a_tr = np.clip(r_a_tr_rm + 0.5 * dr_tr, 0.0, 1.0)  # RG grid

# TR 变量 (单位: MW/m^3)
TR = {
    'PRSUM': df_022['PRSUM'].values,         # 辐射损失
    'PIE': df_022['PIE'].values,             # 离子-电子能量交换
    'PCX': df_022['PCX'].values,             # 电荷交换损失
    'PNF': df_021['PNF'].values,             # 聚变总功率
    'PNFCL_e': df_028.iloc[:, 2].values,     # 聚变对电子加热 (第3列)
    'PNFCL_D': df_028.iloc[:, 3].values,     # 聚变对D加热 (第4列)
    'PNFCL_T': df_028.iloc[:, 4].values,     # 聚变对T加热 (第5列)
    'PNFCL_He': df_028.iloc[:, 5].values if df_028.shape[1] > 5 else np.zeros_like(r_a_tr),
    'PRF_e': df_021.iloc[:, 5].values,       # RF对电子加热 (PRF第一列,NS=1)
}

# 聚变对离子加热 = D + T + He
TR['PNFCL_i'] = TR['PNFCL_D'] + TR['PNFCL_T'] + TR['PNFCL_He']

print("TR data loaded:")
print(f"  r/a range (RM from CSV): {r_a_tr_rm.min():.3f} - {r_a_tr_rm.max():.3f}")
print(f"  r/a range (RG for compare): {r_a_tr.min():.3f} - {r_a_tr.max():.3f}")
for k, v in TR.items():
    print(f"  {k}: max = {v.max():.4f} MW/m^3")

# ============================================================================
# 2. 读取 OMFIT statefile
# ============================================================================
statefile_path = f'{omfit_dir}/statefile_3.000000E+01.nc'
ds = nc4.Dataset(statefile_path)

rho_grid = ds.variables['rho_grid'][:]
rho_omfit = rho_grid / rho_grid.max()

# OMFIT 变量 (单位: W/m^3)
OMFIT = {
    'qfuse': ds.variables['qfuse'][:],      # 聚变对电子加热
    'qfusi': ds.variables['qfusi'][:],      # 聚变对离子加热
    'qdelt': ds.variables['qdelt'][:],      # 离子-电子能量交换
    'qrad': ds.variables['qrad'][:],        # 辐射损失
    'qrfe': ds.variables['qrfe'][:],        # RF对电子加热
    'qcx': ds.variables['qcx'][:],          # 电荷交换损失
}

ds.close()

print("\nOMFIT statefile loaded:")
print(f"  rho range: {rho_omfit.min():.3f} - {rho_omfit.max():.3f}")
for k, v in OMFIT.items():
    print(f"  {k}: max = {v.max():.4e} W/m^3 = {v.max()*1e-6:.4f} MW/m^3")

# ============================================================================
# 3. 单位转换: OMFIT W/m^3 -> MW/m^3
# ============================================================================
for k in OMFIT:
    OMFIT[k] = OMFIT[k] * 1e-6  # W/m^3 -> MW/m^3

# ============================================================================
# 4. 读取 rmin 并转换为 r/a
# ============================================================================
# 从 all_profiles.csv 读取 rmin
df_profiles = pd.read_csv(f'{omfit_dir}/outputs_new/all_profiles.csv')
# rmin_omfit = df_profiles['rmin'].values  # [m]
r_a_omfit = df_profiles['rho']         # r/a = rmin / 2.68

print(f"\nOMFIT r/a (from rmin/{a_minor}):")
print(f"  r/a range: {r_a_omfit.min():.3f} - {r_a_omfit.max():.3f}")

# 将 statefile 变量插值到 all_profiles 的网格上
# statefile 用 rho_omfit, all_profiles 用 rho_profiles
rho_profiles = df_profiles['rho'].values if 'rho' in df_profiles.columns else np.linspace(0, 1, len(rmin_omfit))

for k in OMFIT:
    OMFIT[k] = np.interp(rho_profiles, rho_omfit, OMFIT[k])

# ============================================================================
# 5. 绘图对比
# ============================================================================
fig, axes = plt.subplots(2, 3, figsize=(15, 10))

comparisons = [
    # (OMFIT_key, TR_key, title, sign_omfit, sign_tr)
    ('qfuse', 'PNFCL_e', 'Fusion -> Electron', 1, 1),
    ('qfusi', 'PNFCL_i', 'Fusion -> Ion', 1, 1),
    ('qdelt', 'PIE', 'Ion-Electron Exchange (qdelt/PIE)', 1, 1),
    ('qrad', 'PRSUM', 'Radiation Loss (qrad/PRSUM)', 1, 1),
    ('qrfe', 'PRF_e', 'RF -> Electron', 1, 1),
    ('qcx', 'PCX', 'Charge Exchange Loss', 1, 1),
]

for ax, (omfit_key, tr_key, title, sign_o, sign_t) in zip(axes.flatten(), comparisons):
    # OMFIT
    ax.plot(r_a_omfit, sign_o * OMFIT[omfit_key], 'b-', lw=2, label=f'OMFIT: {omfit_key}')
    # TR
    ax.plot(r_a_tr, sign_t * TR[tr_key], 'r--', lw=2, label=f'TR: {tr_key}')

    ax.set_xlabel('r/a', fontsize=12)
    ax.set_ylabel('MW/m$^3$', fontsize=12)
    ax.set_title(title, fontsize=13)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 1)

fig.suptitle('OMFIT vs TR: Heat Source Comparison', fontsize=14)
fig.tight_layout()
fig.savefig(f'{tr_dir}/compare_omfit_tr_sources.png', dpi=150, bbox_inches='tight')
print(f"\nFigure saved: {tr_dir}/compare_omfit_tr_sources.png")

# ============================================================================
# 6. 保存对比数据
# ============================================================================
# 插值 OMFIT 到 TR 的 r/a 网格
compare_data = {'r_a': r_a_tr}
for omfit_key, tr_key, _, sign_o, sign_t in comparisons:
    omfit_interp = np.interp(r_a_tr, r_a_omfit, sign_o * OMFIT[omfit_key])
    compare_data[f'OMFIT_{omfit_key}'] = omfit_interp
    compare_data[f'TR_{tr_key}'] = sign_t * TR[tr_key]
    compare_data[f'diff_{omfit_key}'] = omfit_interp - sign_t * TR[tr_key]

df_compare = pd.DataFrame(compare_data)
csv_path = f'{tr_dir}/compare_omfit_tr_sources.csv'
df_compare.to_csv(csv_path, index=False)
print(f"Comparison data saved: {csv_path}")
