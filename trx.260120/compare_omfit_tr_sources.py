#!/usr/bin/env python3
"""
对比 OMFIT statefile 热源与 TR 代码输出的对应变量

OMFIT -> TR 对应关系:
  电子:
    qfuse  (聚变对电子加热)    -> PNFCL_e (tr_data_028.csv)
    qdelt  (离子-电子能量交换)  -> QEI     (tr_data_022.csv)
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
tr_dir = '/Users/dengxiaoya/TASK/CFEDR/git/task/trx.260120'

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
csv_021 = find_csv_by_title('POH,PNB,PNF')                  # power vs r (matches both PRSUM and -PRSUM variants)
csv_022 = find_csv_by_title('PRSUM,PRB,PRC,PRL,PCX,PIE')   # radiation/exchange vs r
csv_027 = find_csv_by_title('PNBIN,PNBCL')                 # NBI collisional split vs r
csv_028 = find_csv_by_title('PNFIN,PNFCL')                 # fusion split vs r
csv_037 = find_csv_by_title('@PIN [MW/m$+3$=]')            # PIN by species vs r

print(f"CSV file mapping:")
print(f"  Power sources  -> {csv_021}")
print(f"  Radiation/exch -> {csv_022}")
print(f"  NBI split      -> {csv_027}")
print(f"  Fusion detail  -> {csv_028}")
print(f"  PIN species    -> {csv_037}")

missing = []
if csv_021 is None: missing.append('POH/PNB/PNF/PRSUM/PRF')
if csv_022 is None: missing.append('PRSUM/PRB/PCX/PIE')
if csv_027 is None: missing.append('PNBIN/PNBCL')
if csv_028 is None: missing.append('PNFIN/PNFCL')
if csv_037 is None: missing.append('PIN species')
if missing:
    raise FileNotFoundError(f"Could not find CSV files for: {missing}")

# 读取各个TR CSV文件
df_021 = read_tr_csv(csv_021)
df_022 = read_tr_csv(csv_022)
df_027 = read_tr_csv(csv_027)
df_028 = read_tr_csv(csv_028)
df_037 = read_tr_csv(csv_037)

# TR 径向网格对齐说明：
# tr_data_*.csv 的 X 为 RM=(nr-0.5)*dr；
# 但外部 qrfe 在 TR 中是按 RG=nr*dr 读入并插值（见 tr_prep_prffixed）。
# 为避免与外部 qrfe 对比时出现半网格偏移，统一使用 RG 对齐。
r_a_tr_rm = df_022['X'].values
dr_tr = float(np.median(np.diff(r_a_tr_rm))) if len(r_a_tr_rm) > 1 else 0.02
r_a_tr = np.clip(r_a_tr_rm + 0.5 * dr_tr, 0.0, 1.0)  # RG grid

# TR 变量 (单位: MW/m^3)
TR = {
    'POH': df_021['POH'].values,             # 欧姆加热
    'PRSUM': df_022['PRSUM'].values,         # 辐射损失
    'QEI': df_022['QEI'].values if 'QEI' in df_022.columns else df_022['PIE'].values,  # 电子-离子能量交换(显式诊断)
    'PIE': df_022['PIE'].values,             # 电离损失
    'PCX': df_022['PCX'].values,             # 电荷交换损失
    'PNF': df_021['PNF'].values,             # 聚变总功率
    'PNBCL_e': df_027.iloc[:, 2].values,     # NBI对电子加热
    'PNBCL_D': df_027.iloc[:, 3].values,
    'PNBCL_T': df_027.iloc[:, 4].values,
    'PNBCL_He': df_027.iloc[:, 5].values if df_027.shape[1] > 5 else np.zeros_like(r_a_tr),
    'PNFCL_e': df_028.iloc[:, 2].values,     # 聚变对电子加热 (第3列)
    'PNFCL_D': df_028.iloc[:, 3].values,     # 聚变对D加热 (第4列)
    'PNFCL_T': df_028.iloc[:, 4].values,     # 聚变对T加热 (第5列)
    'PNFCL_He': df_028.iloc[:, 5].values if df_028.shape[1] > 5 else np.zeros_like(r_a_tr),
    'PRF_e': df_021.iloc[:, 5].values,       # RF对电子加热 (PRF第一列,NS=1)
    'PRF_D': df_021.iloc[:, 6].values if df_021.shape[1] > 6 else np.zeros_like(r_a_tr),
    'PRF_T': df_021.iloc[:, 7].values if df_021.shape[1] > 7 else np.zeros_like(r_a_tr),
    'PRF_He': df_021.iloc[:, 8].values if df_021.shape[1] > 8 else np.zeros_like(r_a_tr),
    'PIN_e': df_037['PIN_1'].values,
    'PIN_D': df_037['PIN_2'].values if 'PIN_2' in df_037.columns else np.zeros_like(r_a_tr),
    'PIN_T': df_037['PIN_3'].values if 'PIN_3' in df_037.columns else np.zeros_like(r_a_tr),
    'PIN_He': df_037['PIN_4'].values if 'PIN_4' in df_037.columns else np.zeros_like(r_a_tr),
}

# 聚变对离子加热 = D + T + He
TR['PNFCL_i'] = TR['PNFCL_D'] + TR['PNFCL_T'] + TR['PNFCL_He']
TR['PNBCL_i'] = TR['PNBCL_D'] + TR['PNBCL_T'] + TR['PNBCL_He']
TR['PRF_i'] = TR['PRF_D'] + TR['PRF_T'] + TR['PRF_He']
TR['PIN_i'] = TR['PIN_D'] + TR['PIN_T'] + TR['PIN_He']

# TR净源项：
# - PIN口径（对应TR显式方程）
# - statefile-like口径（额外加入QEI诊断，便于与OMFIT qdelt/qdelt_i比较）
TR['Se_net_PIN'] = TR['PIN_e']
TR['Se_net_state_like'] = TR['PNBCL_e'] + TR['PNFCL_e'] + TR['PRF_e'] + TR['POH'] - TR['PRSUM'] - TR['PIE'] + TR['QEI']
TR['Si_net_PIN'] = TR['PIN_i']
TR['Si_net_state_like'] = TR['PNBCL_i'] + TR['PNFCL_i'] + TR['PRF_i'] - TR['PCX'] - TR['QEI']

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
def get_var(name):
    if name in ds.variables:
        return ds.variables[name][:]
    return np.zeros_like(rho_omfit)

OMFIT = {
    'qfuse': get_var('qfuse'),      # 聚变对电子加热
    'qfusi': get_var('qfusi'),      # 聚变对离子加热
    'qdelt': get_var('qdelt'),      # 离子-电子交换(电子侧)
    'qdelt_i': get_var('qdelt_i'),  # 离子-电子交换(离子侧)
    'qrad': get_var('qrad'),        # 辐射损失
    'qrfe': get_var('qrfe'),        # RF对电子加热
    'qrfi': get_var('qrfi'),        # RF对离子加热
    'qbeame': get_var('qbeame'),    # NBI对电子加热
    'qbeami': get_var('qbeami'),    # NBI对离子加热
    'qe2d': get_var('qe2d'),        # 电子原子过程项
    'qioni': get_var('qioni'),      # 离子原子过程项
    'qmag': get_var('qmag'),        # 欧姆/磁功率项
    'qcx': get_var('qcx'),          # 电荷交换损失
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
# 4. 使用 statefile 原始 rho 网格（不再插值到 all_profiles 网格）
# ============================================================================
r_a_omfit = rho_omfit

print(f"\nOMFIT r/a (statefile rho_grid):")
print(f"  r/a range: {r_a_omfit.min():.3f} - {r_a_omfit.max():.3f}")

# ============================================================================
# 5. 构造电子/离子净源项及分项（TR vs OMFIT）
#    OMFIT口径参考 compare_Se_statefile.py 中 statefile sum
# ============================================================================
TR_e_components = {
    'Fusion': TR['PNFCL_e'],
    'RF': TR['PRF_e'],
    'NBI': TR['PNBCL_e'],
    'Ohmic': TR['POH'],
    'Atomic(-PIE)': -TR['PIE'],
    'e-i exchange(QEI)': TR['QEI'],
    'Radiation(-PRSUM)': -TR['PRSUM'],
}

TR_i_components = {
    'Fusion': TR['PNFCL_i'],
    'RF': TR['PRF_i'],
    'NBI': TR['PNBCL_i'],
    'Atomic(PEX~0)': np.zeros_like(r_a_tr),
    'i-e exchange(-QEI)': -TR['QEI'],
    'CX(-PCX)': -TR['PCX'],
}

OMFIT_e_components = {
    'Fusion': OMFIT['qfuse'],
    'RF': OMFIT['qrfe'],
    'NBI': OMFIT['qbeame'],
    'Ohmic': OMFIT['qmag'],
    'Atomic(qe2d)': OMFIT['qe2d'],
    'e-i exchange(qdelt)': OMFIT['qdelt'],
    'Radiation(-qrad)': -OMFIT['qrad'],
}

OMFIT_i_components = {
    'Fusion': OMFIT['qfusi'],
    'RF': OMFIT['qrfi'],
    'NBI': OMFIT['qbeami'],
    'Atomic(qioni)': OMFIT['qioni'],
    'i-e exchange(qdelt_i)': OMFIT['qdelt_i'],
    'CX(-qcx)': -OMFIT['qcx'],
}

TR['Se_net_components_sum'] = np.sum(np.vstack(list(TR_e_components.values())), axis=0)
TR['Si_net_components_sum'] = np.sum(np.vstack(list(TR_i_components.values())), axis=0)
OMFIT['Se_net_components_sum'] = np.sum(np.vstack(list(OMFIT_e_components.values())), axis=0)
OMFIT['Si_net_components_sum'] = np.sum(np.vstack(list(OMFIT_i_components.values())), axis=0)

# ============================================================================
# 5. 绘图对比（含电子/离子净源项）
# ============================================================================
fig, axes = plt.subplots(3, 3, figsize=(18, 14))

comparisons = [
    # (OMFIT_key, TR_key, title, sign_omfit, sign_tr)
    ('qfuse', 'PNFCL_e', 'Fusion -> Electron', 1, 1),
    ('qfusi', 'PNFCL_i', 'Fusion -> Ion', 1, 1.),
    ('qdelt', 'QEI', 'Ion-Electron Exchange (qdelt/QEI)', 1, 1),
    ('qrad', 'PRSUM', 'Radiation Loss (qrad/PRSUM)', 1, 1),
    ('qrfe', 'PRF_e', 'RF -> Electron', 1, 1),
    ('qrfi', 'PRF_i', 'RF -> Ion', 1, 1),
    ('qcx', 'PCX', 'Charge Exchange Loss', 1, 1),
    ('Se_net_components_sum', 'Se_net_state_like', 'Net Source -> Electron', 1, 1),
    ('Si_net_components_sum', 'Si_net_state_like', 'Net Source -> Ion', 1, 1),
]

for ax, (omfit_key, tr_key, title, sign_o, sign_t) in zip(axes.flatten(), comparisons):
    # OMFIT

    ax.plot(r_a_omfit, sign_o * OMFIT[omfit_key], 'r-', lw=2, label=f'OMFIT: {omfit_key}')
    # TR
    ax.plot(r_a_tr, sign_t * TR[tr_key], 'b--', lw=2, label=f'TR: {tr_key}')

    ax.set_xlabel('rho', fontsize=12)
    ax.set_ylabel('MW/m$^3$', fontsize=12)
    ax.set_title(title, fontsize=13)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 1)
    if omfit_key == 'qrfi':
        ax.set_ylim(-0.03, 0.03)
    elif omfit_key == 'qcx':
        ax.set_ylim(-0.25, 0.6)
    else:
        ax.set_ylim(bottom=-0.25)

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

# ============================================================================
# 7. 电子净源项对比图（总和 + 分项）
# ============================================================================
fig_e, axes_e = plt.subplots(2, 1, figsize=(12, 10), sharex=True)
style_omfit = '-'
style_tr = '--'

# 总净源项
axes_e[0].plot(
    r_a_omfit,
    OMFIT['Se_net_components_sum'],
    style_omfit,
    color='tab:blue',
    lw=2.2,
    label='Electron net source (OMFIT)'
)
axes_e[0].plot(
    r_a_tr,
    TR['Se_net_components_sum'],
    style_tr,
    color='tab:blue',
    lw=2.2,
    label='Electron net source (TR)'
)
axes_e[0].plot(
    r_a_tr,
    TR['Se_net_PIN'],
    ':',
    color='tab:gray',
    lw=1.8,
    label='TR PIN_e (reference)'
)
axes_e[0].set_ylabel('MW/m$^3$', fontsize=12)
axes_e[0].set_title('Electron Net Source: Total Comparison', fontsize=13)
axes_e[0].grid(True, alpha=0.3)
axes_e[0].legend(fontsize=9)
axes_e[0].set_xlim(0, 1)
axes_e[0].set_ylim(bottom=-0.3)

# 分项
e_component_pairs = [
    ('Fusion', OMFIT_e_components['Fusion'], TR_e_components['Fusion']),
    ('RF', OMFIT_e_components['RF'], TR_e_components['RF']),
    ('NBI', OMFIT_e_components['NBI'], TR_e_components['NBI']),
    ('Ohmic/qmag', OMFIT_e_components['Ohmic'], TR_e_components['Ohmic']),
    ('Atomic (qe2d vs -PIE)', OMFIT_e_components['Atomic(qe2d)'], TR_e_components['Atomic(-PIE)']),
    ('e-i exchange', OMFIT_e_components['e-i exchange(qdelt)'], TR_e_components['e-i exchange(QEI)']),
    ('Radiation loss', OMFIT_e_components['Radiation(-qrad)'], TR_e_components['Radiation(-PRSUM)']),
]

e_colors = {
    'Fusion': 'tab:blue',
    'RF': 'tab:orange',
    'NBI': 'tab:green',
    'Ohmic/qmag': 'tab:red',
    'Atomic (qe2d vs -PIE)': 'tab:purple',
    'e-i exchange': 'tab:brown',
    'Radiation loss': 'tab:cyan',
}

for label, omfit_v, tr_v in e_component_pairs:
    c = e_colors[label]
    axes_e[1].plot(r_a_omfit, omfit_v, style_omfit, color=c, lw=1.8, label=f'{label} (OMFIT)')
    axes_e[1].plot(r_a_tr, tr_v, style_tr, color=c, lw=1.6, label=f'{label} (TR)')

axes_e[1].set_xlabel('rho', fontsize=12)
axes_e[1].set_ylabel('MW/m$^3$', fontsize=12)
axes_e[1].set_title('Electron Net Source: Component Breakdown', fontsize=13)
axes_e[1].grid(True, alpha=0.3)
axes_e[1].legend(fontsize=8, ncol=2)
axes_e[1].set_xlim(0, 1)
axes_e[1].set_ylim(bottom=-0.3)

fig_e.tight_layout()
fig_e_path = f'{tr_dir}/compare_omfit_tr_net_source_electron.png'
fig_e.savefig(fig_e_path, dpi=150, bbox_inches='tight')
print(f'Figure saved: {fig_e_path}')

# ============================================================================
# 8. 离子净源项对比图（总和 + 分项）
# ============================================================================
fig_i, axes_i = plt.subplots(2, 1, figsize=(12, 10), sharex=True)

# 总净源项
axes_i[0].plot(
    r_a_omfit,
    OMFIT['Si_net_components_sum'],
    style_omfit,
    color='tab:blue',
    lw=2.2,
    label='Ion net source (OMFIT)'
)
axes_i[0].plot(
    r_a_tr,
    TR['Si_net_components_sum'],
    style_tr,
    color='tab:blue',
    lw=2.2,
    label='Ion net source (TR)'
)
axes_i[0].plot(
    r_a_tr,
    TR['Si_net_PIN'],
    ':',
    color='tab:gray',
    lw=1.8,
    label='TR PIN_i (reference)'
)
axes_i[0].set_ylabel('MW/m$^3$', fontsize=12)
axes_i[0].set_title('Ion Net Source: Total Comparison', fontsize=13)
axes_i[0].grid(True, alpha=0.3)
axes_i[0].legend(fontsize=9)
axes_i[0].set_xlim(0, 1)
axes_i[0].set_ylim(bottom=-0.3)

# 分项
i_component_pairs = [
    ('Fusion', OMFIT_i_components['Fusion'], TR_i_components['Fusion']),
    ('RF', OMFIT_i_components['RF'], TR_i_components['RF']),
    ('NBI', OMFIT_i_components['NBI'], TR_i_components['NBI']),
    ('Atomic (qioni vs ~0)', OMFIT_i_components['Atomic(qioni)'], TR_i_components['Atomic(PEX~0)']),
    ('i-e exchange', OMFIT_i_components['i-e exchange(qdelt_i)'], TR_i_components['i-e exchange(-QEI)']),
    ('CX loss', OMFIT_i_components['CX(-qcx)'], TR_i_components['CX(-PCX)']),
]

i_colors = {
    'Fusion': 'tab:blue',
    'RF': 'tab:orange',
    'NBI': 'tab:green',
    'Atomic (qioni vs ~0)': 'tab:red',
    'i-e exchange': 'tab:purple',
    'CX loss': 'tab:brown',
}

for label, omfit_v, tr_v in i_component_pairs:
    c = i_colors[label]
    axes_i[1].plot(r_a_omfit, omfit_v, style_omfit, color=c, lw=1.8, label=f'{label} (OMFIT)')
    axes_i[1].plot(r_a_tr, tr_v, style_tr, color=c, lw=1.6, label=f'{label} (TR)')

axes_i[1].set_xlabel('rho', fontsize=12)
axes_i[1].set_ylabel('MW/m$^3$', fontsize=12)
axes_i[1].set_title('Ion Net Source: Component Breakdown', fontsize=13)
axes_i[1].grid(True, alpha=0.3)
axes_i[1].legend(fontsize=8, ncol=2)
axes_i[1].set_xlim(0, 1)
axes_i[1].set_ylim(bottom=-0.3)

fig_i.tight_layout()
fig_i_path = f'{tr_dir}/compare_omfit_tr_net_source_ion.png'
fig_i.savefig(fig_i_path, dpi=150, bbox_inches='tight')
print(f'Figure saved: {fig_i_path}')

# ============================================================================
# 9. 保存净源项对比数据（插值到TR网格）
# ============================================================================
net_compare = {'r_a': r_a_tr}

net_compare['OMFIT_Se_net_sum'] = np.interp(r_a_tr, r_a_omfit, OMFIT['Se_net_components_sum'])
net_compare['TR_Se_net_sum'] = TR['Se_net_components_sum']
net_compare['TR_Se_PIN'] = TR['Se_net_PIN']

net_compare['OMFIT_Si_net_sum'] = np.interp(r_a_tr, r_a_omfit, OMFIT['Si_net_components_sum'])
net_compare['TR_Si_net_sum'] = TR['Si_net_components_sum']
net_compare['TR_Si_PIN'] = TR['Si_net_PIN']

for key, val in OMFIT_e_components.items():
    net_compare[f'OMFIT_e_{key}'] = np.interp(r_a_tr, r_a_omfit, val)
for key, val in TR_e_components.items():
    net_compare[f'TR_e_{key}'] = val

for key, val in OMFIT_i_components.items():
    net_compare[f'OMFIT_i_{key}'] = np.interp(r_a_tr, r_a_omfit, val)
for key, val in TR_i_components.items():
    net_compare[f'TR_i_{key}'] = val

net_csv_path = f'{tr_dir}/compare_omfit_tr_net_sources.csv'
pd.DataFrame(net_compare).to_csv(net_csv_path, index=False)
print(f'Net source comparison data saved: {net_csv_path}')
