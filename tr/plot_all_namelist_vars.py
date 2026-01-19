#!/usr/bin/env python3
"""
Plot ALL profiles from CFETR namelist file with variable explanations
"""

import numpy as np
import matplotlib.pyplot as plt

def parse_custom_namelist(filepath):
    """Parse the custom namelist format"""
    data = {}
    current_section = None
    current_key = None
    current_values = []
    
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()
    
    lines = content.split('\n')
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        
        if not line:
            i += 1
            continue
        
        if line.startswith('&'):
            section_name = line[1:].lower()
            current_section = section_name
            data[current_section] = {}
            i += 1
            continue
        
        if line == '/':
            current_section = None
            i += 1
            continue
        
        if '=' in line:
            if current_key and current_values:
                values_str = ' '.join(current_values)
                parsed = parse_values(values_str)
                if current_section:
                    data[current_section][current_key] = parsed
                else:
                    data[current_key] = parsed
            
            parts = line.split('=', 1)
            current_key = parts[0].strip().lower()
            current_values = [parts[1].strip()]
        else:
            current_values.append(line)
        
        i += 1
    
    if current_key and current_values:
        values_str = ' '.join(current_values)
        parsed = parse_values(values_str)
        if current_section:
            data[current_section][current_key] = parsed
        else:
            data[current_key] = parsed
    
    return data

def parse_values(values_str):
    """Parse a string of values into appropriate Python types"""
    values_str = values_str.strip()
    
    if values_str.startswith("'") or values_str.startswith('"'):
        return values_str.strip("'\"")
    
    if '*' in values_str and not values_str.startswith("'"):
        expanded = []
        for part in values_str.split():
            if '*' in part:
                count, val = part.split('*')
                expanded.extend([float(val)] * int(count))
            else:
                try:
                    expanded.append(float(part))
                except:
                    return values_str
        return np.array(expanded)
    
    try:
        values = []
        for v in values_str.split():
            values.append(float(v))
        if len(values) == 1:
            return values[0]
        return np.array(values)
    except:
        return values_str

# Variable explanations
VARIABLE_DESCRIPTIONS = {
    # Top-level variables
    'time_stamp': '时间戳 - 数据生成时间',
    'iteration_step': '迭代步数 - 模拟收敛所需的迭代次数',
    'rho': 'ρ (归一化磁面径向坐标) - √(ψ_tor/ψ_tor_edge)',
    'psi_normalized': '归一化极向磁通 ψ_N',
    'omega': '环向角速度 ω [rad/s] - 等离子体旋转',
    'zeff': '有效电荷数 Zeff',
    'pressure_total': '总压力 [MPa] - 所有粒子种类压力之和',
    'current_density': '总电流密度 <J_tor·R0/R> [A/m²]',
    'bootstrap_current_density_onetwo': '自举电流密度 [A/m²] - 由压力梯度驱动',
    'nbcd_current_density_onetwo': 'NBI电流驱动密度 [A/m²] - 中性束注入驱动',
    'rfcd_current_density_onetwo': 'RF电流驱动密度 [A/m²] - EC/LH波驱动',
    'ohmic_current_density_onetwo': '欧姆电流密度 [A/m²] - 感应电场驱动',
    'r0': '主半径 R0 [m]',
    'b0': '磁轴环向磁场 B0 [T]',
    'description': '数据描述',
    'unit': '单位说明',
    
    # Electron variables
    'electron.name': '粒子名称',
    'electron.charge': '电荷数 (质子电荷单位)',
    'electron.mass': '质量 (质子质量单位)',
    'electron.density': '电子密度 ne [m⁻³]',
    'electron.temperature': '电子温度 Te [keV]',
    'electron.type': '粒子类型',
    'electron.pressure': '电子压力 pe [MPa]',
    'electron.stored_energy_density': '电子储能密度 [MJ/m³]',
    
    # Ion variables
    'ions_1.name': '离子名称 (D=氘)',
    'ions_1.charge': '离子电荷数',
    'ions_1.mass': '离子质量 (AMU)',
    'ions_1.density': '离子密度 ni [m⁻³]',
    'ions_1.temperature': '离子温度 Ti [keV]',
    'ions_1.type': '离子类型',
    'ions_1.pressure': '离子压力 pi [MPa]',
}

# Parse namelist
nml_path = '/Users/dengxiaoya/COREDIV/CFEDR混合运行模式250925版本/CFEDR混合运行模式250925数据/case250925（即case6.2）/profiles_CFETR.namelist'
print("Parsing namelist file...")
data = parse_custom_namelist(nml_path)

# Print all variables
print("\n" + "="*80)
print("NAMELIST 变量列表及说明")
print("="*80)

# Top-level variables
print("\n【顶层变量】")
for key, value in data.items():
    if not isinstance(value, dict):
        desc = VARIABLE_DESCRIPTIONS.get(key, '(无说明)')
        if isinstance(value, np.ndarray):
            print(f"  {key}: 数组[{len(value)}] - {desc}")
        elif isinstance(value, str) and len(value) > 50:
            print(f"  {key}: 字符串 - {desc}")
        else:
            print(f"  {key}: {value} - {desc}")

# Section variables
for section_name, section_data in data.items():
    if isinstance(section_data, dict):
        print(f"\n【{section_name.upper()} 组】")
        for key, value in section_data.items():
            full_key = f"{section_name}.{key}"
            desc = VARIABLE_DESCRIPTIONS.get(full_key, '(无说明)')
            if isinstance(value, np.ndarray):
                print(f"  {key}: 数组[{len(value)}] - {desc}")
            elif isinstance(value, str) and len(value) > 50:
                print(f"  {key}: 字符串 - {desc}")
            else:
                print(f"  {key}: {value} - {desc}")

# Get rho for x-axis
rho = data['rho']

# Count plottable arrays
profile_vars = []
for key, value in data.items():
    if isinstance(value, np.ndarray) and len(value) == len(rho):
        profile_vars.append((key, value, key))
    elif isinstance(value, dict):
        for subkey, subvalue in value.items():
            if isinstance(subvalue, np.ndarray) and len(subvalue) == len(rho):
                profile_vars.append((f"{key}.{subkey}", subvalue, f"{key}:{subkey}"))

print(f"\n共有 {len(profile_vars)} 个可绘制的剖面变量")

# Create comprehensive plot
n_plots = len(profile_vars)
n_cols = 4
n_rows = (n_plots + n_cols - 1) // n_cols

fig, axes = plt.subplots(n_rows, n_cols, figsize=(20, 4*n_rows))
axes = axes.flatten()

for idx, (name, values, label) in enumerate(profile_vars):
    ax = axes[idx]
    ax.plot(rho, values, 'b-', linewidth=1.5)
    ax.set_xlabel('r/a')
    ax.set_title(label, fontsize=9)
    ax.grid(True, alpha=0.3)
    
    # Format y-axis for large numbers
    if np.max(np.abs(values)) > 1e10:
        ax.ticklabel_format(style='scientific', axis='y', scilimits=(0,0))

# Hide unused subplots
for idx in range(len(profile_vars), len(axes)):
    axes[idx].set_visible(False)

plt.tight_layout()
output_path = '/Users/dengxiaoya/TASK/latest/task/tr/namelist_all_profiles.png'
plt.savefig(output_path, dpi=120, bbox_inches='tight')
print(f"\n图片已保存到: {output_path}")

plt.show()
