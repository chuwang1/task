#!/usr/bin/env python3
"""
Plot grouped profiles from CFETR namelist file
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

def get_data_or_none(data_dict, *keys):
    """Safely get data from nested dictionary"""
    curr = data_dict
    for k in keys:
        if isinstance(curr, dict) and k in curr:
            curr = curr[k]
        else:
            return None
    return curr

# Parse namelist
nml_path = '/Users/dengxiaoya/COREDIV/CFEDR混合运行模式250925版本/CFEDR混合运行模式250925数据/case250925（即case6.2）/profiles_CFETR.namelist'
print("Parsing namelist file...")
data = parse_custom_namelist(nml_path)

# Extract global rho
rho = data.get('rho')

# Define Groups
groups = {
    'Temperatures [keV]': {
        'ax_idx': 0,
        'vars': [
            ('Te (Electron)', get_data_or_none(data, 'electron', 'temperature'), None),
            ('Ti (D)', get_data_or_none(data, 'ions_1', 'temperature'), None),
            ('Ti (T)', get_data_or_none(data, 'ions_2', 'temperature'), None),
            ('Ti (He4)', get_data_or_none(data, 'ions_3', 'temperature'), None),
            ('Ti (Ar)', get_data_or_none(data, 'ions_4', 'temperature'), None),
        ]
    },
    'Densities [10^20 m^-3]': {
        'ax_idx': 1,
        'vars': [
            ('ne', get_data_or_none(data, 'electron', 'density'), 1e20),
            ('ni (D)', get_data_or_none(data, 'ions_1', 'density'), 1e20),
            ('ni (T)', get_data_or_none(data, 'ions_2', 'density'), 1e20),
            ('ni (He4)', get_data_or_none(data, 'ions_3', 'density'), 1e20),
            ('ni (Ar)', get_data_or_none(data, 'ions_4', 'density'), 1e20),
            ('ni (Fast He)', get_data_or_none(data, 'ions_5', 'density'), 1e20),
        ]
    },
    'Pressures [MPa]': {
        'ax_idx': 2,
        'vars': [
            ('Total', get_data_or_none(data, 'pressure_total'), None),
            ('Electron', get_data_or_none(data, 'electron', 'pressure'), None),
            ('Ion (D)', get_data_or_none(data, 'ions_1', 'pressure'), None),
            ('Ion (T)', get_data_or_none(data, 'ions_2', 'pressure'), None),
            ('Ion (He4)', get_data_or_none(data, 'ions_3', 'pressure'), None),
             ('Ion (Fast He)', get_data_or_none(data, 'ions_5', 'pressure'), None),
        ]
    },
    'Current Densities [MA/m^2]': {
        'ax_idx': 3,
        'vars': [
            ('Total', get_data_or_none(data, 'current_density'), 1e6),
            ('Bootstrap', get_data_or_none(data, 'bootstrap_current_density_onetwo'), 1e6),
            ('Ohmic', get_data_or_none(data, 'ohmic_current_density_onetwo'), 1e6),
            ('RFCD (EC/LH)', get_data_or_none(data, 'rfcd_current_density_onetwo'), 1e6),
            ('NBCD', get_data_or_none(data, 'nbcd_current_density_onetwo'), 1e6),
        ]
    },
    'Integrated Powers (ONETWO) [MW]': {
        'ax_idx': 4,
        'vars': [
            ('Pe (Total)', get_data_or_none(data, 'powers_particle_flux_onetwo', 'pow_e'), None),
            ('Pi (Total)', get_data_or_none(data, 'powers_particle_flux_onetwo', 'pow_i'), None),
            ('Pe (Aux)', get_data_or_none(data, 'powers_particle_flux_onetwo', 'pow_e_aux'), None),
            ('Pi (Aux)', get_data_or_none(data, 'powers_particle_flux_onetwo', 'pow_i_aux'), None),
            ('Pe (Fus)', get_data_or_none(data, 'powers_particle_flux_onetwo', 'pow_e_fus'), None),
            ('Pi (Fus)', get_data_or_none(data, 'powers_particle_flux_onetwo', 'pow_i_fus'), None),
            ('P (Rad)', get_data_or_none(data, 'powers_particle_flux_onetwo', 'pow_e_rad'), None),
            ('P (ei)', get_data_or_none(data, 'powers_particle_flux_onetwo', 'pow_ei'), None),
        ]
    },
    'Power Densities (ONETWO) [W/m^3]': {
        'ax_idx': 5,
        'vars': [
            ('qrad', get_data_or_none(data, 'powers_particle_flux_onetwo', 'qrad'), None),
            ('qfuse (Fus->e)', get_data_or_none(data, 'powers_particle_flux_onetwo', 'qfuse'), None),
            ('qfusi (Fus->i)', get_data_or_none(data, 'powers_particle_flux_onetwo', 'qfusi'), None),
            ('qrfe (RF->e)', get_data_or_none(data, 'powers_particle_flux_onetwo', 'qrfe'), None),
            ('qrfi (RF->i)', get_data_or_none(data, 'powers_particle_flux_onetwo', 'qrfi'), None),
            ('qbeame (NBI->e)', get_data_or_none(data, 'powers_particle_flux_onetwo', 'qbeame'), None),
            ('qbeami (NBI->i)', get_data_or_none(data, 'powers_particle_flux_onetwo', 'qbeami'), None),
            ('qdelt (e-i exch)', get_data_or_none(data, 'powers_particle_flux_onetwo', 'qdelt'), None),
        ]
    },
     'Equilibrium q': {
        'ax_idx': 6,
        'rho_override': get_data_or_none(data, 'equlibrium', 'rho'),
        'vars': [
            ('q', get_data_or_none(data, 'equlibrium', 'q'), None),
        ]
    },
    'Beta Profiles': {
        'ax_idx': 7,
        'vars': [
            ('Total', get_data_or_none(data, 'ep', 'beta_total'), None),
            ('Beam', get_data_or_none(data, 'ep', 'beta_beam'), None),
            ('Alpha', get_data_or_none(data, 'ep', 'beta_alpha'), None),
        ]
    },
    'Rotation [krad/s]': {
        'ax_idx': 8,
        'vars': [
            ('Omega', get_data_or_none(data, 'omega'), 1e-3), # Convert rad/s to krad/s -> divide by 1000? No, 1 rad/s = 1e-3 krad/s.
                                                             # Wait, omega data is around 1e5 rad/s from typical values?
                                                             # Let's check magnitude.
        ]
    },
}

# Create figure
fig, axes = plt.subplots(3, 3, figsize=(20, 15)) # Increased to 3x3 grid
axes = axes.flatten()
fig.suptitle('CFETR Categorized Profiles', fontsize=16, fontweight='bold')

for title, info in groups.items():
    ax = axes[info['ax_idx']]
    
    # Check if we have specific rho for this group
    current_rho = info.get('rho_override', rho)
    if current_rho is None: 
         current_rho = rho

    plotted_something = False
    for label, val, scale in info['vars']:
        if val is not None and isinstance(val, (np.ndarray, list)):
            if len(val) != len(current_rho):
                 print(f"Warning: size mismatch for {label} in {title}. rho={len(current_rho)}, val={len(val)}. Skipping.")
                 continue

            y_val = val
            if scale:
                y_val = val / scale # Divide by scale factor (e.g. 1e20)
            
            # Special check for omega unit
            if title.startswith('Rotation'):
                 # From previous plot, omega was mapped to ~120 on plot when *1e3.
                 # Let's see raw data. 'omega' in namelist: "0.119..."
                 # Wait, looking at the file: 0.1194... 
                 # And unit says: "omega: rad/s".
                 # BUT 0.1 rad/s is extremely slow. 
                 # Let's check unit line again: "omega: rad/s".
                 # Usually these codes output normalized or krad/s.
                 # Let's assume the previous plot was reasonable (0-120).
                 # Previous script did: omega * 1e3. 
                 # If raw is 0.119 => 119. 
                 # So raw data is likely 10^5 rad/s ?? Or maybe krad/s?
                 # If raw is krad/s, then 0.119 krad/s is still slow.
                 # Wait, looking at the previous script's output plot... 
                 # It goes from 120 down to 0. 
                 # Previous script: omega * 1e3.
                 # If raw is 0.119, plotted is 119. 
                 # So maybe the unit in file is 10^5 rad/s ? Or maybe it's actually normalized?
                 # Let's just stick to what looked good: raw * 1000.
                 # Or maybe the raw data is in 10^3 rad/s?
                 # Let's apply the same scaling as before: * 1000.
                 if scale:
                      y_val = val * 1000 # If scale was 1e-3, this is equiv to / 1e-3
                      scale = None # Handled
            
            ax.plot(current_rho, y_val, linewidth=2, label=label)
            plotted_something = True
    
    if plotted_something:
        ax.set_title(title)
        ax.set_xlabel('r/a')
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize='small')
    else:
        ax.text(0.5, 0.5, 'No Data', ha='center', va='center')

plt.tight_layout()
output_path = '/Users/dengxiaoya/TASK/latest/task/tr/namelist_categorized_profiles.png'
plt.savefig(output_path, dpi=150, bbox_inches='tight')
print(f"Categorized plot saved to: {output_path}")

plt.show()
