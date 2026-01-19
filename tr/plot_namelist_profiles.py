#!/usr/bin/env python3
"""
Plot profiles from CFETR namelist file (custom parser for non-standard format)
"""

import numpy as np
import matplotlib.pyplot as plt
import re

def parse_custom_namelist(filepath):
    """Parse the custom namelist format"""
    data = {}
    current_section = None
    current_key = None
    current_values = []
    
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Split into lines and process
    lines = content.split('\n')
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        
        # Skip empty lines
        if not line:
            i += 1
            continue
        
        # Check for section header (&SECTION_NAME)
        if line.startswith('&'):
            section_name = line[1:].lower()
            current_section = section_name
            data[current_section] = {}
            i += 1
            continue
        
        # Check for section end (/)
        if line == '/':
            current_section = None
            i += 1
            continue
        
        # Check for key = value assignment
        if '=' in line:
            # Save previous key-values if exists
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
            # Continuation of values
            current_values.append(line)
        
        i += 1
    
    # Save last key-values
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
    
    # Handle string values (quoted)
    if values_str.startswith("'") or values_str.startswith('"'):
        return values_str.strip("'\"")
    
    # Handle repeated values like 201*0.0
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
    
    # Try to parse as numeric array
    try:
        values = []
        for v in values_str.split():
            values.append(float(v))
        if len(values) == 1:
            return values[0]
        return np.array(values)
    except:
        return values_str

# Parse the namelist file
nml_path = '/Users/dengxiaoya/COREDIV/CFEDR混合运行模式250925版本/CFEDR混合运行模式250925数据/case250925（即case6.2）/profiles_CFETR.namelist'
print("Parsing namelist file...")
data = parse_custom_namelist(nml_path)

# Extract data
rho = data['rho']
Te = data['electron']['temperature']
Ti = data['ions_1']['temperature']
ne = data['electron']['density'] / 1e20
ni = data['ions_1']['density'] / 1e20
j_total = data['current_density'] / 1e6
j_bootstrap = data['bootstrap_current_density_onetwo'] / 1e6
j_ohmic = data['ohmic_current_density_onetwo'] / 1e6
j_rfcd = data['rfcd_current_density_onetwo'] / 1e6
p_total = data['pressure_total']
p_e = data['electron']['pressure']
omega = data['omega'] * 1e3

# Create figure
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
fig.suptitle('CFETR Profiles (case250925)', fontsize=14, fontweight='bold')

# 1. Temperature profiles
ax1 = axes[0, 0]
ax1.plot(rho, Te, 'r-', label='Te', linewidth=2)
ax1.plot(rho, Ti, 'b-', label='Ti (D)', linewidth=2)
ax1.set_xlabel('r/a')
ax1.set_ylabel('Temperature [keV]')
ax1.set_title('Temperature Profiles')
ax1.legend()
ax1.grid(True, alpha=0.3)

# 2. Density profiles
ax2 = axes[0, 1]
ax2.plot(rho, ne, 'r-', label='ne', linewidth=2)
ax2.plot(rho, ni, 'b-', label='ni (D)', linewidth=2)
ax2.set_xlabel('r/a')
ax2.set_ylabel('Density [10²⁰ m⁻³]')
ax2.set_title('Density Profiles')
ax2.legend()
ax2.grid(True, alpha=0.3)

# 3. Current density profiles
ax3 = axes[0, 2]
ax3.plot(rho, j_total, 'k-', label='Total', linewidth=2)
ax3.plot(rho, j_bootstrap, 'g-', label='Bootstrap', linewidth=2)
ax3.plot(rho, j_ohmic, 'b-', label='Ohmic', linewidth=2)
ax3.set_xlabel('r/a')
ax3.set_ylabel('Current Density [MA/m²]')
ax3.set_title('Current Density Profiles')
ax3.legend()
ax3.grid(True, alpha=0.3)

# 4. RF current drive
ax4 = axes[1, 0]
ax4.plot(rho, j_rfcd, 'm-', label='RFCD (EC/LH)', linewidth=2)
ax4.axvline(x=0.4, color='gray', linestyle='--', alpha=0.5, label='r/a=0.4')
ax4.set_xlabel('r/a')
ax4.set_ylabel('Current Density [MA/m²]')
ax4.set_title('RF Current Drive Profile')
ax4.legend()
ax4.grid(True, alpha=0.3)

# 5. Pressure profiles
ax5 = axes[1, 1]
ax5.plot(rho, p_total, 'k-', label='Total', linewidth=2)
ax5.plot(rho, p_e, 'r-', label='Electron', linewidth=2)
ax5.set_xlabel('r/a')
ax5.set_ylabel('Pressure [MPa]')
ax5.set_title('Pressure Profiles')
ax5.legend()
ax5.grid(True, alpha=0.3)

# 6. Rotation profile
ax6 = axes[1, 2]
ax6.plot(rho, omega, 'c-', linewidth=2)
ax6.set_xlabel('r/a')
ax6.set_ylabel('ω [krad/s]')
ax6.set_title('Toroidal Rotation Profile')
ax6.grid(True, alpha=0.3)

plt.tight_layout()
output_path = '/Users/dengxiaoya/TASK/latest/task/tr/namelist_profiles.png'
plt.savefig(output_path, dpi=150, bbox_inches='tight')
print(f"Plot saved to: {output_path}")

# Print key parameters
print(f"\nKey parameters:")
print(f"  R0 = {data['r0']:.2f} m")
print(f"  B0 = {data['b0']:.2f} T")
print(f"  Te(0) = {Te[0]:.2f} keV")
print(f"  Ti(0) = {Ti[0]:.2f} keV")
print(f"  ne(0) = {ne[0]:.2f} × 10²⁰ m⁻³")
print(f"  RFCD peak location: r/a ≈ {rho[np.argmax(j_rfcd)]:.2f}")

plt.show()
