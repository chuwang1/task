#!/usr/bin/env python3
"""
Plot PNFCL (Fusion Collisional Heating Power) profile from TR output CSV
"""
import numpy as np
import matplotlib.pyplot as plt

def read_tr_csv(filename):
    """Read TR output CSV file"""
    data = []
    with open(filename, 'r') as f:
        lines = f.readlines()
        # Skip header lines
        for line in lines[2:]:
            if line.strip():
                values = [float(x) for x in line.split(',')]
                data.append(values)
    return np.array(data)

# Read data
data = read_tr_csv('tr_data_028.csv')

r_a = data[:, 0]           # r/a
PNFIN = data[:, 1]         # Fusion input power [MW/m³]
PNFCL_e = data[:, 2]       # PNFCL to electrons [MW/m³]
PNFCL_D = data[:, 3]       # PNFCL to D ions [MW/m³]
PNFCL_T = data[:, 4]       # PNFCL to T ions [MW/m³]
PNFCL_He4 = data[:, 5]     # PNFCL to He4 ions [MW/m³]

# Total PNFCL to ions
PNFCL_i = PNFCL_D + PNFCL_T + PNFCL_He4

# Create figure
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Plot 1: All components
ax1 = axes[0]
ax1.plot(r_a, PNFIN, 'k-', linewidth=2, label='PNFIN (Total input)')
ax1.plot(r_a, PNFCL_e, 'r-', linewidth=2, label='PNFCL_e (to electrons)')
ax1.plot(r_a, PNFCL_D, 'b--', linewidth=1.5, label='PNFCL_D (to D)')
ax1.plot(r_a, PNFCL_T, 'g--', linewidth=1.5, label='PNFCL_T (to T)')
ax1.plot(r_a, PNFCL_He4, 'm--', linewidth=1.5, label='PNFCL_He4 (to He4)')
ax1.set_xlabel('r/a', fontsize=12)
ax1.set_ylabel('Power density [MW/m³]', fontsize=12)
ax1.set_title('PNFCL Components vs r/a', fontsize=14)
ax1.legend(loc='upper right')
ax1.grid(True, alpha=0.3)
ax1.set_xlim([0, 1])

# Plot 2: Electron vs Ion heating
ax2 = axes[1]
ax2.plot(r_a, PNFIN, 'k-', linewidth=2, label='PNFIN (Total)')
ax2.plot(r_a, PNFCL_e, 'r-', linewidth=2, label='PNFCL_e (to electrons)')
ax2.plot(r_a, PNFCL_i, 'b-', linewidth=2, label='PNFCL_i (to ions)')
ax2.fill_between(r_a, 0, PNFCL_e, alpha=0.3, color='red', label='_nolegend_')
ax2.fill_between(r_a, PNFCL_e, PNFCL_e+PNFCL_i, alpha=0.3, color='blue', label='_nolegend_')
ax2.set_xlabel('r/a', fontsize=12)
ax2.set_ylabel('Power density [MW/m³]', fontsize=12)
ax2.set_title('Fusion Heating: Electron vs Ion', fontsize=14)
ax2.legend(loc='upper right')
ax2.grid(True, alpha=0.3)
ax2.set_xlim([0, 1])

# Calculate and print fractions
total_e = np.trapz(PNFCL_e * r_a, r_a)
total_i = np.trapz(PNFCL_i * r_a, r_a)
total = total_e + total_i

print("=" * 50)
print("PNFCL Energy Distribution")
print("=" * 50)
print(f"Fraction to electrons: {total_e/total*100:.1f}%")
print(f"Fraction to ions:      {total_i/total*100:.1f}%")
print(f"  - to D:   {np.trapz(PNFCL_D * r_a, r_a)/total*100:.1f}%")
print(f"  - to T:   {np.trapz(PNFCL_T * r_a, r_a)/total*100:.1f}%")
print(f"  - to He4: {np.trapz(PNFCL_He4 * r_a, r_a)/total*100:.1f}%")
print("=" * 50)
print(f"Peak PNFIN at r/a = {r_a[np.argmax(PNFIN)]:.2f}: {np.max(PNFIN):.3f} MW/m³")
print(f"Peak PNFCL_e at r/a = {r_a[np.argmax(PNFCL_e)]:.2f}: {np.max(PNFCL_e):.3f} MW/m³")

plt.tight_layout()
plt.savefig('pnfcl_profile.png', dpi=150)
plt.show()
print("\nSaved: pnfcl_profile.png")
