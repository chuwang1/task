#!/usr/bin/env python3
import os
import sys
import numpy as np
import pandas as pd

# Setup paths
current_dir = os.getcwd()
eq_dir = os.path.abspath(os.path.join(current_dir, '..'))
if eq_dir not in sys.path: sys.path.append(eq_dir)
gfile_dir = os.path.abspath(os.path.join(current_dir, '../../../../gfile'))
if gfile_dir not in sys.path: sys.path.append(gfile_dir)

try:
    from geqdsk import geqdsk
except ImportError:
    print("Could not import geqdsk")
    sys.exit(1)

# 1. Load G-file (Raw)
g_data = None
candidates = [
    os.path.join(current_dir, 'gfile_efit'),
    os.path.join(eq_dir, 'gfile_efit'),
]
for p in candidates:
    if os.path.exists(p):
        g_data = geqdsk(p)
        break

if g_data is None:
    print("G-file not found")
    sys.exit(1)

# PSI_g raw (no transpose yet, let's verify index match)
# geqdsk.psirz is (nw, nh). 
# If R_g is linspace(rleft, ...), R index is 1st dimension.
PSI_g = g_data.psirz 

# 2. Load Eqdata (Raw)
eqdata_psi_path = os.path.join(eq_dir, 'eqdata_PSIRZ.csv')
if not os.path.exists(eqdata_psi_path):
    eqdata_psi_path = os.path.join(current_dir, 'eqdata_PSIRZ.csv')

if not os.path.exists(eqdata_psi_path):
    print("Eqdata not found")
    sys.exit(1)

try:
    PSI_vals = pd.read_csv(eqdata_psi_path, comment='#').values
except:
    PSI_vals = pd.read_csv(eqdata_psi_path, header=None, comment='#').values

# eqdata CSV is usually (NR, NZ) based on my previous discovery (read_eqdata_to_csv saves flattened F-order reshape).
# So indices should match G-file directly if grids are same.
PSI_eqd = PSI_vals

# 3. Compare at specific point
# Center index (assuming 129x129)
idx_r = g_data.nw // 2
idx_z = g_data.nh // 2

val_g = PSI_g[idx_r, idx_z]
val_e = PSI_eqd[idx_r, idx_z]

print(f"Grid Size: {g_data.nw}x{g_data.nh}")
print(f"Comparison at Index ({idx_r}, {idx_z}):")
print(f"  G-file Raw: {val_g:.6e}")
print(f"  Eqdata Raw: {val_e:.6e}")

# Check scaling
if abs(val_g) > 1e-10:
    ratio = val_e / val_g
    print(f"  Ratio (Eqdata/Gfile): {ratio:.6f}")
    print(f"  2*pi: {2*np.pi:.6f}")

# Check offset with 2pi scaling
val_g_scaled = val_g * 2 * np.pi
diff = val_e - val_g_scaled
print(f"  G-file * 2pi: {val_g_scaled:.6e}")
print(f"  Difference (Eqdata - Gfile*2pi): {diff:.6e}")

# Check magnetic axis (peak)
# G-file stores axis index sometimes, or we find max/min
# Assuming plasma, usually min of flux (or max depending on sign)
min_g = PSI_g.min()
min_e = PSI_eqd.min()
print(f"\nMinimum Values (Global):")
print(f"  G-file min: {min_g:.6e}")
print(f"  Eqdata min: {min_e:.6e}")
print(f"  Ratio: {min_e/min_g:.6f}")
print(f"  Diff (min_e - min_g*2pi): {min_e - min_g*2*np.pi:.6e}")

max_g = PSI_g.max()
max_e = PSI_eqd.max()
print(f"\nMaximum Values (Global):")
print(f"  G-file max: {max_g:.6e}")
print(f"  Eqdata max: {max_e:.6e}")

range_g = max_g - min_g
range_e = max_e - min_e
print(f"\nRange:")
print(f"  G-file range: {range_g:.6e}")
print(f"  Eqdata range: {range_e:.6e}")
print(f"  Range Ratio: {range_e/range_g:.6f}")
