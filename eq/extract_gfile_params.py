#!/usr/bin/env python3
"""
Extract parameters from G-file:
- Minor Radius (a)
- Major Radius (R_geo)
- Triangularity (delta)
- Elongation (kappa)
- Plasma Current (Ip)

Based on gfile reading method in in/compare_psi_3way.py
"""

import sys
import os
import numpy as np

# 1. Setup paths to helper scripts
current_dir = os.getcwd()
# Assuming this script is in TASK/latest/task/eq/ or TASK/latest/task/eq/in/
# We need to find the 'gfile' directory which seems to be at ../../../../gfile from in/compare_psi_3way.py
# If we are in eq/, it might be ../../../gfile.
# Let's try to be robust similar to the reference script.

# If this script is in 'eq', then 'in' is a subdir.
# Reference script was in 'eq/in/'.
# Its path logic: 
# eq_dir = ..
# gfile_dir = ../../../../gfile

# Let's adapt. We will try to locate 'gfile' module.
# We will define a search strategy.

script_dir = os.path.dirname(os.path.abspath(__file__))

# Potential locations for gfile module relative to this script
# If script is in eq/ : ../../../gfile
# If script is in eq/in/: ../../../../gfile
possible_gfile_dirs = [
    os.path.abspath(os.path.join(script_dir, '../../../gfile')),
    os.path.abspath(os.path.join(script_dir, '../../../../gfile')),
    os.path.abspath(os.path.join(script_dir, '../gfile')), # Just in case
    '/Users/dengxiaoya/TASK/gfile' # Absolute fallback
]

gfile_module_path = None
for p in possible_gfile_dirs:
    if os.path.exists(os.path.join(p, 'geqdsk.py')):
        gfile_module_path = p
        break

if gfile_module_path:
    if gfile_module_path not in sys.path:
        sys.path.append(gfile_module_path)
    print(f"Found gfile module at: {gfile_module_path}")
else:
    print("Error: Could not find 'gfile' directory containing geqdsk.py")
    sys.exit(1)

try:
    from geqdsk import geqdsk
except ImportError:
    print("Error: Could not import geqdsk class")
    sys.exit(1)

# 2. Load G-file
# Search locations for the gfile itself
candidates = [
    os.path.join(script_dir, 'in', 'gfile_efit'),
    os.path.join(script_dir, 'gfile_efit'),
    os.path.join(script_dir, 'in', 'gfile_efit_bak'),
    os.path.join(script_dir, 'gfile_efit_bak'),
    # Add absolute paths if known or commonly used
]

gfile_path = None
for p in candidates:
    if os.path.exists(p):
        gfile_path = p
        break

if not gfile_path:
    print("Error: G-file (gfile_efit) not found in standard locations.")
    print(f"Searched: {candidates}")
    sys.exit(1)

print(f"Reading G-file: {gfile_path}")
g = geqdsk(gfile_path)

# 3. Calculate Parameters

# Plasma Boundary (LCFS)
rbbbs = np.array(g.rbbbs)
zbbbs = np.array(g.zbbbs)

# Geometric parameters
R_min = np.min(rbbbs)
R_max = np.max(rbbbs)
Z_min = np.min(zbbbs)
Z_max = np.max(zbbbs)

R_geo = (R_max + R_min) / 2.0
a = (R_max - R_min) / 2.0
aspect_ratio = R_geo / a

# Elongation (Kappa)
# Standard definition: kappa = (Z_max - Z_min) / (2*a)
kappa = (Z_max - Z_min) / (2.0 * a)

# Triangularity (Delta)
# Find R at Z_max (top) and Z_min (bottom)
# Since the points are discrete, we find the point with max/min Z
idx_top = np.argmax(zbbbs)
idx_bot = np.argmin(zbbbs)

R_top = rbbbs[idx_top]
R_bot = rbbbs[idx_bot]

delta_u = (R_geo - R_top) / a
delta_l = (R_geo - R_bot) / a
delta = (delta_u + delta_l) / 2.0


# Current
Ip = g.currentA

# Safety Factor on Axis
q_axis = g.qpsi[0]

# --- Internal Inductance (li) Calculation ---
# li = (2 / (R_geo * (mu0 * Ip)**2)) * integral(Bp^2 dV)
# dV = 2 * pi * R * dR * dZ
# Bp^2 = (grad(psi)/R)^2 = ((dpsi/dR)^2 + (dpsi/dZ)^2) / R^2
# Integraion domain: Inside LCFS (psi between psimag and psibdy)

# 1. Setup Grid
nw = g.nw
nh = g.nh
R_1d = np.linspace(g.rleft, g.rleft + g.rdim, nw)
Z_1d = np.linspace(g.zmid - g.zdim/2, g.zmid + g.zdim/2, nh)
R_2d, Z_2d = np.meshgrid(R_1d, Z_1d, indexing='ij') # Shape (nw, nh)

# 2. Gradient of Psi
# g.psirz is (nw, nh)
dpsi_dR, dpsi_dZ = np.gradient(g.psirz, R_1d, Z_1d) # Note: check gradient axis order with indexing='ij'
# With indexing='ij', axis 0 is R, axis 1 is Z.
# np.gradient(f, x, y) returns df/dx (axis 0), df/dy (axis 1)

# 3. Bp Squared
# Bp^2 = ( (dpsi/dR)^2 + (dpsi/dZ)^2 ) / R^2
Bp2 = (dpsi_dR**2 + dpsi_dZ**2) / (R_2d**2)

# 4. Define Mask (Inside LCFS)
# Check polarity of psi
# Usually psimag is min or max.
psi_mag = g.psimag
psi_bdy = g.psibdy

# Determine if we go up or down
is_inside = np.zeros_like(g.psirz, dtype=bool)

if psi_mag < psi_bdy:
    # Psi increases from axis to edge
    is_inside = (g.psirz >= psi_mag) & (g.psirz <= psi_bdy)
else:
    # Psi decreases from axis to edge
    is_inside = (g.psirz <= psi_mag) & (g.psirz >= psi_bdy)

# 5. Volume Integration
# integral(Bp^2 * 2 * pi * R * dR * dZ)
mu0 = 4 * np.pi * 1e-7
dR = R_1d[1] - R_1d[0]
dZ = Z_1d[1] - Z_1d[0]

# Element volume dV = 2 * pi * R * dR * dZ
dV = 2 * np.pi * R_2d * dR * dZ

# Integrate only inside mask
energy_integral = np.sum(Bp2[is_inside] * dV[is_inside])

# Calculate li
# li(3) definition
li = (2 * energy_integral) / (R_geo * (mu0 * Ip)**2)


# Output results
print("\n" + "="*30)
print(" Extracted G-file Parameters")
print("="*30)
print(f"{'Major Radius (R_geo)':<25}: {R_geo:.6f} m")
print(f"{'Minor Radius (a)':<25}: {a:.6f} m")
print(f"{'Aspect Ratio (R/a)':<25}: {aspect_ratio:.6f}")
print(f"{'Elongation (Kappa)':<25}: {kappa:.6f}")
print(f"{'Triangularity (Delta)':<25}: {delta:.6f}")
print(f"  {'Upper Delta':<23}: {delta_u:.6f}")
print(f"  {'Lower Delta':<23}: {delta_l:.6f}")
print(f"{'Plasma Current (Ip)':<25}: {Ip:.6e} A")
print(f"{'Safety Factor (q0)':<25}: {q_axis:.6f}")
print(f"{'Internal Inductance (li)':<25}: {li:.6f}")
print(f"{'Magnetic Axis (R_mag)':<25}: {g.rmaxis:.6f} m")
print(f"{'Magnetic Axis (Z_mag)':<25}: {g.zmaxis:.6f} m")
print("="*30 + "\n")

