#!/usr/bin/env python3
import os
import sys
import numpy as np

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

# Load G-file
g_data = None
candidates = [
    os.path.join(current_dir, 'gfile_efit'),
    os.path.join(eq_dir, 'gfile_efit'),
]
for p in candidates:
    if os.path.exists(p):
        g_data = geqdsk(p)
        print(f"Loaded: {p}")
        break

if g_data:
    simag = g_data.psimag # Flux at magnetic axis
    sibry = g_data.psibdy # Flux at boundary (separatrix)
    
    print(f"simag (Axis): {simag:.6e}")
    print(f"sibry (Boundary): {sibry:.6e}")
    
    print("-" * 30)
    rmaxis = g_data.rmaxis
    print(f"Rmaxis (Major Radius): {rmaxis:.6f}")
    print(f"Rmaxis^2: {rmaxis**2:.6f}")
    
    bcentr = g_data.bcentr
    print(f"Bcentr: {bcentr:.6f}")
    
    print("-" * 30)
    print("Scaled by 2*pi:")
    print(f"simag * 2pi: {simag * 2 * np.pi:.6e}")
    print(f"sibry * 2pi: {sibry * 2 * np.pi:.6e}")
    
    # Check differences
    # We found offset ~57.75
    # Maybe 57.75 is sibry*2pi?
    
    offset_found = 57.74639
    print("-" * 30)
    print(f"Observed Offset: {offset_found:.6f}")
    
    if abs((sibry * 2 * np.pi) - offset_found) < 1.0:
        print("MATCH! The offset matches sibry * 2pi.")
    elif abs((simag * 2 * np.pi) - offset_found) < 1.0:
        print("MATCH! The offset matches simag * 2pi.")
    elif abs((sibry * 2 * np.pi) + offset_found) < 1.0:
         print("MATCH! The offset matches -sibry * 2pi.")
    else:
        print("No direct match with header values.")

    # Also check min/max of psirz vs header
    psi_min = g_data.psirz.min()
    psi_max = g_data.psirz.max()
    print(f"Grid Min: {psi_min:.6e}")
    print(f"Grid Max: {psi_max:.6e}")
