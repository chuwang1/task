#!/usr/bin/env python3
"""
Compare PSI(R, Z) from 3 sources:
1. G-file (using geqdsk.py)
2. eqgs2d CSV (using read_eqgs2d_csv.py)
3. eqdata CSV (eqdata_PSIRZ.csv)
"""

import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# 1. Setup paths to helper scripts
current_dir = os.getcwd()
eq_dir = os.path.abspath(os.path.join(current_dir, '..'))
if eq_dir not in sys.path:
    sys.path.append(eq_dir)

gfile_dir = os.path.abspath(os.path.join(current_dir, '../../../../gfile'))
if gfile_dir not in sys.path:
    sys.path.append(gfile_dir)

print(f"Current Dir: {current_dir}")
print(f"EQ Dir: {eq_dir}")
print(f"GFile Dir: {gfile_dir}")

# 2. Load Data

# --- Source A: G-file ---
g_data = None
try:
    from geqdsk import geqdsk
    # Search locations
    candidates = [
        os.path.join(current_dir, 'gfile_efit'),
        os.path.join(eq_dir, 'gfile_efit'),
        os.path.join(current_dir, 'gfile_efit_bak'),
        os.path.join(eq_dir, 'gfile_efit_bak')
    ]
    
    gfile_path = None
    for p in candidates:
        if os.path.exists(p):
            gfile_path = p
            break
    
    if gfile_path:
        g_data = geqdsk(gfile_path)
        # Construct Grid
        R_g = np.linspace(g_data.rleft, g_data.rleft + g_data.rdim, g_data.nw)
        Z_g = np.linspace(g_data.zmid - 0.5 * g_data.zdim, g_data.zmid + 0.5 * g_data.zdim, g_data.nh)
        PSI_g = g_data.psirz*2*3.14+57.7
        print(f"Loaded G-file: {gfile_path}")
    else:
        print("G-file not found")
except ImportError:
    print("Could not import geqdsk")
except Exception as e:
    print(f"Error loading G-file: {e}")


# --- Source B: eqgs2d CSV ---
PSI_gs2d = None
R_gs2d = None
Z_gs2d = None
try:
    from read_eqgs2d_csv import load_psirz_grid
    # Check parent dir first
    csv_path = os.path.join(eq_dir, 'eqgs2d_01_PSIRZ_grid.csv')
    if not os.path.exists(csv_path):
        csv_path = os.path.join(current_dir, 'eqgs2d_01_PSIRZ_grid.csv')

    if os.path.exists(csv_path):
        R_gs2d, Z_gs2d, PSI_gs2d = load_psirz_grid(csv_path)
        print(f"Loaded eqgs2d CSV: {csv_path}")
    else:
        print("eqgs2d CSV not found")
except ImportError:
    print("Could not import read_eqgs2d_csv")
except Exception as e:
    print(f"Error loading eqgs2d CSV: {e}")


# --- Source C: eqdata CSV ---
PSI_eqd = None
R_eqd = None
Z_eqd = None
try:
    # Check parent then local
    eqdata_psi_path = os.path.join(eq_dir, 'eqdata_PSIRZ.csv')
    if not os.path.exists(eqdata_psi_path):
        eqdata_psi_path = os.path.join(current_dir, 'eqdata_PSIRZ.csv')
        
    eqdata_rg_path = os.path.join(eq_dir, 'eqdata_RG_grid.csv')
    if not os.path.exists(eqdata_rg_path):
        eqdata_rg_path = os.path.join(current_dir, 'eqdata_RG_grid.csv')

    eqdata_zg_path = os.path.join(eq_dir, 'eqdata_ZG_grid.csv')
    if not os.path.exists(eqdata_zg_path):
       eqdata_zg_path = os.path.join(current_dir, 'eqdata_ZG_grid.csv')
    
    if os.path.exists(eqdata_psi_path) and os.path.exists(eqdata_rg_path) and os.path.exists(eqdata_zg_path):
        # Read Grids
        df_rg = pd.read_csv(eqdata_rg_path)
        R_eqd = df_rg['RG'].values
        
        df_zg = pd.read_csv(eqdata_zg_path)
        Z_eqd = df_zg['ZG'].values
        
        # Read Matrix
        try:
            PSI_vals = pd.read_csv(eqdata_psi_path, comment='#').values
        except:
             PSI_vals = pd.read_csv(eqdata_psi_path, header=None, comment='#').values

        PSI_eqd = PSI_vals.T
        
        print(f"Loaded eqdata CSV: {eqdata_psi_path}")
    else:
        print("eqdata CSV files not found")
except Exception as e:
    print(f"Error loading eqdata CSV: {e}")


# 3. Plot Comparison

fig, axes = plt.subplots(1, 3, figsize=(20, 6), constrained_layout=True)

# Common levels
datasets = []
if g_data is not None: datasets.append(PSI_g)
if PSI_gs2d is not None: datasets.append(PSI_gs2d)
if PSI_eqd is not None: datasets.append(PSI_eqd)

if datasets:
    # Print statistics
    if g_data is not None:
        print(f"G-file PSI: min={PSI_g.min():.4e}, max={PSI_g.max():.4e}, range={PSI_g.max()-PSI_g.min():.4e}")
    if PSI_gs2d is not None:
        print(f"Eqgs2d PSI: min={PSI_gs2d.min():.4e}, max={PSI_gs2d.max():.4e}, range={PSI_gs2d.max()-PSI_gs2d.min():.4e}")
    if PSI_eqd is not None:
        print(f"Eqdata PSI: min={PSI_eqd.min():.4e}, max={PSI_eqd.max():.4e}, range={PSI_eqd.max()-PSI_eqd.min():.4e}")

    # Use percentiles to avoid extreme outliers affecting colorbar too much
    vmin = min(np.percentile(d, 1) for d in datasets)
    vmax = max(np.percentile(d, 99) for d in datasets)
    
    # Ensure range isn't zero
    if vmax == vmin:
        vmax += 1.0
        
    levels = np.linspace(vmin, vmax, 40)
else:
    levels = 40

# --- Plot 1: G-file ---
ax = axes[0]
if g_data is not None:
    cf = ax.contourf(R_g, Z_g, PSI_g, levels=levels, cmap='jet')
    ax.set_title('G-file PSI(R,Z)')
    ax.set_xlabel('R (m)')
    ax.set_ylabel('Z (m)')
    fig.colorbar(cf, ax=ax)
    ax.set_aspect('equal')
    # Plot LCFS
    ax.plot(g_data.rbbbs, g_data.zbbbs, 'w--', linewidth=1)
else:
    ax.text(0.5, 0.5, 'G-file data missing', ha='center')

# --- Plot 2: EQGS2D CSV ---
ax = axes[1]
if PSI_gs2d is not None:
    cf = ax.contourf(R_gs2d, Z_gs2d, PSI_gs2d, levels=levels, cmap='jet')
    ax.set_title('EQGS2D CSV PSI(R,Z)')
    ax.set_xlabel('R (m)')
    fig.colorbar(cf, ax=ax)
    ax.set_aspect('equal')
else:
    ax.text(0.5, 0.5, 'EQGS2D data missing', ha='center')

# --- Plot 3: EQDATA CSV ---
ax = axes[2]
if PSI_eqd is not None:
    cf = ax.contourf(R_eqd, Z_eqd, PSI_eqd, levels=levels, cmap='jet')
    ax.set_title('EQDATA CSV PSI(R,Z)')
    ax.set_xlabel('R (m)')
    fig.colorbar(cf, ax=ax)
    ax.set_aspect('equal')
else:
    ax.text(0.5, 0.5, 'EQDATA data missing', ha='center')

plt.suptitle(f'Comparison of PSI(R, Z) from 3 Sources', fontsize=16)
output_img = 'psi_comparison_3way.png'
plt.savefig(output_img, dpi=150)
print(f"Comparison plot saved to: {output_img}")

# --- Additional Plot: Overlay G-file vs Eqdata ---
if g_data is not None and PSI_eqd is not None:
    fig2, ax2 = plt.subplots(figsize=(8, 10))
    
    # Use same levels for direct comparison
    # Plot G-file with solid lines
    cs1 = ax2.contour(R_g, Z_g, PSI_g, levels=levels, colors='blue', linestyles='solid', linewidths=1.5)
    # Add labels to some contours
    # ax2.clabel(cs1, inline=True, fontsize=8, fmt='%.2e')
    
    # Plot Eqdata with dashed lines
    cs2 = ax2.contour(R_eqd, Z_eqd, PSI_eqd, levels=levels, colors='red', linestyles='dashed', linewidths=1.5)
    # ax2.clabel(cs2, inline=True, fontsize=8, fmt='%.2e')
    
    ax2.set_title('Overlay: G-file (Blue/Solid) vs Eqdata (Red/Dashed)')
    ax2.set_xlabel('R (m)')
    ax2.set_ylabel('Z (m)')
    ax2.set_aspect('equal')
    
    # Create custom legend
    from matplotlib.lines import Line2D
    legend_elements = [Line2D([0], [0], color='blue', lw=1.5, label='G-file'),
                       Line2D([0], [0], color='red', lw=1.5, linestyle='--', label='Eqdata')]
    ax2.legend(handles=legend_elements, loc='upper right')
    
    output_overlay = 'gfile_eqdata_overlay.png'
    plt.savefig(output_overlay, dpi=150)
    print(f"Overlay plot saved to: {output_overlay}")

