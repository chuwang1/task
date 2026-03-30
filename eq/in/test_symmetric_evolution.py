import sys
sys.path.append('in')
from solve_gs_from_gfile import read_gfile, solve_gs_fixed_boundary, resolve_profile_csv, load_source_profiles
import numpy as np
import pandas as pd
from scipy.interpolate import RectBivariateSpline
import matplotlib.pyplot as plt

def main():
    g = read_gfile('in/g260206.20000_teq_0114')
    print("Gfile shape:", g.psirz.shape)
    
    psirz_raw = pd.read_csv('eqdata0114Test_PSIRZ.csv', comment="#").to_numpy(dtype=float).T
    rg = pd.read_csv('eqdata0114Test_RG_grid.csv')["RG"].to_numpy(dtype=float)
    zg = pd.read_csv('eqdata0114Test_ZG_grid.csv')["ZG"].to_numpy(dtype=float)
    
    print(f"rg: {len(rg)}, zg: {len(zg)}, psirz: {psirz_raw.shape}")
    
    # Check orientation and interpolate to (g.R, g.Z)
    if psirz_raw.shape == (len(rg), len(zg)):
        spline = RectBivariateSpline(rg, zg, psirz_raw)
        mapped = spline(g.R, g.Z)
    else:
        spline = RectBivariateSpline(zg, rg, psirz_raw.T)
        mapped = spline(g.Z, g.R).T
        
    print(f"Interpolated map shape: {mapped.shape}")
    
    # Overwrite the initial guess: eqdata -> gfile convention
    # eqdata: PSIRZ = 2*pi*(psi_gfile - psibdy)  =>  psi_gfile = PSIRZ / 2pi + psibdy
    g.psirz = mapped / (2.0 * np.pi) + g.psibdy
    
    # Also overwrite initial boundary conditions
    import copy
    orig_psirz = copy.deepcopy(read_gfile('in/g260206.20000_teq_0114').psirz)
    
    # We will keep the boundary as the symmetric one, or replace it?
    # Usually in fixed boundary, the boundary comes from g.psirz directly in the solve function.
    print("Replaced initial guess with symmetric PSIRZ.")

    # Load source profiles (from the target gfile, to drive the evolution)
    prof_csv = resolve_profile_csv('in/g260206.20000_teq_0114', None)
    src_x, pprime, ffprime, _, src_kind = load_source_profiles(g, prof_csv, "rho_tor")
    print(f"Loaded target profiles: {len(src_x)} points")
    
    # Run the solver
    print("Starting solver...")
    
    psi_sol, res_hist, err_hist = solve_gs_fixed_boundary(
        g,
        src_x,
        pprime,
        ffprime,
        source_x_kind=src_kind,
        max_outer=60,
        max_inner=2000,
        omega=0.2,            # Heavy under-relaxation (omega)
        boundary_mode="rectangle",
        mask_mode="both"
    )
    
    print(f"Finished in {len(res_hist)} outer iterations.")
    print("Residual history:", res_hist)
    
    # Extract the original LCFS to compare
    orig_g = read_gfile('in/g260206.20000_teq_0114')
    
    fig, ax = plt.subplots(figsize=(6,8))
    # Original target (red dashed)
    ax.contour(orig_g.R, orig_g.Z, orig_g.psirz.T, levels=40, colors="red", linestyles="dashed", linewidths=1.0)
    # Evolved steady state (blue solid)
    ax.contour(g.R, g.Z, psi_sol.T, levels=40, colors="blue", linestyles="solid", linewidths=1.2)
    # Initial symmetric guess (green dotted)
    ax.contour(g.R, g.Z, g.psirz.T, levels=40, colors="green", linestyles="dotted", linewidths=1.0)
    
    ax.set_aspect('equal')
    ax.set_title("Evolution from Symmetric Guess to Target State")
    from matplotlib.lines import Line2D
    ax.legend(handles=[
        Line2D([0],[0], color="green", ls=":", label="Initial Symmetric Guess"),
        Line2D([0],[0], color="blue", ls="-", label="Evolved Steady State"),
        Line2D([0],[0], color="red", ls="--", label="Target gfile"),
    ])
    fig.savefig('in/symmetric_evolution_test.png', dpi=150)
    print("Saved overlay to in/symmetric_evolution_test.png")
    
if __name__ == "__main__":
    main()
