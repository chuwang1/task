import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import sys
import os

from compare_psi_3way import load_eqdata, load_gfile_profiles_csv, compare_profile, resolve_eqdata_prefix_paths, resolve_eq_q_profile

current_dir = os.getcwd()
eq_dir = os.path.dirname(current_dir)

prefix, paths = resolve_eqdata_prefix_paths(eq_dir, current_dir, 'eqdata_modelg5_rebuild_mdleqf9_qmap_samegeom')
eq_q_override = resolve_eq_q_profile(current_dir, eq_dir, None)
eq = load_eqdata(paths, eq_q_override=eq_q_override)

g_profile = load_gfile_profiles_csv('g260206.20000_teq_0114_profiles_from_gfile.csv')

q_m = compare_profile(g_profile["psi_n_profile"], g_profile["q"], eq["psi_n_q"], eq["q"], psi_max=0.98)
p_m = compare_profile(g_profile["psi_n_profile"], g_profile["p"], eq["psi_n_p"], eq["p"])

print("q_m bounds:")
if q_m:
    print(f"ref: min={np.min(q_m['ref'])}, max={np.max(q_m['ref'])}")
    print(f"cmp: min={np.min(q_m['cmp'])}, max={np.max(q_m['cmp'])}")
    # any nan?
    print(f"nan in ref: {np.isnan(q_m['ref']).sum()}, cmp: {np.isnan(q_m['cmp']).sum()}")
print("p_m bounds:")
if p_m:
    print(f"ref: min={np.min(p_m['ref'])}, max={np.max(p_m['ref'])}")
    print(f"cmp: min={np.min(p_m['cmp'])}, max={np.max(p_m['cmp'])}")

# Check what ax.plot would actually return
fig, axes = plt.subplots(1, 2)
axes[0].plot(q_m["psi"], q_m["ref"], "b-", label="ref")
axes[0].plot(q_m["psi"], q_m["cmp"], "r--", label="cmp")
axes[1].plot(p_m["psi"], p_m["ref"], "b-", label="ref")
axes[1].plot(p_m["psi"], p_m["cmp"], "r--", label="cmp")

# Let's see the extent of the plotted lines
print("x_bound q:", axes[0].get_xlim())
print("y_bound q:", axes[0].get_ylim())
print("x_bound p:", axes[1].get_xlim())
print("y_bound p:", axes[1].get_ylim())

