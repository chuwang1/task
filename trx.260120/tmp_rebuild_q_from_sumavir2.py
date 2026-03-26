import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path

from plot_sumavir2_from_gfile import read_gfile_2d, read_input_profiles_mapping, bilinear_interp

base = '/Users/dengxiaoya/TASK/CFEDR/git/task/trx.260120'
OMFIT_BASE = '/Users/dengxiaoya/CFEDRSW/OMFIT_out'
GFILE = os.path.join(OMFIT_BASE, 'g260206.20000_teq_0114')
INPUT_PROFILES = os.path.join(OMFIT_BASE, 'input.profiles')

lines = open(os.path.join(base, 'eqcalq_dvdpsit_raw.csv')).read().splitlines()
rho = []
qps = []
tts = []
for line in lines[1:]:
    parts = line.split(',')
    vals = [float(x) for x in parts[1].split()]
    rho.append(vals[0])
    qps.append(vals[5])
    tts.append(vals[6])
rho = np.array(rho)
qps = np.array(qps)
tts = np.array(tts)

R, Z, psi_n_2d, psi_p_2d, rmaxis, zmaxis = read_gfile_2d(GFILE)
psi_n_map, rho_map = read_input_profiles_mapping(INPUT_PROFILES)
dpsi_dZ, dpsi_dR = np.gradient(psi_p_2d, Z, R)

fig_dummy, ax = plt.subplots()
sumavir2_g = np.zeros_like(rho_map)
for i, p_target in enumerate(psi_n_map):
    if p_target <= 1e-6:
        sumavir2_g[i] = 0.0
        continue
    cs = ax.contour(R, Z, psi_n_2d, levels=[p_target])
    paths = []
    if len(cs.allsegs) > 0 and len(cs.allsegs[0]) > 0:
        for seg in cs.allsegs[0]:
            paths.append(Path(seg))
    if not paths:
        sumavir2_g[i] = np.nan
        continue
    valid = []
    for p in paths:
        v = p.vertices
        if np.min(v[:, 0]) < rmaxis < np.max(v[:, 0]) and np.min(v[:, 1]) < zmaxis < np.max(v[:, 1]):
            valid.append(p)
    if not valid:
        valid = paths
    path = max(valid, key=lambda x: len(x.vertices))
    v = path.vertices
    if not np.allclose(v[0], v[-1]):
        v = np.vstack([v, v[0]])
    s = 0.0
    for k in range(1, len(v)):
        x1, y1 = v[k - 1]
        x2, y2 = v[k]
        dl = np.hypot(x2 - x1, y2 - y1)
        rm = 0.5 * (x1 + x2)
        zm = 0.5 * (y1 + y2)
        dpr = bilinear_interp(R, Z, dpsi_dR, rm, zm)
        dpz = bilinear_interp(R, Z, dpsi_dZ, rm, zm)
        bp = np.sqrt(dpr * dpr + dpz * dpz) / (2.0 * np.pi * rm)
        if bp <= 1e-30 or rm <= 1e-30:
            continue
        s += dl / (bp * rm * rm)
    sumavir2_g[i] = s
plt.close(fig_dummy)

suma_g_on_tr = np.interp(rho, rho_map, sumavir2_g)
q_rebuilt = suma_g_on_tr * tts / (4.0 * np.pi ** 2)
allp = pd.read_csv('/Users/dengxiaoya/CFEDRSW/OMFIT_out/outputs_new/all_profiles.csv')
q_g = np.interp(rho, allp['rho'].values, np.abs(allp['q'].values))

print('mean rel err rebuilt q vs gfile q =', float(np.mean(np.abs(q_rebuilt - q_g) / np.maximum(np.abs(q_g), 1e-30))))
print('mean rel err TR q vs gfile q      =', float(np.mean(np.abs(qps - q_g) / np.maximum(np.abs(q_g), 1e-30))))
print('mean rel err rebuilt q vs TR q    =', float(np.mean(np.abs(q_rebuilt - qps) / np.maximum(np.abs(qps), 1e-30))))
for idx in [150, 170, 180, 190, 195, 196, 197, 198, 199, 200]:
    print(f'rho={rho[idx]:.4f} q_g={q_g[idx]:.4f} q_rebuilt={q_rebuilt[idx]:.4f} q_tr={qps[idx]:.4f}')

fig, axes = plt.subplots(2, 1, figsize=(8, 8), sharex=True)
axes[0].plot(rho, q_g, '--', lw=2.0, color='tab:red', label='gfile q')
axes[0].plot(rho, q_rebuilt, '-.', lw=2.0, color='tab:green', label='rebuild q = SUMAVIR2_g*TTS/4pi^2')
axes[0].plot(rho, qps, 'o-', ms=3, lw=1.6, color='tab:blue', label='TR QPS')
axes[0].set_ylabel('q')
axes[0].set_title('q rebuilt from gfile SUMAVIR2 and TR TTS')
axes[0].grid(True, alpha=0.3)
axes[0].legend(fontsize=8)

axes[1].plot(rho, q_rebuilt / np.maximum(q_g, 1e-30), '-.', lw=1.8, color='tab:green', label='q_rebuilt/q_gfile')
axes[1].plot(rho, qps / np.maximum(q_g, 1e-30), '--', lw=1.8, color='tab:blue', label='QPS/q_gfile')
axes[1].axhline(1.0, color='k', lw=1, alpha=0.5)
axes[1].set_xlabel('rho')
axes[1].set_ylabel('ratio')
axes[1].grid(True, alpha=0.3)
axes[1].legend(fontsize=8)

out = os.path.join(base, 'q_rebuilt_from_sumavir2_gfile.png')
plt.tight_layout()
plt.savefig(out, dpi=160)
plt.close(fig)
print('Wrote', out)
