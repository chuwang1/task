$ /Users/dengxiaoya/miniforge3/bin/conda run --no-capture-output -n Fusion python - <<'PY'
import re
import numpy as np
import pandas as pd
from pathlib import Path

gfile = Path('in/g260206.20000_teq_0114')
with gfile.open('r', encoding='utf-8') as f:
    head = re.match(r"^\s*(.*)\s+(\d+)\s+(\d+)\s*$", f.readline())
    nw, nh = int(head.group(2)), int(head.group(3))
    def read5():
        line = re.sub(r"([^Ee])-", r"\1 -", f.readline())
        vals = [float(x) for x in re.split(r"\s+", line.strip()) if x]
        return vals[:5]
    r2 = read5()

rdim, zdim, rcentr, rleft, zmid = r2
R = np.linspace(rleft, rleft + rdim, nw)
Z = np.linspace(zmid - 0.5 * zdim, zmid + 0.5 * zdim, nh)
psi = pd.read_csv('in/gs_from_gfile_psirz_solved.csv').values
j0 = int(np.argmin(np.abs(Z - 0.0)))
out = pd.DataFrame({'R': R, 'psi_at_Z0': psi[:, j0]})
out_path = Path('in/gs_from_gfile_outer1_psi_Z0.csv')
out.to_csv(out_path, index=False)
print(f'saved: {out_path.resolve()}')
print(f'Z_index={j0}, Z_value={Z[j0]:.12g}, N={len(R)}')
PY

saved: /Users/dengxiaoya/TASK/latest/task/eq/in/gs_from_gfile_outer1_psi_Z0.csv
Z_index=128, Z_value=0, N=257
$ /Users/dengxiaoya/miniforge3/bin/conda run --no-capture-output -n Fusion python - <<'PY'
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

gfile = Path('in/g260206.20000_teq_0114')
with gfile.open('r', encoding='utf-8') as f:
    head = re.match(r"^\s*(.*)\s+(\d+)\s+(\d+)\s*$", f.readline())
    nw, nh = int(head.group(2)), int(head.group(3))

    def read5():
        line = re.sub(r"([^Ee])-", r"\1 -", f.readline())
        vals = [float(x) for x in re.split(r"\s+", line.strip()) if x]
        return vals[:5]

    r2 = read5()
    read5(); read5(); read5()
    body = f.read().replace('\n', ' ')

body = re.sub(r"(\d)-", r"\1 -", body)
nums = [float(x) for x in re.split(r"\s+", body) if x]

rdim, zdim, rcentr, rleft, zmid = r2
R = np.linspace(rleft, rleft + rdim, nw)
Z = np.linspace(zmid - 0.5 * zdim, zmid + 0.5 * zdim, nh)

ia = 0
ia += nw  # fpol
ia += nw  # pres
ia += nw  # ffprime
ia += nw  # pprime
psirz_g = np.asarray(nums[ia: ia + nw * nh], dtype=float).reshape(nw, nh)  # [R,Z]

psi_sol = pd.read_csv('in/gs_from_gfile_psirz_solved.csv').values  # [R,Z]
j0 = int(np.argmin(np.abs(Z - 0.0)))

psi_g_z0 = psirz_g[:, j0]
psi_s_z0 = psi_sol[:, j0]
delta = psi_s_z0 - psi_g_z0

fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), constrained_layout=True)
axes[0].plot(R, psi_g_z0, lw=2.0, label='gfile psi(Z=0)')
axes[0].plot(R, psi_s_z0, lw=1.8, ls='--', label='solved psi(Z=0), outer=1')
axes[0].set_xlabel('R (m)')
axes[0].set_ylabel('psi')
axes[0].set_title('Z=0 Psi Profile Comparison')
axes[0].grid(alpha=0.3)
axes[0].legend()

axes[1].plot(R, delta, color='tab:red', lw=1.8)
axes[1].axhline(0.0, color='k', lw=1.0, alpha=0.6)
axes[1].set_xlabel('R (m)')
axes[1].set_ylabel('psi_solved - psi_gfile')
axes[1].set_title('Difference at Z=0')
axes[1].grid(alpha=0.3)

out_png = Path('in/gs_from_gfile_outer1_Z0_compare.png')
out_csv = Path('in/gs_from_gfile_outer1_Z0_compare.csv')
fig.savefig(out_png, dpi=160)
plt.close(fig)

pd.DataFrame({
    'R': R,
    'psi_gfile_Z0': psi_g_z0,
    'psi_solved_outer1_Z0': psi_s_z0,
    'delta': delta,
}).to_csv(out_csv, index=False)

print(f'saved: {out_png.resolve()}')
print(f'saved: {out_csv.resolve()}')
print(f'Z_index={j0}, Z_value={Z[j0]:.12g}')
print(f'RMSE_Z0={np.sqrt(np.mean(delta**2)):.6e}, MAXABS_Z0={np.max(np.abs(delta)):.6e}')
PY

saved: /Users/dengxiaoya/TASK/latest/task/eq/in/gs_from_gfile_outer1_Z0_compare.png
saved: /Users/dengxiaoya/TASK/latest/task/eq/in/gs_from_gfile_outer1_Z0_compare.csv
Z_index=128, Z_value=0
RMSE_Z0=9.116668e-01, MAXABS_Z0=4.588506e+00
