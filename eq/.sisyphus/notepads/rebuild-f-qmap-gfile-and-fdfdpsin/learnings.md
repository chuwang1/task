
## Task: Update rebuild_f.sh with gfile-qp and rebuild-csv flags

### Changes Made
1. Added `--use-gfile-qp` flag to `rebuild_f_from_q_eqdata.py` invocation
   - Enables gfile profile data usage in the rebuild process
   - Placed before `--out-prefix` argument

2. Added `--rebuild-csv` argument to `plot_eqipqp_all_variables_big.py` invocation
   - Points to: `in/eqdata_modelg5_rebuild_mdleqf9_qmap_samegeom_f_from_q_eqipqp_stable.csv`
   - This is the CSV output from the rebuild step
   - Prevents plot script from using stale default CSV path

### Key Insight
The plot script's default `--rebuild-csv` would resolve to `in/eqdata_modelg5_rebuild_mdleqf9_qmap_f_from_q_eqipqp.csv` (based on prefix), which differs from the rebuild output. Explicit wiring ensures data consistency between rebuild and plot stages.

### Files Modified
- `/Users/dengxiaoya/TASK/latest/task/eq/rebuild_f.sh` (2 changes, 2 new lines added)

### Evidence
- Diff recorded in: `.sisyphus/evidence/task-1-sh-diff.txt`

## T2/T3: Dual q output columns + F_direct comment (2026-03-03)

### Pattern: Initialize output variables before branching, update inside branch
- `q_gfile_out` initialized to `np.full_like(q, np.nan, dtype=float)` before the `if args.method` branch
- Set to `q_gf.copy()` inside `if args.use_gfile_qp:` (nested inside the eqipqp branch)
- This avoids NameError when method="direct" tries to reference q_gfile_out at DataFrame time

### q semantics in this script
- `q` = q(psi) from radial CSV (dr[c_q]), always the "original" radial q
- `q_use` = q actually fed into _rebuild_f_eqipqp(); may be q from EQIPQP inputs, or q from gfile
- `q_gf` = q loaded from gfile profile CSV when --use-gfile-qp is active
- `q_gfile_out` = q_gf.copy() when gfile q used, NaN array otherwise

### Output CSV new columns
- `q_gfile`: NaN unless --use-gfile-qp; then holds gfile q values on psip grid
- `q_source`: "gfile" or "radial" string per row — machine-readable source flag

### F_direct is always q(psi)-based
- The algebraic formula `F = 4π² q / (dV/dψ · <1/R²>)` uses `q` (radial CSV q), never q_use
- Comment added at line 266 to document this invariant explicitly

### File actual location
- The `.sisyphus/` tree is under `/Users/dengxiaoya/TASK/latest/task/eq/` (NOT eq-rebuild-f-gfile)
- The actual script is at `/Users/dengxiaoya/TASK/latest/task/eq/in/rebuild_f_from_q_eqdata.py`
- Task description referenced eq-rebuild-f-gfile path (incorrect); file found in `eq/in/`
