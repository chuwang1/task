#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

GFILE="${1:-in/g260206.20000_teq_0114}"
PREFIX="${2:-eqdata_modelg5_rebuild_mdleqf9_qmap_samegeom}"
CONDA_ENV="${CONDA_ENV:-Fusion}"
USE_GFILE_LCFS="${USE_GFILE_LCFS:-1}"
INIT_PSI_MODE="${INIT_PSI_MODE:-1}"

echo "== [1/4] Build m5 =="
make m5

echo "== [2/4] Run m5 =="
echo "   gfile : $GFILE"
echo "   prefix: $PREFIX"
echo "   lcfs  : $USE_GFILE_LCFS (1=gfile, 0=EQFBND)"
echo "   init  : $INIT_PSI_MODE (1=gfile PSIRZ, 0=analytic)"
rm -f "$PREFIX"
printf '0\nc\n' | ./m5 "$GFILE" "$PREFIX" "-" "$USE_GFILE_LCFS" "$INIT_PSI_MODE"
if [[ ! -s "$PREFIX" ]]; then
  echo "ERROR: m5 did not produce output file: $PREFIX" >&2
  exit 1
fi

echo "== [3/4] Export eqdata CSV =="
conda run -n "$CONDA_ENV" python3 read_eqdata_to_csv.py "$PREFIX"
conda run -n "$CONDA_ENV" python3 export_2d_from_eqdata.py "$PREFIX"

echo "== [4/4] Plot comparison =="
GBASE="$(basename "$GFILE")"
PROFILE_CSV="$SCRIPT_DIR/in/${GBASE}_profiles_from_gfile.csv"
COMPARE_ARGS=(--prefix "$PREFIX")
if [[ -f "$PROFILE_CSV" ]]; then
  COMPARE_ARGS+=(--gfile-profile-csv "$PROFILE_CSV")
fi

(
  cd "$SCRIPT_DIR/in"
  conda run -n "$CONDA_ENV" python3 compare_psi_3way.py "${COMPARE_ARGS[@]}"
)

echo ""
echo "Done."
echo "  eqdata file : $SCRIPT_DIR/$PREFIX"
echo "  overview fig: $SCRIPT_DIR/in/comparison_overview_big.png"
echo "  profile fig : $SCRIPT_DIR/in/profile_q_p_f_dp_ff_comparison.png"
