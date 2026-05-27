#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_DIR=$(cd -- "${SCRIPT_DIR}/../.." && pwd)
RUN_DIR="${SCRIPT_DIR}/tr_output_3step"
PYTHON_BIN="${PYTHON:-${HOME}/ai-venv/bin/python}"

mkdir -p "${RUN_DIR}/gs" "${SCRIPT_DIR}/results"

ln -sfn "${REPO_DIR}/tr2" "${RUN_DIR}/tr2"
ln -sfn "${REPO_DIR}/density_input0114_rho.csv" "${RUN_DIR}/density_input0114_rho.csv"
ln -sfn "${REPO_DIR}/omfit_prl_for_tr.csv" "${RUN_DIR}/omfit_prl_for_tr.csv"
ln -sfn "${REPO_DIR}/chi_div_grad2_drdrho.dat" "${RUN_DIR}/chi_div_grad2_drdrho.dat"

(
  cd "${RUN_DIR}"
  ./tr2 < tr_onetwo_3step.in > run_3step.log 2>&1
)

"${PYTHON_BIN}" "${SCRIPT_DIR}/compare_tr_output_prc.py" \
  --run-dir "${RUN_DIR}" \
  --gfile "${REPO_DIR}/../eq/in/g260206.20000_teq_0114" \
  --out-prefix "${SCRIPT_DIR}/results/tr_output_3step_prc_comparison"
