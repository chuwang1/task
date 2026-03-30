#!/usr/bin/env bash
set -euo pipefail

myprefix="${1:?Usage: bash rebuild_f.sh <prefix>}"

make eq
./eq < eq_export_samegeom.in

eqipqp_qps_csv="in/${myprefix}_eqipqp_inputs_qps.csv"
conda run -n Fusion python3 -c "import numpy as np, pandas as pd; deq=pd.read_csv('eqgs1d_24_EQIPQP_INPUTS.csv'); dq=pd.read_csv('eqgs1d_04_QPS.csv'); deq['QPSI_ORIG']=deq['QPSI']; deq['QPSI']=np.interp(deq['PSIPV'].to_numpy(float), dq['PSIP'].to_numpy(float), dq['QPS'].to_numpy(float)); deq.to_csv('${eqipqp_qps_csv}', index=False)"

conda run -n Fusion python3 in/rebuild_f_from_q_eqdata.py \
  --prefix "$myprefix" \
  --eqipqp-input-csv "$eqipqp_qps_csv" \
  --method eqipqp \
  --out-prefix "in/${myprefix}_f_from_q_eqipqp_stable"

conda run -n Fusion python3 in/plot_eqipqp_all_variables_big.py \
  --prefix "$myprefix" \
  --eqipqp-input-csv "$eqipqp_qps_csv" \
  --rebuild-csv "in/${myprefix}_f_from_q_eqipqp_stable.csv" \
  --out "in/${myprefix}_eqipqp_all_variables_big.png"
