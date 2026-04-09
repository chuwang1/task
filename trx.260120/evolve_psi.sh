python3 build_gs_profile_csv.py \
  --gfile g260206.20000_teq_0114 \
  --pressure-mode tr_total \
  --f-mode ttrho \
  --output offline_profile_psin_pfast_ttrhoF.csv
python3 "/Users/dengxiaoya/TASK/CFEDR/git/task/eq/solve_gs_from_gfile.py" --gfile "/Users/dengxiaoya/TASK/CFEDR/git/task/trx.260120/g260206.20000_teq_0114" --profile-csv "/Users/dengxiaoya/TASK/CFEDR/git/task/trx.260120/offline_profile_psin_pfast_ttrhoF.csv" --profile-x psi_n --out-prefix "/Users/dengxiaoya/TASK/CFEDR/git/task/trx.260120/offline_currentsolver_pfast_ttrhoF_retest" --write-gfile "/Users/dengxiaoya/TASK/CFEDR/git/task/trx.260120/offline_currentsolver_pfast_ttrhoF_retest.gfile" --freeze-source --max-outer 120
