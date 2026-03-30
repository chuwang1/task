rm -rf ./tr_data_*
./tr2 <$1    
# ./tr2  <tr.iter01.in
/Users/dengxiaoya/miniforge3/envs/Fusion/bin/python  parse_tr_csv.py tr_data_*.csv -n 20
/Users/dengxiaoya/miniforge3/envs/Fusion/bin/python images_to_pdf.py 4
/Users/dengxiaoya/miniforge3/envs/Fusion/bin/python analyze_energy_balance_final.py 
/Users/dengxiaoya/miniforge3/envs/Fusion/bin/python calculate_temperature_profile.py --save-csv 2>&1
/Users/dengxiaoya/miniforge3/envs/Fusion/bin/python compare_temperature_omfit.py 
/Users/dengxiaoya/miniforge3/envs/Fusion/bin/python compare_omfit_tr_sources.py
/Users/dengxiaoya/miniforge3/envs/Fusion/bin/python compare_jtot_jbs_v3.py
/Users/dengxiaoya/miniforge3/envs/Fusion/bin/python verify_energy_balance.py
/Users/dengxiaoya/miniforge3/envs/Fusion/bin/python calculate_cdbm_chi.py
/Users/dengxiaoya/miniforge3/envs/Fusion/bin/python compare_chimix_external_vs_tr.py