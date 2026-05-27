rm -rf out/tr_data_*
./tr2 <$1    
# ./tr2  <tr.iter01.in
/home/chuw/ai-venv/bin/python  parse_tr_csv.py tr_data_*.csv -n 20
/home/chuw/ai-venv/bin/python images_to_pdf.py 4
/home/chuw/ai-venv/bin/python analyze_energy_balance_final.py 
/home/chuw/ai-venv/bin/python calculate_temperature_profile.py --save-csv 2>&1
/home/chuw/ai-venv/bin/python compare_temperature_omfit.py 
/home/chuw/ai-venv/bin/python compare_omfit_tr_sources.py
/home/chuw/ai-venv/bin/python compare_jtot_jbs_v3.py
/home/chuw/ai-venv/bin/python verify_energy_balance.py
/home/chuw/ai-venv/bin/python calculate_cdbm_chi.py
/home/chuw/ai-venv/bin/python compare_chimix_external_vs_tr.py
cp ./tr_data_* out/
# cp out/tr_data_017.csv ./PyMak/data/density.csv 
# cp out/tr_data_019.csv ./PyMak/data/temperature.csv 
# cp out/tr_data_022.csv ./PyMak/data/RADIATION.csv  
# cp out/tr_data_033.csv ./PyMak/data/PNFIN.csv
# cp out/tr_data_033.csv ./PyMak/data/PNFCL.csv
# cp out/tr_data_021.csv ./PyMak/data/power.csv