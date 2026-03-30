./eq <in/eq.CFEDR0114.in
python3 read_eqdata_to_csv.py eqdata0114
python3 export_2d_from_eqdata.py eqdata0114
cd in
python3 compare_psi_3way.py --prefix eqdata0114 --eq-psi-csv eqdata0114_eq_psirz_gfile.csv 
