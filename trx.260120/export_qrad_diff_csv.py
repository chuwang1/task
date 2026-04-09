#!/usr/bin/env python3

import glob
import os

import netCDF4 as nc4
import numpy as np
import pandas as pd


OMFIT_DIR = '/Users/dengxiaoya/CFEDRSW/OMFIT_out'
TR_DIR = '/Users/dengxiaoya/TASK/CFEDR/git/task/trx.260120'
STATEFILE = os.path.join(OMFIT_DIR, 'statefile_3.000000E+01.nc')
OUTPUT_CSV = os.path.join(TR_DIR, 'qrad_diff.csv')


def find_csv_by_title(title_keyword, search_dir=TR_DIR):
    for path in sorted(glob.glob(os.path.join(search_dir, 'tr_data_*.csv'))):
        with open(path, 'r') as fh:
            title = fh.readline().strip()
        if title_keyword in title:
            return path
    return None


def read_tr_csv(path):
    df = pd.read_csv(path, skiprows=1)
    df.columns = df.columns.str.strip()
    return df


def main():
    tr_csv = find_csv_by_title('PRSUM,PRB,PRC,PRL,PCX,PIE,QEI')
    if tr_csv is None:
        raise FileNotFoundError('Could not find TR CSV containing PRSUM/PRB/PRC/PRL/PCX/PIE')

    df_tr = read_tr_csv(tr_csv)
    r_a_tr_rm = df_tr['X'].values
    dr_tr = float(np.median(np.diff(r_a_tr_rm))) if len(r_a_tr_rm) > 1 else 0.02
    r_a_tr = np.clip(r_a_tr_rm + 0.5 * dr_tr, 0.0, 1.0)
    tr_prsum = df_tr['PRSUM'].values

    with nc4.Dataset(STATEFILE) as ds:
        rho_grid = ds.variables['rho_grid'][:]
        r_a_omfit = rho_grid / rho_grid.max()
        qrad = ds.variables['qrad'][:] * 1e-6

    diff_rad = np.interp(r_a_tr, r_a_omfit, qrad) - tr_prsum
    pd.DataFrame({'r_a_omfit': r_a_tr, 'prl': diff_rad}).to_csv(OUTPUT_CSV, index=False)

    print(f'TR csv: {os.path.basename(tr_csv)}')
    print(f'Statefile: {STATEFILE}')
    print(f'Saved: {OUTPUT_CSV}')


if __name__ == '__main__':
    main()
