import numpy as np
import pandas as pd
import f90nml
import os
from scipy.interpolate import interp1d

base = '/Users/dengxiaoya/TASK/latest/task/trx'
os.chdir(base)


def read_tr_csv(filename):
    with open(filename, 'r') as f:
        _ = f.readline()
    df = pd.read_csv(filename, skiprows=1, index_col=False)
    df = df.loc[:, ~df.columns.str.contains('^Unnamed')]
    x = df.iloc[:, 0].values
    data = {col: df[col].values for col in df.columns[1:]}
    return x, data


def read_omfit_profiles(namelist_file='profiles_CFEDR.namelist'):
    nml = f90nml.read(namelist_file)
    pwo = nml['powers_particle_flux_onetwo']
    return {
        'rho': np.array(pwo['rho']),
        'qrad': np.array(pwo['qrad']) / 1e6,
    }


def load_rho_to_ra_mapping(csv_file='~/CFEDRSW/OMFIT_out/rho_rmin.csv'):
    csv_file = os.path.expanduser(csv_file)
    if not os.path.exists(csv_file):
        if os.path.exists('rho_rmin.csv'):
            csv_file = 'rho_rmin.csv'
        else:
            return None
    df = pd.read_csv(csv_file)
    return df['rho'].values, df['rho'].values


def interpolate_to_tr_grid(omfit_data, rho_omfit, r_tr, rho_rmin=None):
    if rho_rmin is not None:
        rho_map, rmin_map = rho_rmin
        a = np.max(rmin_map)
        r_over_a_map = rmin_map / a
        r_over_a_omfit = np.interp(rho_omfit, rho_map, r_over_a_map)
    else:
        r_over_a_omfit = rho_omfit
    f = interp1d(r_over_a_omfit, omfit_data, kind='linear', bounds_error=False, fill_value='extrapolate')
    return f(r_tr)


r_tr, power_src = read_tr_csv('tr_data_159.csv')
_, rad_exch = read_tr_csv('tr_data_160.csv')
prsum_tr = rad_exch.get('PRSUM', power_src.get('-PRSUM', np.zeros_like(r_tr)))
prl_tr = rad_exch.get('PRL', np.zeros_like(r_tr))
prb_tr = rad_exch.get('PRB', np.zeros_like(r_tr))
prc_tr = rad_exch.get('PRC', np.zeros_like(r_tr))

omfit = read_omfit_profiles('profiles_CFEDR.namelist')
rho_rmin = load_rho_to_ra_mapping()
qrad_omfit = interpolate_to_tr_grid(omfit['qrad'], omfit['rho'], r_tr, rho_rmin)

delta = qrad_omfit - prsum_tr
prl_new = prl_tr + delta

pd.DataFrame({'rho': r_tr, 'prl': prl_new}).to_csv('omfit_prl_for_tr.csv', index=False)
pd.DataFrame({
    'rho': r_tr,
    'prl_old': prl_tr,
    'delta_qrad_minus_prsum': delta,
    'prl_new': prl_new,
    'prsum_tr': prsum_tr,
    'prb_tr': prb_tr,
    'prc_tr': prc_tr,
    'qrad_omfit': qrad_omfit,
}).to_csv('omfit_prl_for_tr_debug.csv', index=False)

print('Wrote', os.path.join(base, 'omfit_prl_for_tr.csv'))
print('Debug ', os.path.join(base, 'omfit_prl_for_tr_debug.csv'))
print('PRL old range   :', float(prl_tr.min()), float(prl_tr.max()))
print('Delta range     :', float(delta.min()), float(delta.max()))
print('PRL new range   :', float(prl_new.min()), float(prl_new.max()))
print('Any negative?   :', bool((prl_new < 0).any()))
