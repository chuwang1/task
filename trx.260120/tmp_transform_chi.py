import os
import numpy as np
import pandas as pd

base_tr = '/Users/dengxiaoya/TASK/latest/task/trx'
base_om = '/Users/dengxiaoya/CFEDRSW/OMFIT_out/outputs_new'

chi = pd.read_csv(os.path.join(base_tr, 'chi.dat'), sep=r'\s+')
rho_chi = chi['rho'].values

tr_ar2 = pd.read_csv(os.path.join(base_tr, 'tr_data_137.csv'), skiprows=1)
rho_ar2 = tr_ar2['X'].values
ar2 = tr_ar2['AR2RHO'].values
ar2_on_chi = np.interp(rho_chi, rho_ar2, ar2)

allp = pd.read_csv(os.path.join(base_om, 'all_profiles.csv'))
rho_om = allp['rho'].values
rmin = allp['rmin'].values
drdrho = np.gradient(rmin, rho_om)
drdrho_on_chi = np.interp(rho_chi, rho_om, drdrho)

chi2 = chi.copy()
chi2['chi_e'] = chi['chi_e'].values / (ar2_on_chi * drdrho_on_chi)
chi2['chi_i'] = chi['chi_i'].values / (ar2_on_chi * drdrho_on_chi)

out = os.path.join(base_tr, 'chi_div_grad2_drdrho.dat')
chi2.to_csv(out, sep=' ', index=False, float_format='%.8e')

print('Wrote', out)
print('Original chi_e range:', chi['chi_e'].min(), chi['chi_e'].max())
print('Transformed chi_e range:', chi2['chi_e'].min(), chi2['chi_e'].max())
print('Original chi_i range:', chi['chi_i'].min(), chi['chi_i'].max())
print('Transformed chi_i range:', chi2['chi_i'].min(), chi2['chi_i'].max())

for idx in [0, 20, 40, 80, 120, 160, 199]:
    print(
        f"rho={rho_chi[idx]:.3f} "
        f"ar2={ar2_on_chi[idx]:.4f} drdrho={drdrho_on_chi[idx]:.3f} "
        f"chi_e:{chi.iloc[idx]['chi_e']:.4f}->{chi2.iloc[idx]['chi_e']:.4f} "
        f"chi_i:{chi.iloc[idx]['chi_i']:.4f}->{chi2.iloc[idx]['chi_i']:.4f}"
    )
