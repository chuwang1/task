#!/usr/bin/env python3
"""
Convert density_input.txt to TASK/TR format files
Supports: ne, nD, nT, nHe profiles
"""
import numpy as np
from scipy.optimize import curve_fit

# Read CSV data: rho, ne, nD, nT, nHe (all in 10^20 m^-3)
data = np.loadtxt('../tr/in/density_input.txt', delimiter=',', skiprows=1)
rho = data[:, 0]
ne = data[:, 1]
nD = data[:, 2]
nT = data[:, 3]
nHe = data[:, 4]

# === Output 1: Profile file with all species ===
# Format: rho, ne, nD, nT, nHe (m^-3)
with open('density_prof.dat', 'w') as f:
    f.write('rho  ne(m^-3)  nD(m^-3)  nT(m^-3)  nHe(m^-3)\n')
    for i in range(len(rho)):
        f.write(f'{rho[i]:.6f}  {ne[i]*1e20:.6e}  {nD[i]*1e20:.6e}  {nT[i]*1e20:.6e}  {nHe[i]*1e20:.6e}\n')
print('Created: density_prof.dat (with all species)')

# Calculate density ratios at rho=0 for PN settings
print(f'\n=== Density ratios at rho=0 (for PN settings) ===')
print(f'PN(1) = ne(0) = {ne[0]:.4f} (10^20 m^-3)')
print(f'PN(2) = nD(0) = {nD[0]:.4f}')
print(f'PN(3) = nT(0) = {nT[0]:.4f}')
print(f'PN(4) = nHe(0) = {nHe[0]:.4f}')
print(f'\nRatio nD/ne = {nD[0]/ne[0]:.4f}')
print(f'Ratio nT/ne = {nT[0]/ne[0]:.4f}')
print(f'Ratio nHe/ne = {nHe[0]/ne[0]:.4f}')

# === Output 2: Coefficient file for model_nfixed ===
# Fit to: n(rho) = c0 + c4*(1-rho^2)^c5

def profile_func(rho, c0, c4, c5):
    """Simplified profile function"""
    return c0 + c4 * (1 - rho**2)**c5

# Initial guess
p0 = [ne[-1], ne[0] - ne[-1], 1.0]

# Fit
try:
    popt, _ = curve_fit(profile_func, rho, ne, p0=p0, maxfev=10000)
    c0, c4, c5 = popt

    # Create coefficient file
    # Format: ntime, ncoef, rho_min, rho_max
    #         time
    #         coef(0), coef(1), ..., coef(10)
    with open('nprof_coef_data', 'w') as f:
        f.write('1  10  0.0  1.0\n')  # 1 time point, 10 coefficients
        f.write('0.0\n')  # start time = 0
        # coef: c0, c1=0, c2=0, c3=1, c4, c5, c6=0, c7=0, c8=0, c9=0, c10=1
        f.write(f'{c0:.6e}  0.0  0.0  1.0  {c4:.6e}  {c5:.6e}  0.0  0.0  0.0  0.0  1.0\n')

    print(f'Created: nprof_coef_data (for model_nfixed=1)')
    print(f'  Fitted coefficients: c0={c0:.4f}, c4={c4:.4f}, c5={c5:.4f}')

    # Check fit quality
    ne_fit = profile_func(rho, *popt)
    rms_error = np.sqrt(np.mean((ne - ne_fit)**2))
    print(f'  RMS error: {rms_error:.4f} (10^20 m^-3)')
    print(f'  Relative error: {rms_error/np.mean(ne)*100:.2f}%')

except Exception as e:
    print(f'Fitting failed: {e}')
    print('You may need to manually create the coefficient file')

print('\n=== Usage in tr.iter.in ===')
print('model_prof=11')
print("knam_prof='density_prof.dat'")
print('model_nfixed=1')
print("knam_nfixed='nprof_coef_data'")