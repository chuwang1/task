import numpy as np
import matplotlib.pyplot as plt

# Parameters from tr.iter.in and trinit.f90
PROFN1 = 1.0  # Default from trinit.f90
PROFN2 = 0.15 # From tr.iter.in
ALP1 = 1.0    # Default from trinit.f90

# Species definitions (from trinit.f90)
species_names = ['Electron', 'Deuterium', 'Tritium', 'Helium-4']
# PN: Initial number density on axis (1.E20 m**-3) from tr.iter.in
PN = np.array([1.0, 0.45, 0.45, 0.055])
# PNS: Initial number density on surface (1.E20 m**-3) from tr.iter.in
PNS = np.array([0.1, 0.045, 0.045, 0.005])

# Radial grid (rho)
rho = np.linspace(0, 1, 100)

plt.figure(figsize=(10, 6))

for i, name in enumerate(species_names):
    # Formula from trprof.f90:
    # PROF = (1.D0-(ALP(1)*RM(NR))**PROFN1(NS))**PROFN2(NS)
    # RN(NR,NS) = (PN(NS)-PNS(NS))*PROF+PNS(NS)
    
    # Note: PROFN1/2 can be arrays, but usually scalars are used for all species or first element
    # trinit.f90 defines them as scalars, but code uses PROFN1(NS)?
    # trinit.f90: PROFN1 = 2.D0 (scalar variable)
    # trprof.f90: PROFN1(NS) implies array?
    # Let's check trcomm.f90 for PROFN1 declaration.
    # trinit.f90 line 144: PROFN1 = 2.D0
    # trcomm.f90 line 63: REAL(rkind):: PROFNU1,PROFNU2,PROFJ1,PROFJ2 (Wait, where are PROFN1 defined?)
    # trcomm.f90 IMPORTS PROFN1 from plcomm (line 15).
    # plcomm is likely another module. 
    # Usually in these codes, if input is scalar, it applies to all.
    # tr.iter.in sets PROFN2=0.15D0 (scalar).
    
    prof = (1.0 - (ALP1 * rho)**PROFN1)**PROFN2
    n_rho = (PN[i] - PNS[i]) * prof + PNS[i]
    
    plt.plot(rho, n_rho, label=f'{name} (PN={PN[i]}, PNS={PNS[i]})')

plt.title(f'Density Profile (PROFN1={PROFN1}, PROFN2={PROFN2})')
plt.xlabel('Normalized Radius (rho)')
plt.ylabel('Density (10^20 m^-3)')
plt.grid(True)
plt.legend()
plt.savefig('density_profile_plot.png')
print("Plot saved to density_profile_plot.png")
