#!/usr/bin/env python3
"""
Analyze geometry factors used in TR code for temperature equation.

Based on trmetric.f90:
    DVRHO(NR) = 2*PI*RKAP*RA²*2*PI*RR*RM(NR)   ! dV/dr
    AR2RHO(NR) = 1/(RKAPS*RA)²                  ! <|∇ρ|²>
    AR1RHO(NR) = 1/(RKAPS*RA)                   ! <|∇ρ|>

where RKAPS = sqrt(RKAP), RM(NR) = r/a

The diffusion equation uses:
    FB = DVRHO * AR2RHO / DR²
    DD = FB * chi * DV23
where DV23 = DVRHO^(2/3), DV53 = DVRHO^(5/3)
"""

import numpy as np
import matplotlib.pyplot as plt

# CFEDR parameters
R0 = 8.03        # Major radius [m]
A_MINOR = 2.68   # Minor radius [m]
KAPPA = 1.89     # Elongation

def calculate_geometry_factors(r_a, a=A_MINOR, R0=R0, kappa=KAPPA):
    """
    Calculate TR geometry factors.

    Args:
        r_a: r/a array (normalized radius)
        a: minor radius [m]
        R0: major radius [m]
        kappa: elongation

    Returns:
        dict with geometry factors
    """
    kappa_s = np.sqrt(kappa)  # RKAPS
    r = r_a * a  # r in meters

    # From trmetric.f90
    DVRHO = 2 * np.pi * kappa * a**2 * 2 * np.pi * R0 * r_a
    AR2RHO = 1 / (kappa_s * a)**2
    AR1RHO = 1 / (kappa_s * a)

    # DVRHO * AR2RHO (should simplify)
    DVRHO_AR2 = DVRHO * AR2RHO

    # Simple cylindrical equivalent
    DVRHO_cyl = 4 * np.pi**2 * R0 * r

    # Volume powers
    DV23 = DVRHO**(2/3)
    DV53 = DVRHO**(5/3)

    return {
        'r_a': r_a,
        'r': r,
        'DVRHO': DVRHO,
        'AR2RHO': AR2RHO,
        'AR1RHO': AR1RHO,
        'DVRHO_AR2': DVRHO_AR2,
        'DVRHO_cyl': DVRHO_cyl,
        'DV23': DV23,
        'DV53': DV53,
        'kappa': kappa,
        'kappa_s': kappa_s,
    }


def analyze_diffusion_coefficient(geom, chi, n, DR):
    """
    Analyze effective diffusion coefficient as used in TR.

    DD = FB * chi * DV23
    where FB = DVRHO * AR2RHO / DR²

    Args:
        geom: geometry factors dict
        chi: thermal diffusivity [m²/s]
        n: density [10^20/m³]
        DR: radial step [m]

    Returns:
        dict with diffusion analysis
    """
    FB = geom['DVRHO'] * geom['AR2RHO'] / DR**2
    DD = FB * chi * geom['DV23']

    # Compare with simple cylindrical form
    # Simple: D = (r * n * chi * k) / dr
    # where k = 1.6e4 (unit conversion for n in 10^20, T in keV)
    e_factor = 1.6e4
    D_simple = geom['r'] * n * chi * e_factor / DR

    return {
        'FB': FB,
        'DD': DD,
        'D_simple': D_simple,
        'ratio': DD / (D_simple + 1e-30),
    }


def main():
    print("="*70)
    print("GEOMETRY FACTOR ANALYSIS")
    print("="*70)

    # Create r/a grid
    nr = 50
    r_a = np.linspace(0.01, 0.99, nr)
    DR = (r_a[1] - r_a[0]) * A_MINOR

    # Calculate geometry factors
    geom = calculate_geometry_factors(r_a)

    print(f"\nCFEDR Parameters:")
    print(f"  R0 = {R0:.2f} m")
    print(f"  a  = {A_MINOR:.2f} m")
    print(f"  κ  = {KAPPA:.2f}")
    print(f"  √κ = {geom['kappa_s']:.4f}")
    print(f"  DR = {DR:.4f} m")

    print(f"\nGeometry Factors (at r/a = 0.5):")
    idx = nr // 2
    print(f"  DVRHO     = {geom['DVRHO'][idx]:.4e}")
    print(f"  AR2RHO    = {geom['AR2RHO']:.4e}")
    print(f"  AR1RHO    = {geom['AR1RHO']:.4e}")
    print(f"  DVRHO*AR2 = {geom['DVRHO_AR2'][idx]:.4e}")
    print(f"  DVRHO_cyl = {geom['DVRHO_cyl'][idx]:.4e}")
    print(f"  Ratio     = {geom['DVRHO_AR2'][idx] / geom['DVRHO_cyl'][idx]:.4f}")

    # Analyze effective diffusion
    print(f"\nEffective Diffusion Analysis (chi=1, n=1):")
    diff = analyze_diffusion_coefficient(geom, chi=1.0, n=1.0, DR=DR)
    print(f"  FB (at r/a=0.5)       = {diff['FB'][idx]:.4e}")
    print(f"  DD (at r/a=0.5)       = {diff['DD'][idx]:.4e}")
    print(f"  D_simple (at r/a=0.5) = {diff['D_simple'][idx]:.4e}")
    print(f"  Ratio DD/D_simple     = {diff['ratio'][idx]:.4f}")

    # Geometry correction factor
    # My simple model: q = -n * chi * k * dT/dr
    # TR model: q = -DVRHO * AR2RHO / (DR²) * chi * DV23 * dT
    # The ratio gives the geometry correction

    print(f"\nGeometry Correction Factor:")
    print(f"  Expected from elongation: κ = {KAPPA:.2f}")
    print(f"  Calculated ratio: {diff['ratio'][idx]:.4f}")
    print(f"  This factor modifies the effective chi")

    # Plot
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    ax = axes[0, 0]
    ax.plot(r_a, geom['DVRHO'], 'b-', lw=2, label='DVRHO (elongated)')
    ax.plot(r_a, geom['DVRHO_cyl'], 'r--', lw=2, label='4π²R₀r (cylindrical)')
    ax.set_xlabel('r/a')
    ax.set_ylabel('dV/dr [m²]')
    ax.set_title('Volume Element')
    ax.legend()
    ax.grid(True)

    ax = axes[0, 1]
    ax.plot(r_a, geom['DVRHO_AR2'], 'b-', lw=2, label='DVRHO × AR2RHO')
    ax.plot(r_a, geom['DVRHO_cyl'], 'r--', lw=2, label='4π²R₀r')
    ax.set_xlabel('r/a')
    ax.set_ylabel('[m²/m²]')
    ax.set_title('DVRHO × AR2RHO Comparison')
    ax.legend()
    ax.grid(True)

    ax = axes[1, 0]
    ax.plot(r_a, diff['DD'], 'b-', lw=2, label='DD (TR form)')
    ax.plot(r_a, diff['D_simple'], 'r--', lw=2, label='D_simple')
    ax.set_xlabel('r/a')
    ax.set_ylabel('Diffusion coefficient')
    ax.set_title('Diffusion Coefficient Comparison (chi=1, n=1)')
    ax.legend()
    ax.grid(True)

    ax = axes[1, 1]
    ax.plot(r_a, diff['ratio'], 'b-', lw=2)
    ax.axhline(1.0, color='k', ls='--', lw=1)
    ax.set_xlabel('r/a')
    ax.set_ylabel('DD / D_simple')
    ax.set_title('Geometry Correction Ratio')
    ax.grid(True)

    plt.tight_layout()
    plt.savefig('geometry_factors_analysis.png', dpi=150)
    print(f"\nPlot saved: geometry_factors_analysis.png")

    # Summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"""
The TR code uses geometry factors from elongated flux surface coordinates:
  - DVRHO = 4π²κa²R₀(r/a) - volume element including elongation
  - AR2RHO = 1/(√κ × a)² - gradient correction

When combined: DVRHO × AR2RHO = 4π²R₀r × (1/κ)

For CFEDR (κ={KAPPA}):
  - The elongation reduces effective diffusion by factor ~1/κ ≈ {1/KAPPA:.3f}
  - This means effective chi is LOWER than expected from simple cylindrical model
  - Temperature should be HIGHER with elongation correction

Current discrepancy:
  - Electron: need chi × 1.32 → suggests geometry factor helps
  - Ion: need chi × 0.73 → opposite direction, NOT explained by geometry
""")

    return geom, diff


if __name__ == '__main__':
    main()
