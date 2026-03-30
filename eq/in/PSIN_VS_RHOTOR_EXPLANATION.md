# Radial Coordinates in Tokamak Physics: $\Psi_n$ vs $\rho_{tor}$

In tokamak plasma physics, 1D radial profiles of quantities like temperature, density, and safety factor are mapped to magnetic flux surfaces. Because these flux surfaces are not perfectly circular or concentric, we use flux-surface labels as "radial" coordinates to simplify complex 3D geometry into 1D profiles. 

Two of the most common normalized radial coordinates are **$\Psi_n$** (normalized poloidal flux) and **$\rho_{tor}$** (normalized toroidal flux radius). This document explains their definitions, physical meanings, and specific use cases.

---

## 1. Normalized Poloidal Flux ($\Psi_n$)

The poloidal magnetic flux, $\Psi$, is the total magnetic flux passing through a horizontal ribbon extending from the magnetic axis outward to a given flux surface. 

**$\Psi_n$** (often just called `psi_n` or `psin`) is the poloidal flux normalized such that it equals $0$ at the magnetic axis and $1$ at the Last Closed Flux Surface (LCFS) or separatrix.

### Mathematical Definition
$$ \Psi_n = \frac{\Psi - \Psi_{axis}}{\Psi_{LCFS} - \Psi_{axis}} $$

*Where:*
*   $\Psi$ is the poloidal flux at a given surface.
*   $\Psi_{axis}$ is the poloidal flux at the magnetic axis.
*   $\Psi_{LCFS}$ is the poloidal flux at the plasma boundary.

### Characteristics
*   **Relationship to Radius:** Near the magnetic axis, $\Psi_n$ scales quadratically with the physical minor radius $r$ (i.e., $\Psi_n \propto r^2$).
*   **Primary Use Case:** $\Psi_n$ is the natural coordinate for Magnetohydrodynamic (MHD) equilibrium. The Grad-Shafranov equation, which dictates the macroscopic force balance of the plasma, uses poloidal flux as its fundamental independent variable. Consequently, equilibrium codes (like EFIT) natively output profiles and geometry on a $\Psi_n$ grid.

---

## 2. Normalized Toroidal Flux Radius ($\rho_{tor}$)

The toroidal magnetic flux, $\Phi$, is the total magnetic flux passing perpendicularly through the poloidal cross-sectional area bounded by a specific flux surface. 

**$\rho_{tor}$** (often called `rho_tor` or `rho`) is defined as the square root of the normalized toroidal magnetic flux. Taking the square root gives it the dimension of a length, making it behave similarly to a true geometric radius.

### Mathematical Definition
First, we define the normalized toroidal flux, $\Phi_n$:
$$ \Phi_n = \frac{\Phi - \Phi_{axis}}{\Phi_{LCFS} - \Phi_{axis}} $$

Then, $\rho_{tor}$ is simply the square root of $\Phi_n$:
$$ \rho_{tor} = \sqrt{\Phi_n} = \sqrt{\frac{\Phi - \Phi_{axis}}{\Phi_{LCFS} - \Phi_{axis}}} $$

### Characteristics
*   **Relationship to Radius:** Because toroidal flux $\Phi \propto r^2$ for circular cross-sections, taking the square root means $\rho_{tor}$ scales linearly with the physical minor radius $r$ near the core ($\rho_{tor} \propto r$).
*   **Primary Use Case:** $\rho_{tor}$ is the preferred coordinate for **transport modeling** and **gyrokinetic simulations** (e.g., in codes like TRANSP, GENE, CGYRO, or TGLF). Because it scales linearly with distance, gradients calculated with respect to $\rho_{tor}$ more closely reflect true physical spatial gradients ($\nabla T, \nabla n$), which are the fundamental drivers for turbulence and heat transport.

---

## Comparison Table

| Feature | $\Psi_n$ (Normalized Poloidal Flux) | $\rho_{tor}$ (Normalized Toroidal Radius) |
| :--- | :--- | :--- |
| **Base Quantity** | Poloidal magnetic flux ($\Psi$) | Toroidal magnetic flux ($\Phi$) |
| **Formula** | $\frac{\Psi - \Psi_{axis}}{\Psi_{LCFS} - \Psi_{axis}}$ | $\sqrt{\frac{\Phi - \Phi_{axis}}{\Phi_{LCFS} - \Phi_{axis}}}$ |
| **Scaling near axis** | $\propto r^2$ (Quadratic) | $\propto r$ (Linear) |
| **Physical Analogy** | "Area-like" or "Volume-like" coordinate | "Radius-like" or "Length-like" coordinate |
| **Primary Domain** | MHD, Equilibrium (Grad-Shafranov), EFIT | Transport, Turbulence, Gyrokinetics |
| **Gradient representation** | Skewed near the core (artificially steep) | Accurately represents physical radial gradients |
| **Coordinate link** | $q = \frac{d\Phi}{d\Psi}$ (Safety factor definition) | Safety factor relates $\Psi_n$ and $\rho_{tor}$ grids |

---

## Summary

While $\Psi_n$ and $\rho_{tor}$ are inextricably linked by the safety factor profile $q$ (since $q = d\Phi/d\Psi$), choosing the right coordinate for the right task is essential in fusion research. 

*   **$\Psi_n$** is best when establishing the global force-balance, drawing magnetic geometry, or solving equilibrium problems. 
*   **$\rho_{tor}$** is strictly necessary when calculating gradient-driven heat/particle fluxes and analyzing micro-turbulence, because an "area-like" coordinate like $\Psi_n$ will artificially distort radial gradient calculations near the magnetic axis.

---

## 4. Graphical Verification

The following plot visually verifies the relationship between Normalized Poloidal Flux ($\Psi_n$) and Normalized Toroidal Flux radius ($\rho_{tor}$):

![PSIN vs RHOTOR Mapping](psin_rhotor_mapping_demo.png)