#!/usr/bin/env python3
"""
Read OMFIT profiles namelist and extract species densities and current profiles.
Supports coordinate transformation from rho (sqrt toroidal flux) to r/a.
"""

import f90nml
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os


def build_r_over_a_mapper(rho_tgyro, rmin):
    """
    Build interpolation function to convert rho (sqrt toroidal flux) to r/a.

    Based on TGYRO's rho and rmin (real minor radius) mapping.

    Args:
        rho_tgyro: array of rho values (sqrt normalized toroidal flux)
        rmin: array of rmin values (real minor radius in meters)

    Returns:
        mapper: function that converts rho to r/a
        a: minor radius (max of rmin)
    """
    rho_tgyro = np.array(rho_tgyro, dtype=float)
    rmin = np.array(rmin, dtype=float)

    # Minor radius a
    a = np.max(rmin)

    # Calculate r/a for each TGYRO point
    r_over_a = rmin / a

    def mapper(rho_new):
        """
        Convert rho to r/a via interpolation.

        Args:
            rho_new: rho values to convert (array or scalar)

        Returns:
            r/a values
        """
        rho_new = np.array(rho_new, dtype=float)
        return np.interp(rho_new, rho_tgyro, r_over_a)

    return mapper, a


def load_rho_rmin_mapping(csv_file='~/CFEDRSW/rho_rmin.csv'):
    """Load rho-rmin mapping from CSV file."""
    csv_file = os.path.expanduser(csv_file)
    if not os.path.exists(csv_file):
        print(f"Warning: rho_rmin mapping file not found: {csv_file}")
        return None, None

    df = pd.read_csv(csv_file)
    return df['rho'].values, df['rmin'].values

def read_omfit_profiles(filename='profiles_CFEDR.namelist'):
    """Read OMFIT profiles namelist file."""
    nml = f90nml.read(filename)
    return nml

def read_header_variables(filename='profiles_CFEDR.namelist'):
    """Read variables defined before any namelist block."""
    header_vars = {}
    current_var = None
    current_values = []

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            # Stop when we hit a namelist block
            if line.startswith('&'):
                break

            # Skip empty lines and comments
            if not line or line.startswith('!'):
                continue

            # Check for new variable assignment
            if '=' in line:
                # Save previous variable if exists
                if current_var and current_values:
                    header_vars[current_var] = np.array(current_values)

                parts = line.split('=', 1)
                current_var = parts[0].strip().lower()
                value_str = parts[1].strip()
                current_values = []

                # Parse values
                for val in value_str.split():
                    try:
                        current_values.append(float(val))
                    except ValueError:
                        if val.startswith("'"):
                            header_vars[current_var] = val.strip("'")
                            current_var = None
                            current_values = []
                            break
            else:
                # Continuation of previous variable
                for val in line.split():
                    try:
                        current_values.append(float(val))
                    except ValueError:
                        pass

        # Save last variable
        if current_var and current_values:
            header_vars[current_var] = np.array(current_values)

    return header_vars

def extract_profiles(nml, header_vars=None):
    """Extract key profiles from namelist."""
    profiles = {}

    # Get rho (radial coordinate) - first try header, then namelists
    if header_vars and 'rho' in header_vars:
        profiles['rho'] = header_vars['rho']
    elif 'equlibrium' in nml and 'rho' in nml['equlibrium']:
        profiles['rho'] = np.array(nml['equlibrium']['rho'])

    # Get power and source terms from powers_particle_flux namelists
    power_nml = None
    if 'powers_particle_flux_onetwo' in nml:
        power_nml = nml['powers_particle_flux_onetwo']
    elif 'powers_particle_flux_tgyro' in nml:
        power_nml = nml['powers_particle_flux_tgyro']

    if power_nml:
        # Power flux profiles (MW)
        power_vars = ['pow_e', 'pow_i', 'pow_e_fus', 'pow_i_fus',
                     'pow_e_aux', 'pow_i_aux', 'pow_ei', 'pow_e_rad']
        for var in power_vars:
            if var in power_nml:
                profiles[var] = np.array(power_nml[var])

        # Heating rate profiles (W/m³) - from ONETWO
        heating_vars = ['qrad', 'qfuse', 'qfusi', 'qrfe', 'qrfi',
                       'qbeame', 'qbeami', 'qdelt']
        for var in heating_vars:
            if var in power_nml:
                profiles[var] = np.array(power_nml[var])

        # Particle flux
        flux_vars = ['flow_beam', 'flow_wall', 'flow_beam_and_wall']
        for var in flux_vars:
            if var in power_nml:
                profiles[var] = np.array(power_nml[var])

    # Get electron density
    if 'electron' in nml:
        profiles['ne'] = np.array(nml['electron']['density'])
        profiles['electron_name'] = nml['electron'].get('name', 'e')

    # Get ion densities
    ion_keys = [k for k in nml.keys() if k.startswith('ions_')]
    name_counts = {}  # Track duplicate species names
    for key in sorted(ion_keys):
        ion_data = nml[key]
        name = ion_data.get('name', key)
        density = np.array(ion_data['density'])
        charge = ion_data.get('charge', 1)
        mass = ion_data.get('mass', 1)

        # Handle array-like charge (varies with radius)
        if hasattr(charge, '__len__'):
            charge_avg = np.mean(charge)
        else:
            charge_avg = charge

        # Handle duplicate species names (e.g., thermal and fast alpha)
        if name in name_counts:
            name_counts[name] += 1
            unique_name = f'{name}_{name_counts[name]}'  # e.g., he4_2 for fast alpha
        else:
            name_counts[name] = 1
            unique_name = name

        profiles[f'n_{unique_name}'] = density
        profiles[f'{unique_name}_charge'] = charge_avg
        profiles[f'{unique_name}_mass'] = mass
        print(f"  {key}: name={name} -> n_{unique_name}, mass={mass}, charge_avg={charge_avg:.2f}")

    # Get current densities from header variables
    current_vars = ['current_density', 'bootstrap_current_density_onetwo',
                   'ohmic_current_density_onetwo', 'nbcd_current_density_onetwo',
                   'rfcd_current_density_onetwo']
    for var in current_vars:
        if header_vars and var in header_vars:
            profiles[var] = header_vars[var]
        elif var in nml:
            profiles[var] = np.array(nml[var])

    # Get other parameters from header
    if header_vars:
        if 'r0' in header_vars:
            profiles['R0'] = header_vars['r0'][0] if hasattr(header_vars['r0'], '__len__') else header_vars['r0']
        if 'b0' in header_vars:
            profiles['B0'] = header_vars['b0'][0] if hasattr(header_vars['b0'], '__len__') else header_vars['b0']

    return profiles

def plot_density_profiles(profiles, output_file='omfit_density_profiles.png'):
    """Plot density profiles."""
    # Use r/a if available, otherwise use rho
    if 'r_over_a' in profiles:
        x_coord = profiles['r_over_a']
        x_label = 'r/a'
    else:
        x_coord = profiles.get('rho')
        x_label = 'rho (sqrt normalized toroidal flux)'

    if x_coord is None:
        print("Error: coordinate (rho or r/a) not found")
        return

    # Keep rho for backward compatibility in the function
    rho = x_coord

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # 1. Main species densities (ne, nD, nT)
    ax = axes[0, 0]
    if 'ne' in profiles:
        ax.plot(rho, profiles['ne'] / 1e20, 'b-', label='ne', linewidth=2)
    if 'n_d' in profiles:
        ax.plot(rho, profiles['n_d'] / 1e20, 'r-', label='nD', linewidth=2)
    if 'n_t' in profiles:
        ax.plot(rho, profiles['n_t'] / 1e20, 'g-', label='nT', linewidth=2)
    ax.set_xlabel(x_label)
    ax.set_ylabel('n [10^20 m^-3]')
    ax.set_title('Main Species Density Profiles')
    ax.legend()
    ax.grid(True)
    ax.set_xlim([0, 1])

    # 2. Alpha and impurity densities
    ax = axes[0, 1]
    # Alpha particles (thermal and fast)
    if 'n_he4' in profiles:
        ax.plot(rho, profiles['n_he4'] / 1e18, 'b-', label='nA (thermal)', linewidth=2)
    if 'n_he4_2' in profiles:
        ax.plot(rho, profiles['n_he4_2'] / 1e18, 'c-', label='nA (fast)', linewidth=2)
    if 'n_ar' in profiles:
        ax.plot(rho, profiles['n_ar'] / 1e17, 'm-', label='nAr (x10)', linewidth=2)
    ax.set_xlabel(x_label)
    ax.set_ylabel('n [10^18 m^-3] (Ar: 10^17 m^-3)')
    ax.set_title('Alpha and Impurity Density')
    ax.legend()
    ax.grid(True)
    ax.set_xlim([0, 1])

    # 3. All densities on log scale
    ax = axes[1, 0]
    if 'ne' in profiles:
        ax.semilogy(rho, profiles['ne'], 'b-', label='ne', linewidth=2)
    if 'n_d' in profiles:
        ax.semilogy(rho, profiles['n_d'], 'r-', label='nD', linewidth=2)
    if 'n_t' in profiles:
        ax.semilogy(rho, profiles['n_t'], 'g-', label='nT', linewidth=2)
    if 'n_he4' in profiles:
        ax.semilogy(rho, profiles['n_he4'], 'b--', label='nA (th)', linewidth=2)
    if 'n_he4_2' in profiles:
        ax.semilogy(rho, profiles['n_he4_2'], 'c--', label='nA (fast)', linewidth=2)
    if 'n_ar' in profiles:
        ax.semilogy(rho, profiles['n_ar'], 'm-', label='nAr', linewidth=2)
    ax.set_xlabel(x_label)
    ax.set_ylabel('n [m^-3]')
    ax.set_title('All Densities (log scale)')
    ax.legend()
    ax.grid(True)
    ax.set_xlim([0, 1])

    # 4. Current density profiles
    ax = axes[1, 1]
    def safe_plot_current(ax, rho, data, label, style):
        """Plot current density with dimension check."""
        if len(data) == len(rho):
            ax.plot(rho, data / 1e6, style, label=label, linewidth=2)
        else:
            # Use linear interpolation to rho grid
            x_orig = np.linspace(0, 1, len(data))
            ax.plot(x_orig, data / 1e6, style, label=label, linewidth=2)

    if 'current_density' in profiles:
        safe_plot_current(ax, rho, profiles['current_density'], 'Total J', 'k-')
    if 'bootstrap_current_density_onetwo' in profiles:
        safe_plot_current(ax, rho, profiles['bootstrap_current_density_onetwo'], 'J_BS', 'r--')
    if 'ohmic_current_density_onetwo' in profiles:
        safe_plot_current(ax, rho, profiles['ohmic_current_density_onetwo'], 'J_OH', 'b--')
    if 'nbcd_current_density_onetwo' in profiles:
        j_nb = profiles['nbcd_current_density_onetwo']
        if np.max(np.abs(j_nb)) > 0:
            safe_plot_current(ax, rho, j_nb, 'J_NB', 'g--')
    if 'rfcd_current_density_onetwo' in profiles:
        j_rf = profiles['rfcd_current_density_onetwo']
        if np.max(np.abs(j_rf)) > 0:
            safe_plot_current(ax, rho, j_rf, 'J_RF', 'm--')
    ax.set_xlabel(x_label)
    ax.set_ylabel('J [MA/m^2]')
    ax.set_title('Current Density Profiles')
    if ax.get_legend_handles_labels()[0]:
        ax.legend()
    ax.grid(True)
    ax.set_xlim([0, 1])

    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    print(f"\nPlot saved to {output_file}")

def plot_power_profiles(profiles, output_file='omfit_power_profiles.png'):
    """Plot power flux and heating rate profiles."""
    # Use r/a if available, otherwise use rho
    if 'r_over_a' in profiles:
        x_coord = profiles['r_over_a']
        x_label = 'r/a'
    else:
        x_coord = profiles.get('rho')
        x_label = 'rho'

    if x_coord is None:
        print("Error: coordinate (rho or r/a) not found for power profiles")
        return

    rho = x_coord

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # 1. Power flux profiles (MW)
    ax = axes[0, 0]
    power_plotted = False
    if 'pow_e' in profiles:
        ax.plot(rho, profiles['pow_e'], 'b-', label='pow_e (electron)', linewidth=2)
        power_plotted = True
    if 'pow_i' in profiles:
        ax.plot(rho, profiles['pow_i'], 'r-', label='pow_i (ion)', linewidth=2)
        power_plotted = True
    if 'pow_e_fus' in profiles:
        ax.plot(rho, profiles['pow_e_fus'], 'b--', label='pow_e_fus', linewidth=2)
        power_plotted = True
    if 'pow_i_fus' in profiles:
        ax.plot(rho, profiles['pow_i_fus'], 'r--', label='pow_i_fus', linewidth=2)
        power_plotted = True
    if power_plotted:
        ax.set_xlabel(x_label)
        ax.set_ylabel('Power [MW]')
        ax.set_title('Power Flux Profiles')
        ax.legend()
        ax.grid(True)
        ax.set_xlim([0, 1])
    else:
        ax.text(0.5, 0.5, 'No power flux data', ha='center', va='center', transform=ax.transAxes)

    # 2. Auxiliary and exchange power
    ax = axes[0, 1]
    aux_plotted = False
    if 'pow_e_aux' in profiles:
        ax.plot(rho, profiles['pow_e_aux'], 'b-', label='pow_e_aux', linewidth=2)
        aux_plotted = True
    if 'pow_i_aux' in profiles:
        ax.plot(rho, profiles['pow_i_aux'], 'r-', label='pow_i_aux', linewidth=2)
        aux_plotted = True
    if 'pow_ei' in profiles:
        ax.plot(rho, profiles['pow_ei'], 'g-', label='pow_ei (e-i exchange)', linewidth=2)
        aux_plotted = True
    if 'pow_e_rad' in profiles:
        ax.plot(rho, profiles['pow_e_rad'], 'm-', label='pow_e_rad', linewidth=2)
        aux_plotted = True
    if aux_plotted:
        ax.set_xlabel(x_label)
        ax.set_ylabel('Power [MW]')
        ax.set_title('Auxiliary & Radiation Power')
        ax.legend()
        ax.grid(True)
        ax.set_xlim([0, 1])
    else:
        ax.text(0.5, 0.5, 'No auxiliary power data', ha='center', va='center', transform=ax.transAxes)

    # 3. Heating rate profiles (W/m³)
    ax = axes[1, 0]
    heat_plotted = False
    heating_styles = {
        'qfuse': ('b-', 'qfuse (fusion→e)'),
        'qfusi': ('r-', 'qfusi (fusion→i)'),
        'qrfe': ('c-', 'qrfe (RF→e)'),
        'qrfi': ('m-', 'qrfi (RF→i)'),
        'qbeame': ('g-', 'qbeame (beam→e)'),
        'qbeami': ('y-', 'qbeami (beam→i)'),
    }
    for var, (style, label) in heating_styles.items():
        if var in profiles and np.max(np.abs(profiles[var])) > 0:
            # Convert to MW/m³ for easier reading
            ax.plot(rho, profiles[var] / 1e6, style, label=label, linewidth=2)
            heat_plotted = True
    if heat_plotted:
        ax.set_xlabel(x_label)
        ax.set_ylabel('Heating rate [MW/m³]')
        ax.set_title('Heating Rate Profiles')
        ax.legend()
        ax.grid(True)
        ax.set_xlim([0, 1])
    else:
        ax.text(0.5, 0.5, 'No heating rate data', ha='center', va='center', transform=ax.transAxes)

    # 4. Radiation and energy exchange
    ax = axes[1, 1]
    rad_plotted = False
    if 'qrad' in profiles:
        ax.plot(rho, profiles['qrad'] / 1e6, 'r-', label='qrad', linewidth=2)
        rad_plotted = True
    if 'qdelt' in profiles and np.max(np.abs(profiles['qdelt'])) > 0:
        ax.plot(rho, profiles['qdelt'] / 1e6, 'b-', label='qdelt', linewidth=2)
        rad_plotted = True
    if rad_plotted:
        ax.set_xlabel(x_label)
        ax.set_ylabel('Power density [MW/m³]')
        ax.set_title('Radiation & Other Losses')
        ax.legend()
        ax.grid(True)
        ax.set_xlim([0, 1])
    else:
        ax.text(0.5, 0.5, 'No radiation data', ha='center', va='center', transform=ax.transAxes)

    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    print(f"\nPower profiles plot saved to {output_file}")


def print_summary(profiles):
    """Print summary of profiles."""
    print("\n" + "="*60)
    print("OMFIT Profiles Summary")
    print("="*60)

    if 'R0' in profiles:
        print(f"R0 = {profiles['R0']:.3f} m")
    if 'B0' in profiles:
        print(f"B0 = {profiles['B0']:.3f} T")

    print("\n--- Density Profiles ---")
    if 'ne' in profiles:
        ne = profiles['ne']
        print(f"ne:  center={ne[0]/1e20:.3f}, edge={ne[-1]/1e20:.3f} [10^20 m^-3]")
    if 'n_d' in profiles:
        nd = profiles['n_d']
        print(f"nD:  center={nd[0]/1e20:.3f}, edge={nd[-1]/1e20:.3f} [10^20 m^-3]")
    if 'n_t' in profiles:
        nt = profiles['n_t']
        print(f"nT:  center={nt[0]/1e20:.3f}, edge={nt[-1]/1e20:.3f} [10^20 m^-3]")

    # Alpha particles (thermal and fast)
    if 'n_he4' in profiles:
        na = profiles['n_he4']
        print(f"nA(th): center={na[0]/1e18:.3f}, edge={na[-1]/1e18:.3f} [10^18 m^-3]")
    if 'n_he4_2' in profiles:
        na = profiles['n_he4_2']
        print(f"nA(fast): center={na[0]/1e18:.3f}, edge={na[-1]/1e18:.3f} [10^18 m^-3]")

    if 'n_ar' in profiles:
        nar = profiles['n_ar']
        print(f"nAr: center={nar[0]/1e17:.3f}, edge={nar[-1]/1e17:.3f} [10^17 m^-3]")
        # Calculate Ar concentration
        if 'ne' in profiles:
            ne = profiles['ne']
            nar_frac = nar / ne * 100
            print(f"     Ar concentration: center={nar_frac[0]:.3f}%, edge={nar_frac[-1]:.3f}%")

    print("\n--- Current Density ---")
    if 'current_density' in profiles:
        j = profiles['current_density']
        print(f"J_total: center={j[0]/1e6:.3f}, peak={np.max(j)/1e6:.3f} [MA/m^2]")
    if 'bootstrap_current_density_onetwo' in profiles:
        jbs = profiles['bootstrap_current_density_onetwo']
        print(f"J_BS:    center={jbs[0]/1e3:.3f}, peak={np.max(jbs)/1e3:.3f} [kA/m^2]")

    # Power summary
    print("\n--- Power Sources ---")
    power_vars = {
        'pow_e': 'Total e',
        'pow_i': 'Total i',
        'pow_e_fus': 'Fusion→e',
        'pow_i_fus': 'Fusion→i',
        'pow_e_aux': 'Aux→e',
        'pow_i_aux': 'Aux→i',
        'pow_e_rad': 'Radiation',
    }
    for var, label in power_vars.items():
        if var in profiles:
            p = profiles[var]
            # Find peak value (maximum absolute value) for flux profiles
            peak = np.max(np.abs(p))
            print(f"{label:12s}: peak={peak:.2f} MW")

    # Heating rates
    heating_vars = {
        'qfuse': 'Fusion→e',
        'qfusi': 'Fusion→i',
        'qrfe': 'RF→e',
        'qrfi': 'RF→i',
        'qbeame': 'Beam→e',
        'qbeami': 'Beam→i',
        'qrad': 'Radiation',
    }
    has_heating = any(var in profiles for var in heating_vars)
    if has_heating:
        print("\n--- Heating Rates ---")
        for var, label in heating_vars.items():
            if var in profiles:
                q = profiles[var]
                peak = np.max(np.abs(q))
                if peak > 0:
                    print(f"{label:12s}: peak={peak/1e6:.3f} MW/m³")

def export_profiles_for_tr(profiles, output_file='omfit_profiles_for_tr.csv'):
    """
    Export profiles in CSV format for TR code input.

    The TR code expects density profiles with r/a as the first column.
    """
    # Use r/a if available, otherwise use rho
    if 'r_over_a' in profiles:
        x = profiles['r_over_a']
        x_name = 'r/a'
    else:
        x = profiles.get('rho', np.linspace(0, 1, 201))
        x_name = 'rho'

    # Prepare data dictionary
    data = {x_name: x}

    # Add density profiles (convert to 10^20 m^-3)
    if 'ne' in profiles:
        data['ne'] = profiles['ne'] / 1e20
    if 'n_d' in profiles:
        data['nD'] = profiles['n_d'] / 1e20
    if 'n_t' in profiles:
        data['nT'] = profiles['n_t'] / 1e20

    # Add alpha densities (handle both naming conventions)
    # First alpha (thermal): n_he4
    # Second alpha (fast): n_he4_2
    if 'n_he4' in profiles:
        data['nA_thermal'] = profiles['n_he4'] / 1e20
    if 'n_he4_2' in profiles:
        data['nA_fast'] = profiles['n_he4_2'] / 1e20

    # Add Ar impurity
    if 'n_ar' in profiles:
        data['nAr'] = profiles['n_ar'] / 1e20

    # Add power flux profiles (MW)
    power_vars = ['pow_e', 'pow_i', 'pow_e_fus', 'pow_i_fus',
                  'pow_e_aux', 'pow_i_aux', 'pow_ei', 'pow_e_rad']
    for var in power_vars:
        if var in profiles:
            data[var] = profiles[var]

    # Add heating rate profiles (convert W/m³ to MW/m³)
    heating_vars = ['qrad', 'qfuse', 'qfusi', 'qrfe', 'qrfi',
                    'qbeame', 'qbeami', 'qdelt']
    for var in heating_vars:
        if var in profiles:
            data[var] = profiles[var] / 1e6  # W/m³ -> MW/m³

    # Create DataFrame and export
    df = pd.DataFrame(data)
    df.to_csv(output_file, index=False, float_format='%.6e')
    print(f"\nProfiles exported to {output_file}")
    print(f"  Columns: {list(df.columns)}")
    print(f"  Rows: {len(df)}")

    return df


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Read OMFIT profiles namelist')
    parser.add_argument('file', nargs='?', default='profiles_CFEDR.namelist',
                        help='Input namelist file')
    parser.add_argument('-o', '--output', default='omfit_density_profiles.png',
                        help='Output plot file')
    parser.add_argument('--rho-rmin', default='~/CFEDRSW/rho_rmin.csv',
                        help='CSV file with rho-rmin mapping for coordinate transformation')
    parser.add_argument('--use-r-over-a', action='store_true',
                        help='Use r/a coordinate instead of rho (requires --rho-rmin)')
    parser.add_argument('--export-csv', default=None,
                        help='Export profiles to CSV for TR code input')
    args = parser.parse_args()

    print(f"Reading {args.file}...")
    nml = read_omfit_profiles(args.file)
    header_vars = read_header_variables(args.file)

    print("\nFound namelists:", list(nml.keys()))
    if header_vars:
        print("Header variables:", list(header_vars.keys()))
    print("\nExtracting profiles...")
    profiles = extract_profiles(nml, header_vars)

    # Coordinate transformation: rho -> r/a
    if args.use_r_over_a:
        rho_tgyro, rmin = load_rho_rmin_mapping(args.rho_rmin)
        if rho_tgyro is not None and 'rho' in profiles:
            mapper, a = build_r_over_a_mapper(rho_tgyro, rmin)
            profiles['r_over_a'] = mapper(profiles['rho'])
            profiles['a'] = a
            print(f"\nCoordinate transformation applied:")
            print(f"  Minor radius a = {a:.4f} m")
            print(f"  r/a range: {profiles['r_over_a'][0]:.4f} to {profiles['r_over_a'][-1]:.4f}")
        else:
            print("Warning: Could not apply coordinate transformation")

    print_summary(profiles)
    plot_density_profiles(profiles, args.output)

    # Plot power profiles
    power_output = args.output.replace('.png', '_power.png')
    if power_output == args.output:
        power_output = 'omfit_power_profiles.png'
    plot_power_profiles(profiles, power_output)

    # Export to CSV (default: omfit_profiles_for_tr.csv)
    csv_output = args.export_csv if args.export_csv else 'omfit_profiles_for_tr.csv'
    export_profiles_for_tr(profiles, csv_output)

if __name__ == '__main__':
    main()
