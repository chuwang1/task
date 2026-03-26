#!/usr/bin/env python3
"""
Plot radial profiles of radiation losses: PRSUM, PRB, PRC, PRL, PCX, PIE
Reads data from tr_data_XXX.csv file containing these variables.
"""

import pandas as pd
import matplotlib.pyplot as plt
import glob
import sys

def find_radiation_csv():
    """Find CSV file containing PRSUM,PRB,PRC,PRL,PCX,PIE data"""
    for fname in sorted(glob.glob('tr_data_*.csv')):
        with open(fname, 'r') as f:
            first_line = f.readline()
            if 'PRSUM,PRB,PRC,PRL,PCX,PIE' in first_line and 'vs r' in first_line:
                return fname
    return None

def plot_radiation_loss(csv_file=None):
    # Find the CSV file if not specified
    if csv_file is None:
        csv_file = find_radiation_csv()
        if csv_file is None:
            print("Error: Cannot find CSV file with PRSUM,PRB,PRC,PRL,PCX,PIE data")
            print("Please run TR graphics menu (G R 2) first to generate the data")
            sys.exit(1)

    print(f"Reading: {csv_file}")

    # Read CSV, skip the title line
    df = pd.read_csv(csv_file, skiprows=1)

    # Get column names
    cols = df.columns.tolist()
    print(f"Columns: {cols}")

    # X axis is rho (first column)
    rho = df.iloc[:, 0]

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 7))

    # Plot each variable
    colors = ['black', 'blue', 'red', 'green', 'orange', 'purple']
    labels = ['PRSUM (Total)', 'PRB (Bremsstrahlung)', 'PRC (Cyclotron)',
              'PRL (Line)', 'PCX (CX)', 'PIE (Ionization)']
    linestyles = ['-', '--', '-.', ':', '-', '--']
    linewidths = [2.5, 1.5, 1.5, 1.5, 1.5, 1.5]

    for i, (col, color, label, ls, lw) in enumerate(zip(cols[1:7], colors, labels, linestyles, linewidths)):
        data = df[col]
        # Only plot if data is non-zero somewhere
        if data.abs().max() > 1e-10:
            ax.plot(rho, data, color=color, label=label, linestyle=ls, linewidth=lw)

    ax.set_xlabel(r'$\rho$ (normalized radius)', fontsize=12)
    ax.set_ylabel(r'Power Loss [MW/m$^3$]', fontsize=12)
    ax.set_title('Radial Profiles of Energy Losses', fontsize=14)
    ax.legend(loc='best', fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 1)

    # Save figure
    output_file = csv_file.replace('.csv', '_radiation_loss.png')
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Saved: {output_file}")

    plt.show()

if __name__ == '__main__':
    if len(sys.argv) > 1:
        plot_radiation_loss(sys.argv[1])
    else:
        plot_radiation_loss()
