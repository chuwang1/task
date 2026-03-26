import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# Number of lines to plot for time-space evolution profiles (files 045-064)
# Must be >= 2 to include first and last time points
NUM_EVOLUTION_LINES = 8

def select_evolution_columns(y_cols, num_lines):
    """
    Select columns for time-space evolution plots.
    Always includes first and last, with evenly spaced points in between.

    Args:
        y_cols: list of column names (e.g., ['NE_1', 'NE_2', ..., 'NE_51'])
        num_lines: number of lines to plot (must be >= 2)

    Returns:
        list of selected column names
    """
    n_total = len(y_cols)

    if n_total <= num_lines:
        return list(y_cols)

    if num_lines < 2:
        num_lines = 2

    # Generate evenly spaced indices including first (0) and last (n_total-1)
    indices = np.linspace(0, n_total - 1, num_lines, dtype=int)

    # Ensure unique indices
    indices = sorted(set(indices))

    return [y_cols[i] for i in indices]

def is_evolution_file(csv_file):
    """
    Check if file is a time-space evolution profile (files 045-064).
    These files have many columns like NE_1, NE_2, ..., NE_51
    """
    basename = os.path.basename(csv_file)
    # Check if it's tr_data_045 to tr_data_064
    if 'tr_data_' in basename:
        try:
            num = int(basename.replace('tr_data_', '').replace('.csv', ''))
            return 45 <= num <= 64
        except ValueError:
            pass
    return False

def parse_and_plot(csv_file, num_evolution_lines=None):
    """
    Parse and plot CSV file.

    Args:
        csv_file: path to the CSV file
        num_evolution_lines: number of lines for evolution plots (default: NUM_EVOLUTION_LINES)
    """
    if num_evolution_lines is None:
        num_evolution_lines = NUM_EVOLUTION_LINES

    # Read the first line to get the title (metadata)
    with open(csv_file, 'r') as f:
        title_line = f.readline().strip()

    # Read the data, skipping the first line (metadata), header is on the second line
    try:
        # index_col=False prevents pandas from using the first column as index when
        # there is a trailing comma (more data cols than headers)
        df = pd.read_csv(csv_file, skiprows=1, index_col=False)

        # Drop columns that are unnamed (artifacts of trailing commas)
        df = df.loc[:, ~df.columns.str.contains('^Unnamed')]

        # Skip the first row for tr_data_004.csv to ignore initial transient
        if 'tr_data_004' in os.path.basename(csv_file):
            print(f"Skipping first row for {csv_file} due to initial transient.")
            df = df.iloc[1:]
    except Exception as e:
        print(f"Error reading CSV: {e}")
        return

    # Basic Plotting
    if len(df.columns) > 1:
        x_col = df.columns[0]
        y_cols = list(df.columns[1:])

        # For files with many columns, select subset of columns
        if len(y_cols) > num_evolution_lines:
            selected_cols = select_evolution_columns(y_cols, num_evolution_lines)
            print(f"Selecting {len(selected_cols)} of {len(y_cols)} columns to plot")
        else:
            selected_cols = y_cols

        plt.figure(figsize=(10, 6))
        for y_col in selected_cols:
            plt.plot(df[x_col], df[y_col], label=y_col)

        plt.xlabel(x_col)
        plt.title(title_line)
        plt.legend()
        plt.grid(True)

        output_img = csv_file.replace('.csv', '.png')
        plt.savefig(output_img)
        # print(f"Plot saved to {output_img}")
        # plt.show() # Uncomment if running locally with display
    else:
        print("Not enough columns to plot.")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Parse and plot TR CSV files')
    parser.add_argument('files', nargs='*', default=['tr_data_001.csv'],
                        help='CSV files to process')
    parser.add_argument('-n', '--num-lines', type=int, default=NUM_EVOLUTION_LINES,
                        help=f'Number of lines for evolution plots (default: {NUM_EVOLUTION_LINES})')

    args = parser.parse_args()

    for f in args.files:
        if os.path.exists(f):
            parse_and_plot(f, num_evolution_lines=args.num_lines)
        else:
            print(f"File not found: {f}")
