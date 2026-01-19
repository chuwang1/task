import pandas as pd
import matplotlib.pyplot as plt
import sys
import os

def parse_and_plot(csv_file):
    print(f"Processing {csv_file}...")
    
    # Read the first line to get the title (metadata)
    with open(csv_file, 'r') as f:
        title_line = f.readline().strip()
        
    print(f"Metadata: {title_line}")
    
    # Read the data, skipping the first line (metadata), header is on the second line
    try:
        # index_col=False prevents pandas from using the first column as index when 
        # there is a trailing comma (more data cols than headers)
        df = pd.read_csv(csv_file, skiprows=1, index_col=False)
        
        # Drop columns that are unnamed (artifacts of trailing commas)
        df = df.loc[:, ~df.columns.str.contains('^Unnamed')]
    except Exception as e:
        print(f"Error reading CSV: {e}")
        return

    print("Columns found:", df.columns.tolist())
    print(df.head())
    
    # Basic Plotting
    if len(df.columns) > 1:
        x_col = df.columns[0]
        y_cols = df.columns[1:]
        
        plt.figure(figsize=(10, 6))
        for y_col in y_cols:
            plt.plot(df[x_col], df[y_col], label=y_col)
        
        plt.xlabel(x_col)
        plt.title(title_line)
        plt.legend()
        plt.grid(True)
        
        output_img = csv_file.replace('.csv', '.png')
        plt.savefig(output_img)
        print(f"Plot saved to {output_img}")
        # plt.show() # Uncomment if running locally with display
    else:
        print("Not enough columns to plot.")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        files = sys.argv[1:]
    else:
        # Default to processing a specific file for demonstration
        files = ['tr_data_001.csv'] # Default
        
    for f in files:
        if os.path.exists(f):
            parse_and_plot(f)
        else:
            print(f"File not found: {f}")
