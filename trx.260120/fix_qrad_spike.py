#!/usr/bin/env python3
"""
Fix OMFIT qrad CSV by removing spike at r/a > 1.0
and ensuring smooth extrapolation to r/a = 1.0
"""
import pandas as pd
import numpy as np

# Read original file
df = pd.read_csv('omfit_qrad_for_tr.csv')

print("Original data:")
print(f"  Total points: {len(df)}")
print(f"  r/a range: {df['r_a'].min():.6f} - {df['r_a'].max():.6f}")
print(f"  qrad range: {df['qrad_MW_m3'].min():.6f} - {df['qrad_MW_m3'].max():.6f}")

# Check last few points
print("\nLast 5 points before fix:")
print(df.tail(5).to_string(index=False))

# Remove points with r/a > 1.0
df_fixed = df[df['r_a'] <= 1.0].copy()

# Check if we need to add a point at r/a = 1.0
if df_fixed['r_a'].max() < 1.0:
    # Use the last qrad value for r/a = 1.0
    last_qrad = df_fixed['qrad_MW_m3'].iloc[-1]
    new_row = pd.DataFrame({'r_a': [1.0], 'qrad_MW_m3': [last_qrad]})
    df_fixed = pd.concat([df_fixed, new_row], ignore_index=True)

print(f"\nFixed data:")
print(f"  Total points: {len(df_fixed)}")
print(f"  r/a range: {df_fixed['r_a'].min():.6f} - {df_fixed['r_a'].max():.6f}")
print(f"  qrad range: {df_fixed['qrad_MW_m3'].min():.6f} - {df_fixed['qrad_MW_m3'].max():.6f}")

print("\nLast 5 points after fix:")
print(df_fixed.tail(5).to_string(index=False))

# Save fixed file
output_path = 'omfit_qrad_fixed.csv'
df_fixed.to_csv(output_path, index=False)
print(f"\nSaved: {output_path}")
