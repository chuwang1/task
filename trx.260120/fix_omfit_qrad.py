#!/usr/bin/env python3
"""
Fix the OMFIT qrad CSV file by removing the spike at r/a > 1.0
and also check for any edge anomalies
"""
import pandas as pd
import numpy as np

# Read original file
df = pd.read_csv('omfit_qrad_for_tr.csv')
print("Original data:")
print(f"  Total points: {len(df)}")
print(f"  r/a range: {df['r_a'].min():.6f} - {df['r_a'].max():.6f}")
print(f"  qrad range: {df['qrad_MW_m3'].min():.6f} - {df['qrad_MW_m3'].max():.6f}")

# Check values near edge
print("\nOriginal edge values (last 15 points):")
print(df.tail(15).to_string(index=False))

# Remove points with r/a > 1.0
df_clean = df[df['r_a'] <= 1.0].copy()

# Also check for any sudden jumps in qrad near the edge
# The jump at r/a=0.999 from 0.06 to 0.08 and then spike to 0.56 is problematic
# Let's smooth the last few points or truncate at r/a=0.99

# Option 1: Remove all points with r/a > 0.99
df_truncated = df[df['r_a'] <= 0.99].copy()

# Option 2: Remove the spike but keep up to r/a=0.995
df_no_spike = df[df['r_a'] <= 0.995].copy()

print(f"\nCleaned data (r/a <= 1.0):")
print(f"  Total points: {len(df_clean)}")
print(f"  r/a range: {df_clean['r_a'].min():.6f} - {df_clean['r_a'].max():.6f}")
print(f"  qrad range: {df_clean['qrad_MW_m3'].min():.6f} - {df_clean['qrad_MW_m3'].max():.6f}")

print(f"\nTruncated data (r/a <= 0.99):")
print(f"  Total points: {len(df_truncated)}")
print(f"  r/a range: {df_truncated['r_a'].min():.6f} - {df_truncated['r_a'].max():.6f}")
print(f"  qrad range: {df_truncated['qrad_MW_m3'].min():.6f} - {df_truncated['qrad_MW_m3'].max():.6f}")

# Save cleaned versions
df_clean.to_csv('omfit_qrad_for_tr_clean.csv', index=False)
df_truncated.to_csv('omfit_qrad_for_tr_truncated.csv', index=False)

print("\nSaved cleaned files:")
print("  - omfit_qrad_for_tr_clean.csv (r/a <= 1.0)")
print("  - omfit_qrad_for_tr_truncated.csv (r/a <= 0.99)")

# Show the edge values of the cleaned file
print("\nCleaned edge values (last 10 points):")
print(df_clean.tail(10).to_string(index=False))
