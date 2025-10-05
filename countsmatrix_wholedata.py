#!/usr/bin/env python
# coding: utf-8

import os
import glob
import pandas as pd
import time

# Path where your featureCounts outputs are stored
path = "/home/ddibina_93/RNA_seq/quants"
files = glob.glob(os.path.join(path, "*.txt"))

print("Files found:", files)

all_counts = []

for file in files:
    start_time = time.time()
    # Read the featureCounts file (skip headers starting with "#")
    df = pd.read_csv(file, sep="\t", comment="#")

    # Extract sample name from filename
    sample_name = os.path.basename(file).replace("_featurecounts.txt", "")
    df = df[["Geneid", df.columns[-1]]]  # Keep Geneid + counts column
    df.rename(columns={df.columns[-1]: sample_name}, inplace=True)

    all_counts.append(df)

    elapsed = (time.time() - start_time) / 60
    print(f"✅ Completed {sample_name} | Rows: {df.shape[0]} | Time: {elapsed:.2f} min")

# Merge all samples by Geneid
counts_matrix = all_counts[0]
for df in all_counts[1:]:
    counts_matrix = counts_matrix.merge(df, on="Geneid", how="outer")

# Save merged count matrix
output_file = os.path.join(path, "GSE_counts_matrix.csv")
counts_matrix.to_csv(output_file, index=False)

print("\n🎉 All files processed!")
print("Merged matrix shape:", counts_matrix.shape)
print("Saved to:", output_file)
