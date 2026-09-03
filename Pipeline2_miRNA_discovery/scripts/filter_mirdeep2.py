import sys
import pandas as pd

input_file = sys.argv[1]
output_file = sys.argv[2]

# Find correct header line
with open(input_file) as f:
    lines = f.readlines()

start = None
for i, line in enumerate(lines):
    if line.startswith("provisional id"):
        start = i
        break

# Read table
df = pd.read_csv(input_file, sep="\t", skiprows=start)

# Filter
filtered = df[
    (df["miRDeep2 score"] >= 5) &
    (df["mature read count"] >= 10) &
    (df["star read count"] >= 1)
]

filtered.to_csv(output_file, index=False)
