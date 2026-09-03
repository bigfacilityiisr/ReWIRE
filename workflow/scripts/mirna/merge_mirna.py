import pandas as pd
import glob

# Get all filtered files
files = glob.glob("/media/user/DATADRIVE1/miRNA_Project/mirdeep/*_filtered.csv")

dfs = []
for f in files:
    df = pd.read_csv(f)
    df["sample"] = f.split("/")[-1].split("_")[0]  # add sample name
    dfs.append(df)

# Merge all
merged = pd.concat(dfs, ignore_index=True)

# Remove duplicates based on mature sequence
non_redundant = merged.drop_duplicates(subset=["consensus mature sequence"])

# Save outputs
merged.to_csv("/media/user/DATADRIVE1/miRNA_Project/mirdeep/all_samples_merged.csv", index=False)
non_redundant.to_csv("/media/user/DATADRIVE1/miRNA_Project/mirdeep/non_redundant_miRNAs.csv", index=False)

print("Total miRNAs:", len(merged))
print("Non-redundant miRNAs:", len(non_redundant))

summary = merged.groupby("sample").size().reset_index(name="miRNA_count")
summary.to_csv("/media/user/DATADRIVE1/miRNA_Project/mirdeep/miRNA_summary_per_sample.csv", index=False)
print(summary)
