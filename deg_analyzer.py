import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

input_file = sys.argv[1]

df = pd.read_csv(input_file)

up = df[(df["logFC"] > 1) & (df["pvalue"] < 0.05)]
down = df[(df["logFC"] < -1) & (df["pvalue"] < 0.05)]

print("=" * 50)
print("ARSENIC DEG ANALYSIS REPORT")
print("=" * 50)

print(f"Total Genes        : {len(df)}")
print(f"Upregulated Genes  : {len(up)}")
print(f"Downregulated Genes: {len(down)}")

up.to_csv("upregulated_genes.csv", index=False)
down.to_csv("downregulated_genes.csv", index=False)

with open("summary.txt", "w") as f:
    f.write("ARSENIC DEG ANALYSIS REPORT\n")
    f.write(f"Total Genes: {len(df)}\n")
    f.write(f"Upregulated Genes: {len(up)}\n")
    f.write(f"Downregulated Genes: {len(down)}\n")

df["neglog10p"] = -np.log10(df["pvalue"])

nonsig = df[(df["pvalue"] >= 0.05) |
            ((df["logFC"] <= 1) & (df["logFC"] >= -1))]

up = df[(df["logFC"] > 1) & (df["pvalue"] < 0.05)]

down = df[(df["logFC"] < -1) & (df["pvalue"] < 0.05)]

plt.figure(figsize=(8,6))

plt.scatter(nonsig["logFC"], nonsig["neglog10p"],
            color="gray", label="Not Significant")

plt.scatter(up["logFC"], up["neglog10p"],
            color="red", label="Upregulated")

plt.scatter(down["logFC"], down["neglog10p"],
            color="blue", label="Downregulated")

plt.xlabel("logFC")
plt.ylabel("-log10(p-value)")
plt.title("Arsenic DEG Volcano Plot")

plt.legend()

plt.savefig("volcano_plot.png")
arsenic_genes = [
    "OsNIP2;1",
    "OsABCC1",
    "OsHAC1",
    "OsHAC4",
    "OsPCS1",
    "OsLsi2"
]

print("\nArsenic-related Genes Found")

for gene in arsenic_genes:
    if gene in df["Gene"].values:
        print("-", gene)

print("\nGenerated Files:")
print("- upregulated_genes.csv")
print("- downregulated_genes.csv")
print("- summary.txt")
print("- volcano_plot.png")
