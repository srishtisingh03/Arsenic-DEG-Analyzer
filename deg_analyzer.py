import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_csv("arsenic_deg.csv")

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

plt.figure(figsize=(6,4))
plt.scatter(df["logFC"], df["neglog10p"])
plt.xlabel("logFC")
plt.ylabel("-log10(p-value)")
plt.title("Arsenic DEG Volcano Plot")

plt.savefig("volcano_plot.png")

print("\nGenerated Files:")
print("- upregulated_genes.csv")
print("- downregulated_genes.csv")
print("- summary.txt")
print("- volcano_plot.png")
