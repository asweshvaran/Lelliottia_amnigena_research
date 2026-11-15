import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform

BASE = "/mnt/d/Ewingella_americana_research"
ANI_FILE = f"{BASE}/results/ani/fastani_all_vs_all.tsv"
META_FILE = f"{BASE}/data/metadata_curated.csv"
FIG_DIR = f"{BASE}/results/figures"
TAB_DIR = f"{BASE}/results/tables"
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(TAB_DIR, exist_ok=True)

meta = pd.read_csv(META_FILE)
ref_acc = "GCF_019047465.1"

# Parse FastANI
records = []
with open(ANI_FILE) as f:
    for line in f:
        q, r, ani, frag, total = line.strip().split("\t")
        qb = os.path.basename(q).replace(".fna", "")
        rb = os.path.basename(r).replace(".fna", "")
        records.append((qb, rb, float(ani)))

ani_df = pd.DataFrame(records, columns=["query", "ref", "ani"])
mat = ani_df.pivot(index="query", columns="ref", values="ani")

# Symmetrize matrix and fill missing distant pairs with 75.0
all_accs = sorted(mat.index.tolist())
mat = mat.reindex(index=all_accs, columns=all_accs)
mat = (mat + mat.T) / 2.0
mat = mat.fillna(75.0)
np.fill_diagonal(mat.values, 100.0)

mat.to_csv(f"{TAB_DIR}/fastani_matrix_58x58.csv")

# Assign taxonomic group based on ANI to reference GCF_019047465.1
ani_to_ref = mat[ref_acc].to_dict()
def assign_taxon_group(acc):
    v = ani_to_ref.get(acc, 0.0)
    if v >= 95.0:
        return "E. americana sensu stricto (>=95% ANI)"
    elif 92.0 <= v < 95.0:
        return "Putative novel Ewingella sp. (93.3-93.7% ANI)"
    else:
        return "Misclassified non-Ewingella (<86% ANI)"

meta["ani_to_ref"] = meta["accession"].map(ani_to_ref)
meta["taxon_group"] = meta["accession"].map(assign_taxon_group)
meta["final_include_pangenome"] = (meta["qc_status"] == "PASS") & (meta["taxon_group"] == "E. americana sensu stricto (>=95% ANI)")
meta.to_csv(f"{BASE}/data/metadata_final.csv", index=False)
meta.to_csv(f"{TAB_DIR}/TableS1_Genome_Metadata_QC_ANI.csv", index=False)

print("Taxonomic group counts:")
print(meta["taxon_group"].value_counts())
print("\nFinal included in strict E. americana pan-genome:", meta["final_include_pangenome"].sum())

# --- FIGURE 1: Dataset QC, Taxonomic Verification, & Ecological Distribution ---
plt.rcParams.update({"font.size": 10, "font.family": "DejaVu Sans"})
fig, axes = plt.subplots(1, 3, figsize=(16, 5.2), dpi=300)

# Panel A: ANI to Reference vs GC%
palette_tax = {
    "E. americana sensu stricto (>=95% ANI)": "#1f77b4",
    "Putative novel Ewingella sp. (93.3-93.7% ANI)": "#ff7f0e",
    "Misclassified non-Ewingella (<86% ANI)": "#d62728"
}
sns.scatterplot(
    data=meta, x="ani_to_ref", y="gc_percent", hue="taxon_group",
    style="qc_status", palette=palette_tax, s=75, alpha=0.85, ax=axes[0]
)
axes[0].axvline(95.0, color="black", linestyle="--", linewidth=1, label="95% Species Cutoff")
axes[0].set_title("A. Taxonomic Verification (ANI vs GC%)", fontweight="bold")
axes[0].set_xlabel("ANI to Reference GCF_019047465.1 (%)")
axes[0].set_ylabel("GC Content (%)")
axes[0].legend(fontsize=7.5, loc="lower left", frameon=True)

# Panel B: Genome Size vs Contig Count (QC)
qc_pal = {"PASS": "#2ca02c", "FAIL": "#d62728"}
sns.scatterplot(
    data=meta, x="length_bp", y="contigs", hue="qc_status",
    palette=qc_pal, s=70, alpha=0.85, ax=axes[1]
)
axes[1].set_yscale("log")
axes[1].axhline(300, color="gray", linestyle="--", linewidth=1)
axes[1].set_title("B. Assembly Quality Control", fontweight="bold")
axes[1].set_xlabel("Genome Size (bp)")
axes[1].set_ylabel("Contig Count (log scale)")
axes[1].legend(title="QC Status", fontsize=8.5, title_fontsize=9)

# Panel C: Ecological Niche Distribution (QC-passed E. americana + Putative novel)
sub = meta[meta["qc_status"] == "PASS"]
niche_counts = sub["niche"].value_counts()
sns.barplot(x=niche_counts.values, y=niche_counts.index, palette="viridis", ax=axes[2])
axes[2].set_title("C. Ecological Niches (QC-Passed, n={})".format(len(sub)), fontweight="bold")
axes[2].set_xlabel("Number of Genomes")
axes[2].set_ylabel("")
for i, v in enumerate(niche_counts.values):
    axes[2].text(v + 0.2, i, str(v), va="center", fontsize=9, fontweight="bold")

plt.tight_layout()
plt.savefig(f"{FIG_DIR}/Fig1_Dataset_QC_and_Ecology.png", dpi=300, bbox_inches="tight")
plt.savefig(f"{FIG_DIR}/Fig1_Dataset_QC_and_Ecology.pdf", bbox_inches="tight")
plt.close()

# --- FIGURE 2: All-vs-All FastANI Clustered Heatmap ---
# Map accession -> Niche color & Taxon color
niche_map = meta.set_index("accession")["niche"].to_dict()
tax_map = meta.set_index("accession")["taxon_group"].to_dict()
niche_map["GCF_000009345.1"] = "Outgroup (Y. enterocolitica)"
tax_map["GCF_000009345.1"] = "Outgroup"

niche_colors_lut = {
    "Environment (Soil/Water)": "#2ca02c",
    "Animal": "#ff7f0e",
    "Plant/Rhizosphere": "#8c564b",
    "Human/Clinical": "#d62728",
    "Unknown/Other": "#7f7f7f",
    "Outgroup (Y. enterocolitica)": "#000000"
}
row_colors = pd.Series([niche_colors_lut.get(niche_map.get(a, "Unknown/Other"), "#7f7f7f") for a in mat.index], index=mat.index, name="Ecological Niche")

dist = 100.0 - mat.values
np.fill_diagonal(dist, 0.0)
dist = np.clip(dist, 0, None)
Z = linkage(squareform(dist), method="average")

g = sns.clustermap(
    mat, row_linkage=Z, col_linkage=Z,
    row_colors=row_colors,
    cmap="YlGnBu", vmin=80, vmax=100,
    figsize=(13, 12), xticklabels=False, yticklabels=True,
    cbar_kws={"label": "FastANI (%)"}
)
g.ax_heatmap.tick_params(axis="y", labelsize=6.5)
g.fig.suptitle("All-vs-All FastANI Clustering of E. americana & Related Assemblies (n=58)", y=1.01, fontsize=13, fontweight="bold")
plt.savefig(f"{FIG_DIR}/Fig2_FastANI_Clustered_Heatmap.png", dpi=300, bbox_inches="tight")
plt.savefig(f"{FIG_DIR}/Fig2_FastANI_Clustered_Heatmap.pdf", bbox_inches="tight")
plt.close()
print("Saved Fig1 and Fig2 successfully!")
