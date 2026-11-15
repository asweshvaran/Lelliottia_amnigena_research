import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

BASE = "/mnt/d/Ewingella_americana_research"
ABR_DIR = f"{BASE}/results/abricate"
FIG_DIR = f"{BASE}/results/figures"
TAB_DIR = f"{BASE}/results/tables"

meta = pd.read_csv(f"{BASE}/data/metadata_final.csv")
valid_accs = set(meta[meta["taxon_group"] != "Misclassified non-Ewingella (<86% ANI)"]["accession"])
strict_accs = set(meta[meta["final_include_pangenome"]]["accession"])

for db in ["card", "ncbi", "resfinder", "vfdb", "plasmidfinder"]:
    df = pd.read_csv(f"{ABR_DIR}/abricate_{db}.tsv", sep="\t")
    df["accession"] = df["#FILE"].apply(lambda x: os.path.basename(x).replace(".fna", ""))
    df_strict = df[df["accession"].isin(strict_accs)]
    print(f"=== {db.upper()} | All 57: {len(df)} hits ({df['GENE'].nunique()} genes) | Strict 41 E. americana: {len(df_strict)} hits ({df_strict['GENE'].nunique()} genes) ===")
    if len(df_strict) > 0:
        print(df_strict.groupby("GENE")["accession"].nunique().sort_values(ascending=False).head(15).to_string())
        print()

# Build combined Resistome (CARD + NCBI) + Virulome (VFDB) + Plasmidome (PlasmidFinder) summary table for strict 41 + 3 putative novel (44 Ewingella)
ew_meta = meta[meta["accession"].isin(valid_accs)].copy()
ew_meta = ew_meta.sort_values(["taxon_group", "niche", "accession"])
ordered_accs = ew_meta["accession"].tolist()

def get_matrix(db):
    df = pd.read_csv(f"{ABR_DIR}/abricate_{db}.tsv", sep="\t")
    df["accession"] = df["#FILE"].apply(lambda x: os.path.basename(x).replace(".fna", ""))
    df = df[df["accession"].isin(valid_accs)]
    if len(df) == 0:
        return pd.DataFrame(index=ordered_accs)
    piv = df.pivot_table(index="accession", columns="GENE", values="%IDENTITY", aggfunc="max")
    piv = piv.reindex(ordered_accs).fillna(0.0)
    piv = (piv > 0).astype(int)
    return piv

mat_card = get_matrix("card")
mat_vfdb = get_matrix("vfdb")
mat_plas = get_matrix("plasmidfinder")

mat_card.to_csv(f"{TAB_DIR}/TableS2_CARD_Resistome_Matrix.csv")
mat_vfdb.to_csv(f"{TAB_DIR}/TableS3_VFDB_Virulome_Matrix.csv")
mat_plas.to_csv(f"{TAB_DIR}/TableS4_PlasmidFinder_Matrix.csv")

# Plot Figure 4: Resistome + Virulome + Plasmidome Heatmap across 44 Ewingella strains
fig, axes = plt.subplots(1, 3, figsize=(18, 10), gridspec_kw={"width_ratios": [1.5, 1.8, 0.7]}, sharey=True)

niche_labels = [f"{r['accession']} ({r['niche'][:6]})" for _, r in ew_meta.iterrows()]

sns.heatmap(mat_card, cmap="Reds", cbar=False, linewidths=0.3, linecolor="lightgray", ax=axes[0])
axes[0].set_title("A. Resistome (CARD)", fontweight="bold", fontsize=11)
axes[0].set_yticklabels(niche_labels, fontsize=7.5)
axes[0].set_xticklabels(axes[0].get_xticklabels(), rotation=45, ha="right", fontsize=7.5)
axes[0].set_ylabel("Ewingella Genomes (n=44)", fontsize=10)

if mat_vfdb.shape[1] > 0:
    sns.heatmap(mat_vfdb, cmap="Blues", cbar=False, linewidths=0.3, linecolor="lightgray", ax=axes[1])
axes[1].set_title("B. Virulome (VFDB)", fontweight="bold", fontsize=11)
axes[1].set_xticklabels(axes[1].get_xticklabels(), rotation=45, ha="right", fontsize=7.5)
axes[1].set_ylabel("")

if mat_plas.shape[1] > 0:
    sns.heatmap(mat_plas, cmap="Greens", cbar=False, linewidths=0.3, linecolor="lightgray", ax=axes[2])
axes[2].set_title("C. Plasmid Replicons", fontweight="bold", fontsize=11)
axes[2].set_xticklabels(axes[2].get_xticklabels(), rotation=45, ha="right", fontsize=7.5)
axes[2].set_ylabel("")

plt.suptitle("Resistome, Virulome, and Plasmidome Profiles Across Ewingella Genomes (n=44)", fontsize=13, fontweight="bold", y=0.98)
plt.tight_layout()
plt.savefig(f"{FIG_DIR}/Fig4_Resistome_Virulome_Plasmidome.png", dpi=300, bbox_inches="tight")
plt.savefig(f"{FIG_DIR}/Fig4_Resistome_Virulome_Plasmidome.pdf", bbox_inches="tight")
plt.close()
print("Saved Fig4 and Tables S2-S4!")
