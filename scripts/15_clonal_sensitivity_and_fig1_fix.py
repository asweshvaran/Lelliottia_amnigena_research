import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import curve_fit
from scipy.spatial.distance import pdist, squareform

BASE = "/mnt/d/Ewingella_americana_research"
meta = pd.read_csv(f"{BASE}/data/metadata_final.csv")
ani = pd.read_csv(f"{BASE}/results/tables/fastani_matrix_58x58.csv", index_col=0)
rtab = pd.read_csv(f"{BASE}/results/panaroo_strict_41/gene_presence_absence.Rtab", sep="\t", index_col=0)

FIG_DIR = f"{BASE}/results/figures"
TAB_DIR = f"{BASE}/results/tables"

# --- 1. Dereplicate strict 41 at 99.95% ANI (keep highest N50 / fewest contigs per clonal cluster) ---
strict41 = meta[meta["final_include_pangenome"] == True].copy()
strict41 = strict41.sort_values(["contigs", "contig_n50"], ascending=[True, False])
accs = strict41["accession"].tolist()

clusters = []
assigned = set()
for a in accs:
    if a in assigned:
        continue
    grp = [a]
    assigned.add(a)
    for b in accs:
        if b not in assigned and ani.loc[a, b] >= 99.95:
            grp.append(b)
            assigned.add(b)
    clusters.append(grp)

derep_accs = [c[0] for c in clusters]
print(f"Strict 41 genomes -> {len(derep_accs)} non-clonal clusters (at >=99.95% ANI threshold):")
for c in clusters:
    if len(c) > 1:
        print("  Clonal cluster:", c)

# Sensitivity check: Heaps' law & PERMANOVA on dereplicated set
df_d = rtab[derep_accs]
df_d = df_d[df_d.sum(axis=1) > 0]
mat_d = (df_d.values > 0).astype(int)
n_genes_d, n_gen_d = mat_d.shape
rng = np.random.default_rng(42)
pan_c = np.zeros((100, n_gen_d))
for p in range(100):
    sub = mat_d[:, rng.permutation(n_gen_d)]
    pan_c[p, :] = (np.cumsum(sub, axis=1) > 0).sum(axis=0)
x_d = np.arange(1, n_gen_d + 1)
popt_d, pcov_d = curve_fit(lambda N, k, g: k * (N ** g), x_d, pan_c.mean(axis=0), p0=[3500, 0.25], maxfev=10000)
perr_d = np.sqrt(np.diag(pcov_d))
print(f"Dereplicated (n={n_gen_d}): Total genes={n_genes_d}, Heaps gamma={popt_d[1]:.4f} +/- {perr_d[1]:.4f}")

# PERMANOVA on dereplicated known-niche genomes
meta_idx = meta.set_index("accession")
freqs_d = mat_d.sum(axis=1) / n_gen_d
acc_mat_d = mat_d[(mat_d.sum(axis=1) > 1) & (freqs_d < 0.95), :].T
D_d = squareform(pdist(acc_mat_d, metric="jaccard"))
niches_d = np.array([meta_idx.loc[g, "niche"] for g in derep_accs])
known_m = niches_d != "Unknown/Other"
Dk = D_d[np.ix_(known_m, known_m)]
gk = niches_d[known_m]
print(f"Dereplicated known-niche counts (n={known_m.sum()}):", pd.Series(gk).value_counts().to_dict())

def permanova_pseudo_f(D_mat, grp):
    n_k = len(grp)
    uniq = np.unique(grp)
    ss_total = np.sum(D_mat ** 2) / (2 * n_k)
    ss_within = sum(np.sum(D_mat[np.ix_(grp == u, grp == u)] ** 2) / (2 * (grp == u).sum()) for u in uniq if (grp == u).sum() > 1)
    ss_between = ss_total - ss_within
    f_stat = (ss_between / (len(uniq) - 1)) / (ss_within / (n_k - len(uniq)))
    return f_stat, ss_between / ss_total

obs_f, obs_r2 = permanova_pseudo_f(Dk, gk)
perm_f = [permanova_pseudo_f(Dk, rng.permutation(gk))[0] for _ in range(999)]
p_val_d = (np.sum(np.array(perm_f) >= obs_f) + 1) / 1000.0
print(f"Dereplicated PERMANOVA (n={known_m.sum()}): R2={obs_r2:.3f}, pseudo-F={obs_f:.2f}, p={p_val_d:.4f}")

# --- 2. Regenerate Fig 1 cleanly with exact consistent counts (41 strict + 3 novel = 44 Ewingella QC-passed) ---
plt.rcParams.update({"font.size": 10, "font.family": "DejaVu Sans"})
fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.2), dpi=300)

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
axes[0].set_title("A. Taxonomic Verification (n=57)", fontweight="bold")
axes[0].set_xlabel("ANI to Reference GCF_019047465.1 (%)")
axes[0].set_ylabel("GC Content (%)")
axes[0].legend(fontsize=7.2, loc="lower left", frameon=True)

qc_pal = {"PASS": "#2ca02c", "FAIL": "#d62728"}
sns.scatterplot(
    data=meta, x="length_bp", y="contigs", hue="qc_status",
    palette=qc_pal, s=70, alpha=0.85, ax=axes[1]
)
axes[1].set_yscale("log")
axes[1].axhline(300, color="gray", linestyle="--", linewidth=1)
axes[1].set_title("B. Assembly QC (47 PASS, 10 FAIL)", fontweight="bold")
axes[1].set_xlabel("Genome Size (bp)")
axes[1].set_ylabel("Contig Count (log scale)")
axes[1].legend(title="QC Status", fontsize=8.5, title_fontsize=9)

# Panel C: Strictly the 44 QC-passed Ewingella genomes (41 strict E. americana + 3 putative novel)
ew_qc44 = meta[(meta["qc_status"] == "PASS") & (meta["taxon_group"] != "Misclassified non-Ewingella (<86% ANI)")].copy()
ct_tab = pd.crosstab(ew_qc44["niche"], ew_qc44["taxon_group"])
order_niches = ew_qc44["niche"].value_counts().index.tolist()
ct_tab = ct_tab.reindex(order_niches).fillna(0)
cols_order = [c for c in ["E. americana sensu stricto (>=95% ANI)", "Putative novel Ewingella sp. (93.3-93.7% ANI)"] if c in ct_tab.columns]
ct_tab = ct_tab[cols_order]
ct_tab.plot(kind="barh", stacked=True, color=["#1f77b4", "#ff7f0e"], edgecolor="black", linewidth=0.5, ax=axes[2])
axes[2].invert_yaxis()
axes[2].set_title("C. Ecological Niches (QC-Passed Ewingella, n=44)", fontweight="bold")
axes[2].set_xlabel("Number of Genomes")
axes[2].set_ylabel("")
axes[2].legend(["Strict E. americana (n=41)", "Putative novel sp. (n=3)"], fontsize=7.8, loc="lower right")
for i, n_name in enumerate(order_niches):
    tot = int(ct_tab.loc[n_name].sum())
    axes[2].text(tot + 0.25, i, str(tot), va="center", fontsize=9, fontweight="bold")
axes[2].set_xlim(0, ct_tab.sum(axis=1).max() + 2.2)

plt.tight_layout()
plt.savefig(f"{FIG_DIR}/Fig1_Dataset_QC_and_Ecology.png", dpi=300, bbox_inches="tight")
plt.savefig(f"{FIG_DIR}/Fig1_Dataset_QC_and_Ecology.pdf", bbox_inches="tight")
plt.close()
print("Regenerated Fig1 cleanly with exact n=41 + n=3 = 44 breakdown!")
