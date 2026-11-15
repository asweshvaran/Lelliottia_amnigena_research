import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import curve_fit
from scipy.spatial.distance import pdist, squareform
from scipy.stats import kruskal
from sklearn.decomposition import PCA

BASE = "/mnt/d/Ewingella_americana_research"
RTAB = f"{BASE}/results/panaroo_strict_41/gene_presence_absence.Rtab"
META = f"{BASE}/data/metadata_final.csv"
FIG_DIR = f"{BASE}/results/figures"
TAB_DIR = f"{BASE}/results/tables"
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(TAB_DIR, exist_ok=True)

# Load binary presence/absence matrix (genes x genomes)
df = pd.read_csv(RTAB, sep="\t", index_col=0)
# Ensure binary 0/1
mat = (df.values > 0).astype(int)
n_genes, n_genomes = mat.shape
genomes = df.columns.tolist()
print(f"Loaded Panaroo Rtab: {n_genes} gene clusters across {n_genomes} strict E. americana genomes.")

# Partition counts
freqs = mat.sum(axis=1) / n_genomes
n_strict_core = int((freqs == 1.0).sum())
n_core = int((freqs >= 0.99).sum())
n_soft_core = int(((freqs >= 0.95) & (freqs < 0.99)).sum())
n_shell = int(((freqs >= 0.15) & (freqs < 0.95)).sum())
n_cloud = int((freqs < 0.15).sum())
n_singletons = int((mat.sum(axis=1) == 1).sum())

print(f"Strict Core (100%): {n_strict_core} | Core (>=99%): {n_core} | Soft-core (95-99%): {n_soft_core} | Shell (15-95%): {n_shell} | Cloud (<15%): {n_cloud} | Singletons: {n_singletons}")

# --- 1. Heaps' Law & Core-Genome Decay (100 random permutations) ---
rng = np.random.default_rng(42)
n_perm = 100
pan_curves = np.zeros((n_perm, n_genomes))
core_curves = np.zeros((n_perm, n_genomes))

for p in range(n_perm):
    order = rng.permutation(n_genomes)
    sub = mat[:, order]
    cum_pan = np.cumsum(sub, axis=1) > 0
    cum_core = np.cumprod(sub, axis=1) > 0
    pan_curves[p, :] = cum_pan.sum(axis=0)
    core_curves[p, :] = cum_core.sum(axis=0)

x = np.arange(1, n_genomes + 1)
pan_mean = pan_curves.mean(axis=0)
pan_sd = pan_curves.std(axis=0)
core_mean = core_curves.mean(axis=0)
core_sd = core_curves.std(axis=0)

# Fit Heaps' Law: P(N) = kappa * N^gamma
def heaps(N, kappa, gamma):
    return kappa * (N ** gamma)

popt, pcov = curve_fit(heaps, x, pan_mean, p0=[3500, 0.25], maxfev=10000)
kappa, gamma = popt
perr = np.sqrt(np.diag(pcov))
print(f"Heaps' Law Fit: kappa = {kappa:.2f} +/- {perr[0]:.2f}, gamma = {gamma:.4f} +/- {perr[1]:.4f}")
print("Pan-genome status:", "OPEN (gamma > 0)" if gamma > 0 else "CLOSED")

# Fit Core decay: C(N) = omega + theta * exp(-N / tau)
def core_decay(N, omega, theta, tau):
    return omega + theta * np.exp(-N / tau)

popt_c, _ = curve_fit(core_decay, x, core_mean, p0=[3000, 1000, 10], maxfev=10000)
print(f"Core-genome asymptote (omega) = {popt_c[0]:.1f} genes")

# Save Table 2
tab2 = pd.DataFrame([
    {"Metric": "Total Genomes (Strict E. americana QC-passed)", "Value": n_genomes},
    {"Metric": "Total Pan-Genome Gene Families", "Value": n_genes},
    {"Metric": "Strict Core Genes (100% of strains)", "Value": n_strict_core},
    {"Metric": "Core Genes (99% <= strains <= 100%)", "Value": n_core},
    {"Metric": "Soft-Core Genes (95% <= strains < 99%)", "Value": n_soft_core},
    {"Metric": "Shell Genes (15% <= strains < 95%)", "Value": n_shell},
    {"Metric": "Cloud Genes (<15% of strains)", "Value": n_cloud},
    {"Metric": "Strain-Specific Singletons", "Value": n_singletons},
    {"Metric": "Heaps' Law kappa", "Value": f"{kappa:.2f} ± {perr[0]:.2f}"},
    {"Metric": "Heaps' Law gamma (openness exponent)", "Value": f"{gamma:.4f} ± {perr[1]:.4f}"},
    {"Metric": "Core Asymptote (omega)", "Value": f"{popt_c[0]:.1f}"}
])
tab2.to_csv(f"{TAB_DIR}/Table2_Pangenome_Summary.csv", index=False)

# --- 2. Accessory Genome PCoA / PCA & PERMANOVA by Niche ---
meta = pd.read_csv(META).set_index("accession")
# Accessory genes: present in >1 and <95% of genomes
acc_mask = (mat.sum(axis=1) > 1) & (freqs < 0.95)
acc_mat = mat[acc_mask, :].T  # genomes x accessory genes
print(f"Accessory gene matrix for PCoA: {acc_mat.shape}")

# Jaccard distance + PCoA (classical MDS)
dist_vec = pdist(acc_mat, metric="jaccard")
D = squareform(dist_vec)
n = D.shape[0]
H = np.eye(n) - np.ones((n, n)) / n
B = -0.5 * H @ (D ** 2) @ H
eigvals, eigvecs = np.linalg.eigh(B)
idx = np.argsort(eigvals)[::-1]
eigvals = eigvals[idx]
eigvecs = eigvecs[:, idx]
pos = eigvals > 0
coords = eigvecs[:, :2] * np.sqrt(np.maximum(eigvals[:2], 0))
var_exp = (eigvals[:2] / eigvals[pos].sum()) * 100

pcoa_df = pd.DataFrame(coords, index=genomes, columns=["PCo1", "PCo2"])
pcoa_df["niche"] = [meta.loc[g, "niche"] if g in meta.index else "Unknown" for g in genomes]
pcoa_df["country"] = [meta.loc[g, "country"] if g in meta.index else "Unknown" for g in genomes]
pcoa_df["length_bp"] = [meta.loc[g, "length_bp"] if g in meta.index else np.nan for g in genomes]

# One-way PERMANOVA (pseudo-F with 999 permutations) on known niches
known_mask = pcoa_df["niche"] != "Unknown/Other"
D_known = D[np.ix_(known_mask.values, known_mask.values)]
groups = pcoa_df.loc[known_mask, "niche"].values

def permanova_pseudo_f(D_mat, grp):
    n_k = len(grp)
    uniq = np.unique(grp)
    ss_total = np.sum(D_mat ** 2) / (2 * n_k)
    ss_within = 0.0
    for u in uniq:
        m = (grp == u)
        if m.sum() > 1:
            ss_within += np.sum(D_mat[np.ix_(m, m)] ** 2) / (2 * m.sum())
    ss_between = ss_total - ss_within
    df_b = len(uniq) - 1
    df_w = n_k - len(uniq)
    f_stat = (ss_between / df_b) / (ss_within / df_w)
    r2 = ss_between / ss_total
    return f_stat, r2

obs_f, obs_r2 = permanova_pseudo_f(D_known, groups)
perm_f = []
for _ in range(999):
    perm_grp = rng.permutation(groups)
    f_p, _ = permanova_pseudo_f(D_known, perm_grp)
    perm_f.append(f_p)
p_val = (np.sum(np.array(perm_f) >= obs_f) + 1) / 1000.0
print(f"Accessory Genome PERMANOVA (known niches, n={known_mask.sum()}): R2 = {obs_r2:.3f}, pseudo-F = {obs_f:.2f}, p = {p_val:.4f}")

# --- FIGURE 3: Pan-Genome Architecture & Ecological PCoA ---
plt.rcParams.update({"font.size": 10, "font.family": "DejaVu Sans"})
fig, axes = plt.subplots(1, 3, figsize=(17, 5.2), dpi=300)

# Panel A: Heaps' Law & Core Decay
axes[0].fill_between(x, pan_mean - pan_sd, pan_mean + pan_sd, color="#1f77b4", alpha=0.2)
axes[0].plot(x, pan_mean, color="#1f77b4", lw=2, label="Pan-genome (observed)")
axes[0].plot(x, heaps(x, *popt), color="#0b3c5d", ls="--", lw=1.5, label=f"Heaps' Law (γ = {gamma:.3f})")
axes[0].fill_between(x, core_mean - core_sd, core_mean + core_sd, color="#d62728", alpha=0.2)
axes[0].plot(x, core_mean, color="#d62728", lw=2, label="Core-genome (observed)")
axes[0].set_title("A. Open Pan-Genome & Core Decay", fontweight="bold")
axes[0].set_xlabel("Number of Genomes (N)")
axes[0].set_ylabel("Number of Gene Families")
axes[0].legend(fontsize=8.5, loc="center right")

# Panel B: Pan-Genome Partition Bar Chart
parts = ["Core\n(99-100%)", "Soft-core\n(95-99%)", "Shell\n(15-95%)", "Cloud\n(<15%)"]
counts = [n_core, n_soft_core, n_shell, n_cloud]
colors = ["#2ca02c", "#98df8a", "#ff7f0e", "#1f77b4"]
bars = axes[1].bar(parts, counts, color=colors, edgecolor="black", linewidth=0.7)
axes[1].set_title(f"B. Pan-Genome Partition (Total = {n_genes:,})", fontweight="bold")
axes[1].set_ylabel("Number of Gene Families")
for b, c in zip(bars, counts):
    pct = 100 * c / n_genes
    axes[1].text(b.get_x() + b.get_width()/2, c + 120, f"{c:,}\n({pct:.1f}%)", ha="center", va="bottom", fontsize=8.5, fontweight="bold")
axes[1].set_ylim(0, max(counts) * 1.18)

# Panel C: Accessory Genome PCoA (Jaccard)
niche_pal = {
    "Environment (Soil/Water)": "#2ca02c",
    "Animal": "#ff7f0e",
    "Plant/Rhizosphere": "#8c564b",
    "Human/Clinical": "#d62728",
    "Unknown/Other": "#7f7f7f"
}
sns.scatterplot(data=pcoa_df, x="PCo1", y="PCo2", hue="niche", palette=niche_pal, s=80, alpha=0.88, edgecolor="black", linewidth=0.4, ax=axes[2])
axes[2].set_title(f"C. Accessory Genome PCoA (PERMANOVA R²={obs_r2:.2f}, p={p_val:.3f})", fontweight="bold")
axes[2].set_xlabel(f"PCo1 ({var_exp[0]:.1f}% var)")
axes[2].set_ylabel(f"PCo2 ({var_exp[1]:.1f}% var)")
axes[2].legend(title="Ecological Niche", fontsize=7.5, title_fontsize=8.5, loc="best")

plt.tight_layout()
plt.savefig(f"{FIG_DIR}/Fig3_Pangenome_Dynamics_and_PCoA.png", dpi=300, bbox_inches="tight")
plt.savefig(f"{FIG_DIR}/Fig3_Pangenome_Dynamics_and_PCoA.pdf", bbox_inches="tight")
plt.close()
print("Saved Fig3 and Table 2 successfully!")
