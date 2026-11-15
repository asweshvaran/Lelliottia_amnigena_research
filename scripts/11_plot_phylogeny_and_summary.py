import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from Bio import Phylo

BASE = "/mnt/d/Ewingella_americana_research"
TREE_FILE = f"{BASE}/results/phylogeny/complex_45/ewingella_complex.treefile"
META_FILE = f"{BASE}/data/metadata_final.csv"
RESF_FILE = f"{BASE}/results/abricate/abricate_resfinder.tsv"
FIG_DIR = f"{BASE}/results/figures"
os.makedirs(FIG_DIR, exist_ok=True)

meta = pd.read_csv(META_FILE).set_index("accession")
resf = pd.read_csv(RESF_FILE, sep="\t")
resf["accession"] = resf["#FILE"].apply(lambda x: os.path.basename(x).replace(".fna", ""))

# Summarize key acquired ARGs per accession
arg_summary = {}
for acc, grp in resf.groupby("accession"):
    genes = sorted(set(g.split("_")[0] for g in grp["GENE"] if not g.startswith("Oqx")))
    if genes:
        arg_summary[acc] = ", ".join(genes)

# Load and root tree on Yersinia enterocolitica outgroup
tree = Phylo.read(TREE_FILE, "newick")
outgroup_term = [t for t in tree.get_terminals() if "GCF_000009345.1" in t.name][0]
tree.root_with_outgroup(outgroup_term)
tree.ladderize()

terminals = tree.get_terminals()
n_tips = len(terminals)

niche_colors = {
    "Environment (Soil/Water)": "#2ca02c",
    "Animal": "#ff7f0e",
    "Plant/Rhizosphere": "#8c564b",
    "Human/Clinical": "#d62728",
    "Unknown/Other": "#7f7f7f",
    "Outgroup": "#000000"
}

# Compute x, y coordinates for rectangular phylogram
def get_x_positions(tree):
    depths = tree.depths()
    if not max(depths.values()):
        depths = tree.depths(unit_branch_lengths=True)
    return depths

def get_y_positions(tree):
    maxheight = tree.count_terminals()
    heights = {}
    for i, tip in enumerate(tree.get_terminals()):
        heights[tip] = maxheight - i
    def calc_row(clade):
        for subclade in clade.clades:
            if subclade not in heights:
                calc_row(subclade)
        if clade.clades:
            heights[clade] = (heights[clade.clades[0]] + heights[clade.clades[-1]]) / 2.0
    calc_row(tree.root)
    return heights

x_pos = get_x_positions(tree)
y_pos = get_y_positions(tree)
max_x = max(x_pos.values())

plt.rcParams.update({"font.size": 9, "font.family": "DejaVu Sans"})
fig, ax = plt.subplots(figsize=(15, 13), dpi=300)

def draw_clade(clade):
    x0 = x_pos[clade]
    y0 = y_pos[clade]
    if clade.clades:
        ys = [y_pos[c] for c in clade.clades]
        ax.plot([x0, x0], [min(ys), max(ys)], color="#333333", lw=1.1)
        for c in clade.clades:
            x1 = x_pos[c]
            y1 = y_pos[c]
            ax.plot([x0, x1], [y1, y1], color="#333333", lw=1.1)
            if c.confidence is not None and c.confidence >= 70 and not c.is_terminal():
                ax.text(x0 + 0.002, y1 + 0.18, f"{int(c.confidence)}", fontsize=6.5, color="#555555")
            draw_clade(c)

draw_clade(tree.root)

novel_accs = {"GCF_021498285.1", "GCF_053586105.1", "GCF_053586165.1"}

for tip in terminals:
    acc = tip.name
    xt = x_pos[tip]
    yt = y_pos[tip]
    if acc == "GCF_000009345.1":
        niche = "Outgroup"
        label = "GCF_000009345.1 | Yersinia enterocolitica (Outgroup)"
        col = "#000000"
    else:
        row = meta.loc[acc]
        niche = row["niche"]
        col = niche_colors.get(niche, "#7f7f7f")
        country = row["country"] if pd.notna(row["country"]) else "Unknown"
        tag = " [Putative novel sp.]" if acc in novel_accs else ""
        arg_txt = f" | ARGs: {arg_summary[acc]}" if acc in arg_summary else ""
        label = f"{acc} | {country}{tag}{arg_txt}"

    ax.scatter(xt, yt, color=col, s=45, zorder=5, edgecolor="black", linewidth=0.4)
    ax.plot([xt, max_x * 1.03], [yt, yt], color="#cccccc", ls=":", lw=0.6)
    fw = "bold" if (acc in novel_accs or acc in arg_summary or acc == "GCF_000009345.1") else "normal"
    ax.text(max_x * 1.04, yt, label, va="center", fontsize=7.8, color=col if acc in novel_accs else "#111111", fontweight=fw)

ax.set_xlim(-0.02, max_x * 2.15)
ax.set_ylim(0, n_tips + 1.5)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.set_yticks([])
ax.set_xlabel("Substitutions per site (Core-Genome ML Tree, GTR+F+G4, 1,000 UFBoot)", fontsize=10)
ax.set_title("Rooted Core-Genome Phylogenomics of Ewingella americana & Sister Genomospecies (n=45)", fontsize=12.5, fontweight="bold", pad=12)

handles = [mpatches.Patch(color=c, label=k) for k, c in niche_colors.items()]
ax.legend(handles=handles, title="Ecological Niche", loc="lower left", fontsize=8.5, title_fontsize=9.5, frameon=True)

plt.tight_layout()
plt.savefig(f"{FIG_DIR}/Fig5_Core_Genome_Phylogeny.png", dpi=300, bbox_inches="tight")
plt.savefig(f"{FIG_DIR}/Fig5_Core_Genome_Phylogeny.pdf", bbox_inches="tight")
plt.close()

# --- FIGURE 6: Antitumor Effector Conservation (PMC12710904 link) ---
t5 = pd.read_csv(f"{BASE}/results/tables/TableS5_Antitumor_Candidate_Genes_Pangenome.csv")
selected = [
    ("group_4132", "Putative hemolysin 1 homologue (89 aa)"),
    ("group_6209", "Putative hemolysin 2 homologue (81 aa)"),
    ("group_356", "ShlB/FhaC/HecB hemolysin activator (569 aa)"),
    ("group_6344", "Thermostable hemolysin homologue (226 aa)"),
    ("hha", "Hemolysin modulator Hha (73 aa)"),
    ("ecnB", "Entericidin B toxin homologue (49 aa)"),
    ("cvpA", "Colicin V production protein CvpA (163 aa)"),
    ("group_5994", "T6SS Hcp effector homologue (160 aa)"),
    ("group_6324", "T6SS contractile sheath (172 aa)"),
    ("fnr", "Anaerobic regulator Fnr (251 aa)"),
    ("cheA", "Chemotaxis CheA (678 aa)"),
    ("cheY", "Chemotaxis CheY (130 aa)"),
    ("group_6072", "Flagellin N-term domain homologue (284 aa)"),
]
rows_plot = []
for g, desc in selected:
    m = t5[t5["Gene"] == g]
    if len(m) > 0:
        rows_plot.append({"Label": desc, "Pct": float(m["pct"].values[0]), "Count": int(m["n_genomes"].values[0])})

pdf = pd.DataFrame(rows_plot)
fig, ax = plt.subplots(figsize=(9.5, 5.5), dpi=300)
cols = ["#d62728" if p == 100 else ("#ff7f0e" if p >= 95 else "#1f77b4") for p in pdf["Pct"]]
bars = ax.barh(pdf["Label"][::-1], pdf["Pct"][::-1], color=cols[::-1], edgecolor="black", linewidth=0.6)
ax.set_xlim(0, 115)
ax.set_xlabel("Prevalence Across Strict E. americana Genomes (%, n=41)", fontsize=10)
ax.set_title("Species-Wide Conservation of Candidate Cytotoxic, Hypoxic & PAMP Homologues (n=41)", fontsize=11, fontweight="bold")
for b, cnt, pct in zip(bars, pdf["Count"][::-1], pdf["Pct"][::-1]):
    ax.text(pct + 1.2, b.get_y() + b.get_height()/2, f"{cnt}/41 ({pct:.1f}%)", va="center", fontsize=8.5, fontweight="bold")
plt.tight_layout()
plt.savefig(f"{FIG_DIR}/Fig6_Antitumor_Toolkit_Conservation.png", dpi=300, bbox_inches="tight")
plt.savefig(f"{FIG_DIR}/Fig6_Antitumor_Toolkit_Conservation.pdf", bbox_inches="tight")
plt.close()
print("Saved Fig5 and Fig6 successfully!")
