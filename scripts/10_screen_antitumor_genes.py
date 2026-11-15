import pandas as pd

df = pd.read_csv("/mnt/d/Ewingella_americana_research/results/panaroo_strict_41/gene_presence_absence.csv", low_memory=False)
genome_cols = [c for c in df.columns if c.startswith("GCF_") or c.startswith("GCA_")]
n_tot = len(genome_cols)
df["n_genomes"] = df[genome_cols].notna().sum(axis=1)
df["pct"] = 100.0 * df["n_genomes"] / n_tot

keywords = ["hemolysin", "haemolysin", "cytolysin", "toxin", "colicin", "bacteriocin", "hly", "shlA", "shlB", "tly", "flagellin", "fliC", "cheA", "cheY", "motA", "motB", "type VI", "tss", "vgrG", "hcp", "type III", "fnr", "narG"]
mask = df["Annotation"].fillna("").str.lower().apply(lambda x: any(k.lower() in x for k in keywords)) | \
       df["Gene"].fillna("").str.lower().apply(lambda x: any(k.lower() in x for k in keywords))

sub = df[mask][["Gene", "Annotation", "n_genomes", "pct"]].sort_values(["n_genomes", "Gene"], ascending=[False, True])
print(f"Found {len(sub)} antitumor/cytotoxic/immunomodulatory candidate gene clusters across {n_tot} E. americana genomes:\n")
print(sub.to_string(index=False))
sub.to_csv("/mnt/d/Ewingella_americana_research/results/tables/TableS5_Antitumor_Candidate_Genes_Pangenome.csv", index=False)
