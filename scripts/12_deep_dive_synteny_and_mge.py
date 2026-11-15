import os
import pandas as pd

BASE = "/mnt/d/Ewingella_americana_research"
resf = pd.read_csv(f"{BASE}/results/abricate/abricate_resfinder.tsv", sep="\t")
plas = pd.read_csv(f"{BASE}/results/abricate/abricate_plasmidfinder.tsv", sep="\t")
resf["accession"] = resf["#FILE"].apply(lambda x: os.path.basename(x).replace(".fna", ""))
plas["accession"] = plas["#FILE"].apply(lambda x: os.path.basename(x).replace(".fna", ""))

target_accs = [
    "GCF_056155165.1", "GCF_056157825.1",
    "GCF_023970615.1", "GCF_019048185.1",
    "GCF_053586065.1", "GCF_036487165.1",
    "GCF_036487385.1", "GCF_988867145.1"
]

rows_out = []
print("=== DEEP-DIVE CONTIG CO-LOCALIZATION (ARGs + PLASMIDS) ===")
for acc in target_accs:
    r_sub = resf[(resf["accession"] == acc) & (~resf["GENE"].str.startswith("Oqx"))]
    p_sub = plas[plas["accession"] == acc]
    print(f"\n--- {acc} ---")
    for _, r in r_sub.iterrows():
        print(f"  [ARG] {r['GENE']:15s} | contig={r['SEQUENCE']} | pos={r['START']}-{r['END']} | id={r['%IDENTITY']}% | cov={r['%COVERAGE']}%")
        rows_out.append({"accession": acc, "type": "ARG", "gene": r["GENE"], "contig": r["SEQUENCE"], "start": r["START"], "end": r["END"], "identity": r["%IDENTITY"], "coverage": r["%COVERAGE"]})
    for _, p in p_sub.iterrows():
        print(f"  [REP] {p['GENE']:22s} | contig={p['SEQUENCE']} | pos={p['START']}-{p['END']} | id={p['%IDENTITY']}% | cov={p['%COVERAGE']}%")
        rows_out.append({"accession": acc, "type": "PlasmidReplicon", "gene": p["GENE"], "contig": p["SEQUENCE"], "start": p["START"], "end": p["END"], "identity": p["%IDENTITY"], "coverage": p["%COVERAGE"]})

    # Also check GFF features on the ARG-bearing contigs
    gff_path = f"{BASE}/data/gff_panaroo/{acc}.gff"
    arg_contigs = set(r_sub["SEQUENCE"].tolist())
    if os.path.exists(gff_path) and arg_contigs:
        mge_hits = []
        with open(gff_path) as fg:
            for line in fg:
                if line.startswith("##FASTA"):
                    break
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 9 and parts[0] in arg_contigs and parts[2] == "CDS":
                    attr = parts[8]
                    if any(k in attr.lower() for k in ["transposase", "integrase", "recombinase", "conjug", "type iv", "plasmid", "relaxase", "imir", "oxa", "mcr", "qnr"]):
                        prod = [x.split("=")[1] for x in attr.split(";") if x.startswith("product=")]
                        mge_hits.append((parts[0], parts[3], parts[4], prod[0] if prod else "CDS"))
        if mge_hits:
            print(f"  [MGE/Context on ARG contigs]: {len(mge_hits)} MGE/mobilization CDSs found (e.g., {mge_hits[:6]})")

pd.DataFrame(rows_out).to_csv(f"{BASE}/results/tables/TableS6_ARG_Plasmid_Colocalization.csv", index=False)
print("\nSaved TableS6_ARG_Plasmid_Colocalization.csv")
