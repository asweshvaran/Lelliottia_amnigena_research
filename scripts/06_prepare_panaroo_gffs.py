import csv
import os
import glob

BASE = r"D:\Ewingella_americana_research"
RAW = os.path.join(BASE, "data", "ncbi_dataset", "data")
GFF_DIR = os.path.join(BASE, "data", "gff_panaroo")
os.makedirs(GFF_DIR, exist_ok=True)

with open(os.path.join(BASE, "data", "metadata_final.csv"), encoding="utf-8") as f:
    rows = list(csv.DictReader(f))

strict_41 = []
complex_45 = []

def prepare_gff(acc):
    gffs = glob.glob(os.path.join(RAW, acc, "*.gff"))
    fnas = glob.glob(os.path.join(RAW, acc, "*.fna"))
    if not gffs or not fnas:
        raise FileNotFoundError(f"Missing gff/fna for {acc}")
    out_gff = os.path.join(GFF_DIR, f"{acc}.gff")
    if not os.path.exists(out_gff):
        with open(gffs[0], "r", encoding="utf-8") as fg:
            gff_txt = fg.read()
        with open(out_gff, "w", encoding="utf-8", newline="\n") as fo:
            fo.write(gff_txt)
            if "##FASTA" not in gff_txt:
                if not gff_txt.endswith("\n"):
                    fo.write("\n")
                fo.write("##FASTA\n")
                with open(fnas[0], "r", encoding="utf-8") as ff:
                    fo.write(ff.read())
    return f"/mnt/d/Ewingella_americana_research/data/gff_panaroo/{acc}.gff"

for r in rows:
    acc = r["accession"]
    if r["qc_status"] == "PASS" and r["taxon_group"] != "Misclassified non-Ewingella (<86% ANI)":
        wsl_gff = prepare_gff(acc)
        complex_45.append(wsl_gff)
        if r["final_include_pangenome"] == "True":
            strict_41.append(wsl_gff)

# Add outgroup
complex_45.append(prepare_gff("GCF_000009345.1"))

with open(os.path.join(GFF_DIR, "strict_41_list.txt"), "w", newline="\n") as f:
    f.write("\n".join(strict_41) + "\n")

with open(os.path.join(GFF_DIR, "complex_45_list.txt"), "w", newline="\n") as f:
    f.write("\n".join(complex_45) + "\n")

print(f"Prepared GFFs: strict_41={len(strict_41)}, complex_45={len(complex_45)}")
