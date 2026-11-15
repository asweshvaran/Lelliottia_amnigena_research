import csv
import os
import glob
import shutil

BASE = r"D:\Ewingella_americana_research"
RAW_DATA = os.path.join(BASE, "data", "ncbi_dataset", "data")
FNA_DIR = os.path.join(BASE, "data", "genomes_fna")
OUTGROUP_DIR = os.path.join(BASE, "data", "outgroup_fna")
os.makedirs(FNA_DIR, exist_ok=True)
os.makedirs(OUTGROUP_DIR, exist_ok=True)

def classify_niche(src, host):
    s = (src or "").lower()
    h = (host or "").lower()
    combined = f"{s} {h}".strip()
    if any(k in combined for k in ["homo sapiens", "human", "urine", "faecal", "fecal", "oral", "blood", "clinical"]):
        return "Human/Clinical"
    if any(k in combined for k in ["pig", "equus caballus", "fish", "phocoena", "bombyx", "silkworm", "chicken", "mosquito"]):
        return "Animal"
    if any(k in combined for k in ["potato", "apple", "stem", "vaccinium", "arabidopsis", "cucumis", "rhizosphere", "root/leaf"]):
        return "Plant/Rhizosphere"
    if any(k in combined for k in ["soil", "water", "wastewater", "environment", "chemotactic"]):
        return "Environment (Soil/Water)"
    return "Unknown/Other"

def clean_country(geo):
    if not geo or geo.lower() in ["missing", "not applicable", "not provided", ""]:
        return "Unknown"
    return geo.split(":")[0].strip()

with open(os.path.join(BASE, "data", "metadata_all.csv"), encoding="utf-8") as f:
    rows = list(csv.DictReader(f))

curated = []
qc_pass_count = 0
for r in rows:
    acc = r["accession"]
    fnas = glob.glob(os.path.join(RAW_DATA, acc, "*.fna"))
    if not fnas:
        print(f"WARNING: No .fna found for {acc}")
        continue
    dst = os.path.join(FNA_DIR, f"{acc}.fna")
    if not os.path.exists(dst):
        shutil.copy2(fnas[0], dst)

    comp = float(r["checkm_completeness"]) if r["checkm_completeness"] else None
    contam = float(r["checkm_contamination"]) if r["checkm_contamination"] else 0.0
    ctg = int(r["contigs"]) if r["contigs"] else 9999
    length = int(r["length_bp"]) if r["length_bp"] else 0

    # Strict Q1 QC criteria:
    # Completeness >= 95% (or None if missing, we'll inspect separately, but require contigs <= 300 and length >= 3.8 Mb)
    qc_reasons = []
    if comp is not None and comp < 95.0:
        qc_reasons.append(f"low_comp({comp}%)")
    if comp is None:
        qc_reasons.append("no_checkm")
    if contam > 5.0:
        qc_reasons.append(f"high_contam({contam}%)")
    if ctg > 300:
        qc_reasons.append(f"fragmented({ctg}_contigs)")
    if length < 3800000 or length > 5500000:
        qc_reasons.append(f"atypical_len({length})")

    qc_status = "PASS" if len(qc_reasons) == 0 else "FAIL"
    if qc_status == "PASS":
        qc_pass_count += 1

    r["niche"] = classify_niche(r["isolation_source"], r["host"])
    r["country"] = clean_country(r["geo_loc_name"])
    r["qc_status"] = qc_status
    r["qc_notes"] = ";".join(qc_reasons) if qc_reasons else "OK"
    curated.append(r)

# Copy outgroup
out_fnas = glob.glob(os.path.join(RAW_DATA, "GCF_000009345.1", "*.fna"))
if out_fnas:
    shutil.copy2(out_fnas[0], os.path.join(OUTGROUP_DIR, "GCF_000009345.1.fna"))

out_csv = os.path.join(BASE, "data", "metadata_curated.csv")
with open(out_csv, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=list(curated[0].keys()))
    writer.writeheader()
    writer.writerows(curated)

print(f"Curated {len(curated)} E. americana genomes. Initial QC PASS: {qc_pass_count}, FAIL: {len(curated) - qc_pass_count}")
from collections import Counter
print("Niche distribution (all):", Counter(r["niche"] for r in curated))
print("Niche distribution (QC PASS):", Counter(r["niche"] for r in curated if r["qc_status"] == "PASS"))
