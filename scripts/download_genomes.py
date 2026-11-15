import urllib.request
import json
import csv
import os
import zipfile
import shutil

DATA_DIR = r"D:\Ewingella_americana_research\data"
META_CSV = os.path.join(DATA_DIR, "metadata_all.csv")
ZIP_PATH = os.path.join(DATA_DIR, "ewingella_dataset.zip")
EXTRACT_DIR = os.path.join(DATA_DIR, "ncbi_dataset")

with open(META_CSV, encoding="utf-8") as f:
    rows = list(csv.DictReader(f))

accessions = [r["accession"] for r in rows]
# Add Yersinia enterocolitica subsp. enterocolitica 8081 (GCF_000009345.1) as Yersiniaceae outgroup
outgroup = "GCF_000009345.1"
all_accs = accessions + [outgroup]

print(f"Requesting download of {len(all_accs)} genomes ({len(accessions)} E. americana + 1 outgroup)...")

url = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/download"
payload = json.dumps({
    "accessions": all_accs,
    "include_annotation_type": ["GENOME_FASTA", "GENOME_GFF", "PROT_FASTA", "SEQUENCE_REPORT"]
}).encode("utf-8")

req = urllib.request.Request(
    url,
    data=payload,
    headers={"Content-Type": "application/json", "Accept": "application/zip"},
    method="POST"
)

with urllib.request.urlopen(req, timeout=300) as resp:
    total = 0
    with open(ZIP_PATH, "wb") as out:
        while True:
            chunk = resp.read(1024 * 1024)
            if not chunk:
                break
            out.write(chunk)
            total += len(chunk)
            if total % (20 * 1024 * 1024) < 1024 * 1024:
                print(f"Downloaded {total / (1024*1024):.1f} MB...")

print(f"Download complete: {ZIP_PATH} ({os.path.getsize(ZIP_PATH) / (1024*1024):.2f} MB)")

if os.path.exists(EXTRACT_DIR):
    shutil.rmtree(EXTRACT_DIR)

with zipfile.ZipFile(ZIP_PATH, "r") as zf:
    zf.extractall(DATA_DIR)

print("Extraction complete to", EXTRACT_DIR)
