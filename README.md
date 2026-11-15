# Pan-Genome Architecture, Cryptic Lineages, and One Health Mobilization of Epidemic Carbapenemase and Colistin Resistance Plasmids in *Lelliottia amnigena*

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Bioinformatics](https://img.shields.io/badge/Bioinformatics-Panaroo%20%7C%20FastANI%20%7C%20IQ--TREE-green.svg)](https://github.com/asweshvaran/Lelliottia_amnigena_research)

This repository contains the computational pipeline, analytical scripts, genome accessions, and processed tabular datasets supporting the manuscript:

> **"Pan-Genome Architecture, Cryptic Lineages, and One Health Mobilization of Epidemic Carbapenemase and Colistin Resistance Plasmids in *Lelliottia amnigena*"**  
> **Author:** Asweshvaran R.  
> *Department of Microbiology, Shriram Research Institute, 14–15, Sadarmangala Industrial Area, Whitefield Road, Bengaluru 560048, Karnataka, India.*  
> **Correspondence:** `asweshvaran@gmail.com`  
> *Submitted to: Journal of Applied Microbiology (Oxford University Press)*

---

## Abstract Overview

*Lelliottia amnigena* (formerly *Enterobacter amnigenus*) is an emerging opportunistic pathogen spanning clinical, agricultural, food, and aquatic compartments. This study establishes its species-wide pan-genome architecture, evaluates cryptic sister lineages, and resolves the genomic contexts of acquired resistance determinants:
* **Taxonomic Curation:** 57 publicly available assemblies indexed as *L. amnigena* were triaged against a *Yersinia enterocolitica* outgroup. FastANI identified 50 *L. amnigena* sensu stricto assemblies (>=98.2% ANI), three candidate sister-lineage genomes (93.35%-93.74% ANI), and four misidentified assemblies (<86% ANI).
* **Pan-Genome Dimensions:** Panaroo analysis across 41 QC-passed sensu stricto genomes resolved 9,630 gene families (3,369 strict core [34.98%], 299 soft-core [3.10%], 818 shell [8.49%], and 5,144 cloud [53.42%]). Heaps' law modeling confirms an open pan-genome (gamma = 0.2609), which remains open after clonal dereplication (gamma = 0.2769, n = 32).
* **Clonal Sensitivity in Ecology:** Accessory-genome variation initially appeared significantly associated with ecological niche (R2 = 0.138, p = 0.019, n = 31), but clonal dereplication rendered this association non-significant (R2 = 0.134, p = 0.285, n = 25), demonstrating that niche differentiation signals were driven by clonal oversampling.
* **Mobilome & Resistome:** Species-wide occurrence of intrinsic beta-lactamase `blaLAQ-1` was confirmed in 100% (44/44) of QC-passed *Lelliottia* genomes. Acquired high-risk mobile resistance elements include `blaOXA-181` on a 51.5-kb ColKP3-IncX3 plasmid in soil, `mcr-4.3` on a 123.8-kb ColE10 plasmid in swine, and IncR plasmid pP13-67 (66.8 kb, NZ_CP099512.1) in farm sewage featuring duplicated `tet(A)-floR` modules spanning a 55.3-kb composite outer region.

---

## Repository Structure

```text
Lelliottia_amnigena_research/
├── README.md                           <- Project documentation and replication guide
├── LICENSE                             <- MIT Open Source License
├── environment.yml                     <- Conda environment specification
├── requirements.txt                    <- Python dependency list
├── metadata/
│   ├── 09_Complete_57_Genomes_Catalog.csv         <- Complete metadata catalog for all 57 public assemblies
│   ├── Supplementary_Table_S1_Accounting.csv      <- Assembly accessions, origins, and taxonomic metrics
│   ├── Supplementary_Table_S2_QC_Metrics.csv       <- CheckM completeness, contamination, and QC flags
│   ├── Supplementary_Table_S3_FastANI_Matrix.csv   <- Pairwise FastANI identity matrix
│   ├── Supplementary_Table_S4_Pangenome.csv        <- Pan-genome gene frequency partition accounting
│   ├── Supplementary_Table_S5_Resistome.csv        <- Comprehensive CARD and ResFinder resistance profiles
│   ├── Supplementary_Table_S6_Plasmids_MGE.csv     <- Plasmid replicon typing and mobile genetic elements
│   └── Supplementary_Table_S7_Adaptation_Loci.csv  <- Conserved physiological adaptation and stress operons
└── scripts/
    ├── 01_prepare_genomes_and_metadata.py          <- Assembly download parsing and metadata structuring
    ├── 02_run_fastani.sh                           <- Pairwise average nucleotide identity execution
    ├── 03_run_abricate.sh                          <- Resistome and plasmid replicon screening
    ├── 04_analyze_ani_and_qc.py                    <- QC triage and species-boundary clustering
    ├── 05_analyze_abricate.py                      <- Resistome frequency synthesis and plasmid mapping
    ├── 06_prepare_panaroo_gffs.py                  <- Annotation formatting and GFF preparation
    ├── 07_run_panaroo.sh                           <- Species-core pan-genome graph reconstruction
    ├── 07b_resume_msa_and_complex.sh               <- Multi-locus sequence alignment pipeline
    ├── 08_run_iqtree.sh                            <- Maximum-likelihood phylogenomic tree inference
    ├── 09_analyze_pangenome.py                     <- Heaps' law curve fitting and gene categorization
    ├── 10_screen_antitumor_genes.py                <- Conserved physiological locus screening
    ├── 11_plot_phylogeny_and_summary.py            <- Tree and taxonomic summary figure generation
    ├── 12_deep_dive_synteny_and_mge.py             <- Plasmid pP13-67 structural and syntenic reannotation
    ├── 15_clonal_sensitivity_and_fig1_fix.py       <- Clonal dereplication sensitivity & PERMANOVA ordination
    └── download_genomes.py                         <- Automated NCBI genome retrieval utility
```

---

## Installation & Environment Setup

### 1. Clone the Repository
```bash
git clone https://github.com/asweshvaran/Lelliottia_amnigena_research.git
cd Lelliottia_amnigena_research
```

### 2. Conda Environment
Install all command-line bioinformatics tools and Python libraries via Conda:
```bash
conda env create -f environment.yml
conda activate lelliottia_pangenome
```

Alternatively, install Python dependencies via pip:
```bash
pip install -r requirements.txt
```

---

## Reproduction of Analyses

1. **Genome Retrieval:** Download public GenBank assemblies cataloged in `metadata/09_Complete_57_Genomes_Catalog.csv` using `python scripts/download_genomes.py`.
2. **Taxonomic Triage:** Run `bash scripts/02_run_fastani.sh` followed by `python scripts/04_analyze_ani_and_qc.py` to calculate pairwise ANI and evaluate CheckM quality thresholds.
3. **Resistome & Mobilome Profiling:** Execute `bash scripts/03_run_abricate.sh` to screen CARD, ResFinder, and PlasmidFinder databases; synthesize using `python scripts/05_analyze_abricate.py`.
4. **Pan-Genome Modeling:** Reconstruct gene family networks with `bash scripts/07_run_panaroo.sh` and fit Heaps' law accumulation curves using `python scripts/09_analyze_pangenome.py`.
5. **Core Phylogenomics:** Infer maximum-likelihood phylogenies with IQ-TREE using `bash scripts/08_run_iqtree.sh`.
6. **Ecological Ordination & Clonal Sensitivity:** Run `python scripts/15_clonal_sensitivity_and_fig1_fix.py` to evaluate PERMANOVA across unadjusted vs clonally dereplicated ecological cohorts.
7. **Plasmid Synteny:** Reannotate IncR plasmid pP13-67 and mobile elements using `python scripts/12_deep_dive_synteny_and_mge.py`.

---

## Citation & Contact

If you use these scripts or data in your research, please cite:

> Asweshvaran R. (2026). Pan-Genome Architecture, Cryptic Lineages, and One Health Mobilization of Epidemic Carbapenemase and Colistin Resistance Plasmids in *Lelliottia amnigena*. *Journal of Applied Microbiology* (Under Review).

For inquiries or correspondence, please contact:  
**Asweshvaran R.** — `asweshvaran@gmail.com`
