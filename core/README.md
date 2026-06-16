# Ferro-aging + Aortic Dissection Multi-Omics Analysis

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

Computational pipeline for the study: **"Ferro-aging Drives Vascular Smooth Muscle Cell Senescence in Aortic Dissection: A Multi-Omics Integrative Analysis"**

---

## Overview

This repository contains the complete analysis code for a multi-omics study investigating the role of **ferro-aging** — iron-driven chronic cellular senescence — in aortic dissection (AD). The pipeline integrates single-cell RNA-seq, bulk RNA-seq, drug-gene interaction databases, and cross-species aging transcriptomics.

### Key Findings
- **Ferro-aging is cell-type-specific** in human AD aorta (not global)
- **Macrophage subcluster C4** (99.7% AD-specific) identified as disease-relevant
- **SLC7A11** is the strongest differentially expressed ferro-aging gene (log2FC = +3.05 in macrophages)
- **Bulk RNA-seq** independently validates ferro-aging signature upregulation (p = 0.0043)
- **545 drug candidates** identified; 20 FDA-approved drugs prioritized
- **24 ferro-aging genes** are age-associated in primate arteries; **FOXO3** identified as master upstream regulator

---

## Pipeline Overview

```
┌─────────────────────────────────────────────────────────────────────┐
│  Step 1: Gene Set Construction                                       │
│  └─ build_ferro_aging_geneset.py → 129 genes in 7 categories        │
├─────────────────────────────────────────────────────────────────────┤
│  Step 2: scRNA-seq Analysis (GSE213740: 80,525 cells, 6 AD + 2 Ctrl)│
│  └─ ferro_aging_main_analysis.py  → cell-type atlas + DEGs          │
│  └─ ferro_aging_deep_analysis.py  → subcluster characterization     │
├─────────────────────────────────────────────────────────────────────┤
│  Step 3: Bulk RNA-seq Validation (GSE147026: 4 AD vs 4 Ctrl)        │
│  └─ ferro_aging_bulk_validation.py → cross-platform confirmation    │
├─────────────────────────────────────────────────────────────────────┤
│  Step 4: Drug Repositioning                                          │
│  └─ dgidb_query.py → DGIdb GraphQL API                              │
│  └─ ferro_aging_drug_prediction.py → prioritization + visualization  │
├─────────────────────────────────────────────────────────────────────┤
│  Step 5: Age Validation (Primate Arterial Aging, Nat Commun 2020)    │
│  └─ parse_primate_aging.py → download supplementary data            │
│  └─ ferro_aging_age_validation.py → cross-species cross-validation  │
├─────────────────────────────────────────────────────────────────────┤
│  Step 6: Manuscript Figures                                          │
│  └─ manuscript_figures.py → Figure 1–3 (publication-ready)          │
└─────────────────────────────────────────────────────────────────────┘
```

---

## Data Sources

| Dataset | Accession | Type | Samples | Usage |
|---------|-----------|------|---------|-------|
| Human AD aorta scRNA-seq | [GSE213740](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213740) | 10X scRNA-seq | 6 AD + 2 Normal | Primary discovery |
| Human AD aorta bulk RNA-seq | [GSE147026](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147026) | Bulk mRNA-seq | 4 AD + 4 Control | Independent validation |
| Primate arterial aging scRNA-seq | [Nat Commun 2020](https://doi.org/10.1038/s41467-020-15997-0) | 10X scRNA-seq | Young vs Old | Age cross-validation |
| Drug-Gene Interaction Database | [DGIdb v5.0](https://www.dgidb.org/) | API | — | Drug repositioning |

---

## Installation

```bash
# Clone repository
git clone https://github.com/YOUR_USERNAME/ferro-aging-aad.git
cd ferro-aging-aad

# Create conda environment (recommended)
conda create -n ferro_aging python=3.10 -y
conda activate ferro_aging

# Install dependencies
pip install -r requirements.txt
```

> **Note:** Some scripts require `scanpy` which is best installed via conda:
> ```bash
> conda install -c conda-forge scanpy python-igraph leidenalg -y
> ```

---

## Usage

Scripts should be run in the following order:

```bash
# 1. Build gene set
python scripts/build_ferro_aging_geneset.py

# 2. scRNA-seq primary analysis (requires GSE213740 data in sc_data/)
python scripts/ferro_aging_main_analysis.py

# 3. Deep subcluster analysis
python scripts/ferro_aging_deep_analysis.py

# 4. Bulk RNA-seq validation
python scripts/ferro_aging_bulk_validation.py

# 5. Drug-gene interaction query
python scripts/dgidb_query.py

# 6. Drug prediction + visualization
python scripts/ferro_aging_drug_prediction.py

# 7. Download & parse primate aging data
python scripts/parse_primate_aging.py

# 8. Age validation analysis
python scripts/ferro_aging_age_validation.py

# 9. Generate manuscript figures
python scripts/manuscript_figures.py
```

> ⚠️ **Path Configuration:** Before running, update the `WORK_DIR`, `DATA_DIR`, and `OUT_DIR` variables at the top of each script to match your local directory structure.

---

## Output Files

All results are written to the `results/` directory:

| File | Description |
|------|-------------|
| `ferro_aging_geneset.csv` | 129-gene set with categories & annotations |
| `GSE213740_ferro_aging_processed.h5ad` | Processed scRNA-seq AnnData object |
| `ferro_aging_gene_DE.csv` | Differential expression results |
| `ferro_aging_cell_scores.csv` | Per-cell ferro-aging module scores |
| `ferro_aging_bulk_validation.csv` | Bulk RNA-seq cross-validation |
| `ferro_aging_drug_rankings.csv` | Prioritized drug candidates |
| `ferro_aging_primate_age_genes.csv` | Age-associated ferro-aging genes |
| `ferro_aging_*.png` | Publication-quality figures |

---

## Ferro-aging Gene Set (129 genes, 7 categories)

| Category | Count | Key Genes |
|----------|-------|-----------|
| Core_FerroAging | 7 | ACSL4, GPX4, NCOA4, FTH1, HMOX1, TFRC, SLC40A1 |
| Iron_Metabolism | 35 | TF, FTL, HFE, HAMP, BMP6, STEAP3... |
| Ferroptosis | 28 | SLC7A11, SLC3A2, AIFM2, LPCAT3, PTGS2... |
| Lipid_Peroxidation | 12 | ALOX5, ALOX12, ALOX15, NOX4, CYBB... |
| NRF2_Pathway | 11 | NFE2L2, KEAP1, GCLM, GCLC, NQO1... |
| Senescence_SASP | 18 | CDKN1A, CDKN2A, TP53, CXCL8, IL6, SERPINE1... |
| Vascular_ECM | 18 | MMP2, MMP9, MMP12, FBN1, COL1A1, COL3A1... |

---

## Citation

If you use this code, please cite:

> *[Manuscript under review]*

and the data sources:

- Zhang et al. (2023) *Biomolecules* 13(2):399 — GSE213740
- Liu et al. (2026) *Cell Metabolism* — Ferro-aging concept
- Nat Commun (2020) 11:2628 — Primate arterial aging
- Cannon et al. (2024) *Nucleic Acids Res* 52(D1):D1227 — DGIdb v5.0

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

---

## Contact

For questions, please open an issue or contact the corresponding author.
