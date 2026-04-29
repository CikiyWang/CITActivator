# CITActivator

**Cap-Independent Translation Activators**

A computational toolkit for identifying and characterizing RNA-binding proteins (RBPs) that act as cap-independent translation activators. This repository contains the analysis pipelines, machine-learning models, and visualization scripts used to integrate experimental measurements (translation efficiency, luciferase activity, splicing efficiency, circRNA abundance) with eCLIP-derived binding profiles across 110 RBPs.

---

## Table of Contents

- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Requirements](#requirements)
- [Installation](#installation)
- [Modules](#modules)
  - [1. circos_heatmap](#1-circos_heatmap)
  - [2. SPLSDA_classification](#2-splsda_classification)
  - [3. RBP_eCLIP](#3-rbp_eclip)
  - [4. IRES_overlap](#4-ires_overlap)
- [Data Availability](#data-availability)
- [Citation](#citation)
- [Contact](#contact)

---

## Overview

Cap-independent translation is a critical mechanism that allows ribosomes to initiate protein synthesis under cellular stress, viral infection, or in the context of circular RNAs (circRNAs) — all without the canonical 5' cap structure. This project systematically screens **110 RNA-binding proteins (RBPs)** for their roles as cap-independent translation activators (CITActivators) by integrating four orthogonal experimental readouts with computational modeling of binding behavior.

The pipeline performs four core analyses:

1. **Multi-modal visualization** of experimental phenotypes across all 110 RBPs.
2. **Sparse PLS-DA classification** to identify protein sequence features predictive of CITActivator activity.
3. **Single-nucleotide-resolution peak coverage** from eCLIP data, normalized by transcript expression.
4. **Overlap analysis** between RBP binding peaks and annotated IRES (Internal Ribosome Entry Site) elements.

---

## Repository Structure

```
CITActivator/
├── circos_heatmap/                     # Circos plot visualization
│   └── circos_plot.R
├── SPLSDA_classification/              # sPLS-DA model training & evaluation
│   ├── splsda_train_and_evaluation.R
│   └── 110RBP_amino_seq_feature_matrix_total_length.txt
├── RBP_eCLIP/                          # eCLIP peak coverage analysis
│   └── single_nucleotide_depth_norm_single_threading.py
├── IRES_overlap/                       # IRES vs eCLIP peak overlap
│   └── peak_ires_overlap.py
└── README.md
```

---

## Requirements

### R (≥ 4.0.0)
- `circlize` — circular visualization
- `ComplexHeatmap` — heatmap rendering
- `mixOmics` — sparse PLS-DA implementation
- `caret` — cross-validation utilities
- `ggplot2`, `dplyr`, `tidyr` — data wrangling and plotting

### Python (≥ 3.7)
- `numpy`
- `pandas`
- `pysam` — BAM/SAM parsing
- `pybedtools` — interval operations
- `argparse`


---

## Installation

Clone the repository:

```bash
git clone https://github.com/CikiyWang/CITActivator.git
cd CITActivator
```

No additional compilation is required. Each module is self-contained and can be run independently.

---

## Modules

### 1. circos_heatmap

Generates a multi-track Circos heatmap integrating four experimental phenotypes — **translation efficiency**, **luciferase activity**, **splicing efficiency**, and **circRNA abundance** — measured for each of the 110 RBPs in the screen.

**Input:**
A tab-delimited matrix where rows are RBPs and columns are the four phenotypes (raw or z-score normalized values).

**Run:**

```bash
cd circos_heatmap
Rscript circos_plot.R \
    --input  ../data/110RBP_phenotype_matrix.txt \
    --output ./circos_heatmap.pdf
```

**Output:**
A circular heatmap (`circos_heatmap.pdf`) with one concentric track per phenotype, allowing visual inspection of co-modulated RBPs.

---

### 2. SPLSDA_classification

Trains a **sparse Partial Least Squares Discriminant Analysis (sPLS-DA)** model to classify RBPs as CITActivators vs. non-activators based on amino acid sequence features. The script performs **k-fold cross-validation** with feature selection inside each fold to prevent information leakage.

**Directory:** `./SPLSDA_classification`

**Input:**
- `110RBP_amino_seq_feature_matrix_total_length.txt` — RBP-by-feature matrix derived from full-length protein sequences (amino acid composition, dipeptide frequency, physicochemical descriptors, etc.).

**Run:**

```bash
cd SPLSDA_classification
Rscript splsda_train_and_evaluation.R
```

**Output:**
- Cross-validated classification accuracy, AUC, and confusion matrix.
- Ranked list of features selected in each fold.
- Per-component loading plots and sample projection plots.

---

### 3. RBP_eCLIP

Computes **single-nucleotide-resolution peak coverage** from eCLIP-seq data and normalizes it by transcript expression (TPM/FPKM from RSEM) to correct for abundance bias.

**Run:**

```bash
python single_nucleotide_depth_norm_single_threading.py \
    -depth_dir ./depth_dir/ \
    -rbp       <RBP_NAME> \
    -rep       1 \
    -i_expr    rsem.isoforms.results \
    -o         ./out_dir/
```

**Arguments:**

| Flag         | Description                                                  |
|--------------|--------------------------------------------------------------|
| `-depth_dir` | Directory containing per-nucleotide depth files (`samtools depth` output). |
| `-rbp`       | RBP identifier (e.g., `EIF4G2`).                             |
| `-rep`       | Replicate number (e.g., `1` or `2`).                         |
| `-i_expr`    | RSEM isoform-level expression file for normalization.        |
| `-o`         | Output directory.                                            |

**Output:**
A normalized coverage table (one row per transcript position) suitable for downstream meta-profile analysis.

---

### 4. IRES_overlap

Determines whether eCLIP-derived RBP binding peaks overlap with previously annotated IRES elements. For each RBP, the script reports peak counts intersecting IRES regions within a configurable flanking window.

**Directory:** `./IRES_overlap`

**Run:**

```bash
python peak_ires_overlap.py <RBP_NAME> <REGION>
```

**Arguments:**

| Positional | Description                                              |
|------------|----------------------------------------------------------|
| `RBP_NAME` | RBP identifier matching the eCLIP peak file.             |
| `REGION`   | Flanking window size in nucleotides (e.g., `150`).       |

**Example:**

```bash
python peak_ires_overlap.py LARP4 150
```

**Output:**
`{RBP}_peak_ires_overlap.txt` — a tab-delimited file listing each peak, the overlapping IRES, the overlap length, and the relative position.

---

## Data Availability

- **eCLIP peak files** were obtained from the ENCODE Project (https://www.encodeproject.org/).
- **IRES annotations** were curated from IRESbase and published literature.
- **Phenotype measurements** (translation efficiency, luciferase, splicing, circRNA) were generated in-house and are available upon request.

---

## Citation

If you use CITActivator in your research, please cite:

> Lu, Q., Wang, S., Ye, Y., Yang, Y., & Wang, Z. (2024). Systematic screen of RNA binding proteins that enhance circular RNA translation. bioRxiv, 2024-11.


---

## Contact

For questions, bug reports, or feature requests, please open an [issue](../../issues) or contact the corresponding author at `wangsiqi2019@sibs.ac.cn`.
