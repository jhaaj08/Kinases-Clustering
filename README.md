# Kinase Functional Classification with ESM-2 Layer Selection

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17370925.svg)](https://doi.org/10.5281/zenodo.17370925)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.12+](https://img.shields.io/badge/python-3.12+-blue.svg)](https://www.python.org/downloads/)

> **Key Finding**: Mid-layer averaging (layers 20–30) in ESM-2 improves unsupervised kinase clustering by **+138% ARI** over final-layer representations, while the final layer (33) performs best for supervised classification (79.9% calibrated accuracy).

---

## For Reviewers — Complete Reproduction Guide

### Step 1: Clone the repository

```bash
git clone https://github.com/jhaaj08/Kinases-Clustering.git
cd Kinases-Clustering
```

### Step 2: Install dependencies

```bash
# Option A: Conda (recommended)
conda env create -f environment.yml
conda activate kinase-clustering
conda install -c bioconda hmmer cd-hit

# Option B: pip
python -m venv venv && source venv/bin/activate
pip install -r requirements.txt
brew install hmmer cd-hit        # macOS
# sudo apt-get install hmmer cd-hit  # Ubuntu
```

### Step 3: Download pre-computed embeddings from Zenodo

The ESM-2 embeddings (~7 GB) are not in the repository. Download them from Zenodo:

**DOI**: [10.5281/zenodo.17370925](https://doi.org/10.5281/zenodo.17370925)

```bash
wget https://zenodo.org/records/17370925/files/kinase_data_v1.zip
unzip kinase_data_v1.zip
cp -r kinase_data_v1/embeddings/* embeddings/
```

### Step 4: Run

```bash
make review
```

That's it. Expected runtime: ~25 min (CPU) / ~15 min (GPU).

### Step 5: Open results

```bash
open runs/run_current_data/REPORT.html    # macOS
# xdg-open runs/run_current_data/REPORT.html  # Linux
```

The `REPORT.html` is self-contained — all figures are embedded, all tables are rendered, all key numbers are shown. No internet connection required to view it.

### Step 6: Verify integrity

```bash
make verify RUN_ID=run_current_data
```

Expected output: all checks PASS (SHA256 hashes, manifest counts, split integrity).

---

### What gets created

```
runs/run_current_data/
├── REPORT.html                        ← open this — all figures, tables, numbers
├── results/manuscript_numbers.json   ← every number cited in the manuscript
├── tables/
│   ├── Table1_supervised.csv         ← Table 1: classification performance (8 methods)
│   ├── TableS1.csv                   ← Table S1: clustering layer ablation
│   ├── TableS2.csv                   ← Table S2: baseline comparisons
│   └── Table1.csv                    ← dataset construction summary
├── figures/
│   ├── Figure1_clustering_ari.png
│   ├── Figure2_confusion_matrix.png
│   ├── Figure3_homology_classification.png
│   ├── Figure4_pooling_comparison.png
│   ├── Figure5_calibration.png
│   └── Figure6_retrieval_pr.png
├── MANIFEST.txt                       ← SHA256 hashes for all files
└── data/, results/, models/, ...     ← full run artifacts
```

> **Note**: `runs/` is gitignored — each reviewer creates their own local run.

---

## Expected Key Results

| Metric | Expected Value |
|--------|----------------|
| Best Clustering ARI (layers 20–30) | ~0.305 |
| ARI improvement vs. Layer 33 | ~+138% |
| Supervised accuracy — Layer 33, calibrated, split40 | ~79.9% |
| Calibrated ECE | ~0.095 |
| Retrieval P@1 | ~0.703 |
| Retrieval MRR | ~0.781 |
| Supervised dataset size | 1,362 sequences, 8 classes |

---

## Overview

This repository provides a reproducible pipeline for classifying human protein kinases into functional families using ESM-2 protein language model embeddings. The central finding is that different transformer layers encode qualitatively different information:

- **Clustering**: Intermediate layers (20–30) work best — +138% ARI over the final layer
- **Classification**: Final layer (33) works best — 79.9% calibrated accuracy

### What the pipeline does

1. Creates dataset manifests from pre-extracted kinase domain sequences
2. Links pre-computed ESM-2 embeddings
3. Runs k-means clustering across 4 layer configurations
4. Creates homology-aware train/test splits (40%, 50%, 70% identity thresholds)
5. Trains logistic regression classifiers with Platt scaling calibration
6. Computes baselines (k-NN, MLP, motifs-only LR, random)
7. Computes retrieval metrics (P@k, MRR)
8. Generates manuscript tables and 6 publication figures
9. Produces a self-contained `REPORT.html` with all results

---

## Installation

### Option A: Conda (recommended)

```bash
conda env create -f environment.yml
conda activate kinase-clustering
conda install -c bioconda hmmer cd-hit
```

### Option B: pip

```bash
python -m venv venv && source venv/bin/activate
pip install -r requirements.txt
brew install hmmer cd-hit        # macOS
# sudo apt-get install hmmer cd-hit  # Ubuntu
```

### Verify

```bash
python -c "import torch; from sklearn.cluster import KMeans; print('OK')"
hmmsearch -h | head -1
cd-hit -h | head -1
```

---

## Running the Pipeline

### All Makefile targets

```bash
make review                     # Reviewer shortcut — runs everything as run_current_data
make all                        # Fresh run with auto-generated timestamp ID
make all RUN_ID=my_exp          # Named run
make all RUN_ID=my_exp FORCE=1  # Overwrite existing run
make verify RUN_ID=my_exp       # Integrity check
make zip                        # Create Zenodo ZIP package
make list                       # List all runs
make clean RUN_ID=my_exp        # Remove a specific run
```

### Individual steps (if needed)

```bash
make manifests   RUN_ID=x   # Step 6:  dataset manifests
make embeddings  RUN_ID=x   # Step 8:  link embeddings
make splits      RUN_ID=x   # Step 10: homology-aware splits
make clustering  RUN_ID=x   # Step 9:  k-means clustering
make supervised  RUN_ID=x   # Step 11: logistic regression
make calibration RUN_ID=x   # Step 12: Platt scaling
make baselines   RUN_ID=x   # Step 13: baseline methods
make retrieval   RUN_ID=x   # Step 14: retrieval metrics
make tables      RUN_ID=x   # Step 15: manuscript numbers + tables
make figures     RUN_ID=x   # Step 16: publication figures
make report      RUN_ID=x   # Step 17: REPORT.html
```

---

## Project Structure

```
Kinases-Clustering/
├── pipeline/                    # Reproducible pipeline steps
│   ├── step_06_manifests.py     # Dataset manifests
│   ├── step_08_embeddings.py    # Link embeddings
│   ├── step_09_clustering.py    # K-means clustering
│   ├── step_10_splits.py        # Homology-aware splits
│   ├── step_11_supervised.py    # Logistic regression
│   ├── step_12_calibration.py   # Platt scaling
│   ├── step_13_baselines.py     # Baseline methods
│   ├── step_14_retrieval.py     # Retrieval metrics
│   ├── step_15_build_numbers.py # Manuscript tables + numbers JSON
│   ├── step_16_figures.py       # Publication figures
│   ├── step_17_report.py        # REPORT.html generator
│   ├── generate_manifest.py     # SHA256 manifest
│   ├── membership.py            # Dataset membership validation
│   └── run_manager.py           # Run directory management
│
├── scripts/                     # Data preparation (raw → embeddings)
│   ├── download_uniprot_kinases.py
│   ├── filter_sequences.py
│   ├── deduplicate_sequences.py
│   ├── reduce_redundancy_cdhit.py
│   ├── extract_domains.py
│   ├── assign_labels.py
│   ├── generate_embeddings.py
│   ├── regenerate_embeddings.py
│   └── verify_package.py        # Run integrity checks
│
├── data/                        # Source data (not re-run by default)
│   ├── raw/                     # Original UniProt downloads
│   ├── processed/               # Cleaned sequences + labels
│   ├── domains/                 # Extracted kinase domains (HMMER)
│   ├── manifests/               # Dataset membership (source of truth)
│   ├── splits/                  # Pre-computed train/test splits
│   └── hmm_profiles/            # Pfam HMM files (PF00069, PF07714)
│
├── embeddings/
│   └── esm2_t33_650M/           # Pre-computed ESM-2 embeddings
│       ├── ids.txt
│       ├── embedding_metadata.json
│       └── *.npy                # Embedding arrays
│
├── runs/                        # All experiment outputs (gitignored)
│   ├── current -> run_current_data/   # symlink to latest
│   └── run_current_data/        # Created by `make review`
│       ├── REPORT.html          # Self-contained reviewer report
│       ├── results/manuscript_numbers.json
│       ├── tables/              # Table1_supervised, TableS1, TableS2
│       ├── figures/             # Figure1–Figure6 PNGs
│       └── ...
│
├── figures_output/              # Stable project-level figure copies
│   └── Figure1–Figure6.png
│
├── webapp/                      # Gradio web application
│   ├── app.py
│   └── predictor.py
│
├── docs/
│   └── Simple_English.md        # Plain-language explanation
│
├── MANUSCRIPT.md                # Full manuscript text
├── Makefile
├── requirements.txt
└── environment.yml
```

---

## Data Description

### Dataset Summary

| Stage | Sequences | Classes | Notes |
|-------|-----------|---------|-------|
| UniProt download | 20,262 | — | All reviewed human kinases |
| After cleaning | 6,465 | 11 | Deduplicated, CD-HIT 60% |
| Domain extraction (E < 0.01) | 1,959 | 11 | HMMER + Pfam (PF00069, PF07714) |
| Clustering dataset | 1,387 | 10 | Excluding "Other" |
| Supervised dataset | 1,362 | 8 | Excluding Other, Histidine, RGC |

**Why 8 classes for supervised learning:**
- **Histidine kinases** (7 sequences): Bacterial, mechanistically distinct
- **RGC** (18 sequences): Receptor guanylate cyclases, not true kinases
- Both are retained for clustering (10 classes) but excluded from supervised learning

### Kinase Classes (supervised, n = 1,362)

| Family | Count | Description |
|--------|-------|-------------|
| TK | 490 | Tyrosine kinases |
| CMGC | 240 | CDK, MAPK, GSK3, CLK |
| CAMK | 221 | Calcium/calmodulin-dependent |
| AGC | 139 | PKA, PKG, PKC |
| STE | 130 | MAP kinase cascade |
| TKL | 63 | Tyrosine kinase-like |
| CK1 | 43 | Casein kinase 1 |
| Atypical | 36 | PI3K, mTOR, etc. |

---

## Reproducibility

### Guarantees

- **No stale outputs**: Pipeline aborts if run directory exists (use `FORCE=1` to overwrite)
- **Single source of truth**: All datasets derived from `data/manifests/*.txt`
- **Homology leakage prevention**: CD-HIT ensures no train sequence is >40% similar to any test sequence
- **Fixed seeds**: All experiments use `random_state=42`, `PYTHONHASHSEED=0`

### Reproducibility across runs

| Result | Reproducibility | Notes |
|--------|----------------|-------|
| Clustering (ARI, NMI) | 100% identical | k-means with seed=42 |
| Layer 33 supervised accuracy | 100% identical | Main result |
| Retrieval (P@k, MRR) | 100% identical | Deterministic k-NN |
| Dataset counts | 100% identical | Deterministic filtering |
| Layers 20–30 supervised | ±0.8% | LBFGS convergence variation |
| MLP baseline | ±1.5% | Neural net randomness |
| Random baseline | ±3.5% | Stratified sampling |

All primary claims are 100% reproducible.

### Calibrated vs uncalibrated accuracy

The primary reported accuracy (79.9%) uses Platt scaling via `CalibratedClassifierCV(cv=5)`. This builds an ensemble of 5 sub-models, which improves both accuracy and probability reliability. The uncalibrated logistic regression alone achieves 73.6%.

---

## Figures

| Figure | Content |
|--------|---------|
| Figure 1 | ARI bar chart — 4 ESM-2 layer configurations |
| Figure 2 | 8×8 confusion matrix (split40, Layer 33) |
| Figure 3 | Calibrated accuracy vs. identity threshold (70/50/40%) |
| Figure 4 | Mean pooling vs CLS token — clustering and classification |
| Figure 5 | Reliability diagram before/after Platt scaling |
| Figure 6 | Precision-recall curve at cosine-similarity thresholds |

All figures are generated by `pipeline/step_16_figures.py` at 300 DPI and saved into `{run_dir}/figures/`. They are also embedded in `REPORT.html`.

---

## Expected Runtime

| Step | CPU | GPU |
|------|-----|-----|
| Manifests | <1 min | <1 min |
| Embeddings (link) | <1 min | <1 min |
| Splits (CD-HIT) | ~5 min | ~5 min |
| Clustering | ~2 min | ~1 min |
| Supervised + Calibration | ~5 min | ~2 min |
| Baselines | ~5 min | ~2 min |
| Retrieval | ~2 min | ~1 min |
| Tables + Figures + Report | ~3 min | ~2 min |
| **Total** | **~25 min** | **~15 min** |

---

## Web Application

```bash
python webapp/app.py
# → http://localhost:7860
```

Accepts any kinase sequence: extracts domain (HMMER), generates ESM-2 embedding, predicts family with confidence scores.

---

## Troubleshooting

| Issue | Solution |
|-------|----------|
| `ModuleNotFoundError: esm` | `pip install fair-esm` |
| `cd-hit: command not found` | `conda install -c bioconda cd-hit` or `brew install cd-hit` |
| `hmmsearch: command not found` | `conda install -c bioconda hmmer` or `brew install hmmer` |
| Run directory already exists | `make all RUN_ID=xxx FORCE=1` |
| Results differ slightly | Check Python 3.12+ and `requirements.txt` versions |

---

## Citation

```bibtex
@article{kinase_layer_selection_2025,
  title   = {Layer Selection in Protein Language Models Improves Kinase Functional Classification},
  author  = {[Authors]},
  year    = {2025},
  doi     = {10.5281/zenodo.17370925}
}
```

---

## License

MIT — see [LICENSE](LICENSE).

## Acknowledgments

- **ESM-2**: Meta AI Research ([fair-esm](https://github.com/facebookresearch/esm))
- **Data**: UniProt, Pfam (InterPro)
- **Tools**: HMMER, CD-HIT, scikit-learn, PyTorch
