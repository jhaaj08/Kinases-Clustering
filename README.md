# Kinase Functional Classification with ESM-2 Layer Selection

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17370925.svg)](https://doi.org/10.5281/zenodo.17370925)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.12+](https://img.shields.io/badge/python-3.12+-blue.svg)](https://www.python.org/downloads/)

> **Key Finding**: Mid-layer averaging (layers 20-30) in ESM-2 improves unsupervised kinase clustering by **+138% ARI** over final layer representations, while final layer embeddings perform better for supervised classification.

---

## Table of Contents

- [Overview](#overview)
- [Key Results](#key-results)
- [Quick Start](#quick-start)
- [Installation](#installation)
- [Running the Pipeline](#running-the-pipeline)
- [Project Structure](#project-structure)
- [Data Description](#data-description)
- [Reproducibility](#reproducibility)
- [Web Application](#web-application)
- [Citation](#citation)
- [License](#license)

---

## Overview

This repository contains a complete, reproducible pipeline for classifying protein kinases into functional families using ESM-2 protein language model embeddings. The project investigates how different transformer layers encode functional information, revealing that intermediate layers (20-30) capture better clustering structure while final layers excel at supervised classification.

### What This Project Does

1. **Downloads and cleans** kinase sequences from UniProt (6,465 → 1,959 domain sequences)
2. **Extracts catalytic domains** using HMMER + Pfam profiles (PF00069, PF07714)
3. **Generates ESM-2 embeddings** with multiple layer configurations
4. **Performs unsupervised clustering** (k-means, 10 kinase families)
5. **Trains supervised classifiers** with homology-aware train/test splits
6. **Evaluates** clustering, classification, calibration, and retrieval performance
7. **Provides a web interface** for classifying new kinase sequences

---

## Key Results

| Task | Metric | Value | Details |
|------|--------|-------|---------|
| **Unsupervised Clustering** | ARI | 0.305 | Layers 20-30, k=10, n=1,387 |
| **Clustering Improvement** | vs. Final Layer | +138.4% | 0.128 → 0.305 ARI |
| **Supervised Classification** | Accuracy | 80.2% | Layer 33, 40% identity split |
| **Classification** | Macro-F1 | 0.617 | 8-class problem, n=1,362 |
| **Retrieval** | P@1 | 0.703 | Nearest-neighbor |
| **Retrieval** | MRR | 0.781 | Mean Reciprocal Rank |

### The Layer Selection Paradox

We discovered an interesting finding:
- **Clustering**: Intermediate layers (20-30) work best (+138% ARI improvement)
- **Classification**: Final layer (33) works best (80.2% vs 73.6% accuracy)

This suggests different layers encode different types of information useful for different tasks.

---

## Quick Start

### 30-Second Demo

```bash
# Clone and setup
git clone https://github.com/jhaaj08/Kinases-Clustering.git
cd Kinases-Clustering
pip install -r requirements.txt

# Run a fresh experiment
make all

# Verify results
make verify
```

### What Happens

1. Creates a timestamped run directory: `runs/2025-12-25_120000/`
2. Generates dataset manifests from domain coordinates
3. Links pre-computed embeddings
4. Runs k-means clustering experiments
5. Creates homology-aware train/test splits
6. Trains logistic regression classifiers
7. Computes calibration, baselines, and retrieval metrics
8. Generates manuscript tables
9. **Generates all manuscript figures** (UMAP, clustering, supervised, calibration, retrieval)
10. Creates SHA256 manifest for integrity verification

---

## Installation

### Prerequisites

- Python 3.12+
- Conda (recommended) or pip
- ~50GB disk space (for embeddings)
- 16GB RAM minimum

### Option 1: Conda (Recommended)

```bash
# Create environment
conda env create -f environment.yml
conda activate kinase-clustering

# Install system dependencies
conda install -c bioconda hmmer cd-hit
```

### Option 2: Pip

```bash
# Create virtual environment
python -m venv venv
source venv/bin/activate  # Linux/Mac
# or: venv\Scripts\activate  # Windows

# Install dependencies
pip install -r requirements.txt

# Install system tools (macOS)
brew install hmmer cd-hit

# Install system tools (Ubuntu/Debian)
sudo apt-get install hmmer cd-hit
```

### Verify Installation

```bash
python -c "import torch; import esm; print('PyTorch:', torch.__version__)"
python -c "from sklearn.cluster import KMeans; print('scikit-learn OK')"
hmmsearch -h | head -2
cd-hit -h | head -2
```

---

## Running the Pipeline

### Method 1: Makefile (Recommended)

The Makefile provides dependency-tracked, reproducible builds.

```bash
# Fresh run with auto-generated timestamp ID
make all

# Named run for experiments
make all RUN_ID=experiment_v1

# Force overwrite existing run
make all RUN_ID=experiment_v1 FORCE=1

# Run specific steps
make manifests      # Step 6: Create dataset manifests
make embeddings     # Step 8: Link embeddings
make splits         # Step 10: Create homology-aware splits
make clustering     # Step 9: Run k-means clustering
make supervised     # Step 11: Train classifiers
make calibration    # Step 12: Calibration
make baselines      # Step 13: Baseline comparisons
make retrieval      # Step 14: Retrieval experiment
make tables         # Step 15: Generate manuscript tables
make figures        # Step 16: Generate manuscript figures

# Verify and package
make verify         # Run all integrity checks
make zip            # Create Zenodo package
make list           # List all runs
make clean RUN_ID=xxx  # Remove a specific run
```

### Method 2: Python Scripts

```bash
# Initialize a run directory
python pipeline/run_manager.py --init --run-id my_experiment

# Run individual steps
python pipeline/step_06_manifests.py --run-dir runs/my_experiment
python pipeline/step_08_embeddings.py --run-dir runs/my_experiment
python pipeline/step_10_splits.py --run-dir runs/my_experiment
python pipeline/step_09_clustering.py --run-dir runs/my_experiment
python pipeline/step_11_supervised.py --run-dir runs/my_experiment
python pipeline/step_12_calibration.py --run-dir runs/my_experiment
python pipeline/step_13_baselines.py --run-dir runs/my_experiment
python pipeline/step_14_retrieval.py --run-dir runs/my_experiment
python pipeline/step_15_build_numbers.py --run-dir runs/my_experiment
python pipeline/step_16_figures.py --run-dir runs/my_experiment
python pipeline/generate_manifest.py --run-dir runs/my_experiment

# Verify
python scripts/verify_package.py runs/my_experiment
```

### Understanding Run Directories

Each run creates an isolated directory:

```
runs/2025-12-25_120000/
├── run_config.json           # Run metadata
├── data/
│   ├── manifests/            # Dataset membership (SOURCE OF TRUTH)
│   │   ├── domain_E001.txt   # 1,387 IDs for clustering
│   │   ├── supervised_eligible.txt  # 1,362 IDs for classification
│   │   └── manifest_report.json
│   └── splits/               # Train/test splits
│       ├── split40_train.txt # 1,089 training IDs
│       ├── split40_test.txt  # 273 test IDs
│       └── splits_report.json
├── embeddings/
│   └── esm2_t33_650M/        # Symlinked embeddings
│       ├── ids.txt           # 1,959 embedding IDs
│       └── *.npy → ../../..  # Symlinks to project embeddings
├── results/
│   ├── clustering/           # K-means results
│   ├── supervised/           # Classification results
│   ├── calibration/          # Platt scaling results
│   ├── baselines/            # Comparison methods
│   ├── retrieval/            # NN retrieval results
│   └── manuscript_numbers.json  # All numbers in one place
├── tables/
│   ├── Table1.csv            # Dataset construction
│   ├── TableS1.csv           # Layer ablation
│   └── TableS2.csv           # Baselines comparison
├── figures/                   # Generated figures
│   ├── Fig2_umap_geometry.png
│   ├── Fig3_clustering_metrics.png
│   ├── Fig4_supervised_homology_splits.png
│   ├── Fig5_calibration_reliability.png
│   ├── Fig6_retrieval_metrics.png
│   ├── FigS1_layer_sweep_clustering.png
│   ├── FigS2_dataset_class_distribution.png
│   └── figure_registry.json
├── MANIFEST.txt              # SHA256 hashes
└── MANIFEST.json             # Machine-readable manifest
```

---

## Project Structure

```
Kinases-Clustering/
├── pipeline/                 # NEW: Reproducible pipeline module
│   ├── __init__.py          # Package exports
│   ├── run_manager.py       # Run directory management
│   ├── membership.py        # Dataset membership validation
│   ├── step_06_manifests.py # Create dataset manifests
│   ├── step_08_embeddings.py# Link embeddings
│   ├── step_09_clustering.py# K-means clustering
│   ├── step_10_splits.py    # Homology-aware splits
│   ├── step_11_supervised.py# Logistic regression
│   ├── step_12_calibration.py# Platt scaling
│   ├── step_13_baselines.py # Baseline methods
│   ├── step_14_retrieval.py # NN retrieval
│   ├── step_15_build_numbers.py # Generate tables
│   └── generate_manifest.py # SHA256 manifest
├── scripts/                  # Utility scripts
│   ├── verify_package.py    # Package verification
│   ├── regenerate_embeddings.py # ESM-2 embedding generation
│   └── ...
├── data/                     # Source data
│   ├── raw/                  # Original UniProt data
│   ├── processed/            # Cleaned data + labels
│   ├── domains/              # Extracted kinase domains
│   ├── manifests/            # Dataset membership files
│   ├── splits/               # Train/test splits
│   └── hmm_profiles/         # Pfam HMM files
├── embeddings/               # Pre-computed embeddings
│   └── esm2_t33_650M/
│       ├── domain_E001_layer33_mean.npy
│       ├── domain_E001_layers20_30_mean.npy
│       ├── ids.txt
│       └── embedding_metadata.json
├── runs/                     # Experiment runs (gitignored)
│   ├── .gitkeep
│   ├── current -> 2025-12-25_120000/  # Symlink to latest
│   └── 2025-12-25_120000/    # Timestamped runs
├── webapp/                   # Gradio web application
├── Makefile                  # Build system
├── requirements.txt          # Python dependencies
├── environment.yml           # Conda environment
├── MANUSCRIPT.md             # Full manuscript
└── README.md                 # This file
```

---

## Data Description

### Dataset Summary

| Stage | Sequences | Classes | Notes |
|-------|-----------|---------|-------|
| UniProt download | 20,262 | - | All reviewed kinases |
| After cleaning | 6,465 | 11 | Deduplicated, CD-HIT 60% |
| Domain extraction (E<0.01) | 1,959 | 11 | HMMER + Pfam |
| Clustering dataset | 1,387 | 10 | Excluding "Other" |
| Supervised dataset | 1,362 | 8 | Excluding Other, Histidine, RGC |

### Kinase Families (8 Classes for Supervised Learning)

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

### Why Histidine and RGC are Excluded

- **Histidine kinases** (7 sequences): Bacterial, mechanistically distinct from eukaryotic kinases
- **RGC** (18 sequences): Receptor guanylate cyclases, not true kinases

Both are retained for clustering (10 classes) but excluded from supervised learning (8 classes).

---

## Reproducibility

### Verification

Every run can be verified for integrity:

```bash
# Verify a run directory
python scripts/verify_package.py runs/2025-12-25_120000/

# Verify a ZIP package
python scripts/verify_package.py kinase_data_v1.zip
```

Verification checks:
- ✓ SHA256 hashes match MANIFEST.txt
- ✓ Manifest counts are valid (domain ≥ supervised)
- ✓ Split integrity (train + test = supervised)
- ✓ Tables match manuscript_numbers.json
- ✓ No orphan IDs

### Key Guarantees

1. **No stale outputs**: Pipeline aborts if run directory exists (unless `FORCE=1`)
2. **Single source of truth**: All datasets derived from `data/manifests/*.txt`
3. **Homology leakage prevention**: CD-HIT clustering ensures no sequence in train is >40% similar to test
4. **Fail-fast assertions**: Pipeline stops immediately on integrity violations

### Random Seeds

All experiments use fixed seeds for reproducibility:
- Random seed: 42
- K-means n_init: 10
- Train/test split ratio: 80/20

### Version Pinning

Key package versions for exact reproduction:
- Python: 3.12+
- PyTorch: 2.0+
- fair-esm: 2.0.0
- scikit-learn: 1.7.1
- numpy: 1.26+

All versions are locked in `requirements.txt`. For exact reproduction, use:
```bash
pip install -r requirements.txt --no-deps
```

### Reproducibility Validation

We validated reproducibility across multiple runs. Results:

| Metric Category | Reproducibility | Notes |
|-----------------|-----------------|-------|
| **Clustering (ARI, NMI)** | 100% identical | k-means with seed=42 |
| **Dataset counts** | 100% identical | Deterministic filtering |
| **Train/test splits** | 100% identical | CD-HIT + fixed seed |
| **Layer 33 supervised** | 100% identical | Main result |
| **Retrieval (P@k, MRR)** | 100% identical | k-NN is deterministic |
| **k-NN baseline** | 100% identical | Deterministic |
| Layers 20-30 supervised | ±0.8% | LBFGS convergence variation |
| MLP baseline | ±1.5% | Neural net randomness |
| Random baseline | ±3.5% | Stratified sampling |

**All primary claims are 100% reproducible.**

The Makefile automatically sets environment variables for maximum determinism:
```bash
export PYTHONHASHSEED=0
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
```

---

## Web Application

Launch the interactive web interface:

```bash
python webapp/app.py
```

Open: http://localhost:7860

### Features

- Paste any kinase sequence
- Automatic domain extraction (HMMER)
- ESM-2 embedding generation
- Family prediction with confidence scores
- Motif highlighting (DFG, HRD, APE)
- Similar kinase retrieval

### Example

```python
from webapp.predictor import KinasePredictor

predictor = KinasePredictor(
    model_path='supervised_results/logistic_regression_model.joblib',
    embeddings_path='embeddings/esm2_t33_650M/domain_E001_layer33_mean.npy',
    index_path='embeddings/esm2_t33_650M/ids.txt'
)

result = predictor.predict("MENFQKVEKIGEGTYGVVYKARNKLT...")
print(f"Predicted: {result['top_predictions'][0]}")
print(f"Confidence: {result['confidence_flag']}")
```

---

## For Reviewers

### Complete Reproduction Guide

Follow these exact steps to reproduce all results from scratch:

```bash
# ============================================================
# STEP 1: Clone repository
# ============================================================
git clone https://github.com/jhaaj08/Kinases-Clustering.git
cd Kinases-Clustering

# ============================================================
# STEP 2: Install dependencies
# ============================================================
# Option A: Conda (recommended)
conda env create -f environment.yml
conda activate kinase-clustering
conda install -c bioconda hmmer cd-hit

# Option B: Pip
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
brew install hmmer cd-hit  # macOS
# sudo apt-get install hmmer cd-hit  # Ubuntu

# ============================================================
# STEP 3: Verify installation
# ============================================================
python -c "import torch; import esm; print('PyTorch:', torch.__version__)"
python -c "from sklearn.cluster import KMeans; print('scikit-learn OK')"
hmmsearch -h | head -2
cd-hit -h | head -2

# ============================================================
# STEP 4: Run complete pipeline
# ============================================================
make all RUN_ID=review_run

# ============================================================
# STEP 5: Verify results
# ============================================================
make verify RUN_ID=review_run

# Expected output:
# ============================================================
# VERIFICATION SUMMARY
# ============================================================
# [PASS] ✓ MANIFEST.txt exists with SHA256 hashes
# [PASS] ✓ Manifest counts valid (domain >= supervised)
# [PASS] ✓ Split integrity (train + test = supervised)
# [PASS] ✓ All 8 classes present in supervised dataset
# [PASS] ✓ No orphan IDs in splits
# 
# VERIFICATION PASSED - Package is valid

# ============================================================
# STEP 6: Compare results to expected values
# ============================================================
cat runs/review_run/results/manuscript_numbers.json | python -m json.tool
```

### Expected Key Results

After running the pipeline, verify these numbers match (±0.001 tolerance due to floating point):

| Metric | Expected Value | JSON Path |
|--------|----------------|-----------|
| **Best Clustering ARI** | 0.305 | `clustering.best_ARI` |
| **Layer 33 Clustering ARI** | 0.128 | `clustering.baseline_ARI` |
| **Improvement %** | 138.4% | `clustering.improvement_percent` |
| **Layer 33 Accuracy (split40)** | 0.802 | `supervised.split40.layer33_mean.accuracy` |
| **Layer 33 Macro-F1 (split40)** | 0.617 | `supervised.split40.layer33_mean.macro_f1` |
| **Retrieval P@1** | 0.703 | `retrieval.P@1` |
| **Retrieval MRR** | 0.781 | `retrieval.MRR` |
| **Clustering N** | 1,387 | `dataset.domain_E001_n` |
| **Supervised N** | 1,362 | `dataset.supervised_eligible_n` |
| **Train/Test Split (40%)** | 1,089 / 273 | `splits.split40.n_train/n_test` |

### Quick Verification Script

```bash
# One-liner to check key results
python -c "
import json
with open('runs/review_run/results/manuscript_numbers.json') as f:
    d = json.load(f)
print('=== KEY RESULTS ===')
print(f'Best ARI:       {d[\"clustering\"][\"best_ARI\"]:.3f} (expect 0.305)')
print(f'Improvement:    {d[\"clustering\"][\"improvement_percent\"]:.1f}% (expect 138.4%)')
print(f'Accuracy:       {d[\"supervised\"][\"split40\"][\"layer33_mean\"][\"accuracy\"]:.3f} (expect 0.802)')
print(f'Macro-F1:       {d[\"supervised\"][\"split40\"][\"layer33_mean\"][\"macro_f1\"]:.3f} (expect 0.617)')
print(f'P@1:            {d[\"retrieval\"][\"P@1\"]:.3f} (expect 0.703)')
print(f'MRR:            {d[\"retrieval\"][\"MRR\"]:.3f} (expect 0.781)')
print(f'Clustering N:   {d[\"dataset\"][\"domain_E001_n\"]} (expect 1387)')
print(f'Supervised N:   {d[\"dataset\"][\"supervised_eligible_n\"]} (expect 1362)')
"
```

### Expected Runtime

| Step | Time (CPU) | Time (GPU) |
|------|------------|------------|
| Manifests | <1 min | <1 min |
| Embeddings (link) | <1 min | <1 min |
| Splits (CD-HIT) | ~5 min | ~5 min |
| Clustering | ~2 min | ~1 min |
| Supervised | ~3 min | ~1 min |
| Calibration | ~2 min | ~1 min |
| Baselines | ~5 min | ~2 min |
| Retrieval | ~2 min | ~1 min |
| Tables | <1 min | <1 min |
| Figures | ~3 min | ~2 min |
| **Total** | **~24 min** | **~15 min** |

### Troubleshooting

| Issue | Solution |
|-------|----------|
| `ModuleNotFoundError: esm` | `pip install fair-esm` |
| `cd-hit: command not found` | `conda install -c bioconda cd-hit` or `brew install cd-hit` |
| `hmmsearch: command not found` | `conda install -c bioconda hmmer` or `brew install hmmer` |
| `Run directory already exists` | Use `make all RUN_ID=xxx FORCE=1` to overwrite |
| Results don't match | Ensure Python 3.12+, check `requirements.txt` versions |
| Permission denied on symlinks | Run from local filesystem, not network drive |

### Data Availability

All data and pre-computed embeddings are available on Zenodo:

**DOI**: [10.5281/zenodo.17370925](https://doi.org/10.5281/zenodo.17370925)

Contents:
- Raw UniProt data (6,465 cleaned sequences)
- Processed domain sequences (1,959 domains)
- Pre-computed ESM-2 embeddings (4 configurations, ~7GB)
- Train/test splits at 40/50/70% identity thresholds
- Complete results registries

To use Zenodo data instead of regenerating:
```bash
# Download and extract
wget https://zenodo.org/records/17370925/files/kinase_data_v1.zip
unzip kinase_data_v1.zip

# Copy embeddings to project
cp -r kinase_data_v1/embeddings/* embeddings/

# Run pipeline (will use existing embeddings)
make all RUN_ID=from_zenodo
```

---

## Citation

```bibtex
@article{kinase_layer_selection_2025,
  title={Layer Selection in Protein Language Models Improves 
         Kinase Functional Classification},
  author={[Authors]},
  journal={[Journal]},
  year={2025},
  doi={10.5281/zenodo.17370925}
}
```

---

## License

MIT License - see [LICENSE](LICENSE) file.

---

## Acknowledgments

- **ESM-2**: Meta AI Research ([fair-esm](https://github.com/facebookresearch/esm))
- **Data**: UniProt, Pfam (InterPro)
- **Tools**: HMMER, CD-HIT, scikit-learn, PyTorch

---

## Contact

- **Issues**: https://github.com/jhaaj08/Kinases-Clustering/issues
- **Email**: See git commit history

---

## Additional Documentation

| Document | Description |
|----------|-------------|
| [MANUSCRIPT.md](MANUSCRIPT.md) | Full manuscript text |
| [docs/Simple_English.md](docs/Simple_English.md) | Plain-language explanation |
| [DEPLOYMENT_GUIDE.md](DEPLOYMENT_GUIDE.md) | Deployment instructions |
| [figures_output/FIGURES_PROCESS.md](figures_output/FIGURES_PROCESS.md) | Figure generation details |
