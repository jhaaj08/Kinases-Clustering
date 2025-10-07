# Determinism and Reproducibility

## Overview

This codebase implements comprehensive determinism controls to ensure bit-exact reproducibility across runs and platforms (where possible).

## Deterministic Components

### 1. Random Seeds (All Fixed to 42)

**Python built-in random**:
```python
import random
random.seed(42)
```

**NumPy**:
```python
import numpy as np
np.random.seed(42)
```

**PyTorch**:
```python
import torch
torch.manual_seed(42)
torch.cuda.manual_seed_all(42)  # If GPU used
```

**scikit-learn**:
```python
# All estimators use random_state=42
KMeans(random_state=42)
LogisticRegression(random_state=42)
train_test_split(random_state=42)
```

### 2. Deterministic Algorithms

**PyTorch** (when deterministic=True in configs):
```python
torch.use_deterministic_algorithms(True)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False
```

**K-means**:
```python
KMeans(algorithm='lloyd', random_state=42, n_init=50)
# lloyd algorithm is deterministic with fixed seed
```

### 3. Data Processing Order

**Sorted iteration** (ensures consistent order):
```python
# GOOD: Sorted keys
for uniprot_id in sorted(df['uniprot_id'].unique()):
    process(uniprot_id)

# BAD: Unsorted (order may vary)
for uniprot_id in df['uniprot_id'].unique():
    process(uniprot_id)
```

**Deterministic file reading**:
```python
# Files read in sorted order
files = sorted(Path('dir').glob('*.csv'))
```

### 4. Floating Point Precision

**Controlled precision** (default: fp32):
```python
# ESM-2 embeddings
model.eval()
with torch.no_grad():
    results = model(tokens, repr_layers=layers)
# Default: fp32 on CPU (deterministic)
# GPU: fp16/bf16 may have minor variations (<0.1%)
```

**Reproducible operations**:
```python
# Use stable algorithms
np.mean(embeddings, axis=0)  # Reproducible
embeddings.mean(dim=0)  # PyTorch, reproducible with deterministic mode
```

## Known Sources of Non-Determinism

### 1. External Tools

**CD-HIT** (bioconda):
- Deterministic with fixed parameters
- Output order may vary across versions
- **Mitigation**: We use v4.8.1 (pinned)

**HMMER** (domain search):
- Deterministic for single-threaded (`--cpu 1`)
- Multi-threaded may have minor score variations
- **Mitigation**: We use single-threaded mode

### 2. Platform Differences

**Hardware**:
- CPU (Apple M-series): Fully deterministic
- GPU (NVIDIA): Minor fp16/bf16 variations possible (<0.1% metric difference)
- **Mitigation**: Use fp32 for critical runs, document hardware

**Operating System**:
- File system iteration order can vary
- **Mitigation**: Always sort files/keys before iteration

**Python version**:
- Minor differences in random number generation across Python 3.10/3.11/3.12
- **Mitigation**: We use Python 3.12 (documented in pyproject.toml)

### 3. Parallel Processing

**Multi-threading**:
- scikit-learn uses OpenMP/MKL, may have thread-order non-determinism
- **Mitigation**: Set `n_jobs=1` for critical runs, or accept minor variations

**PyTorch DataLoader**:
- `num_workers > 0` introduces non-determinism
- **Mitigation**: We use `num_workers=0` (single-threaded)

## Reproducibility Checklist

✅ **Seeds set** (Python, NumPy, PyTorch, scikit-learn)  
✅ **Deterministic algorithms** enabled (PyTorch)  
✅ **Sorted iteration** (files, dataframe rows, dictionary keys)  
✅ **Fixed precision** (fp32 default)  
✅ **Pinned versions** (pyproject.toml, environment.yml)  
✅ **Single-threaded** external tools (HMMER, CD-HIT)  
✅ **Documented hardware** (Apple M-series CPU)  
✅ **Saved splits** (data/splits_*.json, exact train/test IDs)  
✅ **Provenance tracking** (data/provenance.json)  
✅ **Configuration files** (configs/config.yaml)  

## Validation

To verify reproducibility:

```bash
# Run pipeline twice
snakemake --cores 1 all --config seed=42
mv results/ results_run1/

snakemake --cores 1 all --config seed=42
mv results/ results_run2/

# Compare outputs (should be identical)
diff results_run1/statistical_comparisons.csv results_run2/statistical_comparisons.csv
```

Expected: **Identical results** on same hardware/OS/Python version.

## Documentation

All experiments log:
- Random seed used
- Hardware (CPU/GPU, model)
- Software versions (Python, PyTorch, scikit-learn, HMMER, CD-HIT)
- Precision (fp32/fp16/bf16)
- Deterministic mode (on/off)

Logs saved to: `logs/*.log`

## References

- PyTorch Reproducibility: https://pytorch.org/docs/stable/notes/randomness.html
- NumPy Random State: https://numpy.org/doc/stable/reference/random/generator.html
- scikit-learn Random State: https://scikit-learn.org/stable/glossary.html#term-random_state

