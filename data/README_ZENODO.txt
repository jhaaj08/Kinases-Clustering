================================================================================
KINASE CLASSIFICATION DATA PACKAGE
Version 1.0
================================================================================

Title: ESM-2 Layer Selection for Kinase Functional Classification
Authors: [Your Name Here]
Date: December 2025
License: CC-BY 4.0

================================================================================
DESCRIPTION
================================================================================

This data package contains the complete dataset for reproducing the kinase
functional classification study using ESM-2 protein language model embeddings.

Key Finding: Different transformer layers are optimal for different tasks:
  - Clustering: Layers 20-30 achieve +137.5% improvement over Layer 33
  - Classification: Layer 33 achieves +3.4% improvement over Layers 20-30

================================================================================
DATA PROVENANCE
================================================================================

Source Database: UniProt SwissProt (reviewed entries only)
UniProt Release: 2025_04
Access Date: September 30, 2025
Query: reviewed:true AND (keyword:KW-0418 OR name:kinase*)
Raw Sequences: 20,262

================================================================================
FILE STRUCTURE
================================================================================

data/
├── raw/                          # Original downloaded data
│   ├── kinases_all.csv           # Full UniProt download (20,262 sequences)
│   ├── kinases_revised.csv       # After CD-HIT 60% clustering (6,465 sequences)
│   ├── uniprot_query.txt         # Exact query used
│   ├── uniprot_release.txt       # Database version and date
│   ├── MANIFEST.txt              # File checksums
│   └── MANIFEST.json             # Machine-readable checksums
│
├── processed/                    # Cleaned and labeled data
│   ├── labels.csv                # Label assignments for all sequences
│   ├── label_policy.json         # Labeling rules
│   ├── label_counts.tsv          # Class distribution
│   └── dataset_manifest_report.json  # Dataset statistics (SOURCE OF TRUTH)
│
├── domains/                      # Extracted kinase domains
│   ├── domains_E001.fasta        # Domain sequences (E-value < 0.01)
│   ├── domain_coords_E001.tsv    # Domain coordinates (1,959 domains)
│   ├── hmmer_domtblout_E001.txt  # Raw HMMER output
│   └── domain_extraction_report.json
│
├── manifests/                    # Dataset definitions
│   ├── whole_seq_excl_other.txt  # Full-length sequences (2,911)
│   ├── domain_E001.txt           # Domain sequences excl. Other (1,392)
│   └── supervised_eligible.txt   # For classification (1,367)
│
├── splits/                       # Homology-aware train/test splits
│   ├── split40_train.txt         # 40% identity (strictest)
│   ├── split40_test.txt
│   ├── split50_train.txt
│   ├── split50_test.txt
│   ├── split70_train.txt
│   ├── split70_test.txt
│   └── splits_report.json        # Full split statistics
│
├── hmm_profiles/                 # Pfam HMM profiles
│   ├── PF00069.hmm               # Pkinase domain
│   └── PF07714.hmm               # Pkinase_Tyr domain
│
results/                          # RESULTS REGISTRIES (SOURCE OF TRUTH)
├── manuscript_numbers.json       # All numbers in the manuscript
├── tables/                       # Generated tables
│   ├── Table1.csv                # Dataset construction
│   ├── TableS1.csv               # Layer ablation
│   └── TableS2.csv               # Baselines comparison
├── clustering/
│   └── clustering_registry.json  # Clustering metrics
├── supervised/
│   ├── lr_split40_metrics.json   # Classification metrics (40%)
│   ├── lr_split50_metrics.json   # Classification metrics (50%)
│   └── lr_split70_metrics.json   # Classification metrics (70%)
├── calibration/
│   └── split40_calibration.json  # Calibration results
├── baselines/
│   └── baselines_split40.csv     # All baseline comparisons
├── retrieval/
│   ├── split40_retrieval.json    # Retrieval experiment
│   └── summary.csv               # Retrieval metrics
└── layer_comparison/
    ├── layer_comparison_results.json  # Layer 33 vs 20-30
    └── layer_comparison_summary.csv

embeddings/esm2_t33_650M/
├── embedding_metadata.json       # Model config and hashes
└── ids.txt                       # Sequence IDs in row order

================================================================================
KEY DATASETS
================================================================================

| Dataset | N | Classes | Description |
|---------|---|---------|-------------|
| Whole-seq (excl. Other) | 2,911 | 10 | Full-length sequences |
| Domain E<0.01 (all) | 1,959 | 11 | All extracted domains |
| Domain E<0.01 (excl. Other) | 1,392 | 10 | For clustering |
| Supervised-eligible | 1,367 | 8 | For classification |

Excluded from classification:
  - Histidine kinases (7 samples): Different catalytic mechanism
  - RGC kinases (18 samples): Not true kinases

================================================================================
RESULTS REGISTRIES (SOURCE OF TRUTH)
================================================================================

All numbers in the manuscript are derived from the results/ directory.
The rule: If a number is not in results/manuscript_numbers.json, it cannot
appear in the manuscript.

Key registries:
  - results/manuscript_numbers.json - Master source of all numbers
  - results/clustering/clustering_registry.json - ARI, NMI, Hungarian accuracy
  - results/supervised/lr_split40_metrics.json - Classification metrics
  - results/baselines/baselines_split40.csv - All baseline comparisons

This prevents "hand-typed drift" where numbers in the text don't match
the actual experimental results.

================================================================================
DATA INTEGRITY NOTES
================================================================================

Known discrepancy (documented for transparency):
  - embeddings/esm2_t33_650M/ids.txt: 1,968 sequences
  - data/domains/domain_coords_E001.tsv: 1,959 domains
  - Difference: 9 sequences are in embeddings but not in domain_coords

This occurred because embeddings were generated from an earlier domain
extraction run. The current domain coords (1,959) is the authoritative
source. Downstream analyses filter to 1,392 (clustering) or 1,367
(classification) anyway, so this discrepancy does not affect results.

================================================================================
REPRODUCIBILITY
================================================================================

All files have SHA-256 checksums in data/raw/MANIFEST.txt. To verify:

  python scripts/generate_raw_manifest.py

This will regenerate hashes and compare against the archived values.

================================================================================
CITATION
================================================================================

If you use this data, please cite:

[Citation will be added after publication]

DOI: [Will be assigned by Zenodo]

================================================================================
CONTACT
================================================================================

For questions about this dataset, please contact:
[Your email here]

================================================================================

