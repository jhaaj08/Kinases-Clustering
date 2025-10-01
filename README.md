# Kinases Clustering Project

This project downloads ALL kinase sequences from UniProt and saves them as a CSV file for analysis and clustering.

## Setup

### Python Dependencies

Install the required dependencies:

```bash
pip install -r requirements.txt
```

### CD-HIT Installation (Required for data cleaning)

CD-HIT is used for sequence clustering at 60% identity to reduce redundancy.

**Install CD-HIT:**

```bash
# macOS (via conda/bioconda)
conda install -c bioconda cd-hit

# Ubuntu/Debian
sudo apt-get install cd-hit

# From source
# https://github.com/weizhongli/cdhit
```

**Verify installation:**
```bash
cd-hit -h
```

## Usage

### Basic Usage

Run the script to download all kinase sequences:

```bash
python download_kinases.py
```

This will create a file `kinases_all.csv` with **all kinase sequences** from UniProt (reviewed/SwissProt entries only).

### CSV Format

The output CSV file contains seven columns:
1. **uniprot_id**: UniProt accession number (e.g., P24941)
2. **protein_name**: Protein name/description (e.g., "Cyclin-dependent kinase 2")
3. **function**: UniProt function annotation text (detailed biological function)
4. **kinome_group_subfamily**: Kinase group/subfamily classification (e.g., "CMGC", "TK", "AGC")
5. **conformation_DFG_aC**: DFG/αC helix conformation states (placeholder - requires KLIFS/PDB data)
6. **inhibitor_class_sensitivity**: Inhibitor class sensitivity (placeholder - requires ChEMBL/experimental data)
7. **sequence**: Full amino acid sequence

Example:
```csv
uniprot_id,protein_name,function,kinome_group_subfamily,conformation_DFG_aC,inhibitor_class_sensitivity,sequence
P24941,Cyclin-dependent kinase 2,FUNCTION: Serine/threonine...,CMGC,Not available,Not available,MENFQKVEK...
```

**Notes:**
- ~82.6% of kinases have function annotations; ~17.4% have "N/A"
- ~98.5% of kinases are classified into kinome groups; ~1.5% are "Unclassified"
- **Conformation and inhibitor sensitivity columns are placeholders**
  - See `EXTERNAL_DATA_INTEGRATION.md` for how to populate these from specialized databases

### Custom Usage

You can use the `download_all_kinases_from_uniprot()` function programmatically:

```python
from download_kinases import download_all_kinases_from_uniprot

# Download all kinases from all organisms
result = download_all_kinases_from_uniprot(
    output_file="kinases_all.csv"
)

# Download only human kinases
result = download_all_kinases_from_uniprot(
    output_file="human_kinases.csv",
    organism="Homo sapiens"
)

# Download mouse kinases
result = download_all_kinases_from_uniprot(
    output_file="mouse_kinases.csv",
    organism="Mus musculus"
)

# Check the results
print(f"Downloaded {result['sequences_downloaded']} sequences")
print(f"Saved to {result['output_file']}")
```

### Function Parameters

- `output_file` (str, optional): Output CSV filename (default: "kinases.csv")
- `organism` (str, optional): Filter by organism name (e.g., "Homo sapiens", "Mus musculus")
  - If `None` (default), downloads kinases from all organisms

### Return Value

The function returns a dictionary containing:
- `sequences_downloaded`: Number of sequences actually downloaded
- `output_file`: Path to the output CSV file
- `query`: UniProt query string used

## Data Statistics

As of the last run:
- **Total kinases downloaded**: 20,262 sequences
- **File size**: ~11 MB
- **Source**: UniProt SwissProt (reviewed entries only)
- **Format**: CSV with 7 columns

**Column Coverage:**
- **Function annotations**: 16,736 kinases (82.6%) have detailed function text
- **Kinome group classification**: 19,966 kinases (98.5%) classified into groups
- **Conformation (DFG/αC)**: Placeholder - requires KLIFS/PDB integration
- **Inhibitor sensitivity**: Placeholder - requires ChEMBL/experimental data integration

**Major Kinase Groups Distribution:**
- Protein kinase superfamily: ~9.2%
- Tyrosine kinases (TK): ~5.7%
- GHMP kinases: ~5.6%
- CAMK group: ~3.5%
- CMGC group: ~3.4%
- Other groups: ~75.6%

See `EXTERNAL_DATA_INTEGRATION.md` for detailed instructions on populating placeholder columns.

## Understanding the Columns

### Kinome Group/Subfamily

This column classifies kinases according to the Manning kinome classification and protein family databases. Major groups include:

- **AGC**: PKA, PKG, PKC families
- **CAMK**: Calcium/calmodulin-dependent kinases
- **CK1**: Casein kinase 1 family
- **CMGC**: CDK, MAPK, GSK3, CLK families
- **STE**: Homologs of yeast Sterile kinases (MAP kinase cascade)
- **TK**: Tyrosine kinases
- **TKL**: Tyrosine kinase-like
- **RGC**: Receptor guanylate cyclases
- **Atypical**: PI3K, mTOR, etc.
- **Histidine**: Histidine kinases (mainly bacterial)

### Conformation (DFG/αC states) - Placeholder

Kinase conformational states are critical for understanding:
- **DFG-in vs DFG-out**: Active vs inactive conformations
- **αC-in vs αC-out**: Helix position affecting ATP binding
- **Drug binding modes**: Type I (DFG-in) vs Type II (DFG-out) inhibitors

To populate: Use KLIFS database (https://klifs.net/) or PDB structural data.

### Inhibitor Class Sensitivity - Placeholder

Classification of kinase sensitivity to different inhibitor types:
- **Type I**: ATP-competitive, bind to active conformation
- **Type II**: Bind to inactive (DFG-out) conformation
- **Type III**: Allosteric, non-ATP competitive
- **Type IV**: Covalent inhibitors
- **Type V**: Bivalent inhibitors

To populate: Use ChEMBL (bioactivity data) or KIDFamMap (inhibitor classification).

## Reading the CSV File

### Using Python pandas:

```python
import pandas as pd

# Read the CSV file
df = pd.read_csv('kinases_all.csv')

# Display basic info
print(f"Total kinases: {len(df)}")
print(df.head())

# Access specific columns
uniprot_ids = df['uniprot_id'].tolist()
protein_names = df['protein_name'].tolist()
functions = df['function'].tolist()
sequences = df['sequence'].tolist()

# Filter kinases by protein name
cdk_kinases = df[df['protein_name'].str.contains('CDK', case=False)]
print(f"CDK kinases: {len(cdk_kinases)}")

# Filter kinases that have function annotations
kinases_with_function = df[df['function'].notna() & (df['function'] != 'N/A')]
print(f"Kinases with function annotations: {len(kinases_with_function)}")
```

### Using Python csv module:

```python
import csv

with open('kinases_all.csv', 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        uniprot_id = row['uniprot_id']
        protein_name = row['protein_name']
        function = row['function']
        sequence = row['sequence']
        print(f"{uniprot_id}: {protein_name}")
        if function != 'N/A':
            print(f"  Function: {function[:100]}...")
```

## Examples

### Example 1: Download Organism-Specific Kinases

```python
from download_kinases import download_all_kinases_from_uniprot

# Human kinases
download_all_kinases_from_uniprot("human_kinases.csv", "Homo sapiens")

# Mouse kinases
download_all_kinases_from_uniprot("mouse_kinases.csv", "Mus musculus")

# E. coli kinases
download_all_kinases_from_uniprot("ecoli_kinases.csv", "Escherichia coli")
```

### Example 2: Analyze Sequence Lengths

```python
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('kinases_all.csv')
df['seq_length'] = df['sequence'].apply(len)

# Plot distribution
df['seq_length'].hist(bins=50)
plt.xlabel('Sequence Length (amino acids)')
plt.ylabel('Frequency')
plt.title('Distribution of Kinase Sequence Lengths')
plt.show()

print(f"Average length: {df['seq_length'].mean():.0f} amino acids")
print(f"Median length: {df['seq_length'].median():.0f} amino acids")
```

## Notes

- The function uses the UniProt REST API (https://rest.uniprot.org)
- Only reviewed (SwissProt) entries are downloaded for high quality
- The query searches for proteins with "kinase" in their name/description
- Network connection required for downloading
- Large downloads may take several minutes
- Be considerate of UniProt's servers when making frequent requests

## Workflow

### Step 1: Download All Kinases

```bash
python download_kinases.py
```

This creates `kinases_all.csv` with 20,262 kinases and 7 columns.

### Step 2: Clean and Reduce Label Cardinality

```bash
python data_clean.py
```

This creates `kinases_revised.csv` with:
- 6,465 representative kinases (CD-HIT clustered at 60% identity)
- 11 major kinase groups (AGC, CAMK, CK1, CMGC, STE, TK, TKL, RGC, Atypical, Histidine, Other)
- 9 columns including sequence, function, and hierarchical classification
- Ready for embedding generation and clustering

**Cleaning pipeline:**
1. Remove exact duplicates: 20,262 → 17,391 (-14.2%)
2. CD-HIT clustering @ 60%: 17,391 → 6,465 (-62.8%)
3. Total reduction: -68.1% (13,797 sequences removed)

**Use `kinases_revised.csv` for clustering and ML tasks.**

## Label Hierarchy in kinases_revised.csv

The cleaned dataset provides 2 main levels of classification:

| Level | Column | Description | Use Case |
|-------|--------|-------------|----------|
| 1 | `kinome_group_major` | 11 major groups | Clustering, high-level classification |
| 2 | `kinome_group_subfamily` | Original subfamily labels | Detailed annotation, literature reference |

### Major Groups (11 classes)
- **Protein Kinases**: AGC, CAMK, CK1, CMGC, STE, TK, TKL, RGC, Atypical
- **Other Kinases**: Histidine, Other (non-eukaryotic & metabolic kinases)

**Distribution in cleaned dataset (6,465 kinases)**:
- Other: 4,536 (70.16%) - heterogeneous group
- TK: 601 (9.30%)
- CMGC: 289 (4.47%)
- CAMK: 287 (4.44%)
- AGC: 185 (2.86%)
- Histidine: 160 (2.47%)
- Atypical: 154 (2.38%)
- STE: 138 (2.13%)
- TKL: 63 (0.97%)
- CK1: 50 (0.77%)
- RGC: 2 (0.03%)

### Cleaning Statistics
- Original: 20,262 kinases
- After duplicate removal: 17,391 kinases (-14.2%)
- After CD-HIT @ 60%: 6,465 kinases (-68.1% total)
- Total removed: 13,797 sequences
  - Exact duplicates: 2,871
  - CD-HIT clustering: 10,926

### Step 3: Generate ESM-2 Embeddings

```bash
python generate_esm2_embeddings.py --input kinases_revised.csv --output-dir kinases_embeddings
```

This generates sequence embeddings using the ESM-2 protein language model (650M parameters):
- **Input**: `kinases_revised.csv` (6,465 kinases)
- **Model**: ESM-2 (esm2_t33_650M_UR50D)
- **Method**: Sliding window (max_len=1022, stride=900) with length-weighted mean pooling
- **Output**: `kinases_embeddings/` folder with:
  - `esm2_embeddings.npy`: (6465, 1280) array - 31.6 MB
  - `esm2_index.csv`: UniProt IDs corresponding to rows
  - `kinases_with_embeddings.csv`: Combined CSV - 95.3 MB

**Processing time**: ~2 hours on CPU (M-series Mac or similar)

**Command-line options**:
```bash
python generate_esm2_embeddings.py \
  --input kinases_revised.csv \
  --output-dir kinases_embeddings \
  --device auto \              # 'auto', 'cpu', or 'cuda'
  --max_len 1022 \            # ESM-2 max tokens per window
  --stride 900                # Sliding window stride
```

### Step 3b: Extract Kinase Domains (HMMER)

Use HMMER with Pfam PF00069 to extract catalytic domains:

```bash
python extract_kinase_domains.py
```

Output: `kinases_domains.csv` with domain-only sequences and boundaries.

### Step 3c: Generate ESM-2 Embeddings (Domains)

```bash
python generate_esm2_embeddings.py --input kinases_domains.csv --output-dir kinases_domains_embeddings
```

This produces domain embeddings `(N, 1280)` in `kinases_domains_embeddings/`.

### Step 4: Motif Features and Clustering

1) Extract motif features from domains:

```bash
python extract_motif_features.py
```

Output: `kinases_domains_with_motifs.csv` with 22 motif features.

2) Fuse domain embeddings + motif features and run K-Means (k=10, excluding "Other"):

```bash
python cluster_with_motifs.py
```

Key results (see `clustering/RESULTS_SUMMARY.md` for full table):
- Domain-only (1280-d): ARI 0.2678, NMI 0.3601, Purity 0.6243, Hungarian 0.4505, Best cluster 87.0%
- Domain + Motifs (1302-d): ARI 0.2741, NMI 0.3658, Purity 0.6251, Hungarian 0.4578, Best cluster 88.2%

Files generated:
- `clustering/kmeans10_domains_assignments.csv`, `clustering/kmeans10_domains_report.txt`
- `clustering/kmeans10_domain_motifs_assignments.csv`, `clustering/kmeans10_domain_motifs_report.txt`
- `clustering/RESULTS_SUMMARY.md`

## Embedding Features

### ESM-2 Embeddings

The ESM-2 model captures:
- **Structural patterns**: Catalytic domain features, secondary structure
- **Functional motifs**: DFG, HRD, APE motifs, P-loop patterns
- **Evolutionary signals**: Sequence conservation, family-specific patterns
- **Substrate specificity**: Tyrosine vs Serine/Threonine kinases separate well

**Key findings**:
- ✅ TK and Histidine kinases show strong separation (>75% purity)
- ✅ Embeddings capture functional similarity beyond phylogeny
- ⚠️ CMGC/CAMK/AGC groups mix (shared Ser/Thr mechanism)

## File Structure

```
Kinases-Clustering/
├── download_kinases.py                 # Download kinases from UniProt
├── data_clean.py                       # Clean data and reduce label cardinality
├── verify_data.py                      # Verify dataset statistics
├── extract_kinase_domains.py           # HMMER-based domain extraction ⭐
├── extract_motif_features.py           # Motif feature extraction ⭐
├── cluster_with_motifs.py              # Fuse features + clustering ⭐
├── generate_esm2_embeddings.py         # ESM-2 embeddings (whole or domain) ⭐
├── requirements.txt                    # Python dependencies
├── README.md                           # This file
├── EXTERNAL_DATA_INTEGRATION.md        # Guide for adding external data
├── kinases_all.csv                     # Raw download (20,262 kinases)
├── kinases_revised.csv                 # Cleaned dataset (6,465 kinases) ⭐
├── kinases_domains.csv                 # Domain-only sequences (HMMER)
├── duplicates_report.csv               # Duplicate sequences report
├── kinases_embeddings/                 # Whole-seq embeddings (optional)
│   ├── esm2_embeddings.npy             # (6465, 1280)
│   ├── esm2_index.csv                  # UniProt ID index
│   └── kinases_with_embeddings.csv     # Combined CSV
├── kinases_domains_embeddings/         # Domain embeddings ⭐
│   ├── esm2_embeddings.npy             # (N, 1280)
│   ├── esm2_index.csv                  # UniProt ID index
│   └── kinases_with_embeddings.csv     # Combined CSV
└── clustering/                         # Clustering analysis results ⭐
    ├── kmeans10_no_other_assignments.csv
    ├── kmeans10_no_other_report.txt
    ├── kmeans10_domains_assignments.csv
    ├── kmeans10_domains_report.txt
    ├── kmeans10_domain_motifs_assignments.csv
    ├── kmeans10_domain_motifs_report.txt
    └── RESULTS_SUMMARY.md
```