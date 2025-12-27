# Kinase Classification Web Application

Interactive web interface for kinase family prediction using ESM-2 embeddings.

## Features

- **Sequence Input**: Paste FASTA or raw sequence
- **Automatic Domain Extraction**: HMMER with Pfam PF00069
- **ESM-2 Embeddings**: Mid-layer averaging (layers 20-33) for optimal performance
- **Top-3 Predictions**: Calibrated probabilities with confidence flagging
- **Motif Analysis**: Visualize DFG, HRD, APE, VAIK, P-loop, gatekeeper
- **Exemplar Retrieval**: Find 5 most similar training sequences
- **JSON Export**: Download complete results

## Quick Start

### Prerequisites

1. **HMMER** (for domain extraction):
```bash
# macOS
conda install -c bioconda hmmer

# Ubuntu/Debian
sudo apt-get install hmmer
```

2. **Python packages**:
```bash
pip install -r webapp/requirements.txt
```

3. **Trained models** (from parent directory):
   - `supervised_results/logistic_regression_model.joblib`
   - `kinases_domains_e0.01_layers_mid/esm2_embeddings.npy`
   - `kinases_domains_e0.01_layers_mid/esm2_index.csv`
   - `kinases_domains_e0.01.csv`
   - `PF00069.hmm`

### Run Locally

```bash
cd /path/to/Kinases-Clustering
python webapp/app.py
```

Then open: **http://localhost:7860**

## Usage

### Input Formats

**FASTA**:
```
>CDK2_HUMAN
MENFQKVEKIGEGTYGVVYKARNKLTGEVVALKKIRLDTETEGVPS...
```

**Raw sequence**:
```
MENFQKVEKIGEGTYGVVYKARNKLTGEVVALKKIRLDTETEGVPS...
```

### Example Sequences

The app includes 3 pre-loaded examples:
- **CDK2** (CMGC family)
- **EGFR** (TK family)  
- **CaMK2** (CAMK family)

### Output

1. **Summary**: Top prediction, confidence level, recommendation
2. **Predictions Tab**: Bar chart with top-3 families and probabilities
3. **Motifs Tab**: Identified motifs with color-coded sequence visualization
4. **Similar Sequences Tab**: Nearest training exemplars with UniProt IDs
5. **Export Tab**: Complete results in JSON format

### Confidence Levels

- 🟢 **High Confidence** (≥70%): Prediction can be used directly
- 🟡 **Needs Review** (<70%): Recommend manual review or experimental validation

## Deployment

### Option 1: Local Server (Recommended for Lab Use)

```bash
# Run on specific port
python webapp/app.py --server-port 7860

# Allow external access
python webapp/app.py --server-name 0.0.0.0
```

### Option 2: Docker Container

```bash
# Build image
docker build -t kinase-classifier webapp/

# Run container
docker run -p 7860:7860 kinase-classifier
```

### Option 3: Hugging Face Spaces

1. Create a new Space on https://huggingface.co/spaces
2. Upload files from `webapp/` directory
3. Include trained models (or download in Space)
4. Space will auto-deploy

**Note**: Large model files (>500 MB) should be stored with Git LFS.

## Architecture

```
webapp/
├── app.py                 # Main Gradio interface
├── predictor.py           # Prediction pipeline
├── motif_highlighter.py   # Sequence visualization
├── requirements.txt       # Python dependencies
├── Dockerfile             # Container configuration
└── README.md              # This file
```

### Prediction Pipeline

```
User Input
    ↓
Domain Extraction (HMMER)
    ↓
ESM-2 Embedding (layers 20-33, mean pooling)
    ↓
Classification (calibrated logistic regression)
    ↓
Top-3 Predictions + Confidence
    ↓
Exemplar Retrieval + Motif Analysis
    ↓
Display Results
```

## Performance

### Speed

- **Domain Extraction**: ~2-5 seconds
- **Embedding Generation**: ~5-10 seconds (CPU), ~1-2 seconds (GPU)
- **Classification**: <1 second
- **Total**: ~10-20 seconds per sequence

### Accuracy

- **Test Accuracy**: 74.9% (homology-aware, 40% identity)
- **Macro F1**: 0.668
- **Top-3 Accuracy**: 94.8%

Best performing families:
- CAMK: F1 = 0.928
- Atypical: F1 = 0.815
- CMGC: F1 = 0.792

## Troubleshooting

### "HMMER not available"

Install HMMER:
```bash
conda install -c bioconda hmmer
# or
sudo apt-get install hmmer
```

### "Model files not found"

Run the pipeline first to generate models:
```bash
cd /path/to/Kinases-Clustering
bash scripts/run_all_pipeline.sh
```

### Out of Memory

The app uses CPU by default. If you have a GPU and want to use it, modify `predictor.py`:

```python
self.device = 'cuda' if torch.cuda.is_available() else 'cpu'
```

### Slow Performance

- Use GPU if available (5-10× faster)
- Consider using smaller ESM-2 model (150M parameters)
- Cache embeddings for frequently analyzed sequences

## Security

**For production deployment**:

- Enable authentication (Gradio supports username/password)
- Rate limiting to prevent abuse
- Input validation (max sequence length, allowed characters)
- HTTPS encryption
- Regular security updates

## Citation

If you use this web app, please cite:

```
[Citation will be added after publication]
```

## License

MIT License - See LICENSE file

## Contact

For issues or questions:
- GitHub Issues: https://github.com/jhaaj08/Kinases-Clustering/issues
- See git commit history for contact information

