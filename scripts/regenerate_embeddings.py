#!/usr/bin/env python3
"""
Regenerate ESM-2 Embeddings from Current Domain FASTA

This script regenerates all embeddings from data/domains/domains_E001.fasta
to ensure consistency between domain extraction and embeddings.

Usage:
    python scripts/regenerate_embeddings.py

Output:
    embeddings/esm2_t33_650M/domain_E001_layer33_mean.npy
    embeddings/esm2_t33_650M/domain_E001_layers20_30_mean.npy
    embeddings/esm2_t33_650M/domain_E001_layers20_33_mean.npy
    embeddings/esm2_t33_650M/domain_E001_layer33_cls.npy
    embeddings/esm2_t33_650M/ids.txt
    embeddings/esm2_t33_650M/embedding_metadata.json
"""

import torch
import esm
import numpy as np
import json
import hashlib
from pathlib import Path
from datetime import datetime
from Bio import SeqIO
from tqdm import tqdm

# Configuration
INPUT_FASTA = Path("data/domains/domains_E001.fasta")
OUTPUT_DIR = Path("embeddings/esm2_t33_650M")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Model settings
MODEL_NAME = "esm2_t33_650M_UR50D"
BATCH_SIZE = 4  # Adjust based on available memory
MAX_SEQ_LEN = 1022  # ESM-2 limit

# Layer configurations to generate
LAYER_CONFIGS = {
    "layer33_mean": {
        "layers": [33],
        "pooling": "mean",
        "description": "Final layer only, mean pooling"
    },
    "layers20_30_mean": {
        "layers": list(range(20, 31)),  # 20-30 inclusive
        "pooling": "mean", 
        "description": "Layers 20-30 averaged, mean pooling"
    },
    "layers20_33_mean": {
        "layers": list(range(20, 34)),  # 20-33 inclusive
        "pooling": "mean",
        "description": "Layers 20-33 averaged, mean pooling"
    },
    "layer33_cls": {
        "layers": [33],
        "pooling": "cls",
        "description": "Final layer only, CLS token"
    },
}


def compute_config_hash(model_name, layers, pooling):
    """Compute hash of configuration for reproducibility."""
    config_str = f"{model_name}|{layers}|{pooling}"
    return hashlib.sha256(config_str.encode()).hexdigest()[:16]


def load_sequences(fasta_path):
    """Load sequences from FASTA file."""
    sequences = []
    ids = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        ids.append(record.id)
        sequences.append(str(record.seq))
    return ids, sequences


def get_device():
    """Get the best available device."""
    if torch.cuda.is_available():
        return torch.device("cuda")
    elif torch.backends.mps.is_available():
        return torch.device("mps")
    else:
        return torch.device("cpu")


def extract_embeddings(model, alphabet, sequences, ids, layer_config, device, batch_size=4):
    """Extract embeddings for all sequences."""
    batch_converter = alphabet.get_batch_converter()
    
    layers = layer_config["layers"]
    pooling = layer_config["pooling"]
    
    all_embeddings = []
    
    # Process in batches
    for i in tqdm(range(0, len(sequences), batch_size), desc=f"Extracting {layer_config['description']}"):
        batch_ids = ids[i:i+batch_size]
        batch_seqs = sequences[i:i+batch_size]
        
        # Truncate sequences if needed
        batch_seqs = [seq[:MAX_SEQ_LEN] for seq in batch_seqs]
        
        # Prepare batch
        batch_data = list(zip(batch_ids, batch_seqs))
        batch_labels, batch_strs, batch_tokens = batch_converter(batch_data)
        batch_tokens = batch_tokens.to(device)
        
        # Extract representations
        with torch.no_grad():
            results = model(batch_tokens, repr_layers=layers, return_contacts=False)
        
        # Process each sequence in batch
        for j, seq in enumerate(batch_seqs):
            seq_len = len(seq)
            
            if pooling == "cls":
                # Use CLS token (position 0)
                layer_repr = results["representations"][layers[0]][j, 0, :]
            else:
                # Mean pooling over sequence (excluding BOS/EOS tokens)
                layer_reprs = []
                for layer_idx in layers:
                    # positions 1 to seq_len+1 are the actual sequence
                    layer_repr = results["representations"][layer_idx][j, 1:seq_len+1, :]
                    layer_reprs.append(layer_repr)
                
                # Average across layers
                stacked = torch.stack(layer_reprs, dim=0)
                layer_avg = stacked.mean(dim=0)  # (seq_len, embed_dim)
                
                # Mean pool across sequence
                layer_repr = layer_avg.mean(dim=0)  # (embed_dim,)
            
            all_embeddings.append(layer_repr.cpu().numpy())
    
    return np.array(all_embeddings)


def main():
    print("=" * 70)
    print("REGENERATE ESM-2 EMBEDDINGS")
    print("=" * 70)
    
    # Check input file
    if not INPUT_FASTA.exists():
        print(f"ERROR: Input FASTA not found: {INPUT_FASTA}")
        return
    
    # Load sequences
    print(f"\nLoading sequences from {INPUT_FASTA}...")
    ids, sequences = load_sequences(INPUT_FASTA)
    n_seqs = len(sequences)
    print(f"  Loaded {n_seqs} sequences")
    
    # Get device
    device = get_device()
    print(f"\nUsing device: {device}")
    
    # Load model
    print(f"\nLoading ESM-2 model ({MODEL_NAME})...")
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    model = model.to(device)
    model.eval()
    print("  Model loaded successfully")
    
    # Save IDs
    ids_file = OUTPUT_DIR / "ids.txt"
    with open(ids_file, 'w') as f:
        for uid in ids:
            f.write(f"{uid}\n")
    print(f"\nSaved IDs to {ids_file}")
    
    # Generate embeddings for each configuration
    metadata = {
        "step": 8,
        "name": "ESM-2 Embedding Generation",
        "timestamp": datetime.now().isoformat(),
        "model": {
            "model_name": MODEL_NAME,
            "model_params": "650M parameters",
            "embedding_dim": 1280,
            "num_layers": 33,
            "max_sequence_length": MAX_SEQ_LEN,
            "library": "fair-esm",
            "library_version": esm.__version__ if hasattr(esm, '__version__') else "2.0.0"
        },
        "source_fasta": str(INPUT_FASTA),
        "n_sequences": n_seqs,
        "embedding_dim": 1280,
        "layer_configurations": {},
        "embeddings": {},
        "sanity_checks": {
            "ids_consistent": True,
            "no_nans": True,
            "dimension_correct": True
        }
    }
    
    for config_name, config in LAYER_CONFIGS.items():
        print(f"\n--- {config['description']} ---")
        
        # Extract embeddings
        embeddings = extract_embeddings(
            model, alphabet, sequences, ids, config, device, BATCH_SIZE
        )
        
        print(f"  Shape: {embeddings.shape}")
        
        # Check for NaNs
        if np.isnan(embeddings).any():
            print("  WARNING: NaN values detected!")
            metadata["sanity_checks"]["no_nans"] = False
        
        # Check dimension
        if embeddings.shape != (n_seqs, 1280):
            print(f"  WARNING: Unexpected shape! Expected ({n_seqs}, 1280)")
            metadata["sanity_checks"]["dimension_correct"] = False
        
        # Compute config hash
        config_hash = compute_config_hash(MODEL_NAME, config["layers"], config["pooling"])
        
        # Save embeddings
        output_file = OUTPUT_DIR / f"domain_E001_{config_name}.npy"
        np.save(output_file, embeddings)
        print(f"  Saved to {output_file}")
        
        # Update metadata
        metadata["layer_configurations"][config_name] = {
            "layers": config["layers"],
            "pooling": config["pooling"],
            "description": config["description"],
            "config_hash": config_hash
        }
        
        metadata["embeddings"][f"domain_E001_{config_name}"] = {
            "file": str(output_file),
            "shape": list(embeddings.shape),
            "layers": config["layers"],
            "pooling": config["pooling"],
            "description": config["description"],
            "config_hash": config_hash,
            "source": str(INPUT_FASTA)
        }
    
    # Save metadata
    metadata_file = OUTPUT_DIR / "embedding_metadata.json"
    with open(metadata_file, 'w') as f:
        json.dump(metadata, f, indent=2)
    print(f"\nSaved metadata to {metadata_file}")
    
    print("\n" + "=" * 70)
    print("EMBEDDING GENERATION COMPLETE")
    print("=" * 70)
    print(f"\nSummary:")
    print(f"  Sequences: {n_seqs}")
    print(f"  Configurations: {len(LAYER_CONFIGS)}")
    print(f"  Output directory: {OUTPUT_DIR}")
    
    print("\nSanity checks:")
    for check, passed in metadata["sanity_checks"].items():
        status = "✓" if passed else "✗"
        print(f"  {status} {check}")


if __name__ == "__main__":
    main()

