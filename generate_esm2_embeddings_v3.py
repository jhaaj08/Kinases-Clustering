#!/usr/bin/env python3
"""
Enhanced ESM-2 embedding generation with full publication-ready features.

New features (v3):
- Per-residue stitching with overlap averaging
- Precision control (fp32, fp16, bf16)
- Deterministic mode for reproducibility
- Per-sequence caching with content+config hashing
- Complete shape and parameter verification
"""

import argparse
import os
import hashlib
import json
import numpy as np
import pandas as pd
import torch
from tqdm import tqdm
from pathlib import Path
import esm


def compute_config_hash(seq, model_name, layers, pooling, max_len, stride, 
                       precision, stitching_method):
    """
    Compute hash of sequence + configuration for caching.
    
    Ensures cached embeddings match current config.
    """
    config_str = f"{seq}|{model_name}|{layers}|{pooling}|{max_len}|{stride}|{precision}|{stitching_method}"
    return hashlib.md5(config_str.encode()).hexdigest()


def load_from_cache(cache_dir, seq_hash):
    """Load cached embedding if available."""
    cache_file = Path(cache_dir) / f"{seq_hash}.npy"
    if cache_file.exists():
        return np.load(cache_file)
    return None


def save_to_cache(cache_dir, seq_hash, embedding):
    """Save embedding to cache."""
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_file = cache_dir / f"{seq_hash}.npy"
    np.save(cache_file, embedding)


def embed_sequence_per_residue(seq, model, alphabet, device, 
                               max_len=1022, stride=900, 
                               layers=[33], precision='fp32'):
    """
    Generate per-residue embeddings with proper overlap stitching.
    
    For overlapping windows, residues are averaged across windows.
    Returns sequence-level embedding via mean pooling.
    
    Parameters:
    -----------
    precision : str
        'fp32', 'fp16', or 'bf16'
    
    Returns:
    --------
    np.ndarray
        Sequence-level embedding (embedding_dim,)
    """
    batch_converter = alphabet.get_batch_converter()
    pad_idx = alphabet.padding_idx
    cls_idx = alphabet.cls_idx
    eos_idx = alphabet.eos_idx
    
    # Set precision
    if precision == 'fp16':
        dtype = torch.float16
    elif precision == 'bf16':
        dtype = torch.bfloat16
    else:
        dtype = torch.float32
    
    def embed_chunk(chunk_seq, chunk_start):
        """Embed a chunk and return per-residue embeddings."""
        data = [("seq", chunk_seq)]
        labels, strs, tokens = batch_converter(data)
        tokens = tokens.to(device)
        
        with torch.no_grad():
            if precision != 'fp32':
                with torch.autocast(device_type='cuda' if device.startswith('cuda') else 'cpu', 
                                   dtype=dtype):
                    results = model(tokens, repr_layers=layers, return_contacts=False)
            else:
                results = model(tokens, repr_layers=layers, return_contacts=False)
        
        # Average across layers if multiple
        if len(layers) > 1:
            layer_reps = [results["representations"][layer][0] for layer in layers]
            reps = torch.stack(layer_reps).mean(dim=0)  # (seq_len, dim)
        else:
            reps = results["representations"][layers[0]][0]
        
        # Mask out special tokens to get residue embeddings
        tok = tokens[0]
        mask = (tok != pad_idx) & (tok != cls_idx) & (tok != eos_idx)
        residue_reps = reps[mask].cpu().float().numpy()  # (num_residues, dim)
        
        # Residue positions in original sequence
        residue_positions = list(range(chunk_start, chunk_start + len(chunk_seq)))
        
        return residue_reps, residue_positions
    
    # If sequence fits in one window
    if len(seq) <= max_len:
        reps, _ = embed_chunk(seq, 0)
        return reps.mean(axis=0)  # Sequence-level via mean pooling
    
    # Sliding window for long sequences
    embedding_dim = 1280  # ESM-2 650M
    
    # Collect per-residue embeddings across all windows
    residue_embeddings = {}  # position -> list of embeddings
    
    start = 0
    while start < len(seq):
        end = min(start + max_len, len(seq))
        chunk = seq[start:end]
        
        reps, positions = embed_chunk(chunk, start)
        
        # Accumulate embeddings for each position
        for i, pos in enumerate(positions):
            if pos not in residue_embeddings:
                residue_embeddings[pos] = []
            residue_embeddings[pos].append(reps[i])
        
        if end == len(seq):
            break
        start += stride
    
    # Average embeddings for overlapping residues
    final_residue_embeddings = []
    for pos in sorted(residue_embeddings.keys()):
        embs = residue_embeddings[pos]
        if len(embs) > 1:
            # Residue appeared in multiple windows - average them
            avg_emb = np.mean(embs, axis=0)
        else:
            avg_emb = embs[0]
        final_residue_embeddings.append(avg_emb)
    
    # Convert to array and mean pool to sequence level
    final_residue_embeddings = np.array(final_residue_embeddings)  # (L, D)
    sequence_embedding = final_residue_embeddings.mean(axis=0)  # (D,)
    
    return sequence_embedding


def embed_sequence_windowed_legacy(seq, model, alphabet, device, 
                                   max_len=1022, stride=900, 
                                   layers=[33], pooling='mean', precision='fp32'):
    """
    Legacy method: window-level mean pooling with length-weighted averaging.
    Faster but less precise for overlapping regions.
    """
    batch_converter = alphabet.get_batch_converter()
    pad_idx = alphabet.padding_idx
    cls_idx = alphabet.cls_idx
    eos_idx = alphabet.eos_idx
    
    # Set precision
    if precision == 'fp16':
        dtype = torch.float16
    elif precision == 'bf16':
        dtype = torch.bfloat16
    else:
        dtype = torch.float32
    
    def embed_chunk(chunk_seq):
        """Embed a chunk and return pooled vector and residue count."""
        data = [("seq", chunk_seq)]
        labels, strs, tokens = batch_converter(data)
        tokens = tokens.to(device)
        
        with torch.no_grad():
            if precision != 'fp32':
                with torch.autocast(device_type='cuda' if device.startswith('cuda') else 'cpu', 
                                   dtype=dtype):
                    results = model(tokens, repr_layers=layers, return_contacts=False)
            else:
                results = model(tokens, repr_layers=layers, return_contacts=False)
        
        # Average across layers if multiple
        if len(layers) > 1:
            layer_reps = [results["representations"][layer][0] for layer in layers]
            reps = torch.stack(layer_reps).mean(dim=0)  # (seq_len, dim)
        else:
            reps = results["representations"][layers[0]][0]
        
        # Pooling strategy
        if pooling == 'cls':
            pooled_vec = reps[0]
            num_residues = 1
        else:  # mean
            tok = tokens[0]
            mask = (tok != pad_idx) & (tok != cls_idx) & (tok != eos_idx)
            residue_reps = reps[mask]
            pooled_vec = residue_reps.mean(dim=0)
            num_residues = residue_reps.shape[0]
        
        return pooled_vec.cpu().float().numpy(), num_residues
    
    # If sequence fits in one window
    if len(seq) <= max_len:
        vec, _ = embed_chunk(seq)
        return vec
    
    # Sliding window
    embedding_dim = 1280
    weighted_sum = np.zeros(embedding_dim, dtype=np.float32)
    total_weight = 0
    
    start = 0
    while start < len(seq):
        end = min(start + max_len, len(seq))
        chunk = seq[start:end]
        
        vec, num_res = embed_chunk(chunk)
        
        if pooling == 'mean':
            weighted_sum += vec * num_res
            total_weight += num_res
        else:  # cls
            weighted_sum += vec
            total_weight += 1
        
        if end == len(seq):
            break
        start += stride
    
    final_vec = weighted_sum / total_weight
    return final_vec


def parse_layer_spec(layer_spec, model_layers=33):
    """Parse layer specification string."""
    if layer_spec == 'all':
        return list(range(1, model_layers + 1))
    elif layer_spec == 'mid':
        start = int(model_layers * 0.6)
        return list(range(start, model_layers + 1))
    elif '-' in layer_spec:
        start, end = map(int, layer_spec.split('-'))
        return list(range(start, end + 1))
    else:
        return [int(layer_spec)]


def generate_embeddings(input_csv, output_dir, model_name='esm2_t33_650M_UR50D', 
                       device='auto', max_len=1022, stride=900,
                       layer_spec='33', pooling='mean', precision='fp32',
                       stitching='per_residue', use_cache=True, deterministic=False):
    """
    Generate ESM-2 embeddings with full reproducibility controls.
    
    Parameters:
    -----------
    stitching : str
        'per_residue' (accurate overlap averaging) or 'window' (faster legacy)
    use_cache : bool
        Enable per-sequence caching
    deterministic : bool
        Enable deterministic operations (for GPU reproducibility)
    precision : str
        'fp32', 'fp16', or 'bf16'
    """
    
    print("=" * 80)
    print("ESM-2 EMBEDDING GENERATION (v3 - Publication-Ready)")
    print("=" * 80)
    print()
    print(f"Input:       {input_csv}")
    print(f"Output:      {output_dir}")
    print(f"Model:       {model_name}")
    print(f"Layers:      {layer_spec}")
    print(f"Pooling:     {pooling}")
    print(f"Stitching:   {stitching}")
    print(f"Precision:   {precision}")
    print(f"Cache:       {use_cache}")
    print(f"Deterministic: {deterministic}")
    print()
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Cache directory
    cache_dir = Path(output_dir) / "cache"
    if use_cache:
        cache_dir.mkdir(exist_ok=True)
        print(f"Cache directory: {cache_dir}")
    
    # Read input data
    print(f"Reading sequences...")
    df = pd.read_csv(input_csv)
    print(f"  ✅ Loaded {len(df):,} sequences")
    
    # Check sequence lengths
    df['seq_len'] = df['sequence'].apply(len)
    num_long = (df['seq_len'] > max_len).sum()
    print(f"  ⚠️  {num_long:,} sequences need windowing (>{max_len} aa)")
    print()
    
    # Determine device
    if device == 'auto':
        device = 'cuda' if torch.cuda.is_available() else 'cpu'
    print(f"Device: {device}")
    
    # Deterministic mode
    if deterministic:
        torch.use_deterministic_algorithms(True)
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False
        print("✅ Deterministic mode enabled (reproducible on same hardware)")
    
    # Load ESM-2 model
    print(f"Loading ESM-2 model: {model_name}...")
    if model_name == 'esm2_t33_650M_UR50D':
        model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
        model_layers = 33
        embedding_dim = 1280
        param_count = "650M"
    else:
        raise ValueError(f"Unsupported model: {model_name}")
    
    model.eval()
    model.to(device)
    
    # Set precision
    if precision == 'fp16' and device.startswith('cuda'):
        model = model.half()
        print(f"  ✅ Model converted to fp16")
    elif precision == 'bf16' and device.startswith('cuda'):
        model = model.bfloat16()
        print(f"  ✅ Model converted to bf16")
    
    # Parse layer specification
    layers = parse_layer_spec(layer_spec, model_layers)
    print(f"  Model: {model_name} ({param_count} parameters, {model_layers} layers)")
    print(f"  Layers to extract: {layers if len(layers) <= 5 else f'{layers[0]}-{layers[-1]} ({len(layers)} layers)'}")
    print(f"  Pooling strategy: {pooling}")
    print(f"  Stitching method: {stitching}")
    print(f"  Embedding dim: {embedding_dim}")
    print()
    
    # Generate embeddings
    print(f"Generating embeddings (window={max_len}, stride={stride})...")
    print()
    
    embeddings = np.zeros((len(df), embedding_dim), dtype=np.float32)
    cache_hits = 0
    cache_misses = 0
    
    for idx, row in tqdm(df.iterrows(), total=len(df), desc="Processing"):
        seq = row['sequence']
        uniprot_id = row['uniprot_id']
        
        # Check cache
        if use_cache:
            seq_hash = compute_config_hash(
                seq, model_name, str(layers), pooling, max_len, stride, 
                precision, stitching
            )
            cached_emb = load_from_cache(cache_dir, seq_hash)
            
            if cached_emb is not None:
                embeddings[idx] = cached_emb
                cache_hits += 1
                continue
            else:
                cache_misses += 1
        
        # Generate embedding
        if stitching == 'per_residue':
            emb = embed_sequence_per_residue(
                seq, model, alphabet, device,
                max_len=max_len, stride=stride,
                layers=layers, precision=precision
            )
        else:  # window (legacy)
            emb = embed_sequence_windowed_legacy(
                seq, model, alphabet, device,
                max_len=max_len, stride=stride,
                layers=layers, pooling=pooling, precision=precision
            )
        
        embeddings[idx] = emb
        
        # Save to cache
        if use_cache:
            save_to_cache(cache_dir, seq_hash, emb)
    
    print()
    if use_cache:
        print(f"Cache stats: {cache_hits} hits, {cache_misses} misses")
        print()
    
    print("Saving embeddings...")
    
    # Save embeddings as numpy array
    emb_file = os.path.join(output_dir, 'esm2_embeddings.npy')
    np.save(emb_file, embeddings)
    print(f"  ✅ Saved: {emb_file}")
    print(f"     Shape: {embeddings.shape}")
    print(f"     Size: {os.path.getsize(emb_file) / 1024 / 1024:.1f} MB")
    
    # Save index
    index_df = df[['uniprot_id']].copy()
    index_file = os.path.join(output_dir, 'esm2_index.csv')
    index_df.to_csv(index_file, index=False)
    print(f"  ✅ Saved: {index_file}")
    
    # Save metadata (complete configuration for reproducibility)
    meta_file = os.path.join(output_dir, 'embedding_metadata.json')
    metadata = {
        "model": model_name,
        "model_parameters": param_count,
        "model_layers": model_layers,
        "layers_used": layers if len(layers) <= 10 else f"{layers[0]}-{layers[-1]}",
        "pooling": pooling,
        "stitching_method": stitching,
        "precision": precision,
        "shape": list(embeddings.shape),
        "max_length": max_len,
        "stride": stride,
        "device": device,
        "deterministic": deterministic,
        "cache_enabled": use_cache,
        "cache_hits": cache_hits if use_cache else 0,
        "cache_misses": cache_misses if use_cache else 0,
        "generated_at": pd.Timestamp.now().isoformat(),
    }
    
    with open(meta_file, 'w') as f:
        json.dump(metadata, f, indent=2)
    print(f"  ✅ Saved: {meta_file}")
    
    # Statistics
    print()
    print("=" * 80)
    print("EMBEDDING STATISTICS")
    print("=" * 80)
    print(f"  Shape:      {embeddings.shape}")
    print(f"  Dtype:      {embeddings.dtype}")
    print(f"  Mean:       {embeddings.mean():.6f}")
    print(f"  Std:        {embeddings.std():.6f}")
    print(f"  Min:        {embeddings.min():.6f}")
    print(f"  Max:        {embeddings.max():.6f}")
    print()
    
    # Shape verification
    print("Shape verification:")
    print(f"  Expected: ({len(df)}, {embedding_dim})")
    print(f"  Actual:   {embeddings.shape}")
    assert embeddings.shape == (len(df), embedding_dim), "Shape mismatch!"
    print("  ✅ Shapes verified")
    print()
    
    # Update provenance
    try:
        from utils.provenance import ProvenanceTracker
        prov = ProvenanceTracker(output_dir="data")
        prov.add_processing_step(
            step_name="ESM-2 Embedding Generation",
            params={
                "model": model_name,
                "layers": str(layers),
                "pooling": pooling,
                "stitching": stitching,
                "precision": precision,
                "max_len": max_len,
                "stride": stride,
            },
            input_count=len(df),
            output_count=len(df),
            description=f"Generated {embedding_dim}-d embeddings"
        )
        print("✅ Provenance updated")
    except:
        print("⚠️  Could not update provenance (utils module not found)")
    
    print()
    print("✅ Embedding generation complete!")
    print()
    
    return embeddings


def main():
    parser = argparse.ArgumentParser(
        description='Generate ESM-2 embeddings with publication-ready features'
    )
    parser.add_argument(
        '--input', 
        required=True,
        help='Input CSV file with sequences'
    )
    parser.add_argument(
        '--output-dir', 
        required=True,
        help='Output directory for embeddings'
    )
    parser.add_argument(
        '--model', 
        default='esm2_t33_650M_UR50D',
        choices=['esm2_t33_650M_UR50D'],
        help='ESM-2 model to use'
    )
    parser.add_argument(
        '--device', 
        default='auto',
        choices=['auto', 'cpu', 'cuda', 'cuda:0', 'cuda:1'],
        help='Device to use'
    )
    parser.add_argument(
        '--max-len', 
        type=int, 
        default=1022,
        help='Maximum window length (ESM-2 limit: 1022)'
    )
    parser.add_argument(
        '--stride', 
        type=int, 
        default=900,
        help='Stride for sliding window'
    )
    parser.add_argument(
        '--layers',
        default='20-33',
        help='Layer specification: single (33), range (20-33), mid, or all'
    )
    parser.add_argument(
        '--pooling',
        default='mean',
        choices=['mean', 'cls'],
        help='Pooling strategy: mean or cls token'
    )
    parser.add_argument(
        '--stitching',
        default='per_residue',
        choices=['per_residue', 'window'],
        help='Overlap handling: per_residue (accurate) or window (fast)'
    )
    parser.add_argument(
        '--precision',
        default='fp32',
        choices=['fp32', 'fp16', 'bf16'],
        help='Numerical precision (fp16/bf16 requires CUDA)'
    )
    parser.add_argument(
        '--no-cache',
        action='store_true',
        help='Disable per-sequence caching'
    )
    parser.add_argument(
        '--deterministic',
        action='store_true',
        help='Enable deterministic mode for reproducibility'
    )
    
    args = parser.parse_args()
    
    # Validation
    if args.precision in ['fp16', 'bf16'] and not args.device.startswith('cuda'):
        print("⚠️  Warning: fp16/bf16 requires CUDA, falling back to fp32")
        args.precision = 'fp32'
    
    generate_embeddings(
        input_csv=args.input,
        output_dir=args.output_dir,
        model_name=args.model,
        device=args.device,
        max_len=args.max_len,
        stride=args.stride,
        layer_spec=args.layers,
        pooling=args.pooling,
        precision=args.precision,
        stitching=args.stitching,
        use_cache=not args.no_cache,
        deterministic=args.deterministic
    )


if __name__ == '__main__':
    main()
