#!/usr/bin/env python3
"""
Enhanced ESM-2 embedding generation with layer probing and pooling options.

New features:
- Layer selection: single layer, mean of layers, or specific range
- Pooling strategy: mean pooling or CLS token
- Command-line arguments for all options
"""

import argparse
import os
import numpy as np
import pandas as pd
import torch
from tqdm import tqdm
import esm


def embed_sequence_windowed(seq, model, alphabet, device, 
                            max_len=1022, stride=900, 
                            layers=[33], pooling='mean'):
    """
    Generate embedding for a single sequence using sliding window approach.
    
    Parameters:
    -----------
    seq : str
        Amino acid sequence
    model : torch.nn.Module
        ESM model
    alphabet : esm.Alphabet
        ESM alphabet
    device : str
        Device to use ('cpu' or 'cuda')
    max_len : int
        Maximum sequence length per window
    stride : int
        Stride for sliding window
    layers : list of int
        Layers to extract (will be averaged if multiple)
    pooling : str
        Pooling strategy: 'mean' or 'cls'
    
    Returns:
    --------
    np.ndarray
        Embedding vector of shape (embedding_dim,)
    """
    batch_converter = alphabet.get_batch_converter()
    pad_idx = alphabet.padding_idx
    cls_idx = alphabet.cls_idx
    eos_idx = alphabet.eos_idx
    
    def embed_chunk(chunk_seq):
        """Embed a single chunk and return pooled vector and residue count."""
        data = [("seq", chunk_seq)]
        labels, strs, tokens = batch_converter(data)
        tokens = tokens.to(device)
        
        with torch.no_grad():
            results = model(tokens, repr_layers=layers, return_contacts=False)
        
        # Average across layers if multiple
        if len(layers) > 1:
            layer_reps = [results["representations"][layer][0] for layer in layers]
            reps = torch.stack(layer_reps).mean(dim=0)  # (seq_len, dim)
        else:
            reps = results["representations"][layers[0]][0]  # (seq_len, dim)
        
        # Pooling strategy
        if pooling == 'cls':
            # Use CLS token (first token)
            pooled_vec = reps[0]  # (dim,)
            num_residues = 1
        else:  # mean
            # Mask out special tokens
            tok = tokens[0]
            mask = (tok != pad_idx) & (tok != cls_idx) & (tok != eos_idx)
            
            # Get only residue embeddings
            residue_reps = reps[mask]  # (num_residues, dim)
            
            # Mean pool
            pooled_vec = residue_reps.mean(dim=0)  # (dim,)
            num_residues = residue_reps.shape[0]
        
        return pooled_vec, num_residues
    
    # If sequence fits in one window
    if len(seq) <= max_len:
        vec, _ = embed_chunk(seq)
        return vec.cpu().numpy()
    
    # Sliding window for long sequences
    embedding_dim = 1280  # ESM-2 650M has 1280 dimensions
    weighted_sum = torch.zeros(embedding_dim, device=device)
    total_weight = 0
    
    start = 0
    while start < len(seq):
        end = min(start + max_len, len(seq))
        chunk = seq[start:end]
        
        vec, num_res = embed_chunk(chunk)
        
        # Weight by number of residues (only for mean pooling)
        if pooling == 'mean':
            weighted_sum += vec * num_res
            total_weight += num_res
        else:  # cls: simple average
            weighted_sum += vec
            total_weight += 1
        
        # Move to next window
        if end == len(seq):
            break
        start += stride
    
    # Weighted average
    final_vec = weighted_sum / total_weight
    
    return final_vec.cpu().numpy()


def parse_layer_spec(layer_spec, model_layers=33):
    """
    Parse layer specification string.
    
    Examples:
    - '33' -> [33] (last layer)
    - '20-30' -> [20, 21, ..., 30] (range)
    - 'mid' -> [20, 21, ..., 30] (middle layers)
    - 'all' -> [1, 2, ..., 33] (all layers)
    """
    if layer_spec == 'all':
        return list(range(1, model_layers + 1))
    elif layer_spec == 'mid':
        start = int(model_layers * 0.6)  # ~20 for 33 layers
        return list(range(start, model_layers + 1))
    elif '-' in layer_spec:
        start, end = map(int, layer_spec.split('-'))
        return list(range(start, end + 1))
    else:
        return [int(layer_spec)]


def generate_embeddings(input_csv, output_dir, model_name='esm2_t33_650M_UR50D', 
                       device='auto', max_len=1022, stride=900,
                       layer_spec='33', pooling='mean'):
    """
    Generate ESM-2 embeddings with specified layer and pooling options.
    """
    
    print("=" * 80)
    print("ENHANCED ESM-2 EMBEDDING GENERATION")
    print("=" * 80)
    print()
    print(f"Input:    {input_csv}")
    print(f"Output:   {output_dir}")
    print(f"Layers:   {layer_spec}")
    print(f"Pooling:  {pooling}")
    print()
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
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
    
    # Load ESM-2 model
    print(f"Loading ESM-2 model: {model_name}...")
    if model_name == 'esm2_t33_650M_UR50D':
        model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
        model_layers = 33
        embedding_dim = 1280
    else:
        raise ValueError(f"Unsupported model: {model_name}")
    
    model.eval()
    model.to(device)
    
    # Parse layer specification
    layers = parse_layer_spec(layer_spec, model_layers)
    print(f"  Layers to extract: {layers if len(layers) <= 5 else f'{layers[0]}-{layers[-1]} ({len(layers)} layers)'}")
    print(f"  Pooling strategy: {pooling}")
    print(f"  Embedding dim: {embedding_dim}")
    print()
    
    # Generate embeddings
    print(f"Generating embeddings (window={max_len}, stride={stride})...")
    print()
    
    embeddings = np.zeros((len(df), embedding_dim), dtype=np.float32)
    
    for idx, row in tqdm(df.iterrows(), total=len(df), desc="Processing"):
        seq = row['sequence']
        emb = embed_sequence_windowed(
            seq, model, alphabet, device, 
            max_len=max_len, stride=stride, 
            layers=layers, pooling=pooling
        )
        embeddings[idx] = emb
    
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
    
    # Save metadata
    meta_file = os.path.join(output_dir, 'embedding_metadata.txt')
    with open(meta_file, 'w') as f:
        f.write(f"Model: {model_name}\n")
        f.write(f"Layers: {layer_spec} -> {layers}\n")
        f.write(f"Pooling: {pooling}\n")
        f.write(f"Shape: {embeddings.shape}\n")
        f.write(f"Max length: {max_len}\n")
        f.write(f"Stride: {stride}\n")
    print(f"  ✅ Saved: {meta_file}")
    
    # Statistics
    print()
    print("=" * 80)
    print("EMBEDDING STATISTICS")
    print("=" * 80)
    print(f"  Shape: {embeddings.shape}")
    print(f"  Mean:  {embeddings.mean():.6f}")
    print(f"  Std:   {embeddings.std():.6f}")
    print(f"  Min:   {embeddings.min():.6f}")
    print(f"  Max:   {embeddings.max():.6f}")
    print()
    print("✅ Embedding generation complete!")
    print()
    
    return embeddings


def main():
    parser = argparse.ArgumentParser(
        description='Generate ESM-2 embeddings with layer and pooling options'
    )
    parser.add_argument(
        '--input', 
        default='kinases_domains.csv',
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
        choices=['auto', 'cpu', 'cuda'],
        help='Device to use'
    )
    parser.add_argument(
        '--max-len', 
        type=int, 
        default=1022,
        help='Maximum window length'
    )
    parser.add_argument(
        '--stride', 
        type=int, 
        default=900,
        help='Stride for sliding window'
    )
    parser.add_argument(
        '--layers',
        default='33',
        help='Layer specification: single (33), range (20-30), mid, or all'
    )
    parser.add_argument(
        '--pooling',
        default='mean',
        choices=['mean', 'cls'],
        help='Pooling strategy: mean or cls token'
    )
    
    args = parser.parse_args()
    
    generate_embeddings(
        input_csv=args.input,
        output_dir=args.output_dir,
        model_name=args.model,
        device=args.device,
        max_len=args.max_len,
        stride=args.stride,
        layer_spec=args.layers,
        pooling=args.pooling
    )


if __name__ == '__main__':
    main()

