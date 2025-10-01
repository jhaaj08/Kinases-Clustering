"""
Generate ESM-2 embeddings for all kinase sequences using sliding window approach.
Uses length-weighted averaging to combine windows for long sequences.
"""

import argparse
import os
import numpy as np
import pandas as pd
import torch
from tqdm import tqdm
import esm


def embed_sequence_windowed(seq, model, alphabet, device, max_len=1022, stride=900, layer=33):
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
        Maximum sequence length per window (default: 1022)
    stride : int
        Stride for sliding window (default: 900)
    layer : int
        Layer to extract representations from (default: 33 for ESM-2)
    
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
        """Embed a single chunk and return mean-pooled vector and residue count."""
        data = [("seq", chunk_seq)]
        labels, strs, tokens = batch_converter(data)
        tokens = tokens.to(device)
        
        with torch.no_grad():
            results = model(tokens, repr_layers=[layer], return_contacts=False)
        
        # Get representations: (batch=1, seq_len, dim)
        reps = results["representations"][layer][0]  # (seq_len, dim)
        
        # Mask out special tokens
        tok = tokens[0]
        mask = (tok != pad_idx) & (tok != cls_idx) & (tok != eos_idx)
        
        # Get only residue embeddings
        residue_reps = reps[mask]  # (num_residues, dim)
        
        # Mean pool
        mean_vec = residue_reps.mean(dim=0)  # (dim,)
        num_residues = residue_reps.shape[0]
        
        return mean_vec, num_residues
    
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
        
        # Weight by number of residues
        weighted_sum += vec * num_res
        total_weight += num_res
        
        # Move to next window
        if end == len(seq):
            break
        start += stride
    
    # Length-weighted average
    final_vec = weighted_sum / total_weight
    
    return final_vec.cpu().numpy()


def generate_embeddings(input_csv, output_dir, model_name='esm2_t33_650M_UR50D', 
                       device='auto', max_len=1022, stride=900):
    """
    Generate ESM-2 embeddings for all sequences in CSV file.
    
    Parameters:
    -----------
    input_csv : str
        Path to input CSV file
    output_dir : str
        Output directory for embeddings
    model_name : str
        ESM model name (default: esm2_t33_650M_UR50D)
    device : str
        Device to use ('auto', 'cpu', or 'cuda')
    max_len : int
        Maximum window length (default: 1022)
    stride : int
        Stride for sliding window (default: 900)
    """
    
    print("=" * 80)
    print("ESM-2 EMBEDDING GENERATION")
    print("=" * 80)
    print()
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Read input data
    print(f"Reading {input_csv}...")
    df = pd.read_csv(input_csv)
    print(f"  Loaded {len(df):,} sequences")
    
    # Check sequence lengths
    df['seq_len'] = df['sequence'].apply(len)
    num_long = (df['seq_len'] > max_len).sum()
    print(f"  {num_long:,} sequences need windowing (>{max_len} aa)")
    print()
    
    # Determine device
    if device == 'auto':
        device = 'cuda' if torch.cuda.is_available() else 'cpu'
    print(f"Device: {device}")
    
    # Load ESM-2 model
    print(f"Loading ESM-2 model: {model_name}...")
    if model_name == 'esm2_t33_650M_UR50D':
        model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
        layer = 33
        embedding_dim = 1280
    else:
        raise ValueError(f"Unsupported model: {model_name}")
    
    model.eval()
    model.to(device)
    print(f"  Model loaded (embedding dim: {embedding_dim})")
    print()
    
    # Generate embeddings
    print(f"Generating embeddings (window={max_len}, stride={stride})...")
    print("This may take 10-30 minutes depending on your hardware...")
    print()
    
    embeddings = np.zeros((len(df), embedding_dim), dtype=np.float32)
    
    for idx, row in tqdm(df.iterrows(), total=len(df), desc="Processing"):
        seq = row['sequence']
        emb = embed_sequence_windowed(
            seq, model, alphabet, device, 
            max_len=max_len, stride=stride, layer=layer
        )
        embeddings[idx] = emb
    
    print()
    print("Saving embeddings...")
    
    # Save embeddings as numpy array
    emb_file = os.path.join(output_dir, 'esm2_embeddings.npy')
    np.save(emb_file, embeddings)
    print(f"  Saved: {emb_file}")
    print(f"  Shape: {embeddings.shape}")
    print(f"  Size: {os.path.getsize(emb_file) / 1024 / 1024:.1f} MB")
    
    # Save index with uniprot_id
    index_df = df[['uniprot_id']].copy()
    index_file = os.path.join(output_dir, 'esm2_index.csv')
    index_df.to_csv(index_file, index=False)
    print(f"  Saved: {index_file}")
    
    # Save combined file with embeddings as columns
    print()
    print("Creating combined CSV with embeddings...")
    emb_df = pd.DataFrame(
        embeddings, 
        columns=[f'emb_{i}' for i in range(embedding_dim)]
    )
    combined_df = pd.concat([df[['uniprot_id']], emb_df], axis=1)
    combined_file = os.path.join(output_dir, 'kinases_with_embeddings.csv')
    combined_df.to_csv(combined_file, index=False)
    print(f"  Saved: {combined_file}")
    print(f"  Columns: uniprot_id + {embedding_dim} embedding dimensions")
    
    # Statistics
    print()
    print("=" * 80)
    print("EMBEDDING STATISTICS")
    print("=" * 80)
    print(f"  Shape: {embeddings.shape}")
    print(f"  Mean: {embeddings.mean():.6f}")
    print(f"  Std:  {embeddings.std():.6f}")
    print(f"  Min:  {embeddings.min():.6f}")
    print(f"  Max:  {embeddings.max():.6f}")
    print()
    
    print("✅ Embedding generation complete!")
    print()
    print("Output files:")
    print(f"  1. {emb_file} - numpy array (6465, 1280)")
    print(f"  2. {index_file} - uniprot_id index")
    print(f"  3. {combined_file} - CSV with uniprot_id + embeddings")
    print()
    
    return embeddings


def main():
    parser = argparse.ArgumentParser(
        description='Generate ESM-2 embeddings for kinase sequences'
    )
    parser.add_argument(
        '--input', 
        default='kinases_revised.csv',
        help='Input CSV file with sequences'
    )
    parser.add_argument(
        '--output-dir', 
        default='kinases_embeddings',
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
        help='Device to use (auto, cpu, or cuda)'
    )
    parser.add_argument(
        '--max-len', 
        type=int, 
        default=1022,
        help='Maximum window length (default: 1022)'
    )
    parser.add_argument(
        '--stride', 
        type=int, 
        default=900,
        help='Stride for sliding window (default: 900)'
    )
    
    args = parser.parse_args()
    
    generate_embeddings(
        input_csv=args.input,
        output_dir=args.output_dir,
        model_name=args.model,
        device=args.device,
        max_len=args.max_len,
        stride=args.stride
    )


if __name__ == '__main__':
    main()
