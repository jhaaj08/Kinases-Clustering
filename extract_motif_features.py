#!/usr/bin/env python3
"""
Extract kinase motif and structural features from domain sequences.

Features extracted:
1. Motif presence/absence (binary):
   - DFG motif ([D][F][G])
   - HRD motif ([H][R][D])
   - APE motif ([A][P][E])
   - P-loop (GxGxxG pattern)
   - β3-Lys motif (VAIK region)
   - αC helix acidic residue (E or D)

2. Quantitative features:
   - Activation loop length (DFG → APE distance)
   - Catalytic loop length (HRD → DFG distance)
   - Gatekeeper residue identity and size
   - Motif positions (normalized by sequence length)

3. Structural features (when available):
   - Gatekeeper hydrophobicity
   - P-loop consensus score
"""

import re
import pandas as pd
import numpy as np
from typing import Dict, Tuple, Optional


# Amino acid properties
AA_HYDROPHOBICITY = {
    'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5, 'F': 2.8,
    'G': -0.4, 'H': -3.2, 'I': 4.5, 'K': -3.9, 'L': 3.8,
    'M': 1.9, 'N': -3.5, 'P': -1.6, 'Q': -3.5, 'R': -4.5,
    'S': -0.8, 'T': -0.7, 'V': 4.2, 'W': -0.9, 'Y': -1.3
}

AA_SIZE = {
    'A': 1, 'C': 2, 'D': 4, 'E': 5, 'F': 8,
    'G': 0, 'H': 6, 'I': 4, 'K': 5, 'L': 4,
    'M': 4, 'N': 4, 'P': 3, 'Q': 5, 'R': 7,
    'S': 2, 'T': 3, 'V': 3, 'W': 10, 'Y': 9
}


def find_motif(sequence: str, pattern: str, motif_name: str) -> Tuple[bool, int]:
    """
    Find a motif using regex pattern.
    
    Returns:
        (found: bool, position: int) - position is -1 if not found
    """
    match = re.search(pattern, sequence)
    if match:
        return True, match.start()
    return False, -1


def find_ploop(sequence: str) -> Tuple[bool, int, float]:
    """
    Find P-loop motif (GxGxxG pattern) and calculate consensus score.
    
    Returns:
        (found: bool, position: int, consensus_score: float)
    """
    # Classic P-loop: G-x-G-x-x-G (positions 0, 2, 5 should be G)
    pattern = r'G.G..G'
    match = re.search(pattern, sequence)
    
    if match:
        motif_seq = match.group()
        # Calculate consensus score (how many Gs are in correct positions)
        score = sum([
            motif_seq[0] == 'G',  # Position 0
            motif_seq[2] == 'G',  # Position 2
            motif_seq[5] == 'G',  # Position 5
        ]) / 3.0
        return True, match.start(), score
    
    return False, -1, 0.0


def find_vaik_region(sequence: str) -> Tuple[bool, int, int]:
    """
    Find β3-Lys motif (VAIK region - characteristic of kinases).
    This is a more flexible search for the conserved lysine region.
    
    Returns:
        (found: bool, position: int, lysine_position: int)
    """
    # Look for VAIK or similar patterns (V/I/L-A/G-I/V-K)
    patterns = [
        r'VAIK',
        r'V[AG][IV]K',
        r'[VIL]A[IV]K',
        r'[VIL][AG][IV]K'
    ]
    
    for pattern in patterns:
        match = re.search(pattern, sequence)
        if match:
            # Position of K in the motif (usually at index 3)
            motif_seq = match.group()
            k_offset = motif_seq.index('K')
            k_pos = match.start() + k_offset
            return True, match.start(), k_pos
    
    return False, -1, -1


def calculate_k_e_salt_bridge_distance(sequence: str, vaik_k_pos: int) -> Tuple[int, bool]:
    """
    Calculate K-E salt bridge distance (sequence-based proxy).
    
    In kinases, the conserved Lys (β3) forms a salt bridge with a Glu in the αC helix.
    This is typically 20-40 residues downstream of the β3-Lys.
    
    Args:
        sequence: Protein sequence
        vaik_k_pos: Position of the conserved Lys from VAIK motif
    
    Returns:
        (distance: int, found: bool) - distance to nearest acidic residue (E/D)
    """
    if vaik_k_pos == -1:
        return -1, False
    
    # Look for acidic residue (E preferred, D acceptable) within 20-50 residues downstream
    search_start = vaik_k_pos + 15
    search_end = min(vaik_k_pos + 55, len(sequence))
    
    if search_start >= len(sequence):
        return -1, False
    
    search_region = sequence[search_start:search_end]
    
    # Find first E or D
    for i, aa in enumerate(search_region):
        if aa in ['E', 'D']:
            distance = search_start - vaik_k_pos + i
            return distance, True
    
    return -1, False


def find_alphac_acidic(sequence: str, search_window: int = 100) -> Tuple[bool, int]:
    """
    Find acidic residue (E or D) in αC helix region.
    We search in the first ~100 residues where αC typically resides.
    """
    search_region = sequence[:search_window]
    
    # Look for conserved glutamate (E) first, then aspartate (D)
    for aa in ['E', 'D']:
        pos = search_region.find(aa)
        if pos != -1:
            return True, pos
    
    return False, -1


def calculate_gatekeeper_features(sequence: str, dfg_pos: int) -> Dict[str, float]:
    """
    Calculate gatekeeper residue features.
    Gatekeeper is typically ~15 residues before DFG motif.
    
    Returns:
        Dictionary with gatekeeper identity, size, and hydrophobicity
    """
    features = {
        'gatekeeper_found': 0,
        'gatekeeper_size': 0,
        'gatekeeper_hydrophobicity': 0,
        'gatekeeper_is_small': 0,  # T, S, A, G
        'gatekeeper_is_large': 0,  # M, F, Y, W
    }
    
    if dfg_pos < 0:
        return features
    
    # Gatekeeper is typically at DFG-15 position
    gk_pos = dfg_pos - 15
    
    if 0 <= gk_pos < len(sequence):
        gk_residue = sequence[gk_pos]
        features['gatekeeper_found'] = 1
        features['gatekeeper_size'] = AA_SIZE.get(gk_residue, 0)
        features['gatekeeper_hydrophobicity'] = AA_HYDROPHOBICITY.get(gk_residue, 0)
        features['gatekeeper_is_small'] = 1 if gk_residue in 'TSAG' else 0
        features['gatekeeper_is_large'] = 1 if gk_residue in 'MFYW' else 0
    
    return features


def extract_motif_features(sequence: str) -> Dict[str, float]:
    """
    Extract all motif and structural features from a kinase domain sequence.
    
    Returns:
        Dictionary of feature name -> value
    """
    features = {}
    
    # 1. Core motifs (binary presence)
    dfg_found, dfg_pos = find_motif(sequence, r'DFG', 'DFG')
    hrd_found, hrd_pos = find_motif(sequence, r'HRD', 'HRD')
    ape_found, ape_pos = find_motif(sequence, r'APE', 'APE')
    
    features['dfg_present'] = int(dfg_found)
    features['hrd_present'] = int(hrd_found)
    features['ape_present'] = int(ape_found)
    
    # 2. P-loop
    ploop_found, ploop_pos, ploop_consensus = find_ploop(sequence)
    features['ploop_present'] = int(ploop_found)
    features['ploop_consensus'] = ploop_consensus
    
    # 3. β3-Lys (VAIK region)
    vaik_found, vaik_pos, vaik_k_pos = find_vaik_region(sequence)
    features['vaik_present'] = int(vaik_found)
    
    # 3b. K-E salt bridge distance (NEW - reviewer requested)
    k_e_distance, k_e_found = calculate_k_e_salt_bridge_distance(sequence, vaik_k_pos)
    features['k_e_distance'] = k_e_distance if k_e_found else -1
    features['k_e_salt_bridge_present'] = int(k_e_found)
    features['k_e_distance_normal'] = int(25 <= k_e_distance <= 40) if k_e_found else 0  # Typical range
    
    # 4. αC helix acidic residue
    alphac_found, alphac_pos = find_alphac_acidic(sequence)
    features['alphac_acidic_present'] = int(alphac_found)
    
    # 5. Loop lengths (quantitative)
    seq_len = len(sequence)
    
    # Activation loop: DFG → APE
    if dfg_found and ape_found and ape_pos > dfg_pos:
        activation_loop_length = ape_pos - dfg_pos
        features['activation_loop_length'] = activation_loop_length
        features['activation_loop_length_norm'] = activation_loop_length / seq_len
    else:
        features['activation_loop_length'] = 0
        features['activation_loop_length_norm'] = 0
    
    # Catalytic loop: HRD → DFG
    if hrd_found and dfg_found and dfg_pos > hrd_pos:
        catalytic_loop_length = dfg_pos - hrd_pos
        features['catalytic_loop_length'] = catalytic_loop_length
        features['catalytic_loop_length_norm'] = catalytic_loop_length / seq_len
    else:
        features['catalytic_loop_length'] = 0
        features['catalytic_loop_length_norm'] = 0
    
    # 6. Motif positions (normalized)
    features['dfg_position_norm'] = dfg_pos / seq_len if dfg_found else -1
    features['hrd_position_norm'] = hrd_pos / seq_len if hrd_found else -1
    features['ape_position_norm'] = ape_pos / seq_len if ape_found else -1
    features['ploop_position_norm'] = ploop_pos / seq_len if ploop_found else -1
    
    # 7. Gatekeeper features
    gk_features = calculate_gatekeeper_features(sequence, dfg_pos)
    features.update(gk_features)
    
    # 8. HRD/DFG state features (NEW - reviewer requested)
    # DFG conformation proxy: check residues around DFG
    if dfg_found and dfg_pos + 5 < seq_len:
        # DFG-in state typically has specific residue patterns
        post_dfg_region = sequence[dfg_pos+3:dfg_pos+6]
        features['dfg_state_hydrophobic'] = sum(1 for aa in post_dfg_region if aa in 'FVILM') / len(post_dfg_region)
    else:
        features['dfg_state_hydrophobic'] = 0
    
    # HRD-DFG spacing (another indicator of active/inactive state)
    if hrd_found and dfg_found and hrd_pos < dfg_pos:
        hrd_dfg_spacing = dfg_pos - hrd_pos
        features['hrd_dfg_spacing'] = hrd_dfg_spacing
        features['hrd_dfg_spacing_normal'] = int(20 <= hrd_dfg_spacing <= 60)  # Typical range
    else:
        features['hrd_dfg_spacing'] = -1
        features['hrd_dfg_spacing_normal'] = 0
    
    # 9. Core motif triad (all three present = canonical kinase)
    features['core_triad_complete'] = int(dfg_found and hrd_found and ape_found)
    
    # Extended motif completeness (including VAIK and K-E)
    extended_motifs_present = sum([
        dfg_found, hrd_found, ape_found, vaik_found, k_e_found
    ])
    features['extended_motif_completeness'] = extended_motifs_present / 5.0
    
    # 10. Overall motif integrity score (for per-sequence reports)
    # This score combines presence and proper spacing of key motifs
    integrity_score = (
        features['core_triad_complete'] * 0.3 +
        features['k_e_salt_bridge_present'] * 0.2 +
        features['ploop_present'] * 0.15 +
        features['vaik_present'] * 0.15 +
        features['hrd_dfg_spacing_normal'] * 0.1 +
        features['k_e_distance_normal'] * 0.1
    )
    features['motif_integrity_score'] = integrity_score
    
    # 11. Sequence-level features
    features['sequence_length'] = seq_len
    
    return features


def extract_all_features(df: pd.DataFrame, sequence_col: str = 'sequence') -> pd.DataFrame:
    """
    Extract motif features for all sequences in dataframe.
    
    Args:
        df: DataFrame with sequence column
        sequence_col: Name of sequence column
    
    Returns:
        DataFrame with original columns + motif feature columns
    """
    print("Extracting motif features from domain sequences...")
    
    feature_list = []
    for idx, row in df.iterrows():
        seq = row[sequence_col]
        features = extract_motif_features(seq)
        features['uniprot_id'] = row['uniprot_id']
        feature_list.append(features)
    
    features_df = pd.DataFrame(feature_list)
    
    # Merge with original dataframe
    result_df = df.merge(features_df, on='uniprot_id', how='left')
    
    print(f"✅ Extracted {len(features_df.columns)-1} motif features")
    print(f"   Total columns: {len(result_df.columns)}")
    
    return result_df


def get_feature_statistics(df: pd.DataFrame) -> None:
    """Print statistics about extracted features."""
    
    print("\n" + "="*80)
    print("MOTIF FEATURE STATISTICS")
    print("="*80)
    
    binary_features = [
        'dfg_present', 'hrd_present', 'ape_present',
        'ploop_present', 'vaik_present', 'alphac_acidic_present',
        'core_triad_complete', 'gatekeeper_found'
    ]
    
    print("\nBinary Features (% present):")
    print("-"*80)
    for feat in binary_features:
        if feat in df.columns:
            pct = df[feat].mean() * 100
            count = df[feat].sum()
            total = len(df)
            print(f"  {feat:<30s}: {count:4d}/{total} ({pct:5.1f}%)")
    
    print("\nQuantitative Features (mean ± std):")
    print("-"*80)
    quant_features = [
        'activation_loop_length', 'catalytic_loop_length',
        'ploop_consensus', 'gatekeeper_size', 'gatekeeper_hydrophobicity',
        'sequence_length'
    ]
    
    for feat in quant_features:
        if feat in df.columns:
            mean = df[feat].mean()
            std = df[feat].std()
            print(f"  {feat:<30s}: {mean:7.2f} ± {std:6.2f}")
    
    print()


def main():
    """Extract motif features and save to CSV."""
    
    print("="*80)
    print("KINASE MOTIF FEATURE EXTRACTION")
    print("="*80)
    print()
    
    # Load domain sequences
    input_file = "kinases_domains.csv"
    output_file = "kinases_domains_with_motifs.csv"
    
    print(f"Loading {input_file}...")
    df = pd.read_csv(input_file)
    print(f"✅ Loaded {len(df)} domain sequences")
    
    # Extract features
    df_with_motifs = extract_all_features(df)
    
    # Save
    df_with_motifs.to_csv(output_file, index=False)
    print(f"\n✅ Saved {len(df_with_motifs)} sequences with motif features to {output_file}")
    
    # Statistics
    get_feature_statistics(df_with_motifs)
    
    # List feature columns
    feature_cols = [col for col in df_with_motifs.columns if col not in df.columns]
    print("Feature columns added:")
    print("-"*80)
    for i, col in enumerate(feature_cols, 1):
        print(f"  {i:2d}. {col}")
    
    print("\n" + "="*80)
    print("✅ MOTIF EXTRACTION COMPLETE!")
    print("="*80)
    print()
    print("Next: Fuse with ESM-2 embeddings and re-cluster")
    print()


if __name__ == "__main__":
    main()

