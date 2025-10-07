#!/usr/bin/env python3
"""
Generate homology-aware train/test splits to prevent data leakage.

Approach:
1. Cluster sequences at 40% identity using CD-HIT
2. Assign cluster IDs to all sequences
3. Use GroupShuffleSplit to ensure no cluster spans train/test
4. Maintain stratification of kinase families
5. Save splits to data/splits.json for reproducibility

This prevents information leakage from homologous sequences.
"""

import os
import sys
import argparse
import subprocess
import tempfile
import pandas as pd
import numpy as np
import json
from pathlib import Path
from collections import Counter
from sklearn.model_selection import StratifiedGroupKFold
from utils.provenance import ProvenanceTracker


def check_cdhit_installed():
    """Check if CD-HIT is installed."""
    try:
        subprocess.run(['cd-hit', '-h'], capture_output=True, timeout=5)
        return True
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


def run_cdhit_clustering(sequences_df, identity=0.4, output_prefix="cdhit_splits"):
    """
    Run CD-HIT to cluster sequences at specified identity.
    
    Args:
        sequences_df: DataFrame with sequences
        identity: Sequence identity threshold (0.0-1.0)
        output_prefix: Prefix for output files
    
    Returns:
        DataFrame with cluster_id column added
    """
    print(f"\nRunning CD-HIT clustering (identity={identity:.0%})...")
    
    # Create temporary FASTA
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        fasta_file = f.name
        for _, row in sequences_df.iterrows():
            f.write(f">{row['uniprot_id']}\n")
            f.write(f"{row['sequence']}\n")
    
    output_file = f"{output_prefix}_{int(identity*100)}.fasta"
    cluster_file = f"{output_file}.clstr"
    
    # Run CD-HIT
    cmd = [
        'cd-hit',
        '-i', fasta_file,
        '-o', output_file,
        '-c', str(identity),
        '-n', '2',  # Word size for ~40% identity
        '-M', '0',  # Unlimited memory
        '-T', '0',  # Use all CPUs
        '-d', '0',  # Full header in output
    ]
    
    print(f"  Command: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        if result.returncode != 0:
            print(f"  ❌ CD-HIT failed: {result.stderr}")
            return None
    except subprocess.TimeoutExpired:
        print("  ❌ CD-HIT timed out")
        return None
    
    print(f"  ✅ CD-HIT complete")
    
    # Parse cluster file
    cluster_map = {}  # uniprot_id -> cluster_id
    current_cluster = -1
    
    with open(cluster_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>Cluster'):
                current_cluster = int(line.split()[1])
            elif '\t' in line and '>' in line and current_cluster >= 0:
                # Extract uniprot_id from cluster line
                # Format: "0	408aa, >P32350... *" or "1	321aa, >Q10156... at 43.30%"
                parts = line.split('>')
                if len(parts) > 1:
                    # Get text after '>' and before '...'
                    uniprot_id = parts[1].split('...')[0]
                    cluster_map[uniprot_id] = current_cluster
    
    # Add cluster IDs to dataframe
    sequences_df['cluster_id'] = sequences_df['uniprot_id'].map(cluster_map)
    
    # Cleanup
    os.unlink(fasta_file)
    
    n_clusters = len(set(cluster_map.values()))
    print(f"  Clustered {len(cluster_map)} sequences into {n_clusters} groups")
    
    # Statistics
    if cluster_map:
        cluster_sizes = Counter(cluster_map.values())
        sizes_list = list(cluster_sizes.values())
        print(f"  Cluster size: mean={np.mean(sizes_list):.1f}, "
              f"median={np.median(sizes_list):.0f}, "
              f"max={max(sizes_list)}")
    else:
        print("  ⚠️  No clusters found!")
    
    return sequences_df


def create_homology_aware_split(df, test_size=0.2, random_state=42, min_class_size=5):
    """
    Create homology-aware train/test split.
    
    Args:
        df: DataFrame with cluster_id and kinome_group_major columns
        test_size: Fraction for test set
        random_state: Random seed
        min_class_size: Minimum samples per class to include
    
    Returns:
        train_ids, test_ids, split_metadata
    """
    print(f"\nCreating homology-aware split...")
    print(f"  Test size: {test_size*100:.0f}%")
    print(f"  Random seed: {random_state}")
    print(f"  Min class size: {min_class_size}")
    
    # Filter out small classes
    class_counts = df['kinome_group_major'].value_counts()
    valid_classes = class_counts[class_counts >= min_class_size].index.tolist()
    
    df_filtered = df[df['kinome_group_major'].isin(valid_classes)].copy()
    
    print(f"  Classes before filtering: {len(class_counts)}")
    print(f"  Classes after filtering: {len(valid_classes)}")
    print(f"  Samples: {len(df)} → {len(df_filtered)}")
    
    # Ensure cluster_id and labels are available
    if 'cluster_id' not in df_filtered.columns:
        print("  ⚠️  No cluster_id found, will use identity as group")
        df_filtered['cluster_id'] = range(len(df_filtered))
    
    # Use StratifiedGroupKFold with n_splits=5, take first split
    # This ensures stratification while respecting cluster boundaries
    sgkf = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=random_state)
    
    X = np.zeros((len(df_filtered), 1))  # Dummy features
    y = df_filtered['kinome_group_major'].values
    groups = df_filtered['cluster_id'].values
    
    # Get first split (approximately 20% test)
    for train_idx, test_idx in sgkf.split(X, y, groups):
        # Use first split
        break
    
    train_ids = df_filtered.iloc[train_idx]['uniprot_id'].tolist()
    test_ids = df_filtered.iloc[test_idx]['uniprot_id'].tolist()
    
    # Verify no cluster spans train and test
    train_clusters = set(df_filtered.iloc[train_idx]['cluster_id'])
    test_clusters = set(df_filtered.iloc[test_idx]['cluster_id'])
    overlap = train_clusters & test_clusters
    
    if overlap:
        print(f"  ⚠️  WARNING: {len(overlap)} clusters span train/test!")
    else:
        print(f"  ✅ No cluster overlap between train/test")
    
    # Class distribution
    train_labels = df_filtered.iloc[train_idx]['kinome_group_major']
    test_labels = df_filtered.iloc[test_idx]['kinome_group_major']
    
    print(f"\n  Train set: {len(train_ids)} samples")
    for cls, count in train_labels.value_counts().items():
        print(f"    {cls:12s}: {count:4d} ({count/len(train_ids)*100:5.1f}%)")
    
    print(f"\n  Test set: {len(test_ids)} samples")
    for cls, count in test_labels.value_counts().items():
        print(f"    {cls:12s}: {count:4d} ({count/len(test_ids)*100:5.1f}%)")
    
    # Metadata
    split_metadata = {
        "method": "StratifiedGroupKFold (n_splits=5, first fold as test)",
        "homology_threshold": "40% identity (CD-HIT)",
        "stratified": True,
        "random_state": random_state,
        "test_size_actual": len(test_ids) / (len(train_ids) + len(test_ids)),
        "test_size_target": test_size,
        "train_n": len(train_ids),
        "test_n": len(test_ids),
        "total_clusters": len(set(groups)),
        "train_clusters": len(train_clusters),
        "test_clusters": len(test_clusters),
        "cluster_overlap": len(overlap),
        "classes_included": valid_classes,
        "classes_excluded": class_counts[class_counts < min_class_size].index.tolist(),
        "train_class_distribution": train_labels.value_counts().to_dict(),
        "test_class_distribution": test_labels.value_counts().to_dict(),
    }
    
    return train_ids, test_ids, split_metadata


def save_splits(train_ids, test_ids, metadata, output_dir="data"):
    """Save split indices and metadata to JSON."""
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    splits_file = output_dir / "splits.json"
    
    splits_data = {
        "created_at": pd.Timestamp.now().isoformat(),
        "train_ids": train_ids,
        "test_ids": test_ids,
        "metadata": metadata,
    }
    
    with open(splits_file, 'w') as f:
        json.dump(splits_data, f, indent=2)
    
    print(f"\n✅ Saved splits to: {splits_file}")
    
    return splits_file


def main():
    parser = argparse.ArgumentParser(
        description='Generate homology-aware train/test splits'
    )
    parser.add_argument(
        '--input',
        default='kinases_domains_e0.01.csv',
        help='Input CSV with sequences and labels'
    )
    parser.add_argument(
        '--identity',
        type=float,
        default=0.4,
        help='CD-HIT identity threshold for clustering (default: 0.4)'
    )
    parser.add_argument(
        '--test-size',
        type=float,
        default=0.2,
        help='Test set fraction (default: 0.2)'
    )
    parser.add_argument(
        '--seed',
        type=int,
        default=42,
        help='Random seed (default: 42)'
    )
    parser.add_argument(
        '--min-class-size',
        type=int,
        default=5,
        help='Minimum samples per class (default: 5)'
    )
    parser.add_argument(
        '--output-dir',
        default='data',
        help='Output directory for splits.json'
    )
    parser.add_argument(
        '--multi-identity',
        action='store_true',
        help='Generate splits at multiple identities (70%%, 50%%, 40%%)'
    )
    
    args = parser.parse_args()
    
    print("="*80)
    print("HOMOLOGY-AWARE SPLIT GENERATION")
    print("="*80)
    print()
    print(f"Input:           {args.input}")
    
    # Determine which identities to process
    if args.multi_identity:
        identities = [0.7, 0.5, 0.4]
        print(f"Mode:            Multi-identity (70%, 50%, 40%)")
    else:
        identities = [args.identity]
        print(f"Identity:        {args.identity}")
    
    print(f"Test size:       {args.test_size}")
    print(f"Random seed:     {args.seed}")
    print(f"Min class size:  {args.min_class_size}")
    print()
    
    # Check CD-HIT
    if not check_cdhit_installed():
        print("❌ CD-HIT not installed!")
        sys.exit(1)
    
    print("✅ CD-HIT is installed")
    
    # Load data
    print(f"\nLoading {args.input}...")
    df = pd.read_csv(args.input)
    print(f"✅ Loaded {len(df):,} sequences")
    
    # Exclude "Other" category
    df_no_other = df[df['kinome_group_major'] != 'Other'].copy()
    print(f"  Excluded 'Other': {len(df_no_other):,} sequences remaining")
    
    # Process each identity threshold
    all_splits = {}
    for identity in identities:
        print("\n" + "="*80)
        print(f"Processing identity threshold: {identity:.0%}")
        print("="*80)
        
        # Run CD-HIT clustering
        df_clustered = run_cdhit_clustering(
            df_no_other, 
            identity=identity,
            output_prefix=f"cdhit_splits_{int(identity*100)}"
        )
        
        if df_clustered is None:
            print(f"❌ CD-HIT clustering failed for identity={identity}")
            continue
        
        # Create splits
        train_ids, test_ids, metadata = create_homology_aware_split(
            df_clustered,
            test_size=args.test_size,
            random_state=args.seed,
            min_class_size=args.min_class_size
        )
        
        # Update metadata with identity
        metadata['identity_threshold'] = identity
        
        # Save splits with identity suffix
        identity_suffix = f"_{int(identity*100)}" if args.multi_identity else ""
        output_file = f"{args.output_dir}/splits{identity_suffix}.json"
        
        os.makedirs(args.output_dir, exist_ok=True)
        with open(output_file, 'w') as f:
            json.dump({
                'train_ids': train_ids,
                'test_ids': test_ids,
                'metadata': metadata
            }, f, indent=2)
        
        all_splits[f"{int(identity*100)}%"] = {
            'file': output_file,
            'train': len(train_ids),
            'test': len(test_ids),
            'metadata': metadata
        }
        
        print(f"\n✅ Splits saved to: {output_file}")
        print(f"  Train: {len(train_ids):,} sequences ({len(train_ids)/(len(train_ids)+len(test_ids))*100:.1f}%)")
        print(f"  Test:  {len(test_ids):,} sequences ({len(test_ids)/(len(train_ids)+len(test_ids))*100:.1f}%)")
        print(f"  Clusters: {metadata.get('n_clusters', 'N/A')}")
        print(f"  No homology overlap (≤{identity:.0%} identity threshold)")
    
    # Update provenance
    print("\n" + "="*80)
    print("Updating provenance...")
    prov = ProvenanceTracker(output_dir=args.output_dir)
    for identity_label, split_info in all_splits.items():
        prov.add_split_info(split_info['metadata'])
    prov.add_cdhit_info(thresholds={"splits": list(identities)})
    print(f"✅ Updated {prov.provenance_file}")
    
    print("\n" + "="*80)
    print("✅ SPLIT GENERATION COMPLETE!")
    print("="*80)
    print()
    print("Summary of all splits:")
    for identity_label, split_info in all_splits.items():
        print(f"\n  {identity_label} identity:")
        print(f"    File:  {split_info['file']}")
        print(f"    Train: {split_info['train']:,} sequences")
        print(f"    Test:  {split_info['test']:,} sequences")
    print()


if __name__ == "__main__":
    main()
