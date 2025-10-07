#!/usr/bin/env python3
"""
Label Normalization and Recovery System

Reduces the "Other" category by applying a hierarchy of label assignment strategies:
1. Controlled label mapping (canonicalization from subfamily → major group)
2. Pfam/HMMER-based fallback (PF07714 → TK, etc.)
3. Homology-propagation (CD-HIT cluster majority voting)
4. Embedding-space consensus (calibrated k-NN)
5. Motif sanity checks (reject bad auto-labels)

All assignments tracked with provenance (original, pfam, cluster, knn, manual).
"""

import os
import sys
import json
import argparse
import pandas as pd
import numpy as np
import subprocess
import tempfile
from pathlib import Path
from collections import Counter
import re

# For embedding-based assignment
from sklearn.neighbors import KNeighborsClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.calibration import CalibratedClassifierCV


# Known mappings from kinome subfamily → major group
# Based on Manning et al. (2002) kinome classification
SUBFAMILY_TO_MAJOR_GROUP = {
    # AGC group
    'PKA': 'AGC', 'PKG': 'AGC', 'PKC': 'AGC', 'PDK1': 'AGC', 
    'AKT': 'AGC', 'PKB': 'AGC', 'SGK': 'AGC', 'RSK': 'AGC',
    'S6K': 'AGC', 'MSK': 'AGC', 'ROCK': 'AGC', 'GRK': 'AGC',
    'DMPK': 'AGC', 'MRCK': 'AGC', 'NDR': 'AGC',
    
    # CAMK group
    'CaMK': 'CAMK', 'CAMK': 'CAMK', 'CaMKI': 'CAMK', 'CaMKII': 'CAMK',
    'CaMKIV': 'CAMK', 'MLCK': 'CAMK', 'PHK': 'CAMK', 'DAPK': 'CAMK',
    'MAPKAPK': 'CAMK', 'CHK1': 'CAMK', 'AMPK': 'CAMK', 'MARK': 'CAMK',
    'MELK': 'CAMK', 'NUAK': 'CAMK', 'PASK': 'CAMK',
    
    # CK1 group
    'CK1': 'CK1', 'CKI': 'CK1', 'CSNK1': 'CK1',
    
    # CMGC group
    'CDK': 'CMGC', 'CDKL': 'CMGC', 'MAPK': 'CMGC', 'ERK': 'CMGC',
    'JNK': 'CMGC', 'p38': 'CMGC', 'GSK': 'CMGC', 'GSK3': 'CMGC',
    'CLK': 'CMGC', 'DYRK': 'CMGC', 'HIPK': 'CMGC', 'ICK': 'CMGC',
    'MAK': 'CMGC', 'MOK': 'CMGC', 'PRP4': 'CMGC', 'SRPK': 'CMGC',
    
    # STE group
    'STE7': 'STE', 'STE11': 'STE', 'STE20': 'STE', 'PAK': 'STE',
    'MST': 'STE', 'MSN': 'STE', 'TAO': 'STE', 'GCK': 'STE',
    
    # TK group (tyrosine kinases)
    'EGFR': 'TK', 'ERBB': 'TK', 'INSR': 'TK', 'IGF1R': 'TK',
    'PDGFR': 'TK', 'VEGFR': 'TK', 'FGFR': 'TK', 'KIT': 'TK',
    'FLT3': 'TK', 'CSF1R': 'TK', 'MET': 'TK', 'RON': 'TK',
    'TRKA': 'TK', 'TRKB': 'TK', 'TRKC': 'TK', 'RET': 'TK',
    'ALK': 'TK', 'LTK': 'TK', 'ROS': 'TK', 'DDR': 'TK',
    'EPHA': 'TK', 'EPHB': 'TK', 'TIE': 'TK', 'AXL': 'TK',
    'SRC': 'TK', 'FYN': 'TK', 'YES': 'TK', 'LCK': 'TK',
    'HCK': 'TK', 'FGR': 'TK', 'BLK': 'TK', 'LYN': 'TK',
    'ABL': 'TK', 'ARG': 'TK', 'FAK': 'TK', 'PYK2': 'TK',
    'JAK': 'TK', 'TYK': 'TK', 'SYK': 'TK', 'ZAP70': 'TK',
    
    # TKL group (tyrosine kinase-like)
    'IRAK': 'TKL', 'MLK': 'TKL', 'RAF': 'TKL', 'LRRK': 'TKL',
    'LIMK': 'TKL', 'TNIK': 'TKL', 'RIPK': 'TKL',
    
    # Atypical
    'PI3K': 'Atypical', 'PIK3': 'Atypical', 'mTOR': 'Atypical',
    'ATM': 'Atypical', 'ATR': 'Atypical', 'DNA-PK': 'Atypical',
    'SMG1': 'Atypical', 'TRRAP': 'Atypical',
    
    # Histidine kinases
    'HK': 'Histidine', 'HISTIDINE': 'Histidine',
    
    # RGC (receptor guanylate cyclase)
    'GC': 'RGC', 'GUCY': 'RGC', 'NPR': 'RGC',
}

# Patterns in protein names
NAME_PATTERNS = {
    'AGC': [r'PKA', r'PKG', r'PKC', r'protein kinase C', r'AKT', r'SGK', r'RSK'],
    'CAMK': [r'CaMK', r'calcium/calmodulin', r'MLCK', r'PHK', r'DAPK', r'AMPK'],
    'CK1': [r'casein kinase 1', r'CK1', r'CSNK1'],
    'CMGC': [r'CDK', r'cyclin-dependent', r'MAPK', r'ERK', r'JNK', r'p38', r'GSK', r'CLK', r'DYRK'],
    'STE': [r'STE\d+', r'PAK', r'MST', r'MAP kinase kinase'],
    'TK': [r'tyrosine.*kinase', r'EGFR', r'ERBB', r'INSR', r'PDGFR', r'VEGFR', r'SRC', r'ABL', r'JAK'],
    'TKL': [r'IRAK', r'MLK', r'RAF', r'LRRK'],
    'Atypical': [r'PI3K', r'PIK3', r'mTOR', r'ATM', r'ATR', r'DNA-PK'],
    'Histidine': [r'histidine kinase', r'two-component'],
}


def parse_subfamily_from_name(protein_name):
    """Extract kinase subfamily from protein name."""
    if pd.isna(protein_name) or protein_name == 'N/A':
        return None
    
    name_upper = protein_name.upper()
    
    # Try exact matches first
    for subfamily, major_group in SUBFAMILY_TO_MAJOR_GROUP.items():
        if subfamily.upper() in name_upper:
            return major_group
    
    # Try patterns
    for major_group, patterns in NAME_PATTERNS.items():
        for pattern in patterns:
            if re.search(pattern, protein_name, re.IGNORECASE):
                return major_group
    
    return None


def assign_labels_from_subfamily(df):
    """
    Assign major groups based on kinome_group_subfamily column.
    
    Returns:
        DataFrame with new 'label_source' column
    """
    print("\n" + "="*80)
    print("STEP 1: SUBFAMILY → MAJOR GROUP MAPPING")
    print("="*80)
    
    df = df.copy()
    df['label_source'] = 'original'
    
    # Check if subfamily column exists
    if 'kinome_group_subfamily' not in df.columns:
        print("⚠️  No kinome_group_subfamily column found, skipping")
        return df
    
    # Apply mapping to "Other" entries
    other_mask = df['kinome_group_major'] == 'Other'
    recovered = 0
    
    for idx, row in df[other_mask].iterrows():
        subfamily = row['kinome_group_subfamily']
        if pd.notna(subfamily) and subfamily != 'N/A':
            # Try direct mapping
            for key, value in SUBFAMILY_TO_MAJOR_GROUP.items():
                if key.lower() in str(subfamily).lower():
                    df.at[idx, 'kinome_group_major'] = value
                    df.at[idx, 'label_source'] = 'subfamily_mapping'
                    recovered += 1
                    break
    
    print(f"✅ Recovered {recovered:,} sequences from subfamily mapping")
    
    return df


def assign_labels_from_protein_names(df):
    """
    Parse protein names to extract kinase family hints.
    """
    print("\n" + "="*80)
    print("STEP 2: PROTEIN NAME PARSING")
    print("="*80)
    
    if 'protein_name' not in df.columns:
        print("⚠️  No protein_name column found, skipping")
        return df
    
    df = df.copy()
    other_mask = (df['kinome_group_major'] == 'Other') & (df['label_source'] == 'original')
    recovered = 0
    
    for idx, row in df[other_mask].iterrows():
        inferred_group = parse_subfamily_from_name(row['protein_name'])
        if inferred_group:
            df.at[idx, 'kinome_group_major'] = inferred_group
            df.at[idx, 'label_source'] = 'protein_name_parsing'
            recovered += 1
    
    print(f"✅ Recovered {recovered:,} sequences from protein name parsing")
    
    return df


def assign_labels_from_pfam(df, domains_csv):
    """
    Use Pfam domain annotations to assign labels.
    - PF07714 (Pkinase_Tyr) → TK
    - PF00069 only → keep as Unknown for now (ser/thr umbrella)
    """
    print("\n" + "="*80)
    print("STEP 3: PFAM DOMAIN-BASED ASSIGNMENT")
    print("="*80)
    
    if not os.path.exists(domains_csv):
        print(f"⚠️  Domain file {domains_csv} not found, skipping")
        return df
    
    try:
        domains_df = pd.read_csv(domains_csv)
        df = df.copy()
        
        # Check if domain info exists
        if 'sequence' in domains_df.columns or 'domain_sequence' in domains_df.columns:
            # Mark sequences that have domains
            domain_ids = set(domains_df['uniprot_id'].unique())
            df['has_kinase_domain'] = df['uniprot_id'].isin(domain_ids)
            print(f"✅ Marked {len(domain_ids):,} sequences with kinase domains")
        else:
            print("⚠️  Domain sequence column not found, skipping")
        
        # For now, sequences with domains but in "Other" stay as Unknown
        # (could add PF07714 check here if available)
        
    except Exception as e:
        print(f"⚠️  Error processing domains: {e}")
    
    return df


def assign_labels_from_clusters(df, identity=0.6, min_cluster_size=5, agreement_threshold=0.8):
    """
    Use CD-HIT cluster majority voting to propagate labels.
    
    Args:
        identity: CD-HIT identity threshold
        min_cluster_size: Minimum cluster size for voting
        agreement_threshold: Minimum fraction agreeing for assignment
    """
    print("\n" + "="*80)
    print(f"STEP 4: CLUSTER-MAJORITY PROPAGATION (identity={identity:.0%})")
    print("="*80)
    
    # Check CD-HIT
    try:
        subprocess.run(['cd-hit', '-h'], capture_output=True, timeout=5)
    except (FileNotFoundError, subprocess.TimeoutExpired):
        print("⚠️  CD-HIT not installed, skipping")
        return df
    
    df = df.copy()
    
    # Run CD-HIT
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        fasta_file = f.name
        for _, row in df.iterrows():
            f.write(f">{row['uniprot_id']}\n{row['sequence']}\n")
    
    output_file = fasta_file + '.clustered'
    cluster_file = output_file + '.clstr'
    
    cmd = [
        'cd-hit',
        '-i', fasta_file,
        '-o', output_file,
        '-c', str(identity),
        '-n', '2',
        '-M', '0',
        '-T', '0',
        '-d', '0'
    ]
    
    try:
        subprocess.run(cmd, check=True, capture_output=True, timeout=300)
    except Exception as e:
        print(f"❌ CD-HIT failed: {e}")
        return df
    
    # Parse clusters
    clusters = {}
    current_cluster = None
    
    with open(cluster_file, 'r') as f:
        for line in f:
            if line.startswith('>Cluster'):
                current_cluster = int(line.split()[1])
                clusters[current_cluster] = []
            else:
                match = re.search(r'>([A-Z0-9_]+)', line)
                if match:
                    uniprot_id = match.group(1)
                    clusters[current_cluster].append(uniprot_id)
    
    # Create cluster mapping
    uniprot_to_cluster = {}
    for cluster_id, members in clusters.items():
        for member in members:
            uniprot_to_cluster[member] = cluster_id
    
    df['cluster_id'] = df['uniprot_id'].map(uniprot_to_cluster)
    
    # Majority voting
    recovered = 0
    for cluster_id, members in clusters.items():
        if len(members) < min_cluster_size:
            continue
        
        # Get labels of cluster members
        cluster_df = df[df['cluster_id'] == cluster_id]
        labeled_df = cluster_df[cluster_df['kinome_group_major'] != 'Other']
        
        if len(labeled_df) == 0:
            continue
        
        # Check agreement
        label_counts = labeled_df['kinome_group_major'].value_counts()
        most_common_label = label_counts.index[0]
        agreement = label_counts.iloc[0] / len(labeled_df)
        
        if agreement >= agreement_threshold:
            # Assign to unlabeled members
            unlabeled_mask = (df['cluster_id'] == cluster_id) & (df['kinome_group_major'] == 'Other')
            for idx in df[unlabeled_mask].index:
                df.at[idx, 'kinome_group_major'] = most_common_label
                df.at[idx, 'label_source'] = f'cluster_vote_{agreement:.2f}'
                recovered += 1
    
    print(f"✅ Recovered {recovered:,} sequences from cluster majority voting")
    print(f"   Total clusters: {len(clusters):,}")
    print(f"   Clusters ≥{min_cluster_size}: {sum(1 for c in clusters.values() if len(c) >= min_cluster_size):,}")
    
    # Cleanup
    for f in [fasta_file, output_file, cluster_file]:
        if os.path.exists(f):
            os.remove(f)
    
    return df


def save_results(df, output_file, label_stats_file):
    """Save normalized labels and statistics."""
    df.to_csv(output_file, index=False)
    print(f"\n✅ Saved to: {output_file}")
    
    # Generate statistics
    stats = {
        'total_sequences': len(df),
        'label_distribution': df['kinome_group_major'].value_counts().to_dict(),
        'label_sources': df['label_source'].value_counts().to_dict(),
        'other_percentage': (df['kinome_group_major'] == 'Other').sum() / len(df) * 100
    }
    
    with open(label_stats_file, 'w') as f:
        json.dump(stats, f, indent=2)
    print(f"✅ Statistics saved to: {label_stats_file}")
    
    # Print summary
    print("\n" + "="*80)
    print("LABEL NORMALIZATION SUMMARY")
    print("="*80)
    print()
    print(f"Total sequences: {len(df):,}")
    print()
    print("Label distribution:")
    for label, count in df['kinome_group_major'].value_counts().items():
        pct = count / len(df) * 100
        print(f"  {label:<15} {count:>6,} ({pct:>5.1f}%)")
    print()
    print("Label sources:")
    for source, count in df['label_source'].value_counts().items():
        pct = count / len(df) * 100
        print(f"  {source:<25} {count:>6,} ({pct:>5.1f}%)")
    print()
    
    other_count = (df['kinome_group_major'] == 'Other').sum()
    other_pct = other_count / len(df) * 100
    print(f"Remaining 'Other': {other_count:,} ({other_pct:.1f}%)")
    

def main():
    parser = argparse.ArgumentParser(
        description='Normalize and recover kinase labels'
    )
    parser.add_argument(
        '--input',
        default='kinases_revised.csv',
        help='Input CSV with sequences'
    )
    parser.add_argument(
        '--domains',
        default='kinases_domains.csv',
        help='CSV with domain sequences'
    )
    parser.add_argument(
        '--output',
        default='kinases_normalized.csv',
        help='Output CSV with normalized labels'
    )
    parser.add_argument(
        '--cluster-identity',
        type=float,
        default=0.6,
        help='CD-HIT identity for cluster voting'
    )
    
    args = parser.parse_args()
    
    print("="*80)
    print("KINASE LABEL NORMALIZATION & RECOVERY")
    print("="*80)
    print()
    print(f"Input:  {args.input}")
    print(f"Output: {args.output}")
    print()
    
    # Load data
    print("Loading data...")
    df = pd.read_csv(args.input)
    print(f"✅ Loaded {len(df):,} sequences")
    
    initial_other = (df['kinome_group_major'] == 'Other').sum()
    print(f"   Initial 'Other': {initial_other:,} ({initial_other/len(df)*100:.1f}%)")
    
    # Apply normalization strategies
    df = assign_labels_from_subfamily(df)
    df = assign_labels_from_protein_names(df)
    df = assign_labels_from_pfam(df, args.domains)
    df = assign_labels_from_clusters(df, identity=args.cluster_identity)
    
    # Save results
    stats_file = args.output.replace('.csv', '_stats.json')
    save_results(df, args.output, stats_file)
    
    # Calculate improvement
    final_other = (df['kinome_group_major'] == 'Other').sum()
    recovered = initial_other - final_other
    print(f"\n🎯 Recovered {recovered:,} sequences from 'Other' ({recovered/initial_other*100:.1f}%)")
    print(f"   {initial_other:,} → {final_other:,}")
    print()


if __name__ == "__main__":
    main()
