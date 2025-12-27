#!/usr/bin/env python3
"""
Step 10: Create Homology-aware Train/Test Splits

This script creates train/test splits that ensure no homologous sequences
span both sets. Uses CD-HIT clustering at different identity thresholds.

Usage:
    python pipeline/step_10_splits.py --run-dir runs/2025-01-01_000000/

Outputs:
    - data/splits/split40_train.txt, split40_test.txt
    - data/splits/split50_train.txt, split50_test.txt
    - data/splits/split70_train.txt, split70_test.txt
    - data/splits/splits_report.json
"""

import argparse
import json
import subprocess
import random
from pathlib import Path
from datetime import datetime
from collections import defaultdict, Counter
from Bio import SeqIO

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from pipeline.membership import load_manifest, assert_split_integrity

# Project root
PROJECT_ROOT = Path(__file__).parent.parent

# Configuration
RANDOM_STATE = 42
TEST_RATIO = 0.2


def parse_cdhit_clusters(clstr_file):
    """Parse CD-HIT cluster file to get cluster assignments."""
    clusters = defaultdict(list)
    current_cluster = None
    
    with open(clstr_file, 'r') as f:
        for line in f:
            if line.startswith('>Cluster'):
                current_cluster = int(line.strip().split()[1])
            else:
                if '>' in line:
                    start = line.index('>') + 1
                    end = line.index('...')
                    seq_id = line[start:end]
                    clusters[current_cluster].append(seq_id)
    
    return clusters


def create_splits_from_clusters(clusters, test_ratio, random_state):
    """Split clusters into train/test, keeping all sequences in a cluster together."""
    random.seed(random_state)
    
    cluster_ids = list(clusters.keys())
    random.shuffle(cluster_ids)
    
    total_seqs = sum(len(seqs) for seqs in clusters.values())
    target_test = int(total_seqs * test_ratio)
    
    test_ids = []
    train_ids = []
    test_size = 0
    
    for cluster_id in cluster_ids:
        seqs = clusters[cluster_id]
        if test_size < target_test:
            test_ids.extend(seqs)
            test_size += len(seqs)
        else:
            train_ids.extend(seqs)
    
    return train_ids, test_ids


def main():
    parser = argparse.ArgumentParser(description="Create homology-aware splits")
    parser.add_argument("--run-dir", type=str, required=True, help="Run directory path")
    args = parser.parse_args()
    
    run_dir = Path(args.run_dir)
    splits_dir = run_dir / "data" / "splits"
    splits_dir.mkdir(parents=True, exist_ok=True)
    temp_dir = run_dir / "data" / "temp_cdhit"
    temp_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 60)
    print("Step 10: Create Homology-aware Splits")
    print("=" * 60)
    
    # Load supervised_eligible manifest
    manifest = load_manifest("supervised_eligible", run_dir)
    print(f"\nLoaded manifest: {len(manifest)} supervised-eligible sequences")
    
    # Load domain FASTA
    domains_fasta = PROJECT_ROOT / "data" / "domains" / "domains_E001.fasta"
    id_to_seq = {}
    for record in SeqIO.parse(domains_fasta, "fasta"):
        uniprot_id = record.id.split('|')[0]
        id_to_seq[uniprot_id] = str(record.seq)
    
    # Filter to manifest IDs
    filtered_seqs = {uid: seq for uid, seq in id_to_seq.items() if uid in manifest}
    print(f"Found {len(filtered_seqs)} sequences matching manifest")
    
    # Create FASTA file for CD-HIT
    fasta_file = temp_dir / "supervised_eligible.fasta"
    with open(fasta_file, 'w') as f:
        for uid, seq in filtered_seqs.items():
            f.write(f">{uid}\n{seq}\n")
    
    # Load labels for class distribution
    labels_file = PROJECT_ROOT / "data" / "processed" / "labels.csv"
    import pandas as pd
    labels_df = pd.read_csv(labels_file)
    id_to_label = dict(zip(labels_df['uniprot_id'], labels_df['label_used_for_experiments']))
    
    # Report structure
    report = {
        "step": 10,
        "name": "Homology-aware Splits",
        "timestamp": datetime.now().isoformat(),
        "method": {
            "algorithm": "CD-HIT clustering + cluster-group split",
            "guarantee": "No cluster spans train/test",
            "test_ratio": TEST_RATIO,
            "random_state": RANDOM_STATE
        },
        "supervised_eligible_n": len(manifest),
        "splits": {}
    }
    
    # Run CD-HIT at each threshold
    thresholds = [70, 50, 40]
    word_sizes = {70: 5, 50: 3, 40: 2}
    
    for threshold in thresholds:
        print(f"\n{'=' * 60}")
        print(f"Processing {threshold}% identity threshold")
        print(f"{'=' * 60}")
        
        identity = threshold / 100.0
        word_size = word_sizes[threshold]
        
        output_prefix = temp_dir / f"cdhit_{threshold}"
        clstr_file = Path(f"{output_prefix}.clstr")
        
        cmd = [
            "cd-hit",
            "-i", str(fasta_file),
            "-o", str(output_prefix),
            "-c", str(identity),
            "-n", str(word_size),
            "-M", "0",
            "-T", "0"
        ]
        
        print(f"  Running CD-HIT with -c {identity}...")
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            print(f"  ✓ CD-HIT completed")
        except subprocess.CalledProcessError as e:
            print(f"  ✗ CD-HIT failed: {e.stderr}")
            continue
        except FileNotFoundError:
            print("  ✗ CD-HIT not found. Install with: conda install -c bioconda cd-hit")
            continue
        
        # Parse clusters
        clusters = parse_cdhit_clusters(clstr_file)
        n_clusters = len(clusters)
        print(f"  Found {n_clusters} clusters")
        
        # Create splits
        train_ids, test_ids = create_splits_from_clusters(clusters, TEST_RATIO, RANDOM_STATE)
        print(f"  Split: {len(train_ids)} train, {len(test_ids)} test")
        
        # Verify disjoint
        overlap = set(train_ids) & set(test_ids)
        if overlap:
            print(f"  ✗ ERROR: {len(overlap)} overlapping IDs!")
        else:
            print(f"  ✓ Train/test are disjoint")
        
        # Get class distributions
        train_labels = [id_to_label.get(uid) for uid in train_ids if uid in id_to_label]
        test_labels = [id_to_label.get(uid) for uid in test_ids if uid in id_to_label]
        
        train_class_dist = dict(Counter(train_labels))
        test_class_dist = dict(Counter(test_labels))
        
        # Save split files
        train_file = splits_dir / f"split{threshold}_train.txt"
        test_file = splits_dir / f"split{threshold}_test.txt"
        
        with open(train_file, 'w') as f:
            for uid in sorted(train_ids):
                f.write(f"{uid}\n")
        
        with open(test_file, 'w') as f:
            for uid in sorted(test_ids):
                f.write(f"{uid}\n")
        
        print(f"  ✓ Saved: {train_file.name}")
        print(f"  ✓ Saved: {test_file.name}")
        
        # Add to report
        report["splits"][f"split{threshold}"] = {
            "identity_threshold": threshold,
            "n_clusters": n_clusters,
            "n_train": len(train_ids),
            "n_test": len(test_ids),
            "n_total": len(train_ids) + len(test_ids),
            "train_class_distribution": train_class_dist,
            "test_class_distribution": test_class_dist,
            "is_disjoint": len(overlap) == 0
        }
    
    # Save report
    report_file = splits_dir / "splits_report.json"
    with open(report_file, 'w') as f:
        json.dump(report, f, indent=2)
    print(f"\n✓ Saved: {report_file}")
    
    # Validate splits using membership module
    print("\n" + "-" * 60)
    print("Validating split integrity...")
    for threshold in thresholds:
        try:
            assert_split_integrity(run_dir, f"split{threshold}")
        except AssertionError as e:
            print(f"  ✗ split{threshold}: {e}")
        except FileNotFoundError:
            print(f"  ⚠ split{threshold}: files not found")
    
    print("\n" + "=" * 60)
    print("Step 10 COMPLETE")
    print("=" * 60)
    print("\n                          Summary")
    print("-" * 60)
    print("Threshold    Clusters     Train      Test       Total")
    print("-" * 60)
    for threshold in thresholds:
        if f"split{threshold}" in report["splits"]:
            s = report["splits"][f"split{threshold}"]
            print(f"{threshold}%         {s['n_clusters']:<12} {s['n_train']:<10} {s['n_test']:<10} {s['n_total']}")
    print("-" * 60)


if __name__ == "__main__":
    main()

