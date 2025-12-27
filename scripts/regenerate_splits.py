#!/usr/bin/env python3
"""
Regenerate homology-aware splits for the corrected supervised_eligible manifest.

This script:
1. Creates a FASTA file from the 1,392 supervised-eligible sequences
2. Runs CD-HIT at 40%, 50%, 70% identity thresholds
3. Creates train/test splits ensuring no cluster spans both sets
4. Saves new split files
"""

import subprocess
import json
import random
import pandas as pd
from pathlib import Path
from collections import defaultdict
from datetime import datetime

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
                # Extract sequence ID from line like: "0	256aa, >A0A075F7E9... *"
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
    
    # Calculate target test size
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
    print("="*60)
    print("Regenerating Homology-aware Splits")
    print("="*60)
    
    # Paths
    manifest_file = Path("data/manifests/supervised_eligible.txt")
    domains_fasta = Path("data/domains/domains_E001.fasta")
    output_dir = Path("data/splits")
    temp_dir = Path("data/temp_cdhit")
    temp_dir.mkdir(exist_ok=True)
    
    # Load manifest IDs
    with open(manifest_file, 'r') as f:
        manifest_ids = set(line.strip() for line in f if line.strip())
    print(f"\nLoaded {len(manifest_ids)} supervised-eligible IDs")
    
    # Load sequences from domain FASTA
    from Bio import SeqIO
    id_to_seq = {}
    for record in SeqIO.parse(domains_fasta, "fasta"):
        # Extract UniProt ID from header (format: Q7TPS0|PF00069|446-703)
        uniprot_id = record.id.split('|')[0]
        id_to_seq[uniprot_id] = str(record.seq)
    
    print(f"Loaded {len(id_to_seq)} sequences from domain FASTA")
    
    # Filter to manifest IDs
    filtered_seqs = {uid: seq for uid, seq in id_to_seq.items() if uid in manifest_ids}
    print(f"Found {len(filtered_seqs)} sequences matching manifest")
    
    # Create FASTA file for CD-HIT
    fasta_file = temp_dir / "supervised_eligible.fasta"
    with open(fasta_file, 'w') as f:
        for uid, seq in filtered_seqs.items():
            f.write(f">{uid}\n")
            f.write(f"{seq}\n")
    print(f"Created FASTA: {fasta_file}")
    
    # Load labels for class distribution
    labels_df = pd.read_csv("data/processed/labels.csv")
    id_to_label = dict(zip(labels_df['uniprot_id'], labels_df['label_used_for_experiments']))
    
    # Report structure
    report = {
        "step": 10,
        "name": "Homology-aware Supervised Splits (Regenerated)",
        "timestamp": datetime.now().isoformat(),
        "method": {
            "algorithm": "CD-HIT clustering then cluster-group split",
            "description": "Sequences clustered at identity threshold. Clusters split into train/test.",
            "guarantee": "No cluster spans train/test - prevents homology leakage",
            "test_ratio": TEST_RATIO,
            "random_state": RANDOM_STATE
        },
        "supervised_eligible_n": len(manifest_ids),
        "splits": {}
    }
    
    # Run CD-HIT at each threshold
    thresholds = [70, 50, 40]
    word_sizes = {70: 5, 50: 3, 40: 2}  # CD-HIT word size requirements
    
    for threshold in thresholds:
        print(f"\n{'='*60}")
        print(f"Processing {threshold}% identity threshold")
        print(f"{'='*60}")
        
        identity = threshold / 100.0
        word_size = word_sizes[threshold]
        
        output_prefix = temp_dir / f"cdhit_{threshold}"
        clstr_file = Path(f"{output_prefix}.clstr")
        
        # Run CD-HIT
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
            print(f"  CD-HIT completed successfully")
        except subprocess.CalledProcessError as e:
            print(f"  CD-HIT failed: {e.stderr}")
            continue
        except FileNotFoundError:
            print("  CD-HIT not found. Install with: conda install -c bioconda cd-hit")
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
            print(f"  ERROR: {len(overlap)} overlapping IDs!")
        else:
            print(f"  ✓ Train/test are disjoint")
        
        # Get class distributions
        train_labels = [id_to_label.get(uid) for uid in train_ids if uid in id_to_label]
        test_labels = [id_to_label.get(uid) for uid in test_ids if uid in id_to_label]
        
        from collections import Counter
        train_class_dist = dict(Counter(train_labels))
        test_class_dist = dict(Counter(test_labels))
        
        # Save split files
        train_file = output_dir / f"split{threshold}_train.txt"
        test_file = output_dir / f"split{threshold}_test.txt"
        
        with open(train_file, 'w') as f:
            for uid in sorted(train_ids):
                f.write(f"{uid}\n")
        
        with open(test_file, 'w') as f:
            for uid in sorted(test_ids):
                f.write(f"{uid}\n")
        
        print(f"  Saved: {train_file}")
        print(f"  Saved: {test_file}")
        
        # Add to report
        report["splits"][f"split{threshold}"] = {
            "identity_threshold": threshold,
            "n_clusters": n_clusters,
            "n_train": len(train_ids),
            "n_test": len(test_ids),
            "n_total": len(train_ids) + len(test_ids),
            "train_test_ratio": round(len(train_ids) / (len(train_ids) + len(test_ids)), 3),
            "train_class_distribution": train_class_dist,
            "test_class_distribution": test_class_dist,
            "is_disjoint": len(overlap) == 0
        }
        
        print(f"\n  Class distribution (train):")
        for cls, count in sorted(train_class_dist.items(), key=lambda x: -x[1]):
            print(f"    {cls:15} {count}")
    
    # Save report
    report_file = output_dir / "splits_report.json"
    with open(report_file, 'w') as f:
        json.dump(report, f, indent=2)
    print(f"\nReport saved: {report_file}")
    
    print(f"\n{'='*60}")
    print("SPLITS REGENERATION COMPLETE")
    print(f"{'='*60}")
    
    print("\n                          Summary                           ")
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

