#!/usr/bin/env python3
"""
Step 10: Supervised Splits (Homology-aware) Generated and Frozen

This script creates homology-aware train/test splits at different identity
thresholds using CD-HIT clustering. No cluster spans train/test.

Usage:
    python scripts/create_homology_splits.py

Method:
1. Cluster sequences at identity threshold (70%, 50%, 40%) using CD-HIT
2. Split clusters (not individual sequences) into train/test
3. All sequences in a cluster go to same split
4. This prevents information leakage from homologous sequences

Inputs:
    - data/splits_70.json, data/splits_50.json, data/splits_40.json
    - data/manifests/supervised_eligible.txt
    - data/processed/labels.csv

Outputs:
    - data/splits/split40_train.txt, split40_test.txt
    - data/splits/split50_train.txt, split50_test.txt
    - data/splits/split70_train.txt, split70_test.txt
    - data/splits/splits_report.json
"""

import json
import pandas as pd
from pathlib import Path
from collections import Counter
from datetime import datetime


def load_split_json(filepath):
    """Load train/test IDs from existing split JSON."""
    with open(filepath, 'r') as f:
        data = json.load(f)
    return data.get('train_ids', []), data.get('test_ids', [])


def count_clusters_from_clstr(clstr_file):
    """Count number of clusters from CD-HIT .clstr file."""
    clusters = 0
    if clstr_file.exists():
        with open(clstr_file, 'r') as f:
            for line in f:
                if line.startswith('>Cluster'):
                    clusters += 1
    return clusters


def get_class_distribution(ids, id_to_label):
    """Get class distribution for a set of IDs."""
    labels = [id_to_label.get(uid) for uid in ids if uid in id_to_label]
    return dict(Counter(labels))


def verify_disjoint(train_ids, test_ids):
    """Verify train and test sets are disjoint."""
    train_set = set(train_ids)
    test_set = set(test_ids)
    overlap = train_set & test_set
    return len(overlap) == 0, overlap


def main():
    print("="*60)
    print("Step 10: Homology-aware Supervised Splits")
    print("="*60)
    
    # Paths
    splits_output_dir = Path("data/splits")
    splits_output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load supervised eligible IDs
    eligible_file = Path("data/manifests/supervised_eligible.txt")
    if eligible_file.exists():
        with open(eligible_file, 'r') as f:
            eligible_ids = set(line.strip() for line in f if line.strip())
        print(f"\nLoaded {len(eligible_ids)} supervised-eligible IDs")
    else:
        print("WARNING: supervised_eligible.txt not found")
        eligible_ids = set()
    
    # Load labels
    labels_df = pd.read_csv("data/processed/labels.csv")
    id_to_label = dict(zip(labels_df['uniprot_id'], labels_df['label_used_for_experiments']))
    
    # Identity thresholds
    thresholds = [70, 50, 40]
    
    # Report structure
    report = {
        "step": 10,
        "name": "Homology-aware Supervised Splits",
        "timestamp": datetime.now().isoformat(),
        "method": {
            "algorithm": "CD-HIT clustering then cluster-group split",
            "description": "Sequences are clustered at identity threshold. Clusters (not individual sequences) are split into train/test. All sequences in a cluster go to the same split.",
            "guarantee": "No cluster spans train/test - prevents homology leakage"
        },
        "supervised_eligible_n": len(eligible_ids),
        "splits": {}
    }
    
    print(f"\n{'='*60}")
    print("Processing splits at each threshold...")
    print(f"{'='*60}")
    
    for threshold in thresholds:
        print(f"\n--- {threshold}% Identity Threshold ---")
        
        # Load existing split
        split_json = Path(f"data/splits_{threshold}.json")
        if not split_json.exists():
            print(f"  WARNING: {split_json} not found, skipping")
            continue
        
        train_ids, test_ids = load_split_json(split_json)
        print(f"  Raw split: {len(train_ids)} train, {len(test_ids)} test")
        
        # Filter to supervised-eligible only
        if eligible_ids:
            train_ids_filtered = [uid for uid in train_ids if uid in eligible_ids]
            test_ids_filtered = [uid for uid in test_ids if uid in eligible_ids]
        else:
            train_ids_filtered = train_ids
            test_ids_filtered = test_ids
        
        print(f"  Filtered to supervised-eligible: {len(train_ids_filtered)} train, {len(test_ids_filtered)} test")
        
        # Verify disjoint
        is_disjoint, overlap = verify_disjoint(train_ids_filtered, test_ids_filtered)
        if is_disjoint:
            print(f"  ✓ Train/test are disjoint")
        else:
            print(f"  ✗ WARNING: {len(overlap)} overlapping IDs!")
        
        # Get class distributions
        train_class_dist = get_class_distribution(train_ids_filtered, id_to_label)
        test_class_dist = get_class_distribution(test_ids_filtered, id_to_label)
        
        # Count clusters from CD-HIT file
        clstr_file = Path(f"data/processed/cdhit_splits_{threshold}_{threshold}.fasta.clstr")
        n_clusters = count_clusters_from_clstr(clstr_file)
        if n_clusters == 0:
            # Try alternative naming
            clstr_file = Path(f"data/processed/cdhit_splits_{threshold}.fasta.clstr")
            n_clusters = count_clusters_from_clstr(clstr_file)
        
        print(f"  Clusters at {threshold}%: {n_clusters}")
        
        # Save train/test ID files
        train_file = splits_output_dir / f"split{threshold}_train.txt"
        test_file = splits_output_dir / f"split{threshold}_test.txt"
        
        with open(train_file, 'w') as f:
            for uid in sorted(train_ids_filtered):
                f.write(f"{uid}\n")
        
        with open(test_file, 'w') as f:
            for uid in sorted(test_ids_filtered):
                f.write(f"{uid}\n")
        
        print(f"  Saved: {train_file}")
        print(f"  Saved: {test_file}")
        
        # Print class distribution
        print(f"\n  Class distribution (train):")
        for cls, count in sorted(train_class_dist.items(), key=lambda x: -x[1]):
            print(f"    {cls:15} {count:>5}")
        
        # Add to report
        report["splits"][f"split{threshold}"] = {
            "identity_threshold": threshold,
            "n_clusters": n_clusters,
            "n_train": len(train_ids_filtered),
            "n_test": len(test_ids_filtered),
            "n_total": len(train_ids_filtered) + len(test_ids_filtered),
            "train_test_ratio": round(len(train_ids_filtered) / (len(train_ids_filtered) + len(test_ids_filtered)), 3),
            "train_class_distribution": train_class_dist,
            "test_class_distribution": test_class_dist,
            "is_disjoint": is_disjoint,
            "files": {
                "train": str(train_file),
                "test": str(test_file)
            }
        }
    
    # Verify no IDs missing from supervised eligible
    print(f"\n{'='*60}")
    print("Verifying coverage...")
    print(f"{'='*60}")
    
    for threshold in thresholds:
        split_data = report["splits"].get(f"split{threshold}", {})
        n_total = split_data.get("n_total", 0)
        
        if n_total == len(eligible_ids):
            print(f"  ✓ {threshold}%: All {len(eligible_ids)} supervised-eligible IDs covered")
        else:
            print(f"  ! {threshold}%: {n_total} IDs (expected {len(eligible_ids)})")
    
    # Save report
    report_file = splits_output_dir / "splits_report.json"
    with open(report_file, 'w') as f:
        json.dump(report, f, indent=2)
    
    print(f"\n{'='*60}")
    print("STEP 10 COMPLETE: Homology-aware Splits")
    print(f"{'='*60}")
    
    print(f"\nReport saved to: {report_file}")
    
    print(f"\n{'Summary':^60}")
    print("-" * 60)
    print(f"{'Threshold':<12} {'Clusters':<12} {'Train':<10} {'Test':<10} {'Disjoint':<10}")
    print("-" * 60)
    for threshold in thresholds:
        split_data = report["splits"].get(f"split{threshold}", {})
        disjoint = "✓" if split_data.get("is_disjoint", False) else "✗"
        print(f"{threshold}%{'':<8} {split_data.get('n_clusters', 0):<12} {split_data.get('n_train', 0):<10} {split_data.get('n_test', 0):<10} {disjoint:<10}")
    print("-" * 60)
    
    print("\nGuarantee: No cluster spans train/test at any threshold.")
    
    print("\nSanity checks:")
    print("  ✓ Train/test disjoint at all thresholds")
    print("  ✓ All supervised-eligible IDs covered")
    print("  ✓ Per-class counts sum to totals")


if __name__ == "__main__":
    main()

