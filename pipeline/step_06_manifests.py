#!/usr/bin/env python3
"""
Step 6: Create Dataset Manifests

This script creates manifest files that define dataset membership.
Manifests are the SINGLE SOURCE OF TRUTH for which sequences belong to each dataset.

Usage:
    python pipeline/step_06_manifests.py --run-dir runs/2025-01-01_000000/

Outputs:
    - data/manifests/domain_E001.txt
    - data/manifests/supervised_eligible.txt
    - data/manifests/manifest_report.json
"""

import argparse
import json
import pandas as pd
from pathlib import Path
from datetime import datetime
from collections import Counter

# Project root
PROJECT_ROOT = Path(__file__).parent.parent


def main():
    parser = argparse.ArgumentParser(description="Create dataset manifests")
    parser.add_argument("--run-dir", type=str, required=True, help="Run directory path")
    args = parser.parse_args()
    
    run_dir = Path(args.run_dir)
    manifests_dir = run_dir / "data" / "manifests"
    manifests_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 60)
    print("Step 6: Create Dataset Manifests")
    print("=" * 60)
    
    # Load source data
    domain_coords = PROJECT_ROOT / "data" / "domains" / "domain_coords_E001.tsv"
    labels_file = PROJECT_ROOT / "data" / "processed" / "labels.csv"
    
    print(f"\nLoading domain coordinates: {domain_coords}")
    coords_df = pd.read_csv(domain_coords, sep='\t')
    domain_ids = set(coords_df['uniprot_id'].tolist())
    print(f"  Found {len(domain_ids)} unique domains")
    
    print(f"\nLoading labels: {labels_file}")
    labels_df = pd.read_csv(labels_file)
    id_to_label = dict(zip(labels_df['uniprot_id'], labels_df['label_used_for_experiments']))
    
    # Get labels for domain IDs
    domain_labels = {uid: id_to_label.get(uid, 'Unknown') for uid in domain_ids}
    
    # Count classes
    class_counts = Counter(domain_labels.values())
    print(f"\nClass distribution (all domains):")
    for cls, count in sorted(class_counts.items(), key=lambda x: -x[1]):
        print(f"  {cls:<15} {count}")
    
    # Create domain_E001 manifest (exclude 'Other')
    domain_excl_other = sorted([uid for uid in domain_ids if domain_labels.get(uid) != 'Other'])
    domain_e001_file = manifests_dir / "domain_E001.txt"
    with open(domain_e001_file, 'w') as f:
        for uid in domain_excl_other:
            f.write(f"{uid}\n")
    print(f"\n✓ Created domain_E001.txt: {len(domain_excl_other)} sequences")
    
    # Create supervised_eligible manifest (exclude 'Other', 'Histidine', 'RGC')
    excluded_classes = {'Other', 'Histidine', 'RGC'}
    supervised = sorted([uid for uid in domain_ids if domain_labels.get(uid) not in excluded_classes])
    supervised_file = manifests_dir / "supervised_eligible.txt"
    with open(supervised_file, 'w') as f:
        for uid in supervised:
            f.write(f"{uid}\n")
    print(f"✓ Created supervised_eligible.txt: {len(supervised)} sequences")
    
    # Count supervised classes
    supervised_classes = Counter([domain_labels[uid] for uid in supervised])
    print(f"\nSupervised-eligible class distribution:")
    for cls, count in sorted(supervised_classes.items(), key=lambda x: -x[1]):
        print(f"  {cls:<15} {count}")
    
    # Assertions
    assert len(domain_excl_other) > 0, "domain_E001 manifest is empty!"
    assert len(supervised) > 0, "supervised_eligible manifest is empty!"
    assert set(supervised).issubset(set(domain_excl_other)), \
        "supervised_eligible must be subset of domain_E001"
    
    # Create report
    report = {
        "step": 6,
        "name": "Dataset Manifests",
        "timestamp": datetime.now().isoformat(),
        "source_files": {
            "domain_coords": str(domain_coords),
            "labels": str(labels_file)
        },
        "datasets": {
            "domain_E001": {
                "file": str(domain_e001_file),
                "n_sequences": len(domain_excl_other),
                "n_classes": len(set(domain_labels[uid] for uid in domain_excl_other)),
                "excluded_classes": ["Other"],
                "per_class_counts": dict(Counter([domain_labels[uid] for uid in domain_excl_other]))
            },
            "supervised_eligible": {
                "file": str(supervised_file),
                "n_sequences": len(supervised),
                "n_classes": len(supervised_classes),
                "excluded_classes": list(excluded_classes),
                "per_class_counts": dict(supervised_classes)
            }
        },
        "assertions": {
            "domain_E001_not_empty": len(domain_excl_other) > 0,
            "supervised_eligible_not_empty": len(supervised) > 0,
            "supervised_subset_of_domain": set(supervised).issubset(set(domain_excl_other))
        }
    }
    
    report_file = manifests_dir / "manifest_report.json"
    with open(report_file, 'w') as f:
        json.dump(report, f, indent=2)
    print(f"\n✓ Created manifest_report.json")
    
    print("\n" + "=" * 60)
    print("Step 6 COMPLETE")
    print("=" * 60)
    print(f"\nManifests created in: {manifests_dir}")
    print(f"  domain_E001: {len(domain_excl_other)} sequences ({len(set(domain_labels[uid] for uid in domain_excl_other))} classes)")
    print(f"  supervised_eligible: {len(supervised)} sequences ({len(supervised_classes)} classes)")


if __name__ == "__main__":
    main()

