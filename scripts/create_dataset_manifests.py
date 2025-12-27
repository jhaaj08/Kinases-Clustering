#!/usr/bin/env python3
"""
Step 6: Define Experiment Datasets

This script creates explicit dataset manifests listing UniProt IDs for each 
experiment subset. These manifests are the authoritative source for which 
sequences are included in each analysis.

Usage:
    python scripts/create_dataset_manifests.py

Inputs:
    - data/processed/labels.csv (from Step 5)
    - data/processed/kinases_domains.csv (E=0.001 domains)
    - data/processed/kinases_domains_e0.01.csv (E=0.01 domains)

Outputs:
    - data/manifests/whole_seq_excl_other.txt
    - data/manifests/domain_E0001.txt
    - data/manifests/domain_E001.txt
    - data/manifests/supervised_eligible.txt
    - data/processed/dataset_manifest_report.json
"""

import pandas as pd
import json
from pathlib import Path
from collections import Counter


def main():
    # Paths
    labels_file = Path("data/processed/labels.csv")
    domains_e0001_file = Path("data/processed/kinases_domains.csv")
    domains_e001_file = Path("data/processed/kinases_domains_e0.01.csv")
    
    manifests_dir = Path("data/manifests")
    manifests_dir.mkdir(parents=True, exist_ok=True)
    
    output_report = Path("data/processed/dataset_manifest_report.json")
    
    # Load data
    print("Loading data...")
    labels_df = pd.read_csv(labels_file)
    
    # Load domain data
    if domains_e0001_file.exists():
        domains_e0001_df = pd.read_csv(domains_e0001_file)
    else:
        print(f"Warning: {domains_e0001_file} not found, using empty dataframe")
        domains_e0001_df = pd.DataFrame({'uniprot_id': []})
    
    if domains_e001_file.exists():
        domains_e001_df = pd.read_csv(domains_e001_file)
    else:
        print(f"Warning: {domains_e001_file} not found, using empty dataframe")
        domains_e001_df = pd.DataFrame({'uniprot_id': []})
    
    # Report data
    report = {
        "step": 6,
        "name": "Dataset Manifest Creation",
        "description": "Defines explicit experiment datasets with UniProt ID lists",
        "policy": {
            "label_column": "label_used_for_experiments",
            "other_excluded": True,
            "min_samples_per_class": 5
        },
        "datasets": {},
        "table_1_counts": {}
    }
    
    # ========== Manifest 1: Whole-seq excluding Other ==========
    print("\n1. Creating whole_seq_excl_other manifest...")
    
    # Get all sequences with labels that are NOT "Other"
    whole_seq_excl_other = labels_df[
        labels_df['label_used_for_experiments'] != 'Other'
    ]['uniprot_id'].tolist()
    whole_seq_excl_other = sorted(set(whole_seq_excl_other))
    
    manifest_path = manifests_dir / "whole_seq_excl_other.txt"
    with open(manifest_path, 'w') as f:
        for uid in whole_seq_excl_other:
            f.write(f"{uid}\n")
    
    # Get per-class counts
    whole_seq_class_counts = labels_df[
        labels_df['label_used_for_experiments'] != 'Other'
    ]['label_used_for_experiments'].value_counts().to_dict()
    
    report["datasets"]["whole_seq_excl_other"] = {
        "manifest_file": str(manifest_path),
        "description": "Full-length sequences excluding 'Other' class",
        "n_sequences": len(whole_seq_excl_other),
        "n_classes": len(whole_seq_class_counts),
        "per_class_counts": {k: int(v) for k, v in whole_seq_class_counts.items()}
    }
    
    print(f"   -> {len(whole_seq_excl_other)} sequences in {len(whole_seq_class_counts)} classes")
    
    # ========== Manifest 2: Domain E=0.001 ==========
    print("\n2. Creating domain_E0001 manifest...")
    
    if len(domains_e0001_df) > 0:
        domain_e0001_ids = sorted(domains_e0001_df['uniprot_id'].unique().tolist())
        
        # Filter to get labels and exclude "Other"
        domain_e0001_with_labels = domains_e0001_df.merge(
            labels_df[['uniprot_id', 'label_used_for_experiments']], 
            on='uniprot_id', 
            how='left'
        )
        domain_e0001_excl_other = domain_e0001_with_labels[
            domain_e0001_with_labels['label_used_for_experiments'] != 'Other'
        ]['uniprot_id'].unique().tolist()
        domain_e0001_excl_other = sorted(domain_e0001_excl_other)
        
        # Get class counts for domains excluding Other
        domain_e0001_class_counts = domain_e0001_with_labels[
            domain_e0001_with_labels['label_used_for_experiments'] != 'Other'
        ].drop_duplicates(subset=['uniprot_id'])['label_used_for_experiments'].value_counts().to_dict()
    else:
        domain_e0001_ids = []
        domain_e0001_excl_other = []
        domain_e0001_class_counts = {}
    
    manifest_path = manifests_dir / "domain_E0001.txt"
    with open(manifest_path, 'w') as f:
        for uid in domain_e0001_excl_other:
            f.write(f"{uid}\n")
    
    report["datasets"]["domain_E0001"] = {
        "manifest_file": str(manifest_path),
        "description": "Domain sequences (E-value < 0.001), excluding 'Other'",
        "evalue_threshold": "0.001",
        "n_sequences": len(domain_e0001_excl_other),
        "n_classes": len(domain_e0001_class_counts),
        "per_class_counts": {k: int(v) for k, v in domain_e0001_class_counts.items()}
    }
    
    print(f"   -> {len(domain_e0001_excl_other)} sequences in {len(domain_e0001_class_counts)} classes")
    
    # ========== Manifest 3: Domain E=0.01 (main dataset) ==========
    print("\n3. Creating domain_E001 manifest (MAIN DATASET)...")
    
    if len(domains_e001_df) > 0:
        domain_e001_ids = sorted(domains_e001_df['uniprot_id'].unique().tolist())
        
        # ALWAYS use labels.csv as the authoritative source for labels
        # Do NOT use kinome_group_major from the domain file - it has messy data
        domain_e001_with_labels = domains_e001_df[['uniprot_id']].drop_duplicates().merge(
            labels_df[['uniprot_id', 'label_used_for_experiments']], 
            on='uniprot_id', 
            how='left'
        )
        
        domain_e001_excl_other = domain_e001_with_labels[
            domain_e001_with_labels['label_used_for_experiments'] != 'Other'
        ]['uniprot_id'].unique().tolist()
        domain_e001_excl_other = sorted(domain_e001_excl_other)
        
        # Get class counts for domains excluding Other
        domain_e001_class_counts = domain_e001_with_labels[
            domain_e001_with_labels['label_used_for_experiments'] != 'Other'
        ].drop_duplicates(subset=['uniprot_id'])['label_used_for_experiments'].value_counts().to_dict()
    else:
        domain_e001_ids = []
        domain_e001_excl_other = []
        domain_e001_class_counts = {}
    
    manifest_path = manifests_dir / "domain_E001.txt"
    with open(manifest_path, 'w') as f:
        for uid in domain_e001_excl_other:
            f.write(f"{uid}\n")
    
    report["datasets"]["domain_E001"] = {
        "manifest_file": str(manifest_path),
        "description": "Domain sequences (E-value < 0.01), excluding 'Other' - MAIN ANALYSIS DATASET",
        "evalue_threshold": "0.01",
        "n_sequences": len(domain_e001_excl_other),
        "n_classes": len(domain_e001_class_counts),
        "per_class_counts": {k: int(v) for k, v in domain_e001_class_counts.items()},
        "is_main_dataset": True
    }
    
    print(f"   -> {len(domain_e001_excl_other)} sequences in {len(domain_e001_class_counts)} classes")
    
    # ========== Manifest 4: Supervised-eligible (min 5 per class) ==========
    print("\n4. Creating supervised_eligible manifest...")
    
    # Use domain E=0.01 as the base for supervised learning
    # Apply minimum sample threshold (n >= 5 per class)
    # ALSO exclude Histidine and RGC for biological/methodological reasons:
    # - Histidine kinases use different catalytic mechanism (HisKA + HATPase)
    # - RGC are receptor guanylate cyclases, not true kinases
    MIN_SAMPLES = 5
    EXCLUDED_CLASSES_BIOLOGICAL = ['Histidine', 'RGC']
    
    if len(domain_e001_class_counts) > 0:
        # Find classes with enough samples AND not biologically excluded
        eligible_classes = [
            cls for cls, count in domain_e001_class_counts.items() 
            if count >= MIN_SAMPLES and cls not in EXCLUDED_CLASSES_BIOLOGICAL
        ]
        
        # Get sequences in eligible classes
        supervised_eligible = domain_e001_with_labels[
            (domain_e001_with_labels['label_used_for_experiments'] != 'Other') &
            (domain_e001_with_labels['label_used_for_experiments'].isin(eligible_classes))
        ]['uniprot_id'].unique().tolist()
        supervised_eligible = sorted(supervised_eligible)
        
        # Get class counts for supervised-eligible
        supervised_class_counts = domain_e001_with_labels[
            (domain_e001_with_labels['label_used_for_experiments'] != 'Other') &
            (domain_e001_with_labels['label_used_for_experiments'].isin(eligible_classes))
        ].drop_duplicates(subset=['uniprot_id'])['label_used_for_experiments'].value_counts().to_dict()
        
        # Track excluded classes
        excluded_classes = []
        for cls, count in domain_e001_class_counts.items():
            if count < MIN_SAMPLES:
                excluded_classes.append({"class": cls, "count": count, "reason": f"n < {MIN_SAMPLES}"})
            elif cls in EXCLUDED_CLASSES_BIOLOGICAL:
                excluded_classes.append({"class": cls, "count": count, "reason": "Different catalytic mechanism/domain architecture"})
    else:
        supervised_eligible = []
        supervised_class_counts = {}
        eligible_classes = []
        excluded_classes = []
    
    manifest_path = manifests_dir / "supervised_eligible.txt"
    with open(manifest_path, 'w') as f:
        for uid in supervised_eligible:
            f.write(f"{uid}\n")
    
    report["datasets"]["supervised_eligible"] = {
        "manifest_file": str(manifest_path),
        "description": f"Sequences eligible for supervised learning (n >= {MIN_SAMPLES} per class)",
        "base_dataset": "domain_E001",
        "min_samples_per_class": MIN_SAMPLES,
        "n_sequences": len(supervised_eligible),
        "n_classes": len(supervised_class_counts),
        "per_class_counts": {k: int(v) for k, v in supervised_class_counts.items()},
        "excluded_classes": excluded_classes
    }
    
    print(f"   -> {len(supervised_eligible)} sequences in {len(supervised_class_counts)} classes")
    if excluded_classes:
        print(f"   -> Excluded {len(excluded_classes)} classes with < {MIN_SAMPLES} samples:")
        for exc in excluded_classes:
            print(f"      - {exc['class']}: {exc['count']} samples")
    
    # ========== Create Table 1 counts summary ==========
    print("\n5. Creating Table 1 summary...")
    
    # Total with labels
    total_labeled = len(labels_df)
    n_other = len(labels_df[labels_df['label_used_for_experiments'] == 'Other'])
    n_non_other = len(labels_df[labels_df['label_used_for_experiments'] != 'Other'])
    
    report["table_1_counts"] = {
        "whole_seq_total": total_labeled,
        "whole_seq_n_other": n_other,
        "whole_seq_n_non_other": n_non_other,
        "whole_seq_excl_other_n": len(whole_seq_excl_other),
        "domain_E0001_n": len(domain_e0001_excl_other),
        "domain_E001_n": len(domain_e001_excl_other),
        "supervised_eligible_n": len(supervised_eligible),
        "supervised_n_classes": len(supervised_class_counts)
    }
    
    # ========== Save report ==========
    print("\n6. Saving manifest report...")
    
    with open(output_report, 'w') as f:
        json.dump(report, f, indent=2)
    
    print(f"\n{'='*60}")
    print("STEP 6 COMPLETE: Dataset Manifests Created")
    print(f"{'='*60}")
    print(f"\nManifests saved to: {manifests_dir}/")
    print(f"Report saved to: {output_report}")
    print(f"\n{'Table 1 Summary':^60}")
    print("-" * 60)
    print(f"{'Dataset':<35} {'N':>10} {'Classes':>10}")
    print("-" * 60)
    print(f"{'Whole-seq (all labels)':<35} {total_labeled:>10}")
    print(f"{'Whole-seq (excl. Other)':<35} {len(whole_seq_excl_other):>10} {len(whole_seq_class_counts):>10}")
    print(f"{'Domains E<0.001 (excl. Other)':<35} {len(domain_e0001_excl_other):>10} {len(domain_e0001_class_counts):>10}")
    print(f"{'Domains E<0.01 (excl. Other)':<35} {len(domain_e001_excl_other):>10} {len(domain_e001_class_counts):>10}")
    print(f"{'Supervised-eligible (n≥5/class)':<35} {len(supervised_eligible):>10} {len(supervised_class_counts):>10}")
    print("-" * 60)
    
    if excluded_classes:
        print(f"\nExcluded from supervised learning (n < {MIN_SAMPLES}):")
        for exc in excluded_classes:
            print(f"  - {exc['class']}: {exc['count']} samples")
    
    print("\nPer-class counts for supervised-eligible dataset:")
    for cls, count in sorted(supervised_class_counts.items(), key=lambda x: -x[1]):
        print(f"  {cls:15} {count:>5}")


if __name__ == "__main__":
    main()

