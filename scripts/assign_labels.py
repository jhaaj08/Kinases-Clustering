#!/usr/bin/env python3
"""
Step 5: Label Assignment Script

This script creates a unified labels file with clear provenance tracking.
It builds labels.csv with columns:
  - label_original: The label from UniProt annotations
  - label_recovered: Labels recovered via parsing (optional)
  - label_used_for_experiments: The only column allowed downstream
  - label_source_tag: How the label was obtained (original/mapping/name_regex)

Usage:
    python scripts/assign_labels.py

Inputs:
    - data/processed/kinases_normalized.csv (from Step 4)

Outputs:
    - data/processed/labels.csv
    - data/processed/label_policy.json
    - data/processed/label_counts.tsv
"""

import pandas as pd
import json
from pathlib import Path
from collections import Counter


def main():
    # Paths
    input_file = Path("data/processed/kinases_normalized.csv")
    output_labels = Path("data/processed/labels.csv")
    output_policy = Path("data/processed/label_policy.json")
    output_counts = Path("data/processed/label_counts.tsv")
    
    print("=" * 60)
    print("Step 5: Label Assignment")
    print("=" * 60)
    
    # Load data
    print(f"\nLoading data from {input_file}...")
    df = pd.read_csv(input_file)
    print(f"Total sequences: {len(df)}")
    
    # Extract label information
    labels_df = pd.DataFrame()
    labels_df["uniprot_id"] = df["uniprot_id"]
    
    # Determine original vs recovered labels
    # - "original" label_source = direct annotation from UniProt
    # - "protein_name_parsing" and "subfamily_mapping" = recovered/inferred labels
    
    # label_original: The raw annotation (for sequences with original source)
    # For recovered labels, this would be what they had before recovery (usually "Other" or missing)
    labels_df["label_original"] = df.apply(
        lambda row: row["kinome_group_major"] if row["label_source"] == "original" else "Other",
        axis=1
    )
    
    # label_recovered: Labels that were inferred/recovered
    labels_df["label_recovered"] = df.apply(
        lambda row: row["kinome_group_major"] if row["label_source"] in ["protein_name_parsing", "subfamily_mapping"] else "",
        axis=1
    )
    
    # label_used_for_experiments: The column used for ALL downstream analyses
    # POLICY DECISION: We use the full label column (kinome_group_major) which includes recovered labels
    # This is documented explicitly in label_policy.json
    labels_df["label_used_for_experiments"] = df["kinome_group_major"]
    
    # label_source_tag: How the label was obtained
    source_mapping = {
        "original": "original",
        "protein_name_parsing": "name_regex",
        "subfamily_mapping": "mapping"
    }
    labels_df["label_source_tag"] = df["label_source"].map(source_mapping)
    
    # Save labels.csv
    print(f"\nSaving labels to {output_labels}...")
    labels_df.to_csv(output_labels, index=False)
    print(f"  Saved {len(labels_df)} label entries")
    
    # Create label policy JSON
    policy = {
        "official_column": "label_used_for_experiments",
        "description": "All analyses use 'label_used_for_experiments' column. This includes both original UniProt annotations and recovered labels.",
        "label_recovery_policy": {
            "method": "Protein name parsing and subfamily mapping were used to recover labels for sequences initially lacking kinome group annotations.",
            "original_source": "Direct UniProt kinome_group field annotation",
            "name_regex_source": "Regular expression matching on protein names (e.g., 'CAMK' in name -> CAMK group)",
            "mapping_source": "Subfamily to major group hierarchical mapping"
        },
        "experiment_usage": "All clustering, supervised classification, and retrieval experiments use 'label_used_for_experiments' column exclusively.",
        "version": "1.0",
        "created_by": "scripts/assign_labels.py"
    }
    
    print(f"Saving label policy to {output_policy}...")
    with open(output_policy, "w") as f:
        json.dump(policy, f, indent=2)
    
    # Create label counts TSV
    print(f"Saving label counts to {output_counts}...")
    
    # Count by label column
    counts_data = []
    
    # Experiment labels (the one that matters)
    exp_counts = Counter(labels_df["label_used_for_experiments"])
    for label, count in sorted(exp_counts.items(), key=lambda x: -x[1]):
        counts_data.append({
            "label_column": "label_used_for_experiments",
            "label_value": label,
            "count": count
        })
    
    # Original labels (for transparency)
    orig_counts = Counter(labels_df["label_original"])
    for label, count in sorted(orig_counts.items(), key=lambda x: -x[1]):
        counts_data.append({
            "label_column": "label_original",
            "label_value": label,
            "count": count
        })
    
    # Source tags
    source_counts = Counter(labels_df["label_source_tag"])
    for source, count in sorted(source_counts.items(), key=lambda x: -x[1]):
        counts_data.append({
            "label_column": "label_source_tag",
            "label_value": source,
            "count": count
        })
    
    counts_df = pd.DataFrame(counts_data)
    counts_df.to_csv(output_counts, sep="\t", index=False)
    
    # Print summary
    print("\n" + "=" * 60)
    print("LABEL ASSIGNMENT SUMMARY")
    print("=" * 60)
    
    # Key counts
    n_other = exp_counts.get("Other", 0)
    n_non_other = sum(c for l, c in exp_counts.items() if l != "Other")
    
    print(f"\nTotal sequences: {len(labels_df)}")
    print(f"N_other (Other class): {n_other}")
    print(f"N_non_other (all other classes): {n_non_other}")
    
    print("\nPer-class counts (label_used_for_experiments):")
    print("-" * 40)
    for label, count in sorted(exp_counts.items(), key=lambda x: -x[1]):
        pct = 100 * count / len(labels_df)
        print(f"  {label:15} : {count:5} ({pct:5.1f}%)")
    
    print("\nLabel sources:")
    print("-" * 40)
    for source, count in sorted(source_counts.items(), key=lambda x: -x[1]):
        pct = 100 * count / len(labels_df)
        print(f"  {source:15} : {count:5} ({pct:5.1f}%)")
    
    print("\n" + "=" * 60)
    print("SANITY CHECKS")
    print("=" * 60)
    
    # Check 1: Only one column used for experiments
    print("\n✓ Official experiment column: 'label_used_for_experiments'")
    
    # Check 2: All sequences have a label
    missing_labels = labels_df["label_used_for_experiments"].isna().sum()
    if missing_labels == 0:
        print("✓ All sequences have labels (no missing values)")
    else:
        print(f"✗ WARNING: {missing_labels} sequences have missing labels!")
    
    # Check 3: No empty labels
    empty_labels = (labels_df["label_used_for_experiments"] == "").sum()
    if empty_labels == 0:
        print("✓ No empty label strings")
    else:
        print(f"✗ WARNING: {empty_labels} sequences have empty labels!")
    
    # Check 4: Count consistency
    total_from_counts = sum(exp_counts.values())
    if total_from_counts == len(labels_df):
        print(f"✓ Label counts sum to total ({total_from_counts} = {len(labels_df)})")
    else:
        print(f"✗ WARNING: Count mismatch: {total_from_counts} != {len(labels_df)}")
    
    print("\n" + "=" * 60)
    print("Step 5 complete!")
    print(f"  - labels.csv: {output_labels}")
    print(f"  - label_policy.json: {output_policy}")
    print(f"  - label_counts.tsv: {output_counts}")
    print("=" * 60)
    
    return {
        "total_sequences": len(labels_df),
        "n_other": n_other,
        "n_non_other": n_non_other,
        "class_counts": dict(exp_counts),
        "source_counts": dict(source_counts)
    }


if __name__ == "__main__":
    main()

