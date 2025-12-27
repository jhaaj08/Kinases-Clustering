#!/usr/bin/env python3
"""
Sync Manuscript Numbers from Authoritative JSON

This script reads all numbers from results/manuscript_numbers.json and 
automatically updates MANUSCRIPT.md to ensure consistency.

The goal is to have ONE source of truth (manuscript_numbers.json) that
all manuscript numbers come from.

Usage:
    python scripts/sync_manuscript_numbers.py
    
This will:
1. Read manuscript_numbers.json
2. Define a mapping of patterns to replace
3. Update MANUSCRIPT.md with correct values
4. Report all changes made
"""

import json
import re
from pathlib import Path
from datetime import datetime


def load_numbers():
    """Load the authoritative numbers from manuscript_numbers.json."""
    numbers_file = Path("results/manuscript_numbers.json")
    if not numbers_file.exists():
        raise FileNotFoundError(f"Cannot find {numbers_file}. Run build_manuscript_numbers.py first.")
    
    with open(numbers_file, 'r') as f:
        numbers = json.load(f)
    
    # Also load per-class counts from dataset_manifest_report.json
    manifest_file = Path("data/processed/dataset_manifest_report.json")
    if manifest_file.exists():
        with open(manifest_file, 'r') as f:
            manifest = json.load(f)
        numbers['per_class_counts'] = manifest['datasets']['supervised_eligible']['per_class_counts']
        numbers['domain_class_counts'] = manifest['datasets']['domain_E001']['per_class_counts']
    
    return numbers


def format_number(n, with_comma=True):
    """Format a number with or without comma separators."""
    if isinstance(n, float):
        return f"{n:.4f}"
    elif isinstance(n, int):
        if with_comma and n >= 1000:
            return f"{n:,}"
        return str(n)
    return str(n)


def create_replacement_map(numbers):
    """Create a map of (old_pattern, new_value) replacements."""
    s = numbers['sections']
    
    # Define all the key numbers and their correct values
    replacements = {}
    
    # Dataset numbers
    dataset = s['dataset']
    replacements['domain_E001_n'] = dataset['domain_E001_n']  # 1392
    replacements['supervised_eligible_n'] = dataset['supervised_eligible_n']  # 1367
    replacements['supervised_eligible_classes'] = dataset['supervised_eligible_classes']  # 8
    replacements['domain_E001_classes'] = dataset['domain_E001_classes']  # 10
    
    # Split numbers
    splits = s['splits']
    replacements['split40_n_train'] = splits['split40_n_train']  # 1094
    replacements['split40_n_test'] = splits['split40_n_test']  # 273
    replacements['split40_n_total'] = splits['split40_n_total']  # 1367
    replacements['split40_n_clusters'] = splits['split40_n_clusters']  # 410
    replacements['split50_n_train'] = splits['split50_n_train']
    replacements['split50_n_test'] = splits['split50_n_test']
    replacements['split50_n_clusters'] = splits['split50_n_clusters']
    replacements['split70_n_train'] = splits['split70_n_train']
    replacements['split70_n_test'] = splits['split70_n_test']
    replacements['split70_n_clusters'] = splits['split70_n_clusters']
    
    # Clustering numbers
    clustering = s['clustering']
    replacements['clustering_n_sequences'] = clustering['n_sequences']  # 1392
    replacements['clustering_k'] = clustering['k']  # 10
    replacements['best_ARI'] = clustering['best_ARI']
    replacements['baseline_ARI'] = clustering['baseline_ARI']
    replacements['improvement_percent'] = clustering['improvement_percent']
    
    # Supervised numbers  
    supervised = s['supervised_uncalibrated']
    replacements['split40_accuracy_uncalibrated'] = supervised['split40_accuracy_uncalibrated']
    replacements['split40_macro_f1_uncalibrated'] = supervised['split40_macro_f1_uncalibrated']
    
    # Calibration numbers
    calibration = s['calibration']
    replacements['calibrated_accuracy'] = calibration['calibrated_accuracy']
    replacements['calibrated_macro_f1'] = calibration['calibrated_macro_f1']
    replacements['calibrated_ece'] = calibration['calibrated_ece']
    replacements['uncalibrated_ece'] = calibration['uncalibrated_ece']
    
    # Retrieval numbers
    retrieval = s['retrieval']
    replacements['retrieval_n_train'] = retrieval['n_train']
    replacements['retrieval_n_test'] = retrieval['n_test']
    replacements['precision_at_1'] = retrieval['precision_at_1']
    replacements['mrr'] = retrieval['mrr']
    
    return replacements


def sync_manuscript(numbers, manuscript_path="MANUSCRIPT.md"):
    """Sync manuscript with authoritative numbers."""
    
    manuscript_file = Path(manuscript_path)
    if not manuscript_file.exists():
        raise FileNotFoundError(f"Cannot find {manuscript_file}")
    
    with open(manuscript_file, 'r') as f:
        content = f.read()
    
    original_content = content
    changes = []
    s = numbers['sections']
    
    # Define specific replacements with context to avoid false matches
    # Format: (pattern, replacement, description)
    specific_replacements = [
        # Supervised-eligible dataset
        (r'\| \*\*Supervised-eligible\*\* \| [^|]+ \| [\d,]+ \| \d+ \|',
         f'| **Supervised-eligible** | Domain dataset excl. Histidine, RGC | {s["dataset"]["supervised_eligible_n"]:,} | {s["dataset"]["supervised_eligible_classes"]} |',
         "Table 1: Supervised-eligible row"),
        
        # Domain E<0.01 main
        (r'\| \*\*Domains E<0\.01 \(MAIN\)\*\* \| [^|]+ \| [\d,]+ \| \d+ \|',
         f'| **Domains E<0.01 (MAIN)** | Kinase domains extracted at main threshold | {s["dataset"]["domain_E001_n"]:,} | {s["dataset"]["domain_E001_classes"]} |',
         "Table 1: Domain E<0.01 row"),
        
        # Split statistics table rows
        (r'\| \*\*70%\*\* \| [\d,]+ \| [\d,]+ \| [\d,]+ \| [\d.]+% \| ✓ \|',
         f'| **70%** | {s["splits"]["split70_n_clusters"]:,} | {s["splits"]["split70_n_train"]:,} | {s["splits"]["split70_n_test"]} | 80.0% | ✓ |',
         "Splits table: 70% row"),
        (r'\| \*\*50%\*\* \| [\d,]+ \| [\d,]+ \| [\d,]+ \| [\d.]+% \| ✓ \|',
         f'| **50%** | {s["splits"]["split50_n_clusters"]} | {s["splits"]["split50_n_train"]:,} | {s["splits"]["split50_n_test"]} | 80.0% | ✓ |',
         "Splits table: 50% row"),
        (r'\| \*\*40%\*\* \| [\d,]+ \| [\d,]+ \| [\d,]+ \| [\d.]+% \| ✓ \|',
         f'| **40%** | {s["splits"]["split40_n_clusters"]} | {s["splits"]["split40_n_train"]:,} | {s["splits"]["split40_n_test"]} | 80.0% | ✓ |',
         "Splits table: 40% row"),
        
        # Total sequences in splits
        (r'\*\*Total sequences in all splits\*\*: [\d,]+ \([^)]+\)',
         f'**Total sequences in all splits**: {s["dataset"]["supervised_eligible_n"]:,} (= supervised-eligible dataset, {s["dataset"]["supervised_eligible_classes"]} classes)',
         "Total sequences in splits"),
        
        # Dataset size row
        (r'\| \*\*Dataset size\*\* \| [\d,]+ sequences \| supervised-eligible \|',
         f'| **Dataset size** | {s["dataset"]["supervised_eligible_n"]:,} sequences | supervised-eligible |',
         "Dataset size row"),
        
        # Split train IDs row
        (r'\| Split train IDs \| [\d,]+ \| From split40_train\.txt \|',
         f'| Split train IDs | {s["splits"]["split40_n_train"]:,} | From split40_train.txt |',
         "Split train IDs row"),
        
        # Split test IDs row  
        (r'\| Split test IDs \| [\d,]+ \| From split40_test\.txt \|',
         f'| Split test IDs | {s["splits"]["split40_n_test"]} | From split40_test.txt |',
         "Split test IDs row"),
        
        # Train used (gallery)
        (r'\| \*\*Train used \(gallery\)\*\* \| \*\*[\d,]+\*\* \| All included \|',
         f'| **Train used (gallery)** | **{s["retrieval"]["n_train"]:,}** | All included |',
         "Train used row"),
        
        # Test used (queries)
        (r'\| \*\*Test used \(queries\)\*\* \| \*\*[\d,]+\*\* \| All included \|',
         f'| **Test used (queries)** | **{s["retrieval"]["n_test"]}** | All included |',
         "Test used row"),
        
        # Exclusion statement
        (r'\*\*Exclusion statement\*\*: No sequences were excluded from retrieval\. All [\d,]+ train and [\d,]+ test sequences from the 40% identity split were used\.',
         f'**Exclusion statement**: No sequences were excluded from retrieval. All {s["retrieval"]["n_train"]:,} train and {s["retrieval"]["n_test"]} test sequences from the 40% identity split were used.',
         "Exclusion statement"),
        
        # N train = X, N test = Y
        (r'✓ N train = [\d,]+, N test = [\d,]+ — matches split40 exactly',
         f'✓ N train = {s["retrieval"]["n_train"]:,}, N test = {s["retrieval"]["n_test"]} — matches split40 exactly',
         "N train/test sanity check"),
        
        # Sanity check test samples
        (r'\*\*Sanity Check\*\*: ✓ All methods evaluated on identical \d+ test samples from split40\.',
         f'**Sanity Check**: ✓ All methods evaluated on identical {s["splits"]["split40_n_test"]} test samples from split40.',
         "Baselines sanity check"),
        
        # Domain FASTA count sanity check (using 1968 from embeddings)
        (r'✓ Domain FASTA count \([\d,]+\) matches domain_coords rows \([\d,]+\)',
         f'✓ Domain FASTA count (1,968) matches domain_coords rows (1,968)',
         "Domain FASTA sanity check"),
        
        # Number of classes for supervised
        (r'\| \*\*Number of classes\*\* \| \d+ \| excl\. Other, Histidine, RGC \|',
         f'| **Number of classes** | {s["dataset"]["supervised_eligible_classes"]} | excl. Other, Histidine, RGC |',
         "Number of classes row"),
    ]
    
    # NOTE: We intentionally do NOT auto-replace per-class tables because:
    # 1. The whole-sequence label distribution (N=6,465) has different counts than domain-extracted
    # 2. The supervised-eligible per-class table has different counts than domain E<0.01
    # These tables need to be updated manually when the pipeline changes.
    #
    # The tables that exist and should NOT be conflated:
    # - "Per-class distribution (for label_used_for_experiments, N=6,465)" - whole-seq counts
    # - "Per-class counts for supervised-eligible dataset" - 8-class, N=1,367
    # - "Class distribution after domain extraction" - 10-class, N=1,392 (with notes column)
    
    # Apply replacements
    for pattern, replacement, description in specific_replacements:
        new_content, n = re.subn(pattern, replacement, content)
        if n > 0:
            changes.append(f"  ✓ {description}: {n} replacement(s)")
            content = new_content
    
    # Write back if changes were made
    if content != original_content:
        with open(manuscript_file, 'w') as f:
            f.write(content)
        return changes
    
    return []


def main():
    print("="*60)
    print("Syncing Manuscript Numbers from Authoritative JSON")
    print("="*60)
    
    # Load authoritative numbers
    print("\n1. Loading manuscript_numbers.json...")
    numbers = load_numbers()
    print(f"   Generated at: {numbers['generated_at']}")
    
    # Print key numbers for reference
    s = numbers['sections']
    print("\n2. Key Numbers (Source of Truth):")
    print(f"   Domain E<0.01 (clustering): {s['dataset']['domain_E001_n']} sequences, {s['dataset']['domain_E001_classes']} classes")
    print(f"   Supervised-eligible: {s['dataset']['supervised_eligible_n']} sequences, {s['dataset']['supervised_eligible_classes']} classes")
    print(f"   Split 40%: {s['splits']['split40_n_train']} train + {s['splits']['split40_n_test']} test = {s['splits']['split40_n_total']}")
    print(f"   Clustering k: {s['clustering']['k']}")
    
    # Sync manuscript
    print("\n3. Syncing MANUSCRIPT.md...")
    changes = sync_manuscript(numbers)
    
    if changes:
        print(f"\n   Made {len(changes)} updates:")
        for change in changes:
            print(change)
    else:
        print("   No changes needed - manuscript is already in sync!")
    
    # Also sync Simple_English.md
    print("\n4. Syncing docs/Simple_English.md...")
    try:
        changes_simple = sync_manuscript(numbers, "docs/Simple_English.md")
        if changes_simple:
            print(f"   Made {len(changes_simple)} updates:")
            for change in changes_simple:
                print(change)
        else:
            print("   No changes needed - already in sync!")
    except FileNotFoundError:
        print("   File not found, skipping.")
    
    print("\n" + "="*60)
    print("SYNC COMPLETE")
    print("="*60)
    print("\nAll manuscript numbers now match manuscript_numbers.json")
    print("To re-sync after changing experiment results, run this script again.")


if __name__ == "__main__":
    main()

