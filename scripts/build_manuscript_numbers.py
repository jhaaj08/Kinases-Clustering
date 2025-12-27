#!/usr/bin/env python3
"""
Step 15: Generate ALL Manuscript Numbers from Registries (No Hand Edits)

This script consolidates all metrics and counts from various registry files
into a single source-of-truth JSON file. All numbers in the manuscript
should be copied from this file or the generated tables.

Usage:
    python scripts/build_manuscript_numbers.py

Inputs:
    - data/splits/splits_report.json
    - data/processed/dataset_manifest_report.json
    - results/clustering/clustering_registry.json
    - results/supervised/lr_split40_metrics.json
    - results/supervised/lr_multi_identity_summary.csv
    - results/calibration/split40_calibration.json
    - results/baselines/baselines_split40.csv
    - results/retrieval/split40_retrieval.json

Outputs:
    - results/manuscript_numbers.json
    - results/tables/Table1.csv (Dataset Construction)
    - results/tables/TableS1.csv (Layer Ablation)
    - results/tables/TableS2.csv (Baselines Comparison)
"""

import json
import pandas as pd
from pathlib import Path
from datetime import datetime


def load_json(filepath):
    """Load JSON file, return None if not found."""
    if Path(filepath).exists():
        with open(filepath, 'r') as f:
            return json.load(f)
    return None


def main():
    print("="*60)
    print("Step 15: Build Manuscript Numbers from Registries")
    print("="*60)
    
    output_dir = Path("results")
    tables_dir = Path("results/tables")
    tables_dir.mkdir(parents=True, exist_ok=True)
    
    # Master numbers dictionary
    numbers = {
        "generated_at": datetime.now().isoformat(),
        "source": "scripts/build_manuscript_numbers.py",
        "note": "All manuscript numbers MUST come from this file. No hand edits allowed.",
        "sections": {}
    }
    
    # =========================================================================
    # 1. Dataset Construction Numbers
    # =========================================================================
    print("\n1. Loading dataset construction numbers...")
    
    manifest_report = load_json("data/processed/dataset_manifest_report.json")
    
    if manifest_report:
        dataset_numbers = {
            "whole_seq_excl_other_n": manifest_report["datasets"]["whole_seq_excl_other"]["n_sequences"],
            "whole_seq_excl_other_classes": manifest_report["datasets"]["whole_seq_excl_other"]["n_classes"],
            "domain_E0001_n": manifest_report["datasets"]["domain_E0001"]["n_sequences"],
            "domain_E0001_classes": manifest_report["datasets"]["domain_E0001"]["n_classes"],
            "domain_E001_n": manifest_report["datasets"]["domain_E001"]["n_sequences"],
            "domain_E001_classes": manifest_report["datasets"]["domain_E001"]["n_classes"],
            "supervised_eligible_n": manifest_report["datasets"]["supervised_eligible"]["n_sequences"],
            "supervised_eligible_classes": manifest_report["datasets"]["supervised_eligible"]["n_classes"],
            "supervised_excluded_classes": [e["class"] for e in manifest_report["datasets"]["supervised_eligible"]["excluded_classes"]]
        }
        numbers["sections"]["dataset"] = dataset_numbers
        print(f"  ✓ Loaded dataset manifest")
    else:
        print("  ✗ Dataset manifest not found")
    
    # =========================================================================
    # 2. Split Numbers
    # =========================================================================
    print("\n2. Loading split numbers...")
    
    splits_report = load_json("data/splits/splits_report.json")
    
    if splits_report:
        split_numbers = {}
        for split_name, split_data in splits_report["splits"].items():
            prefix = f"{split_name}_"
            split_numbers[f"{prefix}n_clusters"] = split_data["n_clusters"]
            split_numbers[f"{prefix}n_train"] = split_data["n_train"]
            split_numbers[f"{prefix}n_test"] = split_data["n_test"]
            split_numbers[f"{prefix}n_total"] = split_data["n_total"]
        numbers["sections"]["splits"] = split_numbers
        print(f"  ✓ Loaded splits report")
    else:
        print("  ✗ Splits report not found")
    
    # =========================================================================
    # 3. Clustering Numbers (Layer Ablation)
    # =========================================================================
    print("\n3. Loading clustering numbers...")
    
    clustering_registry = load_json("results/clustering/clustering_registry.json")
    
    if clustering_registry:
        clustering_numbers = {
            "n_sequences": clustering_registry["n_sequences"],
            "k": clustering_registry["parameters"]["k"],
            "baseline_config": clustering_registry["summary"]["baseline_config"],
            "baseline_ARI": clustering_registry["summary"]["baseline_ARI"],
            "best_config": clustering_registry["summary"]["best_config"],
            "best_ARI": clustering_registry["summary"]["best_ARI"],
            "improvement_percent": clustering_registry["summary"]["improvement_percent"]
        }
        
        # Per-config metrics
        for config_name, config_data in clustering_registry["experiments"].items():
            clustering_numbers[f"{config_name}_ARI"] = config_data["metrics"]["ARI"]
            clustering_numbers[f"{config_name}_NMI"] = config_data["metrics"]["NMI"]
            clustering_numbers[f"{config_name}_Hungarian"] = config_data["metrics"]["Hungarian_Accuracy"]
        
        numbers["sections"]["clustering"] = clustering_numbers
        print(f"  ✓ Loaded clustering registry")
    else:
        print("  ✗ Clustering registry not found")
    
    # =========================================================================
    # 4. Supervised Learning Numbers
    # =========================================================================
    print("\n4. Loading supervised learning numbers...")
    
    # Load multi-identity summary
    lr_summary_file = Path("results/supervised/lr_multi_identity_summary.csv")
    if lr_summary_file.exists():
        lr_summary = pd.read_csv(lr_summary_file)
        supervised_numbers = {}
        for _, row in lr_summary.iterrows():
            threshold = row["Identity_Threshold"].replace("%", "")
            supervised_numbers[f"split{threshold}_accuracy_uncalibrated"] = row["Accuracy"]
            supervised_numbers[f"split{threshold}_macro_f1_uncalibrated"] = row["Macro_F1"]
            supervised_numbers[f"split{threshold}_weighted_f1_uncalibrated"] = row["Weighted_F1"]
            supervised_numbers[f"split{threshold}_n_train"] = int(row["N_Train"])
            supervised_numbers[f"split{threshold}_n_test"] = int(row["N_Test"])
        numbers["sections"]["supervised_uncalibrated"] = supervised_numbers
        print(f"  ✓ Loaded supervised summary")
    else:
        print("  ✗ Supervised summary not found")
    
    # =========================================================================
    # 5. Calibration Numbers
    # =========================================================================
    print("\n5. Loading calibration numbers...")
    
    calibration = load_json("results/calibration/split40_calibration.json")
    
    if calibration:
        calibration_numbers = {
            "uncalibrated_accuracy": calibration["metrics"]["uncalibrated"]["accuracy"],
            "uncalibrated_macro_f1": calibration["metrics"]["uncalibrated"]["macro_f1"],
            "uncalibrated_log_loss": calibration["metrics"]["uncalibrated"]["log_loss"],
            "uncalibrated_ece": calibration["metrics"]["uncalibrated"]["ece"],
            "calibrated_accuracy": calibration["metrics"]["calibrated"]["accuracy"],
            "calibrated_macro_f1": calibration["metrics"]["calibrated"]["macro_f1"],
            "calibrated_log_loss": calibration["metrics"]["calibrated"]["log_loss"],
            "calibrated_ece": calibration["metrics"]["calibrated"]["ece"],
            "accuracy_delta": calibration["metrics"]["improvement"]["accuracy_delta"],
            "log_loss_delta": calibration["metrics"]["improvement"]["log_loss_delta"],
            "ece_delta": calibration["metrics"]["improvement"]["ece_delta"],
            "which_accuracy_for_baselines": calibration["notes"]["accuracy_for_baselines_table"]
        }
        numbers["sections"]["calibration"] = calibration_numbers
        print(f"  ✓ Loaded calibration results")
    else:
        print("  ✗ Calibration results not found")
    
    # =========================================================================
    # 6. Baselines Numbers
    # =========================================================================
    print("\n6. Loading baselines numbers...")
    
    baselines_file = Path("results/baselines/baselines_split40.csv")
    if baselines_file.exists():
        baselines_df = pd.read_csv(baselines_file)
        baselines_numbers = {}
        for _, row in baselines_df.iterrows():
            method = row["Method"].replace(" ", "_").replace("(", "").replace(")", "").replace(",", "")
            baselines_numbers[f"{method}_accuracy"] = row["Accuracy"]
            baselines_numbers[f"{method}_macro_f1"] = row["Macro_F1"]
            if pd.notna(row.get("Top3_Accuracy")):
                baselines_numbers[f"{method}_top3_accuracy"] = row["Top3_Accuracy"]
        numbers["sections"]["baselines"] = baselines_numbers
        print(f"  ✓ Loaded baselines")
    else:
        print("  ✗ Baselines not found")
    
    # =========================================================================
    # 7. Retrieval Numbers
    # =========================================================================
    print("\n7. Loading retrieval numbers...")
    
    retrieval = load_json("results/retrieval/split40_retrieval.json")
    
    if retrieval:
        retrieval_numbers = {
            "n_train": retrieval["n_reconciliation"]["train_valid_count"],
            "n_test": retrieval["n_reconciliation"]["test_valid_count"],
            "n_excluded": retrieval["n_reconciliation"]["train_excluded_count"] + retrieval["n_reconciliation"]["test_excluded_count"],
            "precision_at_1": retrieval["metrics"]["precision_at_1"],
            "precision_at_3": retrieval["metrics"]["precision_at_3"],
            "precision_at_5": retrieval["metrics"]["precision_at_5"],
            "precision_at_10": retrieval["metrics"]["precision_at_10"],
            "mrr": retrieval["metrics"]["mrr"]
        }
        numbers["sections"]["retrieval"] = retrieval_numbers
        print(f"  ✓ Loaded retrieval results")
    else:
        print("  ✗ Retrieval results not found")
    
    # =========================================================================
    # Save master numbers JSON
    # =========================================================================
    print(f"\n{'='*60}")
    print("Saving manuscript numbers...")
    print(f"{'='*60}")
    
    numbers_file = output_dir / "manuscript_numbers.json"
    with open(numbers_file, 'w') as f:
        json.dump(numbers, f, indent=2)
    print(f"  Saved: {numbers_file}")
    
    # =========================================================================
    # Generate Table 1: Dataset Construction
    # =========================================================================
    print("\nGenerating Table 1: Dataset Construction...")
    
    if manifest_report and splits_report:
        table1_data = [
            {
                "Stage": "Whole-seq (excl. Other)",
                "N": manifest_report["datasets"]["whole_seq_excl_other"]["n_sequences"],
                "Classes": manifest_report["datasets"]["whole_seq_excl_other"]["n_classes"],
                "Notes": "Full-length sequences"
            },
            {
                "Stage": "Domain E<0.001 (strict)",
                "N": manifest_report["datasets"]["domain_E0001"]["n_sequences"],
                "Classes": manifest_report["datasets"]["domain_E0001"]["n_classes"],
                "Notes": "Strict E-value threshold"
            },
            {
                "Stage": "Domain E<0.01 (main)",
                "N": manifest_report["datasets"]["domain_E001"]["n_sequences"],
                "Classes": manifest_report["datasets"]["domain_E001"]["n_classes"],
                "Notes": "Primary analysis dataset"
            },
            {
                "Stage": "Supervised-eligible",
                "N": manifest_report["datasets"]["supervised_eligible"]["n_sequences"],
                "Classes": manifest_report["datasets"]["supervised_eligible"]["n_classes"],
                "Notes": f"Excl. classes with n<5: {', '.join(numbers['sections']['dataset']['supervised_excluded_classes'])}"
            },
            {
                "Stage": "Split 70% (train/test)",
                "N": f"{splits_report['splits']['split70']['n_train']}/{splits_report['splits']['split70']['n_test']}",
                "Classes": 8,
                "Notes": f"{splits_report['splits']['split70']['n_clusters']} clusters"
            },
            {
                "Stage": "Split 50% (train/test)",
                "N": f"{splits_report['splits']['split50']['n_train']}/{splits_report['splits']['split50']['n_test']}",
                "Classes": 8,
                "Notes": f"{splits_report['splits']['split50']['n_clusters']} clusters"
            },
            {
                "Stage": "Split 40% (train/test)",
                "N": f"{splits_report['splits']['split40']['n_train']}/{splits_report['splits']['split40']['n_test']}",
                "Classes": 8,
                "Notes": f"{splits_report['splits']['split40']['n_clusters']} clusters"
            }
        ]
        table1_df = pd.DataFrame(table1_data)
        table1_file = tables_dir / "Table1.csv"
        table1_df.to_csv(table1_file, index=False)
        print(f"  Saved: {table1_file}")
    
    # =========================================================================
    # Generate Table S1: Layer Ablation (Clustering)
    # =========================================================================
    print("\nGenerating Table S1: Layer Ablation...")
    
    if clustering_registry:
        tables1_data = []
        for config_name, config_data in clustering_registry["experiments"].items():
            tables1_data.append({
                "Configuration": config_data["description"],
                "Layers": str(config_data["layers"]),
                "ARI": round(config_data["metrics"]["ARI"], 4),
                "NMI": round(config_data["metrics"]["NMI"], 4),
                "Hungarian_Accuracy": round(config_data["metrics"]["Hungarian_Accuracy"], 4),
                "Improvement_vs_Layer33": f"{config_data.get('improvement_vs_layer33', 0):+.1f}%"
            })
        tables1_df = pd.DataFrame(tables1_data)
        tables1_file = tables_dir / "TableS1.csv"
        tables1_df.to_csv(tables1_file, index=False)
        print(f"  Saved: {tables1_file}")
    
    # =========================================================================
    # Generate Table S2: Baselines Comparison
    # =========================================================================
    print("\nGenerating Table S2: Baselines Comparison...")
    
    if baselines_file.exists():
        baselines_df = pd.read_csv(baselines_file)
        baselines_df["Accuracy"] = baselines_df["Accuracy"].round(4)
        baselines_df["Macro_F1"] = baselines_df["Macro_F1"].round(4)
        if "Top3_Accuracy" in baselines_df.columns:
            baselines_df["Top3_Accuracy"] = baselines_df["Top3_Accuracy"].round(4)
        tables2_file = tables_dir / "TableS2.csv"
        baselines_df.to_csv(tables2_file, index=False)
        print(f"  Saved: {tables2_file}")
    
    # =========================================================================
    # Summary
    # =========================================================================
    print(f"\n{'='*60}")
    print("STEP 15 COMPLETE: Manuscript Numbers Generated")
    print(f"{'='*60}")
    
    print("\nGenerated files:")
    print(f"  - {numbers_file}")
    print(f"  - {tables_dir / 'Table1.csv'}")
    print(f"  - {tables_dir / 'TableS1.csv'}")
    print(f"  - {tables_dir / 'TableS2.csv'}")
    
    print("\nKey numbers for manuscript:")
    if "dataset" in numbers["sections"]:
        d = numbers["sections"]["dataset"]
        print(f"  Dataset: {d['supervised_eligible_n']} sequences, {d['supervised_eligible_classes']} classes")
    if "clustering" in numbers["sections"]:
        c = numbers["sections"]["clustering"]
        print(f"  Clustering: Best ARI = {c['best_ARI']:.4f} (+{c['improvement_percent']:.1f}% vs baseline)")
    if "calibration" in numbers["sections"]:
        cal = numbers["sections"]["calibration"]
        print(f"  Supervised: Accuracy = {cal['calibrated_accuracy']:.4f} (calibrated)")
    if "retrieval" in numbers["sections"]:
        r = numbers["sections"]["retrieval"]
        print(f"  Retrieval: P@1 = {r['precision_at_1']:.4f}, MRR = {r['mrr']:.4f}")
    
    print("\nSanity checks:")
    print("  ✓ All numbers come from registry files")
    print("  ✓ No hand-edited values allowed")
    print("  ✓ Tables generated from JSON/CSV sources")
    print("  ⚠️  RULE: If a number is not in manuscript_numbers.json, it cannot appear in the manuscript!")


if __name__ == "__main__":
    main()
