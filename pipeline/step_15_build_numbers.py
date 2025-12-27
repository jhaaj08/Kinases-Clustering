#!/usr/bin/env python3
"""
Step 15: Build Manuscript Numbers and Tables

This script aggregates all results from registries and generates:
- manuscript_numbers.json (all numbers in one place)
- Table1.csv (dataset construction)
- TableS1.csv (layer ablation)
- TableS2.csv (baselines comparison)

Usage:
    python pipeline/step_15_build_numbers.py --run-dir runs/2025-01-01_000000/

Outputs:
    - results/manuscript_numbers.json
    - tables/Table1.csv
    - tables/TableS1.csv
    - tables/TableS2.csv
"""

import argparse
import json
import pandas as pd
from pathlib import Path
from datetime import datetime

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from pipeline.membership import load_manifest, get_dataset_summary

# Project root
PROJECT_ROOT = Path(__file__).parent.parent


def load_json_safe(filepath):
    """Load JSON file, return None if not found."""
    if filepath.exists():
        with open(filepath) as f:
            return json.load(f)
    return None


def main():
    parser = argparse.ArgumentParser(description="Build manuscript numbers")
    parser.add_argument("--run-dir", type=str, required=True, help="Run directory path")
    args = parser.parse_args()
    
    run_dir = Path(args.run_dir)
    results_dir = run_dir / "results"
    tables_dir = run_dir / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 60)
    print("Step 15: Build Manuscript Numbers")
    print("=" * 60)
    
    # Collect all numbers
    numbers = {
        "step": 15,
        "name": "Manuscript Numbers",
        "timestamp": datetime.now().isoformat(),
        "run_dir": str(run_dir),
        "dataset": {},
        "clustering": {},
        "supervised": {},
        "calibration": {},
        "baselines": {},
        "retrieval": {}
    }
    
    # 1. Dataset numbers from manifests
    print("\n1. Loading dataset numbers...")
    manifest_report = load_json_safe(run_dir / "data" / "manifests" / "manifest_report.json")
    if manifest_report:
        numbers["dataset"] = {
            "domain_E001_n": manifest_report["datasets"]["domain_E001"]["n_sequences"],
            "domain_E001_classes": manifest_report["datasets"]["domain_E001"]["n_classes"],
            "supervised_eligible_n": manifest_report["datasets"]["supervised_eligible"]["n_sequences"],
            "supervised_eligible_classes": manifest_report["datasets"]["supervised_eligible"]["n_classes"],
            "per_class_counts": manifest_report["datasets"]["supervised_eligible"]["per_class_counts"]
        }
        print(f"  ✓ Dataset: {numbers['dataset']['supervised_eligible_n']} supervised-eligible")
    
    # 2. Split numbers
    print("\n2. Loading split numbers...")
    splits_report = load_json_safe(run_dir / "data" / "splits" / "splits_report.json")
    if splits_report:
        numbers["splits"] = {}
        for split_name, split_info in splits_report.get("splits", {}).items():
            numbers["splits"][split_name] = {
                "n_train": split_info["n_train"],
                "n_test": split_info["n_test"],
                "n_clusters": split_info["n_clusters"]
            }
        print(f"  ✓ Splits: {len(numbers['splits'])} thresholds")
    
    # 3. Clustering numbers
    print("\n3. Loading clustering numbers...")
    clustering_registry = load_json_safe(results_dir / "clustering" / "clustering_registry.json")
    if clustering_registry:
        numbers["clustering"] = {
            "n_sequences": clustering_registry["n_sequences"],
            "k": clustering_registry["parameters"]["k"],
            "baseline_config": clustering_registry["summary"]["baseline_config"],
            "baseline_ARI": clustering_registry["summary"]["baseline_ARI"],
            "best_config": clustering_registry["summary"]["best_config"],
            "best_ARI": clustering_registry["summary"]["best_ARI"],
            "improvement_percent": clustering_registry["summary"]["improvement_percent"]
        }
        
        # Add per-config metrics
        for config_name, config_data in clustering_registry["experiments"].items():
            numbers["clustering"][f"{config_name}_ARI"] = config_data["metrics"]["ARI"]
            numbers["clustering"][f"{config_name}_NMI"] = config_data["metrics"]["NMI"]
        
        print(f"  ✓ Clustering: best ARI = {numbers['clustering']['best_ARI']:.4f}")
    
    # 4. Supervised numbers
    print("\n4. Loading supervised numbers...")
    supervised_registry = load_json_safe(results_dir / "supervised" / "supervised_registry.json")
    if supervised_registry:
        for split_name, split_results in supervised_registry.get("experiments", {}).items():
            numbers["supervised"][split_name] = {}
            for config_name, config_data in split_results.items():
                numbers["supervised"][split_name][config_name] = config_data["metrics"]
        print(f"  ✓ Supervised: {len(numbers['supervised'])} splits")
    
    # 5. Calibration numbers
    print("\n5. Loading calibration numbers...")
    calibration = load_json_safe(results_dir / "calibration" / "split40_calibration.json")
    if calibration:
        numbers["calibration"] = {
            "uncalibrated_accuracy": calibration["uncalibrated"]["accuracy"],
            "uncalibrated_log_loss": calibration["uncalibrated"]["log_loss"],
            "uncalibrated_ece": calibration["uncalibrated"]["ece"],
            "calibrated_accuracy": calibration["calibrated"]["accuracy"],
            "calibrated_log_loss": calibration["calibrated"]["log_loss"],
            "calibrated_ece": calibration["calibrated"]["ece"]
        }
        print(f"  ✓ Calibration: accuracy {numbers['calibration']['calibrated_accuracy']:.4f}")
    
    # 6. Baselines numbers
    print("\n6. Loading baselines numbers...")
    baselines_dir = results_dir / "baselines"
    for baseline_file in baselines_dir.glob("*_split40.json"):
        baseline_name = baseline_file.stem.replace("_split40", "")
        baseline_data = load_json_safe(baseline_file)
        if baseline_data:
            numbers["baselines"][baseline_name] = baseline_data
    print(f"  ✓ Baselines: {len(numbers['baselines'])} methods")
    
    # 7. Retrieval numbers
    print("\n7. Loading retrieval numbers...")
    retrieval = load_json_safe(results_dir / "retrieval" / "split40_retrieval.json")
    if retrieval:
        numbers["retrieval"] = retrieval.get("metrics", {})
        print(f"  ✓ Retrieval: P@1 = {numbers['retrieval'].get('P@1', 'N/A')}")
    
    # Save manuscript_numbers.json
    numbers_file = results_dir / "manuscript_numbers.json"
    with open(numbers_file, 'w') as f:
        json.dump(numbers, f, indent=2)
    print(f"\n✓ Saved: {numbers_file}")
    
    # Generate Table 1: Dataset Construction
    print("\n8. Generating Table 1...")
    table1_data = []
    if "dataset" in numbers:
        table1_data.append({
            "Stage": "Domains (E < 0.01, excl. Other)",
            "N_sequences": numbers["dataset"].get("domain_E001_n", ""),
            "N_classes": numbers["dataset"].get("domain_E001_classes", "")
        })
        table1_data.append({
            "Stage": "Supervised-eligible",
            "N_sequences": numbers["dataset"].get("supervised_eligible_n", ""),
            "N_classes": numbers["dataset"].get("supervised_eligible_classes", "")
        })
    
    if "splits" in numbers:
        for split_name, split_info in numbers["splits"].items():
            table1_data.append({
                "Stage": f"{split_name} train",
                "N_sequences": split_info["n_train"],
                "N_classes": numbers["dataset"].get("supervised_eligible_classes", "")
            })
            table1_data.append({
                "Stage": f"{split_name} test",
                "N_sequences": split_info["n_test"],
                "N_classes": numbers["dataset"].get("supervised_eligible_classes", "")
            })
    
    table1_df = pd.DataFrame(table1_data)
    table1_file = tables_dir / "Table1.csv"
    table1_df.to_csv(table1_file, index=False)
    print(f"  ✓ Saved: {table1_file}")
    
    # Generate Table S1: Layer Ablation
    print("\n9. Generating Table S1...")
    if clustering_registry:
        tables1_data = []
        for config_name, config_data in clustering_registry["experiments"].items():
            tables1_data.append({
                "Configuration": config_data["description"],
                "Layers": str(config_data.get("layers", [])),
                "ARI": round(config_data["metrics"]["ARI"], 4),
                "NMI": round(config_data["metrics"]["NMI"], 4),
                "Hungarian_Accuracy": round(config_data["metrics"]["Hungarian_Accuracy"], 4),
                "Improvement_vs_Layer33": f"{config_data.get('improvement_vs_layer33', 0):+.1f}%"
            })
        
        tables1_df = pd.DataFrame(tables1_data)
        tables1_file = tables_dir / "TableS1.csv"
        tables1_df.to_csv(tables1_file, index=False)
        print(f"  ✓ Saved: {tables1_file}")
    
    # Generate Table S2: Baselines
    print("\n10. Generating Table S2...")
    tables2_data = []
    for baseline_name, baseline_data in numbers.get("baselines", {}).items():
        tables2_data.append({
            "Method": baseline_data.get("name", baseline_name),
            "Accuracy": round(baseline_data.get("accuracy", 0), 4),
            "Macro_F1": round(baseline_data.get("macro_f1", 0), 4),
            "Log_loss": round(baseline_data.get("log_loss", 0), 4) if baseline_data.get("log_loss") else "N/A"
        })
    
    if tables2_data:
        tables2_df = pd.DataFrame(tables2_data)
        tables2_file = tables_dir / "TableS2.csv"
        tables2_df.to_csv(tables2_file, index=False)
        print(f"  ✓ Saved: {tables2_file}")
    
    # Cross-validate
    print("\n11. Cross-validating...")
    errors = []
    
    # Check dataset numbers match
    if manifest_report and "dataset" in numbers:
        if numbers["dataset"]["supervised_eligible_n"] != manifest_report["datasets"]["supervised_eligible"]["n_sequences"]:
            errors.append("supervised_eligible_n mismatch")
    
    # Check split totals match supervised_eligible
    if "splits" in numbers and "dataset" in numbers:
        for split_name, split_info in numbers["splits"].items():
            total = split_info["n_train"] + split_info["n_test"]
            if total != numbers["dataset"]["supervised_eligible_n"]:
                errors.append(f"{split_name} total ({total}) != supervised_eligible ({numbers['dataset']['supervised_eligible_n']})")
    
    if errors:
        print("  ✗ Validation errors:")
        for e in errors:
            print(f"    - {e}")
    else:
        print("  ✓ All cross-validation checks passed")
    
    print("\n" + "=" * 60)
    print("Step 15 COMPLETE")
    print("=" * 60)
    print(f"\nGenerated files:")
    print(f"  - {numbers_file}")
    print(f"  - {tables_dir}/Table1.csv")
    print(f"  - {tables_dir}/TableS1.csv")
    print(f"  - {tables_dir}/TableS2.csv")
    
    print("\nKey numbers:")
    print(f"  Dataset: {numbers.get('dataset', {}).get('supervised_eligible_n', 'N/A')} sequences")
    print(f"  Clustering: Best ARI = {numbers.get('clustering', {}).get('best_ARI', 'N/A')}")
    print(f"  Supervised: Accuracy = {numbers.get('calibration', {}).get('calibrated_accuracy', 'N/A')}")
    print(f"  Retrieval: P@1 = {numbers.get('retrieval', {}).get('P@1', 'N/A')}")


if __name__ == "__main__":
    main()

