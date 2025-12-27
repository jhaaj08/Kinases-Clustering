#!/usr/bin/env python3
"""
Step 11: Supervised Learning Experiments

This script trains logistic regression classifiers on different layer configurations
and evaluates them across different homology thresholds.

Usage:
    python pipeline/step_11_supervised.py --run-dir runs/2025-01-01_000000/

Outputs:
    - results/supervised/supervised_registry.json
    - results/supervised/lr_split{40,50,70}_metrics.json
"""

import argparse
import json
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score, f1_score, precision_score, recall_score,
    log_loss, classification_report
)
import joblib

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from pipeline.membership import load_manifest, load_split

# Project root
PROJECT_ROOT = Path(__file__).parent.parent

# Configuration
RANDOM_STATE = 42


def compute_metrics(y_true, y_pred, y_proba, classes):
    """Compute classification metrics."""
    return {
        "accuracy": float(accuracy_score(y_true, y_pred)),
        "macro_f1": float(f1_score(y_true, y_pred, average='macro', zero_division=0)),
        "weighted_f1": float(f1_score(y_true, y_pred, average='weighted', zero_division=0)),
        "macro_precision": float(precision_score(y_true, y_pred, average='macro', zero_division=0)),
        "macro_recall": float(recall_score(y_true, y_pred, average='macro', zero_division=0)),
        "log_loss": float(log_loss(y_true, y_proba, labels=classes))
    }


def main():
    parser = argparse.ArgumentParser(description="Run supervised learning experiments")
    parser.add_argument("--run-dir", type=str, required=True, help="Run directory path")
    args = parser.parse_args()
    
    run_dir = Path(args.run_dir)
    results_dir = run_dir / "results" / "supervised"
    results_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 60)
    print("Step 11: Supervised Learning")
    print("=" * 60)
    
    # Load embedding IDs
    ids_file = run_dir / "embeddings" / "esm2_t33_650M" / "ids.txt"
    with open(ids_file) as f:
        embedding_ids = [line.strip() for line in f if line.strip()]
    id_to_idx = {uid: i for i, uid in enumerate(embedding_ids)}
    
    # Load labels
    labels_file = PROJECT_ROOT / "data" / "processed" / "labels.csv"
    labels_df = pd.read_csv(labels_file)
    id_to_label = dict(zip(labels_df['uniprot_id'], labels_df['label_used_for_experiments']))
    
    # Embedding configurations
    configs = {
        "layer33_mean": "domain_E001_layer33_mean.npy",
        "layers20_30_mean": "domain_E001_layers20_30_mean.npy"
    }
    
    # Load embeddings
    embeddings = {}
    for config_name, filename in configs.items():
        emb_file = run_dir / "embeddings" / "esm2_t33_650M" / filename
        if emb_file.exists() or emb_file.is_symlink():
            embeddings[config_name] = np.load(emb_file)
            print(f"Loaded {config_name}: {embeddings[config_name].shape}")
    
    # Results registry
    registry = {
        "step": 11,
        "name": "Supervised Learning",
        "timestamp": datetime.now().isoformat(),
        "experiments": {}
    }
    
    # Run experiments for each threshold
    thresholds = [40, 50, 70]
    
    for threshold in thresholds:
        split_name = f"split{threshold}"
        print(f"\n{'=' * 60}")
        print(f"Processing {threshold}% identity threshold")
        print(f"{'=' * 60}")
        
        try:
            train_ids = load_split(split_name, "train", run_dir)
            test_ids = load_split(split_name, "test", run_dir)
        except FileNotFoundError:
            print(f"  ⚠ Split files not found, skipping")
            continue
        
        print(f"  Train: {len(train_ids)}, Test: {len(test_ids)}")
        
        # Get indices and labels
        train_idx = [id_to_idx[uid] for uid in train_ids if uid in id_to_idx]
        test_idx = [id_to_idx[uid] for uid in test_ids if uid in id_to_idx]
        
        train_labels = [id_to_label[uid] for uid in train_ids if uid in id_to_label]
        test_labels = [id_to_label[uid] for uid in test_ids if uid in id_to_label]
        
        classes = sorted(set(train_labels))
        print(f"  Classes: {len(classes)}")
        
        split_results = {}
        
        for config_name, emb in embeddings.items():
            print(f"\n  {config_name}...")
            
            X_train = emb[train_idx]
            X_test = emb[test_idx]
            
            # Train model
            model = LogisticRegression(
                max_iter=1000,
                random_state=RANDOM_STATE,
                multi_class='multinomial',
                solver='lbfgs'
            )
            model.fit(X_train, train_labels)
            
            # Predict
            y_pred = model.predict(X_test)
            y_proba = model.predict_proba(X_test)
            
            # Compute metrics
            metrics = compute_metrics(test_labels, y_pred, y_proba, model.classes_)
            
            print(f"    Accuracy: {metrics['accuracy']:.4f}")
            print(f"    Macro-F1: {metrics['macro_f1']:.4f}")
            
            split_results[config_name] = {
                "n_train": len(train_idx),
                "n_test": len(test_idx),
                "n_classes": len(classes),
                "metrics": metrics
            }
            
            # Save model for 40% split
            if threshold == 40:
                model_file = results_dir / f"lr_{split_name}_{config_name}.joblib"
                joblib.dump(model, model_file)
        
        # Save split metrics
        metrics_file = results_dir / f"lr_{split_name}_metrics.json"
        with open(metrics_file, 'w') as f:
            json.dump(split_results, f, indent=2)
        print(f"\n  ✓ Saved: {metrics_file.name}")
        
        registry["experiments"][split_name] = split_results
    
    # Save registry
    registry_file = results_dir / "supervised_registry.json"
    with open(registry_file, 'w') as f:
        json.dump(registry, f, indent=2)
    print(f"\n✓ Saved: {registry_file}")
    
    # Create summary table
    summary_data = []
    for split_name, split_results in registry["experiments"].items():
        for config_name, res in split_results.items():
            summary_data.append({
                "Split": split_name,
                "Config": config_name,
                "Accuracy": f"{res['metrics']['accuracy']:.4f}",
                "Macro-F1": f"{res['metrics']['macro_f1']:.4f}",
                "Log-loss": f"{res['metrics']['log_loss']:.4f}"
            })
    
    summary_df = pd.DataFrame(summary_data)
    summary_file = results_dir / "lr_multi_identity_summary.csv"
    summary_df.to_csv(summary_file, index=False)
    print(f"✓ Saved: {summary_file}")
    
    print("\n" + "=" * 60)
    print("Step 11 COMPLETE")
    print("=" * 60)


if __name__ == "__main__":
    main()

