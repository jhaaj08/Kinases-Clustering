#!/usr/bin/env python3
"""
Step 13: Baseline Comparisons

This script runs baseline methods (k-NN, Random, MLP, Motifs) for comparison
with the main logistic regression model.

Usage:
    python pipeline/step_13_baselines.py --run-dir runs/2025-01-01_000000/

Outputs:
    - results/baselines/baselines_split40.csv
    - results/baselines/knn_split40.json
    - results/baselines/random_split40.json
    - results/baselines/mlp_split40.json
"""

import argparse
import json
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.dummy import DummyClassifier
from sklearn.metrics import accuracy_score, f1_score, log_loss

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from pipeline.membership import load_split

# Project root
PROJECT_ROOT = Path(__file__).parent.parent

# Configuration
RANDOM_STATE = 42


def main():
    parser = argparse.ArgumentParser(description="Run baseline comparisons")
    parser.add_argument("--run-dir", type=str, required=True, help="Run directory path")
    args = parser.parse_args()
    
    run_dir = Path(args.run_dir)
    results_dir = run_dir / "results" / "baselines"
    results_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 60)
    print("Step 13: Baseline Comparisons")
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
    
    # Load embeddings (layer33_mean)
    emb_file = run_dir / "embeddings" / "esm2_t33_650M" / "domain_E001_layer33_mean.npy"
    embeddings = np.load(emb_file)
    
    # Use split40
    split_name = "split40"
    train_ids = load_split(split_name, "train", run_dir)
    test_ids = load_split(split_name, "test", run_dir)
    
    print(f"\nUsing {split_name}: {len(train_ids)} train, {len(test_ids)} test")
    
    # Get indices and labels
    train_idx = [id_to_idx[uid] for uid in train_ids if uid in id_to_idx]
    test_idx = [id_to_idx[uid] for uid in test_ids if uid in id_to_idx]
    
    train_labels = [id_to_label[uid] for uid in train_ids if uid in id_to_label]
    test_labels = [id_to_label[uid] for uid in test_ids if uid in id_to_label]
    
    X_train = embeddings[train_idx]
    X_test = embeddings[test_idx]
    
    classes = sorted(set(train_labels))
    
    # Baselines to run
    baselines = {
        "knn": {
            "name": "k-Nearest Neighbors",
            "model": KNeighborsClassifier(n_neighbors=5, metric='cosine')
        },
        "random": {
            "name": "Random (stratified)",
            "model": DummyClassifier(strategy='stratified', random_state=RANDOM_STATE)
        },
        "mlp": {
            "name": "MLP (1 hidden layer)",
            "model": MLPClassifier(
                hidden_layer_sizes=(256,),
                max_iter=500,
                random_state=RANDOM_STATE
            )
        }
    }
    
    results = []
    
    for baseline_key, baseline_info in baselines.items():
        print(f"\n{baseline_info['name']}...")
        
        model = baseline_info['model']
        model.fit(X_train, train_labels)
        
        y_pred = model.predict(X_test)
        
        # Get probabilities if available
        try:
            y_proba = model.predict_proba(X_test)
            ll = log_loss(test_labels, y_proba, labels=model.classes_)
        except AttributeError:
            ll = None
        
        acc = accuracy_score(test_labels, y_pred)
        f1 = f1_score(test_labels, y_pred, average='macro', zero_division=0)
        
        print(f"  Accuracy: {acc:.4f}")
        print(f"  Macro-F1: {f1:.4f}")
        
        baseline_result = {
            "name": baseline_info['name'],
            "accuracy": float(acc),
            "macro_f1": float(f1),
            "log_loss": float(ll) if ll else None
        }
        
        results.append({
            "Method": baseline_info['name'],
            "Accuracy": f"{acc:.4f}",
            "Macro-F1": f"{f1:.4f}",
            "Log-loss": f"{ll:.4f}" if ll else "N/A"
        })
        
        # Save individual result
        result_file = results_dir / f"{baseline_key}_{split_name}.json"
        with open(result_file, 'w') as f:
            json.dump(baseline_result, f, indent=2)
        print(f"  ✓ Saved: {result_file.name}")
    
    # Save summary table
    summary_df = pd.DataFrame(results)
    summary_file = results_dir / f"baselines_{split_name}.csv"
    summary_df.to_csv(summary_file, index=False)
    print(f"\n✓ Saved: {summary_file}")
    
    # Save metadata
    metadata = {
        "step": 13,
        "name": "Baseline Comparisons",
        "timestamp": datetime.now().isoformat(),
        "split": split_name,
        "n_train": len(train_idx),
        "n_test": len(test_idx),
        "n_classes": len(classes)
    }
    
    metadata_file = results_dir / "baselines_metadata.json"
    with open(metadata_file, 'w') as f:
        json.dump(metadata, f, indent=2)
    
    print("\n" + "=" * 60)
    print("Step 13 COMPLETE")
    print("=" * 60)
    print("\n              Baselines Summary")
    print("-" * 60)
    print(f"{'Method':<25} {'Accuracy':<12} {'Macro-F1':<12}")
    print("-" * 60)
    for r in results:
        print(f"{r['Method']:<25} {r['Accuracy']:<12} {r['Macro-F1']:<12}")
    print("-" * 60)


if __name__ == "__main__":
    main()

