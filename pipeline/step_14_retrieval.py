#!/usr/bin/env python3
"""
Step 14: Retrieval Experiment

This script evaluates nearest-neighbor retrieval quality using embeddings.
Computes Precision@k and Mean Reciprocal Rank (MRR).

Usage:
    python pipeline/step_14_retrieval.py --run-dir runs/2025-01-01_000000/

Outputs:
    - results/retrieval/split40_retrieval.json
    - results/retrieval/summary.csv
"""

import argparse
import json
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime
from sklearn.metrics.pairwise import cosine_similarity

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from pipeline.membership import load_split

# Project root
PROJECT_ROOT = Path(__file__).parent.parent


def compute_retrieval_metrics(train_embeddings, test_embeddings, train_labels, test_labels, k_values=[1, 3, 5, 10]):
    """Compute retrieval metrics: Precision@k and MRR."""
    # Compute cosine similarity between test and train
    similarities = cosine_similarity(test_embeddings, train_embeddings)
    
    metrics = {}
    reciprocal_ranks = []
    
    for k in k_values:
        correct_at_k = 0
        
        for i, true_label in enumerate(test_labels):
            # Get top-k neighbors
            top_k_indices = np.argsort(similarities[i])[-k:][::-1]
            top_k_labels = [train_labels[j] for j in top_k_indices]
            
            # Check if true label is in top-k
            if true_label in top_k_labels:
                correct_at_k += 1
            
            # For MRR (only compute once at k=max)
            if k == max(k_values):
                all_indices = np.argsort(similarities[i])[::-1]
                for rank, idx in enumerate(all_indices, 1):
                    if train_labels[idx] == true_label:
                        reciprocal_ranks.append(1.0 / rank)
                        break
                else:
                    reciprocal_ranks.append(0.0)
        
        metrics[f"P@{k}"] = correct_at_k / len(test_labels)
    
    metrics["MRR"] = np.mean(reciprocal_ranks)
    
    return metrics


def main():
    parser = argparse.ArgumentParser(description="Run retrieval experiment")
    parser.add_argument("--run-dir", type=str, required=True, help="Run directory path")
    args = parser.parse_args()
    
    run_dir = Path(args.run_dir)
    results_dir = run_dir / "results" / "retrieval"
    results_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 60)
    print("Step 14: Retrieval Experiment")
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
    
    # Get indices and labels (only for IDs that have embeddings)
    train_idx = []
    train_labels = []
    excluded_train = []
    for uid in train_ids:
        if uid in id_to_idx:
            train_idx.append(id_to_idx[uid])
            train_labels.append(id_to_label[uid])
        else:
            excluded_train.append(uid)
    
    test_idx = []
    test_labels = []
    excluded_test = []
    for uid in test_ids:
        if uid in id_to_idx:
            test_idx.append(id_to_idx[uid])
            test_labels.append(id_to_label[uid])
        else:
            excluded_test.append(uid)
    
    print(f"  Train with embeddings: {len(train_idx)}")
    print(f"  Test with embeddings: {len(test_idx)}")
    
    if excluded_train or excluded_test:
        print(f"  ⚠ Excluded (no embeddings): {len(excluded_train)} train, {len(excluded_test)} test")
        
        # Save excluded IDs
        excluded_file = results_dir / "excluded_ids.txt"
        with open(excluded_file, 'w') as f:
            for uid in excluded_train + excluded_test:
                f.write(f"{uid}\n")
    
    X_train = embeddings[train_idx]
    X_test = embeddings[test_idx]
    
    # Compute retrieval metrics
    print("\nComputing retrieval metrics...")
    k_values = [1, 3, 5, 10]
    metrics = compute_retrieval_metrics(X_train, X_test, train_labels, test_labels, k_values)
    
    for k in k_values:
        print(f"  P@{k}: {metrics[f'P@{k}']:.4f}")
    print(f"  MRR: {metrics['MRR']:.4f}")
    
    # Create report
    report = {
        "step": 14,
        "name": "Retrieval Experiment",
        "timestamp": datetime.now().isoformat(),
        "split": split_name,
        "n_train": len(train_idx),
        "n_test": len(test_idx),
        "n_excluded_train": len(excluded_train),
        "n_excluded_test": len(excluded_test),
        "metrics": {k: float(v) for k, v in metrics.items()}
    }
    
    report_file = results_dir / f"{split_name}_retrieval.json"
    with open(report_file, 'w') as f:
        json.dump(report, f, indent=2)
    print(f"\n✓ Saved: {report_file}")
    
    # Save summary
    summary_data = [{"Metric": k, "Value": f"{v:.4f}"} for k, v in metrics.items()]
    summary_df = pd.DataFrame(summary_data)
    summary_file = results_dir / "summary.csv"
    summary_df.to_csv(summary_file, index=False)
    print(f"✓ Saved: {summary_file}")
    
    print("\n" + "=" * 60)
    print("Step 14 COMPLETE")
    print("=" * 60)
    print("\n              Retrieval Summary")
    print("-" * 40)
    for k, v in metrics.items():
        print(f"  {k:<10} {v:.4f}")
    print("-" * 40)


if __name__ == "__main__":
    main()

