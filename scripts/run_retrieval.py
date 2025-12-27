#!/usr/bin/env python3
"""
Step 14: Retrieval Experiment (kNN Retrieval) with Exact N Reconciliation

This script runs a nearest-neighbor retrieval experiment where we query
with test sequences and retrieve the most similar training sequences,
checking if retrieved sequences share the same functional label.

IMPORTANT: This script explicitly tracks any excluded sequences to prevent
silent N discrepancies that could confuse readers.

Usage:
    python scripts/run_retrieval.py

Inputs:
    - embeddings/esm2_t33_650M/domain_E001_layers20_30_mean.npy
    - embeddings/esm2_t33_650M/ids.txt
    - data/splits/split40_train.txt, split40_test.txt
    - data/processed/labels.csv

Outputs:
    - results/retrieval/split40_retrieval.json
    - results/retrieval/excluded_ids.txt (if any exclusions)
    - results/retrieval/summary.csv
"""

import json
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime
from sklearn.metrics.pairwise import cosine_similarity
from collections import Counter


def load_split_ids(split_file):
    """Load IDs from split file."""
    with open(split_file, 'r') as f:
        return [line.strip() for line in f if line.strip()]


def compute_retrieval_metrics(results, k_values=[1, 3, 5, 10]):
    """Compute precision@k and MRR from retrieval results."""
    metrics = {}
    
    # Precision@k
    for k in k_values:
        correct = 0
        total = 0
        for r in results:
            retrieved_labels = r["retrieved_labels"][:k]
            if r["query_label"] in retrieved_labels:
                correct += 1
            total += 1
        metrics[f"precision_at_{k}"] = correct / total if total > 0 else 0
    
    # Mean Reciprocal Rank (MRR)
    reciprocal_ranks = []
    for r in results:
        query_label = r["query_label"]
        for rank, label in enumerate(r["retrieved_labels"], 1):
            if label == query_label:
                reciprocal_ranks.append(1.0 / rank)
                break
        else:
            reciprocal_ranks.append(0)
    
    metrics["mrr"] = np.mean(reciprocal_ranks) if reciprocal_ranks else 0
    
    return metrics


def main():
    print("="*60)
    print("Step 14: Retrieval Experiment (kNN Retrieval)")
    print("="*60)
    
    # Paths
    embeddings_dir = Path("embeddings/esm2_t33_650M")
    splits_dir = Path("data/splits")
    output_dir = Path("results/retrieval")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load embeddings
    print("\nLoading embeddings...")
    embedding_file = embeddings_dir / "domain_E001_layers20_30_mean.npy"
    embeddings = np.load(embedding_file)
    
    with open(embeddings_dir / "ids.txt", 'r') as f:
        embedding_ids = [line.strip() for line in f]
    
    embedding_id_set = set(embedding_ids)
    id_to_emb_idx = {uid: i for i, uid in enumerate(embedding_ids)}
    
    print(f"  Loaded {len(embedding_ids)} embeddings")
    
    # Load labels
    labels_df = pd.read_csv("data/processed/labels.csv")
    id_to_label = dict(zip(labels_df['uniprot_id'], labels_df['label_used_for_experiments']))
    
    # Load splits
    train_file = splits_dir / "split40_train.txt"
    test_file = splits_dir / "split40_test.txt"
    
    split_train_ids = load_split_ids(train_file)
    split_test_ids = load_split_ids(test_file)
    
    print(f"\n  Split files:")
    print(f"    Train: {len(split_train_ids)} IDs")
    print(f"    Test: {len(split_test_ids)} IDs")
    
    # =========================================================================
    # CRITICAL: Track exclusions explicitly
    # =========================================================================
    print(f"\n{'='*60}")
    print("Reconciling IDs (checking for exclusions)...")
    print(f"{'='*60}")
    
    # Check which train IDs are in embeddings
    train_excluded = []
    train_excluded_reasons = []
    train_valid_ids = []
    
    for uid in split_train_ids:
        if uid not in embedding_id_set:
            train_excluded.append(uid)
            train_excluded_reasons.append("not_in_embeddings")
        elif uid not in id_to_label:
            train_excluded.append(uid)
            train_excluded_reasons.append("no_label")
        elif id_to_label[uid] == 'Other':
            train_excluded.append(uid)
            train_excluded_reasons.append("label_is_Other")
        else:
            train_valid_ids.append(uid)
    
    # Check which test IDs are in embeddings
    test_excluded = []
    test_excluded_reasons = []
    test_valid_ids = []
    
    for uid in split_test_ids:
        if uid not in embedding_id_set:
            test_excluded.append(uid)
            test_excluded_reasons.append("not_in_embeddings")
        elif uid not in id_to_label:
            test_excluded.append(uid)
            test_excluded_reasons.append("no_label")
        elif id_to_label[uid] == 'Other':
            test_excluded.append(uid)
            test_excluded_reasons.append("label_is_Other")
        else:
            test_valid_ids.append(uid)
    
    # Report exclusions
    all_excluded = train_excluded + test_excluded
    all_reasons = train_excluded_reasons + test_excluded_reasons
    
    print(f"\n  Train exclusions: {len(train_excluded)}")
    print(f"  Test exclusions: {len(test_excluded)}")
    print(f"  Total exclusions: {len(all_excluded)}")
    
    if train_excluded:
        reason_counts = Counter(train_excluded_reasons)
        print(f"\n  Train exclusion reasons:")
        for reason, count in reason_counts.items():
            print(f"    - {reason}: {count}")
    
    if test_excluded:
        reason_counts = Counter(test_excluded_reasons)
        print(f"\n  Test exclusion reasons:")
        for reason, count in reason_counts.items():
            print(f"    - {reason}: {count}")
    
    # Save excluded IDs if any
    exclusions_file = output_dir / "excluded_ids.txt"
    if all_excluded:
        with open(exclusions_file, 'w') as f:
            f.write("# Excluded IDs from retrieval experiment\n")
            f.write(f"# Total: {len(all_excluded)} sequences\n")
            f.write("# Format: ID,source,reason\n\n")
            for uid, reason in zip(train_excluded, train_excluded_reasons):
                f.write(f"{uid},train,{reason}\n")
            for uid, reason in zip(test_excluded, test_excluded_reasons):
                f.write(f"{uid},test,{reason}\n")
        print(f"\n  Exclusions saved to: {exclusions_file}")
    else:
        # Create empty file to explicitly document no exclusions
        with open(exclusions_file, 'w') as f:
            f.write("# No sequences were excluded from retrieval experiment\n")
            f.write(f"# Train: {len(train_valid_ids)} (all included)\n")
            f.write(f"# Test: {len(test_valid_ids)} (all included)\n")
        print(f"\n  No exclusions! All split IDs are valid.")
    
    # Final counts
    print(f"\n  Final retrieval set:")
    print(f"    Train (gallery): {len(train_valid_ids)}")
    print(f"    Test (queries): {len(test_valid_ids)}")
    
    # =========================================================================
    # Run retrieval
    # =========================================================================
    print(f"\n{'='*60}")
    print("Running retrieval experiment...")
    print(f"{'='*60}")
    
    # Get embeddings for train (gallery) and test (queries)
    train_indices = [id_to_emb_idx[uid] for uid in train_valid_ids]
    test_indices = [id_to_emb_idx[uid] for uid in test_valid_ids]
    
    X_train = embeddings[train_indices]
    X_test = embeddings[test_indices]
    
    train_labels = [id_to_label[uid] for uid in train_valid_ids]
    test_labels = [id_to_label[uid] for uid in test_valid_ids]
    
    classes = sorted(set(train_labels + test_labels))
    print(f"  Classes: {len(classes)}")
    
    # Compute cosine similarities
    print("  Computing cosine similarities...")
    similarities = cosine_similarity(X_test, X_train)
    
    # For each test sample, get ranked list of training samples
    k_max = 10
    retrieval_results = []
    
    print(f"  Retrieving top-{k_max} neighbors for each query...")
    for i, (query_id, query_label) in enumerate(zip(test_valid_ids, test_labels)):
        # Get indices sorted by similarity (descending)
        ranked_indices = np.argsort(similarities[i])[::-1][:k_max]
        
        retrieved_ids = [train_valid_ids[j] for j in ranked_indices]
        retrieved_labels = [train_labels[j] for j in ranked_indices]
        retrieved_sims = [float(similarities[i, j]) for j in ranked_indices]
        
        retrieval_results.append({
            "query_id": query_id,
            "query_label": query_label,
            "retrieved_ids": retrieved_ids,
            "retrieved_labels": retrieved_labels,
            "retrieved_similarities": retrieved_sims
        })
    
    # Compute metrics
    print("  Computing retrieval metrics...")
    metrics = compute_retrieval_metrics(retrieval_results, k_values=[1, 3, 5, 10])
    
    print(f"\n  Results:")
    print(f"    Precision@1: {metrics['precision_at_1']:.4f}")
    print(f"    Precision@3: {metrics['precision_at_3']:.4f}")
    print(f"    Precision@5: {metrics['precision_at_5']:.4f}")
    print(f"    Precision@10: {metrics['precision_at_10']:.4f}")
    print(f"    MRR: {metrics['mrr']:.4f}")
    
    # =========================================================================
    # Save results
    # =========================================================================
    print(f"\n{'='*60}")
    print("Saving results...")
    print(f"{'='*60}")
    
    # Full results JSON
    result = {
        "step": 14,
        "name": "Retrieval Experiment (kNN)",
        "timestamp": datetime.now().isoformat(),
        "split": "40% identity threshold",
        "n_reconciliation": {
            "split_train_count": len(split_train_ids),
            "split_test_count": len(split_test_ids),
            "train_excluded_count": len(train_excluded),
            "test_excluded_count": len(test_excluded),
            "train_valid_count": len(train_valid_ids),
            "test_valid_count": len(test_valid_ids),
            "exclusion_reasons": dict(Counter(all_reasons)) if all_reasons else {},
            "exclusions_file": str(exclusions_file) if all_excluded else None,
            "statement": f"Used {len(train_valid_ids)} train and {len(test_valid_ids)} test sequences. "
                        f"{len(all_excluded)} excluded: {dict(Counter(all_reasons)) if all_reasons else 'none'}."
        },
        "embedding": {
            "file": str(embedding_file),
            "layers": "20-30 mean"
        },
        "config": {
            "similarity_metric": "cosine",
            "k_values_evaluated": [1, 3, 5, 10]
        },
        "metrics": metrics,
        "n_classes": len(classes),
        "classes": classes,
        "per_query_results": retrieval_results  # Full results for reproducibility
    }
    
    retrieval_file = output_dir / "split40_retrieval.json"
    with open(retrieval_file, 'w') as f:
        json.dump(result, f, indent=2)
    print(f"  Saved: {retrieval_file}")
    
    # Summary CSV
    summary_data = []
    for k in [1, 3, 5, 10]:
        summary_data.append({
            "Metric": f"Precision@{k}",
            "Value": metrics[f"precision_at_{k}"]
        })
    summary_data.append({"Metric": "MRR", "Value": metrics["mrr"]})
    summary_data.append({"Metric": "N_train", "Value": len(train_valid_ids)})
    summary_data.append({"Metric": "N_test", "Value": len(test_valid_ids)})
    summary_data.append({"Metric": "N_excluded", "Value": len(all_excluded)})
    
    summary_df = pd.DataFrame(summary_data)
    summary_file = output_dir / "summary.csv"
    summary_df.to_csv(summary_file, index=False)
    print(f"  Saved: {summary_file}")
    
    # =========================================================================
    # Final summary
    # =========================================================================
    print(f"\n{'='*60}")
    print("STEP 14 COMPLETE: Retrieval Experiment")
    print(f"{'='*60}")
    
    print(f"\n{'N Reconciliation':^60}")
    print("-" * 60)
    print(f"  Split train IDs:     {len(split_train_ids)}")
    print(f"  Split test IDs:      {len(split_test_ids)}")
    print(f"  Train excluded:      {len(train_excluded)}")
    print(f"  Test excluded:       {len(test_excluded)}")
    print(f"  ---")
    print(f"  Train used:          {len(train_valid_ids)}")
    print(f"  Test used:           {len(test_valid_ids)}")
    print("-" * 60)
    
    if all_excluded:
        print(f"\n⚠️  EXCLUSIONS OCCURRED!")
        print(f"    {len(all_excluded)} sequences were excluded")
        print(f"    Reasons: {dict(Counter(all_reasons))}")
        print(f"    List saved to: {exclusions_file}")
    else:
        print(f"\n✓ No exclusions - all split IDs used")
    
    print(f"\n{'Retrieval Metrics':^60}")
    print("-" * 60)
    print(f"  Precision@1:  {metrics['precision_at_1']:.4f}")
    print(f"  Precision@3:  {metrics['precision_at_3']:.4f}")
    print(f"  Precision@5:  {metrics['precision_at_5']:.4f}")
    print(f"  Precision@10: {metrics['precision_at_10']:.4f}")
    print(f"  MRR:          {metrics['mrr']:.4f}")
    print("-" * 60)
    
    print("\nSanity checks:")
    print(f"  ✓ All exclusions explicitly tracked and saved")
    print(f"  ✓ N reconciliation documented in JSON")
    print(f"  ✓ No silent exclusions allowed")


if __name__ == "__main__":
    main()

