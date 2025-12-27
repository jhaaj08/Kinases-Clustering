#!/usr/bin/env python3
"""
Step 9: Clustering (k-means) and Layer Ablation

This script runs k-means clustering on different embedding configurations
and creates a comprehensive registry of results.

Usage:
    python scripts/run_clustering.py

Inputs:
    - embeddings/esm2_t33_650M/*.npy (embedding files)
    - data/processed/labels.csv (ground truth labels)
    - data/processed/dataset_manifest_report.json (dataset info)

Outputs:
    - results/clustering/domain_E001_layer33.json
    - results/clustering/domain_E001_layers20_33.json
    - results/clustering/domain_E001_layers20_30.json
    - results/clustering/domain_E001_layer33_cls.json
    - results/clustering/summary_table.csv (Supplementary Table S1)
    - results/clustering/clustering_registry.json (master registry)
"""

import json
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime
from sklearn.cluster import KMeans
from sklearn.metrics import (
    adjusted_rand_score,
    normalized_mutual_info_score,
    homogeneity_score,
    completeness_score,
    v_measure_score,
    silhouette_score
)
from scipy.optimize import linear_sum_assignment


def compute_hungarian_accuracy(labels_true, labels_pred):
    """
    Compute accuracy using Hungarian algorithm for optimal cluster-to-label mapping.
    """
    from collections import Counter
    
    # Create confusion matrix
    labels_true = np.array(labels_true)
    labels_pred = np.array(labels_pred)
    
    unique_true = np.unique(labels_true)
    unique_pred = np.unique(labels_pred)
    
    # Cost matrix (negative of counts for maximization)
    cost_matrix = np.zeros((len(unique_pred), len(unique_true)))
    
    for i, pred_label in enumerate(unique_pred):
        for j, true_label in enumerate(unique_true):
            cost_matrix[i, j] = -np.sum((labels_pred == pred_label) & (labels_true == true_label))
    
    # Solve assignment problem
    row_ind, col_ind = linear_sum_assignment(cost_matrix)
    
    # Compute accuracy
    correct = 0
    for i, j in zip(row_ind, col_ind):
        pred_label = unique_pred[i]
        true_label = unique_true[j]
        correct += np.sum((labels_pred == pred_label) & (labels_true == true_label))
    
    return correct / len(labels_true)


def compute_purity(labels_true, labels_pred):
    """Compute cluster purity."""
    from collections import Counter
    
    labels_true = np.array(labels_true)
    labels_pred = np.array(labels_pred)
    
    total_correct = 0
    for cluster in np.unique(labels_pred):
        cluster_mask = labels_pred == cluster
        cluster_labels = labels_true[cluster_mask]
        if len(cluster_labels) > 0:
            most_common = Counter(cluster_labels).most_common(1)[0][1]
            total_correct += most_common
    
    return total_correct / len(labels_true)


def run_clustering(embeddings, labels, k, random_state=42):
    """
    Run k-means clustering and compute all metrics.
    """
    # Run k-means
    kmeans = KMeans(n_clusters=k, random_state=random_state, n_init=10)
    cluster_labels = kmeans.fit_predict(embeddings)
    
    # Compute metrics
    metrics = {
        "ARI": float(adjusted_rand_score(labels, cluster_labels)),
        "NMI": float(normalized_mutual_info_score(labels, cluster_labels)),
        "Homogeneity": float(homogeneity_score(labels, cluster_labels)),
        "Completeness": float(completeness_score(labels, cluster_labels)),
        "V_measure": float(v_measure_score(labels, cluster_labels)),
        "Silhouette": float(silhouette_score(embeddings, cluster_labels)),
        "Purity": float(compute_purity(labels, cluster_labels)),
        "Hungarian_Accuracy": float(compute_hungarian_accuracy(labels, cluster_labels))
    }
    
    return metrics, cluster_labels


def main():
    print("="*60)
    print("Step 9: Clustering (k-means) and Layer Ablation")
    print("="*60)
    
    # Paths
    embeddings_dir = Path("embeddings/esm2_t33_650M")
    output_dir = Path("results/clustering")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load labels
    print("\nLoading labels...")
    labels_df = pd.read_csv("data/processed/labels.csv")
    
    # Load embedding metadata
    metadata_file = embeddings_dir / "embedding_metadata.json"
    with open(metadata_file, 'r') as f:
        emb_metadata = json.load(f)
    
    # Load sequence IDs
    ids_file = embeddings_dir / "ids.txt"
    with open(ids_file, 'r') as f:
        embedding_ids = [line.strip() for line in f]
    
    print(f"  Loaded {len(embedding_ids)} embedding IDs")
    
    # Create ID to label mapping
    id_to_label = dict(zip(labels_df['uniprot_id'], labels_df['label_used_for_experiments']))
    
    # Get labels for embedding IDs (exclude "Other")
    valid_indices = []
    valid_labels = []
    for i, uid in enumerate(embedding_ids):
        label = id_to_label.get(uid)
        if label and label != "Other":
            valid_indices.append(i)
            valid_labels.append(label)
    
    print(f"  Valid sequences (excl. Other): {len(valid_labels)}")
    
    # Determine k (number of unique classes)
    unique_classes = sorted(set(valid_labels))
    k = len(unique_classes)
    print(f"  Number of classes (k): {k}")
    print(f"  Classes: {unique_classes}")
    
    # Registry for all results
    registry = {
        "step": 9,
        "name": "Clustering (k-means) Layer Ablation",
        "timestamp": datetime.now().isoformat(),
        "parameters": {
            "algorithm": "k-means",
            "k": k,
            "n_init": 10,
            "random_state": 42,
            "classes": unique_classes
        },
        "n_sequences": len(valid_labels),
        "experiments": {}
    }
    
    # Summary table rows
    summary_rows = []
    
    # Embedding configurations to test
    embedding_configs = [
        ("domain_E001_layer33_mean", "Layer 33 (mean)", [33], "mean"),
        ("domain_E001_layers20_33_mean", "Layers 20-33 (mean)", list(range(20, 34)), "mean"),
        ("domain_E001_layers20_30_mean", "Layers 20-30 (mean)", list(range(20, 31)), "mean"),
        ("domain_E001_layer33_cls", "Layer 33 (CLS)", [33], "cls"),
    ]
    
    print(f"\n{'='*60}")
    print("Running clustering experiments...")
    print(f"{'='*60}")
    
    for config_name, description, layers, pooling in embedding_configs:
        emb_file = embeddings_dir / f"{config_name}.npy"
        
        if not emb_file.exists():
            print(f"\n  Skip: {config_name} (file not found)")
            continue
        
        print(f"\n  {description}...")
        
        # Load embeddings
        embeddings = np.load(emb_file)
        
        # Filter to valid indices
        embeddings_filtered = embeddings[valid_indices]
        
        print(f"    Shape: {embeddings_filtered.shape}")
        
        # Run clustering
        metrics, cluster_labels = run_clustering(
            embeddings_filtered, 
            valid_labels, 
            k=k,
            random_state=42
        )
        
        print(f"    ARI: {metrics['ARI']:.4f}")
        print(f"    NMI: {metrics['NMI']:.4f}")
        print(f"    Hungarian: {metrics['Hungarian_Accuracy']:.4f}")
        
        # Get config hash from metadata
        config_hash = emb_metadata.get("embeddings", {}).get(config_name, {}).get("config_hash", "unknown")
        
        # Create experiment result
        experiment_result = {
            "name": config_name,
            "description": description,
            "embedding_file": str(emb_file),
            "config_hash": config_hash,
            "layers": layers,
            "pooling": pooling,
            "n_samples": len(valid_labels),
            "n_features": embeddings_filtered.shape[1],
            "k": k,
            "metrics": metrics
        }
        
        # Save individual JSON
        result_file = output_dir / f"{config_name}.json"
        with open(result_file, 'w') as f:
            json.dump(experiment_result, f, indent=2)
        
        # Add to registry
        registry["experiments"][config_name] = experiment_result
        
        # Add to summary table
        summary_rows.append({
            "Experiment": description,
            "Layers": f"[{layers[0]}]" if len(layers) == 1 else f"[{layers[0]}-{layers[-1]}]",
            "Pooling": pooling,
            "N": len(valid_labels),
            "k": k,
            "ARI": metrics["ARI"],
            "NMI": metrics["NMI"],
            "Homogeneity": metrics["Homogeneity"],
            "Completeness": metrics["Completeness"],
            "V_measure": metrics["V_measure"],
            "Silhouette": metrics["Silhouette"],
            "Purity": metrics["Purity"],
            "Hungarian": metrics["Hungarian_Accuracy"],
            "Config_Hash": config_hash[:8] + "..."
        })
    
    # Compute improvement percentages
    print(f"\n{'='*60}")
    print("Computing improvements...")
    print(f"{'='*60}")
    
    # Baseline is layer 33 mean
    baseline_ari = registry["experiments"].get("domain_E001_layer33_mean", {}).get("metrics", {}).get("ARI", 0)
    best_ari = 0
    best_config = ""
    
    for config_name, exp in registry["experiments"].items():
        ari = exp["metrics"]["ARI"]
        if ari > best_ari:
            best_ari = ari
            best_config = config_name
        
        # Compute improvement over baseline
        if baseline_ari > 0:
            improvement = (ari - baseline_ari) / baseline_ari * 100
            registry["experiments"][config_name]["improvement_vs_layer33"] = round(improvement, 2)
    
    registry["summary"] = {
        "baseline_config": "domain_E001_layer33_mean",
        "baseline_ARI": baseline_ari,
        "best_config": best_config,
        "best_ARI": best_ari,
        "improvement_percent": round((best_ari - baseline_ari) / baseline_ari * 100, 2) if baseline_ari > 0 else 0
    }
    
    print(f"\n  Baseline (Layer 33): ARI = {baseline_ari:.4f}")
    print(f"  Best ({best_config}): ARI = {best_ari:.4f}")
    print(f"  Improvement: +{registry['summary']['improvement_percent']:.1f}%")
    
    # Save registry
    registry_file = output_dir / "clustering_registry.json"
    with open(registry_file, 'w') as f:
        json.dump(registry, f, indent=2)
    
    # Save summary table
    summary_df = pd.DataFrame(summary_rows)
    summary_file = output_dir / "summary_table.csv"
    summary_df.to_csv(summary_file, index=False)
    
    print(f"\n{'='*60}")
    print("STEP 9 COMPLETE: Clustering Layer Ablation")
    print(f"{'='*60}")
    
    print(f"\nOutput files:")
    print(f"  - {registry_file}")
    print(f"  - {summary_file}")
    for config_name in registry["experiments"]:
        print(f"  - {output_dir / f'{config_name}.json'}")
    
    print(f"\n{'Summary Table (Supplementary Table S1)':^60}")
    print("-" * 80)
    print(f"{'Experiment':<25} {'Layers':<12} {'ARI':>8} {'NMI':>8} {'Hungarian':>10}")
    print("-" * 80)
    for row in summary_rows:
        print(f"{row['Experiment']:<25} {row['Layers']:<12} {row['ARI']:>8.4f} {row['NMI']:>8.4f} {row['Hungarian']:>10.4f}")
    print("-" * 80)
    
    print("\nSanity checks:")
    print(f"  ✓ All metrics computed from same dataset (n={len(valid_labels)})")
    print(f"  ✓ k={k} matches number of classes")
    print(f"  ✓ Registry JSON is source of truth for all reported values")


if __name__ == "__main__":
    main()

