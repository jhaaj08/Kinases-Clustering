#!/usr/bin/env python3
"""
Step 9: Clustering Experiments

This script runs k-means clustering on different embedding configurations
and evaluates performance against ground-truth labels.

Usage:
    python pipeline/step_09_clustering.py --run-dir runs/2025-01-01_000000/

Outputs:
    - results/clustering/clustering_registry.json
    - results/clustering/summary_table.csv
"""

import argparse
import json
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime
from collections import Counter
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from scipy.optimize import linear_sum_assignment

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from pipeline.membership import load_manifest

# Project root
PROJECT_ROOT = Path(__file__).parent.parent

# Clustering parameters
RANDOM_STATE = 42
N_INIT = 10


def hungarian_accuracy(labels_true, labels_pred):
    """Compute Hungarian (optimal) accuracy."""
    confusion = pd.crosstab(pd.Series(labels_true), pd.Series(labels_pred))
    row_ind, col_ind = linear_sum_assignment(-confusion.values)
    return confusion.values[row_ind, col_ind].sum() / len(labels_true)


def main():
    parser = argparse.ArgumentParser(description="Run clustering experiments")
    parser.add_argument("--run-dir", type=str, required=True, help="Run directory path")
    args = parser.parse_args()
    
    run_dir = Path(args.run_dir)
    results_dir = run_dir / "results" / "clustering"
    results_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 60)
    print("Step 9: Clustering Experiments")
    print("=" * 60)
    
    # Load manifest (domain_E001 - includes Histidine and RGC for clustering)
    manifest = load_manifest("domain_E001", run_dir)
    print(f"\nLoaded manifest: {len(manifest)} sequences")
    
    # Load embedding IDs
    ids_file = run_dir / "embeddings" / "esm2_t33_650M" / "ids.txt"
    with open(ids_file) as f:
        embedding_ids = [line.strip() for line in f if line.strip()]
    print(f"Embedding IDs: {len(embedding_ids)}")
    
    # Load labels
    labels_file = PROJECT_ROOT / "data" / "processed" / "labels.csv"
    labels_df = pd.read_csv(labels_file)
    id_to_label = dict(zip(labels_df['uniprot_id'], labels_df['label_used_for_experiments']))
    
    # Filter to manifest IDs that have embeddings
    valid_idx = []
    valid_ids = []
    valid_labels = []
    for i, uid in enumerate(embedding_ids):
        if uid in manifest and uid in id_to_label:
            valid_idx.append(i)
            valid_ids.append(uid)
            valid_labels.append(id_to_label[uid])
    
    unique_classes = sorted(set(valid_labels))
    k = len(unique_classes)
    
    print(f"\nValid sequences: {len(valid_ids)}")
    print(f"Classes (k): {k}")
    print(f"Class distribution:")
    for cls, count in sorted(Counter(valid_labels).items(), key=lambda x: -x[1]):
        print(f"  {cls:<15} {count}")
    
    # Embedding configurations
    configs = {
        "layer33_mean": {
            "file": "domain_E001_layer33_mean.npy",
            "description": "Final layer only (layer 33)",
            "layers": [33]
        },
        "layers20_30_mean": {
            "file": "domain_E001_layers20_30_mean.npy",
            "description": "Layers 20-30 averaged",
            "layers": list(range(20, 31))
        },
        "layers20_33_mean": {
            "file": "domain_E001_layers20_33_mean.npy",
            "description": "Layers 20-33 averaged",
            "layers": list(range(20, 34))
        },
        "layer33_cls": {
            "file": "domain_E001_layer33_cls.npy",
            "description": "Final layer only (CLS token)",
            "layers": [33]
        }
    }
    
    # Run clustering for each configuration
    results = {}
    baseline_ari = None
    
    print("\n" + "-" * 60)
    for config_name, config in configs.items():
        emb_file = run_dir / "embeddings" / "esm2_t33_650M" / config["file"]
        
        if not emb_file.exists():
            # Try symlink resolution
            if emb_file.is_symlink():
                emb_file = emb_file.resolve()
            else:
                print(f"⚠ Skipping {config_name}: file not found")
                continue
        
        embeddings = np.load(emb_file)
        embeddings_filtered = embeddings[valid_idx]
        
        print(f"\n{config['description']}...")
        print(f"  Shape: {embeddings_filtered.shape}")
        
        # Run k-means
        kmeans = KMeans(n_clusters=k, random_state=RANDOM_STATE, n_init=N_INIT)
        pred = kmeans.fit_predict(embeddings_filtered)
        
        # Compute metrics
        ari = adjusted_rand_score(valid_labels, pred)
        nmi = normalized_mutual_info_score(valid_labels, pred)
        hung = hungarian_accuracy(valid_labels, pred)
        
        print(f"  ARI: {ari:.4f}")
        print(f"  NMI: {nmi:.4f}")
        print(f"  Hungarian: {hung:.4f}")
        
        # Track baseline
        if config_name == "layer33_mean":
            baseline_ari = ari
        
        # Compute improvement
        improvement = 0
        if baseline_ari is not None and baseline_ari > 0:
            improvement = ((ari - baseline_ari) / baseline_ari) * 100
        
        results[config_name] = {
            "description": config["description"],
            "layers": config["layers"],
            "n": len(valid_ids),
            "k": k,
            "metrics": {
                "ARI": float(ari),
                "NMI": float(nmi),
                "Hungarian_Accuracy": float(hung)
            },
            "improvement_vs_layer33": float(improvement)
        }
    
    # Find best configuration
    best_config = max(results.keys(), key=lambda x: results[x]["metrics"]["ARI"])
    best_ari = results[best_config]["metrics"]["ARI"]
    
    # Create registry
    registry = {
        "step": 9,
        "name": "Clustering Experiments",
        "timestamp": datetime.now().isoformat(),
        "n_sequences": len(valid_ids),
        "parameters": {
            "k": k,
            "algorithm": "k-means",
            "n_init": N_INIT,
            "random_state": RANDOM_STATE
        },
        "experiments": results,
        "summary": {
            "baseline_config": "layer33_mean",
            "baseline_ARI": float(baseline_ari) if baseline_ari else 0,
            "best_config": best_config,
            "best_ARI": float(best_ari),
            "improvement_percent": float(((best_ari - baseline_ari) / baseline_ari) * 100) if baseline_ari else 0
        }
    }
    
    # Save registry
    registry_file = results_dir / "clustering_registry.json"
    with open(registry_file, 'w') as f:
        json.dump(registry, f, indent=2)
    print(f"\n✓ Saved: {registry_file}")
    
    # Save summary table
    summary_data = []
    for config_name, res in results.items():
        summary_data.append({
            "Configuration": res["description"],
            "Layers": str(res["layers"]),
            "ARI": f"{res['metrics']['ARI']:.4f}",
            "NMI": f"{res['metrics']['NMI']:.4f}",
            "Hungarian": f"{res['metrics']['Hungarian_Accuracy']:.4f}",
            "Improvement": f"{res['improvement_vs_layer33']:+.1f}%"
        })
    
    summary_df = pd.DataFrame(summary_data)
    summary_file = results_dir / "summary_table.csv"
    summary_df.to_csv(summary_file, index=False)
    print(f"✓ Saved: {summary_file}")
    
    print("\n" + "=" * 60)
    print("Step 9 COMPLETE")
    print("=" * 60)
    print(f"\nBest configuration: {best_config}")
    print(f"  ARI: {best_ari:.4f}")
    print(f"  Improvement: +{registry['summary']['improvement_percent']:.1f}% vs baseline")


if __name__ == "__main__":
    main()

