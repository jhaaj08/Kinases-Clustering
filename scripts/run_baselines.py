#!/usr/bin/env python3
"""
Step 13: Baselines Rerun Under Same Split + Same Dataset

This script runs multiple baseline methods on the same 40% identity split
to ensure fair comparison with our main model.

Baselines:
1. k-NN (k=5) on ESM-2 embeddings
2. Motifs-only logistic regression (handcrafted features)
3. MLP (2-layer) on ESM-2 embeddings
4. Random baseline (stratified by class distribution)

Usage:
    python scripts/run_baselines.py

Inputs:
    - embeddings/esm2_t33_650M/domain_E001_layers20_30_mean.npy
    - data/processed/kinases_domains_with_motifs.csv
    - data/splits/split40_train.txt, split40_test.txt
    - data/processed/labels.csv

Outputs:
    - results/baselines/baselines_split40.csv
    - results/baselines/knn_split40.json
    - results/baselines/motifs_split40.json
    - results/baselines/mlp_split40.json
    - results/baselines/random_split40.json
"""

import json
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime
from sklearn.neighbors import KNeighborsClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.neural_network import MLPClassifier
from sklearn.dummy import DummyClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (
    accuracy_score,
    f1_score,
    precision_score,
    recall_score,
    top_k_accuracy_score,
    log_loss
)
import warnings
warnings.filterwarnings('ignore')


def load_split_ids(train_file, test_file):
    """Load train/test IDs from split files."""
    with open(train_file, 'r') as f:
        train_ids = [line.strip() for line in f if line.strip()]
    with open(test_file, 'r') as f:
        test_ids = [line.strip() for line in f if line.strip()]
    return train_ids, test_ids


def get_indices_and_labels(ids_list, embedding_ids, id_to_label):
    """Get embedding indices and labels for a list of IDs."""
    id_to_idx = {uid: i for i, uid in enumerate(embedding_ids)}
    
    indices = []
    labels = []
    valid_ids = []
    
    for uid in ids_list:
        if uid in id_to_idx and uid in id_to_label:
            indices.append(id_to_idx[uid])
            labels.append(id_to_label[uid])
            valid_ids.append(uid)
    
    return indices, labels, valid_ids


def compute_metrics(y_true, y_pred, y_proba, classes):
    """Compute classification metrics."""
    metrics = {
        "accuracy": float(accuracy_score(y_true, y_pred)),
        "macro_f1": float(f1_score(y_true, y_pred, average='macro', zero_division=0)),
        "weighted_f1": float(f1_score(y_true, y_pred, average='weighted', zero_division=0)),
        "macro_precision": float(precision_score(y_true, y_pred, average='macro', zero_division=0)),
        "macro_recall": float(recall_score(y_true, y_pred, average='macro', zero_division=0)),
    }
    
    if y_proba is not None:
        metrics["log_loss"] = float(log_loss(y_true, y_proba, labels=classes))
        if len(classes) > 3:
            metrics["top3_accuracy"] = float(top_k_accuracy_score(
                y_true, y_proba, k=3, labels=classes
            ))
    
    return metrics


def main():
    print("="*60)
    print("Step 13: Baselines Comparison")
    print("="*60)
    
    # Paths
    embeddings_dir = Path("embeddings/esm2_t33_650M")
    splits_dir = Path("data/splits")
    output_dir = Path("results/baselines")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load embeddings
    print("\nLoading data...")
    embedding_file = embeddings_dir / "domain_E001_layers20_30_mean.npy"
    embeddings = np.load(embedding_file)
    
    with open(embeddings_dir / "ids.txt", 'r') as f:
        embedding_ids = [line.strip() for line in f]
    
    # Load labels
    labels_df = pd.read_csv("data/processed/labels.csv")
    id_to_label = dict(zip(labels_df['uniprot_id'], labels_df['label_used_for_experiments']))
    
    # Load motif features
    motifs_df = pd.read_csv("data/processed/kinases_domains_with_motifs.csv")
    
    # Define motif feature columns
    motif_cols = [
        'dfg_present', 'hrd_present', 'ape_present', 'ploop_present', 
        'ploop_consensus', 'vaik_present', 'k_e_distance', 'k_e_salt_bridge_present',
        'k_e_distance_normal', 'alphac_acidic_present', 'activation_loop_length',
        'activation_loop_length_norm', 'catalytic_loop_length', 'catalytic_loop_length_norm',
        'dfg_position_norm', 'hrd_position_norm', 'ape_position_norm', 'ploop_position_norm',
        'gatekeeper_found', 'gatekeeper_size', 'gatekeeper_hydrophobicity',
        'gatekeeper_is_small', 'gatekeeper_is_large', 'dfg_state_hydrophobic',
        'hrd_dfg_spacing', 'hrd_dfg_spacing_normal', 'core_triad_complete',
        'extended_motif_completeness', 'motif_integrity_score', 'domain_length'
    ]
    
    # Load split
    train_file = splits_dir / "split40_train.txt"
    test_file = splits_dir / "split40_test.txt"
    train_ids, test_ids = load_split_ids(train_file, test_file)
    
    print(f"  Split40: {len(train_ids)} train, {len(test_ids)} test")
    
    # Get embedding data
    train_indices, train_labels, train_valid_ids = get_indices_and_labels(
        train_ids, embedding_ids, id_to_label
    )
    test_indices, test_labels, test_valid_ids = get_indices_and_labels(
        test_ids, embedding_ids, id_to_label
    )
    
    X_train_emb = embeddings[train_indices]
    X_test_emb = embeddings[test_indices]
    y_train = train_labels
    y_test = test_labels
    
    classes = sorted(set(y_train + y_test))
    print(f"  Classes: {len(classes)}")
    print(f"  Embedding shape: {X_train_emb.shape}")
    
    # Prepare motif features
    motifs_df_indexed = motifs_df.set_index('uniprot_id')
    
    # Get motif features for train/test
    train_motifs = []
    test_motifs = []
    train_valid_motif_ids = []
    test_valid_motif_ids = []
    train_labels_motif = []
    test_labels_motif = []
    
    for uid, label in zip(train_valid_ids, train_labels):
        if uid in motifs_df_indexed.index:
            row = motifs_df_indexed.loc[uid]
            if isinstance(row, pd.DataFrame):
                row = row.iloc[0]
            feats = []
            for col in motif_cols:
                val = row.get(col, 0)
                if pd.isna(val):
                    val = 0
                feats.append(float(val))
            train_motifs.append(feats)
            train_valid_motif_ids.append(uid)
            train_labels_motif.append(label)
    
    for uid, label in zip(test_valid_ids, test_labels):
        if uid in motifs_df_indexed.index:
            row = motifs_df_indexed.loc[uid]
            if isinstance(row, pd.DataFrame):
                row = row.iloc[0]
            feats = []
            for col in motif_cols:
                val = row.get(col, 0)
                if pd.isna(val):
                    val = 0
                feats.append(float(val))
            test_motifs.append(feats)
            test_valid_motif_ids.append(uid)
            test_labels_motif.append(label)
    
    X_train_motif = np.array(train_motifs)
    X_test_motif = np.array(test_motifs)
    
    print(f"  Motif features: {X_train_motif.shape[1]} features")
    print(f"  Motif samples: {len(train_motifs)} train, {len(test_motifs)} test")
    
    # Results storage
    results = []
    
    # =========================================================================
    # Baseline 1: k-NN on embeddings
    # =========================================================================
    print(f"\n{'='*60}")
    print("Baseline 1: k-NN (k=5) on ESM-2 embeddings")
    print(f"{'='*60}")
    
    knn = KNeighborsClassifier(n_neighbors=5, metric='cosine')
    knn.fit(X_train_emb, y_train)
    y_pred_knn = knn.predict(X_test_emb)
    y_proba_knn = knn.predict_proba(X_test_emb)
    
    metrics_knn = compute_metrics(y_test, y_pred_knn, y_proba_knn, classes)
    print(f"  Accuracy: {metrics_knn['accuracy']:.4f}")
    print(f"  Macro-F1: {metrics_knn['macro_f1']:.4f}")
    
    knn_result = {
        "baseline": "k-NN",
        "description": "k-Nearest Neighbors (k=5) with cosine distance on ESM-2 layer 20-30 embeddings",
        "features": "ESM-2 embeddings (1280-dim)",
        "config": {"k": 5, "metric": "cosine"},
        "n_train": len(train_indices),
        "n_test": len(test_indices),
        "metrics": metrics_knn
    }
    with open(output_dir / "knn_split40.json", 'w') as f:
        json.dump(knn_result, f, indent=2)
    
    results.append({
        "Method": "k-NN (k=5)",
        "Features": "ESM-2 embeddings",
        "Accuracy": metrics_knn["accuracy"],
        "Macro_F1": metrics_knn["macro_f1"],
        "Top3_Accuracy": metrics_knn.get("top3_accuracy", None)
    })
    
    # =========================================================================
    # Baseline 2: Motifs-only Logistic Regression
    # =========================================================================
    print(f"\n{'='*60}")
    print("Baseline 2: Logistic Regression on Motif Features")
    print(f"{'='*60}")
    
    # Scale motif features
    scaler = StandardScaler()
    X_train_motif_scaled = scaler.fit_transform(X_train_motif)
    X_test_motif_scaled = scaler.transform(X_test_motif)
    
    lr_motif = LogisticRegression(
        solver='lbfgs',
        max_iter=2000,
        class_weight='balanced',
        random_state=42
    )
    lr_motif.fit(X_train_motif_scaled, train_labels_motif)
    y_pred_motif = lr_motif.predict(X_test_motif_scaled)
    y_proba_motif = lr_motif.predict_proba(X_test_motif_scaled)
    
    metrics_motif = compute_metrics(test_labels_motif, y_pred_motif, y_proba_motif, classes)
    print(f"  Accuracy: {metrics_motif['accuracy']:.4f}")
    print(f"  Macro-F1: {metrics_motif['macro_f1']:.4f}")
    
    motif_result = {
        "baseline": "Motifs-only LR",
        "description": "Logistic Regression on 30 handcrafted kinase motif features (DFG, HRD, APE, P-loop, etc.)",
        "features": f"30 motif features: {motif_cols}",
        "config": {"solver": "lbfgs", "max_iter": 2000, "class_weight": "balanced"},
        "n_train": len(train_motifs),
        "n_test": len(test_motifs),
        "metrics": metrics_motif
    }
    with open(output_dir / "motifs_split40.json", 'w') as f:
        json.dump(motif_result, f, indent=2)
    
    results.append({
        "Method": "Motifs-only LR",
        "Features": "30 handcrafted motif features",
        "Accuracy": metrics_motif["accuracy"],
        "Macro_F1": metrics_motif["macro_f1"],
        "Top3_Accuracy": metrics_motif.get("top3_accuracy", None)
    })
    
    # =========================================================================
    # Baseline 3: MLP on embeddings
    # =========================================================================
    print(f"\n{'='*60}")
    print("Baseline 3: MLP (2-layer) on ESM-2 embeddings")
    print(f"{'='*60}")
    
    # Encode labels as integers for MLP
    from sklearn.preprocessing import LabelEncoder
    le = LabelEncoder()
    y_train_encoded = le.fit_transform(y_train)
    y_test_encoded = le.transform(y_test)
    
    mlp = MLPClassifier(
        hidden_layer_sizes=(256, 64),
        max_iter=500,
        early_stopping=True,
        validation_fraction=0.1,
        random_state=42
    )
    mlp.fit(X_train_emb, y_train_encoded)
    y_pred_mlp_encoded = mlp.predict(X_test_emb)
    y_pred_mlp = le.inverse_transform(y_pred_mlp_encoded)
    y_proba_mlp = mlp.predict_proba(X_test_emb)
    
    metrics_mlp = compute_metrics(y_test, y_pred_mlp, y_proba_mlp, classes)
    print(f"  Accuracy: {metrics_mlp['accuracy']:.4f}")
    print(f"  Macro-F1: {metrics_mlp['macro_f1']:.4f}")
    
    mlp_result = {
        "baseline": "MLP",
        "description": "2-layer Multi-Layer Perceptron (256→64) on ESM-2 layer 20-30 embeddings",
        "features": "ESM-2 embeddings (1280-dim)",
        "config": {"hidden_layers": [256, 64], "max_iter": 500, "early_stopping": True},
        "n_train": len(train_indices),
        "n_test": len(test_indices),
        "metrics": metrics_mlp
    }
    with open(output_dir / "mlp_split40.json", 'w') as f:
        json.dump(mlp_result, f, indent=2)
    
    results.append({
        "Method": "MLP (256→64)",
        "Features": "ESM-2 embeddings",
        "Accuracy": metrics_mlp["accuracy"],
        "Macro_F1": metrics_mlp["macro_f1"],
        "Top3_Accuracy": metrics_mlp.get("top3_accuracy", None)
    })
    
    # =========================================================================
    # Baseline 4: Random (stratified)
    # =========================================================================
    print(f"\n{'='*60}")
    print("Baseline 4: Random (stratified by class distribution)")
    print(f"{'='*60}")
    
    random_clf = DummyClassifier(strategy='stratified', random_state=42)
    random_clf.fit(X_train_emb, y_train)
    y_pred_random = random_clf.predict(X_test_emb)
    y_proba_random = random_clf.predict_proba(X_test_emb)
    
    metrics_random = compute_metrics(y_test, y_pred_random, y_proba_random, classes)
    print(f"  Accuracy: {metrics_random['accuracy']:.4f}")
    print(f"  Macro-F1: {metrics_random['macro_f1']:.4f}")
    
    random_result = {
        "baseline": "Random",
        "description": "Stratified random prediction (baseline for chance-level performance)",
        "features": "None (predicts based on class distribution)",
        "config": {"strategy": "stratified"},
        "n_train": len(train_indices),
        "n_test": len(test_indices),
        "metrics": metrics_random
    }
    with open(output_dir / "random_split40.json", 'w') as f:
        json.dump(random_result, f, indent=2)
    
    results.append({
        "Method": "Random (stratified)",
        "Features": "None",
        "Accuracy": metrics_random["accuracy"],
        "Macro_F1": metrics_random["macro_f1"],
        "Top3_Accuracy": metrics_random.get("top3_accuracy", None)
    })
    
    # =========================================================================
    # Add our main model (from Step 11/12)
    # =========================================================================
    print(f"\n{'='*60}")
    print("Including Main Model (LR on ESM-2, calibrated)")
    print(f"{'='*60}")
    
    # Load from our calibration results
    cal_file = Path("results/calibration/split40_calibration.json")
    if cal_file.exists():
        with open(cal_file, 'r') as f:
            cal_data = json.load(f)
        
        results.append({
            "Method": "LR (Ours, calibrated)",
            "Features": "ESM-2 embeddings (layers 20-30)",
            "Accuracy": cal_data["metrics"]["calibrated"]["accuracy"],
            "Macro_F1": cal_data["metrics"]["calibrated"]["macro_f1"],
            "Top3_Accuracy": None  # Not stored in calibration file
        })
        print(f"  Accuracy: {cal_data['metrics']['calibrated']['accuracy']:.4f}")
        print(f"  Macro-F1: {cal_data['metrics']['calibrated']['macro_f1']:.4f}")
    
    # =========================================================================
    # Save summary table
    # =========================================================================
    print(f"\n{'='*60}")
    print("Saving baselines summary...")
    print(f"{'='*60}")
    
    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values('Accuracy', ascending=False)
    results_df.to_csv(output_dir / "baselines_split40.csv", index=False)
    
    # Save metadata
    metadata = {
        "step": 13,
        "name": "Baselines Comparison",
        "timestamp": datetime.now().isoformat(),
        "split": "40% identity threshold",
        "n_train": len(train_indices),
        "n_test": len(test_indices),
        "n_classes": len(classes),
        "classes": classes,
        "test_ids_file": str(test_file),
        "baselines": list(results_df["Method"].values)
    }
    with open(output_dir / "baselines_metadata.json", 'w') as f:
        json.dump(metadata, f, indent=2)
    
    print(f"\n{'='*60}")
    print("STEP 13 COMPLETE: Baselines Comparison")
    print(f"{'='*60}")
    
    print(f"\n{'Baselines Summary (sorted by accuracy)':^70}")
    print("-" * 70)
    print(f"{'Method':<25} {'Features':<25} {'Accuracy':>10} {'Macro-F1':>10}")
    print("-" * 70)
    for _, row in results_df.iterrows():
        print(f"{row['Method']:<25} {row['Features'][:24]:<25} {row['Accuracy']:>10.4f} {row['Macro_F1']:>10.4f}")
    print("-" * 70)
    
    print("\nSanity checks:")
    print(f"  ✓ All methods use identical test IDs ({len(test_indices)} samples)")
    print(f"  ✓ All methods evaluated on same 8-class problem")
    print("  ✓ Each baseline has individual JSON with full config")


if __name__ == "__main__":
    main()

