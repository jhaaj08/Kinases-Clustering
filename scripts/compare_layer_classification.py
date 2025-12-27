#!/usr/bin/env python3
"""
Compare Layer 33 vs Layers 20-30 for Classification

This script compares classification performance between:
- Layer 33 (final layer) embeddings
- Layers 20-30 (mid-layer) embeddings

For each configuration, it:
1. Trains logistic regression on homology-aware splits (40%, 50%, 70%)
2. Applies Platt scaling calibration
3. Computes comprehensive metrics

Usage:
    python scripts/compare_layer_classification.py
"""

import json
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime
from sklearn.linear_model import LogisticRegression
from sklearn.calibration import CalibratedClassifierCV
from sklearn.metrics import (
    accuracy_score,
    f1_score,
    log_loss,
    brier_score_loss
)
from sklearn.preprocessing import LabelEncoder
import joblib

# Configuration
EMBEDDING_CONFIGS = {
    "layer33": "domain_E001_layer33_mean.npy",
    "layers20_30": "domain_E001_layers20_30_mean.npy"
}

IDENTITY_THRESHOLDS = [40, 50, 70]
RANDOM_STATE = 42

# Paths
EMBEDDINGS_DIR = Path("embeddings/esm2_t33_650M")
SPLITS_DIR = Path("data/splits")
LABELS_FILE = Path("data/processed/labels.csv")
RESULTS_DIR = Path("results/layer_comparison")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def load_data():
    """Load embeddings, IDs, and labels."""
    # Load IDs
    with open(EMBEDDINGS_DIR / "ids.txt", 'r') as f:
        embedding_ids = [line.strip() for line in f if line.strip()]
    
    # Load labels
    labels_df = pd.read_csv(LABELS_FILE)
    id_to_label = dict(zip(labels_df['uniprot_id'], labels_df['label_used_for_experiments']))
    
    # Load embeddings
    embeddings = {}
    for name, filename in EMBEDDING_CONFIGS.items():
        embeddings[name] = np.load(EMBEDDINGS_DIR / filename)
        print(f"Loaded {name}: {embeddings[name].shape}")
    
    return embeddings, embedding_ids, id_to_label


def load_split(threshold):
    """Load train/test split for a given identity threshold."""
    train_file = SPLITS_DIR / f"split{threshold}_train.txt"
    test_file = SPLITS_DIR / f"split{threshold}_test.txt"
    
    with open(train_file, 'r') as f:
        train_ids = [line.strip() for line in f if line.strip()]
    with open(test_file, 'r') as f:
        test_ids = [line.strip() for line in f if line.strip()]
    
    return train_ids, test_ids


def get_data_for_split(train_ids, test_ids, embeddings, embedding_ids, id_to_label):
    """Get X, y data for a split."""
    id_to_idx = {uid: i for i, uid in enumerate(embedding_ids)}
    
    # Get training data
    train_indices = [id_to_idx[uid] for uid in train_ids if uid in id_to_idx and uid in id_to_label]
    train_labels = [id_to_label[uid] for uid in train_ids if uid in id_to_idx and uid in id_to_label]
    
    # Get test data
    test_indices = [id_to_idx[uid] for uid in test_ids if uid in id_to_idx and uid in id_to_label]
    test_labels = [id_to_label[uid] for uid in test_ids if uid in id_to_idx and uid in id_to_label]
    
    X_train = embeddings[train_indices]
    X_test = embeddings[test_indices]
    
    # Encode labels
    le = LabelEncoder()
    all_labels = sorted(list(set(train_labels + test_labels)))
    le.fit(all_labels)
    
    y_train = le.transform(train_labels)
    y_test = le.transform(test_labels)
    
    return X_train, X_test, y_train, y_test, le, all_labels


def compute_ece(y_true, y_proba, n_bins=10):
    """Compute Expected Calibration Error."""
    confidences = np.max(y_proba, axis=1)
    predictions = np.argmax(y_proba, axis=1)
    accuracies = (predictions == y_true).astype(float)
    
    bin_boundaries = np.linspace(0, 1, n_bins + 1)
    ece = 0.0
    
    for i in range(n_bins):
        in_bin = (confidences > bin_boundaries[i]) & (confidences <= bin_boundaries[i + 1])
        prop_in_bin = in_bin.mean()
        
        if prop_in_bin > 0:
            avg_confidence = confidences[in_bin].mean()
            avg_accuracy = accuracies[in_bin].mean()
            ece += np.abs(avg_confidence - avg_accuracy) * prop_in_bin
    
    return ece


def train_and_evaluate(X_train, X_test, y_train, y_test, classes, calibrate=False):
    """Train logistic regression and evaluate."""
    # Train base model
    model = LogisticRegression(
        solver='lbfgs',
        max_iter=1000,
        class_weight='balanced',
        multi_class='multinomial',
        random_state=RANDOM_STATE,
        C=1.0
    )
    model.fit(X_train, y_train)
    
    if calibrate:
        # Apply Platt scaling calibration
        calibrator = CalibratedClassifierCV(model, method='sigmoid', cv=5)
        calibrator.fit(X_train, y_train)
        y_pred = calibrator.predict(X_test)
        y_proba = calibrator.predict_proba(X_test)
    else:
        y_pred = model.predict(X_test)
        y_proba = model.predict_proba(X_test)
    
    # Compute metrics
    metrics = {
        "accuracy": float(accuracy_score(y_test, y_pred)),
        "macro_f1": float(f1_score(y_test, y_pred, average='macro', zero_division=0)),
        "weighted_f1": float(f1_score(y_test, y_pred, average='weighted', zero_division=0)),
        "log_loss": float(log_loss(y_test, y_proba, labels=range(len(classes)))),
        "ece": float(compute_ece(y_test, y_proba)),
    }
    
    return metrics, model if not calibrate else calibrator


def main():
    print("=" * 70)
    print("LAYER COMPARISON: Classification Performance")
    print("Layer 33 (Final) vs Layers 20-30 (Mid)")
    print("=" * 70)
    
    # Load data
    print("\nLoading data...")
    embeddings, embedding_ids, id_to_label = load_data()
    
    # Results storage
    all_results = []
    
    # Run experiments
    for threshold in IDENTITY_THRESHOLDS:
        print(f"\n{'=' * 70}")
        print(f"Identity Threshold: {threshold}%")
        print("=" * 70)
        
        # Load split
        train_ids, test_ids = load_split(threshold)
        print(f"Split: {len(train_ids)} train, {len(test_ids)} test")
        
        for config_name, embedding_file in EMBEDDING_CONFIGS.items():
            emb = embeddings[config_name]
            
            # Get data
            X_train, X_test, y_train, y_test, le, classes = get_data_for_split(
                train_ids, test_ids, emb, embedding_ids, id_to_label
            )
            
            print(f"\n--- {config_name} ---")
            print(f"  Training: {X_train.shape[0]} samples, {X_train.shape[1]} features")
            print(f"  Testing: {X_test.shape[0]} samples")
            print(f"  Classes: {len(classes)}")
            
            # Uncalibrated
            uncal_metrics, _ = train_and_evaluate(
                X_train, X_test, y_train, y_test, classes, calibrate=False
            )
            
            # Calibrated
            cal_metrics, _ = train_and_evaluate(
                X_train, X_test, y_train, y_test, classes, calibrate=True
            )
            
            print(f"  Uncalibrated: Acc={uncal_metrics['accuracy']:.3f}, F1={uncal_metrics['macro_f1']:.3f}, ECE={uncal_metrics['ece']:.3f}")
            print(f"  Calibrated:   Acc={cal_metrics['accuracy']:.3f}, F1={cal_metrics['macro_f1']:.3f}, ECE={cal_metrics['ece']:.3f}")
            
            # Store results
            all_results.append({
                "threshold": threshold,
                "config": config_name,
                "calibrated": False,
                **uncal_metrics
            })
            all_results.append({
                "threshold": threshold,
                "config": config_name,
                "calibrated": True,
                **cal_metrics
            })
    
    # Create comparison table
    print("\n" + "=" * 70)
    print("COMPARISON SUMMARY")
    print("=" * 70)
    
    df = pd.DataFrame(all_results)
    
    # Pivot for comparison
    print("\n### Uncalibrated Accuracy by Configuration and Threshold ###")
    uncal_df = df[df['calibrated'] == False].pivot(index='threshold', columns='config', values='accuracy')
    uncal_df['delta'] = uncal_df['layers20_30'] - uncal_df['layer33']
    uncal_df['delta_pct'] = (uncal_df['delta'] / uncal_df['layer33'] * 100).round(1)
    print(uncal_df.to_string())
    
    print("\n### Calibrated Accuracy by Configuration and Threshold ###")
    cal_df = df[df['calibrated'] == True].pivot(index='threshold', columns='config', values='accuracy')
    cal_df['delta'] = cal_df['layers20_30'] - cal_df['layer33']
    cal_df['delta_pct'] = (cal_df['delta'] / cal_df['layer33'] * 100).round(1)
    print(cal_df.to_string())
    
    print("\n### Macro-F1 (Calibrated) by Configuration and Threshold ###")
    f1_df = df[df['calibrated'] == True].pivot(index='threshold', columns='config', values='macro_f1')
    f1_df['delta'] = f1_df['layers20_30'] - f1_df['layer33']
    print(f1_df.to_string())
    
    print("\n### ECE (Calibrated) by Configuration and Threshold ###")
    ece_df = df[df['calibrated'] == True].pivot(index='threshold', columns='config', values='ece')
    ece_df['delta'] = ece_df['layers20_30'] - ece_df['layer33']
    print(ece_df.to_string())
    
    # Save full results
    results_file = RESULTS_DIR / "layer_comparison_results.json"
    with open(results_file, 'w') as f:
        json.dump({
            "timestamp": datetime.now().isoformat(),
            "configurations": list(EMBEDDING_CONFIGS.keys()),
            "thresholds": IDENTITY_THRESHOLDS,
            "results": all_results
        }, f, indent=2)
    print(f"\nFull results saved to: {results_file}")
    
    # Save summary CSV
    summary_file = RESULTS_DIR / "layer_comparison_summary.csv"
    df.to_csv(summary_file, index=False)
    print(f"Summary CSV saved to: {summary_file}")
    
    # Print final comparison
    print("\n" + "=" * 70)
    print("FINAL COMPARISON (40% threshold - strictest)")
    print("=" * 70)
    
    layer33_uncal = df[(df['threshold'] == 40) & (df['config'] == 'layer33') & (df['calibrated'] == False)].iloc[0]
    layer33_cal = df[(df['threshold'] == 40) & (df['config'] == 'layer33') & (df['calibrated'] == True)].iloc[0]
    layers2030_uncal = df[(df['threshold'] == 40) & (df['config'] == 'layers20_30') & (df['calibrated'] == False)].iloc[0]
    layers2030_cal = df[(df['threshold'] == 40) & (df['config'] == 'layers20_30') & (df['calibrated'] == True)].iloc[0]
    
    print(f"\n| Metric | Layer 33 | Layers 20-30 | Δ | Δ% |")
    print(f"|--------|----------|--------------|---|-----|")
    
    acc_delta = layers2030_cal['accuracy'] - layer33_cal['accuracy']
    acc_pct = acc_delta / layer33_cal['accuracy'] * 100
    print(f"| Accuracy (cal) | {layer33_cal['accuracy']:.3f} | {layers2030_cal['accuracy']:.3f} | {acc_delta:+.3f} | {acc_pct:+.1f}% |")
    
    f1_delta = layers2030_cal['macro_f1'] - layer33_cal['macro_f1']
    f1_pct = f1_delta / layer33_cal['macro_f1'] * 100
    print(f"| Macro-F1 (cal) | {layer33_cal['macro_f1']:.3f} | {layers2030_cal['macro_f1']:.3f} | {f1_delta:+.3f} | {f1_pct:+.1f}% |")
    
    ll_delta = layers2030_cal['log_loss'] - layer33_cal['log_loss']
    ll_pct = ll_delta / layer33_cal['log_loss'] * 100
    print(f"| Log-loss (cal) | {layer33_cal['log_loss']:.3f} | {layers2030_cal['log_loss']:.3f} | {ll_delta:+.3f} | {ll_pct:+.1f}% |")
    
    ece_delta = layers2030_cal['ece'] - layer33_cal['ece']
    ece_pct = ece_delta / layer33_cal['ece'] * 100
    print(f"| ECE (cal) | {layer33_cal['ece']:.3f} | {layers2030_cal['ece']:.3f} | {ece_delta:+.3f} | {ece_pct:+.1f}% |")
    
    print("\n" + "=" * 70)
    print("EXPERIMENT COMPLETE")
    print("=" * 70)


if __name__ == "__main__":
    main()

