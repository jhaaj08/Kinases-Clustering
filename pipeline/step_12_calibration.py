#!/usr/bin/env python3
"""
Step 12: Model Calibration (Platt Scaling)

This script applies Platt scaling to improve probability calibration
and computes Expected Calibration Error (ECE).

Usage:
    python pipeline/step_12_calibration.py --run-dir runs/2025-01-01_000000/

Outputs:
    - results/calibration/split40_calibration.json
    - models/lr_split40_calibrated.joblib
"""

import argparse
import json
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime
from sklearn.linear_model import LogisticRegression
from sklearn.calibration import CalibratedClassifierCV
from sklearn.metrics import accuracy_score, log_loss
import joblib

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from pipeline.membership import load_split

# Project root
PROJECT_ROOT = Path(__file__).parent.parent

# Configuration
RANDOM_STATE = 42
N_BINS = 10


def compute_ece(y_true, y_proba, n_bins=10):
    """Compute Expected Calibration Error."""
    confidences = np.max(y_proba, axis=1)
    predictions = np.argmax(y_proba, axis=1)
    
    # Convert labels to indices if needed
    if isinstance(y_true[0], str):
        unique_labels = sorted(set(y_true))
        label_to_idx = {l: i for i, l in enumerate(unique_labels)}
        y_true_idx = np.array([label_to_idx[l] for l in y_true])
    else:
        y_true_idx = np.array(y_true)
    
    accuracies = (predictions == y_true_idx)
    
    bin_boundaries = np.linspace(0, 1, n_bins + 1)
    ece = 0.0
    
    for i in range(n_bins):
        in_bin = (confidences > bin_boundaries[i]) & (confidences <= bin_boundaries[i + 1])
        prop_in_bin = in_bin.mean()
        
        if prop_in_bin > 0:
            avg_confidence = confidences[in_bin].mean()
            avg_accuracy = accuracies[in_bin].mean()
            ece += np.abs(avg_accuracy - avg_confidence) * prop_in_bin
    
    return float(ece)


def main():
    parser = argparse.ArgumentParser(description="Apply model calibration")
    parser.add_argument("--run-dir", type=str, required=True, help="Run directory path")
    args = parser.parse_args()
    
    run_dir = Path(args.run_dir)
    results_dir = run_dir / "results" / "calibration"
    results_dir.mkdir(parents=True, exist_ok=True)
    models_dir = run_dir / "models"
    models_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 60)
    print("Step 12: Model Calibration")
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
    
    # Load embeddings (use layer33_mean for calibration)
    emb_file = run_dir / "embeddings" / "esm2_t33_650M" / "domain_E001_layer33_mean.npy"
    embeddings = np.load(emb_file)
    
    # Use split40 for calibration
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
    
    # Train uncalibrated model
    print("\nTraining uncalibrated model...")
    base_model = LogisticRegression(
        max_iter=1000,
        random_state=RANDOM_STATE,
        multi_class='multinomial',
        solver='lbfgs'
    )
    base_model.fit(X_train, train_labels)
    
    # Uncalibrated predictions
    y_pred_uncal = base_model.predict(X_test)
    y_proba_uncal = base_model.predict_proba(X_test)
    
    acc_uncal = accuracy_score(test_labels, y_pred_uncal)
    ll_uncal = log_loss(test_labels, y_proba_uncal, labels=base_model.classes_)
    ece_uncal = compute_ece(test_labels, y_proba_uncal)
    
    print(f"  Accuracy: {acc_uncal:.4f}")
    print(f"  Log-loss: {ll_uncal:.4f}")
    print(f"  ECE: {ece_uncal:.4f}")
    
    # Train calibrated model (Platt scaling via sigmoid)
    print("\nTraining calibrated model (Platt scaling)...")
    calibrated_model = CalibratedClassifierCV(
        base_model,
        method='sigmoid',
        cv=5
    )
    calibrated_model.fit(X_train, train_labels)
    
    # Calibrated predictions
    y_pred_cal = calibrated_model.predict(X_test)
    y_proba_cal = calibrated_model.predict_proba(X_test)
    
    acc_cal = accuracy_score(test_labels, y_pred_cal)
    ll_cal = log_loss(test_labels, y_proba_cal, labels=calibrated_model.classes_)
    ece_cal = compute_ece(test_labels, y_proba_cal)
    
    print(f"  Accuracy: {acc_cal:.4f}")
    print(f"  Log-loss: {ll_cal:.4f}")
    print(f"  ECE: {ece_cal:.4f}")
    
    # Save calibrated model
    model_file = models_dir / f"lr_{split_name}_calibrated.joblib"
    joblib.dump(calibrated_model, model_file)
    print(f"\n✓ Saved model: {model_file}")
    
    # Create calibration report
    report = {
        "step": 12,
        "name": "Model Calibration",
        "timestamp": datetime.now().isoformat(),
        "split": split_name,
        "method": "Platt scaling (sigmoid)",
        "n_cv_folds": 5,
        "n_bins_ece": N_BINS,
        "uncalibrated": {
            "accuracy": float(acc_uncal),
            "log_loss": float(ll_uncal),
            "ece": float(ece_uncal)
        },
        "calibrated": {
            "accuracy": float(acc_cal),
            "log_loss": float(ll_cal),
            "ece": float(ece_cal)
        },
        "improvement": {
            "accuracy_delta": float(acc_cal - acc_uncal),
            "log_loss_reduction_pct": float((ll_uncal - ll_cal) / ll_uncal * 100),
            "ece_reduction_pct": float((ece_uncal - ece_cal) / ece_uncal * 100) if ece_uncal > 0 else 0
        }
    }
    
    report_file = results_dir / f"{split_name}_calibration.json"
    with open(report_file, 'w') as f:
        json.dump(report, f, indent=2)
    print(f"✓ Saved: {report_file}")
    
    print("\n" + "=" * 60)
    print("Step 12 COMPLETE")
    print("=" * 60)
    print("\n                Calibration Summary")
    print("-" * 60)
    print(f"{'Metric':<20} {'Uncalibrated':<15} {'Calibrated':<15} {'Change':<15}")
    print("-" * 60)
    print(f"{'Accuracy':<20} {acc_uncal:<15.4f} {acc_cal:<15.4f} {acc_cal - acc_uncal:+.4f}")
    print(f"{'Log-loss':<20} {ll_uncal:<15.4f} {ll_cal:<15.4f} {ll_cal - ll_uncal:+.4f}")
    print(f"{'ECE':<20} {ece_uncal:<15.4f} {ece_cal:<15.4f} {ece_cal - ece_uncal:+.4f}")
    print("-" * 60)


if __name__ == "__main__":
    main()

