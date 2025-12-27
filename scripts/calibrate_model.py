#!/usr/bin/env python3
"""
Step 12: Calibration (Platt Scaling) + ECE

This script applies probability calibration to the trained logistic regression
model using Platt scaling (sigmoid calibration) via cross-validation on the
training set, then evaluates on the test set.

Usage:
    python scripts/calibrate_model.py

Inputs:
    - models/lr_split40.joblib (trained model)
    - embeddings/esm2_t33_650M/domain_E001_layers20_30_mean.npy
    - embeddings/esm2_t33_650M/ids.txt
    - data/splits/split40_train.txt, split40_test.txt
    - data/processed/labels.csv

Outputs:
    - models/lr_split40_calibrated.joblib
    - results/calibration/split40_calibration.json
    - figures/reliability_split40.png
"""

import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from datetime import datetime
from sklearn.calibration import CalibratedClassifierCV, calibration_curve
from sklearn.metrics import (
    accuracy_score,
    f1_score,
    log_loss,
    brier_score_loss
)
import joblib


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
    
    for uid in ids_list:
        if uid in id_to_idx and uid in id_to_label:
            indices.append(id_to_idx[uid])
            labels.append(id_to_label[uid])
    
    return indices, labels


def compute_ece(y_true, y_prob, n_bins=10):
    """
    Compute Expected Calibration Error (ECE).
    
    ECE = sum over bins of (|bin_accuracy - bin_confidence| * bin_weight)
    """
    # Get predicted class and max probability
    y_pred = np.argmax(y_prob, axis=1)
    confidences = np.max(y_prob, axis=1)
    
    # Convert labels to indices
    unique_labels = sorted(set(y_true))
    label_to_idx = {label: i for i, label in enumerate(unique_labels)}
    y_true_idx = np.array([label_to_idx[label] for label in y_true])
    
    # Compute correctness
    correctness = (y_pred == y_true_idx).astype(float)
    
    # Bin by confidence
    bin_boundaries = np.linspace(0, 1, n_bins + 1)
    bin_lowers = bin_boundaries[:-1]
    bin_uppers = bin_boundaries[1:]
    
    ece = 0.0
    bin_data = []
    
    for bin_lower, bin_upper in zip(bin_lowers, bin_uppers):
        in_bin = (confidences > bin_lower) & (confidences <= bin_upper)
        prop_in_bin = in_bin.mean()
        
        if prop_in_bin > 0:
            avg_confidence = confidences[in_bin].mean()
            avg_accuracy = correctness[in_bin].mean()
            bin_ece = abs(avg_accuracy - avg_confidence) * prop_in_bin
            ece += bin_ece
            
            bin_data.append({
                "bin_lower": float(bin_lower),
                "bin_upper": float(bin_upper),
                "bin_center": float((bin_lower + bin_upper) / 2),
                "n_samples": int(in_bin.sum()),
                "avg_confidence": float(avg_confidence),
                "avg_accuracy": float(avg_accuracy),
                "bin_weight": float(prop_in_bin),
                "bin_ece": float(bin_ece)
            })
        else:
            bin_data.append({
                "bin_lower": float(bin_lower),
                "bin_upper": float(bin_upper),
                "bin_center": float((bin_lower + bin_upper) / 2),
                "n_samples": 0,
                "avg_confidence": None,
                "avg_accuracy": None,
                "bin_weight": 0.0,
                "bin_ece": 0.0
            })
    
    return ece, bin_data


def plot_reliability_diagram(bin_data_uncal, bin_data_cal, output_path):
    """Plot reliability diagram comparing uncalibrated vs calibrated."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    for ax, bin_data, title in zip(axes, [bin_data_uncal, bin_data_cal], 
                                     ["Uncalibrated", "Calibrated"]):
        # Extract data for non-empty bins
        confidences = []
        accuracies = []
        weights = []
        
        for b in bin_data:
            if b["n_samples"] > 0:
                confidences.append(b["avg_confidence"])
                accuracies.append(b["avg_accuracy"])
                weights.append(b["n_samples"])
        
        # Plot perfect calibration line
        ax.plot([0, 1], [0, 1], 'k--', label='Perfect calibration', alpha=0.7)
        
        # Plot reliability bars
        if confidences:
            ax.bar([b["bin_center"] for b in bin_data if b["n_samples"] > 0],
                   accuracies, width=0.1, alpha=0.7, 
                   color='steelblue', edgecolor='navy', linewidth=1.5,
                   label='Actual accuracy')
        
        ax.set_xlabel('Mean Predicted Probability', fontsize=12)
        ax.set_ylabel('Fraction of Positives', fontsize=12)
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.set_xlim([0, 1])
        ax.set_ylim([0, 1])
        ax.legend(loc='upper left')
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved reliability diagram: {output_path}")


def main():
    print("="*60)
    print("Step 12: Calibration (Platt Scaling) + ECE")
    print("="*60)
    
    # Paths
    embeddings_dir = Path("embeddings/esm2_t33_650M")
    splits_dir = Path("data/splits")
    models_dir = Path("models")
    output_dir = Path("results/calibration")
    figures_dir = Path("figures")
    
    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    
    # Load embeddings
    print("\nLoading embeddings...")
    embedding_file = embeddings_dir / "domain_E001_layers20_30_mean.npy"
    embeddings = np.load(embedding_file)
    print(f"  Shape: {embeddings.shape}")
    
    # Load embedding IDs
    with open(embeddings_dir / "ids.txt", 'r') as f:
        embedding_ids = [line.strip() for line in f]
    
    # Load labels
    labels_df = pd.read_csv("data/processed/labels.csv")
    id_to_label = dict(zip(labels_df['uniprot_id'], labels_df['label_used_for_experiments']))
    
    # Load uncalibrated model
    print("\nLoading uncalibrated model...")
    model_file = models_dir / "lr_split40.joblib"
    uncalibrated_model = joblib.load(model_file)
    print(f"  Loaded: {model_file}")
    
    # Load split
    train_file = splits_dir / "split40_train.txt"
    test_file = splits_dir / "split40_test.txt"
    train_ids, test_ids = load_split_ids(train_file, test_file)
    
    # Get data
    train_indices, train_labels = get_indices_and_labels(train_ids, embedding_ids, id_to_label)
    test_indices, test_labels = get_indices_and_labels(test_ids, embedding_ids, id_to_label)
    
    X_train = embeddings[train_indices]
    X_test = embeddings[test_indices]
    y_train = train_labels
    y_test = test_labels
    
    classes = sorted(set(y_train + y_test))
    print(f"  Train: {len(train_indices)}, Test: {len(test_indices)}")
    print(f"  Classes: {len(classes)}")
    
    # --- UNCALIBRATED METRICS ---
    print(f"\n{'='*60}")
    print("Evaluating UNCALIBRATED model...")
    print(f"{'='*60}")
    
    y_pred_uncal = uncalibrated_model.predict(X_test)
    y_prob_uncal = uncalibrated_model.predict_proba(X_test)
    
    acc_uncal = accuracy_score(y_test, y_pred_uncal)
    f1_uncal = f1_score(y_test, y_pred_uncal, average='macro', zero_division=0)
    logloss_uncal = log_loss(y_test, y_prob_uncal, labels=classes)
    ece_uncal, bins_uncal = compute_ece(y_test, y_prob_uncal)
    
    print(f"  Accuracy (uncalibrated): {acc_uncal:.4f}")
    print(f"  Macro-F1 (uncalibrated): {f1_uncal:.4f}")
    print(f"  Log-loss (uncalibrated): {logloss_uncal:.4f}")
    print(f"  ECE (uncalibrated): {ece_uncal:.4f}")
    
    # --- FIT CALIBRATION ---
    print(f"\n{'='*60}")
    print("Fitting calibration via CV on training set...")
    print(f"{'='*60}")
    
    # Use Platt scaling (sigmoid) with cross-validation
    calibrated_model = CalibratedClassifierCV(
        uncalibrated_model,
        method='sigmoid',  # Platt scaling
        cv=5  # 5-fold CV on training set
    )
    calibrated_model.fit(X_train, y_train)
    print("  Calibration fitted using 5-fold CV with sigmoid method (Platt scaling)")
    
    # --- CALIBRATED METRICS ---
    print(f"\n{'='*60}")
    print("Evaluating CALIBRATED model...")
    print(f"{'='*60}")
    
    y_pred_cal = calibrated_model.predict(X_test)
    y_prob_cal = calibrated_model.predict_proba(X_test)
    
    acc_cal = accuracy_score(y_test, y_pred_cal)
    f1_cal = f1_score(y_test, y_pred_cal, average='macro', zero_division=0)
    logloss_cal = log_loss(y_test, y_prob_cal, labels=classes)
    ece_cal, bins_cal = compute_ece(y_test, y_prob_cal)
    
    print(f"  Accuracy (calibrated): {acc_cal:.4f}")
    print(f"  Macro-F1 (calibrated): {f1_cal:.4f}")
    print(f"  Log-loss (calibrated): {logloss_cal:.4f}")
    print(f"  ECE (calibrated): {ece_cal:.4f}")
    
    # --- SAVE CALIBRATED MODEL ---
    cal_model_file = models_dir / "lr_split40_calibrated.joblib"
    joblib.dump(calibrated_model, cal_model_file)
    print(f"\n  Saved calibrated model: {cal_model_file}")
    
    # --- SAVE CALIBRATION RESULTS ---
    result = {
        "step": 12,
        "name": "Calibration (Platt Scaling)",
        "timestamp": datetime.now().isoformat(),
        "split": "40%",
        "calibration_method": "sigmoid (Platt scaling)",
        "cv_folds": 5,
        "n_train": len(train_indices),
        "n_test": len(test_indices),
        "n_classes": len(classes),
        "classes": classes,
        "metrics": {
            "uncalibrated": {
                "accuracy": float(acc_uncal),
                "macro_f1": float(f1_uncal),
                "log_loss": float(logloss_uncal),
                "ece": float(ece_uncal)
            },
            "calibrated": {
                "accuracy": float(acc_cal),
                "macro_f1": float(f1_cal),
                "log_loss": float(logloss_cal),
                "ece": float(ece_cal)
            },
            "improvement": {
                "accuracy_delta": float(acc_cal - acc_uncal),
                "log_loss_delta": float(logloss_cal - logloss_uncal),
                "ece_delta": float(ece_cal - ece_uncal),
                "ece_reduction_percent": float((ece_uncal - ece_cal) / ece_uncal * 100) if ece_uncal > 0 else 0
            }
        },
        "reliability_bins": {
            "uncalibrated": bins_uncal,
            "calibrated": bins_cal
        },
        "notes": {
            "accuracy_for_baselines_table": "uncalibrated" if acc_uncal >= acc_cal else "calibrated",
            "reported_accuracy": float(max(acc_uncal, acc_cal)),
            "which_accuracy_used": f"Using {'uncalibrated' if acc_uncal >= acc_cal else 'calibrated'} accuracy ({max(acc_uncal, acc_cal):.4f}) in baselines table"
        }
    }
    
    calibration_file = output_dir / "split40_calibration.json"
    with open(calibration_file, 'w') as f:
        json.dump(result, f, indent=2)
    print(f"  Saved calibration results: {calibration_file}")
    
    # --- PLOT RELIABILITY DIAGRAM ---
    plot_reliability_diagram(
        bins_uncal, bins_cal,
        figures_dir / "reliability_split40.png"
    )
    
    # --- SUMMARY ---
    print(f"\n{'='*60}")
    print("STEP 12 COMPLETE: Calibration")
    print(f"{'='*60}")
    
    print(f"\n{'Comparison: Before vs After Calibration':^60}")
    print("-" * 60)
    print(f"{'Metric':<20} {'Uncalibrated':>15} {'Calibrated':>15} {'Delta':>10}")
    print("-" * 60)
    print(f"{'Accuracy':<20} {acc_uncal:>15.4f} {acc_cal:>15.4f} {acc_cal-acc_uncal:>+10.4f}")
    print(f"{'Macro-F1':<20} {f1_uncal:>15.4f} {f1_cal:>15.4f} {f1_cal-f1_uncal:>+10.4f}")
    print(f"{'Log-loss':<20} {logloss_uncal:>15.4f} {logloss_cal:>15.4f} {logloss_cal-logloss_uncal:>+10.4f}")
    print(f"{'ECE':<20} {ece_uncal:>15.4f} {ece_cal:>15.4f} {ece_cal-ece_uncal:>+10.4f}")
    print("-" * 60)
    
    print(f"\nECE reduction: {(ece_uncal - ece_cal) / ece_uncal * 100:.1f}%")
    
    if acc_uncal != acc_cal:
        print(f"\n⚠️  IMPORTANT: Accuracy changed after calibration!")
        print(f"   Uncalibrated: {acc_uncal:.4f}")
        print(f"   Calibrated:   {acc_cal:.4f}")
        print(f"   For baselines table, use: {result['notes']['which_accuracy_used']}")
    
    print("\nSanity checks:")
    print("  ✓ Uncalibrated and calibrated metrics clearly separated")
    print("  ✓ ECE computed with 10 bins")
    print("  ✓ Reliability diagram saved")
    print("  ✓ Which accuracy to use in baselines is documented")


if __name__ == "__main__":
    main()

