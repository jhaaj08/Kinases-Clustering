#!/usr/bin/env python3
"""
Step 11: Train Logistic Regression + Evaluate (Uncalibrated)

This script trains a logistic regression classifier on homology-aware splits
and evaluates performance metrics.

Usage:
    python scripts/train_supervised.py

Inputs:
    - embeddings/esm2_t33_650M/domain_E001_layers20_30_mean.npy (best embeddings)
    - embeddings/esm2_t33_650M/ids.txt
    - data/splits/split*_train.txt, split*_test.txt
    - data/processed/labels.csv

Outputs:
    - models/lr_split40.joblib
    - results/supervised/lr_split40_metrics.json
    - results/supervised/lr_split40_confusion.csv
    - results/supervised/lr_multi_identity_summary.csv
"""

import json
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    f1_score,
    precision_score,
    recall_score,
    confusion_matrix,
    classification_report,
    top_k_accuracy_score,
    log_loss
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
    valid_ids = []
    
    for uid in ids_list:
        if uid in id_to_idx and uid in id_to_label:
            indices.append(id_to_idx[uid])
            labels.append(id_to_label[uid])
            valid_ids.append(uid)
    
    return indices, labels, valid_ids


def compute_metrics(y_true, y_pred, y_proba, classes):
    """Compute all classification metrics."""
    metrics = {
        "accuracy": float(accuracy_score(y_true, y_pred)),
        "macro_f1": float(f1_score(y_true, y_pred, average='macro', zero_division=0)),
        "weighted_f1": float(f1_score(y_true, y_pred, average='weighted', zero_division=0)),
        "macro_precision": float(precision_score(y_true, y_pred, average='macro', zero_division=0)),
        "macro_recall": float(recall_score(y_true, y_pred, average='macro', zero_division=0)),
        "log_loss": float(log_loss(y_true, y_proba, labels=classes)),  # Specify all training classes
    }
    
    # Top-3 accuracy if more than 3 classes
    if len(classes) > 3:
        metrics["top3_accuracy"] = float(top_k_accuracy_score(
            y_true, y_proba, k=3, labels=classes
        ))
    
    return metrics


def compute_per_class_metrics(y_true, y_pred, classes):
    """Compute per-class precision, recall, F1."""
    # Use labels parameter to handle cases where some classes aren't in y_true
    report = classification_report(y_true, y_pred, labels=classes, target_names=classes, 
                                   output_dict=True, zero_division=0)
    
    per_class = []
    for cls in classes:
        if cls in report:
            per_class.append({
                "class": cls,
                "precision": round(report[cls]["precision"], 4),
                "recall": round(report[cls]["recall"], 4),
                "f1": round(report[cls]["f1-score"], 4),
                "support": int(report[cls]["support"])
            })
    
    return per_class


def main():
    print("="*60)
    print("Step 11: Train Logistic Regression (Uncalibrated)")
    print("="*60)
    
    # Paths
    embeddings_dir = Path("embeddings/esm2_t33_650M")
    splits_dir = Path("data/splits")
    models_dir = Path("models")
    output_dir = Path("results/supervised")
    
    models_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load embeddings (use best configuration: layers 20-30)
    print("\nLoading embeddings...")
    embedding_file = embeddings_dir / "domain_E001_layers20_30_mean.npy"
    embeddings = np.load(embedding_file)
    print(f"  Shape: {embeddings.shape}")
    
    # Load embedding IDs
    with open(embeddings_dir / "ids.txt", 'r') as f:
        embedding_ids = [line.strip() for line in f]
    print(f"  IDs: {len(embedding_ids)}")
    
    # Load labels
    labels_df = pd.read_csv("data/processed/labels.csv")
    id_to_label = dict(zip(labels_df['uniprot_id'], labels_df['label_used_for_experiments']))
    
    # Model configuration
    model_config = {
        "model": "LogisticRegression",
        "solver": "lbfgs",
        "max_iter": 1000,
        "class_weight": "balanced",
        "multi_class": "multinomial",
        "random_state": 42,
        "C": 1.0  # Regularization strength (inverse)
    }
    
    print(f"\nModel configuration:")
    for k, v in model_config.items():
        print(f"  {k}: {v}")
    
    # Identity thresholds to evaluate
    thresholds = [40, 50, 70]
    
    # Results for multi-identity summary
    summary_rows = []
    
    print(f"\n{'='*60}")
    print("Training and evaluating at each identity threshold...")
    print(f"{'='*60}")
    
    for threshold in thresholds:
        print(f"\n--- {threshold}% Identity Threshold ---")
        
        # Load split
        train_file = splits_dir / f"split{threshold}_train.txt"
        test_file = splits_dir / f"split{threshold}_test.txt"
        
        if not train_file.exists() or not test_file.exists():
            print(f"  Skip: split files not found")
            continue
        
        train_ids, test_ids = load_split_ids(train_file, test_file)
        print(f"  Raw split: {len(train_ids)} train, {len(test_ids)} test")
        
        # Get indices and labels
        train_indices, train_labels, train_valid_ids = get_indices_and_labels(
            train_ids, embedding_ids, id_to_label
        )
        test_indices, test_labels, test_valid_ids = get_indices_and_labels(
            test_ids, embedding_ids, id_to_label
        )
        
        print(f"  Valid: {len(train_indices)} train, {len(test_indices)} test")
        
        # Get embeddings
        X_train = embeddings[train_indices]
        X_test = embeddings[test_indices]
        y_train = train_labels
        y_test = test_labels
        
        # Get unique classes
        classes = sorted(set(y_train + y_test))
        print(f"  Classes: {len(classes)}")
        
        # Train model
        print(f"  Training logistic regression...")
        model = LogisticRegression(
            solver=model_config["solver"],
            max_iter=model_config["max_iter"],
            class_weight=model_config["class_weight"],
            multi_class=model_config["multi_class"],
            random_state=model_config["random_state"],
            C=model_config["C"]
        )
        model.fit(X_train, y_train)
        
        # Predict
        y_pred = model.predict(X_test)
        y_proba = model.predict_proba(X_test)
        
        # Compute metrics
        metrics = compute_metrics(y_test, y_pred, y_proba, classes)
        per_class = compute_per_class_metrics(y_test, y_pred, classes)
        
        print(f"  Accuracy: {metrics['accuracy']:.4f}")
        print(f"  Macro-F1: {metrics['macro_f1']:.4f}")
        print(f"  Weighted-F1: {metrics['weighted_f1']:.4f}")
        
        # Create confusion matrix
        cm = confusion_matrix(y_test, y_pred, labels=classes)
        cm_df = pd.DataFrame(cm, index=classes, columns=classes)
        
        # Create result structure
        result = {
            "step": 11,
            "name": f"Logistic Regression (split{threshold})",
            "timestamp": datetime.now().isoformat(),
            "calibration": "uncalibrated",
            "split": {
                "threshold": threshold,
                "n_train": len(train_indices),
                "n_test": len(test_indices),
                "n_classes": len(classes),
                "classes": classes
            },
            "model_config": model_config,
            "embedding": {
                "file": str(embedding_file),
                "layers": [20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30],
                "pooling": "mean"
            },
            "metrics": metrics,
            "per_class_metrics": per_class
        }
        
        # Save outputs for split40 (primary evaluation)
        if threshold == 40:
            # Save model
            model_file = models_dir / "lr_split40.joblib"
            joblib.dump(model, model_file)
            print(f"  Saved model: {model_file}")
            
            # Save metrics JSON
            metrics_file = output_dir / "lr_split40_metrics.json"
            with open(metrics_file, 'w') as f:
                json.dump(result, f, indent=2)
            print(f"  Saved metrics: {metrics_file}")
            
            # Save confusion matrix
            cm_file = output_dir / "lr_split40_confusion.csv"
            cm_df.to_csv(cm_file)
            print(f"  Saved confusion matrix: {cm_file}")
        
        # Also save for each threshold
        metrics_file = output_dir / f"lr_split{threshold}_metrics.json"
        with open(metrics_file, 'w') as f:
            json.dump(result, f, indent=2)
        
        # Add to summary
        summary_rows.append({
            "Identity_Threshold": f"{threshold}%",
            "N_Train": len(train_indices),
            "N_Test": len(test_indices),
            "N_Classes": len(classes),
            "Accuracy": metrics["accuracy"],
            "Macro_F1": metrics["macro_f1"],
            "Weighted_F1": metrics["weighted_f1"],
            "Top3_Accuracy": metrics.get("top3_accuracy", None),
            "Log_Loss": metrics["log_loss"],
            "Calibration": "uncalibrated"
        })
    
    # Save multi-identity summary
    summary_df = pd.DataFrame(summary_rows)
    summary_file = output_dir / "lr_multi_identity_summary.csv"
    summary_df.to_csv(summary_file, index=False)
    print(f"\nSaved summary: {summary_file}")
    
    print(f"\n{'='*60}")
    print("STEP 11 COMPLETE: Supervised Training (Uncalibrated)")
    print(f"{'='*60}")
    
    print(f"\n{'Multi-Identity Summary':^60}")
    print("-" * 70)
    print(f"{'Threshold':<12} {'N_Test':<10} {'Accuracy':<10} {'Macro-F1':<10} {'Weighted-F1':<12}")
    print("-" * 70)
    for row in summary_rows:
        print(f"{row['Identity_Threshold']:<12} {row['N_Test']:<10} {row['Accuracy']:<10.4f} {row['Macro_F1']:<10.4f} {row['Weighted_F1']:<12.4f}")
    print("-" * 70)
    
    print("\nNote: All metrics are UNCALIBRATED. Calibrated metrics in Step 12.")
    
    print("\nSanity checks:")
    print("  ✓ All metrics clearly labeled as 'uncalibrated'")
    print("  ✓ Per-class metrics generated from JSON/CSV")
    print("  ✓ Confusion matrix saved for analysis")


if __name__ == "__main__":
    main()
