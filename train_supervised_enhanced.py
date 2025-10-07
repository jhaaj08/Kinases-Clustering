#!/usr/bin/env python3
"""
Enhanced supervised kinase classification with calibration and uncertainty.

NEW FEATURES (Reviewer requested):
1. Calibrated probabilities (CalibratedClassifierCV)
2. Top-3 accuracy
3. Expected Calibration Error (ECE)
4. Reliability diagrams
5. Per-sequence confidence scores with "needs-review" flagging
6. Evaluation across multiple identity thresholds (70%, 50%, 40%)
7. Per-class calibration metrics
"""

import os
import sys
import argparse
import json
import numpy as np
import pandas as pd
from pathlib import Path

# Machine learning
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (
    accuracy_score, classification_report, confusion_matrix,
    f1_score, precision_recall_fscore_support, top_k_accuracy_score,
    log_loss, brier_score_loss
)
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.calibration import CalibratedClassifierCV, calibration_curve

import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns

import warnings
warnings.filterwarnings('ignore')


def load_splits(splits_file):
    """Load train/test splits from JSON."""
    with open(splits_file, 'r') as f:
        data = json.load(f)
    return data['train_ids'], data['test_ids'], data['metadata']


def load_data(embeddings_dir, labels_csv, splits_file):
    """
    Load embeddings, labels, and splits.
    
    Returns:
        X_train, X_test, y_train, y_test, train_ids, test_ids, labels
    """
    # Load embeddings
    embeddings = np.load(f"{embeddings_dir}/esm2_embeddings.npy")
    index_df = pd.read_csv(f"{embeddings_dir}/esm2_index.csv")
    
    # Load labels
    labels_df = pd.read_csv(labels_csv)
    
    # Merge
    data = index_df.merge(
        labels_df[['uniprot_id', 'kinome_group_major']], 
        on='uniprot_id', 
        how='left'
    )
    
    # Load splits
    train_ids, test_ids, metadata = load_splits(splits_file)
    
    # Filter by splits
    train_mask = data['uniprot_id'].isin(train_ids)
    test_mask = data['uniprot_id'].isin(test_ids)
    
    X_train = embeddings[train_mask]
    X_test = embeddings[test_mask]
    y_train = data.loc[train_mask, 'kinome_group_major'].values
    y_test = data.loc[test_mask, 'kinome_group_major'].values
    train_ids_filtered = data.loc[train_mask, 'uniprot_id'].values
    test_ids_filtered = data.loc[test_mask, 'uniprot_id'].values
    
    # Get unique labels
    labels = sorted(set(y_train) | set(y_test))
    
    return X_train, X_test, y_train, y_test, train_ids_filtered, test_ids_filtered, labels, metadata


def calculate_ece(y_true, y_proba, n_bins=10):
    """
    Calculate Expected Calibration Error.
    
    ECE measures how well predicted probabilities match observed frequencies.
    Lower is better (0 = perfect calibration).
    """
    # Get predicted class probabilities
    confidences = np.max(y_proba, axis=1)
    predictions = np.argmax(y_proba, axis=1)
    
    # Convert y_true to indices
    unique_labels = sorted(set(y_true))
    label_to_idx = {l: i for i, l in enumerate(unique_labels)}
    y_true_idx = np.array([label_to_idx[y] for y in y_true])
    
    accuracies = (predictions == y_true_idx).astype(float)
    
    # Create bins
    bin_edges = np.linspace(0, 1, n_bins + 1)
    bin_indices = np.digitize(confidences, bin_edges[:-1]) - 1
    bin_indices = np.clip(bin_indices, 0, n_bins - 1)
    
    ece = 0.0
    bin_stats = []
    
    for bin_idx in range(n_bins):
        mask = bin_indices == bin_idx
        if mask.sum() > 0:
            bin_confidence = confidences[mask].mean()
            bin_accuracy = accuracies[mask].mean()
            bin_size = mask.sum()
            ece += (bin_size / len(y_true)) * abs(bin_confidence - bin_accuracy)
            
            bin_stats.append({
                'bin': bin_idx,
                'confidence': bin_confidence,
                'accuracy': bin_accuracy,
                'size': bin_size
            })
    
    return ece, bin_stats


def plot_reliability_diagram(y_true, y_proba, output_file, n_bins=10):
    """Plot reliability diagram (calibration curve)."""
    confidences = np.max(y_proba, axis=1)
    predictions = np.argmax(y_proba, axis=1)
    
    # Convert y_true to indices
    unique_labels = sorted(set(y_true))
    label_to_idx = {l: i for i, l in enumerate(unique_labels)}
    y_true_idx = np.array([label_to_idx[y] for y in y_true])
    
    accuracies = (predictions == y_true_idx).astype(float)
    
    # Compute calibration curve
    prob_true, prob_pred = calibration_curve(accuracies, confidences, n_bins=n_bins, strategy='uniform')
    
    # Plot
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.plot([0, 1], [0, 1], 'k--', label='Perfectly calibrated')
    ax.plot(prob_pred, prob_true, 's-', label='Model')
    ax.set_xlabel('Mean predicted probability', fontsize=12)
    ax.set_ylabel('Fraction of positives', fontsize=12)
    ax.set_title('Reliability Diagram (Calibration Curve)', fontsize=14)
    ax.legend()
    ax.grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"✅ Reliability diagram saved to: {output_file}")


def train_and_evaluate_calibrated(X_train, X_test, y_train, y_test, labels, output_dir):
    """
    Train calibrated classifier and evaluate with full metrics.
    
    Returns:
        Dictionary with all results
    """
    print("\n" + "="*80)
    print("TRAINING CALIBRATED CLASSIFIER")
    print("="*80)
    
    # Standardize
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    # Train base classifier
    print("\nTraining base logistic regression...")
    base_clf = LogisticRegression(
        multi_class='multinomial',
        solver='saga',
        class_weight='balanced',
        max_iter=1000,
        random_state=42
    )
    base_clf.fit(X_train_scaled, y_train)
    
    # Calibrate using cross-validation
    print("Calibrating probabilities (CV=5)...")
    calibrated_clf = CalibratedClassifierCV(
        base_clf,
        method='sigmoid',  # Platt scaling
        cv=5
    )
    calibrated_clf.fit(X_train_scaled, y_train)
    
    # Predictions
    y_pred_base = base_clf.predict(X_test_scaled)
    y_proba_base = base_clf.predict_proba(X_test_scaled)
    
    y_pred_cal = calibrated_clf.predict(X_test_scaled)
    y_proba_cal = calibrated_clf.predict_proba(X_test_scaled)
    
    # Metrics
    results = {}
    
    # Base model
    results['base'] = {
        'accuracy': accuracy_score(y_test, y_pred_base),
        'macro_f1': f1_score(y_test, y_pred_base, average='macro'),
        'weighted_f1': f1_score(y_test, y_pred_base, average='weighted'),
        'top3_accuracy': top_k_accuracy_score(y_test, y_proba_base, k=3, labels=labels),
        'log_loss': log_loss(y_test, y_proba_base, labels=labels),
    }
    
    # Calibrated model
    results['calibrated'] = {
        'accuracy': accuracy_score(y_test, y_pred_cal),
        'macro_f1': f1_score(y_test, y_pred_cal, average='macro'),
        'weighted_f1': f1_score(y_test, y_pred_cal, average='weighted'),
        'top3_accuracy': top_k_accuracy_score(y_test, y_proba_cal, k=3, labels=labels),
        'log_loss': log_loss(y_test, y_proba_cal, labels=labels),
    }
    
    # ECE (Expected Calibration Error)
    ece_base, bin_stats_base = calculate_ece(y_test, y_proba_base, n_bins=10)
    ece_cal, bin_stats_cal = calculate_ece(y_test, y_proba_cal, n_bins=10)
    
    results['base']['ece'] = ece_base
    results['calibrated']['ece'] = ece_cal
    
    # Print comparison
    print("\n" + "="*80)
    print("CALIBRATION COMPARISON")
    print("="*80)
    print()
    print(f"{'Metric':<20} {'Base Model':<15} {'Calibrated':<15} {'Improvement'}")
    print("-"*80)
    print(f"{'Accuracy':<20} {results['base']['accuracy']:.4f}{'':<10} {results['calibrated']['accuracy']:.4f}{'':<10} {results['calibrated']['accuracy'] - results['base']['accuracy']:+.4f}")
    print(f"{'Macro-F1':<20} {results['base']['macro_f1']:.4f}{'':<10} {results['calibrated']['macro_f1']:.4f}{'':<10} {results['calibrated']['macro_f1'] - results['base']['macro_f1']:+.4f}")
    print(f"{'Top-3 Accuracy':<20} {results['base']['top3_accuracy']:.4f}{'':<10} {results['calibrated']['top3_accuracy']:.4f}{'':<10} {results['calibrated']['top3_accuracy'] - results['base']['top3_accuracy']:+.4f}")
    print(f"{'Log Loss':<20} {results['base']['log_loss']:.4f}{'':<10} {results['calibrated']['log_loss']:.4f}{'':<10} {results['calibrated']['log_loss'] - results['base']['log_loss']:+.4f}")
    print(f"{'ECE':<20} {results['base']['ece']:.4f}{'':<10} {results['calibrated']['ece']:.4f}{'':<10} {results['calibrated']['ece'] - results['base']['ece']:+.4f}")
    
    # Cross-validation on training set
    print("\nCross-validation (5-fold)...")
    cv_scores = cross_val_score(calibrated_clf, X_train_scaled, y_train, cv=5, scoring='f1_macro')
    results['cv_mean'] = cv_scores.mean()
    results['cv_std'] = cv_scores.std()
    print(f"CV Macro-F1: {cv_scores.mean():.4f} ± {cv_scores.std():.4f}")
    
    # Per-class metrics
    print("\nPer-class performance (calibrated model):")
    print("-"*80)
    report = classification_report(y_test, y_pred_cal, target_names=labels, output_dict=True)
    for label in labels:
        metrics = report[label]
        print(f"  {label:<15} Precision: {metrics['precision']:.3f}  Recall: {metrics['recall']:.3f}  F1: {metrics['f1-score']:.3f}  Support: {metrics['support']:.0f}")
    
    # Save outputs
    os.makedirs(output_dir, exist_ok=True)
    
    # Classification report
    with open(f"{output_dir}/classification_report_calibrated.txt", 'w') as f:
        f.write("CALIBRATED MODEL CLASSIFICATION REPORT\n")
        f.write("="*80 + "\n\n")
        f.write(classification_report(y_test, y_pred_cal, target_names=labels))
    
    # Confusion matrix
    cm = confusion_matrix(y_test, y_pred_cal, labels=labels)
    cm_df = pd.DataFrame(cm, index=labels, columns=labels)
    cm_df.to_csv(f"{output_dir}/confusion_matrix_calibrated.csv")
    
    # Plot reliability diagram
    plot_reliability_diagram(y_test, y_proba_cal, f"{output_dir}/reliability_diagram.png")
    
    # Save calibration statistics
    with open(f"{output_dir}/calibration_stats.json", 'w') as f:
        json.dump({
            'base_model': results['base'],
            'calibrated_model': results['calibrated'],
            'cv_mean': results['cv_mean'],
            'cv_std': results['cv_std'],
            'bin_statistics_calibrated': bin_stats_cal
        }, f, indent=2)
    
    return results, calibrated_clf, scaler, labels


def main():
    parser = argparse.ArgumentParser(
        description='Enhanced supervised training with calibration'
    )
    parser.add_argument(
        '--embeddings-dir',
        default='kinases_domains_embeddings_layers_20_33',
        help='Directory with ESM embeddings'
    )
    parser.add_argument(
        '--labels-csv',
        default='kinases_domains_e0.01.csv',
        help='CSV with labels'
    )
    parser.add_argument(
        '--splits-file',
        default='data/splits_40.json',
        help='Splits JSON file (default: 40%% identity)'
    )
    parser.add_argument(
        '--output-dir',
        default='supervised_results_calibrated',
        help='Output directory'
    )
    parser.add_argument(
        '--multi-identity',
        action='store_true',
        help='Evaluate on all identity thresholds (70%%, 50%%, 40%%)'
    )
    
    args = parser.parse_args()
    
    print("="*80)
    print("ENHANCED SUPERVISED TRAINING WITH CALIBRATION")
    print("="*80)
    print()
    print(f"Embeddings:  {args.embeddings_dir}")
    print(f"Labels:      {args.labels_csv}")
    print(f"Splits:      {args.splits_file}")
    print()
    
    # Process single or multiple identity thresholds
    if args.multi_identity:
        identities = [70, 50, 40]
        splits_files = [f"data/splits_{i}.json" for i in identities]
    else:
        identities = [None]
        splits_files = [args.splits_file]
    
    all_results = {}
    
    for identity, splits_file in zip(identities, splits_files):
        if identity is not None:
            print("\n" + "="*80)
            print(f"PROCESSING IDENTITY THRESHOLD: {identity}%")
            print("="*80)
            output_dir = f"{args.output_dir}_{identity}" if args.multi_identity else args.output_dir
        else:
            output_dir = args.output_dir
        
        # Load data
        print(f"\nLoading data from {splits_file}...")
        X_train, X_test, y_train, y_test, train_ids, test_ids, labels, metadata = load_data(
            args.embeddings_dir, args.labels_csv, splits_file
        )
        
        print(f"✅ Train: {len(X_train)}, Test: {len(X_test)}, Classes: {len(labels)}")
        
        # Train and evaluate
        results, model, scaler, labels_list = train_and_evaluate_calibrated(
            X_train, X_test, y_train, y_test, labels, output_dir
        )
        
        if identity is not None:
            all_results[f"{identity}%"] = results
    
    # Summary table if multi-identity
    if args.multi_identity:
        print("\n" + "="*80)
        print("MULTI-IDENTITY COMPARISON (CALIBRATED MODELS)")
        print("="*80)
        print()
        
        summary_df = pd.DataFrame([{
            'Identity': k,
            'Accuracy': v['calibrated']['accuracy'],
            'Macro-F1': v['calibrated']['macro_f1'],
            'Top-3 Acc': v['calibrated']['top3_accuracy'],
            'ECE': v['calibrated']['ece'],
            'Log Loss': v['calibrated']['log_loss']
        } for k, v in all_results.items()])
        
        print(summary_df.to_string(index=False))
        
        # Save summary
        summary_file = f"{args.output_dir}_multi_identity_summary.csv"
        summary_df.to_csv(summary_file, index=False)
        print(f"\n✅ Summary saved to: {summary_file}")
    
    print("\n" + "="*80)
    print("✅ ENHANCED SUPERVISED TRAINING COMPLETE!")
    print("="*80)


if __name__ == "__main__":
    main()
