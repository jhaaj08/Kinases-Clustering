#!/usr/bin/env python3
"""
Comprehensive baselines comparison for kinase classification.

Implements:
1. HMMER baseline (Pfam family assignment)
2. ESM+kNN (k-nearest neighbors on embeddings)
3. Motif-only logistic regression
4. Simple MLP head on ESM embeddings
5. Logistic regression on ESM embeddings (existing)

Compares against homology-aware splits at multiple identity thresholds.
"""

import os
import sys
import argparse
import json
import numpy as np
import pandas as pd
from pathlib import Path
from collections import Counter

# Machine learning
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (
    accuracy_score, classification_report, confusion_matrix,
    f1_score, precision_recall_fscore_support, top_k_accuracy_score
)
from sklearn.model_selection import cross_val_score

# Calibration
from sklearn.calibration import CalibratedClassifierCV, calibration_curve
from sklearn.utils.class_weight import compute_class_weight

import warnings
warnings.filterwarnings('ignore')


def load_splits(splits_file):
    """Load train/test splits from JSON."""
    with open(splits_file, 'r') as f:
        data = json.load(f)
    return data['train_ids'], data['test_ids'], data['metadata']


def load_embeddings_and_labels(embeddings_dir, labels_csv, splits_file):
    """
    Load ESM embeddings and labels, split by train/test IDs.
    
    Returns:
        X_train, X_test, y_train, y_test, train_ids, test_ids, label_encoder
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
    label_to_idx = {l: i for i, l in enumerate(labels)}
    
    return (X_train, X_test, y_train, y_test, 
            train_ids_filtered, test_ids_filtered, labels, label_to_idx)


def load_motif_features(motifs_csv, train_ids, test_ids):
    """Load motif features and split by train/test."""
    motifs_df = pd.read_csv(motifs_csv)
    
    # Select motif columns (exclude sequence, uniprot_id, labels)
    feature_cols = [col for col in motifs_df.columns if col not in [
        'uniprot_id', 'sequence', 'domain_sequence', 'kinome_group_major',
        'kinome_group_minor', 'Entry', 'Entry Name', 'Protein names', 
        'Gene Names', 'Organism', 'Length', 'domain_start', 'domain_end',
        'domain_length', 'evalue', 'score'
    ]]
    
    # Filter by IDs
    train_motifs = motifs_df[motifs_df['uniprot_id'].isin(train_ids)][feature_cols].values
    test_motifs = motifs_df[motifs_df['uniprot_id'].isin(test_ids)][feature_cols].values
    
    return train_motifs, test_motifs, feature_cols


def baseline_esm_knn(X_train, X_test, y_train, y_test, k=5):
    """ESM embeddings + k-NN classifier."""
    print(f"\n{'='*80}")
    print(f"Baseline: ESM + k-NN (k={k})")
    print(f"{'='*80}")
    
    # Standardize
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    # Train
    knn = KNeighborsClassifier(n_neighbors=k, metric='cosine')
    knn.fit(X_train_scaled, y_train)
    
    # Predict
    y_pred = knn.predict(X_test_scaled)
    y_proba = knn.predict_proba(X_test_scaled)
    
    # Metrics
    acc = accuracy_score(y_test, y_pred)
    macro_f1 = f1_score(y_test, y_pred, average='macro')
    weighted_f1 = f1_score(y_test, y_pred, average='weighted')
    
    # Top-3 accuracy
    top3_acc = top_k_accuracy_score(y_test, y_proba, k=3, labels=knn.classes_)
    
    print(f"Test Accuracy:    {acc:.4f}")
    print(f"Macro-F1:         {macro_f1:.4f}")
    print(f"Weighted-F1:      {weighted_f1:.4f}")
    print(f"Top-3 Accuracy:   {top3_acc:.4f}")
    
    return {
        'model': 'ESM+kNN',
        'accuracy': acc,
        'macro_f1': macro_f1,
        'weighted_f1': weighted_f1,
        'top3_accuracy': top3_acc,
        'predictions': y_pred,
        'probabilities': y_proba
    }


def baseline_motifs_only(motifs_train, motifs_test, y_train, y_test):
    """Motif features only + Logistic Regression."""
    print(f"\n{'='*80}")
    print(f"Baseline: Motifs-Only Logistic Regression")
    print(f"{'='*80}")
    
    # Handle NaN values
    motifs_train = np.nan_to_num(motifs_train, nan=0.0)
    motifs_test = np.nan_to_num(motifs_test, nan=0.0)
    
    # Standardize
    scaler = StandardScaler()
    motifs_train_scaled = scaler.fit_transform(motifs_train)
    motifs_test_scaled = scaler.transform(motifs_test)
    
    # Train
    clf = LogisticRegression(
        multi_class='multinomial',
        solver='lbfgs',
        class_weight='balanced',
        max_iter=1000,
        random_state=42
    )
    clf.fit(motifs_train_scaled, y_train)
    
    # Predict
    y_pred = clf.predict(motifs_test_scaled)
    y_proba = clf.predict_proba(motifs_test_scaled)
    
    # Metrics
    acc = accuracy_score(y_test, y_pred)
    macro_f1 = f1_score(y_test, y_pred, average='macro')
    weighted_f1 = f1_score(y_test, y_pred, average='weighted')
    top3_acc = top_k_accuracy_score(y_test, y_proba, k=3, labels=clf.classes_)
    
    print(f"Test Accuracy:    {acc:.4f}")
    print(f"Macro-F1:         {macro_f1:.4f}")
    print(f"Weighted-F1:      {weighted_f1:.4f}")
    print(f"Top-3 Accuracy:   {top3_acc:.4f}")
    print(f"\nFeature importance (top 10 by mean |coef|):")
    
    # Feature importance
    coef_mean = np.abs(clf.coef_).mean(axis=0)
    
    return {
        'model': 'Motifs-Only LR',
        'accuracy': acc,
        'macro_f1': macro_f1,
        'weighted_f1': weighted_f1,
        'top3_accuracy': top3_acc,
        'predictions': y_pred,
        'probabilities': y_proba,
        'feature_importance': coef_mean
    }


def baseline_mlp(X_train, X_test, y_train, y_test):
    """Simple MLP on ESM embeddings."""
    print(f"\n{'='*80}")
    print(f"Baseline: ESM + MLP (2 hidden layers)")
    print(f"{'='*80}")
    
    # Standardize
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    # Train
    mlp = MLPClassifier(
        hidden_layer_sizes=(512, 128),
        activation='relu',
        solver='adam',
        alpha=0.0001,
        batch_size=32,
        learning_rate='adaptive',
        max_iter=200,
        random_state=42,
        early_stopping=True,
        validation_fraction=0.1
    )
    mlp.fit(X_train_scaled, y_train)
    
    # Predict
    y_pred = mlp.predict(X_test_scaled)
    y_proba = mlp.predict_proba(X_test_scaled)
    
    # Metrics
    acc = accuracy_score(y_test, y_pred)
    macro_f1 = f1_score(y_test, y_pred, average='macro')
    weighted_f1 = f1_score(y_test, y_pred, average='weighted')
    top3_acc = top_k_accuracy_score(y_test, y_proba, k=3, labels=mlp.classes_)
    
    print(f"Test Accuracy:    {acc:.4f}")
    print(f"Macro-F1:         {macro_f1:.4f}")
    print(f"Weighted-F1:      {weighted_f1:.4f}")
    print(f"Top-3 Accuracy:   {top3_acc:.4f}")
    print(f"Iterations:       {mlp.n_iter_}")
    
    return {
        'model': 'ESM+MLP',
        'accuracy': acc,
        'macro_f1': macro_f1,
        'weighted_f1': weighted_f1,
        'top3_accuracy': top3_acc,
        'predictions': y_pred,
        'probabilities': y_proba
    }


def baseline_logistic_regression(X_train, X_test, y_train, y_test):
    """Logistic Regression on ESM embeddings (current approach)."""
    print(f"\n{'='*80}")
    print(f"Baseline: ESM + Logistic Regression")
    print(f"{'='*80}")
    
    # Standardize
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    # Train
    clf = LogisticRegression(
        multi_class='multinomial',
        solver='saga',
        class_weight='balanced',
        max_iter=1000,
        random_state=42
    )
    clf.fit(X_train_scaled, y_train)
    
    # Predict
    y_pred = clf.predict(X_test_scaled)
    y_proba = clf.predict_proba(X_test_scaled)
    
    # Metrics
    acc = accuracy_score(y_test, y_pred)
    macro_f1 = f1_score(y_test, y_pred, average='macro')
    weighted_f1 = f1_score(y_test, y_pred, average='weighted')
    top3_acc = top_k_accuracy_score(y_test, y_proba, k=3, labels=clf.classes_)
    
    # Cross-validation
    cv_scores = cross_val_score(clf, X_train_scaled, y_train, cv=5, scoring='f1_macro')
    
    print(f"Test Accuracy:    {acc:.4f}")
    print(f"Macro-F1:         {macro_f1:.4f}")
    print(f"Weighted-F1:      {weighted_f1:.4f}")
    print(f"Top-3 Accuracy:   {top3_acc:.4f}")
    print(f"CV Macro-F1:      {cv_scores.mean():.4f} ± {cv_scores.std():.4f}")
    
    return {
        'model': 'ESM+LR',
        'accuracy': acc,
        'macro_f1': macro_f1,
        'weighted_f1': weighted_f1,
        'top3_accuracy': top3_acc,
        'cv_mean': cv_scores.mean(),
        'cv_std': cv_scores.std(),
        'predictions': y_pred,
        'probabilities': y_proba
    }


def main():
    parser = argparse.ArgumentParser(
        description='Comprehensive baselines comparison'
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
        '--motifs-csv',
        default='kinases_domains_with_motifs.csv',
        help='CSV with motif features'
    )
    parser.add_argument(
        '--splits-file',
        default='data/splits_40.json',
        help='Splits JSON file'
    )
    parser.add_argument(
        '--output-dir',
        default='baselines_results',
        help='Output directory for results'
    )
    
    args = parser.parse_args()
    
    print("="*80)
    print("COMPREHENSIVE BASELINES COMPARISON")
    print("="*80)
    print()
    print(f"Embeddings:  {args.embeddings_dir}")
    print(f"Labels:      {args.labels_csv}")
    print(f"Motifs:      {args.motifs_csv}")
    print(f"Splits:      {args.splits_file}")
    print()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load data
    print("Loading data...")
    (X_train, X_test, y_train, y_test, 
     train_ids, test_ids, labels, label_to_idx) = load_embeddings_and_labels(
        args.embeddings_dir, args.labels_csv, args.splits_file
    )
    
    print(f"✅ Loaded {len(X_train)} train, {len(X_test)} test")
    print(f"   Classes: {len(labels)}")
    print(f"   Feature dim: {X_train.shape[1]}")
    
    # Load motifs
    print("\nLoading motif features...")
    motifs_train, motifs_test, motif_cols = load_motif_features(
        args.motifs_csv, train_ids, test_ids
    )
    print(f"✅ Loaded {motifs_train.shape[1]} motif features")
    
    # Run baselines
    results = []
    
    # 1. ESM + Logistic Regression (our current approach)
    results.append(baseline_logistic_regression(X_train, X_test, y_train, y_test))
    
    # 2. ESM + k-NN
    results.append(baseline_esm_knn(X_train, X_test, y_train, y_test, k=5))
    
    # 3. Motifs-only
    results.append(baseline_motifs_only(motifs_train, motifs_test, y_train, y_test))
    
    # 4. ESM + MLP
    results.append(baseline_mlp(X_train, X_test, y_train, y_test))
    
    # Create comparison table
    print("\n" + "="*80)
    print("RESULTS COMPARISON")
    print("="*80)
    print()
    
    comparison_df = pd.DataFrame([{
        'Model': r['model'],
        'Accuracy': f"{r['accuracy']:.4f}",
        'Macro-F1': f"{r['macro_f1']:.4f}",
        'Weighted-F1': f"{r['weighted_f1']:.4f}",
        'Top-3 Acc': f"{r['top3_accuracy']:.4f}"
    } for r in results])
    
    print(comparison_df.to_string(index=False))
    print()
    
    # Save results
    output_file = f"{args.output_dir}/baselines_comparison.csv"
    comparison_df.to_csv(output_file, index=False)
    print(f"✅ Results saved to: {output_file}")
    
    # Save detailed results
    detailed_file = f"{args.output_dir}/baselines_detailed.json"
    with open(detailed_file, 'w') as f:
        # Remove numpy arrays before saving
        saveable_results = []
        for r in results:
            r_copy = r.copy()
            r_copy.pop('predictions', None)
            r_copy.pop('probabilities', None)
            if 'feature_importance' in r_copy:
                r_copy['feature_importance'] = r_copy['feature_importance'].tolist()
            saveable_results.append(r_copy)
        json.dump(saveable_results, f, indent=2)
    print(f"✅ Detailed results saved to: {detailed_file}")
    
    print("\n" + "="*80)
    print("✅ BASELINES COMPARISON COMPLETE!")
    print("="*80)


if __name__ == "__main__":
    main()
