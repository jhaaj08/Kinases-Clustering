#!/usr/bin/env python3
"""
Compute calibration curves for Final vs Mid layer configurations.

Generates data for Figure 6: Calibration and Retrieval Quality.
"""

import os
import json
import numpy as np
import pandas as pd
from pathlib import Path

from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.calibration import calibration_curve

import warnings
warnings.filterwarnings('ignore')


def load_data(embeddings_dir, labels_csv, splits_file):
    """Load embeddings, labels, and splits."""
    embeddings = np.load(f"{embeddings_dir}/esm2_embeddings.npy")
    index_df = pd.read_csv(f"{embeddings_dir}/esm2_index.csv")
    labels_df = pd.read_csv(labels_csv)
    
    data = index_df.merge(
        labels_df[['uniprot_id', 'kinome_group_major']], 
        on='uniprot_id', 
        how='left'
    )
    
    with open(splits_file, 'r') as f:
        splits = json.load(f)
    
    train_mask = data['uniprot_id'].isin(splits['train_ids'])
    test_mask = data['uniprot_id'].isin(splits['test_ids'])
    
    X_train = embeddings[train_mask]
    X_test = embeddings[test_mask]
    y_train = data.loc[train_mask, 'kinome_group_major'].values
    y_test = data.loc[test_mask, 'kinome_group_major'].values
    
    valid_train = pd.notna(y_train) & (y_train != 'Other')
    valid_test = pd.notna(y_test) & (y_test != 'Other')
    
    return (X_train[valid_train], X_test[valid_test], 
            y_train[valid_train], y_test[valid_test])


def calculate_ece(y_true, y_proba, n_bins=10):
    """Calculate Expected Calibration Error."""
    confidences = np.max(y_proba, axis=1)
    predictions = np.argmax(y_proba, axis=1)
    
    unique_labels = sorted(set(y_true))
    label_to_idx = {l: i for i, l in enumerate(unique_labels)}
    y_true_idx = np.array([label_to_idx[y] for y in y_true])
    
    accuracies = (predictions == y_true_idx).astype(float)
    
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
                'confidence': float(bin_confidence),
                'accuracy': float(bin_accuracy),
                'size': int(bin_size)
            })
    
    return ece, bin_stats


def get_calibration_curve(y_true, y_proba, n_bins=10):
    """Compute calibration curve data."""
    confidences = np.max(y_proba, axis=1)
    predictions = np.argmax(y_proba, axis=1)
    
    unique_labels = sorted(set(y_true))
    label_to_idx = {l: i for i, l in enumerate(unique_labels)}
    y_true_idx = np.array([label_to_idx[y] for y in y_true])
    
    correct = (predictions == y_true_idx).astype(float)
    
    prob_true, prob_pred = calibration_curve(correct, confidences, n_bins=n_bins, strategy='uniform')
    
    return prob_pred, prob_true


def train_and_calibrate(X_train, X_test, y_train, y_test):
    """Train model and return calibration data."""
    scaler = StandardScaler()
    X_train_s = scaler.fit_transform(X_train)
    X_test_s = scaler.transform(X_test)
    
    clf = LogisticRegression(
        multi_class='multinomial',
        solver='saga',
        class_weight='balanced',
        max_iter=1000,
        random_state=42
    )
    clf.fit(X_train_s, y_train)
    
    y_proba = clf.predict_proba(X_test_s)
    
    ece, bin_stats = calculate_ece(y_test, y_proba)
    prob_pred, prob_true = get_calibration_curve(y_test, y_proba)
    
    return {
        'ece': ece,
        'bin_stats': bin_stats,
        'calibration_curve': {
            'prob_pred': prob_pred.tolist(),
            'prob_true': prob_true.tolist()
        }
    }


def main():
    print("="*80)
    print("CALIBRATION COMPARISON: FINAL VS MID LAYERS")
    print("="*80)
    print()
    
    layer_configs = {
        'Final (Layer 33)': 'kinases_domains_e0.01_embeddings',
        'Mid (Layers 19-33)': 'kinases_domains_e0.01_layers_mid',
    }
    
    labels_csv = 'data/processed/kinases_domains.csv'
    if not Path(labels_csv).exists():
        labels_csv = 'data/processed/kinases_domains_e0.01.csv'
    
    splits_file = 'data/splits_40.json'
    output_dir = 'calibration_comparison_results'
    os.makedirs(output_dir, exist_ok=True)
    
    all_results = {}
    
    for config_name, emb_dir in layer_configs.items():
        print(f"\n{config_name}:")
        
        if not Path(emb_dir).exists():
            print(f"  ⚠️ Not found: {emb_dir}")
            continue
        
        X_train, X_test, y_train, y_test = load_data(emb_dir, labels_csv, splits_file)
        print(f"  Train: {len(X_train)}, Test: {len(X_test)}")
        
        results = train_and_calibrate(X_train, X_test, y_train, y_test)
        print(f"  ECE: {results['ece']:.4f}")
        
        all_results[config_name] = results
    
    # Save results
    with open(f"{output_dir}/calibration_comparison.json", 'w') as f:
        json.dump(all_results, f, indent=2)
    
    # Create summary DataFrame
    summary_data = []
    for config_name, results in all_results.items():
        summary_data.append({
            'Configuration': config_name,
            'ECE': results['ece']
        })
    
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(f"{output_dir}/ece_comparison.csv", index=False)
    
    print(f"\n✅ Results saved to: {output_dir}/")
    print(f"   - calibration_comparison.json")
    print(f"   - ece_comparison.csv")
    
    print("\n" + "="*80)
    print("✅ CALIBRATION COMPARISON COMPLETE!")
    print("="*80)


if __name__ == "__main__":
    main()

