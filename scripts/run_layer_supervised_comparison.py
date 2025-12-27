#!/usr/bin/env python3
"""
Run supervised classification experiments across different layer configurations.

Compares:
1. Final layer only (layer 33)
2. Layers 20-30 (mid)
3. Layers 19-33 (extended mid)

Outputs results for Figure 4: Supervised Classification Performance.
"""

import os
import json
import numpy as np
import pandas as pd
from pathlib import Path

from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (
    accuracy_score, f1_score, confusion_matrix,
    classification_report, precision_recall_fscore_support
)
from sklearn.model_selection import cross_val_score

import warnings
warnings.filterwarnings('ignore')


def load_data(embeddings_dir, labels_csv, splits_file):
    """Load embeddings, labels, and splits."""
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
    with open(splits_file, 'r') as f:
        splits = json.load(f)
    train_ids = splits['train_ids']
    test_ids = splits['test_ids']
    
    # Filter by splits
    train_mask = data['uniprot_id'].isin(train_ids)
    test_mask = data['uniprot_id'].isin(test_ids)
    
    X_train = embeddings[train_mask]
    X_test = embeddings[test_mask]
    y_train = data.loc[train_mask, 'kinome_group_major'].values
    y_test = data.loc[test_mask, 'kinome_group_major'].values
    
    # Get unique labels (excluding Other and NaN)
    valid_mask_train = pd.notna(y_train) & (y_train != 'Other')
    valid_mask_test = pd.notna(y_test) & (y_test != 'Other')
    
    X_train = X_train[valid_mask_train]
    y_train = y_train[valid_mask_train]
    X_test = X_test[valid_mask_test]
    y_test = y_test[valid_mask_test]
    
    labels = sorted(set(y_train) | set(y_test))
    
    return X_train, X_test, y_train, y_test, labels


def train_and_evaluate(X_train, X_test, y_train, y_test, labels):
    """Train logistic regression and return metrics."""
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
    
    # Metrics
    accuracy = accuracy_score(y_test, y_pred)
    macro_f1 = f1_score(y_test, y_pred, average='macro')
    weighted_f1 = f1_score(y_test, y_pred, average='weighted')
    
    # Cross-validation
    cv_scores = cross_val_score(clf, X_train_scaled, y_train, cv=5, scoring='f1_macro')
    
    # Confusion matrix
    cm = confusion_matrix(y_test, y_pred, labels=labels)
    
    # Per-class metrics
    precision, recall, f1, support = precision_recall_fscore_support(
        y_test, y_pred, labels=labels, average=None
    )
    
    per_class = {}
    for i, label in enumerate(labels):
        per_class[label] = {
            'precision': precision[i],
            'recall': recall[i],
            'f1': f1[i],
            'support': int(support[i])
        }
    
    return {
        'accuracy': accuracy,
        'macro_f1': macro_f1,
        'weighted_f1': weighted_f1,
        'cv_mean': cv_scores.mean(),
        'cv_std': cv_scores.std(),
        'confusion_matrix': cm,
        'per_class': per_class,
        'labels': labels,
        'train_size': len(X_train),
        'test_size': len(X_test)
    }


def main():
    print("="*80)
    print("LAYER CONFIGURATION SUPERVISED CLASSIFICATION COMPARISON")
    print("="*80)
    print()
    
    # Configuration
    layer_configs = {
        'Final (Layer 33)': 'kinases_domains_e0.01_embeddings',
        'Mid (Layers 20-30)': 'kinases_domains_e0.01_layers_20_30',
        'Extended Mid (Layers 19-33)': 'kinases_domains_e0.01_layers_mid',
    }
    
    labels_csv = 'data/processed/kinases_domains.csv'
    splits_file = 'data/splits_40.json'
    output_dir = 'supervised_results_layer_comparison'
    
    # Check files exist
    if not Path(labels_csv).exists():
        labels_csv = 'data/processed/kinases_domains_e0.01.csv'
    
    os.makedirs(output_dir, exist_ok=True)
    
    all_results = {}
    all_confusion_matrices = {}
    
    for config_name, emb_dir in layer_configs.items():
        print(f"\n{'='*60}")
        print(f"Configuration: {config_name}")
        print(f"Embeddings: {emb_dir}")
        print(f"{'='*60}")
        
        if not Path(emb_dir).exists():
            print(f"  ⚠️ Directory not found: {emb_dir}")
            continue
        
        # Load and train
        X_train, X_test, y_train, y_test, labels = load_data(
            emb_dir, labels_csv, splits_file
        )
        print(f"  Train: {len(X_train)}, Test: {len(X_test)}, Classes: {len(labels)}")
        
        results = train_and_evaluate(X_train, X_test, y_train, y_test, labels)
        
        print(f"\n  Results:")
        print(f"    Accuracy:    {results['accuracy']:.4f}")
        print(f"    Macro-F1:    {results['macro_f1']:.4f}")
        print(f"    Weighted-F1: {results['weighted_f1']:.4f}")
        print(f"    CV Macro-F1: {results['cv_mean']:.4f} ± {results['cv_std']:.4f}")
        
        all_results[config_name] = results
        all_confusion_matrices[config_name] = results['confusion_matrix']
    
    # Summary table
    print("\n" + "="*80)
    print("SUMMARY COMPARISON")
    print("="*80)
    print()
    
    summary_data = []
    for config_name, results in all_results.items():
        summary_data.append({
            'Configuration': config_name,
            'Accuracy': results['accuracy'],
            'Macro-F1': results['macro_f1'],
            'Weighted-F1': results['weighted_f1'],
            'CV Macro-F1': f"{results['cv_mean']:.4f} ± {results['cv_std']:.4f}",
            'Train Size': results['train_size'],
            'Test Size': results['test_size']
        })
    
    summary_df = pd.DataFrame(summary_data)
    print(summary_df.to_string(index=False))
    
    # Calculate improvements
    if 'Final (Layer 33)' in all_results and 'Extended Mid (Layers 19-33)' in all_results:
        baseline = all_results['Final (Layer 33)']
        best = all_results['Extended Mid (Layers 19-33)']
        
        acc_improvement = (best['accuracy'] - baseline['accuracy']) / baseline['accuracy'] * 100
        f1_improvement = (best['macro_f1'] - baseline['macro_f1']) / baseline['macro_f1'] * 100
        
        print(f"\n  Improvement (Extended Mid vs Final):")
        print(f"    Accuracy:  +{acc_improvement:.1f}%")
        print(f"    Macro-F1:  +{f1_improvement:.1f}%")
    
    # Save results
    summary_df.to_csv(f"{output_dir}/layer_comparison_summary.csv", index=False)
    
    # Save per-class metrics
    per_class_data = []
    for config_name, results in all_results.items():
        for label, metrics in results['per_class'].items():
            per_class_data.append({
                'Configuration': config_name,
                'Class': label,
                'Precision': metrics['precision'],
                'Recall': metrics['recall'],
                'F1': metrics['f1'],
                'Support': metrics['support']
            })
    
    per_class_df = pd.DataFrame(per_class_data)
    per_class_df.to_csv(f"{output_dir}/per_class_metrics.csv", index=False)
    
    # Save confusion matrices
    for config_name, cm in all_confusion_matrices.items():
        labels = all_results[config_name]['labels']
        cm_df = pd.DataFrame(cm, index=labels, columns=labels)
        safe_name = config_name.replace(' ', '_').replace('(', '').replace(')', '').replace('-', '_')
        cm_df.to_csv(f"{output_dir}/confusion_matrix_{safe_name}.csv")
    
    # Save full results as JSON
    json_results = {}
    for config_name, results in all_results.items():
        json_results[config_name] = {
            'accuracy': results['accuracy'],
            'macro_f1': results['macro_f1'],
            'weighted_f1': results['weighted_f1'],
            'cv_mean': results['cv_mean'],
            'cv_std': results['cv_std'],
            'per_class': results['per_class'],
            'train_size': results['train_size'],
            'test_size': results['test_size']
        }
    
    with open(f"{output_dir}/full_results.json", 'w') as f:
        json.dump(json_results, f, indent=2)
    
    print(f"\n✅ Results saved to: {output_dir}/")
    print(f"   - layer_comparison_summary.csv")
    print(f"   - per_class_metrics.csv")
    print(f"   - confusion_matrix_*.csv")
    print(f"   - full_results.json")
    
    print("\n" + "="*80)
    print("✅ LAYER COMPARISON COMPLETE!")
    print("="*80)


if __name__ == "__main__":
    main()



