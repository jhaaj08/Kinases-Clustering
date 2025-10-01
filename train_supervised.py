#!/usr/bin/env python3
"""
Supervised training: Multinomial logistic regression on ESM-2 embeddings.

This script demonstrates the "upper bound" of classification performance
using the same embeddings that produced clustering results.

Key comparisons:
- Unsupervised (clustering): ARI, NMI, purity
- Supervised (this script): Macro-F1, per-class F1, accuracy
"""

import os
import argparse
import numpy as np
import pandas as pd
import joblib
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_score
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.metrics import (
    classification_report,
    confusion_matrix,
    accuracy_score,
    f1_score,
)


def load_data(emb_dir, labels_file, exclude_other=True):
    """
    Load embeddings and align with labels.
    
    Returns:
        X: embeddings (N, 1280)
        y: labels (N,)
        ids: uniprot_ids (N,)
        label_names: original label strings
    """
    print("="*80)
    print("LOADING DATA")
    print("="*80)
    print()
    
    # Load embeddings
    emb_file = os.path.join(emb_dir, "esm2_embeddings.npy")
    idx_file = os.path.join(emb_dir, "esm2_index.csv")
    
    X = np.load(emb_file)
    ids_df = pd.read_csv(idx_file)
    ids = ids_df['uniprot_id'].astype(str).values
    
    print(f"✅ Loaded embeddings: {X.shape}")
    
    # Load labels
    labels_df = pd.read_csv(labels_file)
    labels_df['uniprot_id'] = labels_df['uniprot_id'].astype(str)
    
    # Merge
    emb_df = pd.DataFrame({'uniprot_id': ids, 'emb_idx': range(len(ids))})
    merged = emb_df.merge(
        labels_df[['uniprot_id', 'kinome_group_major']],
        on='uniprot_id',
        how='left'
    )
    
    print(f"✅ Loaded labels from: {labels_file}")
    
    # Filter
    if exclude_other:
        mask = merged['kinome_group_major'] != 'Other'
        merged = merged[mask]
        print(f"✅ Excluded 'Other' category")
    
    # Remove classes with too few samples (< 5) for stratified split
    class_counts = merged['kinome_group_major'].value_counts()
    small_classes = class_counts[class_counts < 5].index.tolist()
    
    if small_classes:
        print(f"⚠️  Removing classes with < 5 samples: {', '.join(small_classes)}")
        mask = ~merged['kinome_group_major'].isin(small_classes)
        merged = merged[mask]
    
    # Extract aligned data
    emb_indices = merged['emb_idx'].values
    X_aligned = X[emb_indices]
    label_names = merged['kinome_group_major'].values
    ids_aligned = merged['uniprot_id'].values
    
    # Encode labels
    le = LabelEncoder()
    y = le.fit_transform(label_names)
    
    print()
    print(f"Final dataset:")
    print(f"  Samples: {len(X_aligned):,}")
    print(f"  Features: {X_aligned.shape[1]}")
    print(f"  Classes: {len(le.classes_)}")
    print()
    
    print("Class distribution:")
    for i, cls in enumerate(le.classes_):
        count = (y == i).sum()
        pct = count / len(y) * 100
        print(f"  {cls:12s}: {count:4d} ({pct:5.1f}%)")
    print()
    
    return X_aligned, y, ids_aligned, label_names, le


def train_and_evaluate(X, y, label_encoder, random_state=42):
    """
    Train logistic regression with cross-validation and test evaluation.
    
    Returns:
        Dictionary with trained model and metrics
    """
    print("="*80)
    print("TRAINING & EVALUATION")
    print("="*80)
    print()
    
    # Stratified split
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, stratify=y, random_state=random_state
    )
    
    print(f"Train set: {len(X_train):,} samples")
    print(f"Test set:  {len(X_test):,} samples")
    print()
    
    # Standardize
    print("Standardizing features...")
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    print("✅ Standardization complete")
    print()
    
    # Train with class weighting
    print("Training Multinomial Logistic Regression...")
    print("  Solver: saga (for multinomial)")
    print("  Penalty: L2")
    print("  Class weight: balanced (handles imbalance)")
    print("  Max iter: 1000")
    print()
    
    clf = LogisticRegression(
        multi_class='multinomial',
        solver='saga',
        penalty='l2',
        C=1.0,
        class_weight='balanced',
        max_iter=1000,
        random_state=random_state,
        n_jobs=-1
    )
    
    clf.fit(X_train_scaled, y_train)
    print("✅ Training complete")
    print()
    
    # Cross-validation on train set
    print("Running 5-fold stratified cross-validation on train set...")
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=random_state)
    cv_scores = cross_val_score(clf, X_train_scaled, y_train, cv=cv, 
                                scoring='f1_macro', n_jobs=-1)
    
    print(f"  CV Macro-F1: {cv_scores.mean():.4f} ± {cv_scores.std():.4f}")
    print(f"  Individual folds: {[f'{s:.4f}' for s in cv_scores]}")
    print()
    
    # Test set evaluation
    print("Evaluating on held-out test set...")
    y_pred = clf.predict(X_test_scaled)
    y_proba = clf.predict_proba(X_test_scaled)
    
    test_acc = accuracy_score(y_test, y_pred)
    test_f1_macro = f1_score(y_test, y_pred, average='macro')
    test_f1_weighted = f1_score(y_test, y_pred, average='weighted')
    
    print(f"  Accuracy:       {test_acc:.4f}")
    print(f"  Macro-F1:       {test_f1_macro:.4f}")
    print(f"  Weighted-F1:    {test_f1_weighted:.4f}")
    print()
    
    # Per-class metrics
    print("Per-class performance:")
    print("-"*80)
    report = classification_report(
        y_test, y_pred,
        target_names=label_encoder.classes_,
        output_dict=True
    )
    
    for cls in label_encoder.classes_:
        metrics = report[cls]
        print(f"  {cls:12s}  Precision: {metrics['precision']:.3f}  "
              f"Recall: {metrics['recall']:.3f}  F1: {metrics['f1-score']:.3f}  "
              f"Support: {int(metrics['support'])}")
    print()
    
    # Confusion matrix
    cm = confusion_matrix(y_test, y_pred)
    
    return {
        'model': clf,
        'scaler': scaler,
        'label_encoder': label_encoder,
        'X_train': X_train_scaled,
        'X_test': X_test_scaled,
        'y_train': y_train,
        'y_test': y_test,
        'y_pred': y_pred,
        'y_proba': y_proba,
        'test_accuracy': test_acc,
        'test_f1_macro': test_f1_macro,
        'test_f1_weighted': test_f1_weighted,
        'cv_scores': cv_scores,
        'confusion_matrix': cm,
        'classification_report': report,
    }


def save_results(results, output_dir="supervised_results"):
    """Save trained model and results."""
    
    os.makedirs(output_dir, exist_ok=True)
    
    print("="*80)
    print("SAVING RESULTS")
    print("="*80)
    print()
    
    # Save model
    model_file = os.path.join(output_dir, "logistic_regression_model.joblib")
    joblib.dump({
        'model': results['model'],
        'scaler': results['scaler'],
        'label_encoder': results['label_encoder']
    }, model_file)
    print(f"✅ Saved model: {model_file}")
    
    # Save confusion matrix
    cm_df = pd.DataFrame(
        results['confusion_matrix'],
        index=results['label_encoder'].classes_,
        columns=results['label_encoder'].classes_
    )
    cm_file = os.path.join(output_dir, "confusion_matrix.csv")
    cm_df.to_csv(cm_file)
    print(f"✅ Saved confusion matrix: {cm_file}")
    
    # Save classification report
    report_file = os.path.join(output_dir, "classification_report.txt")
    with open(report_file, 'w') as f:
        f.write("SUPERVISED CLASSIFICATION RESULTS\n")
        f.write("="*80 + "\n\n")
        
        f.write("Model: Multinomial Logistic Regression (saga, L2, balanced)\n")
        f.write(f"Train size: {len(results['X_train']):,} samples\n")
        f.write(f"Test size:  {len(results['X_test']):,} samples\n")
        f.write("\n")
        
        f.write("Cross-Validation Results (5-fold on train):\n")
        f.write(f"  Macro-F1: {results['cv_scores'].mean():.4f} ± {results['cv_scores'].std():.4f}\n")
        f.write(f"  Folds: {', '.join([f'{s:.4f}' for s in results['cv_scores']])}\n")
        f.write("\n")
        
        f.write("Test Set Results:\n")
        f.write(f"  Accuracy:    {results['test_accuracy']:.4f}\n")
        f.write(f"  Macro-F1:    {results['test_f1_macro']:.4f}\n")
        f.write(f"  Weighted-F1: {results['test_f1_weighted']:.4f}\n")
        f.write("\n")
        
        f.write("Per-Class Performance:\n")
        f.write("-"*80 + "\n")
        f.write(f"{'Class':<12} {'Precision':>10} {'Recall':>10} {'F1-Score':>10} {'Support':>10}\n")
        f.write("-"*80 + "\n")
        
        report = results['classification_report']
        for cls in results['label_encoder'].classes_:
            metrics = report[cls]
            f.write(f"{cls:<12} {metrics['precision']:>10.3f} {metrics['recall']:>10.3f} "
                   f"{metrics['f1-score']:>10.3f} {int(metrics['support']):>10d}\n")
        
        f.write("\n")
        f.write(f"Macro avg    {report['macro avg']['precision']:>10.3f} "
               f"{report['macro avg']['recall']:>10.3f} "
               f"{report['macro avg']['f1-score']:>10.3f} "
               f"{int(report['macro avg']['support']):>10d}\n")
        f.write(f"Weighted avg {report['weighted avg']['precision']:>10.3f} "
               f"{report['weighted avg']['recall']:>10.3f} "
               f"{report['weighted avg']['f1-score']:>10.3f} "
               f"{int(report['weighted avg']['support']):>10d}\n")
    
    print(f"✅ Saved report: {report_file}")
    print()


def compare_with_clustering(supervised_results, output_dir="supervised_results"):
    """Generate comparison with clustering results."""
    
    print("="*80)
    print("SUPERVISED vs UNSUPERVISED COMPARISON")
    print("="*80)
    print()
    
    # Load clustering results if available
    clustering_file = "clustering/systematic_experiments_results.csv"
    
    comparison_text = []
    comparison_text.append("SUPERVISED vs UNSUPERVISED COMPARISON\n")
    comparison_text.append("="*80 + "\n\n")
    
    if os.path.exists(clustering_file):
        clustering_df = pd.read_csv(clustering_file)
        best_clustering = clustering_df.iloc[0]  # Best result (sorted by ARI)
        
        comparison_text.append(f"Best Unsupervised (Clustering):\n")
        comparison_text.append(f"  Configuration: {best_clustering['experiment']}\n")
        comparison_text.append(f"  ARI:           {best_clustering['ARI']:.4f}\n")
        comparison_text.append(f"  NMI:           {best_clustering['NMI']:.4f}\n")
        comparison_text.append(f"  Hungarian Acc: {best_clustering['Hungarian']:.4f}\n")
        comparison_text.append(f"  Purity:        {best_clustering['Purity']:.4f}\n")
        comparison_text.append("\n")
        
        print("Best Unsupervised (K-Means Clustering):")
        print(f"  ARI:           {best_clustering['ARI']:.4f}")
        print(f"  NMI:           {best_clustering['NMI']:.4f}")
        print(f"  Hungarian Acc: {best_clustering['Hungarian']:.4f}")
        print()
    
    comparison_text.append(f"Supervised (Logistic Regression):\n")
    comparison_text.append(f"  Test Accuracy: {supervised_results['test_accuracy']:.4f}\n")
    comparison_text.append(f"  Macro-F1:      {supervised_results['test_f1_macro']:.4f}\n")
    comparison_text.append(f"  Weighted-F1:   {supervised_results['test_f1_weighted']:.4f}\n")
    comparison_text.append(f"  CV Macro-F1:   {supervised_results['cv_scores'].mean():.4f} ± {supervised_results['cv_scores'].std():.4f}\n")
    comparison_text.append("\n")
    
    print("Supervised (Logistic Regression on same embeddings):")
    print(f"  Test Accuracy: {supervised_results['test_accuracy']:.4f}")
    print(f"  Macro-F1:      {supervised_results['test_f1_macro']:.4f}")
    print(f"  Weighted-F1:   {supervised_results['test_f1_weighted']:.4f}")
    print()
    
    comparison_text.append("Key Insights:\n")
    comparison_text.append("-"*80 + "\n")
    comparison_text.append("1. Supervised learning uses labels directly, so it achieves higher\n")
    comparison_text.append("   accuracy than unsupervised clustering's Hungarian matching.\n")
    comparison_text.append("\n")
    comparison_text.append("2. However, clustering validated the embedding quality WITHOUT labels,\n")
    comparison_text.append("   which guided feature engineering (domain extraction, layer selection).\n")
    comparison_text.append("\n")
    comparison_text.append("3. Both use the same embeddings (domain ESM-2, layers 20-33 averaged),\n")
    comparison_text.append("   demonstrating that better unsupervised separability predicts better\n")
    comparison_text.append("   supervised performance.\n")
    comparison_text.append("\n")
    
    print("Key Finding:")
    print("  Supervised accuracy is higher (as expected with labels), but")
    print("  clustering guided us to the right embedding space first!")
    print()
    
    # Save comparison
    comparison_file = os.path.join(output_dir, "supervised_vs_clustering.txt")
    with open(comparison_file, 'w') as f:
        f.writelines(comparison_text)
    
    print(f"✅ Saved comparison: {comparison_file}")
    print()


def main():
    parser = argparse.ArgumentParser(
        description='Train supervised classifier on ESM-2 embeddings'
    )
    parser.add_argument(
        '--emb-dir',
        default='kinases_domains_e0.01_layers_mid',
        help='Directory containing embeddings'
    )
    parser.add_argument(
        '--labels',
        default='kinases_domains_e0.01.csv',
        help='CSV file with labels'
    )
    parser.add_argument(
        '--output-dir',
        default='supervised_results',
        help='Output directory for results'
    )
    parser.add_argument(
        '--seed',
        type=int,
        default=42,
        help='Random seed for reproducibility'
    )
    
    args = parser.parse_args()
    
    print()
    print("╔" + "="*78 + "╗")
    print("║" + " "*20 + "SUPERVISED LEARNING PIPELINE" + " "*29 + "║")
    print("╚" + "="*78 + "╝")
    print()
    
    print(f"Embeddings:  {args.emb_dir}")
    print(f"Labels:      {args.labels}")
    print(f"Output:      {args.output_dir}")
    print(f"Random seed: {args.seed}")
    print()
    
    # Load data
    X, y, ids, label_names, label_encoder = load_data(
        args.emb_dir, args.labels, exclude_other=True
    )
    
    # Train and evaluate
    results = train_and_evaluate(X, y, label_encoder, random_state=args.seed)
    
    # Save results
    save_results(results, output_dir=args.output_dir)
    
    # Compare with clustering
    compare_with_clustering(results, output_dir=args.output_dir)
    
    print("="*80)
    print("✅ SUPERVISED TRAINING COMPLETE!")
    print("="*80)
    print()
    print(f"Key Results:")
    print(f"  Test Accuracy: {results['test_accuracy']:.1%}")
    print(f"  Macro-F1:      {results['test_f1_macro']:.4f}")
    print(f"  CV Macro-F1:   {results['cv_scores'].mean():.4f} ± {results['cv_scores'].std():.4f}")
    print()
    print(f"All results saved to: {args.output_dir}/")
    print()


if __name__ == "__main__":
    main()

