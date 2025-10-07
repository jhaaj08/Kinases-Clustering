#!/usr/bin/env python3
"""
Comprehensive statistical analysis for unsupervised clustering.

Implements reviewer requirements:
1. Bootstrapped confidence intervals for all metrics
2. Permutation tests for comparisons (domain vs full-length)
3. Effect size calculations (Cohen's d)
4. Outlier detection (cluster-flipping sequences)
5. Statistical significance testing
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from collections import Counter

from sklearn.metrics import (
    adjusted_rand_score,
    normalized_mutual_info_score,
    homogeneity_completeness_v_measure,
    silhouette_score
)
from sklearn.utils import resample
from scipy import stats
from scipy.optimize import linear_sum_assignment
import warnings
warnings.filterwarnings('ignore')


def calculate_purity(labels_true, labels_pred):
    """Calculate cluster purity."""
    contingency = pd.crosstab(labels_pred, labels_true)
    purity = contingency.max(axis=1).sum() / len(labels_true)
    return purity


def calculate_hungarian_accuracy(labels_true, labels_pred):
    """Calculate Hungarian (optimal assignment) accuracy."""
    contingency = pd.crosstab(labels_pred, labels_true).values
    row_ind, col_ind = linear_sum_assignment(-contingency)
    accuracy = contingency[row_ind, col_ind].sum() / len(labels_true)
    return accuracy


def calculate_all_metrics(labels_true, labels_pred, embeddings=None):
    """Calculate all clustering metrics."""
    metrics = {
        'ari': adjusted_rand_score(labels_true, labels_pred),
        'nmi': normalized_mutual_info_score(labels_true, labels_pred),
        'purity': calculate_purity(labels_true, labels_pred),
        'hungarian_acc': calculate_hungarian_accuracy(labels_true, labels_pred),
    }
    
    # Add V-measure components
    h, c, v = homogeneity_completeness_v_measure(labels_true, labels_pred)
    metrics['homogeneity'] = h
    metrics['completeness'] = c
    metrics['v_measure'] = v
    
    # Add silhouette if embeddings provided
    if embeddings is not None:
        try:
            metrics['silhouette'] = silhouette_score(embeddings, labels_pred, metric='cosine')
        except:
            metrics['silhouette'] = np.nan
    
    return metrics


def bootstrap_confidence_interval(labels_true, labels_pred, embeddings=None, 
                                  n_bootstrap=1000, confidence=0.95):
    """
    Calculate bootstrapped confidence intervals for all metrics.
    
    Args:
        labels_true: Ground truth labels
        labels_pred: Predicted cluster labels
        embeddings: Optional embeddings for silhouette score
        n_bootstrap: Number of bootstrap samples
        confidence: Confidence level (default 95%)
    
    Returns:
        Dictionary with mean, CI lower, CI upper for each metric
    """
    print(f"\nComputing bootstrapped confidence intervals ({n_bootstrap} samples)...")
    
    n_samples = len(labels_true)
    bootstrap_results = {
        'ari': [], 'nmi': [], 'purity': [], 'hungarian_acc': [],
        'homogeneity': [], 'completeness': [], 'v_measure': []
    }
    if embeddings is not None:
        bootstrap_results['silhouette'] = []
    
    for i in range(n_bootstrap):
        if (i + 1) % 200 == 0:
            print(f"  Bootstrap iteration {i+1}/{n_bootstrap}")
        
        # Resample with replacement
        indices = resample(range(n_samples), n_samples=n_samples, random_state=i)
        
        y_true_boot = labels_true[indices]
        y_pred_boot = labels_pred[indices]
        emb_boot = embeddings[indices] if embeddings is not None else None
        
        # Calculate metrics
        metrics = calculate_all_metrics(y_true_boot, y_pred_boot, emb_boot)
        
        for key, value in metrics.items():
            if not np.isnan(value):
                bootstrap_results[key].append(value)
    
    # Calculate confidence intervals
    alpha = 1 - confidence
    results = {}
    
    for metric, values in bootstrap_results.items():
        if len(values) > 0:
            values = np.array(values)
            mean = values.mean()
            ci_lower = np.percentile(values, 100 * alpha / 2)
            ci_upper = np.percentile(values, 100 * (1 - alpha / 2))
            
            results[metric] = {
                'mean': mean,
                'ci_lower': ci_lower,
                'ci_upper': ci_upper,
                'std': values.std()
            }
        else:
            results[metric] = {
                'mean': np.nan,
                'ci_lower': np.nan,
                'ci_upper': np.nan,
                'std': np.nan
            }
    
    return results


def permutation_test(labels_true, labels_pred_1, labels_pred_2, 
                     embeddings_1=None, embeddings_2=None,
                     n_permutations=10000, metric='ari'):
    """
    Permutation test to compare two clustering results.
    
    Args:
        labels_true: Ground truth labels
        labels_pred_1: Predictions from method 1
        labels_pred_2: Predictions from method 2
        embeddings_1: Embeddings for method 1 (optional)
        embeddings_2: Embeddings for method 2 (optional)
        n_permutations: Number of permutations
        metric: Metric to compare ('ari', 'nmi', 'purity', etc.)
    
    Returns:
        Dictionary with observed difference, p-value, and effect size
    """
    print(f"\nRunning permutation test ({n_permutations} permutations)...")
    print(f"Comparing: Method 1 vs Method 2 on metric '{metric}'")
    
    # Calculate observed difference
    metrics_1 = calculate_all_metrics(labels_true, labels_pred_1, embeddings_1)
    metrics_2 = calculate_all_metrics(labels_true, labels_pred_2, embeddings_2)
    
    observed_diff = metrics_1[metric] - metrics_2[metric]
    
    print(f"  Method 1 {metric}: {metrics_1[metric]:.4f}")
    print(f"  Method 2 {metric}: {metrics_2[metric]:.4f}")
    print(f"  Observed difference: {observed_diff:.4f}")
    
    # Permutation test
    n_samples = len(labels_true)
    perm_diffs = []
    
    # Combine predictions
    combined = np.column_stack([labels_pred_1, labels_pred_2])
    
    for i in range(n_permutations):
        if (i + 1) % 2000 == 0:
            print(f"  Permutation {i+1}/{n_permutations}")
        
        # Randomly assign to method 1 or 2
        np.random.seed(i)
        swap = np.random.randint(0, 2, size=n_samples)
        
        perm_pred_1 = np.where(swap == 0, combined[:, 0], combined[:, 1])
        perm_pred_2 = np.where(swap == 0, combined[:, 1], combined[:, 0])
        
        # Calculate metrics
        perm_metrics_1 = calculate_all_metrics(labels_true, perm_pred_1)
        perm_metrics_2 = calculate_all_metrics(labels_true, perm_pred_2)
        
        perm_diff = perm_metrics_1[metric] - perm_metrics_2[metric]
        perm_diffs.append(perm_diff)
    
    perm_diffs = np.array(perm_diffs)
    
    # Calculate p-value (two-tailed)
    p_value = np.mean(np.abs(perm_diffs) >= np.abs(observed_diff))
    
    # Calculate effect size (Cohen's d)
    pooled_std = np.std(perm_diffs)
    if pooled_std > 0:
        cohens_d = observed_diff / pooled_std
    else:
        cohens_d = np.nan
    
    results = {
        'metric': metric,
        'method_1_score': metrics_1[metric],
        'method_2_score': metrics_2[metric],
        'observed_difference': observed_diff,
        'p_value': p_value,
        'cohens_d': cohens_d,
        'significant': p_value < 0.05,
        'permutation_diffs_mean': perm_diffs.mean(),
        'permutation_diffs_std': perm_diffs.std()
    }
    
    return results


def identify_cluster_flippers(labels_pred_1, labels_pred_2, sequence_ids, 
                               top_n=50, output_file=None):
    """
    Identify sequences that switch clusters between two conditions.
    
    Args:
        labels_pred_1: Cluster assignments from condition 1
        labels_pred_2: Cluster assignments from condition 2
        sequence_ids: Sequence identifiers
        top_n: Number of flippers to report
        output_file: Optional file to save results
    
    Returns:
        DataFrame with cluster flippers
    """
    print(f"\nIdentifying cluster-flipping sequences...")
    
    # Find sequences that changed clusters
    flipped = labels_pred_1 != labels_pred_2
    flip_indices = np.where(flipped)[0]
    
    print(f"  Total sequences: {len(labels_pred_1):,}")
    print(f"  Cluster flippers: {len(flip_indices):,} ({len(flip_indices)/len(labels_pred_1)*100:.1f}%)")
    
    # Create DataFrame
    flippers_df = pd.DataFrame({
        'sequence_id': sequence_ids[flip_indices],
        'cluster_condition_1': labels_pred_1[flip_indices],
        'cluster_condition_2': labels_pred_2[flip_indices],
        'flip_index': flip_indices
    })
    
    # Sort by frequency of cluster pairs (most common flips first)
    flippers_df['flip_pair'] = flippers_df.apply(
        lambda x: f"{x['cluster_condition_1']}→{x['cluster_condition_2']}", 
        axis=1
    )
    
    flip_counts = flippers_df['flip_pair'].value_counts()
    flippers_df['flip_pair_frequency'] = flippers_df['flip_pair'].map(flip_counts)
    flippers_df = flippers_df.sort_values('flip_pair_frequency', ascending=False)
    
    # Report top flip patterns
    print(f"\n  Top {min(10, len(flip_counts))} cluster flip patterns:")
    for pattern, count in flip_counts.head(10).items():
        print(f"    {pattern}: {count:,} sequences")
    
    # Save if requested
    if output_file:
        flippers_df.head(top_n).to_csv(output_file, index=False)
        print(f"\n✅ Top {top_n} flippers saved to: {output_file}")
    
    return flippers_df


def main():
    parser = argparse.ArgumentParser(
        description='Statistical analysis for clustering results'
    )
    parser.add_argument(
        '--clustering-file',
        required=True,
        help='CSV with clustering assignments'
    )
    parser.add_argument(
        '--compare-file',
        help='Optional second clustering file for comparison'
    )
    parser.add_argument(
        '--embeddings-file',
        help='Optional embeddings file (.npy) for silhouette score'
    )
    parser.add_argument(
        '--n-bootstrap',
        type=int,
        default=1000,
        help='Number of bootstrap samples (default: 1000)'
    )
    parser.add_argument(
        '--n-permutations',
        type=int,
        default=10000,
        help='Number of permutations for test (default: 10000)'
    )
    parser.add_argument(
        '--output-dir',
        default='clustering_statistics',
        help='Output directory for results'
    )
    
    args = parser.parse_args()
    
    print("="*80)
    print("CLUSTERING STATISTICAL ANALYSIS")
    print("="*80)
    print()
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load clustering results
    print(f"Loading clustering results from {args.clustering_file}...")
    df = pd.read_csv(args.clustering_file)
    
    required_cols = ['uniprot_id', 'kinome_group_major', 'cluster']
    if not all(col in df.columns for col in required_cols):
        print(f"❌ Missing required columns: {required_cols}")
        sys.exit(1)
    
    labels_true = df['kinome_group_major'].values
    labels_pred_1 = df['cluster'].values
    sequence_ids = df['uniprot_id'].values
    
    print(f"✅ Loaded {len(df):,} sequences")
    print(f"   Classes: {len(set(labels_true)):,}")
    print(f"   Clusters: {len(set(labels_pred_1)):,}")
    
    # Load embeddings if available
    embeddings_1 = None
    if args.embeddings_file and os.path.exists(args.embeddings_file):
        print(f"\nLoading embeddings from {args.embeddings_file}...")
        embeddings_1 = np.load(args.embeddings_file)
        print(f"✅ Loaded embeddings: {embeddings_1.shape}")
    
    # 1. Bootstrapped confidence intervals
    print("\n" + "="*80)
    print("1. BOOTSTRAPPED CONFIDENCE INTERVALS")
    print("="*80)
    
    ci_results = bootstrap_confidence_interval(
        labels_true, labels_pred_1, embeddings_1,
        n_bootstrap=args.n_bootstrap
    )
    
    print("\nResults (95% confidence intervals):")
    print("-"*80)
    for metric, values in ci_results.items():
        print(f"{metric.upper():15} {values['mean']:.4f} [{values['ci_lower']:.4f}, {values['ci_upper']:.4f}] ±{values['std']:.4f}")
    
    # Save CI results
    ci_df = pd.DataFrame(ci_results).T
    ci_file = f"{args.output_dir}/confidence_intervals.csv"
    ci_df.to_csv(ci_file)
    print(f"\n✅ Confidence intervals saved to: {ci_file}")
    
    # 2. Permutation test (if comparison file provided)
    if args.compare_file:
        print("\n" + "="*80)
        print("2. PERMUTATION TEST")
        print("="*80)
        
        print(f"\nLoading comparison clustering from {args.compare_file}...")
        df_compare = pd.read_csv(args.compare_file)
        labels_pred_2 = df_compare['cluster'].values
        
        # Load comparison embeddings if available
        embeddings_2 = None
        compare_emb_file = args.compare_file.replace('.csv', '_embeddings.npy')
        if os.path.exists(compare_emb_file):
            embeddings_2 = np.load(compare_emb_file)
        
        # Run permutation tests for all metrics
        perm_results = []
        for metric in ['ari', 'nmi', 'purity', 'hungarian_acc']:
            result = permutation_test(
                labels_true, labels_pred_1, labels_pred_2,
                embeddings_1, embeddings_2,
                n_permutations=args.n_permutations,
                metric=metric
            )
            perm_results.append(result)
        
        # Print results
        print("\nPermutation test results:")
        print("-"*80)
        print(f"{'Metric':<15} {'Method 1':<10} {'Method 2':<10} {'Diff':<10} {'p-value':<10} {'Cohen d':<10} {'Sig?'}")
        print("-"*80)
        for result in perm_results:
            sig_marker = "***" if result['p_value'] < 0.001 else "**" if result['p_value'] < 0.01 else "*" if result['p_value'] < 0.05 else "ns"
            print(f"{result['metric'].upper():<15} {result['method_1_score']:<10.4f} {result['method_2_score']:<10.4f} "
                  f"{result['observed_difference']:<10.4f} {result['p_value']:<10.4f} {result['cohens_d']:<10.2f} {sig_marker}")
        
        # Save permutation results
        perm_df = pd.DataFrame(perm_results)
        perm_file = f"{args.output_dir}/permutation_test.csv"
        perm_df.to_csv(perm_file, index=False)
        print(f"\n✅ Permutation test results saved to: {perm_file}")
        
        # 3. Identify cluster flippers
        print("\n" + "="*80)
        print("3. CLUSTER-FLIPPING SEQUENCES")
        print("="*80)
        
        flippers_file = f"{args.output_dir}/cluster_flippers.csv"
        flippers_df = identify_cluster_flippers(
            labels_pred_1, labels_pred_2, sequence_ids,
            top_n=50, output_file=flippers_file
        )
    
    print("\n" + "="*80)
    print("✅ STATISTICAL ANALYSIS COMPLETE!")
    print("="*80)
    print()


if __name__ == "__main__":
    main()

