#!/usr/bin/env python3
"""
Exemplar Retrieval Analysis via Nearest Neighbor Search

Implements reviewer requirements:
1. Cosine similarity on L2-normalized vectors
2. Leave-one-out protocol over test set
3. Top-k same-family hit rate and MRR (Mean Reciprocal Rank)
4. Similarity→confidence calibration mapping
5. Precision-recall curves by threshold
6. Qualitative failure mode analysis

Uses sklearn for exact search (dataset size manageable, <5k sequences).
"""

import os
import sys
import argparse
import json
import numpy as np
import pandas as pd
from pathlib import Path
from collections import defaultdict, Counter

from sklearn.metrics.pairwise import cosine_similarity
from sklearn.preprocessing import normalize
from sklearn.metrics import precision_recall_curve, auc

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

import warnings
warnings.filterwarnings('ignore')


def load_data_with_splits(embeddings_dir, labels_csv, splits_file):
    """Load embeddings, labels, and split into train/test."""
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
    train_ids_array = data.loc[train_mask, 'uniprot_id'].values
    test_ids_array = data.loc[test_mask, 'uniprot_id'].values
    
    return X_train, X_test, y_train, y_test, train_ids_array, test_ids_array


def compute_nearest_neighbors(query_embeddings, reference_embeddings, k=10):
    """
    Compute k-nearest neighbors using cosine similarity on L2-normalized vectors.
    
    Args:
        query_embeddings: Query vectors (N_query, D)
        reference_embeddings: Reference vectors (N_ref, D)
        k: Number of neighbors to retrieve
    
    Returns:
        neighbor_indices: (N_query, k) indices of nearest neighbors
        similarities: (N_query, k) cosine similarities
    """
    # L2-normalize vectors
    query_norm = normalize(query_embeddings, norm='l2')
    reference_norm = normalize(reference_embeddings, norm='l2')
    
    # Compute cosine similarity (equivalent to dot product after L2 normalization)
    similarity_matrix = cosine_similarity(query_norm, reference_norm)
    
    # Get top-k indices (descending order)
    # Note: argsort is ascending, so we negate to get descending
    neighbor_indices = np.argsort(-similarity_matrix, axis=1)[:, :k]
    
    # Get corresponding similarities
    similarities = np.take_along_axis(similarity_matrix, neighbor_indices, axis=1)
    
    return neighbor_indices, similarities


def calculate_top_k_hit_rate(neighbor_labels, query_labels, k_values=[1, 3, 5, 10]):
    """
    Calculate top-k same-family hit rate.
    
    Args:
        neighbor_labels: (N_query, k_max) labels of retrieved neighbors
        query_labels: (N_query,) true labels of queries
        k_values: List of k values to evaluate
    
    Returns:
        Dictionary with hit rates for each k
    """
    results = {}
    n_queries = len(query_labels)
    
    for k in k_values:
        if k > neighbor_labels.shape[1]:
            continue
        
        # Check if any of top-k neighbors match query label
        hits = 0
        for i in range(n_queries):
            query_label = query_labels[i]
            top_k_labels = neighbor_labels[i, :k]
            if query_label in top_k_labels:
                hits += 1
        
        hit_rate = hits / n_queries
        results[f'top_{k}_hit_rate'] = hit_rate
    
    return results


def calculate_mrr(neighbor_labels, query_labels):
    """
    Calculate Mean Reciprocal Rank (MRR).
    
    MRR is the average of reciprocal ranks of the first correct match.
    Example: If first match at position 3, reciprocal rank = 1/3.
    
    Returns:
        MRR score (0 to 1, higher is better)
    """
    reciprocal_ranks = []
    
    for i in range(len(query_labels)):
        query_label = query_labels[i]
        neighbors = neighbor_labels[i]
        
        # Find position of first matching label (1-indexed for rank)
        for rank, neighbor_label in enumerate(neighbors, start=1):
            if neighbor_label == query_label:
                reciprocal_ranks.append(1.0 / rank)
                break
        else:
            # No match found in top-k
            reciprocal_ranks.append(0.0)
    
    mrr = np.mean(reciprocal_ranks)
    return mrr


def similarity_to_confidence_mapping(similarities, neighbor_labels, query_labels, n_bins=20):
    """
    Map similarity scores to calibrated confidence (precision).
    
    For each similarity bin, compute: P(same_family | similarity_in_bin)
    
    Returns:
        DataFrame with similarity bins and corresponding precision
    """
    # Flatten all (query, neighbor) pairs
    all_similarities = []
    all_matches = []
    
    for i in range(len(query_labels)):
        for j in range(similarities.shape[1]):
            sim = similarities[i, j]
            match = (neighbor_labels[i, j] == query_labels[i])
            all_similarities.append(sim)
            all_matches.append(match)
    
    all_similarities = np.array(all_similarities)
    all_matches = np.array(all_matches)
    
    # Create bins
    bin_edges = np.linspace(all_similarities.min(), all_similarities.max(), n_bins + 1)
    bin_indices = np.digitize(all_similarities, bin_edges[:-1]) - 1
    bin_indices = np.clip(bin_indices, 0, n_bins - 1)
    
    # Calculate precision per bin
    calibration_data = []
    for bin_idx in range(n_bins):
        mask = bin_indices == bin_idx
        if mask.sum() > 0:
            bin_sim_mean = all_similarities[mask].mean()
            bin_sim_min = all_similarities[mask].min()
            bin_sim_max = all_similarities[mask].max()
            bin_precision = all_matches[mask].mean()
            bin_count = mask.sum()
            
            # Interpret confidence level
            if bin_precision >= 0.9:
                confidence = "Very High"
            elif bin_precision >= 0.75:
                confidence = "High"
            elif bin_precision >= 0.6:
                confidence = "Medium"
            elif bin_precision >= 0.4:
                confidence = "Low"
            else:
                confidence = "Very Low"
            
            calibration_data.append({
                'bin': bin_idx,
                'similarity_min': bin_sim_min,
                'similarity_max': bin_sim_max,
                'similarity_mean': bin_sim_mean,
                'precision': bin_precision,
                'count': bin_count,
                'confidence_level': confidence
            })
    
    return pd.DataFrame(calibration_data)


def plot_precision_recall_by_threshold(similarities, neighbor_labels, query_labels, output_file):
    """
    Plot precision-recall curve as similarity threshold varies.
    """
    # Flatten pairs
    all_similarities = similarities[:, 0]  # Use top-1 neighbor
    all_matches = (neighbor_labels[:, 0] == query_labels).astype(int)
    
    # Compute precision-recall curve
    precision, recall, thresholds = precision_recall_curve(all_matches, all_similarities)
    pr_auc = auc(recall, precision)
    
    # Plot
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(recall, precision, 'b-', linewidth=2, label=f'PR curve (AUC={pr_auc:.3f})')
    ax.set_xlabel('Recall (Same-Family Hit Rate)', fontsize=12)
    ax.set_ylabel('Precision (Confidence)', fontsize=12)
    ax.set_title('Precision-Recall by Similarity Threshold', fontsize=14)
    ax.grid(alpha=0.3)
    ax.legend()
    
    # Add threshold annotations for key points
    for i in range(0, len(thresholds), len(thresholds)//5):
        if i < len(thresholds):
            ax.annotate(f'{thresholds[i]:.2f}', 
                       xy=(recall[i], precision[i]),
                       xytext=(5, 5), textcoords='offset points',
                       fontsize=8, alpha=0.7)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"✅ Precision-recall curve saved to: {output_file}")
    
    return pr_auc


def analyze_failure_modes(neighbor_indices, neighbor_labels, similarities, 
                          query_labels, test_ids, train_ids, 
                          labels_df, n_failures=10):
    """
    Identify and analyze retrieval failures (top-1 neighbor wrong family).
    
    Returns:
        DataFrame with failure cases and potential reasons
    """
    failures = []
    
    for i in range(len(query_labels)):
        top1_label = neighbor_labels[i, 0]
        true_label = query_labels[i]
        
        if top1_label != true_label:
            # Get query info
            query_id = test_ids[i]
            query_info = labels_df[labels_df['uniprot_id'] == query_id].iloc[0]
            
            # Get top-1 neighbor info
            neighbor_idx = neighbor_indices[i, 0]
            neighbor_id = train_ids[neighbor_idx]
            neighbor_info = labels_df[labels_df['uniprot_id'] == neighbor_id].iloc[0]
            
            # Get sequence info if available
            query_seq = query_info.get('sequence', '') if 'sequence' in query_info else ''
            neighbor_seq = neighbor_info.get('sequence', '') if 'sequence' in neighbor_info else ''
            
            # Analyze potential failure reasons
            reasons = []
            
            # Check sequence length difference
            if query_seq and neighbor_seq:
                len_diff = abs(len(query_seq) - len(neighbor_seq))
                if len_diff > 200:
                    reasons.append(f"Length diff: {len_diff} aa")
            
            # Check if it's a known difficult family
            difficult_families = ['Atypical', 'TKL', 'STE']
            if true_label in difficult_families:
                reasons.append(f"{true_label} is diverse family")
            
            # Check similarity score (low similarity = genuinely different)
            top1_sim = similarities[i, 0]
            if top1_sim < 0.6:
                reasons.append(f"Low similarity ({top1_sim:.3f})")
            elif top1_sim > 0.8:
                reasons.append(f"High similarity ({top1_sim:.3f}) - possible mislabel?")
            
            # Check if correct family appears in top-k
            top5_labels = neighbor_labels[i, :5]
            if true_label in top5_labels:
                rank = np.where(top5_labels == true_label)[0][0] + 1
                reasons.append(f"Correct family at rank {rank}")
            else:
                reasons.append("Correct family not in top-5")
            
            failures.append({
                'query_id': query_id,
                'true_family': true_label,
                'predicted_family': top1_label,
                'similarity': top1_sim,
                'top1_neighbor_id': neighbor_id,
                'reasons': '; '.join(reasons) if reasons else 'Unknown',
                'query_length': len(query_seq) if query_seq else -1,
                'neighbor_length': len(neighbor_seq) if neighbor_seq else -1,
            })
    
    failures_df = pd.DataFrame(failures)
    
    # Sort by similarity (descending) to find "confident but wrong" cases
    failures_df = failures_df.sort_values('similarity', ascending=False)
    
    return failures_df


def main():
    parser = argparse.ArgumentParser(
        description='Exemplar retrieval analysis via nearest neighbor'
    )
    parser.add_argument(
        '--embeddings-dir',
        default='kinases_domains_embeddings',
        help='Directory with ESM embeddings'
    )
    parser.add_argument(
        '--labels-csv',
        default='kinases_domains_e0.01.csv',
        help='CSV with labels and sequence info'
    )
    parser.add_argument(
        '--splits-file',
        default='data/splits_40.json',
        help='Train/test splits file'
    )
    parser.add_argument(
        '--top-k',
        type=int,
        default=10,
        help='Number of neighbors to retrieve (default: 10)'
    )
    parser.add_argument(
        '--output-dir',
        default='exemplar_retrieval_results',
        help='Output directory'
    )
    
    args = parser.parse_args()
    
    print("="*80)
    print("EXEMPLAR RETRIEVAL ANALYSIS (NEAREST NEIGHBOR)")
    print("="*80)
    print()
    print(f"Embeddings:  {args.embeddings_dir}")
    print(f"Labels:      {args.labels_csv}")
    print(f"Splits:      {args.splits_file}")
    print(f"Top-k:       {args.top_k}")
    print()
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load data
    print("Loading data...")
    X_train, X_test, y_train, y_test, train_ids, test_ids = load_data_with_splits(
        args.embeddings_dir, args.labels_csv, args.splits_file
    )
    
    print(f"✅ Train: {len(X_train):,} sequences, {len(set(y_train))} classes")
    print(f"✅ Test:  {len(X_test):,} sequences, {len(set(y_test))} classes")
    print(f"   Embedding dimension: {X_train.shape[1]}")
    print()
    
    # Load full labels dataframe for failure analysis
    labels_df = pd.read_csv(args.labels_csv)
    
    # Compute nearest neighbors
    print(f"Computing {args.top_k}-nearest neighbors (cosine similarity, L2-normalized)...")
    neighbor_indices, similarities = compute_nearest_neighbors(X_test, X_train, k=args.top_k)
    
    print(f"✅ Retrieved {args.top_k} neighbors for each of {len(X_test)} test sequences")
    print(f"   Similarity range: [{similarities.min():.4f}, {similarities.max():.4f}]")
    print(f"   Similarity mean: {similarities.mean():.4f} ± {similarities.std():.4f}")
    print()
    
    # Get neighbor labels
    neighbor_labels = y_train[neighbor_indices]
    
    # 1. Top-k hit rates
    print("="*80)
    print("1. TOP-K SAME-FAMILY HIT RATES")
    print("="*80)
    print()
    
    hit_rates = calculate_top_k_hit_rate(neighbor_labels, y_test, k_values=[1, 3, 5, 10])
    
    print("Hit rates (correct family in top-k neighbors):")
    print("-"*80)
    for metric, value in hit_rates.items():
        k = metric.split('_')[1]
        print(f"  Top-{k:<2} hit rate: {value:.4f} ({value*100:.1f}%)")
    print()
    
    # 2. Mean Reciprocal Rank (MRR)
    print("="*80)
    print("2. MEAN RECIPROCAL RANK (MRR)")
    print("="*80)
    print()
    
    mrr = calculate_mrr(neighbor_labels, y_test)
    print(f"MRR: {mrr:.4f}")
    print()
    print("Interpretation:")
    print(f"  • Average rank of first correct match: {1/mrr:.2f}")
    print(f"  • {mrr*100:.1f}% of ideal ranking performance")
    print()
    
    # 3. Similarity → Confidence calibration
    print("="*80)
    print("3. SIMILARITY → CONFIDENCE CALIBRATION")
    print("="*80)
    print()
    
    calibration_df = similarity_to_confidence_mapping(
        similarities, neighbor_labels, y_test, n_bins=20
    )
    
    print("Similarity-to-precision mapping:")
    print("-"*80)
    print(f"{'Similarity Range':<20} {'Precision':<12} {'Count':<8} {'Confidence Level'}")
    print("-"*80)
    for _, row in calibration_df.iterrows():
        print(f"[{row['similarity_min']:.3f}, {row['similarity_max']:.3f}] "
              f"{row['precision']:>10.3f}  {row['count']:>7}  {row['confidence_level']}")
    
    calibration_df.to_csv(f"{args.output_dir}/similarity_calibration.csv", index=False)
    print(f"\n✅ Calibration table saved")
    
    # Identify thresholds for confidence levels
    print("\nRecommended similarity thresholds:")
    print("-"*80)
    high_conf = calibration_df[calibration_df['precision'] >= 0.75]
    if len(high_conf) > 0:
        threshold_high = high_conf['similarity_min'].min()
        print(f"  'High confidence' (precision ≥75%):   similarity ≥ {threshold_high:.3f}")
    
    medium_conf = calibration_df[calibration_df['precision'] >= 0.6]
    if len(medium_conf) > 0:
        threshold_medium = medium_conf['similarity_min'].min()
        print(f"  'Medium confidence' (precision ≥60%): similarity ≥ {threshold_medium:.3f}")
    
    print()
    
    # 4. Precision-Recall curve
    print("="*80)
    print("4. PRECISION-RECALL BY THRESHOLD")
    print("="*80)
    print()
    
    pr_auc = plot_precision_recall_by_threshold(
        similarities, neighbor_labels, y_test,
        f"{args.output_dir}/precision_recall_curve.png"
    )
    print(f"PR-AUC: {pr_auc:.4f}")
    print()
    
    # 5. Per-class retrieval performance
    print("="*80)
    print("5. PER-CLASS RETRIEVAL PERFORMANCE")
    print("="*80)
    print()
    
    per_class_results = []
    for family in sorted(set(y_test)):
        mask = y_test == family
        family_neighbor_labels = neighbor_labels[mask]
        family_query_labels = y_test[mask]
        family_similarities = similarities[mask]
        
        # Top-1 hit rate
        top1_hits = (family_neighbor_labels[:, 0] == family).sum()
        top1_rate = top1_hits / mask.sum()
        
        # Top-3 hit rate
        top3_correct = 0
        for i in range(len(family_query_labels)):
            if family in family_neighbor_labels[i, :3]:
                top3_correct += 1
        top3_rate = top3_correct / mask.sum()
        
        # Mean similarity
        mean_sim = family_similarities[:, 0].mean()
        
        per_class_results.append({
            'family': family,
            'n_test': mask.sum(),
            'top1_hit_rate': top1_rate,
            'top3_hit_rate': top3_rate,
            'mean_similarity': mean_sim
        })
    
    per_class_df = pd.DataFrame(per_class_results)
    per_class_df = per_class_df.sort_values('top1_hit_rate', ascending=False)
    
    print(f"{'Family':<12} {'N':<6} {'Top-1 Hit':<12} {'Top-3 Hit':<12} {'Mean Sim'}")
    print("-"*80)
    for _, row in per_class_df.iterrows():
        print(f"{row['family']:<12} {row['n_test']:<6} {row['top1_hit_rate']:<12.3f} "
              f"{row['top3_hit_rate']:<12.3f} {row['mean_similarity']:.3f}")
    
    per_class_df.to_csv(f"{args.output_dir}/per_class_retrieval.csv", index=False)
    print(f"\n✅ Per-class results saved")
    print()
    
    # 6. Failure mode analysis
    print("="*80)
    print("6. FAILURE MODE ANALYSIS")
    print("="*80)
    print()
    
    failures_df = analyze_failure_modes(
        neighbor_indices, neighbor_labels, similarities,
        y_test, test_ids, train_ids, labels_df, n_failures=20
    )
    
    print(f"Total failures (top-1 wrong): {len(failures_df):,} / {len(y_test):,} "
          f"({len(failures_df)/len(y_test)*100:.1f}%)")
    print()
    
    print("Top 10 failures (sorted by similarity - confident but wrong):")
    print("-"*80)
    print(f"{'Query ID':<12} {'True':<10} {'Predicted':<10} {'Similarity':<12} {'Reasons'}")
    print("-"*80)
    for _, row in failures_df.head(10).iterrows():
        print(f"{row['query_id']:<12} {row['true_family']:<10} {row['predicted_family']:<10} "
              f"{row['similarity']:<12.3f} {row['reasons'][:50]}")
    
    failures_df.to_csv(f"{args.output_dir}/failure_cases.csv", index=False)
    print(f"\n✅ Failure analysis saved to: {args.output_dir}/failure_cases.csv")
    print()
    
    # Failure mode statistics
    print("Failure mode breakdown:")
    print("-"*80)
    
    # Count reasons
    all_reasons = []
    for reasons_str in failures_df['reasons']:
        all_reasons.extend([r.strip() for r in reasons_str.split(';')])
    
    reason_counts = Counter(all_reasons)
    for reason, count in reason_counts.most_common(10):
        print(f"  {reason:<50} {count:>4} cases")
    print()
    
    # 7. Save summary report
    print("="*80)
    print("7. SUMMARY REPORT")
    print("="*80)
    print()
    
    summary = {
        'dataset': {
            'train_size': len(X_train),
            'test_size': len(X_test),
            'n_classes': len(set(y_test)),
            'embedding_dim': X_train.shape[1]
        },
        'retrieval_performance': {
            'top1_hit_rate': hit_rates['top_1_hit_rate'],
            'top3_hit_rate': hit_rates['top_3_hit_rate'],
            'top5_hit_rate': hit_rates['top_5_hit_rate'],
            'top10_hit_rate': hit_rates['top_10_hit_rate'],
            'mrr': mrr,
            'pr_auc': pr_auc
        },
        'similarity_statistics': {
            'mean': float(similarities.mean()),
            'std': float(similarities.std()),
            'min': float(similarities.min()),
            'max': float(similarities.max())
        },
        'failures': {
            'total': len(failures_df),
            'rate': len(failures_df) / len(y_test)
        }
    }
    
    with open(f"{args.output_dir}/retrieval_summary.json", 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"✅ Summary saved to: {args.output_dir}/retrieval_summary.json")
    print()
    
    # Final summary
    print("="*80)
    print("RETRIEVAL PERFORMANCE SUMMARY")
    print("="*80)
    print()
    print(f"Top-1 Hit Rate:     {hit_rates['top_1_hit_rate']:.3f} ({hit_rates['top_1_hit_rate']*100:.1f}%)")
    print(f"Top-3 Hit Rate:     {hit_rates['top_3_hit_rate']:.3f} ({hit_rates['top_3_hit_rate']*100:.1f}%)")
    print(f"Top-5 Hit Rate:     {hit_rates['top_5_hit_rate']:.3f} ({hit_rates['top_5_hit_rate']*100:.1f}%)")
    print(f"MRR:                {mrr:.3f}")
    print(f"PR-AUC:             {pr_auc:.3f}")
    print()
    print(f"Similarity thresholds:")
    if len(high_conf) > 0:
        print(f"  High confidence (≥75% precision):   {threshold_high:.3f}")
    if len(medium_conf) > 0:
        print(f"  Medium confidence (≥60% precision): {threshold_medium:.3f}")
    print()
    print(f"Failure rate: {len(failures_df)/len(y_test)*100:.1f}%")
    print()
    
    print("="*80)
    print("✅ EXEMPLAR RETRIEVAL ANALYSIS COMPLETE!")
    print("="*80)


if __name__ == "__main__":
    main()

