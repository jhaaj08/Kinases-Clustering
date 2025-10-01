#!/usr/bin/env python3
"""
Fuse domain ESM-2 embeddings with motif features and evaluate clustering.

Pipeline:
1. Load domain ESM-2 embeddings (1280-d)
2. Load motif features (22-d)
3. Concatenate features (early fusion)
4. Standardize combined features
5. Run K-Means clustering (k=10, excluding "Other")
6. Compare metrics with domain-only embeddings
"""

import os
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (
    adjusted_rand_score,
    normalized_mutual_info_score,
    homogeneity_completeness_v_measure,
    confusion_matrix,
    silhouette_score,
)
from scipy.optimize import linear_sum_assignment


def load_embeddings(emb_dir="kinases_domains_embeddings"):
    """Load ESM-2 domain embeddings."""
    emb_file = os.path.join(emb_dir, "esm2_embeddings.npy")
    idx_file = os.path.join(emb_dir, "esm2_index.csv")
    
    X = np.load(emb_file)
    ids = pd.read_csv(idx_file)["uniprot_id"].astype(str).values
    
    return X, ids


def load_motif_features(motif_file="kinases_domains_with_motifs.csv"):
    """Load motif features from CSV."""
    df = pd.read_csv(motif_file)
    
    # Select motif feature columns (exclude metadata and sequence)
    exclude_cols = [
        'uniprot_id', 'protein_name', 'function', 'kinome_group_subfamily',
        'kinome_group_major', 'family_slim', 'conformation_DFG_aC',
        'inhibitor_class_sensitivity', 'sequence', 'env_from', 'env_to',
        'domain_length', 'evalue', 'score'
    ]
    
    motif_cols = [col for col in df.columns if col not in exclude_cols]
    
    return df, motif_cols


def fuse_features(emb_ids, embeddings, motif_df, motif_cols):
    """
    Fuse embeddings with motif features.
    
    Returns:
        X_fused: Combined feature matrix
        ids_aligned: Aligned UniProt IDs
        y_labels: Kinome group labels
    """
    print("Fusing embeddings with motif features...")
    
    # Create embedding dataframe
    emb_df = pd.DataFrame({
        'uniprot_id': emb_ids,
        'emb_idx': range(len(emb_ids))
    })
    
    # Merge with motif features
    merged = emb_df.merge(
        motif_df[['uniprot_id'] + motif_cols + ['kinome_group_major']],
        on='uniprot_id',
        how='inner'
    )
    
    print(f"  Embeddings: {len(emb_ids)} sequences")
    print(f"  Motif data: {len(motif_df)} sequences")
    print(f"  Merged: {len(merged)} sequences")
    
    # Get embeddings for matched sequences
    emb_indices = merged['emb_idx'].values
    X_emb = embeddings[emb_indices]  # (N, 1280)
    
    # Get motif features
    X_motif = merged[motif_cols].values  # (N, 22)
    
    # Concatenate (early fusion)
    X_fused = np.concatenate([X_emb, X_motif], axis=1)  # (N, 1302)
    
    ids_aligned = merged['uniprot_id'].values
    y_labels = merged['kinome_group_major'].values
    
    print(f"  Fused feature shape: {X_fused.shape}")
    print(f"    - ESM-2 embeddings: {X_emb.shape[1]} dims")
    print(f"    - Motif features: {X_motif.shape[1]} dims")
    print(f"    - Total: {X_fused.shape[1]} dims")
    
    return X_fused, ids_aligned, y_labels


def run_clustering_experiment(X, ids, y_labels, experiment_name, k=10):
    """
    Run K-Means clustering and compute all metrics.
    
    Returns:
        Dictionary of metrics and predictions
    """
    print(f"\nRunning clustering experiment: {experiment_name}")
    print("-"*80)
    
    # Exclude "Other" category
    mask = y_labels != "Other"
    X_clean = X[mask]
    ids_clean = ids[mask]
    y_clean = y_labels[mask]
    
    print(f"  Dataset: {len(X_clean):,} samples (excluded {(~mask).sum()} 'Other')")
    print(f"  Features: {X_clean.shape[1]} dimensions")
    print(f"  Labels: {len(np.unique(y_clean))} kinase groups")
    
    # Standardize features
    print(f"  Standardizing features...")
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_clean)
    
    # K-Means
    print(f"  Running K-Means (k={k}, n_init=50)...")
    km = KMeans(n_clusters=k, random_state=42, n_init=50, max_iter=500)
    pred = km.fit_predict(X_scaled)
    
    # Metrics
    print(f"  Computing metrics...")
    ari = adjusted_rand_score(y_clean, pred)
    nmi = normalized_mutual_info_score(y_clean, pred)
    h, c, v = homogeneity_completeness_v_measure(y_clean, pred)
    sil = silhouette_score(X_scaled, pred, sample_size=min(10000, len(X_scaled)), random_state=42)
    
    # Confusion matrix and purity
    labels = sorted(np.unique(y_clean))
    label_to_idx = {lbl: i for i, lbl in enumerate(labels)}
    y_idx = np.array([label_to_idx[lbl] for lbl in y_clean])
    cm = confusion_matrix(y_idx, pred)
    
    purity = cm.max(axis=0).sum() / cm.sum()
    
    # Hungarian
    cost = cm.max() - cm
    row_ind, col_ind = linear_sum_assignment(cost)
    hungarian_acc = cm[row_ind, col_ind].sum() / cm.sum()
    
    # Per-cluster analysis
    cluster_analysis = []
    for c_id in range(k):
        mask_c = pred == c_id
        counts = pd.Series(y_clean[mask_c]).value_counts()
        
        if len(counts) > 0:
            top_label = counts.index[0]
            top_count = counts.iloc[0]
            total = mask_c.sum()
            
            top3 = ", ".join([f"{lbl}({cnt})" for lbl, cnt in counts.head(3).items()])
            
            cluster_analysis.append({
                'cluster': c_id,
                'size': int(total),
                'majority_label': top_label,
                'majority_count': int(top_count),
                'majority_pct': round(top_count / total * 100, 1),
                'top3': top3
            })
    
    cluster_df = pd.DataFrame(cluster_analysis).sort_values('majority_pct', ascending=False)
    
    results = {
        'experiment': experiment_name,
        'n_samples': len(X_clean),
        'n_features': X_clean.shape[1],
        'ARI': ari,
        'NMI': nmi,
        'Homogeneity': h,
        'Completeness': c,
        'V-measure': v,
        'Silhouette': sil,
        'Purity': purity,
        'Hungarian': hungarian_acc,
        'predictions': pred,
        'ids': ids_clean,
        'labels': y_clean,
        'cluster_analysis': cluster_df,
        'best_cluster_purity': cluster_df['majority_pct'].max() / 100.0
    }
    
    print(f"\n  ✅ Clustering complete!")
    
    return results


def print_results(results):
    """Print results in a nice format."""
    
    print("\n" + "="*80)
    print(f"RESULTS: {results['experiment']}")
    print("="*80)
    
    print(f"\nDataset: {results['n_samples']:,} samples, {results['n_features']} features")
    print()
    
    print("Metrics:")
    print("-"*80)
    print(f"  ARI (Adjusted Rand Index):      {results['ARI']:.4f}")
    print(f"  NMI (Normalized Mutual Info):   {results['NMI']:.4f}")
    print(f"  Homogeneity:                     {results['Homogeneity']:.4f}")
    print(f"  Completeness:                    {results['Completeness']:.4f}")
    print(f"  V-measure:                       {results['V-measure']:.4f}")
    print(f"  Silhouette Score:                {results['Silhouette']:.4f}")
    print(f"  Purity:                          {results['Purity']:.4f}")
    print(f"  Hungarian Accuracy:              {results['Hungarian']:.4f}")
    print(f"  Best Cluster Purity:             {results['best_cluster_purity']:.1%}")
    print()
    
    print("Top 5 Clusters by Purity:")
    print("-"*80)
    top5 = results['cluster_analysis'].head(5)
    for _, row in top5.iterrows():
        print(f"  Cluster {row['cluster']:2d} (n={row['size']:3d}): {row['majority_label']:12s} {row['majority_pct']:5.1f}%  |  {row['top3']}")
    print()


def compare_results(baseline, enhanced):
    """Compare two sets of results."""
    
    print("\n" + "╔" + "="*78 + "╗")
    print("║" + " "*20 + "DOMAIN-ONLY vs DOMAIN+MOTIFS" + " "*30 + "║")
    print("╚" + "="*78 + "╝")
    print()
    
    metrics = ['ARI', 'NMI', 'Homogeneity', 'Completeness', 'V-measure', 
               'Silhouette', 'Purity', 'Hungarian', 'best_cluster_purity']
    
    print(f"{'Metric':<25} {'Domain-Only':<15} {'Domain+Motifs':<15} {'Change':<20}")
    print("-"*80)
    
    for metric in metrics:
        base_val = baseline[metric]
        enh_val = enhanced[metric]
        change = enh_val - base_val
        pct_change = (change / base_val * 100) if base_val != 0 else 0
        
        sign = "+" if change > 0 else ""
        emoji = "✅" if change > 0 else ("⚠️" if change < 0 else "➖")
        
        print(f"{metric:<25} {base_val:<15.4f} {enh_val:<15.4f} {sign}{change:.4f} ({sign}{pct_change:+.1f}%) {emoji}")
    
    print()


def save_results(results, output_dir="clustering"):
    """Save clustering results to files."""
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Save assignments
    assign_df = pd.DataFrame({
        'uniprot_id': results['ids'],
        'cluster': results['predictions'],
        'kinome_group_major': results['labels']
    })
    
    assign_file = os.path.join(output_dir, "kmeans10_domain_motifs_assignments.csv")
    assign_df.to_csv(assign_file, index=False)
    
    # Save report
    report_file = os.path.join(output_dir, "kmeans10_domain_motifs_report.txt")
    with open(report_file, 'w') as f:
        f.write(f"K-Means Clustering: {results['experiment']}\n")
        f.write("="*80 + "\n\n")
        f.write(f"Dataset: {results['n_samples']:,} samples\n")
        f.write(f"Features: {results['n_features']} dimensions\n\n")
        
        f.write("Metrics:\n")
        f.write("-"*80 + "\n")
        for metric in ['ARI', 'NMI', 'Homogeneity', 'Completeness', 'V-measure', 
                       'Silhouette', 'Purity', 'Hungarian']:
            f.write(f"  {metric:<25s}: {results[metric]:.4f}\n")
        
        f.write(f"\nBest Cluster Purity: {results['best_cluster_purity']:.1%}\n\n")
        
        f.write("Per-Cluster Analysis:\n")
        f.write("-"*80 + "\n")
        f.write(results['cluster_analysis'][['cluster', 'size', 'majority_label', 
                                             'majority_count', 'majority_pct']].to_string(index=False))
    
    print(f"  Saved: {assign_file}")
    print(f"  Saved: {report_file}")
    
    return assign_file, report_file


def main():
    """Main pipeline for motif-fused clustering."""
    
    print("="*80)
    print("DOMAIN + MOTIF FEATURES: CLUSTERING EXPERIMENT")
    print("="*80)
    print()
    
    # Load embeddings
    print("Loading domain ESM-2 embeddings...")
    X_emb, emb_ids = load_embeddings()
    print(f"✅ Loaded {len(X_emb):,} embeddings ({X_emb.shape[1]} dims)")
    
    # Load motif features
    print("\nLoading motif features...")
    motif_df, motif_cols = load_motif_features()
    print(f"✅ Loaded {len(motif_df):,} sequences with {len(motif_cols)} motif features")
    
    # Fuse features
    print()
    X_fused, ids_aligned, y_labels = fuse_features(emb_ids, X_emb, motif_df, motif_cols)
    
    # Run clustering with fused features
    print()
    results_fused = run_clustering_experiment(
        X_fused, ids_aligned, y_labels,
        experiment_name="Domain ESM-2 + Motif Features"
    )
    
    print_results(results_fused)
    
    # Also run domain-only for comparison (same dataset)
    print("\n" + "="*80)
    print("Re-running domain-only clustering for fair comparison...")
    print("="*80)
    
    # Extract just embeddings for the same IDs
    emb_df = pd.DataFrame({'uniprot_id': emb_ids, 'emb_idx': range(len(emb_ids))})
    merged_ids = emb_df[emb_df['uniprot_id'].isin(ids_aligned)]
    X_emb_aligned = X_emb[merged_ids['emb_idx'].values]
    
    results_domain_only = run_clustering_experiment(
        X_emb_aligned, ids_aligned, y_labels,
        experiment_name="Domain ESM-2 Only"
    )
    
    print_results(results_domain_only)
    
    # Compare results
    compare_results(results_domain_only, results_fused)
    
    # Save results
    print("\nSaving results...")
    save_results(results_fused)
    
    # Final summary
    print("\n" + "="*80)
    print("EXPERIMENT COMPLETE!")
    print("="*80)
    
    # Check if motifs helped
    ari_improvement = results_fused['ARI'] - results_domain_only['ARI']
    purity_improvement = results_fused['best_cluster_purity'] - results_domain_only['best_cluster_purity']
    
    print()
    if ari_improvement > 0:
        print("✅ MOTIF FEATURES IMPROVED CLUSTERING!")
        print(f"   ARI improved by {ari_improvement:.4f} ({ari_improvement/results_domain_only['ARI']*100:+.1f}%)")
    elif ari_improvement < -0.005:
        print("⚠️  MOTIF FEATURES DECREASED PERFORMANCE")
        print(f"   ARI decreased by {ari_improvement:.4f} ({ari_improvement/results_domain_only['ARI']*100:+.1f}%)")
    else:
        print("➖ MOTIF FEATURES HAD MINIMAL IMPACT")
        print(f"   ARI change: {ari_improvement:.4f}")
    
    print()
    print("Key metrics comparison:")
    print(f"  ARI:              {results_domain_only['ARI']:.4f} → {results_fused['ARI']:.4f}")
    print(f"  Purity:           {results_domain_only['Purity']:.4f} → {results_fused['Purity']:.4f}")
    print(f"  Best TK cluster:  {results_domain_only['best_cluster_purity']:.1%} → {results_fused['best_cluster_purity']:.1%}")
    print()
    
    return results_domain_only, results_fused


if __name__ == "__main__":
    main()

