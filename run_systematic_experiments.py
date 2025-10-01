#!/usr/bin/env python3
"""
Systematic experiments to improve clustering accuracy.

Experiments:
1. Domain coverage: Different E-values + multiple HMMs
2. Layer probing: Single layer vs mean of mid-layers
3. Pooling strategy: Mean vs CLS token

This script orchestrates all experiments and generates a comprehensive comparison.
"""

import os
import subprocess
import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (
    adjusted_rand_score,
    normalized_mutual_info_score,
    homogeneity_completeness_v_measure,
    silhouette_score,
)
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import confusion_matrix


def run_command(cmd, desc):
    """Run a shell command and print status."""
    print(f"\n{'='*80}")
    print(f"Running: {desc}")
    print(f"Command: {cmd}")
    print(f"{'='*80}\n")
    
    result = subprocess.run(cmd, shell=True, capture_output=False, text=True)
    
    if result.returncode != 0:
        print(f"⚠️  Command failed with exit code {result.returncode}")
        return False
    
    print(f"\n✅ {desc} complete\n")
    return True


def run_clustering(emb_file, idx_file, labels_file, k=10, exclude_other=True):
    """
    Run K-Means clustering and compute metrics.
    
    Returns dict with metrics.
    """
    # Load embeddings
    X = np.load(emb_file)
    ids = pd.read_csv(idx_file)['uniprot_id'].astype(str).values
    
    # Load labels (from original kinases_domains.csv or similar)
    if os.path.exists(labels_file):
        labels_df = pd.read_csv(labels_file)
    else:
        # Try kinases_revised.csv
        labels_df = pd.read_csv('kinases_revised.csv')
    
    # Merge to get labels
    emb_df = pd.DataFrame({'uniprot_id': ids})
    merged = emb_df.merge(
        labels_df[['uniprot_id', 'kinome_group_major']],
        on='uniprot_id',
        how='left'
    )
    
    y_labels = merged['kinome_group_major'].values
    
    # Exclude "Other" if requested
    if exclude_other:
        mask = y_labels != "Other"
        X_clean = X[mask]
        y_clean = y_labels[mask]
    else:
        X_clean = X
        y_clean = y_labels
    
    if len(X_clean) < k:
        return None
    
    # Standardize
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_clean)
    
    # K-Means
    km = KMeans(n_clusters=k, random_state=42, n_init=50, max_iter=500)
    pred = km.fit_predict(X_scaled)
    
    # Metrics
    ari = adjusted_rand_score(y_clean, pred)
    nmi = normalized_mutual_info_score(y_clean, pred)
    h, c, v = homogeneity_completeness_v_measure(y_clean, pred)
    sil = silhouette_score(X_scaled, pred, sample_size=min(10000, len(X_scaled)), random_state=42)
    
    # Purity
    labels = sorted(np.unique(y_clean))
    label_to_idx = {lbl: i for i, lbl in enumerate(labels)}
    y_idx = np.array([label_to_idx[lbl] for lbl in y_clean])
    cm = confusion_matrix(y_idx, pred)
    purity = cm.max(axis=0).sum() / cm.sum()
    
    # Hungarian
    cost = cm.max() - cm
    row_ind, col_ind = linear_sum_assignment(cost)
    hungarian_acc = cm[row_ind, col_ind].sum() / cm.sum()
    
    return {
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
    }


def main():
    print("="*80)
    print("SYSTEMATIC EXPERIMENTS: DOMAIN COVERAGE + LAYER PROBING")
    print("="*80)
    print()
    
    results = []
    
    # Baseline: Current best (domain-only, E=0.001, PF00069, layer 33, mean pooling)
    print("\n" + "╔"+"="*78+"╗")
    print("║" + " "*25 + "BASELINE (CURRENT BEST)" + " "*28 + "║")
    print("╚"+"="*78+"╝")
    
    baseline_emb = "kinases_domains_embeddings/esm2_embeddings.npy"
    baseline_idx = "kinases_domains_embeddings/esm2_index.csv"
    baseline_labels = "kinases_domains.csv"
    
    if os.path.exists(baseline_emb):
        print("Using existing baseline embeddings...")
        baseline_metrics = run_clustering(baseline_emb, baseline_idx, baseline_labels)
        if baseline_metrics:
            baseline_metrics['experiment'] = 'Baseline: E=0.001, PF00069, Layer=33, Pool=mean'
            results.append(baseline_metrics)
            print(f"  ARI: {baseline_metrics['ARI']:.4f}")
            print(f"  NMI: {baseline_metrics['NMI']:.4f}")
    
    # Experiment 1: Improved domain coverage with relaxed E-values
    print("\n" + "╔"+"="*78+"╗")
    print("║" + " "*20 + "EXPERIMENT 1: DOMAIN COVERAGE" + " "*29 + "║")
    print("╚"+"="*78+"╝")
    
    for evalue in [0.01, 0.1]:
        exp_name = f"domains_e{evalue}"
        domain_file = f"kinases_domains_e{evalue}.csv"
        emb_dir = f"kinases_domains_e{evalue}_embeddings"
        
        # Extract domains
        if not os.path.exists(domain_file):
            cmd = (f"python extract_kinase_domains_v2.py "
                   f"--input kinases_revised.csv "
                   f"--output {domain_file} "
                   f"--hmms PF00069 PF07714 "
                   f"--evalue {evalue}")
            run_command(cmd, f"Extract domains (E={evalue})")
        
        # Generate embeddings (default: layer 33, mean pooling)
        if not os.path.exists(emb_dir):
            cmd = (f"python generate_esm2_embeddings_v2.py "
                   f"--input {domain_file} "
                   f"--output-dir {emb_dir} "
                   f"--layers 33 "
                   f"--pooling mean")
            run_command(cmd, f"Generate embeddings (E={evalue})")
        
        # Clustering
        emb_file = f"{emb_dir}/esm2_embeddings.npy"
        idx_file = f"{emb_dir}/esm2_index.csv"
        
        if os.path.exists(emb_file):
            print(f"\nClustering: E={evalue}...")
            metrics = run_clustering(emb_file, idx_file, domain_file)
            if metrics:
                metrics['experiment'] = f'E={evalue}, PF00069+PF07714, Layer=33, Pool=mean'
                results.append(metrics)
                print(f"  ARI: {metrics['ARI']:.4f} (Δ={metrics['ARI']-baseline_metrics['ARI']:+.4f})")
                print(f"  NMI: {metrics['NMI']:.4f} (Δ={metrics['NMI']-baseline_metrics['NMI']:+.4f})")
    
    # Experiment 2: Layer probing (use best domain file from above)
    print("\n" + "╔"+"="*78+"╗")
    print("║" + " "*22 + "EXPERIMENT 2: LAYER PROBING" + " "*29 + "║")
    print("╚"+"="*78+"╝")
    
    # Use E=0.01 domains for layer experiments
    domain_file = "kinases_domains_e0.01.csv"
    
    if not os.path.exists(domain_file):
        print(f"⚠️  {domain_file} not found, using kinases_domains.csv")
        domain_file = "kinases_domains.csv"
    
    for layer_spec in ['mid', '20-30']:
        exp_name = f"layers_{layer_spec.replace('-', '_')}"
        emb_dir = f"kinases_domains_e0.01_{exp_name}"
        
        # Generate embeddings with different layers
        if not os.path.exists(emb_dir):
            cmd = (f"python generate_esm2_embeddings_v2.py "
                   f"--input {domain_file} "
                   f"--output-dir {emb_dir} "
                   f"--layers {layer_spec} "
                   f"--pooling mean")
            run_command(cmd, f"Generate embeddings (layers={layer_spec})")
        
        # Clustering
        emb_file = f"{emb_dir}/esm2_embeddings.npy"
        idx_file = f"{emb_dir}/esm2_index.csv"
        
        if os.path.exists(emb_file):
            print(f"\nClustering: Layers={layer_spec}...")
            metrics = run_clustering(emb_file, idx_file, domain_file)
            if metrics:
                metrics['experiment'] = f'E=0.01, Layers={layer_spec}, Pool=mean'
                results.append(metrics)
                print(f"  ARI: {metrics['ARI']:.4f} (Δ={metrics['ARI']-baseline_metrics['ARI']:+.4f})")
                print(f"  NMI: {metrics['NMI']:.4f} (Δ={metrics['NMI']-baseline_metrics['NMI']:+.4f})")
    
    # Experiment 3: Pooling strategy (CLS vs mean)
    print("\n" + "╔"+"="*78+"╗")
    print("║" + " "*20 + "EXPERIMENT 3: POOLING STRATEGY" + " "*27 + "║")
    print("╚"+"="*78+"╝")
    
    emb_dir_cls = "kinases_domains_e0.01_cls"
    
    if not os.path.exists(emb_dir_cls):
        cmd = (f"python generate_esm2_embeddings_v2.py "
               f"--input {domain_file} "
               f"--output-dir {emb_dir_cls} "
               f"--layers 33 "
               f"--pooling cls")
        run_command(cmd, "Generate embeddings (CLS pooling)")
    
    # Clustering
    emb_file = f"{emb_dir_cls}/esm2_embeddings.npy"
    idx_file = f"{emb_dir_cls}/esm2_index.csv"
    
    if os.path.exists(emb_file):
        print(f"\nClustering: CLS pooling...")
        metrics = run_clustering(emb_file, idx_file, domain_file)
        if metrics:
            metrics['experiment'] = 'E=0.01, Layer=33, Pool=CLS'
            results.append(metrics)
            print(f"  ARI: {metrics['ARI']:.4f} (Δ={metrics['ARI']-baseline_metrics['ARI']:+.4f})")
            print(f"  NMI: {metrics['NMI']:.4f} (Δ={metrics['NMI']-baseline_metrics['NMI']:+.4f})")
    
    # Generate comparison report
    print("\n" + "="*80)
    print("GENERATING COMPARISON REPORT")
    print("="*80)
    
    if results:
        results_df = pd.DataFrame(results)
        
        # Sort by ARI (best first)
        results_df = results_df.sort_values('ARI', ascending=False)
        
        # Save detailed results
        results_df.to_csv('clustering/systematic_experiments_results.csv', index=False)
        print("✅ Saved: clustering/systematic_experiments_results.csv")
        
        # Print summary
        print("\n" + "="*80)
        print("RESULTS SUMMARY (sorted by ARI)")
        print("="*80)
        print()
        
        # Print table
        print(f"{'Experiment':<50} {'ARI':>8} {'NMI':>8} {'Purity':>8} {'Hungarian':>10}")
        print("-"*80)
        
        for _, row in results_df.iterrows():
            print(f"{row['experiment']:<50} {row['ARI']:>8.4f} {row['NMI']:>8.4f} {row['Purity']:>8.4f} {row['Hungarian']:>10.4f}")
        
        print()
        
        # Best experiment
        best = results_df.iloc[0]
        print("🏆 BEST CONFIGURATION:")
        print(f"   {best['experiment']}")
        print(f"   ARI: {best['ARI']:.4f}")
        print(f"   NMI: {best['NMI']:.4f}")
        print(f"   Purity: {best['Purity']:.4f}")
        print(f"   Hungarian: {best['Hungarian']:.4f}")
        print()
        
        # Improvement over baseline
        if len(results_df) > 1:
            baseline_ari = baseline_metrics['ARI']
            best_ari = best['ARI']
            improvement = best_ari - baseline_ari
            pct_improvement = (improvement / baseline_ari) * 100
            
            print(f"📈 IMPROVEMENT OVER BASELINE:")
            print(f"   ARI: {baseline_ari:.4f} → {best_ari:.4f} ({improvement:+.4f}, {pct_improvement:+.1f}%)")
            print()
    
    print("="*80)
    print("✅ SYSTEMATIC EXPERIMENTS COMPLETE!")
    print("="*80)


if __name__ == "__main__":
    main()

