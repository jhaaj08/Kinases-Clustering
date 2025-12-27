#!/usr/bin/env python3
"""
Run supervised classification across homology thresholds for different layer configs.

Tests Final layer vs Mid-layer at 70%, 50%, and 40% identity thresholds
to demonstrate generalization under homology constraints.

Output: Data for Figure 5 - Generalization Under Homology Constraints
"""

import os
import json
import numpy as np
import pandas as pd
from pathlib import Path

from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score, f1_score

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
    train_ids = splits['train_ids']
    test_ids = splits['test_ids']
    
    train_mask = data['uniprot_id'].isin(train_ids)
    test_mask = data['uniprot_id'].isin(test_ids)
    
    X_train = embeddings[train_mask]
    X_test = embeddings[test_mask]
    y_train = data.loc[train_mask, 'kinome_group_major'].values
    y_test = data.loc[test_mask, 'kinome_group_major'].values
    
    # Filter out Other and NaN
    valid_train = pd.notna(y_train) & (y_train != 'Other')
    valid_test = pd.notna(y_test) & (y_test != 'Other')
    
    return (X_train[valid_train], X_test[valid_test], 
            y_train[valid_train], y_test[valid_test])


def train_evaluate(X_train, X_test, y_train, y_test):
    """Train and evaluate classifier."""
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
    
    y_pred = clf.predict(X_test_s)
    
    return {
        'accuracy': accuracy_score(y_test, y_pred),
        'macro_f1': f1_score(y_test, y_pred, average='macro'),
        'train_size': len(X_train),
        'test_size': len(X_test)
    }


def main():
    print("="*80)
    print("HOMOLOGY GENERALIZATION EXPERIMENTS")
    print("="*80)
    print()
    
    # Configurations
    layer_configs = {
        'Final (Layer 33)': 'kinases_domains_e0.01_embeddings',
        'Mid (Layers 19-33)': 'kinases_domains_e0.01_layers_mid',
    }
    
    identity_thresholds = [70, 50, 40]
    
    labels_csv = 'data/processed/kinases_domains.csv'
    if not Path(labels_csv).exists():
        labels_csv = 'data/processed/kinases_domains_e0.01.csv'
    
    output_dir = 'supervised_results_homology'
    os.makedirs(output_dir, exist_ok=True)
    
    # Results storage
    all_results = []
    
    for identity in identity_thresholds:
        splits_file = f'data/splits_{identity}.json'
        
        print(f"\n{'='*60}")
        print(f"IDENTITY THRESHOLD: {identity}%")
        print(f"{'='*60}")
        
        if not Path(splits_file).exists():
            print(f"  ⚠️ Splits file not found: {splits_file}")
            continue
        
        for config_name, emb_dir in layer_configs.items():
            print(f"\n  {config_name}:")
            
            if not Path(emb_dir).exists():
                print(f"    ⚠️ Embeddings not found: {emb_dir}")
                continue
            
            X_train, X_test, y_train, y_test = load_data(
                emb_dir, labels_csv, splits_file
            )
            
            results = train_evaluate(X_train, X_test, y_train, y_test)
            
            print(f"    Accuracy: {results['accuracy']:.4f}")
            print(f"    Macro-F1: {results['macro_f1']:.4f}")
            print(f"    Train/Test: {results['train_size']}/{results['test_size']}")
            
            all_results.append({
                'Identity_Threshold': identity,
                'Layer_Config': config_name,
                'Accuracy': results['accuracy'],
                'Macro_F1': results['macro_f1'],
                'Train_Size': results['train_size'],
                'Test_Size': results['test_size']
            })
    
    # Create results DataFrame
    results_df = pd.DataFrame(all_results)
    
    # Pivot for easier analysis
    print("\n" + "="*80)
    print("SUMMARY: ACCURACY BY IDENTITY × LAYER CONFIG")
    print("="*80)
    
    acc_pivot = results_df.pivot(
        index='Identity_Threshold', 
        columns='Layer_Config', 
        values='Accuracy'
    )
    print(acc_pivot.to_string())
    
    print("\n" + "="*80)
    print("SUMMARY: MACRO-F1 BY IDENTITY × LAYER CONFIG")
    print("="*80)
    
    f1_pivot = results_df.pivot(
        index='Identity_Threshold', 
        columns='Layer_Config', 
        values='Macro_F1'
    )
    print(f1_pivot.to_string())
    
    # Calculate performance gap (Final - Mid)
    print("\n" + "="*80)
    print("PERFORMANCE GAP (Final - Mid)")
    print("="*80)
    
    gap_data = []
    for identity in identity_thresholds:
        final_row = results_df[(results_df['Identity_Threshold'] == identity) & 
                               (results_df['Layer_Config'] == 'Final (Layer 33)')]
        mid_row = results_df[(results_df['Identity_Threshold'] == identity) & 
                             (results_df['Layer_Config'] == 'Mid (Layers 19-33)')]
        
        if len(final_row) > 0 and len(mid_row) > 0:
            acc_gap = final_row['Accuracy'].values[0] - mid_row['Accuracy'].values[0]
            f1_gap = final_row['Macro_F1'].values[0] - mid_row['Macro_F1'].values[0]
            
            gap_data.append({
                'Identity_Threshold': identity,
                'Accuracy_Gap': acc_gap,
                'Macro_F1_Gap': f1_gap
            })
            
            print(f"  {identity}% identity:")
            print(f"    Accuracy gap: {acc_gap:+.4f} ({acc_gap*100:+.2f}%)")
            print(f"    Macro-F1 gap: {f1_gap:+.4f}")
    
    gap_df = pd.DataFrame(gap_data)
    
    # Save results
    results_df.to_csv(f"{output_dir}/homology_results.csv", index=False)
    acc_pivot.to_csv(f"{output_dir}/accuracy_by_identity_layer.csv")
    f1_pivot.to_csv(f"{output_dir}/macro_f1_by_identity_layer.csv")
    gap_df.to_csv(f"{output_dir}/performance_gap.csv", index=False)
    
    print(f"\n✅ Results saved to: {output_dir}/")
    print(f"   - homology_results.csv")
    print(f"   - accuracy_by_identity_layer.csv")
    print(f"   - macro_f1_by_identity_layer.csv")
    print(f"   - performance_gap.csv")
    
    print("\n" + "="*80)
    print("✅ HOMOLOGY GENERALIZATION EXPERIMENTS COMPLETE!")
    print("="*80)


if __name__ == "__main__":
    main()


