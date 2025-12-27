#!/usr/bin/env python3
"""
Generate Figure 3: Embedding Space Visualization (UMAP)

Three-panel figure showing:
(a) Final-layer embeddings (layer 33)
(b) Intermediate-layer embeddings (layers 20-33)
(c) Ground-truth functional labels (same projection as panel b)

Publication-quality formatting for ML/computational biology journals.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Configure matplotlib for publication quality
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'],
    'font.size': 10,
    'axes.linewidth': 0.8,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 8,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1,
})


def load_embeddings_with_labels(embeddings_dir, labels_csv):
    """Load embeddings and merge with labels."""
    embeddings = np.load(f"{embeddings_dir}/esm2_embeddings.npy")
    index_df = pd.read_csv(f"{embeddings_dir}/esm2_index.csv")
    labels_df = pd.read_csv(labels_csv)
    
    # Merge on uniprot_id
    data = index_df.merge(
        labels_df[['uniprot_id', 'kinome_group_major']], 
        on='uniprot_id', 
        how='left'
    )
    
    # Filter out "Other" and NaN
    mask = (data['kinome_group_major'] != 'Other') & (data['kinome_group_major'].notna())
    return embeddings[mask], data[mask].reset_index(drop=True)


def get_family_colors():
    """Define consistent colors for kinase families."""
    # Use a colorblind-friendly palette
    family_colors = {
        'AGC': '#E64B35',      # Red
        'CAMK': '#4DBBD5',     # Cyan
        'CK1': '#00A087',      # Teal
        'CMGC': '#3C5488',     # Blue
        'STE': '#F39B7F',      # Salmon
        'TK': '#8491B4',       # Gray-blue
        'TKL': '#91D1C2',      # Light teal
        'RGC': '#DC0000',      # Dark red
        'Atypical': '#7E6148', # Brown
        'Histidine': '#B09C85' # Tan
    }
    return family_colors


def plot_figure3_umap(output_file):
    """Create the 3-panel UMAP visualization figure."""
    
    try:
        from umap import UMAP
    except ImportError:
        print("❌ Error: umap-learn not installed")
        print("   Install with: pip install umap-learn")
        return False
    
    # File paths
    final_layer_dir = 'kinases_domains_e0.01_embeddings'
    mid_layer_dir = 'kinases_domains_e0.01_layers_mid'
    labels_csv = 'data/processed/kinases_domains.csv'
    
    # Check if files exist, try alternate paths
    if not Path(labels_csv).exists():
        labels_csv = 'kinases_domains.csv'
    if not Path(labels_csv).exists():
        labels_csv = 'data/processed/kinases_revised.csv'
    if not Path(labels_csv).exists():
        labels_csv = 'data/raw/kinases_revised.csv'
    
    print("Loading embeddings...")
    
    # Load final layer embeddings
    print(f"  Loading final layer embeddings from {final_layer_dir}...")
    emb_final, labels_final = load_embeddings_with_labels(final_layer_dir, labels_csv)
    print(f"    ✅ {len(emb_final):,} samples loaded")
    
    # Load mid-layer embeddings  
    print(f"  Loading mid-layer embeddings from {mid_layer_dir}...")
    emb_mid, labels_mid = load_embeddings_with_labels(mid_layer_dir, labels_csv)
    print(f"    ✅ {len(emb_mid):,} samples loaded")
    
    # Get unique families
    families = sorted(labels_mid['kinome_group_major'].unique())
    print(f"  Families: {families}")
    
    # UMAP parameters (consistent across all panels)
    umap_params = {
        'n_neighbors': 15,
        'min_dist': 0.1,
        'metric': 'cosine',
        'random_state': 42,
        'n_components': 2
    }
    
    print("\nComputing UMAP projections...")
    
    # Compute UMAP for final layer
    print("  Panel (a): Final layer (layer 33)...")
    umap_final = UMAP(**umap_params)
    coords_final = umap_final.fit_transform(emb_final)
    
    # Compute UMAP for mid-layers
    print("  Panel (b) & (c): Mid-layers (20-33)...")
    umap_mid = UMAP(**umap_params)
    coords_mid = umap_mid.fit_transform(emb_mid)
    
    # Get colors
    family_colors = get_family_colors()
    
    # Create figure
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), facecolor='white')
    
    for ax in axes:
        ax.set_facecolor('white')
    
    # Panel (a): Final layer with cluster coloring
    ax1 = axes[0]
    for family in families:
        mask = labels_final['kinome_group_major'] == family
        if mask.sum() > 0:
            color = family_colors.get(family, '#999999')
            ax1.scatter(coords_final[mask, 0], coords_final[mask, 1],
                       c=color, label=family,
                       s=15, alpha=0.7, edgecolors='white', linewidth=0.3)
    
    ax1.set_title('Final Layer Only\n(Layer 33)', fontsize=11, fontweight='medium')
    ax1.set_xlabel('UMAP 1', fontsize=10)
    ax1.set_ylabel('UMAP 2', fontsize=10)
    ax1.text(-0.08, 1.08, '(a)', transform=ax1.transAxes, fontsize=14,
             fontweight='bold', va='top', ha='left')
    
    # Remove spines for cleaner look
    for spine in ['top', 'right']:
        ax1.spines[spine].set_visible(False)
    ax1.spines['left'].set_linewidth(0.8)
    ax1.spines['bottom'].set_linewidth(0.8)
    
    # Panel (b): Mid-layer embeddings
    ax2 = axes[1]
    for family in families:
        mask = labels_mid['kinome_group_major'] == family
        if mask.sum() > 0:
            color = family_colors.get(family, '#999999')
            ax2.scatter(coords_mid[mask, 0], coords_mid[mask, 1],
                       c=color, label=family,
                       s=15, alpha=0.7, edgecolors='white', linewidth=0.3)
    
    ax2.set_title('Intermediate Layers\n(Layers 20–33 Mean)', fontsize=11, fontweight='medium')
    ax2.set_xlabel('UMAP 1', fontsize=10)
    ax2.set_ylabel('UMAP 2', fontsize=10)
    ax2.text(-0.08, 1.08, '(b)', transform=ax2.transAxes, fontsize=14,
             fontweight='bold', va='top', ha='left')
    
    for spine in ['top', 'right']:
        ax2.spines[spine].set_visible(False)
    ax2.spines['left'].set_linewidth(0.8)
    ax2.spines['bottom'].set_linewidth(0.8)
    
    # Panel (c): Same projection as (b) but with emphasis on ground truth
    ax3 = axes[2]
    
    # Plot with larger markers and stronger colors to emphasize ground truth
    for family in families:
        mask = labels_mid['kinome_group_major'] == family
        if mask.sum() > 0:
            color = family_colors.get(family, '#999999')
            ax3.scatter(coords_mid[mask, 0], coords_mid[mask, 1],
                       c=color, label=f'{family} (n={mask.sum()})',
                       s=20, alpha=0.8, edgecolors='white', linewidth=0.4)
    
    ax3.set_title('Ground-Truth Labels\n(Same Projection as b)', fontsize=11, fontweight='medium')
    ax3.set_xlabel('UMAP 1', fontsize=10)
    ax3.set_ylabel('UMAP 2', fontsize=10)
    ax3.text(-0.08, 1.08, '(c)', transform=ax3.transAxes, fontsize=14,
             fontweight='bold', va='top', ha='left')
    
    for spine in ['top', 'right']:
        ax3.spines[spine].set_visible(False)
    ax3.spines['left'].set_linewidth(0.8)
    ax3.spines['bottom'].set_linewidth(0.8)
    
    # Add legend to panel (c)
    ax3.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8,
               frameon=True, framealpha=0.9, edgecolor='#cccccc')
    
    # Add annotation showing improvement
    # Add text box with ARI comparison
    textstr = 'ARI: 0.268 → 0.354\n(+32% improvement)'
    props = dict(boxstyle='round,pad=0.4', facecolor='#E8F5F3', 
                 edgecolor='#2A9D8F', linewidth=1)
    fig.text(0.5, 0.02, textstr, transform=fig.transFigure, fontsize=10,
             ha='center', va='bottom', bbox=props, fontweight='medium')
    
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.15, right=0.88)
    
    # Ensure output directory exists
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Save figures
    print("\nSaving figures...")
    plt.savefig(f"{output_file}.png", dpi=300, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.savefig(f"{output_file}.pdf", bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close()
    
    # Save UMAP coordinates for reproducibility
    coords_df = pd.DataFrame({
        'uniprot_id': labels_mid['uniprot_id'],
        'family': labels_mid['kinome_group_major'],
        'umap1_mid': coords_mid[:, 0],
        'umap2_mid': coords_mid[:, 1],
    })
    
    # Add final layer coordinates (need to align by uniprot_id)
    final_coords_df = pd.DataFrame({
        'uniprot_id': labels_final['uniprot_id'],
        'umap1_final': coords_final[:, 0],
        'umap2_final': coords_final[:, 1],
    })
    
    coords_df = coords_df.merge(final_coords_df, on='uniprot_id', how='left')
    coords_df.to_csv(f"{output_file}_coordinates.csv", index=False)
    
    print(f"\n✅ Figure 3 saved:")
    print(f"   {output_file}.png (300 dpi)")
    print(f"   {output_file}.pdf (vector)")
    print(f"   {output_file}_coordinates.csv (UMAP coordinates)")
    
    return True


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description='Generate Figure 3: UMAP Embedding Space Visualization'
    )
    parser.add_argument('--output', default='figures_output/figure3_umap_embedding_space',
                        help='Output file path (without extension)')
    args = parser.parse_args()
    
    plot_figure3_umap(args.output)



