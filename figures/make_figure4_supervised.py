#!/usr/bin/env python3
"""
Generate Figure 4: Supervised Classification Performance

Three-panel figure showing:
(a) Accuracy vs Layer Configuration
(b) Macro-F1 vs Layer Configuration  
(c) Confusion matrices: Final layer vs Mid-layer (side-by-side)

Uses actual experimental results from layer comparison.
"""

import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

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
    'legend.fontsize': 9,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1,
})


def load_results(results_dir='supervised_results_layer_comparison'):
    """Load experimental results."""
    with open(f"{results_dir}/full_results.json", 'r') as f:
        results = json.load(f)
    
    summary_df = pd.read_csv(f"{results_dir}/layer_comparison_summary.csv")
    per_class_df = pd.read_csv(f"{results_dir}/per_class_metrics.csv")
    
    # Load confusion matrices
    confusion_matrices = {}
    for config in results.keys():
        safe_name = config.replace(' ', '_').replace('(', '').replace(')', '').replace('-', '_')
        cm_file = f"{results_dir}/confusion_matrix_{safe_name}.csv"
        if Path(cm_file).exists():
            confusion_matrices[config] = pd.read_csv(cm_file, index_col=0)
    
    return results, summary_df, per_class_df, confusion_matrices


def plot_figure4(output_file):
    """Create Figure 4: Supervised Classification Performance."""
    
    # Load results
    results_dir = 'supervised_results_layer_comparison'
    
    if not Path(results_dir).exists():
        print(f"❌ Results directory not found: {results_dir}")
        print("   Run: python scripts/run_layer_supervised_comparison.py first")
        return False
    
    results, summary_df, per_class_df, confusion_matrices = load_results(results_dir)
    
    # Create figure with 3 panels
    fig = plt.figure(figsize=(14, 10), facecolor='white')
    
    # Layout: 2 bar charts on top, 2 confusion matrices on bottom
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.2], hspace=0.35, wspace=0.3)
    
    ax1 = fig.add_subplot(gs[0, 0])  # Accuracy
    ax2 = fig.add_subplot(gs[0, 1])  # Macro-F1
    ax3 = fig.add_subplot(gs[1, 0])  # Confusion matrix - Final
    ax4 = fig.add_subplot(gs[1, 1])  # Confusion matrix - Mid
    
    for ax in [ax1, ax2, ax3, ax4]:
        ax.set_facecolor('white')
    
    # Data for bar charts
    configs = ['Final\n(Layer 33)', 'Mid\n(Layers 20-30)', 'Extended Mid\n(Layers 19-33)']
    config_keys = ['Final (Layer 33)', 'Mid (Layers 20-30)', 'Extended Mid (Layers 19-33)']
    
    accuracies = [results[k]['accuracy'] for k in config_keys]
    macro_f1s = [results[k]['macro_f1'] for k in config_keys]
    
    # Colors
    colors = ['#E63946', '#6C757D', '#2A9D8F']
    
    # Panel (a): Accuracy
    bars1 = ax1.bar(configs, accuracies, color=colors, edgecolor='#333333', linewidth=1.2)
    ax1.set_ylabel('Accuracy', fontsize=11, fontweight='medium')
    ax1.set_ylim(0.65, 0.80)
    ax1.set_title('Test Accuracy by Layer Configuration', fontsize=11, fontweight='medium')
    
    # Add value labels
    for bar, val in zip(bars1, accuracies):
        ax1.text(bar.get_x() + bar.get_width()/2, val + 0.005,
                f'{val:.1%}', ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    # Highlight best
    best_idx = np.argmax(accuracies)
    bars1[best_idx].set_edgecolor('#FFD700')
    bars1[best_idx].set_linewidth(3)
    
    ax1.text(-0.12, 1.08, '(a)', transform=ax1.transAxes, fontsize=14,
             fontweight='bold', va='top', ha='left')
    
    for spine in ['top', 'right']:
        ax1.spines[spine].set_visible(False)
    ax1.yaxis.grid(True, linestyle='-', alpha=0.2, linewidth=0.5)
    
    # Panel (b): Macro-F1
    bars2 = ax2.bar(configs, macro_f1s, color=colors, edgecolor='#333333', linewidth=1.2)
    ax2.set_ylabel('Macro-F1', fontsize=11, fontweight='medium')
    ax2.set_ylim(0.55, 0.70)
    ax2.set_title('Macro-F1 by Layer Configuration', fontsize=11, fontweight='medium')
    
    for bar, val in zip(bars2, macro_f1s):
        ax2.text(bar.get_x() + bar.get_width()/2, val + 0.005,
                f'{val:.3f}', ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    best_idx = np.argmax(macro_f1s)
    bars2[best_idx].set_edgecolor('#FFD700')
    bars2[best_idx].set_linewidth(3)
    
    ax2.text(-0.12, 1.08, '(b)', transform=ax2.transAxes, fontsize=14,
             fontweight='bold', va='top', ha='left')
    
    for spine in ['top', 'right']:
        ax2.spines[spine].set_visible(False)
    ax2.yaxis.grid(True, linestyle='-', alpha=0.2, linewidth=0.5)
    
    # Panel (c): Confusion matrix - Final Layer
    if 'Final (Layer 33)' in confusion_matrices:
        cm_final = confusion_matrices['Final (Layer 33)']
        
        # Normalize for display
        cm_norm = cm_final.div(cm_final.sum(axis=1), axis=0)
        
        sns.heatmap(cm_norm, annot=cm_final.values, fmt='d', cmap='Blues',
                   ax=ax3, cbar=False, linewidths=0.5, linecolor='white',
                   annot_kws={'size': 8})
        ax3.set_title('Final Layer (Layer 33)\nAccuracy: 75.1%', fontsize=10, fontweight='medium')
        ax3.set_xlabel('Predicted', fontsize=10)
        ax3.set_ylabel('True', fontsize=10)
        ax3.tick_params(axis='both', labelsize=8)
        
        ax3.text(-0.15, 1.12, '(c)', transform=ax3.transAxes, fontsize=14,
                fontweight='bold', va='top', ha='left')
    
    # Panel (d): Confusion matrix - Extended Mid
    if 'Extended Mid (Layers 19-33)' in confusion_matrices:
        cm_mid = confusion_matrices['Extended Mid (Layers 19-33)']
        
        cm_norm = cm_mid.div(cm_mid.sum(axis=1), axis=0)
        
        sns.heatmap(cm_norm, annot=cm_mid.values, fmt='d', cmap='Greens',
                   ax=ax4, cbar=False, linewidths=0.5, linecolor='white',
                   annot_kws={'size': 8})
        ax4.set_title('Extended Mid (Layers 19-33)\nAccuracy: 73.8%', fontsize=10, fontweight='medium')
        ax4.set_xlabel('Predicted', fontsize=10)
        ax4.set_ylabel('True', fontsize=10)
        ax4.tick_params(axis='both', labelsize=8)
        
        ax4.text(-0.15, 1.12, '(d)', transform=ax4.transAxes, fontsize=14,
                fontweight='bold', va='top', ha='left')
    
    # Add note about unsupervised vs supervised finding
    note_text = ("Note: For supervised classification, final layer performs best.\n"
                 "This contrasts with unsupervised clustering where mid-layers excel (+32% ARI).")
    fig.text(0.5, 0.01, note_text, ha='center', va='bottom', fontsize=9,
             style='italic', color='#555555',
             bbox=dict(boxstyle='round,pad=0.4', facecolor='#F5F5F5', 
                      edgecolor='#CCCCCC', linewidth=0.5))
    
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.08)
    
    # Ensure output directory exists
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Save figures
    plt.savefig(f"{output_file}.png", dpi=300, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.savefig(f"{output_file}.pdf", bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close()
    
    # Save underlying data
    summary_df.to_csv(f"{output_file}_summary.csv", index=False)
    
    print(f"\n✅ Figure 4 saved:")
    print(f"   {output_file}.png (300 dpi)")
    print(f"   {output_file}.pdf (vector)")
    print(f"   {output_file}_summary.csv (underlying data)")
    
    return True


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description='Generate Figure 4: Supervised Classification Performance'
    )
    parser.add_argument('--output', default='figures_output/figure4_supervised_performance',
                        help='Output file path (without extension)')
    args = parser.parse_args()
    
    plot_figure4(args.output)


