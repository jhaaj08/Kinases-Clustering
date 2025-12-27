#!/usr/bin/env python3
"""
Generate Figure 6: Calibration and Retrieval Quality

Three-panel figure showing:
(a) Calibration curves (Final vs Mid layer)
(b) Expected Calibration Error (ECE) bar plot
(c) Nearest-neighbor retrieval precision@k

Purpose: Shows reliability + interpretability, not just accuracy.
"""

import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
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


def load_calibration_data():
    """Load calibration comparison results."""
    with open('calibration_comparison_results/calibration_comparison.json', 'r') as f:
        return json.load(f)


def load_retrieval_data():
    """Load retrieval results."""
    with open('exemplar_retrieval_results/retrieval_summary.json', 'r') as f:
        retrieval = json.load(f)
    
    per_class = pd.read_csv('exemplar_retrieval_results/per_class_retrieval.csv')
    similarity_cal = pd.read_csv('exemplar_retrieval_results/similarity_calibration.csv')
    
    return retrieval, per_class, similarity_cal


def plot_figure6(output_file):
    """Create Figure 6: Calibration and Retrieval Quality."""
    
    # Check data exists
    if not Path('calibration_comparison_results/calibration_comparison.json').exists():
        print("Error: Run scripts/run_calibration_comparison.py first")
        return False
    
    if not Path('exemplar_retrieval_results/retrieval_summary.json').exists():
        print("Error: Retrieval results not found")
        return False
    
    # Load data
    cal_data = load_calibration_data()
    retrieval, per_class, similarity_cal = load_retrieval_data()
    
    # Colors
    final_color = '#E63946'
    mid_color = '#2A9D8F'
    
    # Create figure
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), facecolor='white')
    
    for ax in axes:
        ax.set_facecolor('white')
    
    # Panel (a): Calibration Curves
    ax1 = axes[0]
    
    # Perfect calibration line
    ax1.plot([0, 1], [0, 1], 'k--', linewidth=1.5, label='Perfect calibration', alpha=0.7)
    
    # Extract calibration curves from bin statistics
    for config_name, config_data in cal_data.items():
        bin_stats = config_data['bin_stats']
        
        if len(bin_stats) > 0:
            confs = [b['confidence'] for b in bin_stats]
            accs = [b['accuracy'] for b in bin_stats]
            sizes = [b['size'] for b in bin_stats]
            
            color = final_color if 'Final' in config_name else mid_color
            label = config_name.split(' (')[0]  # Short name
            
            ax1.plot(confs, accs, 'o-', color=color, linewidth=2, 
                    markersize=6, label=f'{label} (ECE={config_data["ece"]:.3f})')
    
    ax1.set_xlabel('Mean Predicted Probability', fontsize=11, fontweight='medium')
    ax1.set_ylabel('Fraction Correct', fontsize=11, fontweight='medium')
    ax1.set_title('Reliability Diagram', fontsize=11, fontweight='medium')
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1)
    ax1.legend(loc='lower right', frameon=True, framealpha=0.9)
    ax1.set_aspect('equal')
    
    ax1.text(-0.15, 1.08, '(a)', transform=ax1.transAxes, fontsize=14,
             fontweight='bold', va='top', ha='left')
    
    for spine in ['top', 'right']:
        ax1.spines[spine].set_visible(False)
    ax1.grid(True, linestyle='-', alpha=0.2, linewidth=0.5)
    
    # Panel (b): ECE Bar Plot
    ax2 = axes[1]
    
    configs = ['Final\n(Layer 33)', 'Mid\n(Layers 19-33)']
    eces = [
        cal_data['Final (Layer 33)']['ece'],
        cal_data['Mid (Layers 19-33)']['ece']
    ]
    colors = [final_color, mid_color]
    
    bars = ax2.bar(configs, eces, color=colors, edgecolor='#333333', linewidth=1.2)
    
    ax2.set_ylabel('Expected Calibration Error (ECE)', fontsize=11, fontweight='medium')
    ax2.set_title('Calibration Quality\n(Lower is Better)', fontsize=11, fontweight='medium')
    ax2.set_ylim(0, 0.25)
    
    # Add value labels
    for bar, ece in zip(bars, eces):
        ax2.text(bar.get_x() + bar.get_width()/2, ece + 0.008,
                f'{ece:.3f}', ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    # Highlight better (lower ECE)
    best_idx = np.argmin(eces)
    bars[best_idx].set_edgecolor('#FFD700')
    bars[best_idx].set_linewidth(3)
    
    # Add improvement annotation
    improvement = (eces[0] - eces[1]) / eces[0] * 100
    ax2.annotate(f'{improvement:.0f}% better\ncalibration',
                xy=(1, eces[1]), xytext=(1.3, 0.10),
                fontsize=9, ha='left', va='center',
                arrowprops=dict(arrowstyle='->', color=mid_color, linewidth=1.5),
                color=mid_color, fontweight='bold')
    
    ax2.text(-0.15, 1.08, '(b)', transform=ax2.transAxes, fontsize=14,
             fontweight='bold', va='top', ha='left')
    
    for spine in ['top', 'right']:
        ax2.spines[spine].set_visible(False)
    ax2.yaxis.grid(True, linestyle='-', alpha=0.2, linewidth=0.5)
    
    # Panel (c): Retrieval Precision@k
    ax3 = axes[2]
    
    k_values = [1, 3, 5, 10]
    precision_at_k = [
        retrieval['retrieval_performance']['top1_hit_rate'],
        retrieval['retrieval_performance']['top3_hit_rate'],
        retrieval['retrieval_performance']['top5_hit_rate'],
        retrieval['retrieval_performance']['top10_hit_rate']
    ]
    
    ax3.plot(k_values, precision_at_k, 'o-', color='#3C5488', 
             linewidth=2.5, markersize=10, markeredgecolor='white', 
             markeredgewidth=2)
    
    # Add value labels
    for k, p in zip(k_values, precision_at_k):
        offset = 0.03 if k != 10 else -0.03
        va = 'bottom' if k != 10 else 'top'
        ax3.text(k, p + offset, f'{p:.1%}', ha='center', va=va, 
                fontsize=9, fontweight='bold')
    
    # Add MRR annotation
    mrr = retrieval['retrieval_performance']['mrr']
    ax3.axhline(y=mrr, color='#F39B7F', linestyle='--', linewidth=1.5, alpha=0.8)
    ax3.text(8, mrr + 0.02, f'MRR = {mrr:.3f}', fontsize=9, 
             color='#F39B7F', fontweight='bold')
    
    ax3.set_xlabel('k (Number of Retrieved Neighbors)', fontsize=11, fontweight='medium')
    ax3.set_ylabel('Precision@k (Hit Rate)', fontsize=11, fontweight='medium')
    ax3.set_title('Nearest-Neighbor Retrieval', fontsize=11, fontweight='medium')
    ax3.set_xlim(0, 12)
    ax3.set_ylim(0.6, 1.0)
    ax3.set_xticks(k_values)
    
    ax3.text(-0.15, 1.08, '(c)', transform=ax3.transAxes, fontsize=14,
             fontweight='bold', va='top', ha='left')
    
    for spine in ['top', 'right']:
        ax3.spines[spine].set_visible(False)
    ax3.yaxis.grid(True, linestyle='-', alpha=0.2, linewidth=0.5)
    
    plt.tight_layout()
    
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
    summary_data = pd.DataFrame({
        'Configuration': ['Final (Layer 33)', 'Mid (Layers 19-33)'],
        'ECE': eces
    })
    
    retrieval_data = pd.DataFrame({
        'k': k_values,
        'precision_at_k': precision_at_k
    })
    
    summary_data.to_csv(f"{output_file}_ece_data.csv", index=False)
    retrieval_data.to_csv(f"{output_file}_retrieval_data.csv", index=False)
    
    print(f"\n✅ Figure 6 saved:")
    print(f"   {output_file}.png (300 dpi)")
    print(f"   {output_file}.pdf (vector)")
    print(f"   {output_file}_ece_data.csv")
    print(f"   {output_file}_retrieval_data.csv")
    
    print("\n📊 Summary:")
    print(f"   Final layer ECE: {eces[0]:.3f}")
    print(f"   Mid layers ECE:  {eces[1]:.3f} ({improvement:.0f}% better)")
    print(f"   Retrieval MRR:   {mrr:.3f}")
    print(f"   Precision@1:     {precision_at_k[0]:.1%}")
    print(f"   Precision@5:     {precision_at_k[2]:.1%}")
    
    return True


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description='Generate Figure 6: Calibration and Retrieval Quality'
    )
    parser.add_argument('--output', default='figures_output/figure6_calibration_retrieval',
                        help='Output file path (without extension)')
    args = parser.parse_args()
    
    plot_figure6(args.output)


