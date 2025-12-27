#!/usr/bin/env python3
"""
Generate Panel (c): Layer-wise clustering performance for ESM-2.

Shows ARI performance for different ESM-2 layer configurations based on
actual experimental results.

Publication-quality formatting suitable for ML/computational biology journals.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
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


def plot_layer_performance(output_file):
    """Create panel (c): layer configuration clustering performance plot."""
    
    # Actual experimental data from documented results
    configurations = [
        'Layer 33\n(final only)',
        'All layers\n(1–33 mean)',
        'Layers 20–30\n(mean)',
        'Layers 20–33\n(mean)'
    ]
    ari_values = [0.268, 0.312, 0.353, 0.354]
    
    # Colors: highlight the best (layers 20-33) and baseline (layer 33)
    colors = ['#E63946', '#6C757D', '#2A9D8F', '#2A9D8F']
    edge_colors = ['#B71C1C', '#495057', '#1B7D71', '#1B7D71']
    
    # Create figure with white background
    fig, ax = plt.subplots(figsize=(6, 4.5), facecolor='white')
    ax.set_facecolor('white')
    
    # Create bar positions
    x_pos = np.arange(len(configurations))
    bar_width = 0.65
    
    # Plot bars
    bars = ax.bar(x_pos, ari_values, width=bar_width, color=colors,
                  edgecolor=edge_colors, linewidth=1.5, zorder=3)
    
    # Add value labels on top of bars
    for i, (bar, val) in enumerate(zip(bars, ari_values)):
        height = bar.get_height()
        # Add ARI value
        ax.text(bar.get_x() + bar.get_width()/2, height + 0.008,
                f'{val:.3f}', ha='center', va='bottom', fontsize=10,
                fontweight='bold', color='#1a1a2e')
    
    # Add improvement annotation for the best configuration
    # Arrow from layer 33 to layers 20-33
    ax.annotate('',
                xy=(3, 0.354), xytext=(0, 0.268),
                arrowprops=dict(arrowstyle='->', color='#2A9D8F',
                               connectionstyle='arc3,rad=-0.3', 
                               linewidth=2, linestyle='--'),
                zorder=2)
    
    # Add +32% label
    ax.text(1.5, 0.38, '+32%', fontsize=12, fontweight='bold',
            color='#2A9D8F', ha='center', va='bottom',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#E8F5F3',
                     edgecolor='#2A9D8F', linewidth=1))
    
    # Add baseline indicator
    ax.axhline(y=0.268, color='#E63946', linestyle=':', linewidth=1.2,
               alpha=0.7, zorder=1)
    ax.text(3.45, 0.268, 'baseline', fontsize=8, color='#E63946',
            va='center', ha='left', style='italic')
    
    # Configure axes
    ax.set_xlabel('Layer Configuration', fontsize=11, fontweight='medium')
    ax.set_ylabel('Adjusted Rand Index (ARI)', fontsize=11, fontweight='medium')
    
    # Set x-axis ticks
    ax.set_xticks(x_pos)
    ax.set_xticklabels(configurations, fontsize=9)
    
    # Set y-axis limits and ticks
    ax.set_ylim(0, 0.42)
    ax.set_yticks([0.0, 0.1, 0.2, 0.3, 0.4])
    
    # Thin axis lines, minimal styling
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    for spine in ['bottom', 'left']:
        ax.spines[spine].set_linewidth(0.8)
        ax.spines[spine].set_color('#333333')
    
    # Light grid (y-axis only, very subtle)
    ax.yaxis.grid(True, linestyle='-', alpha=0.2, linewidth=0.5, zorder=0)
    ax.xaxis.grid(False)
    
    # Add tick marks
    ax.tick_params(axis='both', which='major', length=4, width=0.8,
                   color='#333333', direction='out')
    ax.tick_params(axis='x', length=0)  # No ticks on x-axis for cleaner look
    
    # Add panel label "(c)" in top-left corner
    ax.text(-0.12, 1.05, '(c)', transform=ax.transAxes, fontsize=14,
            fontweight='bold', va='top', ha='left')
    
    # Add subtitle explaining the finding
    ax.set_title('ESM-2 Layer Selection Impact on Clustering',
                 fontsize=11, fontweight='medium', pad=10, color='#333333')
    
    # Tight layout
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
    data_df = pd.DataFrame({
        'configuration': ['Layer 33 (final)', 'All layers (1-33 mean)',
                          'Layers 20-30 (mean)', 'Layers 20-33 (mean)'],
        'layers': ['33', '1-33', '20-30', '20-33'],
        'ari': ari_values,
        'relative_improvement': ['baseline (0%)', '+16.4%', '+31.7%', '+32.1%']
    })
    data_df.to_csv(f"{output_file}_data.csv", index=False)
    
    print(f"✅ Layer performance figure saved:")
    print(f"   {output_file}.png (300 dpi)")
    print(f"   {output_file}.pdf (vector)")
    print(f"   {output_file}_data.csv (underlying data)")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description='Generate Panel (c): ESM-2 layer-wise clustering performance'
    )
    parser.add_argument('--output', default='figures_output/figure_panel_c_layer_performance',
                        help='Output file path (without extension)')
    args = parser.parse_args()
    
    plot_layer_performance(args.output)

