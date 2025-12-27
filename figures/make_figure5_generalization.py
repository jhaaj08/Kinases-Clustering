#!/usr/bin/env python3
"""
Generate Figure 5: Generalization Under Homology Constraints

Three-panel figure showing:
(a) Accuracy vs Homology Threshold (70% / 50% / 40%)
(b) Macro-F1 vs Homology Threshold
(c) Relative performance gap (Final − Mid) across thresholds

Purpose: Demonstrates robust abstraction, not memorization.
"""

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


def load_results(results_dir='supervised_results_homology'):
    """Load homology generalization results."""
    results_df = pd.read_csv(f"{results_dir}/homology_results.csv")
    gap_df = pd.read_csv(f"{results_dir}/performance_gap.csv")
    return results_df, gap_df


def plot_figure5(output_file):
    """Create Figure 5: Generalization Under Homology Constraints."""
    
    results_dir = 'supervised_results_homology'
    
    if not Path(results_dir).exists():
        print(f"❌ Results directory not found: {results_dir}")
        print("   Run: python scripts/run_homology_generalization.py first")
        return False
    
    results_df, gap_df = load_results(results_dir)
    
    # Extract data
    identities = [70, 50, 40]
    identity_labels = ['70%', '50%', '40%']
    
    # Get values for each config
    final_acc = []
    mid_acc = []
    final_f1 = []
    mid_f1 = []
    
    for identity in identities:
        final_row = results_df[(results_df['Identity_Threshold'] == identity) & 
                               (results_df['Layer_Config'] == 'Final (Layer 33)')]
        mid_row = results_df[(results_df['Identity_Threshold'] == identity) & 
                             (results_df['Layer_Config'] == 'Mid (Layers 19-33)')]
        
        final_acc.append(final_row['Accuracy'].values[0])
        mid_acc.append(mid_row['Accuracy'].values[0])
        final_f1.append(final_row['Macro_F1'].values[0])
        mid_f1.append(mid_row['Macro_F1'].values[0])
    
    # Calculate gaps
    acc_gaps = [f - m for f, m in zip(final_acc, mid_acc)]
    f1_gaps = [f - m for f, m in zip(final_f1, mid_f1)]
    
    # Colors
    final_color = '#E63946'
    mid_color = '#2A9D8F'
    
    # Create figure
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), facecolor='white')
    
    for ax in axes:
        ax.set_facecolor('white')
    
    x = np.arange(len(identities))
    width = 0.35
    
    # Panel (a): Accuracy vs Homology Threshold
    ax1 = axes[0]
    bars1_final = ax1.bar(x - width/2, final_acc, width, label='Final (Layer 33)', 
                          color=final_color, edgecolor='#333333', linewidth=1)
    bars1_mid = ax1.bar(x + width/2, mid_acc, width, label='Mid (Layers 19-33)',
                        color=mid_color, edgecolor='#333333', linewidth=1)
    
    ax1.set_ylabel('Accuracy', fontsize=11, fontweight='medium')
    ax1.set_xlabel('Sequence Identity Threshold', fontsize=11, fontweight='medium')
    ax1.set_title('Test Accuracy vs Homology Constraint', fontsize=11, fontweight='medium')
    ax1.set_xticks(x)
    ax1.set_xticklabels(identity_labels)
    ax1.set_ylim(0.65, 0.90)
    ax1.legend(loc='upper right', frameon=True, framealpha=0.9)
    
    # Add value labels
    for bar in bars1_final:
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.008,
                f'{bar.get_height():.1%}', ha='center', va='bottom', fontsize=8)
    for bar in bars1_mid:
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.008,
                f'{bar.get_height():.1%}', ha='center', va='bottom', fontsize=8)
    
    ax1.text(-0.12, 1.08, '(a)', transform=ax1.transAxes, fontsize=14,
             fontweight='bold', va='top', ha='left')
    
    for spine in ['top', 'right']:
        ax1.spines[spine].set_visible(False)
    ax1.yaxis.grid(True, linestyle='-', alpha=0.2, linewidth=0.5)
    
    # Add arrow showing degradation direction
    ax1.annotate('', xy=(2.3, 0.74), xytext=(0.3, 0.82),
                arrowprops=dict(arrowstyle='->', color='#666666', 
                               linewidth=1.5, linestyle='--'))
    ax1.text(1.3, 0.70, 'More stringent\n(harder test)', fontsize=8, 
             ha='center', color='#666666', style='italic')
    
    # Panel (b): Macro-F1 vs Homology Threshold
    ax2 = axes[1]
    bars2_final = ax2.bar(x - width/2, final_f1, width, label='Final (Layer 33)',
                          color=final_color, edgecolor='#333333', linewidth=1)
    bars2_mid = ax2.bar(x + width/2, mid_f1, width, label='Mid (Layers 19-33)',
                        color=mid_color, edgecolor='#333333', linewidth=1)
    
    ax2.set_ylabel('Macro-F1', fontsize=11, fontweight='medium')
    ax2.set_xlabel('Sequence Identity Threshold', fontsize=11, fontweight='medium')
    ax2.set_title('Macro-F1 vs Homology Constraint', fontsize=11, fontweight='medium')
    ax2.set_xticks(x)
    ax2.set_xticklabels(identity_labels)
    ax2.set_ylim(0.55, 0.85)
    ax2.legend(loc='upper right', frameon=True, framealpha=0.9)
    
    for bar in bars2_final:
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.008,
                f'{bar.get_height():.3f}', ha='center', va='bottom', fontsize=8)
    for bar in bars2_mid:
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.008,
                f'{bar.get_height():.3f}', ha='center', va='bottom', fontsize=8)
    
    ax2.text(-0.12, 1.08, '(b)', transform=ax2.transAxes, fontsize=14,
             fontweight='bold', va='top', ha='left')
    
    for spine in ['top', 'right']:
        ax2.spines[spine].set_visible(False)
    ax2.yaxis.grid(True, linestyle='-', alpha=0.2, linewidth=0.5)
    
    # Panel (c): Performance Gap (Final - Mid)
    ax3 = axes[2]
    
    # Bar colors based on sign (positive = final better, negative = mid better)
    gap_colors = ['#E63946' if g > 0 else '#2A9D8F' for g in acc_gaps]
    
    bars3 = ax3.bar(x, [g * 100 for g in acc_gaps], width=0.5, 
                    color=gap_colors, edgecolor='#333333', linewidth=1)
    
    ax3.axhline(y=0, color='#333333', linewidth=1, linestyle='-')
    
    ax3.set_ylabel('Accuracy Gap (percentage points)', fontsize=11, fontweight='medium')
    ax3.set_xlabel('Sequence Identity Threshold', fontsize=11, fontweight='medium')
    ax3.set_title('Performance Gap\n(Final − Mid)', fontsize=11, fontweight='medium')
    ax3.set_xticks(x)
    ax3.set_xticklabels(identity_labels)
    ax3.set_ylim(-5, 5)
    
    # Add value labels
    for bar, gap in zip(bars3, acc_gaps):
        offset = 0.3 if gap >= 0 else -0.5
        ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + offset,
                f'{gap*100:+.1f}%', ha='center', va='bottom' if gap >= 0 else 'top',
                fontsize=9, fontweight='bold')
    
    # Add annotations
    ax3.text(0.5, 4.2, 'Final layer better', fontsize=8, color='#E63946', 
             ha='center', style='italic')
    ax3.text(0.5, -4.2, 'Mid-layers better', fontsize=8, color='#2A9D8F',
             ha='center', style='italic')
    
    ax3.text(-0.12, 1.08, '(c)', transform=ax3.transAxes, fontsize=14,
             fontweight='bold', va='top', ha='left')
    
    for spine in ['top', 'right']:
        ax3.spines[spine].set_visible(False)
    ax3.yaxis.grid(True, linestyle='-', alpha=0.2, linewidth=0.5)
    
    # Add key takeaway annotation
    takeaway = ("Key finding: Performance degrades gracefully with stricter homology constraints,\n"
                "demonstrating learned functional abstractions rather than sequence memorization.")
    fig.text(0.5, -0.02, takeaway, ha='center', va='top', fontsize=10,
             style='italic', color='#333333',
             bbox=dict(boxstyle='round,pad=0.4', facecolor='#F0F8FF',
                      edgecolor='#4472C4', linewidth=1))
    
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.15)
    
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
        'Identity_Threshold': identity_labels,
        'Final_Accuracy': final_acc,
        'Mid_Accuracy': mid_acc,
        'Final_Macro_F1': final_f1,
        'Mid_Macro_F1': mid_f1,
        'Accuracy_Gap': acc_gaps,
        'Macro_F1_Gap': f1_gaps
    })
    summary_data.to_csv(f"{output_file}_data.csv", index=False)
    
    print(f"\n✅ Figure 5 saved:")
    print(f"   {output_file}.png (300 dpi)")
    print(f"   {output_file}.pdf (vector)")
    print(f"   {output_file}_data.csv (underlying data)")
    
    # Print summary
    print("\n📊 Summary:")
    print(f"   70% identity: Final {final_acc[0]:.1%} vs Mid {mid_acc[0]:.1%} (gap: {acc_gaps[0]*100:+.1f}%)")
    print(f"   50% identity: Final {final_acc[1]:.1%} vs Mid {mid_acc[1]:.1%} (gap: {acc_gaps[1]*100:+.1f}%)")
    print(f"   40% identity: Final {final_acc[2]:.1%} vs Mid {mid_acc[2]:.1%} (gap: {acc_gaps[2]*100:+.1f}%)")
    
    return True


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description='Generate Figure 5: Generalization Under Homology Constraints'
    )
    parser.add_argument('--output', default='figures_output/figure5_homology_generalization',
                        help='Output file path (without extension)')
    args = parser.parse_args()
    
    plot_figure5(args.output)

