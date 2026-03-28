#!/usr/bin/env python3
"""
Step 16: Figure Generation

Generates 6 publication-quality figures saved to:
  - {run_dir}/figures/          (per-run archive)
  - figures_output/             (project-level stable copy)

Figures:
  Figure1_clustering_ari.png        ARI bar chart across 4 layer configs
  Figure2_confusion_matrix.png      8×8 confusion matrix (split40, layer33)
  Figure3_homology_classification.png  Accuracy vs identity threshold (line chart)
  Figure4_pooling_comparison.png    Mean vs CLS pooling comparison
  Figure5_calibration.png           Reliability diagram before/after Platt scaling
  Figure6_retrieval_pr.png          Precision-recall curve at cosine-sim thresholds

Usage:
    python pipeline/step_16_figures.py --run-dir runs/2025-01-01_000000/
"""

import argparse
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
from pathlib import Path
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

# Project root
PROJECT_ROOT = Path(__file__).parent.parent

# ── Style ──────────────────────────────────────────────────────────────────────
COLORS = {
    'layer33_mean':     '#2E86AB',   # blue
    'layers20_30_mean': '#A23B72',   # magenta
    'layers20_33_mean': '#F18F01',   # orange
    'layer33_cls':      '#C73E1D',   # red
    'accent':           '#4CAF50',   # green accent
}

LAYER_LABELS = {
    'layer33_mean':     'Layer 33\n(mean)',
    'layers20_30_mean': 'Layers 20–30\n(mean)',
    'layers20_33_mean': 'Layers 20–33\n(mean)',
    'layer33_cls':      'Layer 33\n(CLS)',
}

CLASS_COLORS = {
    'TK':       '#e41a1c',
    'CMGC':     '#377eb8',
    'CAMK':     '#4daf4a',
    'AGC':      '#984ea3',
    'STE':      '#ff7f00',
    'TKL':      '#ffff33',
    'CK1':      '#a65628',
    'Atypical': '#f781bf',
}

plt.rcParams.update({
    'font.family':      'sans-serif',
    'font.size':        11,
    'axes.titlesize':   12,
    'axes.labelsize':   11,
    'xtick.labelsize':  10,
    'ytick.labelsize':  10,
    'legend.fontsize':  10,
    'figure.dpi':       150,
    'savefig.dpi':      300,
    'savefig.bbox':     'tight',
    'axes.spines.top':  False,
    'axes.spines.right': False,
})


def load_json(path):
    if Path(path).exists():
        with open(path) as f:
            return json.load(f)
    return None


# ── Figure 1: Clustering ARI bar chart ────────────────────────────────────────
def generate_figure1(run_dir: Path, out_dir: Path):
    print("  Generating Figure1_clustering_ari.png...")

    registry = load_json(run_dir / "results" / "clustering" / "clustering_registry.json")
    if not registry:
        print("    ⚠ clustering_registry.json not found — skipping Figure 1")
        return

    configs   = ['layer33_mean', 'layer33_cls', 'layers20_30_mean', 'layers20_33_mean']
    ari_vals  = []
    nmi_vals  = []
    labels    = []
    colors    = []
    improv    = []

    baseline_ari = registry["experiments"].get("layer33_mean", {}).get("metrics", {}).get("ARI", 0)

    for cfg in configs:
        exp = registry["experiments"].get(cfg)
        if exp is None:
            continue
        ari = exp["metrics"]["ARI"]
        nmi = exp["metrics"]["NMI"]
        ari_vals.append(ari)
        nmi_vals.append(nmi)
        labels.append(LAYER_LABELS.get(cfg, cfg))
        colors.append(COLORS.get(cfg, '#888888'))
        imp = ((ari - baseline_ari) / baseline_ari * 100) if baseline_ari > 0 else 0
        improv.append(imp)

    fig, ax = plt.subplots(figsize=(8, 5))
    x     = np.arange(len(labels))
    bars  = ax.bar(x, ari_vals, color=colors, edgecolor='black', linewidth=0.6, width=0.55)

    # Value + improvement labels
    for bar, val, imp in zip(bars, ari_vals, improv):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.008,
                f'{val:.3f}',
                ha='center', va='bottom', fontsize=10, fontweight='bold')
        if abs(imp) > 0.1:
            sign = '+' if imp > 0 else ''
            ax.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() / 2,
                    f'{sign}{imp:.0f}%',
                    ha='center', va='center', fontsize=9,
                    color='white', fontweight='bold')

    ax.set_xticks(x)
    ax.set_xticklabels(labels, ha='center')
    ax.set_ylabel('Adjusted Rand Index (ARI)')
    ax.set_title('Figure 1. Clustering Performance (ARI) Across ESM-2 Layer Configurations',
                 pad=12)
    ax.set_ylim(0, max(ari_vals) * 1.25)
    ax.axhline(baseline_ari, color='grey', linestyle='--', linewidth=0.8,
               label=f'Baseline (Layer 33 mean) = {baseline_ari:.3f}')
    ax.legend(loc='upper left', frameon=False)

    plt.tight_layout()
    outpath = out_dir / "Figure1_clustering_ari.png"
    plt.savefig(outpath)
    plt.close()
    print(f"    ✓ {outpath.name}")


# ── Figure 2: Confusion matrix ─────────────────────────────────────────────────
def generate_figure2(run_dir: Path, out_dir: Path):
    print("  Generating Figure2_confusion_matrix.png...")

    cm_file = run_dir / "results" / "supervised" / "lr_split40_confusion_layer33_mean.json"
    if not cm_file.exists():
        print("    ⚠ confusion matrix file not found — skipping Figure 2")
        return

    data    = load_json(cm_file)
    cm      = np.array(data["confusion_matrix"])
    classes = data["classes"]

    # Row-normalise (recall per class)
    row_sums = cm.sum(axis=1, keepdims=True)
    cm_norm  = np.where(row_sums > 0, cm / row_sums, 0)

    fig, ax = plt.subplots(figsize=(8, 7))
    im = ax.imshow(cm_norm, cmap='Blues', vmin=0, vmax=1, aspect='auto')
    plt.colorbar(im, ax=ax, fraction=0.046, label='Recall (row-normalised)')

    ax.set_xticks(range(len(classes)))
    ax.set_yticks(range(len(classes)))
    ax.set_xticklabels(classes, rotation=45, ha='right')
    ax.set_yticklabels(classes)
    ax.set_xlabel('Predicted Label')
    ax.set_ylabel('True Label')
    ax.set_title('Figure 2. Confusion Matrix — Supervised Classification\n'
                 '(Layer 33 mean, 40% Identity Split, 8 Kinase Families)',
                 pad=10)

    # Cell annotations
    for i in range(len(classes)):
        for j in range(len(classes)):
            val = cm_norm[i, j]
            raw = int(cm[i, j])
            color = 'white' if val > 0.55 else 'black'
            ax.text(j, i, f'{val:.2f}\n({raw})',
                    ha='center', va='center', fontsize=7.5, color=color)

    plt.tight_layout()
    outpath = out_dir / "Figure2_confusion_matrix.png"
    plt.savefig(outpath)
    plt.close()
    print(f"    ✓ {outpath.name}")


# ── Figure 3: Homology classification line chart ──────────────────────────────
def generate_figure3(run_dir: Path, out_dir: Path):
    print("  Generating Figure3_homology_classification.png...")

    # Prefer calibrated accuracy; fall back to supervised registry
    configs_plotted = ['layer33_mean', 'layers20_30_mean']
    thresholds      = [70, 50, 40]
    accuracy_data   = {cfg: [] for cfg in configs_plotted}

    for thresh in thresholds:
        split = f"split{thresh}"
        for cfg in configs_plotted:
            # Try per-config calibration file first
            cal_path = run_dir / "results" / "calibration" / f"{split}_calibration_{cfg}.json"
            cal      = load_json(cal_path)
            if cal:
                acc = cal["calibrated"]["accuracy"]
            else:
                # Fall back to supervised registry (uncalibrated)
                sup_path = run_dir / "results" / "supervised" / f"lr_{split}_metrics.json"
                sup      = load_json(sup_path)
                if sup and cfg in sup:
                    entry = sup[cfg]
                    # Handle both flat and nested {'metrics': {...}} structures
                    if "metrics" in entry:
                        acc = entry["metrics"].get("accuracy")
                    else:
                        acc = entry.get("accuracy")
                else:
                    acc = None
            accuracy_data[cfg].append(acc)

    fig, ax = plt.subplots(figsize=(7, 5))

    style = {
        'layer33_mean':     dict(color=COLORS['layer33_mean'],     marker='o', linewidth=2, label='Layer 33 (mean)'),
        'layers20_30_mean': dict(color=COLORS['layers20_30_mean'], marker='s', linewidth=2, label='Layers 20–30 (mean)'),
    }

    for cfg in configs_plotted:
        vals = accuracy_data[cfg]
        valid_x = [thresholds[i] for i, v in enumerate(vals) if v is not None]
        valid_y = [v for v in vals if v is not None]
        if valid_y:
            ax.plot(valid_x, valid_y, **style[cfg])
            for xi, yi in zip(valid_x, valid_y):
                ax.annotate(f'{yi:.3f}', (xi, yi),
                            textcoords='offset points', xytext=(0, 8),
                            ha='center', fontsize=9)

    ax.set_xlabel('Homology Identity Threshold (%)')
    ax.set_ylabel('Calibrated Accuracy')
    ax.set_title('Figure 3. Classification Performance Across\nHomology Identity Thresholds', pad=10)
    ax.set_xticks(thresholds)
    ax.set_xticklabels(['70%', '50%', '40%'])
    ax.set_xlim(75, 35)   # reversed: strict (40%) on right
    ax.set_ylim(0.5, 1.0)
    ax.legend(frameon=False)
    ax.grid(axis='y', alpha=0.3)
    ax.annotate('← Stricter split', xy=(0.98, 0.04), xycoords='axes fraction',
                ha='right', fontsize=9, color='grey')

    plt.tight_layout()
    outpath = out_dir / "Figure3_homology_classification.png"
    plt.savefig(outpath)
    plt.close()
    print(f"    ✓ {outpath.name}")


# ── Figure 4: Pooling strategy comparison ─────────────────────────────────────
def generate_figure4(run_dir: Path, out_dir: Path):
    print("  Generating Figure4_pooling_comparison.png...")

    registry   = load_json(run_dir / "results" / "clustering" / "clustering_registry.json")
    cal_33     = load_json(run_dir / "results" / "calibration" / "split40_calibration_layer33_mean.json")

    if not registry:
        print("    ⚠ clustering_registry.json not found — skipping Figure 4")
        return

    # Clustering ARI
    ari_mean = registry["experiments"].get("layer33_mean", {}).get("metrics", {}).get("ARI", 0)
    ari_cls  = registry["experiments"].get("layer33_cls",  {}).get("metrics", {}).get("ARI", 0)
    ari_mid  = registry["experiments"].get("layers20_30_mean", {}).get("metrics", {}).get("ARI", 0)

    # Supervised accuracy (calibrated if available, else from supervised registry)
    def get_acc(cfg):
        cal = load_json(run_dir / "results" / "calibration" / f"split40_calibration_{cfg}.json")
        if cal:
            return cal["calibrated"]["accuracy"]
        sup = load_json(run_dir / "results" / "supervised" / "lr_split40_metrics.json")
        if sup and cfg in sup:
            entry = sup[cfg]
            if "metrics" in entry:
                return entry["metrics"].get("accuracy")
            return entry.get("accuracy")
        return None

    acc_mean = get_acc("layer33_mean")
    acc_mid  = get_acc("layers20_30_mean")

    # Layout: 2 groups (Clustering ARI | Classification Accuracy)
    # Each group: 2–3 bars
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))

    # Left: Clustering ARI - mean vs CLS vs mid-layer
    ax1 = axes[0]
    bars_cfg = ['layer33_mean', 'layer33_cls', 'layers20_30_mean']
    bars_ari = [ari_mean, ari_cls, ari_mid]
    bars_lbl = ['Layer 33\n(mean)', 'Layer 33\n(CLS)', 'Layers 20–30\n(mean)']
    bars_col = [COLORS['layer33_mean'], COLORS['layer33_cls'], COLORS['layers20_30_mean']]

    b1 = ax1.bar(bars_lbl, bars_ari, color=bars_col, edgecolor='black', linewidth=0.6, width=0.5)
    for bar, val in zip(b1, bars_ari):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.005,
                 f'{val:.3f}', ha='center', va='bottom', fontsize=10, fontweight='bold')
    ax1.set_ylabel('Adjusted Rand Index (ARI)')
    ax1.set_title('(a) Clustering Performance')
    ax1.set_ylim(0, max(bars_ari) * 1.3)

    # Right: Classification Accuracy - mean vs mid-layer (CLS not run for classification)
    ax2 = axes[1]
    acc_vals = [v for v in [acc_mean, acc_mid] if v is not None]
    acc_lbls = ['Layer 33\n(mean pooling)', 'Layers 20–30\n(mean pooling)'][:len(acc_vals)]
    acc_cols = [COLORS['layer33_mean'], COLORS['layers20_30_mean']][:len(acc_vals)]

    if acc_vals:
        b2 = ax2.bar(acc_lbls, acc_vals, color=acc_cols, edgecolor='black', linewidth=0.6, width=0.4)
        for bar, val in zip(b2, acc_vals):
            ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.005,
                     f'{val:.3f}', ha='center', va='bottom', fontsize=10, fontweight='bold')
        ax2.set_ylabel('Calibrated Accuracy')
        ax2.set_title('(b) Supervised Classification')
        ax2.set_ylim(0.5, min(1.0, max(acc_vals) * 1.15))

    fig.suptitle('Figure 4. Effect of Pooling Strategy on Performance\n'
                 '(40% Identity Split)', y=1.02)
    plt.tight_layout()
    outpath = out_dir / "Figure4_pooling_comparison.png"
    plt.savefig(outpath)
    plt.close()
    print(f"    ✓ {outpath.name}")


# ── Figure 5: Calibration reliability diagram ─────────────────────────────────
def generate_figure5(run_dir: Path, out_dir: Path):
    print("  Generating Figure5_calibration.png...")

    cal = load_json(run_dir / "results" / "calibration" / "split40_calibration_layer33_mean.json")
    if not cal:
        # Try compat file
        cal_compat = load_json(run_dir / "results" / "calibration" / "split40_calibration.json")
        if not cal_compat:
            print("    ⚠ calibration file not found — skipping Figure 5")
            return
        print("    ⚠ using compat calibration file (no per-bin data)")
        # Build a minimal synthetic reliability plot
        fig, ax = plt.subplots(figsize=(6, 5))
        ece_uncal = cal_compat.get("uncalibrated_ece", 0)
        ece_cal   = cal_compat.get("calibrated_ece", 0)
        ax.plot([0, 1], [0, 1], 'k--', label='Perfect calibration')
        ax.set_xlabel('Mean Predicted Probability')
        ax.set_ylabel('Fraction of Positives')
        ax.set_title(f'Figure 5. Calibration\nECE: {ece_uncal:.3f} → {ece_cal:.3f} (Platt scaling)')
        ax.legend(frameon=False)
        plt.tight_layout()
        plt.savefig(out_dir / "Figure5_calibration.png")
        plt.close()
        return

    uncal_bins = cal["uncalibrated"].get("reliability_bins", [])
    cal_bins   = cal["calibrated"].get("reliability_bins", [])
    ece_uncal  = cal["uncalibrated"]["ece"]
    ece_cal    = cal["calibrated"]["ece"]
    acc_uncal  = cal["uncalibrated"]["accuracy"]
    acc_cal    = cal["calibrated"]["accuracy"]

    fig, ax = plt.subplots(figsize=(6, 5.5))

    # Perfect calibration diagonal
    ax.plot([0, 1], [0, 1], 'k--', linewidth=1.2, label='Perfect calibration', zorder=5)

    def plot_bins(bins, color, label, linestyle='-', alpha=0.8):
        if not bins:
            return
        conf = [b["mean_confidence"]  for b in bins if b["count"] > 0]
        frac = [b["fraction_positive"] for b in bins if b["count"] > 0]
        ax.plot(conf, frac, marker='o', markersize=5, color=color,
                linestyle=linestyle, linewidth=1.8, label=label, alpha=alpha)
        # Fill between curve and diagonal
        ax.fill_betweenx(frac, conf, conf,  alpha=0.05, color=color)

    plot_bins(uncal_bins, COLORS['layer33_mean'],     f'Before Platt scaling (ECE={ece_uncal:.3f})', '--')
    plot_bins(cal_bins,   COLORS['layers20_30_mean'], f'After Platt scaling  (ECE={ece_cal:.3f})',   '-')

    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel('Mean Predicted Probability')
    ax.set_ylabel('Fraction of Positives')
    ax.set_title('Figure 5. Calibration Before and After Platt Scaling\n'
                 f'(Layer 33, split40)  Accuracy: {acc_uncal:.3f} → {acc_cal:.3f}',
                 pad=10)
    ax.legend(frameon=False, loc='upper left')
    ax.grid(alpha=0.25)

    plt.tight_layout()
    outpath = out_dir / "Figure5_calibration.png"
    plt.savefig(outpath)
    plt.close()
    print(f"    ✓ {outpath.name}")


# ── Figure 6: Retrieval PR curve ──────────────────────────────────────────────
def generate_figure6(run_dir: Path, out_dir: Path):
    print("  Generating Figure6_retrieval_pr.png...")

    pr_data     = load_json(run_dir / "results" / "retrieval" / "split40_pr_curve.json")
    basic_data  = load_json(run_dir / "results" / "retrieval" / "split40_retrieval.json")

    fig, ax = plt.subplots(figsize=(7, 5))

    if pr_data:
        thresholds = pr_data["thresholds"]
        precision  = pr_data["precision"]
        recall     = pr_data["recall"]
        p_at_1     = pr_data.get("p_at_1", None)
        p_at_3     = pr_data.get("p_at_3", None)

        # Sort by recall ascending for a clean curve
        pairs = sorted(zip(recall, precision))
        rec_sorted  = [p[0] for p in pairs]
        prec_sorted = [p[1] for p in pairs]

        ax.plot(rec_sorted, prec_sorted, color=COLORS['layer33_mean'],
                linewidth=2, label='Precision–Recall curve (cosine-sim thresholds)')

        # Reference lines for P@1, P@3
        if p_at_1:
            ax.axhline(p_at_1, color=COLORS['layers20_30_mean'], linestyle='--',
                       linewidth=1.2, label=f'P@1 = {p_at_1:.3f}')
        if p_at_3:
            ax.axhline(p_at_3, color=COLORS['layer33_cls'], linestyle=':',
                       linewidth=1.2, label=f'P@3 = {p_at_3:.3f}')

        ax.set_xlabel('Recall')
        ax.set_ylabel('Precision')

    elif basic_data:
        # Fallback: bar chart of P@k metrics
        metrics = basic_data.get("metrics", {})
        ks   = [k for k in ['P@1', 'P@3', 'P@5', 'P@10', 'MRR'] if k in metrics]
        vals = [metrics[k] for k in ks]
        ax.bar(ks, vals, color=COLORS['layer33_mean'], edgecolor='black', linewidth=0.6)
        for xi, val in enumerate(vals):
            ax.text(xi, val + 0.01, f'{val:.3f}', ha='center', fontsize=10, fontweight='bold')
        ax.set_ylabel('Score')
        ax.set_ylim(0, 1.1)
    else:
        print("    ⚠ No retrieval data found — skipping Figure 6")
        plt.close()
        return

    if basic_data:
        m    = basic_data.get("metrics", {})
        mrr  = m.get("MRR", None)
        info = f'MRR = {mrr:.3f}' if mrr else ''
        ax.set_title(f'Figure 6. Exemplar-Based Retrieval\n'
                     f'(Layer 33, split40)  {info}', pad=10)

    ax.legend(frameon=False)
    ax.grid(alpha=0.25)

    plt.tight_layout()
    outpath = out_dir / "Figure6_retrieval_pr.png"
    plt.savefig(outpath)
    plt.close()
    print(f"    ✓ {outpath.name}")


# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description='Generate manuscript figures')
    parser.add_argument('--run-dir', type=str, required=True)
    args    = parser.parse_args()
    run_dir = Path(args.run_dir)

    print("=" * 60)
    print("Step 16: Figure Generation")
    print("=" * 60)
    print(f"Run directory: {run_dir}")

    # Per-run figures dir
    figures_dir = run_dir / "figures"
    figures_dir.mkdir(exist_ok=True)

    # Project-level stable output
    project_out = PROJECT_ROOT / "figures_output"
    project_out.mkdir(exist_ok=True)

    print("\nGenerating figures...")

    generators = [
        generate_figure1,
        generate_figure2,
        generate_figure3,
        generate_figure4,
        generate_figure5,
        generate_figure6,
    ]

    generated = []
    for gen in generators:
        try:
            gen(run_dir, figures_dir)
        except Exception as e:
            print(f"    ✗ {gen.__name__} failed: {e}")

    # Copy PNGs to project-level figures_output/
    print("\nCopying to figures_output/...")
    import shutil
    for png in sorted(figures_dir.glob("Figure*.png")):
        dest = project_out / png.name
        shutil.copy2(png, dest)
        generated.append(png.name)
        print(f"  ✓ figures_output/{png.name}")

    # Write figure registry
    registry = {
        "step": 16,
        "name": "Figure Generation",
        "timestamp": datetime.now().isoformat(),
        "output_dirs": [str(figures_dir), str(project_out)],
        "figures": {
            "Figure1_clustering_ari.png": {
                "title": "Clustering Performance (ARI) Across ESM-2 Layer Configurations",
                "cited_in": "Results §3.1",
                "source": "results/clustering/clustering_registry.json"
            },
            "Figure2_confusion_matrix.png": {
                "title": "Confusion Matrix — Supervised Classification (Layer 33, split40)",
                "cited_in": "Results §3.2",
                "source": "results/supervised/lr_split40_confusion_layer33_mean.json"
            },
            "Figure3_homology_classification.png": {
                "title": "Classification Accuracy Across Homology Identity Thresholds",
                "cited_in": "Results §3.3",
                "source": "results/calibration/split40_calibration_*.json"
            },
            "Figure4_pooling_comparison.png": {
                "title": "Effect of Pooling Strategy on Clustering and Classification",
                "cited_in": "Results §3.4",
                "source": "results/clustering/clustering_registry.json + calibration"
            },
            "Figure5_calibration.png": {
                "title": "Calibration Before and After Platt Scaling",
                "cited_in": "Results §3.5",
                "source": "results/calibration/split40_calibration_layer33_mean.json"
            },
            "Figure6_retrieval_pr.png": {
                "title": "Exemplar-Based Retrieval: Precision-Recall at Cosine Similarity Thresholds",
                "cited_in": "Results §3.6",
                "source": "results/retrieval/split40_pr_curve.json"
            }
        }
    }

    registry_file = figures_dir / "figure_registry.json"
    with open(registry_file, 'w') as f:
        json.dump(registry, f, indent=2)
    print(f"\n✓ Registry: {registry_file}")

    # Also copy registry to project-level
    import shutil
    shutil.copy2(registry_file, project_out / "figure_registry.json")

    print("\n" + "=" * 60)
    print("Step 16 COMPLETE")
    print("=" * 60)
    print(f"\nFigures saved to:")
    print(f"  {figures_dir}")
    print(f"  {project_out}")
    for png in sorted(figures_dir.glob("Figure*.png")):
        print(f"    ✓ {png.name}")


if __name__ == "__main__":
    main()
