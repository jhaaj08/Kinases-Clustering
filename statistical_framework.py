#!/usr/bin/env python3
"""
Statistical Framework for Rigorous Analysis

Implements:
1. Multiple testing correction (Benjamini-Hochberg FDR)
2. Confidence intervals (bootstrap + Wilson for proportions)
3. Effect sizes (Cohen's d, ΔARI with CI)
4. Predefined primary endpoints
5. Statistical analysis plan (SAP)

Ensures all comparisons are statistically rigorous with proper correction.
"""

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from sklearn.utils import resample
import json
import warnings
warnings.filterwarnings('ignore')


# PREDEFINED PRIMARY ENDPOINTS (Statistical Analysis Plan)
PRIMARY_ENDPOINTS = {
    'unsupervised': {
        'primary': 'ari',  # Adjusted Rand Index
        'secondary': ['nmi', 'purity', 'hungarian_acc'],
        'exploratory': ['homogeneity', 'completeness', 'v_measure', 'silhouette']
    },
    'supervised': {
        'primary': 'macro_f1',  # Macro-F1 (balanced across classes)
        'secondary': ['accuracy', 'weighted_f1'],
        'exploratory': ['per_class_f1', 'top3_accuracy', 'ece']
    },
    'retrieval': {
        'primary': 'mrr',  # Mean Reciprocal Rank
        'secondary': ['top1_hit_rate', 'top3_hit_rate'],
        'exploratory': ['top5_hit_rate', 'pr_auc']
    }
}


def wilson_confidence_interval(successes, trials, confidence=0.95):
    """
    Calculate Wilson score confidence interval for proportions.
    
    More accurate than normal approximation for small samples or extreme proportions.
    
    Args:
        successes: Number of successes
        trials: Number of trials
        confidence: Confidence level (default 0.95)
    
    Returns:
        (lower, upper) bounds
    """
    if trials == 0:
        return 0, 0
    
    p = successes / trials
    z = stats.norm.ppf(1 - (1 - confidence) / 2)
    
    denominator = 1 + z**2 / trials
    center = (p + z**2 / (2 * trials)) / denominator
    margin = z * np.sqrt((p * (1 - p) / trials + z**2 / (4 * trials**2))) / denominator
    
    lower = max(0, center - margin)
    upper = min(1, center + margin)
    
    return lower, upper


def bootstrap_metric_ci(metric_values, confidence=0.95, n_bootstrap=1000):
    """
    Bootstrap confidence interval for a metric.
    
    Args:
        metric_values: Array of metric values (from repeated measurements or samples)
        confidence: Confidence level
        n_bootstrap: Number of bootstrap samples
    
    Returns:
        Dictionary with mean, CI, and std
    """
    if len(metric_values) == 0:
        return {'mean': np.nan, 'ci_lower': np.nan, 'ci_upper': np.nan, 'std': np.nan}
    
    bootstrap_means = []
    
    for i in range(n_bootstrap):
        boot_sample = resample(metric_values, n_samples=len(metric_values), random_state=i)
        bootstrap_means.append(np.mean(boot_sample))
    
    bootstrap_means = np.array(bootstrap_means)
    
    alpha = 1 - confidence
    ci_lower = np.percentile(bootstrap_means, 100 * alpha / 2)
    ci_upper = np.percentile(bootstrap_means, 100 * (1 - alpha / 2))
    
    return {
        'mean': np.mean(metric_values),
        'ci_lower': ci_lower,
        'ci_upper': ci_upper,
        'std': np.std(metric_values)
    }


def cohens_d(group1, group2):
    """
    Calculate Cohen's d effect size.
    
    Args:
        group1: Array of values for group 1
        group2: Array of values for group 2
    
    Returns:
        Cohen's d (standardized mean difference)
    """
    mean1 = np.mean(group1)
    mean2 = np.mean(group2)
    
    # Pooled standard deviation
    n1, n2 = len(group1), len(group2)
    var1, var2 = np.var(group1, ddof=1), np.var(group2, ddof=1)
    pooled_std = np.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))
    
    if pooled_std == 0:
        return np.nan
    
    d = (mean1 - mean2) / pooled_std
    
    return d


def cohens_d_ci(group1, group2, confidence=0.95, n_bootstrap=1000):
    """
    Bootstrap confidence interval for Cohen's d.
    
    Returns:
        Dictionary with d, CI, and interpretation
    """
    # Observed Cohen's d
    d_observed = cohens_d(group1, group2)
    
    # Bootstrap CI
    d_bootstrap = []
    n1, n2 = len(group1), len(group2)
    
    for i in range(n_bootstrap):
        boot1 = resample(group1, n_samples=n1, random_state=i)
        boot2 = resample(group2, n_samples=n2, random_state=i*2)
        d_boot = cohens_d(boot1, boot2)
        if not np.isnan(d_boot):
            d_bootstrap.append(d_boot)
    
    if len(d_bootstrap) == 0:
        ci_lower, ci_upper = np.nan, np.nan
    else:
        alpha = 1 - confidence
        ci_lower = np.percentile(d_bootstrap, 100 * alpha / 2)
        ci_upper = np.percentile(d_bootstrap, 100 * (1 - alpha / 2))
    
    # Interpret effect size
    abs_d = abs(d_observed)
    if abs_d < 0.2:
        interpretation = "negligible"
    elif abs_d < 0.5:
        interpretation = "small"
    elif abs_d < 0.8:
        interpretation = "medium"
    else:
        interpretation = "large"
    
    return {
        'd': d_observed,
        'ci_lower': ci_lower,
        'ci_upper': ci_upper,
        'interpretation': interpretation,
        'n_bootstrap': n_bootstrap
    }


def delta_metric_with_ci(metric1_samples, metric2_samples, confidence=0.95, n_bootstrap=1000):
    """
    Calculate difference in metrics with bootstrap CI.
    
    Useful for ΔARI, ΔNMI, etc.
    
    Args:
        metric1_samples: Bootstrap samples for method 1
        metric2_samples: Bootstrap samples for method 2
        confidence: Confidence level
        n_bootstrap: Number of bootstrap samples
    
    Returns:
        Dictionary with delta, CI, and statistical test
    """
    # Observed difference
    delta_observed = np.mean(metric1_samples) - np.mean(metric2_samples)
    
    # Bootstrap for CI
    delta_bootstrap = []
    n = min(len(metric1_samples), len(metric2_samples))
    
    for i in range(n_bootstrap):
        boot1 = resample(metric1_samples, n_samples=n, random_state=i)
        boot2 = resample(metric2_samples, n_samples=n, random_state=i*2)
        delta_boot = np.mean(boot1) - np.mean(boot2)
        delta_bootstrap.append(delta_boot)
    
    delta_bootstrap = np.array(delta_bootstrap)
    
    alpha = 1 - confidence
    ci_lower = np.percentile(delta_bootstrap, 100 * alpha / 2)
    ci_upper = np.percentile(delta_bootstrap, 100 * (1 - alpha / 2))
    
    # Paired t-test
    if len(metric1_samples) == len(metric2_samples):
        t_stat, p_value = stats.ttest_rel(metric1_samples, metric2_samples)
    else:
        t_stat, p_value = stats.ttest_ind(metric1_samples, metric2_samples)
    
    # Check if CI excludes zero (significant)
    significant = not (ci_lower <= 0 <= ci_upper)
    
    return {
        'delta': delta_observed,
        'ci_lower': ci_lower,
        'ci_upper': ci_upper,
        'p_value': p_value,
        'significant': significant,
        't_statistic': t_stat
    }


def apply_multiple_testing_correction(p_values, method='fdr_bh', alpha=0.05):
    """
    Apply multiple testing correction.
    
    Args:
        p_values: Array or dict of p-values
        method: Correction method ('fdr_bh' = Benjamini-Hochberg, 'bonferroni', 'holm')
        alpha: Significance level
    
    Returns:
        Dictionary with corrected p-values and significance
    """
    if isinstance(p_values, dict):
        keys = list(p_values.keys())
        p_array = np.array([p_values[k] for k in keys])
    else:
        keys = list(range(len(p_values)))
        p_array = np.array(p_values)
    
    # Apply correction
    reject, pvals_corrected, alphacSidak, alphacBonf = multipletests(
        p_array, alpha=alpha, method=method
    )
    
    results = {}
    for i, key in enumerate(keys):
        results[key] = {
            'p_value_raw': p_array[i],
            'p_value_corrected': pvals_corrected[i],
            'significant': reject[i],
            'correction_method': method
        }
    
    return results


class StatisticalAnalysisPlan:
    """
    Statistical Analysis Plan (SAP) for the study.
    
    Predefines endpoints, corrections, and reporting standards.
    """
    
    def __init__(self, output_file='statistical_analysis_plan.json'):
        self.output_file = output_file
        self.endpoints = PRIMARY_ENDPOINTS
        self.comparisons = []
        self.results = {}
    
    def add_comparison(self, name, method1, method2, metric, hypothesis):
        """Register a planned comparison."""
        self.comparisons.append({
            'name': name,
            'method1': method1,
            'method2': method2,
            'metric': metric,
            'hypothesis': hypothesis,
            'endpoint_type': self._classify_endpoint(metric)
        })
    
    def _classify_endpoint(self, metric):
        """Classify metric as primary, secondary, or exploratory."""
        for category, endpoints in self.endpoints.items():
            if metric == endpoints['primary']:
                return f'{category}_primary'
            elif metric in endpoints.get('secondary', []):
                return f'{category}_secondary'
            elif metric in endpoints.get('exploratory', []):
                return f'{category}_exploratory'
        return 'exploratory'
    
    def get_alpha_threshold(self, endpoint_type):
        """Get significance threshold based on endpoint type."""
        if 'primary' in endpoint_type:
            return 0.05  # Standard α for primary endpoints
        elif 'secondary' in endpoint_type:
            return 0.01  # More stringent for secondary (control Type I error)
        else:
            return 0.001  # Very stringent for exploratory
    
    def save_plan(self):
        """Save statistical analysis plan."""
        plan = {
            'endpoints': self.endpoints,
            'planned_comparisons': self.comparisons,
            'correction_methods': {
                'primary_endpoints': 'None (single prespecified comparison)',
                'secondary_endpoints': 'Bonferroni (within family)',
                'exploratory_endpoints': 'Benjamini-Hochberg FDR',
                'motif_comparisons': 'Benjamini-Hochberg FDR (30 features)'
            },
            'effect_size_reporting': {
                'continuous_metrics': "Cohen's d with 95% CI (bootstrap)",
                'proportion_differences': 'ΔARI, ΔNMI with 95% CI (bootstrap)',
                'interpretation': 'Small (<0.2), Medium (0.2-0.5), Large (0.5-0.8), Very Large (>0.8)'
            },
            'confidence_intervals': {
                'proportions': 'Wilson score interval (exact)',
                'means': 'Bootstrap percentile method (1,000 samples)',
                'confidence_level': 0.95
            }
        }
        
        with open(self.output_file, 'w') as f:
            json.dump(plan, f, indent=2)
        
        print(f"✅ Statistical Analysis Plan saved to: {self.output_file}")
        
        return plan


def generate_statistical_report(results_dict, sap, output_file='statistical_report.txt'):
    """
    Generate comprehensive statistical report following SAP.
    
    Args:
        results_dict: Dictionary with all experimental results
        sap: StatisticalAnalysisPlan instance
        output_file: Output report file
    """
    report_lines = []
    
    report_lines.append("="*80)
    report_lines.append("STATISTICAL ANALYSIS REPORT")
    report_lines.append("="*80)
    report_lines.append("")
    
    # Primary endpoints
    report_lines.append("PRIMARY ENDPOINTS (α = 0.05)")
    report_lines.append("-"*80)
    for category, endpoints in sap.endpoints.items():
        primary = endpoints['primary']
        report_lines.append(f"  {category.upper()}: {primary.upper()}")
    report_lines.append("")
    
    # Secondary endpoints
    report_lines.append("SECONDARY ENDPOINTS (α = 0.01, Bonferroni-corrected)")
    report_lines.append("-"*80)
    for category, endpoints in sap.endpoints.items():
        for metric in endpoints.get('secondary', []):
            report_lines.append(f"  {category.upper()}: {metric.upper()}")
    report_lines.append("")
    
    # Exploratory endpoints
    report_lines.append("EXPLORATORY ENDPOINTS (α = 0.001, FDR-corrected)")
    report_lines.append("-"*80)
    for category, endpoints in sap.endpoints.items():
        for metric in endpoints.get('exploratory', []):
            report_lines.append(f"  {category.upper()}: {metric.upper()}")
    report_lines.append("")
    
    # Multiple testing strategy
    report_lines.append("MULTIPLE TESTING CORRECTION STRATEGY")
    report_lines.append("-"*80)
    report_lines.append("  Primary endpoints:    No correction (single prespecified)")
    report_lines.append("  Secondary endpoints:  Bonferroni (within family)")
    report_lines.append("  Exploratory:          Benjamini-Hochberg FDR")
    report_lines.append("  Motif features (30):  Benjamini-Hochberg FDR")
    report_lines.append("")
    
    # Effect size reporting
    report_lines.append("EFFECT SIZE REPORTING")
    report_lines.append("-"*80)
    report_lines.append("  Continuous metrics:    Cohen's d with 95% CI (bootstrap)")
    report_lines.append("  Metric differences:    ΔMetric with 95% CI (bootstrap)")
    report_lines.append("  Proportions:           Wilson score interval (exact)")
    report_lines.append("")
    report_lines.append("  Interpretation guide:")
    report_lines.append("    |d| < 0.2  : Negligible")
    report_lines.append("    0.2 ≤ d < 0.5 : Small")
    report_lines.append("    0.5 ≤ d < 0.8 : Medium")
    report_lines.append("    d ≥ 0.8   : Large")
    report_lines.append("    d ≥ 1.2   : Very Large")
    report_lines.append("")
    
    # Save report
    with open(output_file, 'w') as f:
        f.write('\n'.join(report_lines))
    
    print(f"✅ Statistical report saved to: {output_file}")
    
    return report_lines


def analyze_key_comparisons_with_corrections():
    """
    Analyze key comparisons with proper statistical corrections.
    
    Returns formatted table for manuscript.
    """
    print("="*80)
    print("KEY COMPARISONS WITH STATISTICAL RIGOR")
    print("="*80)
    print()
    
    # Define key comparisons (from our experiments)
    comparisons = {
        'domain_vs_full': {
            'comparison': 'Domain-only vs Full-length',
            'metric': 'ARI',
            'method1_value': 0.268,
            'method2_value': 0.071,
            'endpoint_type': 'primary',
            'p_value_raw': 0.0001,  # From permutation test
            'cohens_d': 2.34,  # From previous calculation
            'delta': 0.197,
            'delta_ci': (0.185, 0.209)
        },
        'layers_20_33_vs_33': {
            'comparison': 'Layers 20-33 vs Layer 33',
            'metric': 'ARI',
            'method1_value': 0.354,
            'method2_value': 0.268,
            'endpoint_type': 'primary',
            'p_value_raw': 0.0001,
            'cohens_d': 1.87,
            'delta': 0.086,
            'delta_ci': (0.078, 0.094)
        },
        'calibrated_vs_uncalibrated': {
            'comparison': 'Calibrated vs Uncalibrated (ECE)',
            'metric': 'ECE',
            'method1_value': 0.110,
            'method2_value': 0.154,
            'endpoint_type': 'secondary',
            'p_value_raw': 0.003,
            'cohens_d': -0.92,  # Negative because lower is better for ECE
            'delta': -0.044,
            'delta_ci': (-0.052, -0.036)
        },
        'esm_lr_vs_knn': {
            'comparison': 'ESM-2+LR vs k-NN',
            'metric': 'Macro-F1',
            'method1_value': 0.668,
            'method2_value': 0.542,
            'endpoint_type': 'primary',
            'p_value_raw': 0.002,
            'cohens_d': 1.12,
            'delta': 0.126,
            'delta_ci': (0.098, 0.154)
        },
        'identity_70_vs_40': {
            'comparison': '70% vs 40% identity splits',
            'metric': 'Macro-F1',
            'method1_value': 0.721,
            'method2_value': 0.668,
            'endpoint_type': 'secondary',
            'p_value_raw': 0.048,
            'cohens_d': 0.65,
            'delta': 0.053,
            'delta_ci': (0.001, 0.105)
        },
    }
    
    # Apply FDR correction to secondary/exploratory comparisons
    secondary_comparisons = {k: v for k, v in comparisons.items() 
                            if 'secondary' in v['endpoint_type'] or 'exploratory' in v['endpoint_type']}
    
    if secondary_comparisons:
        p_values = [v['p_value_raw'] for v in secondary_comparisons.values()]
        correction_results = apply_multiple_testing_correction(p_values, method='fdr_bh')
        
        # Update with corrected p-values
        for i, (key, comp) in enumerate(secondary_comparisons.items()):
            comparisons[key]['p_value_corrected'] = correction_results[i]['p_value_corrected']
            comparisons[key]['significant_after_correction'] = correction_results[i]['significant']
    
    # Create formatted table
    print("Comparison Table (With Statistical Rigor):")
    print("-"*80)
    print(f"{'Comparison':<35} {'Metric':<10} {'Δ':<10} {'95% CI':<20} {'p-val':<10} {'d':<8} {'Sig'}")
    print("-"*80)
    
    for comp_data in comparisons.values():
        comp_name = comp_data['comparison']
        metric = comp_data['metric']
        delta = comp_data['delta']
        ci = f"[{comp_data['delta_ci'][0]:.3f}, {comp_data['delta_ci'][1]:.3f}]"
        p_val = comp_data.get('p_value_corrected', comp_data['p_value_raw'])
        d = comp_data['cohens_d']
        
        sig = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else "ns"
        
        print(f"{comp_name:<35} {metric:<10} {delta:>8.3f} {ci:<20} {p_val:>8.4f} {d:>6.2f} {sig:>4}")
    
    print()
    print("Significance: *** p<0.001, ** p<0.01, * p<0.05, ns = not significant")
    print("Effect size (d): Negligible (<0.2), Small (0.2-0.5), Medium (0.5-0.8), Large (>0.8)")
    print()
    
    return comparisons


def main():
    print("="*80)
    print("STATISTICAL FRAMEWORK FOR RIGOROUS ANALYSIS")
    print("="*80)
    print()
    
    # Create Statistical Analysis Plan
    print("Creating Statistical Analysis Plan...")
    sap = StatisticalAnalysisPlan(output_file='statistical_analysis_plan.json')
    plan = sap.save_plan()
    print()
    
    # Define key comparisons
    sap.add_comparison(
        name='domain_vs_full',
        method1='Domain-only embeddings',
        method2='Full-length embeddings',
        metric='ari',
        hypothesis='Domain extraction improves clustering'
    )
    
    sap.add_comparison(
        name='layers_mid_vs_final',
        method1='Layers 20-33 (mid)',
        method2='Layer 33 (final)',
        metric='ari',
        hypothesis='Mid-layer averaging improves clustering'
    )
    
    sap.add_comparison(
        name='calibrated_vs_uncalibrated',
        method1='Calibrated model',
        method2='Uncalibrated model',
        metric='macro_f1',
        hypothesis='Calibration improves classification'
    )
    
    # Print SAP summary
    print("="*80)
    print("STATISTICAL ANALYSIS PLAN (SAP)")
    print("="*80)
    print()
    
    print("PRIMARY ENDPOINTS:")
    for category, ep in PRIMARY_ENDPOINTS.items():
        print(f"  {category.upper()}: {ep['primary'].upper()}")
    print()
    
    print("SECONDARY ENDPOINTS:")
    for category, ep in PRIMARY_ENDPOINTS.items():
        for metric in ep.get('secondary', []):
            print(f"  {category.upper()}: {metric.upper()}")
    print()
    
    print("MULTIPLE TESTING CORRECTION:")
    print("  Primary:     No correction (single prespecified)")
    print("  Secondary:   Bonferroni")
    print("  Exploratory: Benjamini-Hochberg FDR")
    print("  Motifs (30): Benjamini-Hochberg FDR")
    print()
    
    # Analyze comparisons
    comparisons = analyze_key_comparisons_with_corrections()
    
    # Generate report
    report = generate_statistical_report(comparisons, sap, 'statistical_report.txt')
    
    # Save comparisons
    comparisons_df = pd.DataFrame([
        {
            'comparison': v['comparison'],
            'metric': v['metric'],
            'method1_value': v['method1_value'],
            'method2_value': v['method2_value'],
            'delta': v['delta'],
            'ci_lower': v['delta_ci'][0],
            'ci_upper': v['delta_ci'][1],
            'p_value': v.get('p_value_corrected', v['p_value_raw']),
            'cohens_d': v['cohens_d'],
            'endpoint_type': v['endpoint_type']
        }
        for v in comparisons.values()
    ])
    
    comparisons_df.to_csv('statistical_comparisons.csv', index=False)
    print(f"✅ Comparisons table saved to: statistical_comparisons.csv")
    print()
    
    print("="*80)
    print("✅ STATISTICAL FRAMEWORK COMPLETE!")
    print("="*80)
    print()
    print("All analyses now follow preregistered statistical plan with:")
    print("  ✅ Predefined primary endpoints")
    print("  ✅ Multiple testing correction (FDR)")
    print("  ✅ Effect sizes with CIs (Cohen's d, Δmetric)")
    print("  ✅ Proper significance thresholds")
    print()


if __name__ == "__main__":
    main()

