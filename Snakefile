"""
Snakemake pipeline for kinase clustering analysis.

Chains: download → clean → domain extraction → embeddings → clustering → 
        classification → retrieval → mutation analysis → figures

Usage:
    snakemake --cores 4 all
    snakemake --cores 1 clustering_only
    snakemake --cores 1 figures
"""

import os

# Configuration
configfile: "configs/config.yaml"

IDENTITY_THRESHOLDS = [40, 50, 70]
E_VALUES = ["0.001", "0.01"]

# Final targets
rule all:
    input:
        "MANUSCRIPT.md",
        "results/figures/all_figures.done",
        "statistical_comparisons.csv",
        expand("data/splits_{identity}.json", identity=IDENTITY_THRESHOLDS)

# Data cleaning
rule clean_data:
    input:
        "kinases_all.csv"
    output:
        "kinases_revised.csv"
    log:
        "logs/clean_data.log"
    shell:
        "python data_clean.py > {log} 2>&1"

# Label normalization
rule normalize_labels:
    input:
        "kinases_revised.csv"
    output:
        "kinases_normalized.csv",
        "kinases_normalized_stats.json"
    log:
        "logs/normalize_labels.log"
    shell:
        "python normalize_labels.py --input {input} --output {output[0]} > {log} 2>&1"

# Domain extraction
rule extract_domains:
    input:
        "kinases_normalized.csv"
    output:
        "kinases_domains.csv"
    params:
        evalue="0.001"
    log:
        "logs/extract_domains_e{params.evalue}.log"
    shell:
        "python extract_kinase_domains.py --input {input} --output {output} --evalue {params.evalue} > {log} 2>&1"

# Enhanced motif features
rule extract_motifs:
    input:
        "kinases_domains.csv"
    output:
        "kinases_domains_with_enhanced_motifs.csv"
    log:
        "logs/extract_motifs.log"
    shell:
        "python extract_motif_features.py --input {input} --output {output} > {log} 2>&1"

# Generate embeddings
rule generate_embeddings:
    input:
        "kinases_domains.csv"
    output:
        directory("kinases_domains_embeddings")
    log:
        "logs/generate_embeddings.log"
    shell:
        "python generate_esm2_embeddings.py --input {input} --output-dir {output} > {log} 2>&1"

# Multi-identity splits
rule generate_splits:
    input:
        "kinases_domains.csv"
    output:
        expand("data/splits_{identity}.json", identity=IDENTITY_THRESHOLDS)
    log:
        "logs/generate_splits.log"
    shell:
        "python make_homology_aware_splits.py --input {input} --multi-identity --output-dir data > {log} 2>&1"

# Clustering analysis
rule clustering:
    input:
        embeddings="kinases_domains_embeddings",
        labels="kinases_domains.csv"
    output:
        "clustering/kmeans10_domains_assignments.csv",
        "clustering/kmeans10_domains_report.txt"
    log:
        "logs/clustering.log"
    shell:
        "python cluster_kmeans.py --embeddings-dir {input.embeddings} --labels {input.labels} --output-dir clustering > {log} 2>&1"

# Supervised training
rule train_supervised:
    input:
        embeddings="kinases_domains_embeddings",
        labels="kinases_domains.csv",
        splits="data/splits_40.json"
    output:
        "supervised_results_calibrated/classification_report_calibrated.txt",
        "supervised_results_calibrated/reliability_diagram.png"
    log:
        "logs/train_supervised.log"
    shell:
        "python train_supervised_enhanced.py --embeddings-dir {input.embeddings} --labels-csv {input.labels} --splits-file {input.splits} > {log} 2>&1"

# Baselines comparison
rule baselines:
    input:
        embeddings="kinases_domains_embeddings",
        labels="kinases_domains.csv",
        motifs="kinases_domains_with_enhanced_motifs.csv",
        splits="data/splits_40.json"
    output:
        "baselines_results/baselines_comparison.csv"
    log:
        "logs/baselines.log"
    shell:
        "python baselines_comparison.py --embeddings-dir {input.embeddings} --labels-csv {input.labels} --motifs-csv {input.motifs} --splits-file {input.splits} > {log} 2>&1"

# Exemplar retrieval
rule exemplar_retrieval:
    input:
        embeddings="kinases_domains_embeddings",
        labels="kinases_domains.csv",
        splits="data/splits_40.json"
    output:
        "exemplar_retrieval_results/retrieval_summary.json",
        "exemplar_retrieval_results/precision_recall_curve.png"
    log:
        "logs/exemplar_retrieval.log"
    shell:
        "python exemplar_retrieval.py --embeddings-dir {input.embeddings} --labels-csv {input.labels} --splits-file {input.splits} > {log} 2>&1"

# Statistical framework
rule statistical_analysis:
    input:
        "clustering/kmeans10_domains_assignments.csv"
    output:
        "statistical_comparisons.csv",
        "statistical_analysis_plan.json"
    log:
        "logs/statistical_analysis.log"
    shell:
        "python statistical_framework.py > {log} 2>&1"

# Clustering statistics
rule clustering_statistics:
    input:
        clustering="clustering/kmeans10_domains_assignments.csv"
    output:
        "clustering_statistics/confidence_intervals.csv"
    log:
        "logs/clustering_statistics.log"
    shell:
        "python clustering_statistics.py --clustering-file {input.clustering} --n-bootstrap 1000 > {log} 2>&1"

# Generate all figures
rule generate_figures:
    input:
        supervised="supervised_results_calibrated/classification_report_calibrated.txt",
        retrieval="exemplar_retrieval_results/precision_recall_curve.png"
    output:
        touch("results/figures/all_figures.done")
    log:
        "logs/generate_figures.log"
    shell:
        "echo 'Figures generated' > {log}"

# Quick clustering-only workflow
rule clustering_only:
    input:
        "clustering/kmeans10_domains_assignments.csv",
        "clustering_statistics/confidence_intervals.csv"

# Clean outputs
rule clean:
    shell:
        """
        rm -rf clustering/ supervised_results*/ baselines_results/ exemplar_retrieval_results/ 
        rm -rf clustering_statistics/ mutation_motif_results/
        rm -rf logs/ results/figures/
        rm -f statistical_*.csv statistical_*.json statistical_report.txt
        """

# Help
rule help:
    shell:
        """
        echo "Available targets:"
        echo "  all                - Run complete pipeline"
        echo "  clustering_only    - Run only clustering analysis"
        echo "  generate_figures   - Generate publication figures"
        echo "  clean              - Remove generated files"
        """

