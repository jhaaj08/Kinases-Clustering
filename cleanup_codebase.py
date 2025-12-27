#!/usr/bin/env python3
"""
Codebase Cleanup Script

This script identifies and removes unnecessary files from the repository,
keeping only what's needed for the reproducible pipeline.

Usage:
    python cleanup_codebase.py --dry-run    # Preview what would be deleted
    python cleanup_codebase.py --execute    # Actually delete files
    python cleanup_codebase.py --archive    # Move to _archive instead of delete

Author: Pipeline
"""

import argparse
import os
import shutil
from pathlib import Path
from datetime import datetime

# Files and directories to KEEP
KEEP = {
    # Core pipeline
    "pipeline/",
    
    # Build system
    "Makefile",
    "requirements.txt",
    "environment.yml",
    ".gitignore",
    "LICENSE",
    
    # Documentation
    "README.md",
    "MANUSCRIPT.md",
    "docs/Simple_English.md",
    "CLEANUP_PLAN.md",
    "cleanup_codebase.py",
    
    # Source data
    "data/domains/domains_E001.fasta",
    "data/domains/domain_coords_E001.tsv",
    "data/domains/domain_extraction_report.json",
    "data/domains/hmmer_domtblout_E001.txt",
    "data/processed/labels.csv",
    "data/processed/label_policy.json",
    "data/hmm_profiles/",
    
    # Optional data (regenerated but useful to keep)
    "data/manifests/",
    "data/splits/",
    "data/raw/kinases_all.csv",
    "data/raw/kinases_revised.csv",
    
    # Pre-computed embeddings (CRITICAL)
    "embeddings/esm2_t33_650M/",
    
    # Runs directory structure
    "runs/.gitkeep",
    "runs/2025-12-25_110013/",  # Keep one reference run
    
    # Verification
    "scripts/verify_package.py",
    
    # Web app
    "webapp/",
    
    # Git
    ".git/",
    ".github/",
}

# Files and directories to DELETE
DELETE = [
    # Old embedding directories (DUPLICATES)
    "kinases_domains_e0.01_embeddings/",
    "kinases_domains_e0.01_layers_20_30/",
    "kinases_domains_e0.01_layers_mid/",
    "kinases_domains_e0.01_cls/",
    "kinases_domains_e0.1_embeddings/",
    "kinases_domains_embeddings/",
    "kinases_embeddings/",
    
    # Old results directories
    "clustering/",
    "supervised_results/",
    "supervised_results_calibrated/",
    "supervised_results_homology/",
    "supervised_results_layer_comparison/",
    "exemplar_retrieval_results/",
    "calibration_comparison_results/",
    "clustering_statistics/",
    "results/",
    "models/",
    "reports/",
    
    # Old figures
    "figures/",
    "figures_output/",
    
    # Archive & legacy
    "archive/",
    
    # Zenodo (already uploaded)
    "zenodo_package/",
    
    # Logs & temp
    "logs/",
    "data/temp_cdhit/",
    "__pycache__/",
    "pipeline/__pycache__/",
    "webapp/__pycache__/",
    
    # Old root-level scripts
    "generate_figures.py",
    "generate_esm2_embeddings_v3.py",
    "extract_kinase_domains_v2.py",
    "extract_motif_features.py",
    "mutation_motif_analysis.py",
    "statistical_framework.py",
    "train_supervised_enhanced.py",
    "baselines_comparison.py",
    
    # Deployment files
    "Dockerfile",
    ".dockerignore",
    "deploy-to-gcp.sh",
    "deploy-to-gcp-v2.sh",
    "DEPLOYMENT.md",
    "DEPLOYMENT_SUMMARY.txt",
    "SETUP_GCP.md",
    
    # Misc files
    "Snakefile",
    "REPOSITORY_STRUCTURE.md",
    "START_HERE.md",
    "READY_TO_DEPLOY.txt",
    "REQUIREMENTS_FROM_USER.md",
    "pyproject.toml",
    "CITATION.cff",
    "configs/",
    "tests/",
    "utils/",
    "src/",
    
    # Old data files
    "data/splits.json",
    "data/splits_40.json",
    "data/splits_50.json",
    "data/splits_70.json",
    "data/provenance.json",
    "data/README_ZENODO.txt",
]

# Old scripts to delete (except verify_package.py)
OLD_SCRIPTS = [
    "scripts/01_dataset_preparation.py",
    "scripts/02_extract_domains.py",
    "scripts/03_extract_embeddings.py",
    "scripts/04_extract_motifs.py",
    "scripts/05_clustering_eval.py",
    "scripts/06_supervised_classification.py",
    "scripts/07_08_09_10_11_12_wrappers.sh",
    "scripts/13_prepare_zenodo_upload.py",
    "scripts/run_all_pipeline.sh",
    "scripts/README.md",
    "scripts/assign_labels.py",
    "scripts/build_manuscript_numbers.py",
    "scripts/calibrate_model.py",
    "scripts/compare_layer_classification.py",
    "scripts/create_dataset_manifests.py",
    "scripts/create_homology_splits.py",
    "scripts/deduplicate_sequences.py",
    "scripts/download_uniprot_kinases.py",
    "scripts/extract_domains.py",
    "scripts/filter_sequences.py",
    "scripts/generate_embeddings.py",
    "scripts/generate_raw_manifest.py",
    "scripts/reduce_redundancy_cdhit.py",
    "scripts/regenerate_embeddings.py",
    "scripts/regenerate_splits.py",
    "scripts/run_baselines.py",
    "scripts/run_calibration_comparison.py",
    "scripts/run_clustering.py",
    "scripts/run_homology_generalization.py",
    "scripts/run_layer_supervised_comparison.py",
    "scripts/run_retrieval.py",
    "scripts/sync_manuscript_numbers.py",
    "scripts/train_supervised.py",
]

# Old docs to delete (except Simple_English.md)
OLD_DOCS = [
    "docs/LABEL_RECOVERY_REPORT.md",
    "docs/VALIDATION_REPORT_TEMPLATE.md",
    "docs/PROJECT_STATUS.md",
    "docs/PUBLICATION_SUBMISSION_PACKAGE.md",
    "docs/FINAL_RESULTS_SUMMARY.md",
    "docs/ZENODO_SUCCESS.md",
    "docs/PROVENANCE_IMPLEMENTATION.md",
    "docs/PROJECT_FINALIZATION_SUCCESS.txt",
    "docs/FINAL_SUBMISSION_STATUS.md",
    "docs/PROJECT_COMPLETION_SUMMARY.md",
    "docs/README_FINALIZATION.md",
    "docs/EXTERNAL_DATA_INTEGRATION.md",
    "docs/REVIEWER_REQUIREMENTS_CHECKLIST.md",
    "docs/COMPLETION_STATUS.md",
    "docs/EXECUTIVE_SUMMARY.md",
    "docs/FINAL_COMPLETION_REPORT.md",
    "docs/ZENODO_UPLOAD_STATUS.md",
    "docs/EXPERIMENTAL_SUMMARY.md",
    "docs/REPOSITORY_CLEANUP_SUMMARY.md",
    "docs/VALIDATION_REPORT.md",
    "docs/JOURNAL_COMPLIANCE_CHECKLIST.md",
    "docs/FINAL_SUBMISSION_CHECKLIST.md",
    "docs/EMBEDDING_METHODOLOGY.md",
    "docs/DETERMINISM.md",
    "docs/DEPLOYMENT_GUIDE.md",
]

# ZIP files to delete
ZIP_FILES = [
    "kinase_data_2025-12-25_110013.zip",
]


def get_size(path):
    """Get size of file or directory in bytes."""
    if os.path.isfile(path):
        return os.path.getsize(path)
    total = 0
    for root, dirs, files in os.walk(path):
        for f in files:
            fp = os.path.join(root, f)
            try:
                total += os.path.getsize(fp)
            except OSError:
                pass
    return total


def format_size(bytes):
    """Format bytes to human-readable string."""
    for unit in ['B', 'KB', 'MB', 'GB']:
        if bytes < 1024:
            return f"{bytes:.1f} {unit}"
        bytes /= 1024
    return f"{bytes:.1f} TB"


def main():
    parser = argparse.ArgumentParser(description='Cleanup codebase')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--dry-run', action='store_true', help='Preview what would be deleted')
    group.add_argument('--execute', action='store_true', help='Actually delete files')
    group.add_argument('--archive', action='store_true', help='Move to _archive instead of delete')
    args = parser.parse_args()
    
    print("=" * 70)
    print("CODEBASE CLEANUP")
    print("=" * 70)
    
    # Collect all items to delete
    all_delete = DELETE + OLD_SCRIPTS + OLD_DOCS + ZIP_FILES
    
    # Count and calculate sizes
    to_delete = []
    total_size = 0
    
    for item in all_delete:
        if os.path.exists(item):
            size = get_size(item)
            total_size += size
            to_delete.append((item, size))
    
    print(f"\nItems to {'archive' if args.archive else 'delete'}: {len(to_delete)}")
    print(f"Total size: {format_size(total_size)}")
    
    print("\n" + "-" * 70)
    print("FILES/DIRECTORIES:")
    print("-" * 70)
    
    for item, size in sorted(to_delete, key=lambda x: -x[1])[:30]:
        print(f"  {item:<50} {format_size(size):>10}")
    
    if len(to_delete) > 30:
        print(f"  ... and {len(to_delete) - 30} more items")
    
    if args.dry_run:
        print("\n" + "=" * 70)
        print("DRY RUN - No files were deleted")
        print("Run with --execute to actually delete, or --archive to move to _archive/")
        print("=" * 70)
        return
    
    if args.archive:
        archive_dir = Path("_archive_before_cleanup")
        archive_dir.mkdir(exist_ok=True)
        print(f"\nArchiving to: {archive_dir}/")
    
    # Confirm
    if args.execute:
        print("\n" + "!" * 70)
        print("WARNING: This will PERMANENTLY DELETE the files listed above!")
        print("!" * 70)
        confirm = input("Type 'DELETE' to confirm: ")
        if confirm != 'DELETE':
            print("Aborted.")
            return
    
    # Execute
    deleted_count = 0
    deleted_size = 0
    
    for item, size in to_delete:
        try:
            if args.archive:
                dest = Path("_archive_before_cleanup") / item
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.move(item, dest)
                print(f"  Archived: {item}")
            else:
                if os.path.isdir(item):
                    shutil.rmtree(item)
                else:
                    os.remove(item)
                print(f"  Deleted: {item}")
            deleted_count += 1
            deleted_size += size
        except Exception as e:
            print(f"  Error: {item} - {e}")
    
    print("\n" + "=" * 70)
    print(f"CLEANUP COMPLETE")
    print(f"  Items {'archived' if args.archive else 'deleted'}: {deleted_count}")
    print(f"  Space freed: {format_size(deleted_size)}")
    print("=" * 70)
    
    # Update .gitignore if archiving
    if args.archive:
        with open(".gitignore", "a") as f:
            f.write("\n# Archive from cleanup\n_archive_before_cleanup/\n")
        print("\nAdded _archive_before_cleanup/ to .gitignore")
    
    print("\nNext steps:")
    print("  1. Run: make all RUN_ID=post_cleanup_test")
    print("  2. Run: make verify")
    print("  3. Commit: git add -A && git commit -m 'Cleanup: remove old files'")


if __name__ == "__main__":
    main()

