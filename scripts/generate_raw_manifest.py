#!/usr/bin/env python3
"""
Generate Raw Data Manifest with SHA-256 Hashes

This script creates a manifest file documenting all raw data files
with their checksums, sizes, and line counts. This is essential for:
1. Verifying data integrity after transfer
2. Zenodo upload preparation
3. Reproducibility verification

Usage:
    python scripts/generate_raw_manifest.py

Output:
    data/raw/MANIFEST.txt - Human-readable manifest
    data/raw/MANIFEST.json - Machine-readable manifest
"""

import hashlib
import json
from pathlib import Path
from datetime import datetime
import os

# Raw files to include in manifest
RAW_FILES = [
    "data/raw/kinases_all.csv",
    "data/raw/kinases_revised.csv", 
    "data/raw/uniprot_query.txt",
    "data/raw/uniprot_release.txt",
]

# Additional key processed files for the data package
PROCESSED_FILES = [
    "data/processed/labels.csv",
    "data/processed/label_policy.json",
    "data/processed/dataset_manifest_report.json",
    "data/domains/domains_E001.fasta",
    "data/domains/domain_coords_E001.tsv",
    "data/domains/domain_extraction_report.json",
    "data/domains/hmmer_domtblout_E001.txt",
    "data/manifests/supervised_eligible.txt",
    "data/manifests/domain_E001.txt",
    "data/splits/split40_train.txt",
    "data/splits/split40_test.txt",
    "data/splits/split50_train.txt",
    "data/splits/split50_test.txt",
    "data/splits/split70_train.txt",
    "data/splits/split70_test.txt",
    "data/splits/splits_report.json",
]

# HMM profiles (important for reproducibility)
HMM_FILES = [
    "data/hmm_profiles/PF00069.hmm",
    "data/hmm_profiles/PF07714.hmm",
]

# Results registries (SOURCE OF TRUTH - prevents hand-typed drift)
RESULTS_FILES = [
    "results/manuscript_numbers.json",
    "results/tables/Table1.csv",
    "results/tables/TableS1.csv",
    "results/tables/TableS2.csv",
    "results/clustering/clustering_registry.json",
    "results/clustering/summary_table.csv",
    "results/supervised/lr_split40_metrics.json",
    "results/supervised/lr_split50_metrics.json",
    "results/supervised/lr_split70_metrics.json",
    "results/supervised/lr_multi_identity_summary.csv",
    "results/supervised/supervised_registry.json",
    "results/calibration/split40_calibration.json",
    "results/baselines/baselines_split40.csv",
    "results/baselines/knn_split40.json",
    "results/baselines/mlp_split40.json",
    "results/baselines/motifs_split40.json",
    "results/baselines/random_split40.json",
    "results/retrieval/split40_retrieval.json",
    "results/retrieval/summary.csv",
    "results/layer_comparison/layer_comparison_results.json",
    "results/layer_comparison/layer_comparison_summary.csv",
]

# Embeddings metadata (not the actual embeddings - too large)
EMBEDDINGS_META = [
    "embeddings/esm2_t33_650M/embedding_metadata.json",
    "embeddings/esm2_t33_650M/ids.txt",
]


def compute_sha256(filepath):
    """Compute SHA-256 hash of a file."""
    sha256_hash = hashlib.sha256()
    with open(filepath, "rb") as f:
        for byte_block in iter(lambda: f.read(4096), b""):
            sha256_hash.update(byte_block)
    return sha256_hash.hexdigest()


def count_lines(filepath):
    """Count lines in a file."""
    try:
        with open(filepath, 'r', errors='ignore') as f:
            return sum(1 for _ in f)
    except:
        return -1


def get_file_info(filepath):
    """Get comprehensive file information."""
    path = Path(filepath)
    if not path.exists():
        return {
            "path": filepath,
            "exists": False,
            "status": "MISSING"
        }
    
    stat = path.stat()
    return {
        "path": filepath,
        "exists": True,
        "size_bytes": stat.st_size,
        "size_human": format_size(stat.st_size),
        "lines": count_lines(filepath),
        "sha256": compute_sha256(filepath),
        "modified": datetime.fromtimestamp(stat.st_mtime).isoformat(),
        "status": "OK"
    }


def format_size(size_bytes):
    """Format bytes to human readable string."""
    for unit in ['B', 'KB', 'MB', 'GB']:
        if size_bytes < 1024:
            return f"{size_bytes:.1f} {unit}"
        size_bytes /= 1024
    return f"{size_bytes:.1f} TB"


def main():
    print("=" * 70)
    print("RAW DATA MANIFEST GENERATOR")
    print("=" * 70)
    
    manifest = {
        "generated": datetime.now().isoformat(),
        "generator": "scripts/generate_raw_manifest.py",
        "purpose": "Zenodo upload and reproducibility verification",
        "categories": {}
    }
    
    all_files = []
    
    # Process each category
    categories = [
        ("raw_data", RAW_FILES, "Core raw data files from UniProt"),
        ("processed_key", PROCESSED_FILES, "Key processed files for analysis"),
        ("hmm_profiles", HMM_FILES, "Pfam HMM profiles for domain extraction"),
        ("results_registries", RESULTS_FILES, "Results registries (SOURCE OF TRUTH)"),
        ("embeddings_meta", EMBEDDINGS_META, "Embedding metadata and IDs"),
    ]
    
    for category_name, file_list, description in categories:
        print(f"\n### {category_name.upper()} ###")
        print(f"Description: {description}")
        print("-" * 60)
        
        category_files = []
        for filepath in file_list:
            info = get_file_info(filepath)
            category_files.append(info)
            all_files.append(info)
            
            if info["exists"]:
                print(f"✓ {info['path']}")
                print(f"    Size: {info['size_human']} | Lines: {info['lines']}")
                print(f"    SHA256: {info['sha256'][:16]}...")
            else:
                print(f"✗ {info['path']} - MISSING")
        
        manifest["categories"][category_name] = {
            "description": description,
            "files": category_files
        }
    
    # Summary statistics
    total_size = sum(f.get("size_bytes", 0) for f in all_files if f["exists"])
    total_files = len([f for f in all_files if f["exists"]])
    missing_files = len([f for f in all_files if not f["exists"]])
    
    manifest["summary"] = {
        "total_files": total_files,
        "missing_files": missing_files,
        "total_size_bytes": total_size,
        "total_size_human": format_size(total_size)
    }
    
    # Write JSON manifest
    json_path = Path("data/raw/MANIFEST.json")
    with open(json_path, 'w') as f:
        json.dump(manifest, f, indent=2)
    print(f"\n✓ JSON manifest saved: {json_path}")
    
    # Write human-readable manifest
    txt_path = Path("data/raw/MANIFEST.txt")
    with open(txt_path, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("KINASE CLASSIFICATION DATA PACKAGE - FILE MANIFEST\n")
        f.write("=" * 70 + "\n\n")
        f.write(f"Generated: {manifest['generated']}\n")
        f.write(f"Purpose: Zenodo upload and reproducibility verification\n\n")
        
        f.write("-" * 70 + "\n")
        f.write("SUMMARY\n")
        f.write("-" * 70 + "\n")
        f.write(f"Total files: {total_files}\n")
        f.write(f"Missing files: {missing_files}\n")
        f.write(f"Total size: {format_size(total_size)}\n\n")
        
        for category_name, cat_data in manifest["categories"].items():
            f.write("-" * 70 + "\n")
            f.write(f"{category_name.upper()}: {cat_data['description']}\n")
            f.write("-" * 70 + "\n")
            for file_info in cat_data["files"]:
                if file_info["exists"]:
                    f.write(f"\nFile: {file_info['path']}\n")
                    f.write(f"  Size: {file_info['size_human']}\n")
                    f.write(f"  Lines: {file_info['lines']}\n")
                    f.write(f"  SHA256: {file_info['sha256']}\n")
                else:
                    f.write(f"\nFile: {file_info['path']} - MISSING\n")
            f.write("\n")
        
        f.write("=" * 70 + "\n")
        f.write("END OF MANIFEST\n")
        f.write("=" * 70 + "\n")
    
    print(f"✓ Text manifest saved: {txt_path}")
    
    # Print final summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Total files tracked: {total_files}")
    print(f"Missing files: {missing_files}")
    print(f"Total size: {format_size(total_size)}")
    
    if missing_files > 0:
        print("\n⚠️  WARNING: Some files are missing!")
        print("   Run the data processing pipeline first.")
    else:
        print("\n✓ All files present and hashed.")
        print("✓ Ready for Zenodo upload!")
    
    print("\nNext steps for Zenodo:")
    print("  1. Create a ZIP: zip -r kinase_data_v1.zip data/raw/ data/processed/ data/domains/ data/manifests/ data/splits/ data/hmm_profiles/")
    print("  2. Go to https://zenodo.org and create a new upload")
    print("  3. Upload the ZIP file")
    print("  4. Add the DOI to your MANUSCRIPT.md")


if __name__ == "__main__":
    main()

