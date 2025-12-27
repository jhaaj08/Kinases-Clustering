#!/usr/bin/env python3
"""
Generate MANIFEST.txt with SHA256 hashes for all files in the run.

This script creates a manifest file with SHA256 hashes for verification
of the data package integrity (e.g., for Zenodo upload).

Usage:
    python pipeline/generate_manifest.py --run-dir runs/2025-01-01_000000/

Outputs:
    - MANIFEST.txt
    - MANIFEST.json
"""

import argparse
import hashlib
import json
from pathlib import Path
from datetime import datetime


def compute_sha256(filepath):
    """Compute SHA256 hash of a file."""
    sha256 = hashlib.sha256()
    with open(filepath, 'rb') as f:
        for chunk in iter(lambda: f.read(8192), b''):
            sha256.update(chunk)
    return sha256.hexdigest()


def main():
    parser = argparse.ArgumentParser(description="Generate MANIFEST.txt")
    parser.add_argument("--run-dir", type=str, required=True, help="Run directory path")
    args = parser.parse_args()
    
    run_dir = Path(args.run_dir)
    
    print("=" * 60)
    print("Generating MANIFEST.txt")
    print("=" * 60)
    
    # Define important files to track
    important_patterns = [
        "run_config.json",
        "data/manifests/*.txt",
        "data/manifests/*.json",
        "data/splits/*.txt",
        "data/splits/*.json",
        "embeddings/esm2_t33_650M/ids.txt",
        "embeddings/esm2_t33_650M/embedding_metadata.json",
        "results/**/*.json",
        "results/**/*.csv",
        "tables/*.csv",
    ]
    
    files_to_hash = []
    
    # Collect files
    for pattern in important_patterns:
        if '*' in pattern:
            # Glob pattern
            for f in run_dir.glob(pattern):
                if f.is_file() and not f.is_symlink():
                    files_to_hash.append(f)
        else:
            # Exact file
            f = run_dir / pattern
            if f.is_file() and not f.is_symlink():
                files_to_hash.append(f)
    
    # Sort for consistency
    files_to_hash = sorted(set(files_to_hash))
    
    # Compute hashes
    manifest_data = {
        "generated_at": datetime.now().isoformat(),
        "run_dir": str(run_dir),
        "total_files": len(files_to_hash),
        "files": []
    }
    
    manifest_lines = [
        f"# MANIFEST for {run_dir.name}",
        f"# Generated: {datetime.now().isoformat()}",
        f"# Files: {len(files_to_hash)}",
        "#",
        "# Format: SHA256  <relative_path>",
        "#",
    ]
    
    print(f"\nHashing {len(files_to_hash)} files...")
    
    for filepath in files_to_hash:
        rel_path = filepath.relative_to(run_dir)
        file_hash = compute_sha256(filepath)
        file_size = filepath.stat().st_size
        
        manifest_data["files"].append({
            "path": str(rel_path),
            "sha256": file_hash,
            "size_bytes": file_size
        })
        
        manifest_lines.append(f"{file_hash}  {rel_path}")
        print(f"  ✓ {rel_path}")
    
    # Write MANIFEST.txt
    manifest_txt = run_dir / "MANIFEST.txt"
    with open(manifest_txt, 'w') as f:
        f.write('\n'.join(manifest_lines) + '\n')
    print(f"\n✓ Saved: {manifest_txt}")
    
    # Write MANIFEST.json
    manifest_json = run_dir / "MANIFEST.json"
    with open(manifest_json, 'w') as f:
        json.dump(manifest_data, f, indent=2)
    print(f"✓ Saved: {manifest_json}")
    
    print("\n" + "=" * 60)
    print("MANIFEST GENERATION COMPLETE")
    print("=" * 60)
    print(f"\nTotal files: {len(files_to_hash)}")
    print(f"Output: {manifest_txt}")


if __name__ == "__main__":
    main()

