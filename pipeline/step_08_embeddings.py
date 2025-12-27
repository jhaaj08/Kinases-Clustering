#!/usr/bin/env python3
"""
Step 8: Link/Copy Embeddings to Run Directory

This script links or copies pre-computed embeddings to the run directory.
Since embeddings are expensive to compute, we reuse them from the project-level
embeddings/ directory.

Usage:
    python pipeline/step_08_embeddings.py --run-dir runs/2025-01-01_000000/

Outputs:
    - embeddings/esm2_t33_650M/ids.txt (copied)
    - embeddings/esm2_t33_650M/*.npy (symlinked)
    - embeddings/esm2_t33_650M/embedding_metadata.json (copied)
"""

import argparse
import json
import shutil
from pathlib import Path
from datetime import datetime

# Project root
PROJECT_ROOT = Path(__file__).parent.parent


def main():
    parser = argparse.ArgumentParser(description="Link embeddings to run directory")
    parser.add_argument("--run-dir", type=str, required=True, help="Run directory path")
    parser.add_argument("--copy", action="store_true", help="Copy instead of symlink")
    args = parser.parse_args()
    
    run_dir = Path(args.run_dir)
    run_emb_dir = run_dir / "embeddings" / "esm2_t33_650M"
    run_emb_dir.mkdir(parents=True, exist_ok=True)
    
    src_emb_dir = PROJECT_ROOT / "embeddings" / "esm2_t33_650M"
    
    print("=" * 60)
    print("Step 8: Link Embeddings")
    print("=" * 60)
    
    if not src_emb_dir.exists():
        raise FileNotFoundError(f"Source embeddings not found: {src_emb_dir}")
    
    # Files to copy (small text files)
    copy_files = ["ids.txt", "embedding_metadata.json"]
    
    # Files to symlink (large numpy files)
    npy_files = list(src_emb_dir.glob("*.npy"))
    
    print(f"\nSource: {src_emb_dir}")
    print(f"Destination: {run_emb_dir}")
    
    # Copy small files
    for fname in copy_files:
        src = src_emb_dir / fname
        dst = run_emb_dir / fname
        if src.exists():
            shutil.copy2(src, dst)
            print(f"  ✓ Copied: {fname}")
        else:
            print(f"  ⚠ Not found: {fname}")
    
    # Symlink or copy large files
    for src in npy_files:
        dst = run_emb_dir / src.name
        if dst.exists():
            dst.unlink()
        
        if args.copy:
            shutil.copy2(src, dst)
            print(f"  ✓ Copied: {src.name}")
        else:
            # Create relative symlink for portability
            rel_path = Path("../../../..") / src.relative_to(PROJECT_ROOT)
            dst.symlink_to(rel_path)
            print(f"  ✓ Linked: {src.name}")
    
    # Verify ids.txt
    ids_file = run_emb_dir / "ids.txt"
    if ids_file.exists():
        with open(ids_file) as f:
            ids = [line.strip() for line in f if line.strip()]
        print(f"\n  IDs count: {len(ids)}")
    
    # Update metadata with run info
    meta_file = run_emb_dir / "embedding_metadata.json"
    if meta_file.exists():
        with open(meta_file) as f:
            metadata = json.load(f)
        
        metadata["run_info"] = {
            "linked_at": datetime.now().isoformat(),
            "source_dir": str(src_emb_dir),
            "method": "copy" if args.copy else "symlink"
        }
        
        with open(meta_file, 'w') as f:
            json.dump(metadata, f, indent=2)
    
    print("\n" + "=" * 60)
    print("Step 8 COMPLETE")
    print("=" * 60)
    print(f"\nEmbeddings ready in: {run_emb_dir}")
    print(f"  NPY files: {len(npy_files)}")
    print(f"  IDs: {len(ids) if ids_file.exists() else 'N/A'}")


if __name__ == "__main__":
    main()

