#!/usr/bin/env python3
"""
download_uniprot_kinases.py - Download kinase sequences from UniProt

This script downloads kinase sequences from UniProt using the REST API,
producing both FASTA and TSV files with full provenance tracking.

Usage:
    python scripts/download_uniprot_kinases.py [--dry-run]

Outputs:
    data/raw/uniprot_kinases.fasta    - Protein sequences
    data/raw/uniprot_kinases.tsv      - Metadata table
    data/raw/uniprot_query.txt        - Exact query used
    data/raw/uniprot_release.txt      - Release version and access date

Author: Kinases-Clustering Project
Date: 2025
"""

import argparse
import json
import os
import sys
from datetime import datetime
from pathlib import Path
from urllib.parse import urlencode
import requests
import time

# ============================================================================
# CONFIGURATION - SINGLE SOURCE OF TRUTH
# ============================================================================

# The exact UniProt query string
UNIPROT_QUERY = "reviewed:true AND (keyword:KW-0418 OR name:kinase*)"

# UniProt REST API base URL
UNIPROT_API_BASE = "https://rest.uniprot.org/uniprotkb"

# Fields to download in TSV format
TSV_FIELDS = [
    "accession",
    "id",                      # Entry name
    "protein_name",
    "gene_names",
    "organism_name",
    "organism_id",
    "length",
    "keyword",
    "ft_domain",               # Domain annotations
    "xref_pfam",               # Pfam cross-references
    "cc_function",             # Function annotation
    "cc_catalytic_activity",   # Catalytic activity
    "sequence",
]

# Output directory
OUTPUT_DIR = Path("data/raw")

# Request settings
REQUEST_TIMEOUT = 300  # 5 minutes
RETRY_ATTEMPTS = 3
RETRY_DELAY = 10  # seconds


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def get_uniprot_release():
    """Get current UniProt release information."""
    try:
        response = requests.head(
            "https://rest.uniprot.org/uniprotkb/search",
            timeout=30
        )
        # Extract release from headers if available
        release = response.headers.get("X-UniProt-Release", "unknown")
        return release
    except Exception as e:
        print(f"Warning: Could not fetch release info: {e}")
        return "unknown"


def download_with_retry(url, params, timeout=REQUEST_TIMEOUT, retries=RETRY_ATTEMPTS):
    """Download with retry logic for robustness."""
    for attempt in range(retries):
        try:
            response = requests.get(url, params=params, timeout=timeout, stream=True)
            response.raise_for_status()
            return response
        except requests.exceptions.RequestException as e:
            if attempt < retries - 1:
                print(f"  Attempt {attempt + 1} failed: {e}. Retrying in {RETRY_DELAY}s...")
                time.sleep(RETRY_DELAY)
            else:
                raise


def count_fasta_sequences(fasta_path):
    """Count sequences in a FASTA file."""
    count = 0
    with open(fasta_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                count += 1
    return count


def count_tsv_rows(tsv_path):
    """Count data rows in a TSV file (excluding header)."""
    with open(tsv_path, 'r') as f:
        lines = f.readlines()
    # Subtract 1 for header
    return len(lines) - 1 if lines else 0


def extract_unique_accessions(tsv_path):
    """Extract unique accession IDs from TSV."""
    accessions = set()
    with open(tsv_path, 'r') as f:
        header = f.readline().strip().split('\t')
        acc_idx = header.index('Entry') if 'Entry' in header else 0
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) > acc_idx:
                accessions.add(parts[acc_idx])
    return accessions


# ============================================================================
# MAIN DOWNLOAD FUNCTIONS
# ============================================================================

def download_fasta(query, output_path, dry_run=False):
    """Download sequences in FASTA format."""
    print(f"\n📥 Downloading FASTA sequences...")
    
    url = f"{UNIPROT_API_BASE}/stream"
    params = {
        "query": query,
        "format": "fasta",
        "compressed": "false",
    }
    
    if dry_run:
        print(f"  [DRY RUN] Would download from: {url}")
        print(f"  [DRY RUN] Query: {query}")
        return 0
    
    response = download_with_retry(url, params)
    
    # Write to file
    with open(output_path, 'wb') as f:
        for chunk in response.iter_content(chunk_size=8192):
            f.write(chunk)
    
    # Count sequences
    n_sequences = count_fasta_sequences(output_path)
    print(f"  ✓ Downloaded {n_sequences:,} sequences to {output_path}")
    return n_sequences


def download_tsv(query, output_path, fields, dry_run=False):
    """Download metadata in TSV format."""
    print(f"\n📥 Downloading TSV metadata...")
    
    url = f"{UNIPROT_API_BASE}/stream"
    params = {
        "query": query,
        "format": "tsv",
        "fields": ",".join(fields),
        "compressed": "false",
    }
    
    if dry_run:
        print(f"  [DRY RUN] Would download from: {url}")
        print(f"  [DRY RUN] Fields: {', '.join(fields)}")
        return 0
    
    response = download_with_retry(url, params)
    
    # Write to file
    with open(output_path, 'wb') as f:
        for chunk in response.iter_content(chunk_size=8192):
            f.write(chunk)
    
    # Count rows
    n_rows = count_tsv_rows(output_path)
    print(f"  ✓ Downloaded {n_rows:,} entries to {output_path}")
    return n_rows


def save_query_file(query, output_path):
    """Save the exact query string for reproducibility."""
    with open(output_path, 'w') as f:
        f.write(f"# UniProt Query String\n")
        f.write(f"# This is the exact query used to download kinase sequences\n")
        f.write(f"# Do not modify - used for reproducibility\n\n")
        f.write(f"{query}\n")
    print(f"  ✓ Saved query to {output_path}")


def save_release_file(release, access_date, output_path, n_fasta, n_tsv, n_unique):
    """Save release and provenance information."""
    with open(output_path, 'w') as f:
        f.write(f"# UniProt Download Provenance\n")
        f.write(f"# Generated by: scripts/download_uniprot_kinases.py\n\n")
        f.write(f"uniprot_release: {release}\n")
        f.write(f"access_date: {access_date}\n")
        f.write(f"download_timestamp: {datetime.now().isoformat()}\n")
        f.write(f"database: SwissProt (reviewed entries only)\n")
        f.write(f"n_fasta_sequences: {n_fasta}\n")
        f.write(f"n_tsv_rows: {n_tsv}\n")
        f.write(f"n_unique_accessions: {n_unique}\n")
        f.write(f"fasta_tsv_match: {n_fasta == n_tsv}\n")
    print(f"  ✓ Saved release info to {output_path}")


def update_provenance_json(release, access_date, n_sequences, n_unique):
    """Update the main provenance.json file."""
    provenance_path = Path("data/provenance.json")
    
    if provenance_path.exists():
        with open(provenance_path, 'r') as f:
            provenance = json.load(f)
    else:
        provenance = {}
    
    # Update UniProt section
    provenance["uniprot"] = {
        "source": "UniProt SwissProt (reviewed entries only)",
        "release": release,
        "access_date": access_date,
        "download_timestamp": datetime.now().isoformat(),
        "query": UNIPROT_QUERY,
        "download_format": "FASTA + TSV",
        "records_retrieved": n_sequences,
        "unique_accessions": n_unique,
        "fields_downloaded": TSV_FIELDS,
        "isoforms": "canonical only (UniProt default)",
        "url": "https://www.uniprot.org/",
        "download_script": "scripts/download_uniprot_kinases.py",
        "raw_files": {
            "fasta": "data/raw/uniprot_kinases.fasta",
            "tsv": "data/raw/uniprot_kinases.tsv",
            "query": "data/raw/uniprot_query.txt",
            "release": "data/raw/uniprot_release.txt"
        }
    }
    
    provenance["updated_at"] = datetime.now().isoformat()
    
    with open(provenance_path, 'w') as f:
        json.dump(provenance, f, indent=2)
    
    print(f"  ✓ Updated {provenance_path}")


# ============================================================================
# MAIN
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Download kinase sequences from UniProt",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python scripts/download_uniprot_kinases.py           # Full download
    python scripts/download_uniprot_kinases.py --dry-run # Preview only
        """
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Preview download without actually fetching data"
    )
    args = parser.parse_args()
    
    print("=" * 70)
    print("UniProt Kinase Download Script")
    print("=" * 70)
    
    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # Define output paths
    fasta_path = OUTPUT_DIR / "uniprot_kinases.fasta"
    tsv_path = OUTPUT_DIR / "uniprot_kinases.tsv"
    query_path = OUTPUT_DIR / "uniprot_query.txt"
    release_path = OUTPUT_DIR / "uniprot_release.txt"
    
    # Check for existing files
    if fasta_path.exists() or tsv_path.exists():
        print("\n⚠️  WARNING: Raw data files already exist!")
        print(f"    {fasta_path}: {'EXISTS' if fasta_path.exists() else 'not found'}")
        print(f"    {tsv_path}: {'EXISTS' if tsv_path.exists() else 'not found'}")
        print("\n    To re-download, manually delete these files first.")
        print("    This prevents accidental overwriting of original data.")
        if not args.dry_run:
            sys.exit(1)
    
    # Print configuration
    print(f"\n📋 Configuration:")
    print(f"    Query: {UNIPROT_QUERY}")
    print(f"    Database: SwissProt (reviewed only)")
    print(f"    Output: {OUTPUT_DIR}/")
    
    # Get release info
    print(f"\n🔍 Checking UniProt release...")
    release = get_uniprot_release()
    access_date = datetime.now().strftime("%Y-%m-%d")
    print(f"    Release: {release}")
    print(f"    Access date: {access_date}")
    
    # Download FASTA
    n_fasta = download_fasta(UNIPROT_QUERY, fasta_path, dry_run=args.dry_run)
    
    # Download TSV
    n_tsv = download_tsv(UNIPROT_QUERY, tsv_path, TSV_FIELDS, dry_run=args.dry_run)
    
    if not args.dry_run:
        # Save query file
        print(f"\n📝 Saving provenance files...")
        save_query_file(UNIPROT_QUERY, query_path)
        
        # Count unique accessions
        unique_accessions = extract_unique_accessions(tsv_path)
        n_unique = len(unique_accessions)
        
        # Save release file
        save_release_file(release, access_date, release_path, n_fasta, n_tsv, n_unique)
        
        # Update provenance.json
        update_provenance_json(release, access_date, n_fasta, n_unique)
        
        # Sanity checks
        print(f"\n✅ Sanity Checks:")
        print(f"    FASTA sequences: {n_fasta:,}")
        print(f"    TSV rows: {n_tsv:,}")
        print(f"    Unique accessions: {n_unique:,}")
        
        if n_fasta == n_tsv:
            print(f"    ✓ FASTA count matches TSV count")
        else:
            print(f"    ⚠️  MISMATCH: FASTA ({n_fasta}) != TSV ({n_tsv})")
            print(f"       This may indicate isoforms or download issues")
        
        if n_fasta == n_unique:
            print(f"    ✓ All accessions are unique")
        else:
            print(f"    ⚠️  Duplicate accessions detected: {n_fasta - n_unique}")
        
        print(f"\n🎉 Download complete!")
        print(f"\n📁 Output files:")
        print(f"    {fasta_path}")
        print(f"    {tsv_path}")
        print(f"    {query_path}")
        print(f"    {release_path}")
        
        # Summary for manuscript
        print(f"\n📝 For MANUSCRIPT.md:")
        print(f"    N_raw_sequences = {n_fasta:,}")
        print(f"    N_raw_unique_accessions = {n_unique:,}")
    else:
        print(f"\n[DRY RUN] No files were downloaded.")


if __name__ == "__main__":
    main()

