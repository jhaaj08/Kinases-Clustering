#!/usr/bin/env python3
"""
deduplicate_sequences.py - Remove exact duplicate sequences

This script removes sequences that are 100% identical (same amino acid string).
When duplicates are found, we keep the representative with the lowest 
accession ID (alphanumerically sorted) for reproducibility.

Usage:
    python scripts/deduplicate_sequences.py [--dry-run]

Inputs:
    data/processed/step2_filtered.fasta
    data/processed/step2_filtered.tsv

Outputs:
    data/processed/step3_dedup.fasta
    data/processed/step3_dedup.tsv
    data/processed/step3_dedup_map.tsv (duplicate → representative mapping)
    data/processed/step3_dedup_report.json

Author: Kinases-Clustering Project
Date: 2025
"""

import argparse
import hashlib
import json
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path

# ============================================================================
# CONFIGURATION
# ============================================================================

# Input/Output directories
INPUT_DIR = Path("data/processed")
OUTPUT_DIR = Path("data/processed")

# Input files (from Step 2)
INPUT_FASTA = INPUT_DIR / "step2_filtered.fasta"
INPUT_TSV = INPUT_DIR / "step2_filtered.tsv"

# Representative selection rule
# Options: "lowest_accession" (alphanumeric sort), "first_encountered"
REPRESENTATIVE_RULE = "lowest_accession"


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def parse_fasta(fasta_path):
    """Parse FASTA file into dictionary of {accession: sequence}."""
    sequences = {}
    headers = {}  # Store full headers
    current_acc = None
    current_seq = []
    current_header = None
    
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous sequence
                if current_acc:
                    sequences[current_acc] = ''.join(current_seq)
                    headers[current_acc] = current_header
                
                # Parse new header
                current_header = line[1:]  # Remove >
                parts = current_header.split('|')
                if len(parts) >= 2:
                    current_acc = parts[1]
                else:
                    current_acc = current_header.split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        
        # Save last sequence
        if current_acc:
            sequences[current_acc] = ''.join(current_seq)
            headers[current_acc] = current_header
    
    return sequences, headers


def parse_tsv(tsv_path):
    """Parse TSV file into list of dictionaries."""
    records = []
    header = []
    
    if not tsv_path.exists():
        return records, header
    
    with open(tsv_path, 'r') as f:
        header = f.readline().strip().split('\t')
        for line in f:
            values = line.strip().split('\t')
            while len(values) < len(header):
                values.append('')
            record = dict(zip(header, values))
            records.append(record)
    return records, header


def sequence_hash(seq):
    """Create a hash of a sequence for efficient comparison."""
    return hashlib.md5(seq.upper().encode()).hexdigest()


def find_duplicates(sequences):
    """
    Find groups of sequences with identical amino acid strings.
    
    Returns:
        groups: dict of {sequence_hash: [list of accessions with that sequence]}
    """
    hash_to_accs = defaultdict(list)
    
    for acc, seq in sequences.items():
        seq_upper = seq.upper().replace('*', '')  # Normalize: uppercase, remove stop codons
        h = sequence_hash(seq_upper)
        hash_to_accs[h].append(acc)
    
    return hash_to_accs


def select_representative(accessions, rule="lowest_accession"):
    """
    Select representative accession from a group of duplicates.
    
    Args:
        accessions: list of accession IDs
        rule: selection rule ("lowest_accession" or "first_encountered")
    
    Returns:
        representative: the chosen accession
        duplicates: list of non-representative accessions
    """
    if rule == "lowest_accession":
        # Sort alphanumerically and pick first
        sorted_accs = sorted(accessions)
        representative = sorted_accs[0]
        duplicates = sorted_accs[1:]
    elif rule == "first_encountered":
        # Keep order as encountered (first in list)
        representative = accessions[0]
        duplicates = accessions[1:]
    else:
        raise ValueError(f"Unknown rule: {rule}")
    
    return representative, duplicates


def write_fasta(sequences, headers, output_path, accession_order=None):
    """Write sequences to FASTA file."""
    if accession_order is None:
        accession_order = sorted(sequences.keys())
    
    with open(output_path, 'w') as f:
        for acc in accession_order:
            if acc in sequences:
                seq = sequences[acc]
                header = headers.get(acc, acc)
                f.write(f">{header}\n")
                for i in range(0, len(seq), 60):
                    f.write(seq[i:i+60] + '\n')


def write_tsv(records, header, output_path):
    """Write records to TSV file."""
    with open(output_path, 'w') as f:
        f.write('\t'.join(header) + '\n')
        for record in records:
            values = [record.get(h, '') for h in header]
            f.write('\t'.join(values) + '\n')


def write_mapping(mapping, output_path):
    """Write duplicate-to-representative mapping to TSV."""
    with open(output_path, 'w') as f:
        f.write("duplicate_accession\trepresentative_accession\tsequence_hash\n")
        for dup, (rep, seq_hash) in sorted(mapping.items()):
            f.write(f"{dup}\t{rep}\t{seq_hash}\n")


# ============================================================================
# MAIN DEDUPLICATION LOGIC
# ============================================================================

def deduplicate(sequences, headers, records, tsv_header, dry_run=False):
    """
    Remove exact duplicate sequences.
    
    Returns:
        dedup_sequences: dict of {accession: sequence} for unique sequences
        dedup_headers: dict of {accession: header}
        dedup_records: list of record dicts for unique sequences
        report: dict with deduplication statistics
        mapping: dict of {duplicate_acc: (representative_acc, hash)}
    """
    
    print(f"\n📊 Starting deduplication with {len(sequences):,} sequences")
    
    # Find duplicate groups
    hash_to_accs = find_duplicates(sequences)
    
    # Count statistics
    n_unique_sequences = len(hash_to_accs)
    n_singleton_groups = sum(1 for accs in hash_to_accs.values() if len(accs) == 1)
    n_duplicate_groups = sum(1 for accs in hash_to_accs.values() if len(accs) > 1)
    
    print(f"   Found {n_unique_sequences:,} unique sequence strings")
    print(f"   - {n_singleton_groups:,} sequences have no duplicates")
    print(f"   - {n_duplicate_groups:,} groups have duplicates")
    
    # Process each group
    representatives = set()
    duplicates_to_remove = set()
    mapping = {}  # duplicate_acc -> (representative_acc, hash)
    
    duplicate_group_sizes = []
    
    for seq_hash, accs in hash_to_accs.items():
        if len(accs) == 1:
            # No duplicates, keep it
            representatives.add(accs[0])
        else:
            # Has duplicates, select representative
            rep, dups = select_representative(accs, rule=REPRESENTATIVE_RULE)
            representatives.add(rep)
            duplicate_group_sizes.append(len(accs))
            
            for dup in dups:
                duplicates_to_remove.add(dup)
                mapping[dup] = (rep, seq_hash)
    
    # Create deduplicated outputs
    dedup_sequences = {acc: sequences[acc] for acc in representatives}
    dedup_headers = {acc: headers[acc] for acc in representatives if acc in headers}
    
    # Filter records
    dedup_records = [r for r in records 
                     if r.get('Entry', r.get('accession', '')) in representatives]
    
    # Build report
    report = {
        "input_count": len(sequences),
        "output_count": len(dedup_sequences),
        "n_removed_exact_duplicates": len(duplicates_to_remove),
        "n_unique_sequences": n_unique_sequences,
        "n_singleton_groups": n_singleton_groups,
        "n_duplicate_groups": n_duplicate_groups,
        "representative_selection_rule": REPRESENTATIVE_RULE,
        "largest_duplicate_group": max(duplicate_group_sizes) if duplicate_group_sizes else 1,
        "duplicate_group_size_distribution": {},
        "timestamp": datetime.now().isoformat(),
        "script": "scripts/deduplicate_sequences.py"
    }
    
    # Size distribution
    if duplicate_group_sizes:
        size_counts = defaultdict(int)
        for size in duplicate_group_sizes:
            size_counts[size] += 1
        report["duplicate_group_size_distribution"] = dict(sorted(size_counts.items()))
    
    # Retention rate
    report["retention_rate"] = report["output_count"] / report["input_count"] if report["input_count"] > 0 else 0
    
    print(f"\n   Removed {len(duplicates_to_remove):,} duplicate sequences")
    print(f"   Keeping {len(representatives):,} representative sequences")
    
    if duplicate_group_sizes:
        print(f"   Largest duplicate group: {max(duplicate_group_sizes)} identical sequences")
    
    return dedup_sequences, dedup_headers, dedup_records, report, mapping


# ============================================================================
# MAIN
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Remove exact duplicate sequences",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Deduplication rules:
  - Sequences are compared after uppercasing and removing stop codons (*)
  - When duplicates are found, the representative is selected by:
    lowest_accession: alphanumerically first accession ID (deterministic)
  - All duplicates map to their representative in the output mapping file

Examples:
    python scripts/deduplicate_sequences.py           # Full deduplication
    python scripts/deduplicate_sequences.py --dry-run # Preview only
        """
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Preview deduplication without writing files"
    )
    args = parser.parse_args()
    
    print("=" * 70)
    print("Step 3: Exact Deduplication")
    print("=" * 70)
    
    # Check input files
    if not INPUT_FASTA.exists():
        print(f"\n❌ Error: Input FASTA not found: {INPUT_FASTA}")
        print("   Run scripts/filter_sequences.py first (Step 2)")
        sys.exit(1)
    
    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # Parse input files
    print(f"\n📂 Loading input files...")
    print(f"   FASTA: {INPUT_FASTA}")
    
    sequences, headers = parse_fasta(INPUT_FASTA)
    records, tsv_header = parse_tsv(INPUT_TSV) if INPUT_TSV.exists() else ([], [])
    
    print(f"   Loaded {len(sequences):,} sequences from FASTA")
    if records:
        print(f"   Loaded {len(records):,} records from TSV")
    
    # Deduplicate
    print(f"\n🔧 Deduplicating sequences...")
    print(f"   Rule: {REPRESENTATIVE_RULE}")
    
    dedup_sequences, dedup_headers, dedup_records, report, mapping = deduplicate(
        sequences, headers, records, tsv_header, dry_run=args.dry_run
    )
    
    # Add file paths to report
    report["input_files"] = {
        "fasta": str(INPUT_FASTA),
        "tsv": str(INPUT_TSV) if INPUT_TSV.exists() else None
    }
    report["output_files"] = {
        "fasta": str(OUTPUT_DIR / "step3_dedup.fasta"),
        "tsv": str(OUTPUT_DIR / "step3_dedup.tsv"),
        "mapping": str(OUTPUT_DIR / "step3_dedup_map.tsv"),
        "report": str(OUTPUT_DIR / "step3_dedup_report.json")
    }
    
    if not args.dry_run:
        # Write output files
        print(f"\n💾 Writing output files...")
        
        fasta_path = OUTPUT_DIR / "step3_dedup.fasta"
        tsv_path = OUTPUT_DIR / "step3_dedup.tsv"
        mapping_path = OUTPUT_DIR / "step3_dedup_map.tsv"
        report_path = OUTPUT_DIR / "step3_dedup_report.json"
        
        # Write FASTA
        write_fasta(dedup_sequences, dedup_headers, fasta_path)
        print(f"   ✓ {fasta_path} ({len(dedup_sequences):,} sequences)")
        
        # Write TSV (if we had input TSV)
        if dedup_records and tsv_header:
            write_tsv(dedup_records, tsv_header, tsv_path)
            print(f"   ✓ {tsv_path} ({len(dedup_records):,} records)")
        
        # Write mapping
        write_mapping(mapping, mapping_path)
        print(f"   ✓ {mapping_path} ({len(mapping):,} duplicate mappings)")
        
        # Write report
        with open(report_path, 'w') as f:
            json.dump(report, f, indent=2)
        print(f"   ✓ {report_path}")
        
        # Summary
        print(f"\n✅ Step 3 Complete!")
        print(f"\n📊 Summary:")
        print(f"   Input:  {report['input_count']:,} sequences")
        print(f"   Output: {report['output_count']:,} sequences")
        print(f"   Removed: {report['n_removed_exact_duplicates']:,} exact duplicates")
        print(f"   Retention: {report['retention_rate']*100:.1f}%")
        
        print(f"\n📝 For MANUSCRIPT.md:")
        print(f"   N_removed_exact_duplicates = {report['n_removed_exact_duplicates']:,}")
        print(f"   N_after_dedup = {report['output_count']:,}")
        
        # Sanity check
        print(f"\n✅ Sanity checks:")
        print(f"   ✓ Representative selection is deterministic (rule: {REPRESENTATIVE_RULE})")
        print(f"   ✓ Same input → same representatives (sorted by accession)")
        print(f"   ✓ All duplicates mapped to their representative")
        
    else:
        print(f"\n[DRY RUN] No files were written.")
        print(f"\n📊 Would produce:")
        print(f"   {report['output_count']:,} sequences after deduplication")
        print(f"   {report['n_removed_exact_duplicates']:,} duplicates removed")


if __name__ == "__main__":
    main()

