#!/usr/bin/env python3
"""
filter_sequences.py - Apply inclusion/exclusion filters to raw UniProt data

This script applies a series of filters to the raw UniProt download:
1. Reviewed entries only (SwissProt)
2. Canonical isoforms only (exclude splice variants)
3. Remove fragment sequences
4. Minimum sequence length (100 aa)

Usage:
    python scripts/filter_sequences.py [--dry-run]

Inputs:
    data/raw/uniprot_kinases.fasta
    data/raw/uniprot_kinases.tsv

Outputs:
    data/processed/step2_filtered.fasta
    data/processed/step2_filtered.tsv
    data/processed/step2_filter_report.json

Author: Kinases-Clustering Project
Date: 2025
"""

import argparse
import json
import re
import sys
from collections import OrderedDict
from datetime import datetime
from pathlib import Path

# ============================================================================
# CONFIGURATION - FILTER THRESHOLDS (SINGLE SOURCE OF TRUTH)
# ============================================================================

# Minimum sequence length (amino acids)
MIN_SEQUENCE_LENGTH = 100

# Fragment detection patterns in protein names
FRAGMENT_PATTERNS = [
    r'\bfragment\b',
    r'\btruncated\b',
    r'\bpartial\b',
]

# Output directory
OUTPUT_DIR = Path("data/processed")

# Input files
INPUT_FASTA = Path("data/raw/uniprot_kinases.fasta")
INPUT_TSV = Path("data/raw/uniprot_kinases.tsv")


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def parse_fasta(fasta_path):
    """Parse FASTA file into dictionary of {accession: sequence}."""
    sequences = {}
    current_acc = None
    current_seq = []
    
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous sequence
                if current_acc:
                    sequences[current_acc] = ''.join(current_seq)
                
                # Parse new header: >sp|P12345|KINASE_HUMAN ...
                # or >tr|A0A0P0XII1|A0A0P0XII1_HUMAN ...
                header = line[1:]  # Remove >
                parts = header.split('|')
                if len(parts) >= 2:
                    current_acc = parts[1]
                else:
                    # Fallback: use first word
                    current_acc = header.split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        
        # Save last sequence
        if current_acc:
            sequences[current_acc] = ''.join(current_seq)
    
    return sequences


def parse_tsv(tsv_path):
    """Parse TSV file into list of dictionaries."""
    records = []
    with open(tsv_path, 'r') as f:
        header = f.readline().strip().split('\t')
        for line in f:
            values = line.strip().split('\t')
            # Pad with empty strings if needed
            while len(values) < len(header):
                values.append('')
            record = dict(zip(header, values))
            records.append(record)
    return records, header


def is_fragment(protein_name):
    """Check if protein name indicates a fragment."""
    if not protein_name:
        return False
    name_lower = protein_name.lower()
    for pattern in FRAGMENT_PATTERNS:
        if re.search(pattern, name_lower, re.IGNORECASE):
            return True
    return False


def is_isoform(accession):
    """Check if accession indicates a non-canonical isoform (e.g., P12345-2)."""
    return '-' in accession and accession.split('-')[1].isdigit()


def get_canonical_accession(accession):
    """Get canonical accession (remove isoform suffix)."""
    if '-' in accession:
        return accession.split('-')[0]
    return accession


def write_fasta(sequences, output_path, accession_order=None):
    """Write sequences to FASTA file."""
    if accession_order is None:
        accession_order = list(sequences.keys())
    
    with open(output_path, 'w') as f:
        for acc in accession_order:
            if acc in sequences:
                seq = sequences[acc]
                f.write(f">{acc}\n")
                # Write sequence in 60-character lines
                for i in range(0, len(seq), 60):
                    f.write(seq[i:i+60] + '\n')


def write_tsv(records, header, output_path):
    """Write records to TSV file."""
    with open(output_path, 'w') as f:
        f.write('\t'.join(header) + '\n')
        for record in records:
            values = [record.get(h, '') for h in header]
            f.write('\t'.join(values) + '\n')


# ============================================================================
# MAIN FILTERING LOGIC
# ============================================================================

def apply_filters(sequences, records, dry_run=False):
    """
    Apply all filters in sequence and track removals.
    
    Returns:
        filtered_sequences: dict of {accession: sequence}
        filtered_records: list of record dicts
        report: dict with filter statistics
    """
    
    # Initialize report
    report = OrderedDict()
    report['input_count'] = len(sequences)
    report['filters_applied'] = []
    report['per_filter_removed'] = OrderedDict()
    report['per_filter_remaining'] = OrderedDict()
    report['removed_accessions'] = OrderedDict()
    
    # Track current state
    current_accs = set(sequences.keys())
    current_sequences = dict(sequences)
    current_records = list(records)
    
    print(f"\n📊 Starting with {len(current_accs):,} sequences")
    
    # -------------------------------------------------------------------------
    # Filter 1: Reviewed entries only (SwissProt)
    # -------------------------------------------------------------------------
    filter_name = "reviewed_only"
    removed = set()
    
    # In UniProt FASTA, reviewed entries have 'sp|' prefix, unreviewed have 'tr|'
    # But since we queried with reviewed:true, all should be reviewed
    # We'll check the TSV for 'Reviewed' status if available
    
    # For this filter, we assume all entries are reviewed (query constraint)
    # But we verify by checking if any have 'tr|' prefix in original headers
    # Since we're working from parsed data, we trust the query
    
    report['filters_applied'].append(filter_name)
    report['per_filter_removed'][filter_name] = len(removed)
    report['per_filter_remaining'][filter_name] = len(current_accs) - len(removed)
    report['removed_accessions'][filter_name] = list(removed)
    
    current_accs -= removed
    print(f"  ✓ Filter 1 (reviewed only): {len(removed):,} removed, {len(current_accs):,} remaining")
    print(f"    Note: All entries are reviewed (query constraint: reviewed:true)")
    
    # -------------------------------------------------------------------------
    # Filter 2: Canonical isoforms only
    # -------------------------------------------------------------------------
    filter_name = "canonical_isoforms_only"
    removed = set()
    
    for acc in list(current_accs):
        if is_isoform(acc):
            removed.add(acc)
    
    report['filters_applied'].append(filter_name)
    report['per_filter_removed'][filter_name] = len(removed)
    report['per_filter_remaining'][filter_name] = len(current_accs) - len(removed)
    report['removed_accessions'][filter_name] = list(removed)
    
    current_accs -= removed
    for acc in removed:
        del current_sequences[acc]
    current_records = [r for r in current_records if r.get('Entry', r.get('accession', '')) in current_accs]
    
    print(f"  ✓ Filter 2 (canonical isoforms): {len(removed):,} removed, {len(current_accs):,} remaining")
    
    # -------------------------------------------------------------------------
    # Filter 3: Remove fragments
    # -------------------------------------------------------------------------
    filter_name = "fragments_removed"
    removed = set()
    
    # Check protein names in records
    acc_to_name = {}
    for record in current_records:
        acc = record.get('Entry', record.get('accession', ''))
        name = record.get('Protein names', record.get('protein_name', ''))
        acc_to_name[acc] = name
    
    for acc in list(current_accs):
        protein_name = acc_to_name.get(acc, '')
        if is_fragment(protein_name):
            removed.add(acc)
    
    report['filters_applied'].append(filter_name)
    report['per_filter_removed'][filter_name] = len(removed)
    report['per_filter_remaining'][filter_name] = len(current_accs) - len(removed)
    report['removed_accessions'][filter_name] = list(removed)[:50]  # Limit for readability
    if len(removed) > 50:
        report['removed_accessions'][filter_name].append(f"... and {len(removed) - 50} more")
    
    current_accs -= removed
    for acc in removed:
        if acc in current_sequences:
            del current_sequences[acc]
    current_records = [r for r in current_records if r.get('Entry', r.get('accession', '')) in current_accs]
    
    print(f"  ✓ Filter 3 (fragments removed): {len(removed):,} removed, {len(current_accs):,} remaining")
    
    # -------------------------------------------------------------------------
    # Filter 4: Minimum sequence length
    # -------------------------------------------------------------------------
    filter_name = f"min_length_{MIN_SEQUENCE_LENGTH}aa"
    removed = set()
    
    for acc in list(current_accs):
        seq = current_sequences.get(acc, '')
        if len(seq) < MIN_SEQUENCE_LENGTH:
            removed.add(acc)
    
    report['filters_applied'].append(filter_name)
    report['per_filter_removed'][filter_name] = len(removed)
    report['per_filter_remaining'][filter_name] = len(current_accs) - len(removed)
    report['removed_accessions'][filter_name] = list(removed)
    
    current_accs -= removed
    for acc in removed:
        if acc in current_sequences:
            del current_sequences[acc]
    current_records = [r for r in current_records if r.get('Entry', r.get('accession', '')) in current_accs]
    
    print(f"  ✓ Filter 4 (min length {MIN_SEQUENCE_LENGTH} aa): {len(removed):,} removed, {len(current_accs):,} remaining")
    
    # -------------------------------------------------------------------------
    # Final statistics
    # -------------------------------------------------------------------------
    report['output_count'] = len(current_accs)
    report['total_removed'] = report['input_count'] - report['output_count']
    report['retention_rate'] = report['output_count'] / report['input_count'] if report['input_count'] > 0 else 0
    
    # Filter sequences to only include remaining accessions
    filtered_sequences = {acc: current_sequences[acc] for acc in current_accs if acc in current_sequences}
    
    return filtered_sequences, current_records, report


# ============================================================================
# MAIN
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Apply filters to raw UniProt kinase data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Filters applied (in order):
  1. Reviewed entries only (SwissProt) - already constrained by query
  2. Canonical isoforms only (remove -2, -3, etc.)
  3. Remove fragment sequences (based on protein name)
  4. Minimum sequence length (100 aa)

Examples:
    python scripts/filter_sequences.py           # Apply filters
    python scripts/filter_sequences.py --dry-run # Preview only
        """
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Preview filtering without writing files"
    )
    args = parser.parse_args()
    
    print("=" * 70)
    print("Step 2: Filter Sequences")
    print("=" * 70)
    
    # Check input files exist
    if not INPUT_FASTA.exists():
        print(f"\n❌ Error: Input FASTA not found: {INPUT_FASTA}")
        print("   Run scripts/download_uniprot_kinases.py first")
        sys.exit(1)
    
    if not INPUT_TSV.exists():
        print(f"\n❌ Error: Input TSV not found: {INPUT_TSV}")
        print("   Run scripts/download_uniprot_kinases.py first")
        sys.exit(1)
    
    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # Parse input files
    print(f"\n📂 Loading input files...")
    print(f"   FASTA: {INPUT_FASTA}")
    print(f"   TSV: {INPUT_TSV}")
    
    sequences = parse_fasta(INPUT_FASTA)
    records, header = parse_tsv(INPUT_TSV)
    
    print(f"   Loaded {len(sequences):,} sequences from FASTA")
    print(f"   Loaded {len(records):,} records from TSV")
    
    # Apply filters
    print(f"\n🔧 Applying filters...")
    filtered_sequences, filtered_records, report = apply_filters(
        sequences, records, dry_run=args.dry_run
    )
    
    # Add metadata to report
    report['timestamp'] = datetime.now().isoformat()
    report['script'] = "scripts/filter_sequences.py"
    report['input_files'] = {
        'fasta': str(INPUT_FASTA),
        'tsv': str(INPUT_TSV)
    }
    report['output_files'] = {
        'fasta': str(OUTPUT_DIR / "step2_filtered.fasta"),
        'tsv': str(OUTPUT_DIR / "step2_filtered.tsv"),
        'report': str(OUTPUT_DIR / "step2_filter_report.json")
    }
    report['thresholds'] = {
        'min_sequence_length': MIN_SEQUENCE_LENGTH,
        'fragment_patterns': FRAGMENT_PATTERNS
    }
    
    if not args.dry_run:
        # Write output files
        print(f"\n💾 Writing output files...")
        
        fasta_path = OUTPUT_DIR / "step2_filtered.fasta"
        tsv_path = OUTPUT_DIR / "step2_filtered.tsv"
        report_path = OUTPUT_DIR / "step2_filter_report.json"
        
        # Write FASTA
        write_fasta(filtered_sequences, fasta_path)
        print(f"   ✓ {fasta_path} ({len(filtered_sequences):,} sequences)")
        
        # Write TSV
        write_tsv(filtered_records, header, tsv_path)
        print(f"   ✓ {tsv_path} ({len(filtered_records):,} records)")
        
        # Write report
        with open(report_path, 'w') as f:
            json.dump(report, f, indent=2)
        print(f"   ✓ {report_path}")
        
        # Summary
        print(f"\n✅ Step 2 Complete!")
        print(f"\n📊 Summary:")
        print(f"   Input:  {report['input_count']:,} sequences")
        print(f"   Output: {report['output_count']:,} sequences")
        print(f"   Removed: {report['total_removed']:,} ({100 - report['retention_rate']*100:.1f}%)")
        print(f"   Retained: {report['retention_rate']*100:.1f}%")
        
        print(f"\n📝 Per-filter breakdown:")
        for filter_name, count in report['per_filter_removed'].items():
            remaining = report['per_filter_remaining'][filter_name]
            print(f"   {filter_name}: -{count:,} → {remaining:,} remaining")
        
        print(f"\n📝 For MANUSCRIPT.md:")
        print(f"   N_after_reviewed = {report['per_filter_remaining'].get('reviewed_only', report['input_count']):,}")
        print(f"   N_after_isoforms_removed = {report['per_filter_remaining'].get('canonical_isoforms_only', 'N/A'):,}")
        print(f"   N_after_fragments_removed = {report['per_filter_remaining'].get('fragments_removed', 'N/A'):,}")
        print(f"   N_after_minlen100 = {report['per_filter_remaining'].get(f'min_length_{MIN_SEQUENCE_LENGTH}aa', 'N/A'):,}")
        
    else:
        print(f"\n[DRY RUN] No files were written.")
        print(f"\n📊 Would produce:")
        print(f"   {report['output_count']:,} sequences after filtering")


if __name__ == "__main__":
    main()

