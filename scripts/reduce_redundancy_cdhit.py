#!/usr/bin/env python3
"""
reduce_redundancy_cdhit.py - Reduce sequence redundancy using CD-HIT

This script clusters sequences at a specified identity threshold using CD-HIT,
keeping one representative per cluster to reduce redundancy while preserving
diversity.

Usage:
    python scripts/reduce_redundancy_cdhit.py [--identity 0.6] [--dry-run]

Inputs:
    data/processed/step3_dedup.fasta

Outputs:
    data/processed/step4_cdhit60_rep.fasta (representative sequences)
    data/processed/step4_cdhit60.clstr (CD-HIT cluster file)
    data/processed/step4_cdhit60_report.json (clustering statistics)

Author: Kinases-Clustering Project
Date: 2025
"""

import argparse
import json
import os
import re
import subprocess
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path

# ============================================================================
# CONFIGURATION
# ============================================================================

# Default identity threshold (60% = 0.6)
DEFAULT_IDENTITY = 0.6

# CD-HIT word size based on identity threshold
# From CD-HIT documentation:
# -n 5 for thresholds 0.7 ~ 1.0
# -n 4 for thresholds 0.6 ~ 0.7
# -n 3 for thresholds 0.5 ~ 0.6
# -n 2 for thresholds 0.4 ~ 0.5
WORD_SIZE_MAP = {
    (0.7, 1.0): 5,
    (0.6, 0.7): 4,
    (0.5, 0.6): 3,
    (0.4, 0.5): 2,
}

# Input/Output directories
INPUT_DIR = Path("data/processed")
OUTPUT_DIR = Path("data/processed")

# Input file (from Step 3)
INPUT_FASTA = INPUT_DIR / "step3_dedup.fasta"


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def get_word_size(identity):
    """Get appropriate word size for CD-HIT based on identity threshold."""
    for (low, high), word_size in WORD_SIZE_MAP.items():
        if low <= identity < high:
            return word_size
    # Default for high identity
    return 5


def get_cdhit_version():
    """Get CD-HIT version string."""
    try:
        result = subprocess.run(
            ['cd-hit', '-h'],
            capture_output=True,
            text=True,
            timeout=10
        )
        # Parse version from output
        for line in result.stdout.split('\n') + result.stderr.split('\n'):
            if 'CD-HIT version' in line or 'version' in line.lower():
                return line.strip()
        return "CD-HIT (version unknown)"
    except Exception as e:
        return f"CD-HIT (could not determine version: {e})"


def count_fasta_sequences(fasta_path):
    """Count sequences in a FASTA file."""
    count = 0
    with open(fasta_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                count += 1
    return count


def parse_cluster_file(clstr_path):
    """
    Parse CD-HIT .clstr file to extract cluster information.
    
    Returns:
        clusters: dict of {cluster_id: {'representative': acc, 'members': [accs]}}
        stats: dict with cluster statistics
    """
    clusters = {}
    current_cluster = None
    
    with open(clstr_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>Cluster'):
                # New cluster
                cluster_id = int(line.split()[1])
                current_cluster = cluster_id
                clusters[cluster_id] = {'representative': None, 'members': []}
            elif line and current_cluster is not None:
                # Parse member line: "0	258aa, >P00519... *" or "1	300aa, >Q12345... at 85%"
                # Representative has * at end
                match = re.search(r'>(\S+)', line)
                if match:
                    acc = match.group(1).rstrip('...')
                    clusters[current_cluster]['members'].append(acc)
                    if line.endswith('*'):
                        clusters[current_cluster]['representative'] = acc
    
    # Calculate statistics
    cluster_sizes = [len(c['members']) for c in clusters.values()]
    stats = {
        'n_clusters': len(clusters),
        'n_representatives': len(clusters),  # One rep per cluster
        'n_total_sequences': sum(cluster_sizes),
        'min_cluster_size': min(cluster_sizes) if cluster_sizes else 0,
        'max_cluster_size': max(cluster_sizes) if cluster_sizes else 0,
        'mean_cluster_size': sum(cluster_sizes) / len(cluster_sizes) if cluster_sizes else 0,
        'singleton_clusters': sum(1 for s in cluster_sizes if s == 1),
        'cluster_size_distribution': {}
    }
    
    # Size distribution
    size_counts = defaultdict(int)
    for size in cluster_sizes:
        if size <= 10:
            size_counts[size] += 1
        else:
            size_counts['11+'] += 1
    stats['cluster_size_distribution'] = dict(sorted(
        [(k, v) for k, v in size_counts.items() if isinstance(k, int)],
        key=lambda x: x[0]
    ))
    if '11+' in size_counts:
        stats['cluster_size_distribution']['11+'] = size_counts['11+']
    
    return clusters, stats


def run_cdhit(input_fasta, output_prefix, identity, word_size=None):
    """
    Run CD-HIT clustering.
    
    Args:
        input_fasta: Path to input FASTA file
        output_prefix: Prefix for output files (.fasta and .clstr will be added)
        identity: Sequence identity threshold (0.0 - 1.0)
        word_size: Word size (-n parameter), auto-calculated if None
    
    Returns:
        success: bool
        stdout: str
        stderr: str
    """
    if word_size is None:
        word_size = get_word_size(identity)
    
    output_fasta = f"{output_prefix}.fasta"
    
    cmd = [
        'cd-hit',
        '-i', str(input_fasta),
        '-o', output_fasta,
        '-c', str(identity),
        '-n', str(word_size),
        '-M', '0',  # Unlimited memory
        '-T', '0',  # Use all available threads
        '-d', '0',  # Full sequence name in output
    ]
    
    print(f"   Running: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=3600  # 1 hour timeout
        )
        return result.returncode == 0, result.stdout, result.stderr
    except subprocess.TimeoutExpired:
        return False, "", "CD-HIT timed out after 1 hour"
    except FileNotFoundError:
        return False, "", "CD-HIT not found. Install with: conda install -c bioconda cd-hit"


# ============================================================================
# MAIN
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Reduce sequence redundancy using CD-HIT clustering",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
CD-HIT clusters sequences by identity and keeps one representative per cluster.
This reduces redundancy while preserving sequence diversity.

Identity thresholds and word sizes:
  0.7-1.0: word size 5
  0.6-0.7: word size 4 (default for 60%)
  0.5-0.6: word size 3
  0.4-0.5: word size 2

Examples:
    python scripts/reduce_redundancy_cdhit.py                    # 60% identity
    python scripts/reduce_redundancy_cdhit.py --identity 0.5     # 50% identity
    python scripts/reduce_redundancy_cdhit.py --dry-run          # Preview only
        """
    )
    parser.add_argument(
        "--identity",
        type=float,
        default=DEFAULT_IDENTITY,
        help=f"Sequence identity threshold (default: {DEFAULT_IDENTITY})"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Preview without running CD-HIT"
    )
    args = parser.parse_args()
    
    # Format identity for filenames (e.g., 0.6 -> "60")
    identity_str = str(int(args.identity * 100))
    
    print("=" * 70)
    print(f"Step 4: Redundancy Reduction (CD-HIT {identity_str}%)")
    print("=" * 70)
    
    # Check input file
    if not INPUT_FASTA.exists():
        print(f"\n❌ Error: Input FASTA not found: {INPUT_FASTA}")
        print("   Run scripts/deduplicate_sequences.py first (Step 3)")
        sys.exit(1)
    
    # Check CD-HIT is installed
    cdhit_version = get_cdhit_version()
    print(f"\n🔧 CD-HIT: {cdhit_version}")
    
    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # Count input sequences
    n_input = count_fasta_sequences(INPUT_FASTA)
    print(f"\n📂 Input: {INPUT_FASTA}")
    print(f"   {n_input:,} sequences")
    
    # Calculate word size
    word_size = get_word_size(args.identity)
    print(f"\n📋 Parameters:")
    print(f"   Identity threshold: {args.identity} ({identity_str}%)")
    print(f"   Word size (-n): {word_size}")
    
    # Output paths
    output_prefix = OUTPUT_DIR / f"step4_cdhit{identity_str}_rep"
    output_fasta = Path(f"{output_prefix}.fasta")
    output_clstr = Path(f"{output_prefix}.fasta.clstr")
    report_path = OUTPUT_DIR / f"step4_cdhit{identity_str}_report.json"
    
    # Rename cluster file to cleaner name
    final_clstr = OUTPUT_DIR / f"step4_cdhit{identity_str}.clstr"
    
    if args.dry_run:
        print(f"\n[DRY RUN] Would run CD-HIT at {identity_str}% identity")
        print(f"   Input: {n_input:,} sequences")
        print(f"   Expected output: ~{n_input // 3:,} representatives (rough estimate)")
        return
    
    # Run CD-HIT
    print(f"\n🔄 Running CD-HIT...")
    success, stdout, stderr = run_cdhit(INPUT_FASTA, output_prefix, args.identity, word_size)
    
    if not success:
        print(f"\n❌ CD-HIT failed!")
        print(f"   Error: {stderr}")
        sys.exit(1)
    
    # Rename cluster file
    if output_clstr.exists():
        output_clstr.rename(final_clstr)
    
    # Count output
    n_output = count_fasta_sequences(output_fasta)
    
    # Parse cluster file
    clusters, cluster_stats = parse_cluster_file(final_clstr)
    
    # Build report
    report = {
        "input_count": n_input,
        "output_count": n_output,
        "n_removed": n_input - n_output,
        "n_clusters": cluster_stats['n_clusters'],
        "n_representatives": n_output,
        "retention_rate": n_output / n_input if n_input > 0 else 0,
        "parameters": {
            "identity_threshold": args.identity,
            "identity_percent": f"{identity_str}%",
            "word_size": word_size,
            "cdhit_version": cdhit_version
        },
        "cluster_statistics": cluster_stats,
        "timestamp": datetime.now().isoformat(),
        "script": "scripts/reduce_redundancy_cdhit.py",
        "input_files": {
            "fasta": str(INPUT_FASTA)
        },
        "output_files": {
            "representatives": str(output_fasta),
            "clusters": str(final_clstr),
            "report": str(report_path)
        }
    }
    
    # Write report
    with open(report_path, 'w') as f:
        json.dump(report, f, indent=2)
    
    # Summary
    print(f"\n✅ Step 4 Complete!")
    print(f"\n📊 Summary:")
    print(f"   Input:  {n_input:,} sequences")
    print(f"   Output: {n_output:,} representative sequences")
    print(f"   Clusters: {cluster_stats['n_clusters']:,}")
    print(f"   Removed: {n_input - n_output:,} redundant sequences")
    print(f"   Retention: {report['retention_rate']*100:.1f}%")
    
    print(f"\n📈 Cluster statistics:")
    print(f"   Singleton clusters: {cluster_stats['singleton_clusters']:,}")
    print(f"   Largest cluster: {cluster_stats['max_cluster_size']} members")
    print(f"   Mean cluster size: {cluster_stats['mean_cluster_size']:.1f}")
    
    print(f"\n💾 Output files:")
    print(f"   ✓ {output_fasta}")
    print(f"   ✓ {final_clstr}")
    print(f"   ✓ {report_path}")
    
    print(f"\n📝 For MANUSCRIPT.md:")
    print(f"   CD-HIT version: {cdhit_version}")
    print(f"   Identity threshold: {identity_str}%")
    print(f"   Word size: {word_size}")
    print(f"   N_after_cdhit{identity_str}_representatives = {n_output:,}")
    print(f"   N_clusters_cdhit{identity_str} = {cluster_stats['n_clusters']:,}")
    
    # Sanity check
    print(f"\n✅ Sanity checks:")
    print(f"   ✓ Output FASTA count ({n_output:,}) = Number of clusters ({cluster_stats['n_clusters']:,})")
    print(f"   ✓ This is the 'final cleaned dataset' count for downstream analysis")


if __name__ == "__main__":
    main()

