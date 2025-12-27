#!/usr/bin/env python3
"""
Step 7: Domain Extraction (HMMER + Pfam)

This script extracts kinase catalytic domains using HMMER with Pfam profiles.
It creates all required artifacts and generates a comprehensive report.

Usage:
    python scripts/extract_domains.py

Inputs:
    - data/processed/kinases_normalized.csv (from earlier steps)
    - data/processed/labels.csv (from Step 5)
    - data/hmm_profiles/PF00069.hmm (or downloads from Pfam)
    - data/hmm_profiles/PF07714.hmm (or downloads from Pfam)

Outputs:
    - data/domains/hmmer_domtblout_E0001.txt
    - data/domains/hmmer_domtblout_E001.txt
    - data/domains/domains_E0001.fasta
    - data/domains/domains_E001.fasta
    - data/domains/domain_coords_E0001.tsv
    - data/domains/domain_coords_E001.tsv
    - data/domains/domain_extraction_report.json
"""

import os
import sys
import json
import subprocess
import tempfile
from pathlib import Path
from collections import Counter
from datetime import datetime

import pandas as pd


# Configuration
PFAM_PROFILES = {
    'PF00069': 'Pkinase',        # Protein kinase domain (Ser/Thr/Tyr)
    'PF07714': 'Pkinase_Tyr',    # Protein tyrosine kinase
}

E_VALUE_THRESHOLDS = [0.001, 0.01]


def check_hmmer_version():
    """Get HMMER version."""
    try:
        result = subprocess.run(
            ['hmmsearch', '-h'],
            capture_output=True,
            text=True,
            timeout=5
        )
        # Parse version from help output
        for line in result.stdout.split('\n'):
            if 'HMMER' in line and ('hmmer' in line.lower() or 'version' in line.lower()):
                return line.strip()
        return "HMMER (version unknown)"
    except Exception:
        return None


def get_pfam_version(hmm_file):
    """Extract Pfam version from HMM file."""
    try:
        with open(hmm_file, 'r') as f:
            for line in f:
                if line.startswith('ACC'):
                    return line.split()[1].strip()
                if line.startswith('NAME'):
                    name = line.split()[1].strip()
        return "Unknown"
    except Exception:
        return "Unknown"


def parse_domtblout(domtblout_file, e_threshold):
    """Parse HMMER domtblout file and filter by E-value threshold."""
    hits = []
    
    if not os.path.exists(domtblout_file):
        return hits
    
    with open(domtblout_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.split()
            if len(parts) < 23:
                continue
            
            target_name = parts[0]
            pfam_name = parts[3]
            pfam_acc = parts[4].split('.')[0] if '.' in parts[4] else parts[4]
            
            # Domain-specific values
            i_evalue = float(parts[12])  # Independent E-value
            score = float(parts[13])
            
            # Envelope coordinates (conservative boundaries)
            env_from = int(parts[19])
            env_to = int(parts[20])
            
            if i_evalue <= e_threshold:
                hits.append({
                    'uniprot_id': target_name,
                    'pfam_id': pfam_acc,
                    'pfam_name': pfam_name,
                    'env_from': env_from,
                    'env_to': env_to,
                    'domain_length': env_to - env_from + 1,
                    'evalue': i_evalue,
                    'bitscore': score
                })
    
    return hits


def get_best_domain_per_protein(hits):
    """Keep only the best-scoring domain per protein."""
    # Group by uniprot_id, keep highest score
    best = {}
    for hit in hits:
        uid = hit['uniprot_id']
        if uid not in best or hit['bitscore'] > best[uid]['bitscore']:
            best[uid] = hit
    return list(best.values())


def extract_domain_sequence(full_sequence, env_from, env_to):
    """Extract domain sequence using envelope coordinates (1-indexed)."""
    # Convert to 0-indexed Python slicing
    return full_sequence[env_from - 1:env_to]


def main():
    # Paths
    processed_dir = Path("data/processed")
    domains_dir = Path("data/domains")
    hmm_dir = Path("data/hmm_profiles")
    
    domains_dir.mkdir(parents=True, exist_ok=True)
    
    # Load input data
    print("Loading input data...")
    sequences_file = processed_dir / "kinases_normalized.csv"
    labels_file = processed_dir / "labels.csv"
    
    if not sequences_file.exists():
        print(f"ERROR: {sequences_file} not found")
        sys.exit(1)
    
    sequences_df = pd.read_csv(sequences_file)
    labels_df = pd.read_csv(labels_file) if labels_file.exists() else None
    
    print(f"  Loaded {len(sequences_df)} sequences")
    
    # Create sequence lookup
    seq_lookup = dict(zip(sequences_df['uniprot_id'], sequences_df['sequence']))
    
    # Get HMMER version
    hmmer_version = check_hmmer_version()
    if not hmmer_version:
        print("WARNING: HMMER not found. Using existing domtblout files.")
    else:
        print(f"  HMMER: {hmmer_version}")
    
    # Get Pfam versions
    pfam_versions = {}
    for pfam_id in PFAM_PROFILES:
        hmm_file = hmm_dir / f"{pfam_id}.hmm"
        if hmm_file.exists():
            pfam_versions[pfam_id] = get_pfam_version(hmm_file)
    
    # Report structure
    report = {
        "step": 7,
        "name": "Domain Extraction",
        "timestamp": datetime.now().isoformat(),
        "tools": {
            "hmmer": hmmer_version or "Not available (using cached results)",
            "pfam_profiles": PFAM_PROFILES,
            "pfam_versions": pfam_versions
        },
        "parameters": {
            "extraction_rule": "envelope coordinates (env_from, env_to)",
            "selection_rule": "best-scoring domain per protein (highest bitscore)",
            "e_value_thresholds": E_VALUE_THRESHOLDS
        },
        "results": {}
    }
    
    # Process each E-value threshold
    for e_thresh in E_VALUE_THRESHOLDS:
        e_str = str(e_thresh).replace('.', '')
        print(f"\n{'='*60}")
        print(f"Processing E-value threshold: {e_thresh}")
        print(f"{'='*60}")
        
        # Collect hits from all Pfam profiles
        all_hits = []
        
        for pfam_id in PFAM_PROFILES:
            # Look for existing domtblout file
            domtblout_patterns = [
                hmm_dir / f"hmmer_{pfam_id}_e{e_thresh}.domtblout",
                hmm_dir / f"hmmer_{pfam_id}_e{e_str}.domtblout",
                hmm_dir / "hmmer_results.domtblout"
            ]
            
            domtblout_file = None
            for pattern in domtblout_patterns:
                if pattern.exists():
                    domtblout_file = pattern
                    break
            
            if domtblout_file:
                print(f"  Parsing {domtblout_file}...")
                hits = parse_domtblout(domtblout_file, e_thresh)
                print(f"    Found {len(hits)} hits for {pfam_id}")
                all_hits.extend(hits)
            else:
                print(f"  WARNING: No domtblout file found for {pfam_id}")
        
        # Get best domain per protein
        best_hits = get_best_domain_per_protein(all_hits)
        print(f"\n  Best domain per protein: {len(best_hits)} sequences")
        
        # Create domain coordinates TSV
        coords_file = domains_dir / f"domain_coords_E{e_str}.tsv"
        coords_df = pd.DataFrame(best_hits)
        coords_df.to_csv(coords_file, sep='\t', index=False)
        print(f"  Saved: {coords_file}")
        
        # Extract domain sequences and create FASTA
        fasta_file = domains_dir / f"domains_E{e_str}.fasta"
        with open(fasta_file, 'w') as f:
            for hit in best_hits:
                uid = hit['uniprot_id']
                if uid in seq_lookup:
                    domain_seq = extract_domain_sequence(
                        seq_lookup[uid],
                        hit['env_from'],
                        hit['env_to']
                    )
                    f.write(f">{uid}|{hit['pfam_id']}|{hit['env_from']}-{hit['env_to']}\n")
                    f.write(f"{domain_seq}\n")
        
        print(f"  Saved: {fasta_file}")
        
        # Copy/merge domtblout files
        merged_domtblout = domains_dir / f"hmmer_domtblout_E{e_str}.txt"
        with open(merged_domtblout, 'w') as outf:
            outf.write(f"# HMMER domain table output\n")
            outf.write(f"# E-value threshold: {e_thresh}\n")
            outf.write(f"# Pfam profiles: {', '.join(PFAM_PROFILES.keys())}\n")
            outf.write(f"# Generated: {datetime.now().isoformat()}\n")
            outf.write("#\n")
            
            for pfam_id in PFAM_PROFILES:
                domtblout_patterns = [
                    hmm_dir / f"hmmer_{pfam_id}_e{e_thresh}.domtblout",
                    hmm_dir / "hmmer_results.domtblout"
                ]
                
                for pattern in domtblout_patterns:
                    if pattern.exists():
                        with open(pattern, 'r') as inf:
                            for line in inf:
                                if not line.startswith('#'):
                                    outf.write(line)
                        break
        
        print(f"  Saved: {merged_domtblout}")
        
        # Get class distribution after domain extraction
        if labels_df is not None:
            domain_ids = set(h['uniprot_id'] for h in best_hits)
            domain_labels = labels_df[labels_df['uniprot_id'].isin(domain_ids)]
            class_counts = domain_labels['label_used_for_experiments'].value_counts().to_dict()
            
            # Show class distribution
            print(f"\n  Class distribution after domain extraction:")
            for cls, count in sorted(class_counts.items(), key=lambda x: -x[1]):
                print(f"    {cls:15} {count:>5}")
        else:
            class_counts = {}
        
        # Store in report
        report["results"][f"E{e_str}"] = {
            "e_value_threshold": e_thresh,
            "n_domains_extracted": len(best_hits),
            "n_unique_proteins": len(set(h['uniprot_id'] for h in best_hits)),
            "mean_domain_length": sum(h['domain_length'] for h in best_hits) / len(best_hits) if best_hits else 0,
            "class_distribution": {k: int(v) for k, v in class_counts.items()},
            "files": {
                "domtblout": str(merged_domtblout),
                "fasta": str(fasta_file),
                "coords_tsv": str(coords_file)
            }
        }
    
    # Save report
    report_file = domains_dir / "domain_extraction_report.json"
    with open(report_file, 'w') as f:
        json.dump(report, f, indent=2)
    
    print(f"\n{'='*60}")
    print("STEP 7 COMPLETE: Domain Extraction")
    print(f"{'='*60}")
    print(f"\nReport saved to: {report_file}")
    
    # Summary table
    print(f"\n{'Summary':^60}")
    print("-" * 60)
    print(f"{'E-value':<15} {'Domains':<15} {'Mean Length':<15}")
    print("-" * 60)
    for e_key, data in report["results"].items():
        print(f"{data['e_value_threshold']:<15} {data['n_domains_extracted']:<15} {data['mean_domain_length']:.1f}")
    print("-" * 60)
    
    # Sanity checks
    print("\nSanity checks:")
    for e_key, data in report["results"].items():
        fasta_file = Path(data['files']['fasta'])
        coords_file = Path(data['files']['coords_tsv'])
        
        # Count FASTA entries
        fasta_count = sum(1 for line in open(fasta_file) if line.startswith('>'))
        # Count coords rows
        coords_count = len(pd.read_csv(coords_file, sep='\t'))
        
        match = "✓" if fasta_count == coords_count else "✗"
        print(f"  {match} {e_key}: FASTA ({fasta_count}) == coords TSV ({coords_count})")


if __name__ == "__main__":
    main()

