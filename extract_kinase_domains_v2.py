#!/usr/bin/env python3
"""
Enhanced kinase domain extraction with multiple HMMs and E-value options.

New features:
- Support for multiple Pfam profiles (PF00069 + PF07714)
- Configurable E-value thresholds
- Command-line arguments
- Merging results from multiple HMMs
"""

import os
import sys
import argparse
import subprocess
import tempfile
import pandas as pd
import requests
import gzip
from pathlib import Path


PFAM_PROFILES = {
    'PF00069': 'Pkinase',        # Protein kinase domain (Ser/Thr/Tyr)
    'PF07714': 'Pkinase_Tyr',    # Protein tyrosine kinase
}


def download_pfam_hmm(pfam_id, output_file):
    """Download Pfam HMM profile."""
    url = f"https://www.ebi.ac.uk/interpro/api/entry/pfam/{pfam_id}?annotation=hmm"
    
    print(f"Downloading {pfam_id} ({PFAM_PROFILES.get(pfam_id, 'unknown')})...")
    
    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        
        # Decompress gzip content
        hmm_content = gzip.decompress(response.content).decode('utf-8')
        
        with open(output_file, 'w') as f:
            f.write(hmm_content)
        
        print(f"  ✅ Downloaded to {output_file}")
        return output_file
    
    except Exception as e:
        print(f"  ❌ Error: {e}")
        return None


def check_hmmer_installed():
    """Check if HMMER is installed."""
    try:
        subprocess.run(['hmmsearch', '-h'], 
                      capture_output=True, 
                      timeout=5)
        return True
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


def run_hmmsearch(hmm_file, fasta_file, output_file, e_value=0.001):
    """Run HMMER hmmsearch."""
    cmd = [
        'hmmsearch',
        '--domtblout', output_file,
        '-E', str(e_value),
        hmm_file,
        fasta_file
    ]
    
    try:
        result = subprocess.run(cmd, 
                              capture_output=True, 
                              text=True,
                              timeout=600)
        
        if result.returncode != 0:
            print(f"    ❌ HMMER error: {result.stderr}")
            return None
        
        return output_file
    
    except Exception as e:
        print(f"    ❌ Error: {e}")
        return None


def parse_hmmer_output(hmmer_output, pfam_id):
    """Parse HMMER domtblout format."""
    domains = []
    
    if not os.path.exists(hmmer_output):
        return pd.DataFrame()
    
    with open(hmmer_output, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.split()
            if len(fields) < 23:
                continue
            
            domains.append({
                'uniprot_id': fields[0],
                'pfam_id': pfam_id,
                'pfam_name': PFAM_PROFILES.get(pfam_id, 'unknown'),
                'evalue': float(fields[12]),
                'score': float(fields[13]),
                'ali_from': int(fields[17]),
                'ali_to': int(fields[18]),
                'env_from': int(fields[19]),
                'env_to': int(fields[20]),
            })
    
    return pd.DataFrame(domains)


def merge_hmm_results(results_list):
    """
    Merge results from multiple HMM searches.
    Keep best domain per sequence (lowest E-value).
    """
    if not results_list:
        return pd.DataFrame()
    
    # Concatenate all results
    all_domains = pd.concat(results_list, ignore_index=True)
    
    if len(all_domains) == 0:
        return all_domains
    
    # Sort by E-value and score
    all_domains = all_domains.sort_values(['uniprot_id', 'evalue', 'score'], 
                                         ascending=[True, True, False])
    
    # Keep best domain per sequence
    best_domains = all_domains.groupby('uniprot_id').first().reset_index()
    
    return best_domains


def extract_domain_sequences(kinases_df, domains_df, use_envelope=True):
    """Extract domain sequences from full sequences."""
    # Merge to get sequences
    merged = domains_df.merge(kinases_df[['uniprot_id', 'sequence']], 
                             on='uniprot_id', 
                             how='left')
    
    # Extract domain subsequence
    def extract_subseq(row):
        if use_envelope:
            start = row['env_from'] - 1
            end = row['env_to']
        else:
            start = row['ali_from'] - 1
            end = row['ali_to']
        
        return row['sequence'][start:end]
    
    merged['domain_sequence'] = merged.apply(extract_subseq, axis=1)
    merged['domain_length'] = merged['domain_sequence'].apply(len)
    
    return merged


def save_domain_dataset(domains_df, kinases_df, output_file):
    """Save domain-only dataset."""
    # Start with original data
    domain_data = kinases_df[kinases_df['uniprot_id'].isin(domains_df['uniprot_id'])].copy()
    
    # Replace full sequence with domain sequence
    domain_map = dict(zip(domains_df['uniprot_id'], domains_df['domain_sequence']))
    domain_data['sequence'] = domain_data['uniprot_id'].map(domain_map)
    
    # Add domain metadata
    meta_cols = ['uniprot_id', 'pfam_id', 'pfam_name', 'env_from', 'env_to', 
                 'domain_length', 'evalue', 'score']
    domain_meta = domains_df[meta_cols]
    domain_data = domain_data.merge(domain_meta, on='uniprot_id', how='left')
    
    # Save
    domain_data.to_csv(output_file, index=False)
    
    return domain_data


def main():
    parser = argparse.ArgumentParser(
        description='Extract kinase domains with multiple HMMs and E-value options'
    )
    parser.add_argument(
        '--input',
        default='kinases_revised.csv',
        help='Input CSV file with kinase sequences'
    )
    parser.add_argument(
        '--output',
        default='kinases_domains_v2.csv',
        help='Output CSV file for domain sequences'
    )
    parser.add_argument(
        '--hmms',
        nargs='+',
        default=['PF00069', 'PF07714'],
        choices=['PF00069', 'PF07714'],
        help='Pfam HMM profiles to use (default: PF00069 PF07714)'
    )
    parser.add_argument(
        '--evalue',
        type=float,
        default=0.01,
        help='E-value threshold (default: 0.01)'
    )
    parser.add_argument(
        '--envelope',
        action='store_true',
        default=True,
        help='Use envelope boundaries (default: True)'
    )
    
    args = parser.parse_args()
    
    print("="*80)
    print("ENHANCED KINASE DOMAIN EXTRACTION")
    print("="*80)
    print()
    print(f"Input:    {args.input}")
    print(f"Output:   {args.output}")
    print(f"HMMs:     {', '.join(args.hmms)}")
    print(f"E-value:  {args.evalue}")
    print(f"Envelope: {args.envelope}")
    print()
    
    # Check HMMER
    if not check_hmmer_installed():
        print("❌ HMMER not installed!")
        sys.exit(1)
    print("✅ HMMER is installed")
    
    # Load kinases
    print(f"\nLoading {args.input}...")
    df = pd.read_csv(args.input)
    print(f"✅ Loaded {len(df):,} kinases")
    
    # Download HMMs
    hmm_files = {}
    for pfam_id in args.hmms:
        hmm_file = f"{pfam_id}.hmm"
        if not os.path.exists(hmm_file):
            result = download_pfam_hmm(pfam_id, hmm_file)
            if result:
                hmm_files[pfam_id] = hmm_file
        else:
            print(f"✅ Using existing HMM: {hmm_file}")
            hmm_files[pfam_id] = hmm_file
    
    if not hmm_files:
        print("❌ No HMM files available!")
        sys.exit(1)
    
    # Create temporary FASTA
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        fasta_file = f.name
        for _, row in df.iterrows():
            f.write(f">{row['uniprot_id']}\n")
            f.write(f"{row['sequence']}\n")
    
    print(f"✅ Created FASTA: {fasta_file}")
    
    # Run HMMER for each HMM
    print(f"\nRunning HMMER searches (E-value ≤ {args.evalue})...")
    all_results = []
    
    for pfam_id, hmm_file in hmm_files.items():
        print(f"\n  {pfam_id} ({PFAM_PROFILES[pfam_id]})...")
        output_file = f"hmmer_{pfam_id}_e{args.evalue}.domtblout"
        
        result = run_hmmsearch(hmm_file, fasta_file, output_file, args.evalue)
        
        if result:
            domains_df = parse_hmmer_output(output_file, pfam_id)
            if len(domains_df) > 0:
                print(f"    ✅ Found {len(domains_df)} hits in {domains_df['uniprot_id'].nunique()} sequences")
                all_results.append(domains_df)
            else:
                print(f"    ⚠️  No hits found")
    
    # Merge results
    print(f"\nMerging results from {len(all_results)} HMM searches...")
    merged_domains = merge_hmm_results(all_results)
    
    if len(merged_domains) == 0:
        print("❌ No domains found!")
        sys.exit(1)
    
    print(f"✅ Total unique domains: {len(merged_domains):,}")
    
    # Show HMM distribution
    print("\nDomain source distribution:")
    for pfam_id, count in merged_domains['pfam_id'].value_counts().items():
        pct = count / len(merged_domains) * 100
        print(f"  {pfam_id} ({PFAM_PROFILES[pfam_id]}): {count:,} ({pct:.1f}%)")
    
    # Extract sequences
    print("\nExtracting domain sequences...")
    domains_with_seq = extract_domain_sequences(df, merged_domains, args.envelope)
    print(f"✅ Extracted {len(domains_with_seq):,} domain sequences")
    print(f"   Length: {domains_with_seq['domain_length'].mean():.0f} ± {domains_with_seq['domain_length'].std():.0f} aa")
    
    # Save
    print(f"\nSaving to {args.output}...")
    domain_data = save_domain_dataset(domains_with_seq, df, args.output)
    print(f"✅ Saved {len(domain_data):,} sequences")
    
    # Cleanup
    os.unlink(fasta_file)
    
    # Statistics
    print("\n" + "="*80)
    print("EXTRACTION SUMMARY")
    print("="*80)
    print(f"Input sequences:       {len(df):,}")
    print(f"Domains found:         {len(merged_domains):,}")
    print(f"Success rate:          {len(merged_domains)/len(df)*100:.1f}%")
    print(f"E-value threshold:     {args.evalue}")
    print(f"Output file:           {args.output}")
    print()
    print(f"Mean domain length:    {domains_with_seq['domain_length'].mean():.0f} aa")
    print(f"Mean full length:      {df['sequence'].apply(len).mean():.0f} aa")
    print(f"Reduction:             {(1 - domains_with_seq['domain_length'].mean()/df['sequence'].apply(len).mean())*100:.1f}%")
    print()
    
    # Distribution by kinome group
    if 'kinome_group_major' in domain_data.columns:
        print("Domains per kinome group:")
        group_counts = domain_data['kinome_group_major'].value_counts()
        for group, count in group_counts.items():
            pct = count / len(domain_data) * 100
            print(f"  {group:12s}: {count:4d} ({pct:5.1f}%)")
    
    print("\n" + "="*80)
    print("✅ EXTRACTION COMPLETE!")
    print("="*80)
    print()


if __name__ == "__main__":
    main()

