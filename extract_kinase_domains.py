#!/usr/bin/env python3
"""
Extract kinase catalytic domains using HMMER and Pfam HMM profile.

This script:
1. Downloads the Pfam Protein Kinase domain HMM (PF00069)
2. Runs HMMER to identify domain boundaries
3. Extracts domain sequences
4. Saves domain-only sequences for re-embedding
"""

import os
import sys
import subprocess
import tempfile
import pandas as pd
import requests
from pathlib import Path


def download_pfam_hmm(pfam_id="PF00069", output_file="Pkinase.hmm"):
    """
    Download Pfam HMM profile for protein kinase domain.
    
    Args:
        pfam_id: Pfam accession (default: PF00069 - Protein kinase domain)
        output_file: Output HMM file path
    
    Returns:
        Path to downloaded HMM file
    """
    import gzip
    
    url = f"https://www.ebi.ac.uk/interpro/api/entry/pfam/{pfam_id}?annotation=hmm"
    
    print(f"Downloading Pfam HMM for {pfam_id}...")
    print(f"URL: {url}")
    
    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        
        # Decompress gzip content
        hmm_content = gzip.decompress(response.content).decode('utf-8')
        
        with open(output_file, 'w') as f:
            f.write(hmm_content)
        
        print(f"✅ Downloaded and decompressed HMM to {output_file}")
        return output_file
    
    except Exception as e:
        print(f"❌ Error downloading HMM: {e}")
        print("\nAlternative: Download manually from:")
        print(f"https://www.ebi.ac.uk/interpro/entry/pfam/{pfam_id}/")
        sys.exit(1)


def check_hmmer_installed():
    """Check if HMMER is installed."""
    try:
        result = subprocess.run(['hmmsearch', '-h'], 
                              capture_output=True, 
                              timeout=5)
        return True
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


def run_hmmsearch(hmm_file, fasta_file, output_file, e_value=0.001):
    """
    Run HMMER hmmsearch to find kinase domains.
    
    Args:
        hmm_file: Path to HMM profile
        fasta_file: Path to input FASTA file
        output_file: Path to output file
        e_value: E-value threshold (default: 0.001)
    
    Returns:
        Path to output file
    """
    print(f"\nRunning HMMER hmmsearch...")
    print(f"  HMM: {hmm_file}")
    print(f"  Input: {fasta_file}")
    print(f"  E-value threshold: {e_value}")
    
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
            print(f"❌ HMMER error: {result.stderr}")
            sys.exit(1)
        
        print(f"✅ HMMER search complete")
        return output_file
    
    except subprocess.TimeoutExpired:
        print("❌ HMMER timed out (>10 min)")
        sys.exit(1)
    except Exception as e:
        print(f"❌ Error running HMMER: {e}")
        sys.exit(1)


def parse_hmmer_output(hmmer_output):
    """
    Parse HMMER domtblout format to extract domain boundaries.
    
    Returns:
        DataFrame with columns: uniprot_id, ali_from, ali_to, env_from, env_to, score, evalue
    """
    print("\nParsing HMMER results...")
    
    domains = []
    with open(hmmer_output, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.split()
            if len(fields) < 23:
                continue
            
            # Parse domtblout format
            # Columns: target, acc, tlen, query, acc, qlen, E-value, score, bias,
            #          #, of, c-Evalue, i-Evalue, score, bias, hmm from, hmm to,
            #          ali from, ali to, env from, env to, acc, description
            
            domains.append({
                'uniprot_id': fields[0],
                'evalue': float(fields[12]),
                'score': float(fields[13]),
                'ali_from': int(fields[17]),  # Alignment start in sequence
                'ali_to': int(fields[18]),    # Alignment end in sequence
                'env_from': int(fields[19]),  # Envelope start
                'env_to': int(fields[20]),    # Envelope end
            })
    
    df = pd.DataFrame(domains)
    
    if len(df) == 0:
        print("⚠️  No domains found!")
        return df
    
    print(f"✅ Found {len(df)} domain hits in {df['uniprot_id'].nunique()} sequences")
    
    # Keep best domain per sequence (lowest E-value, highest score)
    df = df.sort_values(['uniprot_id', 'evalue', 'score'], 
                       ascending=[True, True, False])
    df = df.groupby('uniprot_id').first().reset_index()
    
    print(f"   Kept {len(df)} best domains (one per sequence)")
    
    return df


def extract_domains(kinases_df, domains_df, use_envelope=True):
    """
    Extract domain sequences from full sequences.
    
    Args:
        kinases_df: DataFrame with 'uniprot_id' and 'sequence' columns
        domains_df: DataFrame from parse_hmmer_output
        use_envelope: Use envelope boundaries (more conservative) vs alignment boundaries
    
    Returns:
        DataFrame with domain-only sequences
    """
    print("\nExtracting domain sequences...")
    
    # Merge to get sequences
    merged = domains_df.merge(kinases_df[['uniprot_id', 'sequence']], 
                             on='uniprot_id', 
                             how='left')
    
    # Extract domain subsequence
    def extract_subseq(row):
        if use_envelope:
            start = row['env_from'] - 1  # Convert to 0-based
            end = row['env_to']
        else:
            start = row['ali_from'] - 1
            end = row['ali_to']
        
        return row['sequence'][start:end]
    
    merged['domain_sequence'] = merged.apply(extract_subseq, axis=1)
    merged['domain_length'] = merged['domain_sequence'].apply(len)
    
    print(f"✅ Extracted {len(merged)} domain sequences")
    print(f"   Domain length: {merged['domain_length'].mean():.0f} ± {merged['domain_length'].std():.0f} aa")
    print(f"   Range: {merged['domain_length'].min()}-{merged['domain_length'].max()} aa")
    
    return merged


def save_domain_sequences(domains_df, kinases_df, output_file):
    """
    Save domain-only dataset matching the format of kinases_revised.csv
    
    Args:
        domains_df: DataFrame with extracted domains
        kinases_df: Original kinases_revised.csv DataFrame
        output_file: Output CSV file path
    """
    print(f"\nCreating domain-only dataset...")
    
    # Start with original data
    domain_data = kinases_df[kinases_df['uniprot_id'].isin(domains_df['uniprot_id'])].copy()
    
    # Replace full sequence with domain sequence
    domain_map = dict(zip(domains_df['uniprot_id'], domains_df['domain_sequence']))
    domain_data['sequence'] = domain_data['uniprot_id'].map(domain_map)
    
    # Add domain metadata
    domain_meta = domains_df[['uniprot_id', 'env_from', 'env_to', 'domain_length', 'evalue', 'score']]
    domain_data = domain_data.merge(domain_meta, on='uniprot_id', how='left')
    
    # Save
    domain_data.to_csv(output_file, index=False)
    
    print(f"✅ Saved {len(domain_data)} domain sequences to {output_file}")
    print(f"\nColumns in output:")
    for col in domain_data.columns:
        print(f"  - {col}")
    
    return domain_data


def main():
    """Main workflow for domain extraction."""
    
    print("="*80)
    print("KINASE DOMAIN EXTRACTION PIPELINE")
    print("="*80)
    print()
    
    # Paths
    input_csv = "kinases_revised.csv"
    output_csv = "kinases_domains.csv"
    hmm_file = "Pkinase.hmm"
    
    # Check prerequisites
    if not check_hmmer_installed():
        print("❌ HMMER not installed!")
        print("   Install: conda install -c bioconda hmmer")
        sys.exit(1)
    
    print("✅ HMMER is installed")
    
    # Load kinases
    print(f"\nLoading {input_csv}...")
    df = pd.read_csv(input_csv)
    print(f"✅ Loaded {len(df)} kinases")
    
    # Download HMM if not present
    if not os.path.exists(hmm_file):
        download_pfam_hmm(output_file=hmm_file)
    else:
        print(f"✅ Using existing HMM: {hmm_file}")
    
    # Create temporary FASTA file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        fasta_file = f.name
        for _, row in df.iterrows():
            f.write(f">{row['uniprot_id']}\n")
            f.write(f"{row['sequence']}\n")
    
    print(f"✅ Created temporary FASTA: {fasta_file}")
    
    # Run HMMER
    hmmer_output = "hmmer_results.domtblout"
    run_hmmsearch(hmm_file, fasta_file, hmmer_output)
    
    # Parse results
    domains_df = parse_hmmer_output(hmmer_output)
    
    if len(domains_df) == 0:
        print("❌ No domains found! Check HMM file and input sequences.")
        sys.exit(1)
    
    # Extract domain sequences
    domains_with_seq = extract_domains(df, domains_df, use_envelope=True)
    
    # Save domain-only dataset
    domain_data = save_domain_sequences(domains_with_seq, df, output_csv)
    
    # Cleanup temporary files
    os.unlink(fasta_file)
    
    # Statistics
    print("\n" + "="*80)
    print("EXTRACTION SUMMARY")
    print("="*80)
    print(f"Input sequences:       {len(df):,}")
    print(f"Domains found:         {len(domains_df):,}")
    print(f"Success rate:          {len(domains_df)/len(df)*100:.1f}%")
    print(f"Output file:           {output_csv}")
    print()
    print(f"Mean domain length:    {domains_with_seq['domain_length'].mean():.0f} aa")
    print(f"Mean full length:      {df['sequence'].apply(len).mean():.0f} aa")
    print(f"Reduction:             {(1 - domains_with_seq['domain_length'].mean()/df['sequence'].apply(len).mean())*100:.1f}%")
    print()
    
    # Show distribution by kinome group
    if 'kinome_group_major' in domain_data.columns:
        print("Domains found per kinome group:")
        group_counts = domain_data['kinome_group_major'].value_counts()
        for group, count in group_counts.items():
            pct = count / len(domain_data) * 100
            print(f"  {group:12s}: {count:4d} ({pct:5.1f}%)")
    
    print("\n" + "="*80)
    print("✅ DOMAIN EXTRACTION COMPLETE!")
    print("="*80)
    print()
    print("Next steps:")
    print("  1. Re-embed domains: python generate_esm2_embeddings.py --input kinases_domains.csv --outdir kinases_domains_embeddings")
    print("  2. Re-run clustering and compare results")
    print()


if __name__ == "__main__":
    main()

