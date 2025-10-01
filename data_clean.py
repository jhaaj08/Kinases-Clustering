"""
Data cleaning script for kinases dataset.
Reduces label cardinality, handles duplicates, and applies CD-HIT clustering.
"""

import pandas as pd
import numpy as np
import subprocess
import os
import tempfile
from pathlib import Path


def map_to_major_group(subfamily: str) -> str:
    """
    Map kinome_group_subfamily to major kinase groups based on Manning classification.
    
    Parameters:
    -----------
    subfamily : str
        The subfamily classification from UniProt
    
    Returns:
    --------
    str
        Major kinase group (AGC, CAMK, CK1, CMGC, STE, TK, TKL, RGC, Atypical, Histidine, Other)
    """
    if pd.isna(subfamily) or subfamily == 'Unclassified':
        return 'Other'
    
    s = str(subfamily).upper()
    
    # Manning kinase groups
    if 'AGC' in s or any(x in s for x in ['PKA', 'PKC', 'PKG', 'AKT', 'SGK', 'RSK', 'ROCK', 'GRK']):
        return 'AGC'
    
    if 'CAMK' in s or any(x in s for x in ['CALCIUM', 'CALMODULIN', 'AMPK', 'MARK', 'MELK', 'DAPK']):
        return 'CAMK'
    
    if 'CK1' in s or 'CASEIN KINASE 1' in s or 'CSNK1' in s:
        return 'CK1'
    
    if 'CMGC' in s or any(x in s for x in ['CDK', 'MAPK', 'GSK', 'CLK', 'DYRK', 'ERK', 'JNK']):
        return 'CMGC'
    
    if 'STE' in s or any(x in s for x in ['MAP2K', 'MAP3K', 'PAK', 'STE7', 'STE11', 'STE20']):
        return 'STE'
    
    if s.startswith('TK ') or s == 'TK' or any(x in s for x in ['TYROSINE KINASE', 'RECEPTOR TYROSINE', 'SRC', 'ABL', 'EGFR', 'PDGFR']):
        return 'TK'
    
    if 'TKL' in s or any(x in s for x in ['TGF-BETA', 'TGFBR', 'IRAK', 'LRRK', 'MLK']):
        return 'TKL'
    
    if 'RGC' in s or any(x in s for x in ['GUANYLATE CYCLASE', 'GUANYLYL CYCLASE']):
        return 'RGC'
    
    if 'ATYPICAL' in s or any(x in s for x in ['PI3K', 'PIKK', 'ALPHA-KINASE', 'RIO', 'MTOR']):
        return 'Atypical'
    
    if 'HISTIDINE' in s or 'TWO-COMPONENT' in s:
        return 'Histidine'
    
    # Metabolic and other kinases
    return 'Other'


def create_family_slim(df: pd.DataFrame, min_count: int = 50) -> pd.Series:
    """
    Create a slim family classification by grouping rare families.
    
    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame with kinome_group_subfamily and kinome_group_major columns
    min_count : int
        Minimum count for a family to be kept separate (default: 50)
    
    Returns:
    --------
    pd.Series
        Slim family labels
    """
    counts = df['kinome_group_subfamily'].value_counts()
    keep_families = set(counts[counts >= min_count].index)
    
    def slim_mapper(row):
        subfamily = row['kinome_group_subfamily']
        major = row['kinome_group_major']
        
        if subfamily in keep_families:
            return subfamily
        else:
            return f"{major}:Other"
    
    return df.apply(slim_mapper, axis=1)


def check_duplicates(df: pd.DataFrame) -> tuple:
    """
    Check for duplicate sequences in the dataset.
    
    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame with sequence column
    
    Returns:
    --------
    tuple
        (has_duplicates, duplicate_df, num_duplicate_groups)
    """
    # Normalize sequences for comparison
    seq_normalized = df['sequence'].astype(str).str.strip().str.upper()
    
    # Find duplicates
    is_duplicate = seq_normalized.duplicated(keep=False)
    duplicate_df = df[is_duplicate].copy()
    
    # Count unique duplicate sequences
    num_duplicate_groups = seq_normalized[seq_normalized.duplicated()].nunique()
    
    has_duplicates = len(duplicate_df) > 0
    
    return has_duplicates, duplicate_df, num_duplicate_groups


def check_cdhit_installed() -> bool:
    """
    Check if CD-HIT is installed and accessible.
    
    Returns:
    --------
    bool
        True if CD-HIT is available, False otherwise
    """
    try:
        result = subprocess.run(['cd-hit', '-h'], 
                              capture_output=True, 
                              timeout=5)
        return True
    except (subprocess.TimeoutExpired, FileNotFoundError, subprocess.CalledProcessError):
        return False


def run_cdhit_clustering(df: pd.DataFrame, identity_threshold: float = 0.6) -> pd.DataFrame:
    """
    Cluster sequences using CD-HIT at specified identity threshold.
    
    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame with sequences
    identity_threshold : float
        Sequence identity threshold (0.0 to 1.0), default 0.6 (60%)
    
    Returns:
    --------
    pd.DataFrame
        DataFrame with representative sequences only
    """
    print(f"Running CD-HIT clustering at {identity_threshold*100:.0f}% identity...")
    
    # Create temporary directory for CD-HIT files
    with tempfile.TemporaryDirectory() as tmpdir:
        input_fasta = os.path.join(tmpdir, 'input.fasta')
        output_prefix = os.path.join(tmpdir, 'output')
        
        # Write sequences to FASTA file
        print("  Writing sequences to FASTA...")
        with open(input_fasta, 'w') as f:
            for idx, row in df.iterrows():
                # Use uniprot_id as sequence identifier
                f.write(f">{row['uniprot_id']}\n")
                f.write(f"{row['sequence']}\n")
        
        # Run CD-HIT
        print(f"  Running CD-HIT (this may take a few minutes)...")
        cmd = [
            'cd-hit',
            '-i', input_fasta,
            '-o', output_prefix,
            '-c', str(identity_threshold),  # sequence identity threshold
            '-n', '3',                       # word size (3 for ~60% identity)
            '-M', '0',                       # unlimited memory
            '-T', '0',                       # use all CPUs
            '-d', '0'                        # full sequence description
        ]
        
        try:
            result = subprocess.run(cmd, 
                                  capture_output=True, 
                                  text=True, 
                                  timeout=600)  # 10 min timeout
            
            if result.returncode != 0:
                print(f"  ⚠️  CD-HIT error: {result.stderr}")
                print("  Falling back to duplicate removal only...")
                return df
            
            # Read cluster file to get representatives
            cluster_file = f"{output_prefix}.clstr"
            representatives = set()
            
            with open(cluster_file, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        continue
                    if '*' in line:  # Representative sequence
                        # Extract UniProt ID from line like: "0	813aa, >A0A075F7E9... *"
                        parts = line.split('>')
                        if len(parts) > 1:
                            uniprot_id = parts[1].split('.')[0]
                            representatives.add(uniprot_id)
            
            # Filter dataframe to keep only representatives
            df_filtered = df[df['uniprot_id'].isin(representatives)].copy()
            
            num_clusters = len(representatives)
            reduction = len(df) - num_clusters
            
            print(f"  ✅ CD-HIT clustering complete")
            print(f"     Input sequences:  {len(df):,}")
            print(f"     Clusters found:   {num_clusters:,}")
            print(f"     Sequences removed: {reduction:,} ({reduction/len(df)*100:.1f}%)")
            
            return df_filtered
            
        except subprocess.TimeoutExpired:
            print("  ⚠️  CD-HIT timed out. Falling back to duplicate removal only...")
            return df
        except Exception as e:
            print(f"  ⚠️  CD-HIT error: {str(e)}")
            print("  Falling back to duplicate removal only...")
            return df


def data_clean(input_file: str = 'kinases_all.csv', 
               output_file: str = 'kinases_revised.csv',
               min_family_count: int = 50,
               remove_duplicates: bool = True,
               use_cdhit: bool = True,
               cdhit_identity: float = 0.6) -> dict:
    """
    Clean kinases dataset by reducing label cardinality and handling duplicates.
    
    Parameters:
    -----------
    input_file : str
        Input CSV file path (default: 'kinases_all.csv')
    output_file : str
        Output CSV file path (default: 'kinases_revised.csv')
    min_family_count : int
        Minimum count for families in family_slim (default: 50)
    remove_duplicates : bool
        Whether to remove duplicate sequences (default: True)
    use_cdhit : bool
        Whether to use CD-HIT for sequence clustering (default: True)
    cdhit_identity : float
        CD-HIT sequence identity threshold 0.0-1.0 (default: 0.6 for 60%)
    
    Returns:
    --------
    dict
        Statistics about the cleaning process
    """
    
    print("=" * 80)
    print("KINASES DATA CLEANING")
    print("=" * 80)
    print()
    
    # Read input file
    print(f"Reading {input_file}...")
    df = pd.read_csv(input_file)
    original_count = len(df)
    print(f"  Loaded {original_count:,} kinases")
    print()
    
    # Check for duplicates
    print("Checking for duplicate sequences...")
    has_dups, dup_df, num_dup_groups = check_duplicates(df)
    
    if has_dups:
        print(f"  ⚠️  Found {len(dup_df):,} rows with duplicate sequences")
        print(f"  ⚠️  {num_dup_groups:,} unique sequences are duplicated")
        
        # Save duplicates for review
        dup_file = 'duplicates_report.csv'
        dup_df.to_csv(dup_file, index=False)
        print(f"  Saved duplicates to {dup_file}")
        
        if remove_duplicates:
            # Keep first occurrence of each sequence
            seq_normalized = df['sequence'].astype(str).str.strip().str.upper()
            df = df[~seq_normalized.duplicated(keep='first')].copy()
            print(f"  ✅ Removed duplicates, keeping first occurrence")
            print(f"  New count: {len(df):,} kinases")
    else:
        print("  ✅ No duplicate sequences found")
    print()
    
    # Apply CD-HIT clustering if requested
    count_before_cdhit = len(df)
    if use_cdhit:
        print(f"CD-HIT Clustering (identity threshold: {cdhit_identity*100:.0f}%)...")
        
        # Check if CD-HIT is installed
        if check_cdhit_installed():
            print("  ✅ CD-HIT is installed")
            df = run_cdhit_clustering(df, identity_threshold=cdhit_identity)
        else:
            print("  ⚠️  CD-HIT not found. Please install CD-HIT:")
            print("     macOS:   brew install cd-hit")
            print("     Ubuntu:  sudo apt-get install cd-hit")
            print("     conda:   conda install -c bioconda cd-hit")
            print("  Skipping CD-HIT clustering...")
        print()
    
    # Add major group classification
    print("Creating major group classification...")
    df['kinome_group_major'] = df['kinome_group_subfamily'].apply(map_to_major_group)
    
    major_group_counts = df['kinome_group_major'].value_counts()
    print(f"  ✅ Mapped to {len(major_group_counts)} major groups:")
    for group, count in major_group_counts.items():
        print(f"     {group:15s}: {count:5,} ({count/len(df)*100:5.1f}%)")
    print()
    
    # Create slim family classification
    print(f"Creating slim family classification (min_count={min_family_count})...")
    original_families = df['kinome_group_subfamily'].nunique()
    
    df['family_slim'] = create_family_slim(df, min_count=min_family_count)
    
    slim_families = df['family_slim'].nunique()
    print(f"  ✅ Reduced from {original_families} to {slim_families} families")
    print(f"  Reduction: {(1 - slim_families/original_families)*100:.1f}%")
    print()
    
    # Show distribution of slim families
    print("Top 15 slim families:")
    family_counts = df['family_slim'].value_counts().head(15)
    for family, count in family_counts.items():
        print(f"  {family:50s}: {count:5,} ({count/len(df)*100:5.1f}%)")
    print()
    
    # Reorder columns for better readability
    column_order = [
        'uniprot_id',
        'protein_name',
        'function',
        'kinome_group_subfamily',
        'kinome_group_major',
        'family_slim',
        'conformation_DFG_aC',
        'inhibitor_class_sensitivity',
        'sequence'
    ]
    df = df[column_order]
    
    # Save cleaned data
    print(f"Saving cleaned data to {output_file}...")
    df.to_csv(output_file, index=False)
    print(f"  ✅ Saved {len(df):,} kinases")
    print()
    
    # Summary statistics
    cdhit_removed = count_before_cdhit - len(df) if use_cdhit else 0
    total_removed = original_count - len(df)
    duplicates_only = total_removed - cdhit_removed
    
    stats = {
        'original_count': original_count,
        'final_count': len(df),
        'duplicates_removed': duplicates_only,
        'cdhit_removed': cdhit_removed,
        'total_removed': total_removed,
        'original_families': original_families,
        'slim_families': slim_families,
        'major_groups': len(major_group_counts),
        'cdhit_used': use_cdhit and check_cdhit_installed(),
        'cdhit_identity': cdhit_identity if use_cdhit else None,
        'output_file': output_file
    }
    
    print("=" * 80)
    print("CLEANING SUMMARY")
    print("=" * 80)
    print(f"  Original kinases:        {stats['original_count']:,}")
    print(f"  Final kinases:           {stats['final_count']:,}")
    print(f"  Total removed:           {stats['total_removed']:,} ({stats['total_removed']/stats['original_count']*100:.1f}%)")
    print(f"    - Exact duplicates:    {stats['duplicates_removed']:,}")
    if stats['cdhit_used']:
        print(f"    - CD-HIT clustering:   {stats['cdhit_removed']:,} (at {cdhit_identity*100:.0f}% identity)")
    print(f"  Original families:       {stats['original_families']:,}")
    print(f"  Slim families:           {stats['slim_families']:,}")
    print(f"  Major groups:            {stats['major_groups']:,}")
    print(f"  Output file:             {stats['output_file']}")
    print("=" * 80)
    print()
    
    return stats


def main():
    """
    Main function to run data cleaning with CD-HIT at 60% identity.
    """
    stats = data_clean(
        input_file='kinases_all.csv',
        output_file='kinases_revised.csv',
        min_family_count=50,
        remove_duplicates=True,
        use_cdhit=True,
        cdhit_identity=0.6  # 60% sequence identity threshold
    )
    
    print("✅ Data cleaning complete!")
    print(f"   Use '{stats['output_file']}' for downstream analysis")
    if stats['cdhit_used']:
        print(f"   CD-HIT clustering applied at {stats['cdhit_identity']*100:.0f}% identity")
    print()


if __name__ == "__main__":
    main()
