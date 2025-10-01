"""
Quick script to verify and show statistics about the downloaded kinase data.
"""

import csv
import pandas as pd


def verify_kinases_csv(csv_file: str = "kinases_all.csv"):
    """
    Verify and display statistics about the kinases CSV file.
    
    Parameters:
    -----------
    csv_file : str
        Path to the CSV file to verify
    """
    
    print("=" * 80)
    print("Kinase Data Verification")
    print("=" * 80)
    print()
    
    # Read CSV
    df = pd.read_csv(csv_file)
    
    # Basic statistics
    print(f"📊 Basic Statistics:")
    print(f"   Total kinases: {len(df):,}")
    print(f"   Columns: {list(df.columns)}")
    print()
    
    # Sequence length statistics
    df['seq_length'] = df['sequence'].apply(len)
    print(f"📏 Sequence Length Statistics:")
    print(f"   Minimum length: {df['seq_length'].min()} amino acids")
    print(f"   Maximum length: {df['seq_length'].max()} amino acids")
    print(f"   Average length: {df['seq_length'].mean():.1f} amino acids")
    print(f"   Median length: {df['seq_length'].median():.1f} amino acids")
    print()
    
    # Function annotation statistics
    has_function = df['function'].notna() & (df['function'] != 'N/A')
    print(f"🔬 Function Annotation Coverage:")
    print(f"   With annotations: {has_function.sum():,} ({has_function.sum()/len(df)*100:.1f}%)")
    print(f"   Without annotations (N/A): {(~has_function).sum():,} ({(~has_function).sum()/len(df)*100:.1f}%)")
    print()
    
    # Kinome group classification statistics
    is_classified = df['kinome_group_subfamily'] != 'Unclassified'
    print(f"🧬 Kinome Group Classification:")
    print(f"   Classified: {is_classified.sum():,} ({is_classified.sum()/len(df)*100:.1f}%)")
    print(f"   Unclassified: {(~is_classified).sum():,} ({(~is_classified).sum()/len(df)*100:.1f}%)")
    print()
    
    # Top kinase groups
    df['main_group'] = df['kinome_group_subfamily'].str.split(' / ').str[0].str.split().str[0]
    main_groups = df['main_group'].value_counts().head(10)
    print(f"   Top 10 Kinase Groups:")
    for group, count in main_groups.items():
        print(f"     {group:20s}: {count:5d} ({count/len(df)*100:4.1f}%)")
    print()
    
    # Conformation and inhibitor data availability
    has_conformation = df['conformation_DFG_aC'] != 'Not available'
    has_inhibitor = df['inhibitor_class_sensitivity'] != 'Not available'
    print(f"⚠️  External Data Fields (Placeholders):")
    print(f"   Conformation data available: {has_conformation.sum():,} ({has_conformation.sum()/len(df)*100:.1f}%)")
    print(f"   Inhibitor data available: {has_inhibitor.sum():,} ({has_inhibitor.sum()/len(df)*100:.1f}%)")
    print(f"   See EXTERNAL_DATA_INTEGRATION.md for population instructions")
    print()
    
    # Sample kinases
    print(f"📋 Sample Kinases (first 5):")
    for idx, row in df.head(5).iterrows():
        protein_name = row['protein_name']
        if len(protein_name) > 60:
            protein_name = protein_name[:60] + "..."
        print(f"   {idx+1}. {row['uniprot_id']:12s} - {protein_name}")
        print(f"      Group: {row['kinome_group_subfamily'][:40]}")
    
    print()
    print(f"✅ Data verification complete!")
    print(f"   File: {csv_file}")
    print(f"   Status: Valid CSV with {len(df):,} kinase sequences")
    print()


if __name__ == "__main__":
    verify_kinases_csv("kinases_all.csv")
