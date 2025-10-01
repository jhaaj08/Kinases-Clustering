"""
Script to download all kinase sequences from UniProt and save as CSV.
"""

import requests
import csv
import time
from typing import Optional
import re


def extract_kinome_group(protein_families: str, keywords: str, protein_name: str) -> str:
    """
    Extract kinome group/subfamily information from UniProt fields.
    
    This is a basic extraction from available UniProt data. For complete
    kinome classification, integrate with Manning kinome classification or
    specialized kinase databases.
    """
    # Combine available information
    combined = f"{protein_families} {keywords} {protein_name}".upper()
    
    # Common kinase group patterns (based on Manning classification)
    kinase_groups = {
        'AGC': ['PKA', 'PKC', 'PKG', 'AKT', 'SGK', 'RSK', 'ROCK', 'GRK'],
        'CAMK': ['CAMK', 'CALCIUM', 'CALMODULIN', 'AMPK', 'MARK', 'MELK'],
        'CK1': ['CK1', 'CASEIN KINASE 1', 'CSNK1'],
        'CMGC': ['CDK', 'MAPK', 'GSK', 'CLK', 'DYRK', 'ERK', 'JNK'],
        'STE': ['STE7', 'STE11', 'STE20', 'MAP2K', 'MAP3K', 'PAK'],
        'TK': ['TYROSINE KINASE', 'RECEPTOR TYROSINE', 'SRC', 'ABL', 'EGFR', 'PDGFR'],
        'TKL': ['TGF-BETA', 'TGFBR', 'IRAK', 'LRRK', 'MLK'],
        'RGC': ['GUANYLATE CYCLASE', 'GUANYLYL CYCLASE'],
        'Atypical': ['PI3K', 'PIKK', 'ALPHA-KINASE', 'RIO'],
        'Histidine': ['HISTIDINE KINASE', 'TWO-COMPONENT']
    }
    
    # Try to match kinase groups
    for group, patterns in kinase_groups.items():
        for pattern in patterns:
            if pattern in combined:
                # Try to get more specific subfamily info
                if 'SUBFAMILY' in protein_families or 'FAMILY' in protein_families:
                    return f"{group} / {protein_families[:50]}"
                return f"{group}"
    
    # If no specific group found, return family info
    if protein_families and protein_families != "N/A":
        return protein_families[:100]
    
    return "Unclassified"


def download_all_kinases_from_uniprot(output_file: str = "kinases.csv", 
                                       organism: Optional[str] = None) -> dict:
    """
    Download ALL kinase sequences from UniProt and save as CSV.
    
    CSV Columns:
    - uniprot_id: UniProt accession number
    - protein_name: Protein name/description
    - function: Function annotation text from UniProt
    - kinome_group_subfamily: Kinase group/subfamily (extracted from UniProt fields)
    - conformation_DFG_aC: DFG/αC conformation state (placeholder - requires KLIFS/PDB data)
    - inhibitor_class_sensitivity: Inhibitor class sensitivity (placeholder - requires experimental data)
    - sequence: Amino acid sequence
    
    Note: Conformation and inhibitor sensitivity fields are placeholders.
    For complete data, integrate with specialized databases like KLIFS, ChEMBL, or KIDFamMap.
    
    Parameters:
    -----------
    output_file : str
        Output filename for the CSV (default: "kinases.csv")
    organism : str, optional
        Organism name to filter results (e.g., "Homo sapiens", "human")
        If None, downloads from all organisms
    
    Returns:
    --------
    dict
        Dictionary with download statistics including:
        - 'sequences_downloaded': number of sequences downloaded
        - 'output_file': path to output file
        - 'query': UniProt query used
    """
    
    base_url = "https://rest.uniprot.org/uniprotkb/stream"
    
    # Construct query for kinases
    # Search for proteins with "kinase" in their name/description
    query = "(protein_name:kinase) AND (reviewed:true)"
    
    # Add organism filter if specified
    if organism:
        query += f" AND (organism_name:{organism})"
    
    # Parameters for the API request
    # Using TSV format to easily parse the data
    params = {
        'query': query,
        'format': 'tsv',
        'fields': 'accession,protein_name,cc_function,protein_families,keyword,sequence'
    }
    
    print(f"Querying UniProt for ALL kinase sequences...")
    print(f"Query: {query}")
    print("This may take a while...")
    
    try:
        # Make the request - using stream endpoint to get all results
        response = requests.get(base_url, params=params, stream=True)
        response.raise_for_status()
        
        # Parse TSV data
        lines = response.text.strip().split('\n')
        
        if len(lines) <= 1:
            print("Warning: No sequences found matching the query.")
            return {
                'sequences_downloaded': 0,
                'output_file': output_file,
                'query': query
            }
        
        # Write to CSV
        csv_data = []
        header = lines[0].split('\t')
        
        print(f"Processing {len(lines)-1} kinase sequences...")
        
        for i, line in enumerate(lines[1:], 1):
            if i % 100 == 0:
                print(f"  Processed {i} sequences...")
            
            parts = line.split('\t')
            if len(parts) >= 6:
                uniprot_id = parts[0]
                protein_name = parts[1]
                function = parts[2] if parts[2] else "N/A"
                protein_families = parts[3] if parts[3] else "N/A"
                keywords = parts[4] if parts[4] else "N/A"
                sequence = parts[5]
                
                # Extract kinome group/subfamily from protein_families and keywords
                kinome_group = extract_kinome_group(protein_families, keywords, protein_name)
                
                # Placeholders for data that requires external databases
                # These would need to be populated from KLIFS, PDB, ChEMBL, etc.
                conformation = "Not available"  # Requires structural data
                inhibitor_sensitivity = "Not available"  # Requires experimental data
                
                csv_data.append({
                    'uniprot_id': uniprot_id,
                    'protein_name': protein_name,
                    'function': function,
                    'kinome_group_subfamily': kinome_group,
                    'conformation_DFG_aC': conformation,
                    'inhibitor_class_sensitivity': inhibitor_sensitivity,
                    'sequence': sequence
                })
        
        # Save to CSV file
        with open(output_file, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=[
                'uniprot_id', 
                'protein_name', 
                'function', 
                'kinome_group_subfamily',
                'conformation_DFG_aC',
                'inhibitor_class_sensitivity',
                'sequence'
            ])
            writer.writeheader()
            writer.writerows(csv_data)
        
        num_sequences = len(csv_data)
        print(f"\nSuccessfully downloaded {num_sequences} kinase sequences")
        print(f"Saved to: {output_file}")
        
        return {
            'sequences_downloaded': num_sequences,
            'output_file': output_file,
            'query': query
        }
        
    except requests.exceptions.RequestException as e:
        print(f"Error downloading from UniProt: {e}")
        return {
            'sequences_downloaded': 0,
            'output_file': output_file,
            'query': query,
            'error': str(e)
        }


def main():
    """
    Main function to download all kinases.
    """
    print("=" * 70)
    print("Kinase Sequence Downloader - CSV Format")
    print("=" * 70)
    print()
    
    # Download all kinases
    result = download_all_kinases_from_uniprot(
        output_file="kinases_all.csv",
        organism=None  # Set to "Homo sapiens" to filter for human kinases only
    )
    
    print()
    print("Download Summary:")
    print(f"  - Downloaded: {result['sequences_downloaded']} sequences")
    print(f"  - Output file: {result['output_file']}")
    print(f"  - Format: CSV (uniprot_id, protein_name, function, kinome_group_subfamily,")
    print(f"              conformation_DFG_aC, inhibitor_class_sensitivity, sequence)")
    print()
    print("Note: conformation_DFG_aC and inhibitor_class_sensitivity are placeholders.")
    print("      For complete data, integrate with KLIFS, PDB, ChEMBL, or KIDFamMap databases.")
    
    # Example: Download only human kinases
    # result2 = download_all_kinases_from_uniprot(
    #     output_file="human_kinases_all.csv",
    #     organism="Homo sapiens"
    # )


if __name__ == "__main__":
    main()