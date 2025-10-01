# External Data Integration Guide

This guide explains how to populate the **Conformation (DFG/αC states)** and **Inhibitor class sensitivity** columns with data from specialized kinase databases.

## Current CSV Columns

1. ✅ **uniprot_id** - Populated from UniProt
2. ✅ **protein_name** - Populated from UniProt
3. ✅ **function** - Populated from UniProt
4. ✅ **kinome_group_subfamily** - Extracted from UniProt fields (98.5% coverage)
5. ⚠️ **conformation_DFG_aC** - Placeholder (requires external data)
6. ⚠️ **inhibitor_class_sensitivity** - Placeholder (requires external data)
7. ✅ **sequence** - Populated from UniProt

## How to Populate Conformation Data (DFG/αC states)

### Background

Kinase conformations, particularly the **DFG (Asp-Phe-Gly) motif** and **αC helix** states, are crucial for:
- Understanding kinase activation mechanisms
- Drug design and inhibitor binding
- Structural biology studies

Common conformations:
- **DFG-in**: Active conformation, ATP binding competent
- **DFG-out**: Inactive conformation, opens hydrophobic pocket
- **αC-in/out**: Helix position affecting ATP binding

### Data Sources

#### 1. KLIFS Database (Recommended)
- **Website**: https://klifs.net/
- **API**: https://klifs.net/swagger/
- **Coverage**: ~7,000 kinase-ligand complexes
- **Data**: DFG conformation, αC helix position, detailed structural annotations

**Example API Usage:**
```python
import requests

def get_klifs_conformation(uniprot_id):
    """Get conformation data from KLIFS."""
    url = f"https://klifs.net/api/structures_pdb_list?kinase-uniprot={uniprot_id}"
    response = requests.get(url)
    
    if response.status_code == 200:
        data = response.json()
        if data:
            # Extract DFG and αC conformation from first structure
            structure = data[0]
            dfg = structure.get('dfg', 'Unknown')
            ac = structure.get('ac_helix', 'Unknown')
            return f"DFG-{dfg}/αC-{ac}"
    
    return "Not available"

# Example
conformation = get_klifs_conformation("P24941")  # CDK2
print(conformation)
```

#### 2. PDB (Protein Data Bank)
- **Website**: https://www.rcsb.org/
- **API**: https://data.rcsb.org/
- **Coverage**: ~200,000 structures (not all kinases)
- **Data**: Full 3D structures (requires manual/automated annotation)

#### 3. KinaMetrix
- **Website**: https://eidogen-sertanty.com/kinametrix.php
- **Coverage**: Commercial database with detailed kinase conformations
- **Data**: DFG, αC, and other structural features

### Integration Script Example

```python
import pandas as pd
import requests
import time

def update_conformation_data(csv_file='kinases_all.csv', output_file='kinases_with_conformation.csv'):
    """Update CSV with conformation data from KLIFS."""
    
    df = pd.read_csv(csv_file)
    
    conformations = []
    for idx, row in df.iterrows():
        if idx % 100 == 0:
            print(f"Processing {idx}/{len(df)}...")
        
        # Try to get from KLIFS
        conformation = get_klifs_conformation(row['uniprot_id'])
        conformations.append(conformation)
        
        # Be respectful of API rate limits
        time.sleep(0.1)
    
    df['conformation_DFG_aC'] = conformations
    df.to_csv(output_file, index=False)
    
    print(f"Updated {len(df)} kinases")
    print(f"With conformation data: {len(df[df['conformation_DFG_aC'] != 'Not available'])}")
```

## How to Populate Inhibitor Class Sensitivity

### Background

Kinase inhibitors are classified by their binding mode:

- **Type I**: Bind to active (DFG-in) conformation, compete with ATP
- **Type II**: Bind to inactive (DFG-out) conformation, access hydrophobic pocket
- **Type III**: Bind to allosteric sites, not ATP-competitive
- **Type IV**: Covalent inhibitors
- **Type V**: Bivalent inhibitors

### Data Sources

#### 1. ChEMBL Database (Recommended)
- **Website**: https://www.ebi.ac.uk/chembl/
- **API**: https://www.ebi.ac.uk/chembl/api/data/docs
- **Coverage**: >2 million bioactivity data points
- **Data**: Inhibitor potency, selectivity, binding modes

**Example API Usage:**
```python
import requests

def get_chembl_inhibitors(uniprot_id):
    """Get inhibitor data from ChEMBL."""
    url = f"https://www.ebi.ac.uk/chembl/api/data/activity.json?target_uniprot={uniprot_id}&pchembl_value__gte=6"
    response = requests.get(url)
    
    if response.status_code == 200:
        data = response.json()
        activities = data.get('activities', [])
        
        # Analyze inhibitor types from activities
        inhibitor_types = set()
        for activity in activities[:10]:  # Sample first 10
            # You would need additional logic to classify inhibitor types
            # This often requires cross-referencing with literature or databases
            pass
        
        if inhibitor_types:
            return ', '.join(inhibitor_types)
    
    return "Not available"
```

#### 2. KIDFamMap
- **Website**: http://mobiosd-hub.com/kidfammap/
- **Coverage**: ~300 kinases with detailed inhibitor classification
- **Data**: Type I, II, III classification based on binding modes

#### 3. DrugBank
- **Website**: https://www.drugbank.ca/
- **Coverage**: FDA-approved and experimental kinase inhibitors
- **Data**: Drug-target interactions, mechanism of action

#### 4. Literature-Based Curation
- Use PubMed API to search for inhibitor studies
- Extract inhibitor types from published papers
- Cross-reference with clinical trials data

### Integration Script Example

```python
import pandas as pd
import requests

def get_inhibitor_sensitivity(uniprot_id, gene_name):
    """Get inhibitor sensitivity from multiple sources."""
    
    # Try ChEMBL
    chembl_url = f"https://www.ebi.ac.uk/chembl/api/data/activity.json?target_uniprot={uniprot_id}&pchembl_value__gte=7"
    response = requests.get(chembl_url)
    
    if response.status_code == 200:
        data = response.json()
        if data.get('activities'):
            # Simple classification based on number of inhibitors
            count = len(data['activities'])
            if count > 100:
                return "High sensitivity (Type I/II)"
            elif count > 20:
                return "Moderate sensitivity"
            elif count > 0:
                return "Low sensitivity"
    
    return "Not available"

def update_inhibitor_data(csv_file='kinases_all.csv', output_file='kinases_with_inhibitors.csv'):
    """Update CSV with inhibitor sensitivity data."""
    
    df = pd.read_csv(csv_file)
    
    sensitivities = []
    for idx, row in df.iterrows():
        if idx % 50 == 0:
            print(f"Processing {idx}/{len(df)}...")
        
        # Extract gene name from protein_name
        gene_name = row['protein_name'].split()[0]
        sensitivity = get_inhibitor_sensitivity(row['uniprot_id'], gene_name)
        sensitivities.append(sensitivity)
        
        time.sleep(0.2)  # Rate limiting
    
    df['inhibitor_class_sensitivity'] = sensitivities
    df.to_csv(output_file, index=False)
    
    print(f"Updated {len(df)} kinases")
```

## Recommended Workflow

1. **Start with KLIFS** for conformation data (best coverage for kinases)
2. **Use ChEMBL** for inhibitor activity data
3. **Cross-reference with literature** for specific inhibitor types
4. **Validate** against experimental data when available

## Important Considerations

### Data Quality
- Not all kinases have structural data (PDB coverage ~30-40%)
- Inhibitor classification requires careful curation
- Multiple conformations may exist for same kinase

### Rate Limiting
- Respect API rate limits (typically 10-20 requests/second)
- Cache results to avoid repeated queries
- Consider bulk downloads when available

### Data Updates
- Databases are regularly updated with new structures and inhibitors
- Plan for periodic re-downloads to keep data current
- Version your datasets

## Example Complete Integration Script

```python
import pandas as pd
import requests
import time
from typing import Optional

def integrate_external_data(
    input_csv='kinases_all.csv',
    output_csv='kinases_complete.csv',
    use_klifs=True,
    use_chembl=True
):
    """
    Integrate external data sources to populate conformation and inhibitor columns.
    """
    
    df = pd.read_csv(input_csv)
    
    print(f"Processing {len(df)} kinases...")
    print(f"Fetching from: {'KLIFS ' if use_klifs else ''}{'ChEMBL' if use_chembl else ''}")
    
    for idx, row in df.iterrows():
        if idx % 100 == 0:
            print(f"Progress: {idx}/{len(df)} ({idx/len(df)*100:.1f}%)")
        
        uniprot_id = row['uniprot_id']
        
        # Update conformation from KLIFS
        if use_klifs and row['conformation_DFG_aC'] == 'Not available':
            conformation = get_klifs_conformation(uniprot_id)
            df.at[idx, 'conformation_DFG_aC'] = conformation
            time.sleep(0.1)
        
        # Update inhibitor sensitivity from ChEMBL
        if use_chembl and row['inhibitor_class_sensitivity'] == 'Not available':
            sensitivity = get_inhibitor_sensitivity(uniprot_id, "")
            df.at[idx, 'inhibitor_class_sensitivity'] = sensitivity
            time.sleep(0.2)
    
    # Save updated dataset
    df.to_csv(output_csv, index=False)
    
    # Print statistics
    print(f"\nIntegration Complete!")
    print(f"Conformations populated: {len(df[df['conformation_DFG_aC'] != 'Not available'])} ({len(df[df['conformation_DFG_aC'] != 'Not available'])/len(df)*100:.1f}%)")
    print(f"Inhibitor data populated: {len(df[df['inhibitor_class_sensitivity'] != 'Not available'])} ({len(df[df['inhibitor_class_sensitivity'] != 'Not available'])/len(df)*100:.1f}%)")
    
    return df

# Run integration
# integrated_df = integrate_external_data()
```

## Further Resources

- **KLIFS Documentation**: https://klifs.net/swagger/
- **ChEMBL Web Services**: https://chembl.gitbook.io/chembl-interface-documentation/web-services
- **Manning Kinome Classification**: http://kinase.com/web/current/kinbase/
- **PDB API**: https://www.rcsb.org/docs/programmatic-access/web-services-overview
- **Review Article on Kinase Inhibitor Types**: Nature Reviews Drug Discovery (search for "kinase inhibitor classification")

## Support

For questions or issues with data integration:
1. Check database documentation and API limits
2. Validate UniProt IDs are current
3. Consider organism-specific data availability
4. Contact database support teams for bulk access options
