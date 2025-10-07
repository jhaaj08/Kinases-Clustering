#!/usr/bin/env python3
"""
Initialize provenance tracking for the kinase clustering project.

Captures:
- Tool versions (HMMER, CD-HIT, Python packages)
- Data sources (UniProt, Pfam)
- Inclusion/exclusion rules
- Processing parameters

Run this once after setup to document the computational environment.
"""

import os
import sys
import subprocess
import datetime
from pathlib import Path
from utils.provenance import ProvenanceTracker, get_tool_version


def main():
    print("="*80)
    print("INITIALIZING PROVENANCE TRACKING")
    print("="*80)
    print()
    
    # Create tracker
    prov = ProvenanceTracker(output_dir="data")
    
    # Capture environment
    print("Capturing environment...")
    prov.add_environment_info()
    print("✅ Environment captured")
    
    # Capture tool versions
    print("\nCapturing tool versions...")
    
    # HMMER
    hmmer_version = get_tool_version('hmmsearch')
    if hmmer_version:
        print(f"  HMMER: {hmmer_version}")
        prov.add_hmmer_info(
            version=hmmer_version,
            params={
                "evalue_threshold": 0.001,
                "boundary_type": "envelope",
                "domain_selection": "best (lowest E-value, then highest score)",
                "multi_domain_handling": "keep best domain per sequence",
                "coordinate_system": "1-based HMMER → 0-based Python via slice [start-1:end]",
                "fallback_strategy": "PF00069 → PF07714 → drop if none found",
            }
        )
    else:
        print("  ⚠️  HMMER not found")
    
    # CD-HIT
    cdhit_version = get_tool_version('cd-hit')
    if cdhit_version:
        print(f"  CD-HIT: {cdhit_version}")
        prov.add_cdhit_info(
            version=cdhit_version,
            thresholds={
                "data_cleaning": 0.60,
                "homology_aware_splits": 0.40,
                "rationale": "60% for redundancy reduction, 40% for preventing test leakage"
            }
        )
    else:
        print("  ⚠️  CD-HIT not found")
    
    print("✅ Tool versions captured")
    
    # Document inclusion/exclusion rules
    print("\nDocumenting inclusion/exclusion rules...")
    prov.add_inclusion_exclusion_rules({
        "data_source": "UniProt SwissProt (reviewed entries only)",
        "query": "reviewed:true AND (keyword:KW-0418 OR name:kinase*)",
        "isoforms": "canonical isoform only (default from UniProt)",
        "fragments": "removed if flagged as fragment by UniProt",
        "minimum_sequence_length": 100,
        "minimum_domain_length": 50,
        "domain_required": "Yes - sequences without PF00069/PF07714 excluded from domain analysis",
        "multi_domain_handling": "keep best-scoring domain per sequence (lowest E-value, then highest bit score)",
        "controlled_vocabulary_labels": [
            "AGC", "CAMK", "CK1", "CMGC", "STE", "TK", "TKL", 
            "RGC", "Atypical", "Histidine", "Other"
        ],
        "label_source": "UniProt kinome_group annotations + Manning classification",
        "missing_labels": "classified as 'Other'",
        "minimum_class_size_for_training": 5,
    })
    print("✅ Rules documented")
    
    # Document Pfam sources
    print("\nDocumenting Pfam sources...")
    prov.add_pfam_info(
        pfam_ids=['PF00069', 'PF07714'],
        urls={
            'PF00069': 'https://www.ebi.ac.uk/interpro/api/entry/pfam/PF00069?annotation=hmm',
            'PF07714': 'https://www.ebi.ac.uk/interpro/api/entry/pfam/PF07714?annotation=hmm'
        },
        response_metadata={
            'note': 'HMM files downloaded via InterPro API (gzip-compressed)',
            'access_date': datetime.datetime.now().isoformat(),
        }
    )
    print("✅ Pfam sources documented")
    
    # Save summary
    print()
    print(prov.get_summary())
    
    print()
    print("="*80)
    print("✅ PROVENANCE INITIALIZATION COMPLETE!")
    print("="*80)
    print(f"\nProvenance file: {prov.provenance_file}")
    print()
    print("This file contains:")
    print("  • Tool versions (HMMER, CD-HIT, Python packages)")
    print("  • Data sources (UniProt, Pfam)")
    print("  • Inclusion/exclusion rules")
    print("  • Processing parameters")
    print("  • Coordinate systems and conventions")
    print()
    print("Use this file to:")
    print("  • Document methods in publications")
    print("  • Reproduce analysis with exact versions")
    print("  • Verify data integrity")
    print()


if __name__ == "__main__":
    main()
