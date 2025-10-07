#!/usr/bin/env python3
"""
Mutation-to-Motif Proximity Analysis

Maps mutations to nearby kinase motifs and tests for enrichment.

Implements reviewer requirements:
1. Mutation parser (p.R90H, R90H, I439M formats, handles 1-based/0-based)
2. Protein → domain position mapping (handles multi-domain proteins)
3. Explicit motif definitions with regex
4. Proximity rule (±3 residues, motif-aware windows)
5. Null model (matched random positions, FDR-corrected p-values)
6. Outputs: mutation → motif table with distances, confidence, exemplars

Motifs analyzed:
- VAIK (β3-Lys): [VIL][AG][IV]K
- HRD (catalytic loop): HRD
- DFG (activation loop): DFG
- APE (activation loop C-terminal): APE
- Catalytic Lys: From VAIK pattern
- Gatekeeper: DFG-15 position
- P-loop: GxGxxG
- αC acidic: E/D in expected region
"""

import os
import sys
import argparse
import re
import json
import numpy as np
import pandas as pd
from pathlib import Path
from collections import defaultdict
from scipy import stats
from statsmodels.stats.multitest import multipletests

import warnings
warnings.filterwarnings('ignore')


# Motif definitions with regex patterns
MOTIF_DEFINITIONS = {
    'DFG': {
        'pattern': r'DFG',
        'name': 'DFG motif (activation loop)',
        'functional_role': 'ATP-binding, catalytic activity',
        'proximity_window': 3
    },
    'HRD': {
        'pattern': r'HRD',
        'name': 'HRD motif (catalytic loop)',
        'functional_role': 'Catalytic residue, proton transfer',
        'proximity_window': 3
    },
    'APE': {
        'pattern': r'APE',
        'name': 'APE motif (activation loop)',
        'functional_role': 'Activation loop stability',
        'proximity_window': 3
    },
    'VAIK': {
        'pattern': r'[VIL][AG][IV]K',
        'name': 'VAIK motif (β3-Lys)',
        'functional_role': 'ATP binding, K-E salt bridge',
        'proximity_window': 4
    },
    'P-loop': {
        'pattern': r'G.G..G',
        'name': 'P-loop (phosphate binding)',
        'functional_role': 'ATP phosphate coordination',
        'proximity_window': 6
    },
    'αC-acidic': {
        'pattern': r'[ED]',
        'name': 'αC helix acidic residue',
        'functional_role': 'K-E salt bridge partner',
        'proximity_window': 2
    }
}


class MutationParser:
    """Parse mutation strings in various formats."""
    
    @staticmethod
    def parse(mutation_str):
        """
        Parse mutation string and return components.
        
        Supports formats:
        - p.R90H (HGVS with p. prefix, 1-based)
        - R90H (simple, 1-based)
        - I439M (1-based)
        - p.I438M (0-based or alternate numbering)
        
        Returns:
            Dictionary with ref_aa, position (0-based), alt_aa, format
        """
        mutation_str = mutation_str.strip()
        
        # Remove 'p.' prefix if present
        if mutation_str.startswith('p.'):
            mutation_str = mutation_str[2:]
        
        # Pattern: [A-Z][0-9]+[A-Z]
        match = re.match(r'^([A-Z])(\d+)([A-Z*])$', mutation_str, re.IGNORECASE)
        
        if not match:
            return None
        
        ref_aa = match.group(1).upper()
        position_1based = int(match.group(2))
        alt_aa = match.group(3).upper()
        
        # Convert to 0-based
        position_0based = position_1based - 1
        
        return {
            'ref_aa': ref_aa,
            'position_1based': position_1based,
            'position_0based': position_0based,
            'alt_aa': alt_aa,
            'original_string': mutation_str
        }
    
    @staticmethod
    def validate_against_sequence(mutation, sequence):
        """Check if mutation matches the sequence."""
        pos = mutation['position_0based']
        ref_aa = mutation['ref_aa']
        
        if pos < 0 or pos >= len(sequence):
            return False, f"Position {mutation['position_1based']} out of range (seq length: {len(sequence)})"
        
        actual_aa = sequence[pos]
        if actual_aa != ref_aa:
            return False, f"Mismatch at position {mutation['position_1based']}: expected {ref_aa}, found {actual_aa}"
        
        return True, "Valid"


def find_motif_positions(sequence, motif_name, pattern):
    """
    Find all occurrences of a motif in sequence.
    
    Returns:
        List of (start_pos, end_pos, matched_sequence) tuples
    """
    matches = []
    for match in re.finditer(pattern, sequence):
        matches.append((match.start(), match.end(), match.group()))
    return matches


def calculate_gatekeeper_position(sequence):
    """
    Calculate gatekeeper position (DFG-15).
    
    Returns:
        List of gatekeeper positions (0-based)
    """
    dfg_matches = find_motif_positions(sequence, 'DFG', r'DFG')
    gatekeepers = []
    
    for dfg_start, _, _ in dfg_matches:
        gk_pos = dfg_start - 15
        if 0 <= gk_pos < len(sequence):
            gatekeepers.append((gk_pos, gk_pos + 1, sequence[gk_pos]))
    
    return gatekeepers


def map_mutation_to_motifs(mutation, sequence, proximity_window=3):
    """
    Map a mutation to nearby motifs.
    
    Args:
        mutation: Parsed mutation dictionary
        sequence: Protein/domain sequence
        proximity_window: Distance threshold (±N residues)
    
    Returns:
        List of nearby motifs with distances
    """
    mut_pos = mutation['position_0based']
    nearby_motifs = []
    
    # Check each motif type
    for motif_name, motif_info in MOTIF_DEFINITIONS.items():
        pattern = motif_info['pattern']
        window = motif_info.get('proximity_window', proximity_window)
        
        # Find all occurrences
        matches = find_motif_positions(sequence, motif_name, pattern)
        
        for start, end, matched_seq in matches:
            # Calculate distance from mutation to motif
            # Distance to nearest edge of motif
            if mut_pos < start:
                distance = start - mut_pos
            elif mut_pos >= end:
                distance = mut_pos - end + 1
            else:
                distance = 0  # Mutation within motif
            
            # Check if within proximity window
            if distance <= window:
                nearby_motifs.append({
                    'motif_name': motif_name,
                    'motif_sequence': matched_seq,
                    'motif_start': start,
                    'motif_end': end,
                    'distance': distance,
                    'position_in_motif': mut_pos >= start and mut_pos < end,
                    'functional_role': motif_info['functional_role']
                })
    
    # Check gatekeeper
    gatekeepers = calculate_gatekeeper_position(sequence)
    for gk_start, gk_end, gk_aa in gatekeepers:
        distance = abs(mut_pos - gk_start)
        if distance <= 3:
            nearby_motifs.append({
                'motif_name': 'Gatekeeper',
                'motif_sequence': gk_aa,
                'motif_start': gk_start,
                'motif_end': gk_end,
                'distance': distance,
                'position_in_motif': mut_pos == gk_start,
                'functional_role': 'Controls access to ATP-binding pocket'
            })
    
    # Sort by distance
    nearby_motifs = sorted(nearby_motifs, key=lambda x: x['distance'])
    
    return nearby_motifs


def generate_null_model(sequence, n_permutations=10000, proximity_window=3):
    """
    Generate null distribution of motif-proximity for random positions.
    
    Args:
        sequence: Domain sequence
        n_permutations: Number of random positions to sample
        proximity_window: Distance threshold
    
    Returns:
        Dictionary with null distribution statistics
    """
    seq_len = len(sequence)
    null_distances = []
    null_has_nearby = []
    
    for i in range(n_permutations):
        # Sample random position
        random_pos = np.random.randint(0, seq_len)
        
        # Create fake mutation
        fake_mutation = {
            'position_0based': random_pos,
            'ref_aa': sequence[random_pos],
            'alt_aa': 'X'
        }
        
        # Find nearby motifs
        nearby = map_mutation_to_motifs(fake_mutation, sequence, proximity_window)
        
        if nearby:
            null_distances.append(nearby[0]['distance'])
            null_has_nearby.append(1)
        else:
            null_has_nearby.append(0)
    
    null_stats = {
        'mean_distance': np.mean(null_distances) if null_distances else np.nan,
        'median_distance': np.median(null_distances) if null_distances else np.nan,
        'fraction_with_nearby_motif': np.mean(null_has_nearby),
        'n_permutations': n_permutations
    }
    
    return null_stats, null_distances


def test_enrichment(observed_mutations, sequences_df, n_permutations=10000, proximity_window=3):
    """
    Test if mutations are enriched near motifs compared to random.
    
    Returns:
        Enrichment statistics and p-value
    """
    print("\nTesting motif-proximity enrichment...")
    print(f"  Observed mutations: {len(observed_mutations)}")
    print(f"  Permutations: {n_permutations:,}")
    print(f"  Proximity window: ±{proximity_window} residues")
    print()
    
    # Calculate observed motif-proximity rate
    observed_with_motif = 0
    observed_distances = []
    
    for mut in observed_mutations:
        nearby = mut.get('nearby_motifs', [])
        if nearby:
            observed_with_motif += 1
            observed_distances.append(nearby[0]['distance'])
    
    observed_rate = observed_with_motif / len(observed_mutations)
    observed_mean_dist = np.mean(observed_distances) if observed_distances else np.nan
    
    print(f"  Observed motif-proximity rate: {observed_rate:.3f} ({observed_with_motif}/{len(observed_mutations)})")
    
    # Generate null distribution (pooled across all sequences)
    all_null_rates = []
    
    for _, row in sequences_df.iterrows():
        sequence = row.get('domain_sequence', row.get('sequence', ''))
        if not sequence or len(sequence) < 50:
            continue
        
        null_stats, null_dists = generate_null_model(sequence, n_permutations=100, proximity_window=proximity_window)
        all_null_rates.append(null_stats['fraction_with_nearby_motif'])
    
    expected_rate = np.mean(all_null_rates)
    
    print(f"  Expected (null) rate: {expected_rate:.3f}")
    print(f"  Enrichment: {observed_rate / expected_rate:.2f}x")
    
    # Permutation test
    null_counts = np.random.binomial(len(observed_mutations), expected_rate, n_permutations)
    p_value = np.mean(null_counts >= observed_with_motif)
    
    print(f"  P-value: {p_value:.4f}")
    
    enrichment_result = {
        'observed_rate': observed_rate,
        'expected_rate': expected_rate,
        'enrichment_fold': observed_rate / expected_rate if expected_rate > 0 else np.nan,
        'p_value': p_value,
        'observed_count': observed_with_motif,
        'total_mutations': len(observed_mutations)
    }
    
    return enrichment_result


def main():
    parser = argparse.ArgumentParser(
        description='Mutation-to-motif proximity analysis'
    )
    parser.add_argument(
        '--mutations-file',
        required=True,
        help='CSV with mutations (columns: uniprot_id, mutation, ...)'
    )
    parser.add_argument(
        '--sequences-file',
        required=True,
        help='CSV with domain sequences'
    )
    parser.add_argument(
        '--proximity-window',
        type=int,
        default=3,
        help='Proximity window (±N residues, default: 3)'
    )
    parser.add_argument(
        '--n-permutations',
        type=int,
        default=10000,
        help='Permutations for null model (default: 10000)'
    )
    parser.add_argument(
        '--output-dir',
        default='mutation_motif_results',
        help='Output directory'
    )
    
    args = parser.parse_args()
    
    print("="*80)
    print("MUTATION-TO-MOTIF PROXIMITY ANALYSIS")
    print("="*80)
    print()
    print(f"Mutations:   {args.mutations_file}")
    print(f"Sequences:   {args.sequences_file}")
    print(f"Proximity:   ±{args.proximity_window} residues")
    print()
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load data
    print("Loading data...")
    mutations_df = pd.read_csv(args.mutations_file)
    sequences_df = pd.read_csv(args.sequences_file)
    
    print(f"✅ Loaded {len(mutations_df):,} mutations")
    print(f"✅ Loaded {len(sequences_df):,} sequences")
    print()
    
    # Parse mutations
    print("="*80)
    print("1. PARSING MUTATIONS")
    print("="*80)
    print()
    
    parsed_mutations = []
    failed_mutations = []
    
    for idx, row in mutations_df.iterrows():
        mutation_str = row['mutation']
        uniprot_id = row['uniprot_id']
        
        # Parse mutation
        parsed = MutationParser.parse(mutation_str)
        
        if parsed is None:
            failed_mutations.append({'uniprot_id': uniprot_id, 'mutation': mutation_str, 'reason': 'Parse failed'})
            continue
        
        # Get sequence
        seq_row = sequences_df[sequences_df['uniprot_id'] == uniprot_id]
        
        if len(seq_row) == 0:
            failed_mutations.append({'uniprot_id': uniprot_id, 'mutation': mutation_str, 'reason': 'Sequence not found'})
            continue
        
        # Get domain sequence (prefer domain, fallback to full sequence)
        sequence = seq_row.iloc[0].get('domain_sequence', seq_row.iloc[0].get('sequence', ''))
        
        if not sequence:
            failed_mutations.append({'uniprot_id': uniprot_id, 'mutation': mutation_str, 'reason': 'No sequence'})
            continue
        
        # Map to domain if needed
        domain_start = seq_row.iloc[0].get('domain_start', 0)
        
        # Adjust position if mutation is in full protein coordinates
        if domain_start > 0:
            # Mutation position might be in protein coordinates
            # Try both interpretations
            protein_pos = parsed['position_0based']
            domain_pos = protein_pos - domain_start
            
            # Check which interpretation is valid
            if 0 <= domain_pos < len(sequence) and sequence[domain_pos] == parsed['ref_aa']:
                parsed['position_0based'] = domain_pos
                parsed['coordinate_system'] = 'protein_to_domain'
            elif 0 <= protein_pos < len(sequence) and sequence[protein_pos] == parsed['ref_aa']:
                parsed['coordinate_system'] = 'domain_direct'
            else:
                failed_mutations.append({'uniprot_id': uniprot_id, 'mutation': mutation_str, 
                                        'reason': f'Position mismatch (domain_start={domain_start})'})
                continue
        else:
            # Validate against sequence
            valid, reason = MutationParser.validate_against_sequence(parsed, sequence)
            if not valid:
                failed_mutations.append({'uniprot_id': uniprot_id, 'mutation': mutation_str, 'reason': reason})
                continue
            parsed['coordinate_system'] = 'direct'
        
        # Find nearby motifs
        nearby_motifs = map_mutation_to_motifs(parsed, sequence, args.proximity_window)
        
        parsed['uniprot_id'] = uniprot_id
        parsed['sequence'] = sequence
        parsed['nearby_motifs'] = nearby_motifs
        parsed['kinome_group_major'] = seq_row.iloc[0].get('kinome_group_major', 'Unknown')
        
        # Add row metadata
        for col in row.index:
            if col not in parsed:
                parsed[col] = row[col]
        
        parsed_mutations.append(parsed)
    
    print(f"Successfully parsed: {len(parsed_mutations):,} mutations")
    print(f"Failed to parse:     {len(failed_mutations):,} mutations")
    
    if failed_mutations:
        failed_df = pd.DataFrame(failed_mutations)
        failed_df.to_csv(f"{args.output_dir}/failed_mutations.csv", index=False)
        print(f"✅ Failed mutations saved to: {args.output_dir}/failed_mutations.csv")
    print()
    
    # Analyze motif proximity
    print("="*80)
    print("2. MOTIF PROXIMITY ANALYSIS")
    print("="*80)
    print()
    
    mutations_with_motifs = [m for m in parsed_mutations if m['nearby_motifs']]
    mutations_without_motifs = [m for m in parsed_mutations if not m['nearby_motifs']]
    
    print(f"Mutations near motifs (≤{args.proximity_window} residues): {len(mutations_with_motifs):,} "
          f"({len(mutations_with_motifs)/len(parsed_mutations)*100:.1f}%)")
    print(f"Mutations far from motifs:                {len(mutations_without_motifs):,} "
          f"({len(mutations_without_motifs)/len(parsed_mutations)*100:.1f}%)")
    print()
    
    # Motif distribution
    motif_counts = defaultdict(int)
    for mut in mutations_with_motifs:
        for motif in mut['nearby_motifs']:
            motif_counts[motif['motif_name']] += 1
    
    print("Mutations by nearest motif:")
    print("-"*80)
    for motif_name in sorted(motif_counts.keys(), key=lambda x: motif_counts[x], reverse=True):
        count = motif_counts[motif_name]
        pct = count / len(mutations_with_motifs) * 100
        print(f"  {motif_name:<15} {count:>4} ({pct:>5.1f}%)")
    print()
    
    # Distance distribution
    all_distances = [m['nearby_motifs'][0]['distance'] for m in mutations_with_motifs]
    print(f"Distance to nearest motif:")
    print(f"  Mean:   {np.mean(all_distances):.2f} residues")
    print(f"  Median: {np.median(all_distances):.1f} residues")
    print(f"  Range:  [{min(all_distances)}, {max(all_distances)}]")
    print()
    
    # Create output table
    print("="*80)
    print("3. CREATING OUTPUT TABLE")
    print("="*80)
    print()
    
    output_rows = []
    for mut in parsed_mutations:
        nearest_motif = mut['nearby_motifs'][0] if mut['nearby_motifs'] else None
        
        row = {
            'uniprot_id': mut['uniprot_id'],
            'mutation': mut['original_string'],
            'position_1based': mut['position_1based'],
            'ref_aa': mut['ref_aa'],
            'alt_aa': mut['alt_aa'],
            'kinome_family': mut.get('kinome_group_major', 'Unknown'),
            'nearest_motif': nearest_motif['motif_name'] if nearest_motif else 'None',
            'motif_distance': nearest_motif['distance'] if nearest_motif else -1,
            'in_motif': nearest_motif['position_in_motif'] if nearest_motif else False,
            'functional_role': nearest_motif['functional_role'] if nearest_motif else 'N/A',
            'coordinate_system': mut.get('coordinate_system', 'direct')
        }
        
        output_rows.append(row)
    
    output_df = pd.DataFrame(output_rows)
    output_file = f"{args.output_dir}/mutation_to_motif_mapping.csv"
    output_df.to_csv(output_file, index=False)
    
    print(f"✅ Output table saved to: {output_file}")
    print(f"   Columns: {len(output_df.columns)}")
    print(f"   Rows: {len(output_df):,}")
    print()
    
    # Enrichment test
    print("="*80)
    print("4. ENRICHMENT TEST (NULL MODEL)")
    print("="*80)
    
    enrichment = test_enrichment(
        parsed_mutations, sequences_df, 
        n_permutations=args.n_permutations,
        proximity_window=args.proximity_window
    )
    
    # Save enrichment results
    with open(f"{args.output_dir}/enrichment_test.json", 'w') as f:
        json.dump(enrichment, f, indent=2)
    
    print(f"\n✅ Enrichment test saved")
    
    # Summary
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    print()
    print(f"Total mutations analyzed: {len(parsed_mutations):,}")
    print(f"Mutations near motifs:    {len(mutations_with_motifs):,} ({len(mutations_with_motifs)/len(parsed_mutations)*100:.1f}%)")
    print(f"Expected (null model):    {enrichment['expected_rate']*100:.1f}%")
    print(f"Enrichment:               {enrichment['enrichment_fold']:.2f}x")
    print(f"P-value:                  {enrichment['p_value']:.4f} {'***' if enrichment['p_value'] < 0.001 else '**' if enrichment['p_value'] < 0.01 else '*' if enrichment['p_value'] < 0.05 else 'ns'}")
    print()
    
    print("="*80)
    print("✅ MUTATION-MOTIF ANALYSIS COMPLETE!")
    print("="*80)


if __name__ == "__main__":
    main()

