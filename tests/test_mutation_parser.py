#!/usr/bin/env python3
"""
Unit tests for mutation parser.

Tests various mutation formats and edge cases.
"""

import pytest
from mutation_motif_analysis import MutationParser


def test_parse_hgvs_format():
    """Test HGVS format (p.R90H)."""
    result = MutationParser.parse("p.R90H")
    
    assert result is not None
    assert result['ref_aa'] == 'R'
    assert result['position_1based'] == 90
    assert result['position_0based'] == 89
    assert result['alt_aa'] == 'H'


def test_parse_simple_format():
    """Test simple format (R90H)."""
    result = MutationParser.parse("R90H")
    
    assert result is not None
    assert result['ref_aa'] == 'R'
    assert result['position_1based'] == 90
    assert result['alt_aa'] == 'H'


def test_parse_lowercase():
    """Test lowercase input."""
    result = MutationParser.parse("r90h")
    
    assert result is not None
    assert result['ref_aa'] == 'R'
    assert result['alt_aa'] == 'H'


def test_parse_nonsense_mutation():
    """Test nonsense mutation (to stop codon)."""
    result = MutationParser.parse("R90*")
    
    assert result is not None
    assert result['ref_aa'] == 'R'
    assert result['alt_aa'] == '*'


def test_parse_invalid_format():
    """Test invalid mutation string."""
    result = MutationParser.parse("invalid")
    assert result is None
    
    result = MutationParser.parse("R-90-H")
    assert result is None
    
    result = MutationParser.parse("90RH")  # Missing ref
    assert result is None


def test_validate_against_sequence():
    """Test mutation validation against sequence."""
    mutation = {
        'ref_aa': 'K',
        'position_0based': 5,
        'position_1based': 6,
        'alt_aa': 'R'
    }
    
    # Valid mutation
    sequence = "ABCDEKFGH"
    valid, msg = MutationParser.validate_against_sequence(mutation, sequence)
    assert valid, f"Should be valid: {msg}"
    
    # Invalid: wrong amino acid
    sequence = "ABCDERGH"  # Position 5 is R, not K
    valid, msg = MutationParser.validate_against_sequence(mutation, sequence)
    assert not valid, "Should detect mismatch"
    assert "Mismatch" in msg


def test_validate_out_of_range():
    """Test mutation outside sequence range."""
    mutation = {
        'ref_aa': 'K',
        'position_0based': 100,
        'position_1based': 101,
        'alt_aa': 'R'
    }
    
    sequence = "ABCDEFGH"  # Only 8 residues
    valid, msg = MutationParser.validate_against_sequence(mutation, sequence)
    assert not valid, "Should detect out of range"
    assert "out of range" in msg.lower()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

