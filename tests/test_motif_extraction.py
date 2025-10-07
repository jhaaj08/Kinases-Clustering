#!/usr/bin/env python3
"""
Unit tests for motif extraction.

Tests motif finding with known sequences.
"""

import pytest
import sys
sys.path.insert(0, '.')
from extract_motif_features import find_motif, find_ploop, find_vaik_region


def test_find_dfg_motif():
    """Test DFG motif finding."""
    sequence = "ABCDEFGIJK"  # DFG at position 4-6
    found, pos = find_motif(sequence, r'DFG', 'DFG')
    
    assert found, "DFG should be found"
    assert pos == 4, f"Expected position 4, got {pos}"


def test_find_hrd_motif():
    """Test HRD motif finding."""
    sequence = "XXXXXHRDXXX"  # HRD at position 5-7
    found, pos = find_motif(sequence, r'HRD', 'HRD')
    
    assert found, "HRD should be found"
    assert pos == 5, f"Expected position 5, got {pos}"


def test_motif_not_found():
    """Test motif not found."""
    sequence = "ABCDEFGHIJK"  # No HRD
    found, pos = find_motif(sequence, r'HRD', 'HRD')
    
    assert not found, "HRD should not be found"
    assert pos == -1, "Position should be -1 when not found"


def test_find_ploop():
    """Test P-loop (GxGxxG) finding."""
    sequence = "XXXXXGAGAAGXXX"  # GAGAAG matches GxGxxG at position 5
    found, pos, score = find_ploop(sequence)
    
    assert found, "P-loop should be found"
    assert pos == 5, f"Expected position 5, got {pos}"
    assert score == 1.0, f"Perfect match should have score 1.0, got {score}"


def test_find_ploop_imperfect():
    """Test P-loop with imperfect match."""
    sequence = "XXXXXGAGAAXXX"  # GAGAA (5th position not G, score 2/3)
    found, pos, score = find_ploop(sequence)
    
    assert found, "P-loop should be found"
    assert abs(score - 2.0/3.0) < 1e-10, f"Expected score 0.667, got {score}"


def test_find_vaik():
    """Test VAIK motif finding."""
    sequence = "XXXXXVAIKXXX"  # Exact VAIK at position 5
    found, pos, k_pos = find_vaik_region(sequence)
    
    assert found, "VAIK should be found"
    assert pos == 5, f"Expected position 5, got {pos}"
    assert k_pos == 8, f"Expected K position 8, got {k_pos}"


def test_find_vaik_variant():
    """Test VAIK variant (LAIK)."""
    sequence = "XXXXXLAIKXXX"  # LAIK at position 5
    found, pos, k_pos = find_vaik_region(sequence)
    
    assert found, "LAIK variant should be found"
    assert pos == 5, f"Expected position 5, got {pos}"


def test_multiple_motifs():
    """Test sequence with multiple DFG motifs."""
    sequence = "DFGXXXDFGXXXDFG"  # Three DFG motifs
    found, pos = find_motif(sequence, r'DFG', 'DFG')
    
    assert found, "First DFG should be found"
    assert pos == 0, f"Should find first occurrence at position 0, got {pos}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

