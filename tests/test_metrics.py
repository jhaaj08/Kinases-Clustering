#!/usr/bin/env python3
"""
Unit tests for clustering and classification metrics.

Tests metrics on toy data with known ground truth answers.
"""

import numpy as np
import pytest
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score


def test_ari_perfect_clustering():
    """Test ARI with perfect clustering."""
    labels_true = np.array([0, 0, 0, 1, 1, 1])
    labels_pred = np.array([0, 0, 0, 1, 1, 1])
    
    ari = adjusted_rand_score(labels_true, labels_pred)
    assert abs(ari - 1.0) < 1e-10, f"Expected ARI=1.0, got {ari}"


def test_ari_random_clustering():
    """Test ARI with random clustering (should be near 0)."""
    labels_true = np.array([0, 0, 0, 1, 1, 1])
    labels_pred = np.array([0, 1, 0, 1, 0, 1])
    
    ari = adjusted_rand_score(labels_true, labels_pred)
    # Random should be near 0 (within reasonable bounds)
    assert -0.2 < ari < 0.2, f"Expected ARI≈0, got {ari}"


def test_nmi_perfect_clustering():
    """Test NMI with perfect clustering."""
    labels_true = np.array([0, 0, 0, 1, 1, 1])
    labels_pred = np.array([0, 0, 0, 1, 1, 1])
    
    nmi = normalized_mutual_info_score(labels_true, labels_pred)
    assert abs(nmi - 1.0) < 1e-10, f"Expected NMI=1.0, got {nmi}"


def test_purity_perfect():
    """Test purity calculation with perfect clustering."""
    from clustering_statistics import calculate_purity
    
    labels_true = np.array([0, 0, 0, 1, 1, 1])
    labels_pred = np.array([0, 0, 0, 1, 1, 1])
    
    purity = calculate_purity(labels_true, labels_pred)
    assert abs(purity - 1.0) < 1e-10, f"Expected purity=1.0, got {purity}"


def test_purity_known_value():
    """Test purity with known ground truth."""
    from clustering_statistics import calculate_purity
    
    # Cluster 0: 3 of class A, 0 of class B
    # Cluster 1: 1 of class A, 2 of class B
    # Purity = (3 + 2) / 6 = 5/6 ≈ 0.833
    labels_true = np.array([0, 0, 0, 1, 1, 1])
    labels_pred = np.array([0, 0, 0, 1, 0, 1])
    
    purity = calculate_purity(labels_true, labels_pred)
    expected = 5.0 / 6.0
    assert abs(purity - expected) < 1e-10, f"Expected purity={expected:.3f}, got {purity:.3f}"


def test_hungarian_accuracy():
    """Test Hungarian accuracy calculation."""
    from clustering_statistics import calculate_hungarian_accuracy
    
    # Perfect clustering (just relabeled)
    labels_true = np.array([0, 0, 0, 1, 1, 1])
    labels_pred = np.array([1, 1, 1, 0, 0, 0])  # Swapped labels
    
    hungarian_acc = calculate_hungarian_accuracy(labels_true, labels_pred)
    assert abs(hungarian_acc - 1.0) < 1e-10, f"Expected hungarian_acc=1.0, got {hungarian_acc}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

