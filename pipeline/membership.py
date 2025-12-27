#!/usr/bin/env python3
"""
Membership Module - Dataset membership validation and assertions.

This module enforces that manifest files are the SINGLE SOURCE OF TRUTH
for dataset membership. All downstream steps must derive from manifests.

Usage:
    from pipeline.membership import load_manifest, assert_split_integrity
    
    manifest_ids = load_manifest("supervised_eligible", run_dir)
    assert_split_integrity(run_dir, "split40")
"""

from pathlib import Path
from typing import Set, Optional, Dict, List


def load_manifest(name: str, run_dir: Path) -> Set[str]:
    """
    Load manifest IDs - the ONLY way to get dataset membership.
    
    Args:
        name: Manifest name (e.g., "supervised_eligible", "domain_E001")
        run_dir: Path to the run directory
    
    Returns:
        Set of UniProt IDs in the manifest.
    
    Raises:
        FileNotFoundError: If manifest file doesn't exist.
        ValueError: If manifest is empty.
    """
    manifest_file = run_dir / "data" / "manifests" / f"{name}.txt"
    
    if not manifest_file.exists():
        raise FileNotFoundError(f"Manifest not found: {manifest_file}")
    
    with open(manifest_file) as f:
        ids = {line.strip() for line in f if line.strip()}
    
    if not ids:
        raise ValueError(f"Manifest is empty: {manifest_file}")
    
    return ids


def load_split(split_name: str, split_type: str, run_dir: Path) -> Set[str]:
    """
    Load split IDs (train or test).
    
    Args:
        split_name: Split name (e.g., "split40", "split50", "split70")
        split_type: Either "train" or "test"
        run_dir: Path to the run directory
    
    Returns:
        Set of UniProt IDs in the split.
    """
    split_file = run_dir / "data" / "splits" / f"{split_name}_{split_type}.txt"
    
    if not split_file.exists():
        raise FileNotFoundError(f"Split file not found: {split_file}")
    
    with open(split_file) as f:
        return {line.strip() for line in f if line.strip()}


def load_embedding_ids(run_dir: Path) -> Set[str]:
    """
    Load embedding IDs from the ids.txt file.
    
    Args:
        run_dir: Path to the run directory
    
    Returns:
        Set of UniProt IDs with embeddings.
    """
    ids_file = run_dir / "embeddings" / "esm2_t33_650M" / "ids.txt"
    
    if not ids_file.exists():
        raise FileNotFoundError(f"Embedding IDs not found: {ids_file}")
    
    with open(ids_file) as f:
        return {line.strip() for line in f if line.strip()}


def assert_split_integrity(run_dir: Path, split_name: str) -> None:
    """
    Assert that a split is valid and matches the supervised_eligible manifest.
    
    Checks:
    1. len(train) + len(test) == len(supervised_eligible)
    2. train ∩ test == ∅ (disjoint)
    3. train ∪ test == supervised_eligible (full coverage)
    
    Args:
        run_dir: Path to the run directory
        split_name: Split name (e.g., "split40")
    
    Raises:
        AssertionError: If any check fails.
    """
    manifest = load_manifest("supervised_eligible", run_dir)
    train = load_split(split_name, "train", run_dir)
    test = load_split(split_name, "test", run_dir)
    
    # Check 1: Count matches
    expected_total = len(manifest)
    actual_total = len(train) + len(test)
    assert actual_total == expected_total, (
        f"Split count mismatch for {split_name}: "
        f"{actual_total} != {expected_total} (expected)"
    )
    
    # Check 2: Train and test are disjoint
    overlap = train & test
    assert len(overlap) == 0, (
        f"Train/test overlap in {split_name}: {len(overlap)} IDs overlap. "
        f"Sample: {list(overlap)[:5]}"
    )
    
    # Check 3: Full coverage
    combined = train | test
    missing_from_manifest = manifest - combined
    extra_in_split = combined - manifest
    
    assert combined == manifest, (
        f"Split coverage mismatch for {split_name}: "
        f"{len(missing_from_manifest)} IDs missing from split, "
        f"{len(extra_in_split)} extra IDs in split"
    )
    
    print(f"[membership] ✓ {split_name} integrity verified: "
          f"{len(train)} train + {len(test)} test = {len(manifest)} total")


def assert_embedding_coverage(run_dir: Path, manifest_name: str = "domain_E001") -> None:
    """
    Assert that all manifest IDs have corresponding embeddings.
    
    Args:
        run_dir: Path to the run directory
        manifest_name: Manifest to check against (default: domain_E001)
    
    Raises:
        AssertionError: If IDs are missing from embeddings.
    """
    manifest = load_manifest(manifest_name, run_dir)
    embedding_ids = load_embedding_ids(run_dir)
    
    missing = manifest - embedding_ids
    
    assert len(missing) == 0, (
        f"Embedding coverage error: {len(missing)} IDs from {manifest_name} "
        f"are missing from embeddings. Sample: {list(missing)[:5]}"
    )
    
    print(f"[membership] ✓ Embedding coverage verified: "
          f"all {len(manifest)} {manifest_name} IDs have embeddings")


def assert_no_orphans(run_dir: Path) -> None:
    """
    Assert that no orphan IDs exist (IDs in embeddings but not in manifests).
    
    Args:
        run_dir: Path to the run directory
    
    Raises:
        AssertionError: If orphan IDs exist.
    """
    embedding_ids = load_embedding_ids(run_dir)
    
    # Load domain_E001 manifest (should contain all embedded sequences)
    try:
        # Try to load from run dir first
        domain_coords = run_dir / "data" / "domains" / "domain_coords_E001.tsv"
        if domain_coords.exists():
            import pandas as pd
            coords_df = pd.read_csv(domain_coords, sep='\t')
            domain_ids = set(coords_df['uniprot_id'].tolist())
        else:
            # Fallback to project-level domain coords
            from .run_manager import PROJECT_ROOT
            domain_coords = PROJECT_ROOT / "data" / "domains" / "domain_coords_E001.tsv"
            if domain_coords.exists():
                import pandas as pd
                coords_df = pd.read_csv(domain_coords, sep='\t')
                domain_ids = set(coords_df['uniprot_id'].tolist())
            else:
                print("[membership] ⚠ Cannot check orphans: domain_coords not found")
                return
    except Exception as e:
        print(f"[membership] ⚠ Cannot check orphans: {e}")
        return
    
    orphans = embedding_ids - domain_ids
    
    assert len(orphans) == 0, (
        f"Orphan IDs detected: {len(orphans)} IDs in embeddings "
        f"but not in domain_coords. Sample: {list(orphans)[:5]}"
    )
    
    print(f"[membership] ✓ No orphan IDs: all {len(embedding_ids)} "
          f"embedding IDs exist in domain_coords")


def assert_manifest_hierarchy(run_dir: Path) -> None:
    """
    Assert that manifest hierarchy is valid:
    - domain_E001 >= supervised_eligible
    - supervised_eligible == sum of all splits
    
    Args:
        run_dir: Path to the run directory
    
    Raises:
        AssertionError: If hierarchy is invalid.
    """
    try:
        domain_E001 = load_manifest("domain_E001", run_dir)
    except FileNotFoundError:
        print("[membership] ⚠ domain_E001 manifest not found, skipping hierarchy check")
        return
    
    try:
        supervised = load_manifest("supervised_eligible", run_dir)
    except FileNotFoundError:
        print("[membership] ⚠ supervised_eligible manifest not found, skipping hierarchy check")
        return
    
    # supervised_eligible should be a subset of domain_E001
    assert supervised.issubset(domain_E001), (
        f"Manifest hierarchy error: supervised_eligible contains "
        f"{len(supervised - domain_E001)} IDs not in domain_E001"
    )
    
    print(f"[membership] ✓ Manifest hierarchy verified: "
          f"domain_E001 ({len(domain_E001)}) >= supervised_eligible ({len(supervised)})")


def get_dataset_summary(run_dir: Path) -> Dict[str, int]:
    """
    Get a summary of dataset sizes from manifests.
    
    Args:
        run_dir: Path to the run directory
    
    Returns:
        Dictionary with dataset counts.
    """
    summary = {}
    
    manifest_names = [
        "domain_E001",
        "supervised_eligible",
    ]
    
    for name in manifest_names:
        try:
            ids = load_manifest(name, run_dir)
            summary[name] = len(ids)
        except FileNotFoundError:
            summary[name] = None
    
    # Add split counts
    for threshold in [40, 50, 70]:
        split_name = f"split{threshold}"
        try:
            train = load_split(split_name, "train", run_dir)
            test = load_split(split_name, "test", run_dir)
            summary[f"{split_name}_train"] = len(train)
            summary[f"{split_name}_test"] = len(test)
        except FileNotFoundError:
            summary[f"{split_name}_train"] = None
            summary[f"{split_name}_test"] = None
    
    return summary


def validate_all(run_dir: Path, verbose: bool = True) -> List[str]:
    """
    Run all membership validations.
    
    Args:
        run_dir: Path to the run directory
        verbose: Print progress messages
    
    Returns:
        List of error messages (empty if all pass).
    """
    errors = []
    
    if verbose:
        print(f"\n[membership] Validating membership for: {run_dir}")
        print("=" * 60)
    
    # 1. Check manifest hierarchy
    try:
        assert_manifest_hierarchy(run_dir)
    except AssertionError as e:
        errors.append(str(e))
    except FileNotFoundError as e:
        if verbose:
            print(f"[membership] ⚠ Skipping hierarchy check: {e}")
    
    # 2. Check split integrity for all splits
    for threshold in [40, 50, 70]:
        try:
            assert_split_integrity(run_dir, f"split{threshold}")
        except AssertionError as e:
            errors.append(str(e))
        except FileNotFoundError as e:
            if verbose:
                print(f"[membership] ⚠ Skipping split{threshold} check: {e}")
    
    # 3. Check embedding coverage
    try:
        assert_embedding_coverage(run_dir)
    except AssertionError as e:
        errors.append(str(e))
    except FileNotFoundError as e:
        if verbose:
            print(f"[membership] ⚠ Skipping embedding coverage check: {e}")
    
    # 4. Check for orphans
    try:
        assert_no_orphans(run_dir)
    except AssertionError as e:
        errors.append(str(e))
    
    if verbose:
        print("=" * 60)
        if errors:
            print(f"[membership] ✗ FAILED with {len(errors)} error(s)")
            for err in errors:
                print(f"  - {err}")
        else:
            print("[membership] ✓ All membership checks passed")
    
    return errors


if __name__ == "__main__":
    import argparse
    from .run_manager import get_run_dir, get_current_run
    
    parser = argparse.ArgumentParser(description="Validate dataset membership")
    parser.add_argument("--run-dir", type=str, help="Run directory path or ID")
    parser.add_argument("--summary", action="store_true", help="Print dataset summary")
    
    args = parser.parse_args()
    
    if args.run_dir:
        run_dir = Path(args.run_dir)
        if not run_dir.exists():
            run_dir = get_run_dir(args.run_dir)
    else:
        run_dir = get_current_run()
    
    if args.summary:
        summary = get_dataset_summary(run_dir)
        print("\nDataset Summary:")
        for key, value in summary.items():
            print(f"  {key}: {value}")
    else:
        errors = validate_all(run_dir)
        exit(1 if errors else 0)

