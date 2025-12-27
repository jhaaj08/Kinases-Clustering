#!/usr/bin/env python3
"""
Run Manager - Create and validate run directories.

This module enforces the "no stale output" rule by requiring either:
1. A fresh run directory (new RUN_ID)
2. Explicit --force flag to overwrite

Usage:
    from pipeline.run_manager import init_run, get_run_dir
    
    run_dir = init_run()  # Auto-generated timestamp
    run_dir = init_run(run_id="experiment_v2")  # Named run
    run_dir = init_run(run_id="experiment_v2", force=True)  # Overwrite
"""

import os
import json
import shutil
from datetime import datetime
from pathlib import Path
from typing import Optional


# Project root directory
PROJECT_ROOT = Path(__file__).parent.parent

# Runs directory
RUNS_DIR = PROJECT_ROOT / "runs"


def generate_run_id() -> str:
    """Generate a timestamp-based run ID."""
    return datetime.now().strftime("%Y-%m-%d_%H%M%S")


def init_run(run_id: Optional[str] = None, force: bool = False) -> Path:
    """
    Initialize a new run directory.
    
    Args:
        run_id: Optional run identifier. If None, generates timestamp-based ID.
        force: If True, overwrites existing run directory. Default False.
    
    Returns:
        Path to the run directory.
    
    Raises:
        RuntimeError: If run directory exists and force=False.
    """
    if run_id is None:
        run_id = generate_run_id()
    
    run_dir = RUNS_DIR / run_id
    
    # Enforce "no stale output" rule
    if run_dir.exists() and not force:
        raise RuntimeError(
            f"Run directory '{run_dir}' already exists.\n"
            f"Options:\n"
            f"  1. Use --force to overwrite\n"
            f"  2. Specify a new run ID\n"
            f"  3. Delete the existing directory manually"
        )
    
    # Clean up if forcing overwrite
    if force and run_dir.exists():
        print(f"[run_manager] Removing existing run directory: {run_dir}")
        shutil.rmtree(run_dir)
    
    # Create directory structure
    run_dir.mkdir(parents=True, exist_ok=True)
    
    subdirs = [
        "data/manifests",
        "data/splits",
        "data/domains",
        "embeddings/esm2_t33_650M",
        "results/clustering",
        "results/supervised",
        "results/calibration",
        "results/baselines",
        "results/retrieval",
        "results/layer_comparison",
        "tables",
    ]
    
    for subdir in subdirs:
        (run_dir / subdir).mkdir(parents=True, exist_ok=True)
    
    # Save run configuration
    config = {
        "run_id": run_id,
        "created_at": datetime.now().isoformat(),
        "project_root": str(PROJECT_ROOT),
        "python_version": f"{os.sys.version_info.major}.{os.sys.version_info.minor}.{os.sys.version_info.micro}",
    }
    
    config_file = run_dir / "run_config.json"
    with open(config_file, "w") as f:
        json.dump(config, f, indent=2)
    
    # Update 'current' symlink
    _update_current_symlink(run_dir)
    
    print(f"[run_manager] Initialized run directory: {run_dir}")
    return run_dir


def _update_current_symlink(run_dir: Path) -> None:
    """Update the 'current' symlink to point to the latest run."""
    current_link = RUNS_DIR / "current"
    
    # Remove existing symlink if it exists
    if current_link.is_symlink():
        current_link.unlink()
    elif current_link.exists():
        # It's a real directory, not a symlink - remove it
        shutil.rmtree(current_link)
    
    # Create new symlink (relative path for portability)
    current_link.symlink_to(run_dir.name)
    print(f"[run_manager] Updated symlink: runs/current -> {run_dir.name}")


def get_run_dir(run_id: str) -> Path:
    """
    Get the path to an existing run directory.
    
    Args:
        run_id: The run identifier.
    
    Returns:
        Path to the run directory.
    
    Raises:
        FileNotFoundError: If run directory doesn't exist.
    """
    run_dir = RUNS_DIR / run_id
    
    if not run_dir.exists():
        raise FileNotFoundError(f"Run directory not found: {run_dir}")
    
    return run_dir


def get_current_run() -> Path:
    """
    Get the path to the current (most recent) run directory.
    
    Returns:
        Path to the current run directory.
    
    Raises:
        FileNotFoundError: If no current run exists.
    """
    current_link = RUNS_DIR / "current"
    
    if not current_link.exists():
        raise FileNotFoundError(
            "No current run found. Run 'make all' or 'init_run()' first."
        )
    
    return current_link.resolve()


def list_runs() -> list:
    """List all available run directories."""
    if not RUNS_DIR.exists():
        return []
    
    runs = []
    for item in RUNS_DIR.iterdir():
        if item.is_dir() and not item.is_symlink() and item.name != ".gitkeep":
            config_file = item / "run_config.json"
            if config_file.exists():
                with open(config_file) as f:
                    config = json.load(f)
                runs.append({
                    "run_id": item.name,
                    "created_at": config.get("created_at"),
                    "path": str(item)
                })
    
    # Sort by creation time (newest first)
    runs.sort(key=lambda x: x["created_at"] or "", reverse=True)
    return runs


def validate_run_dir(run_dir: Path) -> dict:
    """
    Validate that a run directory has all required components.
    
    Returns:
        Dictionary with validation results.
    """
    required_files = [
        "run_config.json",
        "data/manifests/supervised_eligible.txt",
        "data/splits/split40_train.txt",
        "data/splits/split40_test.txt",
        "embeddings/esm2_t33_650M/ids.txt",
        "results/manuscript_numbers.json",
    ]
    
    results = {"valid": True, "missing": [], "present": []}
    
    for rel_path in required_files:
        full_path = run_dir / rel_path
        if full_path.exists():
            results["present"].append(rel_path)
        else:
            results["missing"].append(rel_path)
            results["valid"] = False
    
    return results


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Run directory manager")
    parser.add_argument("--init", action="store_true", help="Initialize a new run")
    parser.add_argument("--run-id", type=str, help="Run ID (default: auto-generated)")
    parser.add_argument("--force", action="store_true", help="Force overwrite existing run")
    parser.add_argument("--list", action="store_true", help="List all runs")
    parser.add_argument("--current", action="store_true", help="Show current run")
    
    args = parser.parse_args()
    
    if args.init:
        run_dir = init_run(run_id=args.run_id, force=args.force)
        print(f"Run directory: {run_dir}")
    elif args.list:
        runs = list_runs()
        if runs:
            print("Available runs:")
            for run in runs:
                print(f"  {run['run_id']} (created: {run['created_at']})")
        else:
            print("No runs found.")
    elif args.current:
        try:
            current = get_current_run()
            print(f"Current run: {current}")
        except FileNotFoundError as e:
            print(str(e))
    else:
        parser.print_help()

