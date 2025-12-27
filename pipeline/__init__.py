"""
Pipeline module for reproducible kinase classification experiments.

This module provides:
- Run directory management (run_manager)
- Dataset membership validation (membership)
- Step execution utilities
"""

from .run_manager import init_run, get_run_dir, get_current_run
from .membership import (
    load_manifest,
    assert_split_integrity,
    assert_embedding_coverage,
    assert_no_orphans
)

__version__ = "1.0.0"

