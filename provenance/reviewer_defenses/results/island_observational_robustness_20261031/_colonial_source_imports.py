"""Load exactly the source-locked GloNAF/Lenzner plant parser from PR #120.

Source-only helper: uses no GBIF butterfly responses. Never duplicate or
reinterpret raw GloNAF taxon species bit positions.
"""
from __future__ import annotations
from pathlib import Path
import sys

SCRIPT_DIR=Path("_external/colonial/scripts")


def _prepare():
    if not SCRIPT_DIR.exists():
        raise RuntimeError("Frozen original GloNAF source scripts not checked out")
    sys.path.insert(0,str(SCRIPT_DIR))


def load_glonaf_helpers():
    _prepare()
    import analyze_matched_glonaf_consumer_history as science
    return science


def load_glonaf_index(original):
    _prepare()
    from analyze_documented_host_specificity import glonaf_index
    return glonaf_index(original)
