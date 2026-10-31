"""Historical country-bilateral LSBCI source/selection/inference gates.

The 42 MB full UNCTAD raw data are downloaded in a separate source-acquisition
workflow, never required in the normal GEB review/CI suite.
"""
import importlib.util
import json
import math
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
DIR = ROOT / "provenance/reviewer_defenses/results/island_observational_robustness_20261031"
FILE = DIR / "full_unctad_lsbci_2010.py"
AUDIT = DIR / "full_lsbci2010_audit.json"
PAIRS = DIR / "full_lsbci2010_island_pairs.csv"


def module():
    spec = importlib.util.spec_from_file_location("full_unctad_lsbci_2010", FILE)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_country_assignment_preserves_multicountry_island():
    m = module()
    islands, mixed = m.fixed_islands()
    assert len(islands) == 50
    assert len(mixed) == 1
    assert sorted(mixed["9"]) == ["BRN", "IDN", "MYS"]
    assert "9" not in islands
    assert m.iso3("United Kingdom of Great Britain and Northern Ireland") == "GBR"
    assert m.iso3("China, Taiwan Province of") == "TWN"
    assert m.iso3("Viet Nam") == "VNM"
    assert m.iso3("Brunei Darussalam") == "BRN"
    assert m.iso3("unmapped fictional country") is None


def test_no_invented_route_or_absence():
    m = module()
    assert "bilateral" in m.RAW_URL.lower() or "LSBCI" in m.RAW_URL
    assert m.PINNED_BLOBS["port_island"] == "0a01831056c2cd59a4bd41bd1dbd4c84d77cf922"
    assert m.RAW_SHA256 == "edab41de33b9a6b4fc568bf454544dedfcbf32117e2bad7ccc77203e1ee6c33e"
    assert math.isclose(m.partial_rank([
        dict(x=1, y=4), dict(x=2, y=1), dict(x=3, y=3), dict(x=4, y=2)
    ], "x", "y"), -0.4, abs_tol=1e-12)


def test_completed_source_receipt_if_present():
    if "_anonymous_stage" in ROOT.parts:
        pytest.skip("Island/port source files not part of anonymous GEB review")
    if not AUDIT.exists() or not PAIRS.exists():
        pytest.skip("Full 2010 acquisition is an independent one-time workflow")
    m = module()
    d = json.loads(AUDIT.read_text(encoding="utf8"))
    pairs = m.read_csv(PAIRS)
    assert d["status"] == "SOURCE_VERIFIED_COUNTRY_BILATERAL_COHORT_NOT_PORT_OD"
    assert d["source_sha256"] == m.RAW_SHA256
    assert len(pairs) == 1275
    eligible = [x for x in pairs if x["status"] == "INDEX_2010_OBSERVED"]
    assert len(eligible) == d["pair_status_counts"]["INDEX_2010_OBSERVED"]
    assert d["sample_size"]["island_pairs"] == len(eligible)
    assert all(x["lsbci_2010"] != "" for x in eligible)
    assert all(x["lsbci_2010"] == "" for x in pairs
               if x["status"] != "INDEX_2010_OBSERVED")
    assert d["eligible_islands"] == 51
    assert d["single_country_islands"] == 50
    # Source outputs are never promoted to port-pair routes by type.
    assert "COUNTRY_BILATERAL" in d["status"]
