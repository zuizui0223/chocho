"""Synthetic + frozen-margin checks for 2026 country recording-turnover diagnostic."""
from __future__ import annotations

import importlib.util
import json
from fractions import Fraction
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "scripts" / "analyze_gbif_country_record_turnover_20261020.py"
FROZEN = ROOT / "provenance/reviewer_defenses/results/gbif_country_20261020/gbif_country_primary_strata.json"
COUNTS = ROOT / "data/exploratory/gbif_country_20261020/gbif_country_counts.csv"
MANIFEST = ROOT / "data/exploratory/gbif_country_20261020/original_239_taxon_manifest.json"
spec = importlib.util.spec_from_file_location("country_record_turnover", SRC)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def test_conditional_reference_matches_hand_computation():
    rows = [
        {"N": 4, "A": 2, "K": 1, "O": 1, "E": Fraction(1, 2),
         "genus": "Example", "species": "Example alpha",
         "positive_countries": ["AA"], "added_positive_countries": ["AA"]}
    ]
    r = mod.result_for_direction(rows, "toy", draws=19)
    assert r["positive_record_cells"] == 1
    assert r["observed_added_positive"] == 1
    assert r["conditional_expected_added"] == .5
    assert r["excess_percentage_points"] == 50
    assert r["exact_conditional_geographic_upper_tail"] == .5
    assert r["status"] == "BELOW_PROTOCOL_INFORMATION_THRESHOLD"


def test_reporting_background_excludes_target_and_uses_both_windows():
    from collections import defaultdict
    states = [{"species": "A a", "genus": "A", "stratum": "8",
               "cells": [{"country": "CU", "added": 1, "positive": 1},
                         {"country": "DO", "added": 0, "positive": 1}]}]
    b = defaultdict(set)
    for w in mod.WINDOWS:
        b["core", w, "CU"] = {"A a", "B b", "C c", "D d"}
    b["core", "2010_2017", "DO"] = {"A a", "B b", "C c", "D d"}
    b["core", "2018_2025", "DO"] = {"A a", "B b"}
    vals = {
        ("core", "2018_2025", "A a", "CU"): 2,
        ("core", "2010_2017", "A a", "DO"): 1,
    }
    out, audit = mod.states_for(states, vals, b, "core", 3, 1)
    assert [x["country"] for x in out[0]["cells"]] == ["CU"]
    assert out[0]["cells"][0]["state"] == "01"
    assert audit["paired_effort_candidate_cells"] == 1


def test_directional_conditioning_recomputes_exposure_margin():
    data = [{"species": "A a", "genus": "A", "stratum": "8",
             "cells": [
                 {"country": "AA", "added": 1, "early": 0, "late": 1},
                 {"country": "BB", "added": 0, "early": 0, "late": 0},
                 {"country": "CC", "added": 1, "early": 1, "late": 0},
                 {"country": "DD", "added": 0, "early": 1, "late": 1},
             ]}]
    forward = mod.build_direction_rows(data, "late_only")
    backwards = mod.build_direction_rows(data, "early_only")
    assert len(forward) == 1
    assert (forward[0]["N"], forward[0]["A"], forward[0]["K"], forward[0]["O"]) == (2, 1, 1, 1)
    assert len(backwards) == 1
    assert (backwards[0]["N"], backwards[0]["A"], backwards[0]["K"], backwards[0]["O"]) == (2, 1, 1, 1)


def test_source_and_primary_margin_frozen_without_fitting():
    assert FROZEN.exists() and COUNTS.exists() and MANIFEST.exists()
    strata, counts, background, source_info = mod.read_sources(COUNTS, MANIFEST, FROZEN)
    assert len(strata) == 255
    assert source_info["source_taxa"] == 229
    assert source_info["country_rows_including_excluded_aliases"] == 23433
    assert source_info["source_record_totals"]["all_accepted_bases"] == 33817891
    assert sum(x["K"] for x in strata) == 48
    paired, audit = mod.states_for(strata, counts, background, "all_accepted_bases", 3, 1)
    assert audit["paired_effort_candidate_cells"] <= 3188
    assert all(c["state"] in {"00", "01", "10", "11"}
               for row in paired for c in row["cells"])
