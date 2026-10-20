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
         "positive_countries": ["AA"], "added_positive_countries": ["AA"],
         "positive_records": [{"species": "Example alpha", "country": "AA",
                              "botanical_level1": "8", "added_host_country": 1,
                              "early_record_count": 0, "late_record_count": 1,
                              "other_taxa_early": 3, "other_taxa_late": 3}]}
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
                 {"country": "AA", "added": 1, "early": 0, "late": 1, "early_count": 0, "late_count": 1, "other_early": 3, "other_late": 3},
                 {"country": "BB", "added": 0, "early": 0, "late": 0, "early_count": 0, "late_count": 0, "other_early": 3, "other_late": 3},
                 {"country": "CC", "added": 1, "early": 1, "late": 0, "early_count": 1, "late_count": 0, "other_early": 3, "other_late": 3},
                 {"country": "DD", "added": 0, "early": 1, "late": 1, "early_count": 1, "late_count": 1, "other_early": 3, "other_late": 3},
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


def test_common_random_numbers_make_identical_scopes_comparable():
    from fractions import Fraction
    row = {"N": 4, "A": 2, "K": 1, "O": 1, "E": Fraction(1, 2),
           "genus": "Example", "species": "Example alpha",
           "positive_countries": ["AA"], "added_positive_countries": ["AA"],
         "positive_records": [{"species": "Example alpha", "country": "AA",
                              "botanical_level1": "8", "added_host_country": 1,
                              "early_record_count": 0, "late_record_count": 1,
                              "other_taxa_early": 3, "other_taxa_late": 3}]}
    a = mod.result_for_direction([row], "accepted/late_only", draws=99)
    b = mod.result_for_direction([row], "core/late_only", draws=99)
    assert a["genus_bootstrap_95_interval_pp"] == b["genus_bootstrap_95_interval_pp"]



def test_published_record_turnover_summary_has_no_unearned_inference():
    path = ROOT / "provenance/reviewer_defenses/results/gbif_country_20261020/gbif_country_record_turnover_summary_20261020.json"
    report = json.loads(path.read_text(encoding="utf-8"))
    assert report["source"]["source_taxa"] == 229
    assert report["source"]["country_rows_including_excluded_aliases"] == 23433
    assert len(report["all_settings"]) == 12
    primary = next(s for s in report["all_settings"] if s["setting"] == report["primary_focal"])
    states = primary["states"]
    assert sum(states["states_added"].values()) + sum(states["states_control"].values()) == 2684
    assert primary["late_only"]["positive_record_cells"] == 16
    assert primary["late_only"]["observed_added_positive"] == 9
    assert primary["late_only"]["informative_strata"] == 10
    assert primary["early_only"]["positive_record_cells"] == 7
    assert all(s[key]["status"] == "BELOW_PROTOCOL_INFORMATION_THRESHOLD"
               for s in report["all_settings"] for key in ("late_only","early_only"))
    assert report["regional_posthoc_diagnostics"]["exclude_caribbean"]["late_only"]["excess_percentage_points"] < 0



def test_positive_transition_details_preserve_species_country_and_record_counts():
    state = [{"species": "A a", "genus": "A", "stratum": "8", "cells": [
        {"country": "AA", "added": 1, "early": 0, "late": 1,
         "early_count": 0, "late_count": 3, "other_early": 4, "other_late": 4},
        {"country": "BB", "added": 0, "early": 0, "late": 0,
         "early_count": 0, "late_count": 0, "other_early": 4, "other_late": 4}
    ]}]
    r = mod.result_for_direction(mod.build_direction_rows(state, "late_only"),
                                 "accepted/late_only", draws=19)
    assert r["positive_record_cells"] == len(r["positive_species_country_details"]) == 1
    cell = r["positive_species_country_details"][0]
    assert (cell["species"], cell["country"], cell["added_host_country"]) == ("A a", "AA", 1)
    assert cell["late_record_count"] == 3
    assert r["positive_window_single_record_cells"] == 0


def test_posthoc_strict_late_record_depth_does_not_reclassify_earlier_records():
    states = [{"species": "A a", "genus": "A", "stratum": "8", "cells": [
        {"country": "AA", "added": 1, "early": 0, "late": 1,
         "early_count": 0, "late_count": 3, "other_early": 4, "other_late": 4},
        {"country": "BB", "added": 0, "early": 0, "late": 0,
         "early_count": 0, "late_count": 0, "other_early": 4, "other_late": 4},
        {"country": "CC", "added": 1, "early": 1, "late": 1,
         "early_count": 1, "late_count": 6, "other_early": 4, "other_late": 4},
        {"country": "DD", "added": 0, "early": 0, "late": 1,
         "early_count": 0, "late_count": 1, "other_early": 4, "other_late": 4},
    ]}]
    depth2 = mod.strict_late_record_depth_rows(states, min_late_records=2)
    depth5 = mod.strict_late_record_depth_rows(states, min_late_records=5)
    assert (depth2[0]["N"], depth2[0]["A"], depth2[0]["K"], depth2[0]["O"]) == (3, 1, 1, 1)
    assert (depth5[0]["N"], depth5[0]["A"], depth5[0]["K"], depth5[0]["O"]) == (3, 1, 0, 0)
    assert "CC" not in depth2[0]["positive_countries"]
