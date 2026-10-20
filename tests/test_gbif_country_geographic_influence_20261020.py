"""Recompute all posthoc region effects from frozen cells, never from readout text."""
from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "provenance/reviewer_defenses/results/gbif_country_20261020/gbif_country_primary_strata.json"
RESULT = ROOT / "provenance/reviewer_defenses/results/gbif_country_20261020/geographic_dependence_posthoc_20261020.json"
SPEC = importlib.util.spec_from_file_location(
    "gbif_geo_influence", ROOT / "scripts/analyze_gbif_country_geographic_influence_20261020.py")
mod = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mod)


def test_frozen_output_rebuild_matches_exactly():
    actual = mod.build(json.loads(SOURCE.read_text(encoding="utf-8")))
    archived = json.loads(RESULT.read_text(encoding="utf-8"))
    assert actual == archived
    x = actual["results"]["baseline"]
    assert (x["candidate_strata"], x["candidate_cells"], x["added_cells"],
            x["positive_cells"], x["observed_added_positive"],
            x["informative_strata"]) == (255, 3188, 1527, 48, 26, 21)
    assert math.isclose(x["excess_pp"], 3.912639308, abs_tol=1e-8)
    assert math.isclose(x["conditional_geographic_upper_tail"], 0.283249574739, abs_tol=1e-10)


def test_posthoc_region_sign_and_information_boundaries():
    x = mod.build(json.loads(SOURCE.read_text(encoding="utf-8")))["results"]
    assert x["exclude_caribbean"]["positive_cells"] == 33
    assert x["exclude_caribbean"]["excess_pp"] < 0
    assert x["exclude_caribbean_and_central_america"]["excess_pp"] < 0
    assert x["exclude_level1_8"]["excess_pp"] < 0
    assert x["only_level1_8"]["excess_pp"] > 0
    assert all(row["informative_strata"] < 50 for row in x.values())


def test_exact_hypergeometric_toy():
    assert mod.exact_upper_tail([{"N": 3, "A": 1, "K": 1}], 1) == pytest.approx(1/3)
    assert mod.exact_upper_tail([{"N": 3, "A": 1, "K": 1}], 0) == pytest.approx(1)


def test_validate_rejects_duplicate_country_in_frozen_row():
    frozen = json.loads(SOURCE.read_text(encoding="utf-8"))
    frozen[0]["cells"].append(dict(frozen[0]["cells"][0]))
    with pytest.raises(ValueError):
        mod.validate(frozen)


def test_leave_one_country_and_genus_include_all_groups():
    z = json.loads(RESULT.read_text(encoding="utf-8"))
    frozen = json.loads(SOURCE.read_text(encoding="utf-8"))
    assert len(z["leave_one_country_out"]) == len({c["country"] for x in frozen for c in x["cells"]})
    assert len(z["leave_one_genus_out"]) == len({x["genus"] for x in frozen})
    assert sum(x["omitted_positive_cells"] for x in z["leave_one_country_out"]) == 48
    assert sum(x["omitted_added_positive_cells"] for x in z["leave_one_country_out"]) == 26

def test_exact_rational_rounding_regression():
    from fractions import Fraction
    assert mod.round_nonnegative_fraction(Fraction(1, 2 * 10**12), 12) == 0
    assert mod.round_nonnegative_fraction(Fraction(3, 2 * 10**12), 12) == 2e-12
    frozen = json.loads(SOURCE.read_text(encoding="utf-8"))
    ie = mod.summarize(frozen, omit_countries={"IE"})
    assert ie["expected_added_positive"] == 23.577761617467
