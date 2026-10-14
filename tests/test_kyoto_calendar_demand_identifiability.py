"""Offline tests use fake records solely to guard parsers and inference boundaries."""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

PATH = Path(__file__).resolve().parents[1] / "scripts/audit_kyoto_calendar_demand_identifiability.py"
SPEC = importlib.util.spec_from_file_location("kyoto_demand_calendar", PATH)
assert SPEC is not None and SPEC.loader is not None
mod = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mod)


def make_cages():
    rows = []
    for cage in range(1, 31):
        mixed = cage <= 20
        a_present = cage <= 25
        s_present = mixed or cage >= 26
        for d in range(23):
            row = {k: "0" for k in mod.EXPECTED_LARVAL_COLUMNS}
            row.update({
                "plot.no": str(cage), "treat.no": str(cage),
                "s.density": str(4 if s_present else 0),
                "a.density": str(4 if a_present else 0),
                "date": f"2012-09-{d+1:02d}", "day": str(d),
                "a1": str(4 if a_present and d < 3 else 0),
                "a3": str(4 if a_present and 3 <= d <= 10 else 0),
                "s1": str(4 if s_present and d < 5 else 0),
                "s3": str(4 if s_present and 5 <= d <= 12 else 0),
            })
            rows.append(row)
    return list(mod.EXPECTED_LARVAL_COLUMNS), rows


def test_original_cage_unit_and_species_late_instar_lag_are_distinct():
    cols, rows = make_cages()
    a = mod.summarize_stages(cols, rows)
    assert a["n_source_cages"] == 30
    assert a["n_initially_mixed_cages"] == 20
    assert a["n_mixed_cages_both_species_observed_late"] == 20
    assert a["observed_pairwise_S_minus_A_first_late_day_median"] == 2
    assert a["paired_late_onset_sign_counts"] == {
        "S_first_late_earlier": 0,
        "same_first_observed_day": 0,
        "A_first_late_earlier": 20,
    }
    assert a["mixed_cage_census_overlap_counts"]["both_late"] == 20 * 6
    assert a["n_cage_by_initial_species_presence"] == {
        "initial_A_0_S_1": 5, "initial_A_1_S_0": 5, "initial_A_1_S_1": 20
    }
    assert "not food intake" in a["boundary"]


def test_original_stage_audit_rejects_duplicate_dates_per_cage():
    cols, rows = make_cages()
    rows[1]["day"] = rows[0]["day"]
    with pytest.raises(ValueError, match="duplicate census days"):
        mod.summarize_stages(cols, rows)


def test_original_stage_audit_rejects_unknown_late_stage_count():
    cols, rows = make_cages()
    rows[7]["a3"] = ""
    with pytest.raises(ValueError, match="invalid 3rd"):
        mod.summarize_stages(cols, rows)


def test_food_dates_and_trial_ids_are_not_mistaken_for_true_dry_food_intake():
    cols = list(mod.EXPECTED_FOOD_COLUMNS)
    rows = [
        {"id": "1", "date": "2012-06-01", "taken": "3", "given": "7", "mass": "4", "no.individuals": "10", "memo": ""},
        {"id": "1", "date": "2012-06-02", "taken": "2", "given": "6", "mass": "6", "no.individuals": "10", "memo": "redundant"},
        {"id": "2", "date": "2012-06-01", "taken": "1", "given": "6", "mass": "3", "no.individuals": "1", "memo": ""},
    ]
    out = mod.summarize_food(cols, rows)
    assert out["raw_rows"] == 3
    assert out["unique_food_trial_ids"] == 2
    assert out["trial_ids_with_repeated_dates"] == 1
    assert out["group_size_counts"] == {"1": 1, "10": 2}
    assert out["memo_counts"]["redundant"] == 1
    assert "NOT interpretable" in out["boundary"]


def test_food_parse_fails_closed_for_missing_sample_identity():
    cols = list(mod.EXPECTED_FOOD_COLUMNS)
    row = {k: "1" for k in cols}
    row["id"] = ""
    with pytest.raises(ValueError, match="without larva/group ID"):
        mod.summarize_food(cols, [row])


def test_food_parse_fails_for_missing_required_field():
    cols = ["id", "date", "given"]
    with pytest.raises(ValueError, match="required published"):
        mod.summarize_food(cols, [{k: "1" for k in cols}])
