"""The frozen examples are artificial mass-ledger mechanics, never insect results."""
from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
PATH = ROOT / "scripts/reconstruct_kyoto_late_instar_intake.py"
SPEC = importlib.util.spec_from_file_location("kyoto_intake", PATH)
assert SPEC is not None and SPEC.loader is not None
kyoto = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(kyoto)


def _row(ident, day, *, stage, given="", taken="", n=1, memo=""):
    return {
        "id": ident, "date": f"2016/7/{day}", "stage": stage,
        "given": str(given), "taken": str(taken),
        "mass": "0.02", "no.individuals": str(n), "memo": memo,
    }


def synthetic():
    rows = [
        {"id": "j", "date": "2016/6/27", "stage": "1",
         "given": "0.3", "taken": "", "mass": "0.002", "no.individuals": "10", "memo": ""},
        {"id": "j", "date": "2016/6/29", "stage": "2",
         "given": "0.4", "taken": "0.2", "mass": "0.02", "no.individuals": "10", "memo": ""},
        _row("j1", 3, stage="3", given="0.7"),
        _row("j1", 5, stage="4", given="0.6", taken="0.3"),
        _row("j1", 7, stage="5", given="0.5", taken="0.1"),
        _row("j1", 9, stage="pp", taken="0.2"),
        _row("j1", 11, stage="p"),
        _row("j2", 3, stage="3", given="0.1"),
        _row("j2", 5, stage="4", given="0.2", taken="0.3"),  # impossible; skip
        _row("j2", 7, stage="5", given="0.2", taken="0.1"),
        _row("j2", 9, stage="p", taken="0.1"),
    ]
    return rows


def test_dry_calibration_preserves_units_and_zero():
    assert kyoto.fresh_g_to_dry_mg(0) == 0.0
    assert kyoto.fresh_g_to_dry_mg(0.5) == pytest.approx(
        math.exp(-1.48 + 0.90 * math.log(500))
    )
    with pytest.raises(ValueError):
        kyoto.fresh_g_to_dry_mg(-1)
    with pytest.raises(ValueError):
        kyoto.fresh_g_to_dry_mg(float("nan"))


def test_exact_two_day_previous_given_current_taken_with_negative_guard():
    x = kyoto.reconstruct(synthetic(), species="A")
    assert x["founder_neonate_group_start"] == "2016-06-27"
    assert x["n_single_larva_ids"] == 2
    assert x["n_pupa_coded_ids"] == 2
    assert x["n_intervals_included"] == 5
    assert x["skipped_adjacent_pairs"]["negative_dry_consumption"] == 1
    assert x["median_day_since_hatch_half_observed_late_intake"] is not None
    expected_first = kyoto.fresh_g_to_dry_mg(0.7) - kyoto.fresh_g_to_dry_mg(0.3)
    assert x["calendar_series"][0]["day_since_founder_group_start"] == 8
    assert x["calendar_series"][0]["sum_dry_mg_across_tracked_individuals"] == pytest.approx(expected_first, abs=1e-5)
    assert x["status"] == "SOURCE_DESCRIBED_LAB_CONSUMPTION_RECONSTRUCTION_NOT_COMPETITION"


def test_redundant_row_excludes_both_adjacent_intervals():
    rows = synthetic()
    for row in rows:
        if row["id"] == "j1" and row["date"] == "2016/7/5":
            row["memo"] = "redundant"
    x = kyoto.reconstruct(rows, species="A")
    assert x["n_intervals_included"] == 3
    assert x["skipped_adjacent_pairs"]["memo_redundant_adjacent"] == 2


def test_duplicate_day_must_fail():
    rows = synthetic()
    rows[-1]["date"] = rows[-2]["date"]
    with pytest.raises(ValueError, match="duplicate identity/calendar-day"):
        kyoto.reconstruct(rows, species="A")


def test_bad_stage_must_fail():
    rows = synthetic()
    rows[0]["stage"] = "sixth"
    with pytest.raises(ValueError, match="unrecognized"):
        kyoto.reconstruct(rows, species="A")


def test_no_group_or_fabricated_origin_fails():
    rows = synthetic()
    rows[0]["no.individuals"] = "1"
    with pytest.raises(ValueError, match="10-neonate founder"):
        kyoto.reconstruct(rows, species="A")
