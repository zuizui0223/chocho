"""Synthetic fixtures: mathematical shape diagnostics, not real butterflies."""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest

BASE = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(BASE))
SPEC = importlib.util.spec_from_file_location("matched_kyoto", BASE / "analyze_kyoto_matched_age_intake_shape.py")
assert SPEC is not None and SPEC.loader is not None
mod = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mod)


def toy_species(masses_by_id):
    events = [
        {"original_larva_id": larva, "end_day_since_founder_hatch": day,
         "consumed_leaf_dry_mg": value}
        for larva, masses in masses_by_id.items()
        for day, value in zip(mod.PRIMARY_AGES, masses)
    ]
    return {"interval_events": events, "n_single_larva_ids": len(masses_by_id),
            "table1_dry_intake_within_2pct": True}


def toy():
    return {
        "A": toy_species({"a1": [3, 2, 1, 1, 1, 1],
                          "a2": [3, 3, 1, 1, 1, 1],
                          "a3": [3, 2, 2, 1, 1, 1]}),
        "S": toy_species({"s1": [1, 1, 1, 1, 2, 3],
                          "s2": [1, 1, 1, 1, 3, 3],
                          "s3": [1, 1, 1, 2, 2, 3]}),
    }


def test_primary_matched_complete_case_detects_synthetic_temporal_difference():
    r = mod.primary_with_leave_one_out(toy())
    assert r["by_species"]["A"]["included_complete_matched_age_ids"] == 3
    assert r["by_species"]["S"]["included_complete_matched_age_ids"] == 3
    assert sum(r["by_species"]["A"]["each_larva_equal_weight_normalized_age_profile"].values()) == pytest.approx(1)
    assert r["temporal_demand_profile_overlap"] < 1
    assert r["total_variation_between_normalized_profiles"] == pytest.approx(1 - r["temporal_demand_profile_overlap"])
    assert r["A_minus_S_early_share"] > 0
    assert r["leave_one_individual_out"]["sign_consistent_with_full_sample"] is True
    assert r["hypothetical_reversal_condition"]["possible_under_illustrative_scarcity_allocation"] is True
    lo, hi = r["hypothetical_reversal_condition"]["k_A_over_k_S_open_interval"]
    assert lo < 1 < hi


def test_identical_temporal_shapes_yield_zero_tv_and_no_reversal():
    a = {"a1": [2, 2, 2, 2, 2, 2], "a2": [1, 1, 1, 1, 1, 1]}
    s = {"s1": [3, 3, 3, 3, 3, 3], "s2": [2, 2, 2, 2, 2, 2]}
    r = mod.primary_with_leave_one_out({"A": toy_species(a), "S": toy_species(s)})
    assert r["total_variation_between_normalized_profiles"] == pytest.approx(0)
    assert r["A_minus_S_early_share"] == pytest.approx(0)
    assert r["hypothetical_reversal_condition"]["possible_under_illustrative_scarcity_allocation"] is False


def test_complete_case_exclusions_cannot_be_silent():
    source = toy()
    source["A"]["interval_events"] = [
        e for e in source["A"]["interval_events"]
        if not (e["original_larva_id"] == "a1" and e["end_day_since_founder_hatch"] == 18)
    ]
    p = mod.profile_comparison(source, ages=mod.PRIMARY_AGES)
    assert p["by_species"]["A"]["included_complete_matched_age_ids"] == 2
    reason = p["by_species"]["A"]["excluded_ids_and_reason"]
    assert len(reason) == 1
    assert reason[0]["original_id"] == "a1"
    assert reason[0]["missing_ages"] == [18]


def test_duplicate_source_ages_rejected():
    source = toy()
    source["A"]["interval_events"].append(dict(source["A"]["interval_events"][0]))
    with pytest.raises(ValueError, match="duplicate original larva"):
        mod.profile_comparison(source, ages=mod.PRIMARY_AGES)


def test_too_few_complete_case_profiles_fails_explicitly():
    source = toy()
    source["A"]["interval_events"] = [
        e for e in source["A"]["interval_events"]
        if e["original_larva_id"] == "a1"
        or (e["original_larva_id"] in {"a2", "a3"} and e["end_day_since_founder_hatch"] != 18)
    ]
    with pytest.raises(ValueError, match="fewer than two"):
        mod.profile_comparison(source, ages=mod.PRIMARY_AGES)


def test_secondary_windows_are_declared_not_optimized():
    for window in mod.ALT_WINDOWS:
        r = mod.profile_comparison(toy(), ages=window)
        assert r["age_bin_end_days_since_source_species_founder_start"] == list(window)
        assert r["temporal_demand_profile_overlap"] <= 1


def test_no_source_gate_no_promotion():
    with pytest.raises(ValueError, match="no verified"):
        mod.analyze({"status": "SOURCE_BLOCKED_OR_INVALID"})
    x = toy()
    x["A"]["table1_dry_intake_within_2pct"] = False
    with pytest.raises(ValueError, match="Table 1"):
        mod.analyze({"status": "SOURCE_RECONSTRUCTED_DESCRIPTIVE_ONLY",
                     "per_species": x, "source_identity": {}})


def test_published_source_reproduction_receipt_arithmetic_and_claim_limits():
    receipt_path = BASE.parent / "provenance/reviewer_defenses/results/kyoto_matched_age_intake_shape_20261014.json"
    x = json.loads(receipt_path.read_text(encoding="utf-8"))
    c = x["comparison"]
    av = list(c["A_equal_larva_profile_by_age"].values())
    sv = list(c["S_equal_larva_profile_by_age"].values())
    assert c["complete_case_A"] == 9 and c["complete_case_S"] == 9
    assert sum(av) == pytest.approx(1.0)
    assert sum(sv) == pytest.approx(1.0)
    overlap = sum(min(a, b) for a, b in zip(av, sv))
    assert overlap == pytest.approx(c["temporal_overlap"])
    assert 1 - overlap == pytest.approx(c["total_variation"])
    assert sum(av[:3]) - sum(sv[:3]) == pytest.approx(c["early_share_A_minus_S"])
    assert c["toy_k_A_over_k_S_open_interval"][0] == pytest.approx(
        c["S_early_share"] / c["A_early_share"]
    )
    assert c["toy_k_A_over_k_S_open_interval"][1] == pytest.approx(
        (1 - c["S_early_share"]) / (1 - c["A_early_share"])
    )
    assert x["source"]["individual_a_source_md5"] == "08c2bcfb9a896f7eb89bebd74153f769"
    assert x["source"]["individual_s_source_md5"] == "174f0bd6ba056bcd49b348427da9f749"
    assert x["status"].startswith("SOURCE_REPRODUCED_EXPLORATORY")
    assert x["claim_boundary"] and any("no" in str(z).lower() for z in x["claim_boundary"])
