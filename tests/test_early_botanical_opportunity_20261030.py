"""Exact support reconciliation for archived 51-island botanical recording comparison."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import math
import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "provenance/reviewer_defenses/results/early_botanical_observation_support_20261030/audit_early_botanical_opportunity.py"
spec = importlib.util.spec_from_file_location("early_botanical_opportunity", SCRIPT)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def test_coverage_max_never_sums_overlapping_ports():
    rows = [
        {"port_id":"1","landmass_id":"7","early_plant_records":"3","all_plant_records":"15"},
        {"port_id":"2","landmass_id":"7","early_plant_records":"5","all_plant_records":"25"},
        {"port_id":"3","landmass_id":"8","early_plant_records":"0","all_plant_records":"10"}
    ]
    x=mod.coverage_by_island(rows)
    assert x["7"]["early_max_per_port"]==5
    assert x["7"]["early_port_record_sum"]==8
    assert x["8"]["early_max_per_port"]==0


def test_fixed_thresholds_preserve_zero_as_detection_nondetection():
    pairs=[
        {"island_a":"7","island_b":"8","min_of_island_early_max":0,
         "shared_new_species":2,"early_matched_butterfly_species":0},
        {"island_a":"7","island_b":"9","min_of_island_early_max":12,
         "shared_new_species":4,"early_matched_butterfly_species":1},
    ]
    assert mod.staged_counts(pairs,0)["shared_new_butterfly_species_island_pairs"]==6
    assert mod.staged_counts(pairs,1)["shared_new_butterfly_species_island_pairs"]==4
    assert mod.staged_counts(pairs,10)["historically_recorded_matching_host_species_island_pairs"]==1
    assert mod.staged_counts(pairs,100)["match_fraction_when_denominator_positive"] is None
    with pytest.raises(ValueError):
        mod.coverage_by_island([{"port_id":"1","landmass_id":"7","early_plant_records":"10","all_plant_records":"5"}])
    with pytest.raises(ValueError):
        mod.unordered_pair("8","8")


def test_full_fixed_51_island_source_reconciliation():
    # The GEB anonymous-stage deliberately omits exploratory source CSVs;
    # ordinary paper CI must instead execute all source-pinned assertions.
    if "_anonymous_stage" in ROOT.parts:
        pytest.skip("Anonymous manuscript export intentionally excludes exploratory island sources")
    r=mod.audit(ROOT)
    saved_path=ROOT / "provenance/reviewer_defenses/results/early_botanical_observation_support_20261030/early_botanical_opportunity_summary.json"
    saved=json.loads(saved_path.read_text(encoding="utf-8"))
    assert saved==r
    assert (r["original_islands"],r["original_ports"],r["original_island_pairs"],
            r["original_shared_new_species_pair_denominator"]) == (51,194,1275,1163)
    s=r["fixed_thresholds"]
    assert [v["island_pairs"] for v in s]==[1275,946,666,231]
    assert [v["shared_new_butterfly_species_island_pairs"] for v in s]==[1163,800,586,189]
    assert [v["historically_recorded_matching_host_species_island_pairs"] for v in s]==[5,5,5,4]
    c=r["ports_early_all_aggregate"]
    assert c["islands_with_no_early_plant_indexed_records"]==7
    assert math.isclose(c["largest_island_fraction_of_early_query_sum"],117356/155002)
    assert len(r["five_early_host_link_coverage"])==5
    assert r["unobserved_plant_not_proven_absent"]
