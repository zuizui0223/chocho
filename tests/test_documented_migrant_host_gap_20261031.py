"""Recompute documented-migrant island gap diagnostics from pinned public sources.

IMPORTANT: Not-listed is unknown; GBIF record gaps are not observed flight paths.
"""
from __future__ import annotations
import importlib.util
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
D=ROOT/"provenance/reviewer_defenses/results/butterfly_maritime_trait_gap_20261031"
SPEC=importlib.util.spec_from_file_location("chocho_migrant_host_gap",D/"analyze_documented_migrant_subset.py")
assert SPEC is not None and SPEC.loader is not None
MOD=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MOD)

def close(a,b):
    if a is None or b is None:
        assert a is None and b is None
    else:
        assert math.isclose(a,b,rel_tol=0,abs_tol=1e-12),(a,b)

def test_migrant_public_193_screen_and_1500_bootstrap_reproduces():
    saved=json.loads((D/"documented_migrant_host_gap_results.json").read_text(encoding="utf8"))
    got=MOD.analyze(draws=1500)
    assert got["status"]==saved["status"]=="DESCRIPTIVE_MIGRANT_SUBSET_WEAK_UNSTABLE_HOST_PROXY_ASSOCIATION"
    assert got["n_source_225"]==saved["n_source_225"]==225
    assert got["n_source_documented_candidate"]==saved["n_source_documented_candidate"]==111
    assert got["documented_subset"]["species"]==saved["documented_subset"]["species"]==106
    assert got["documented_subset"]["prior_record_gaps_ge500km"]==saved["documented_subset"]["prior_record_gaps_ge500km"]==119
    assert got["all_193_screened"]==saved["all_193_screened"]==193
    for k in ("documented_subset","non_listed_unknown_coverage_NOT_A_CONTROL"):
        assert got[k]["species"]==saved[k]["species"]
        for metric, value in saved[k]["correlations"].items():
            close(value,got[k]["correlations"][metric])
    for n,section in saved["min_events_sensitivity"].items():
        assert section["species"]==got["min_events_sensitivity"][n]["species"]
        for k,value in section["correlations"].items():
            close(value,got["min_events_sensitivity"][n]["correlations"][k])
    for k,info in saved["species_resampled_uncertainty"].items():
        actual=got["species_resampled_uncertainty"][k]
        assert info["successful_draws"]==actual["successful_draws"]==1500
        for a,b in zip(info["species_bootstrap_95"],actual["species_bootstrap_95"]):
            close(a,b)
    a=saved["retrospective_5_year_followup"]
    b=got["retrospective_5_year_followup"]
    assert a["same_migrant_species_alien_vs_native_comparable_species"]==b["same_migrant_species_alien_vs_native_comparable_species"]==12
    close(a["same_species_equal_mean_reindex_fraction_difference"],b["same_species_equal_mean_reindex_fraction_difference"])
    assert a["same_species_direction"]==b["same_species_direction"]

def test_no_false_nonmigrant_or_earlier_host_claims():
    result=json.loads((D/"documented_migrant_host_gap_results.json").read_text(encoding="utf8"))
    x=json.loads((D/"migration_gbif_predicate_audit.json").read_text(encoding="utf8"))
    assert x["status"]=="SOURCE_RECOVERED_INDIRECT_MIGRANT_PROXY"
    assert x["GBIF_raw_predicate_taxon_values"]==564
    assert x["GBIF_duplicate_predicate_values"]==2
    assert x["GBIF_predicate_taxon_keys"]==562
    assert x["candidate_225_exact_matched"]==111
    assert result["migration_classification_warning"].startswith("not listed = migration status unknown")
    assert "NOT_A_CONTROL" in result
    assert any("postdate" in s for s in result["limits"])
    assert result["documented_subset"]["correlations"]["host_fraction_vs_500km_gap"]<0.2
    assert result["min_events_sensitivity"]["5"]["correlations"]["host_fraction_vs_500km_gap"]<0.05
