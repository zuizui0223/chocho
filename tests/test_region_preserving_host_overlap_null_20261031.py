"""Frozen source-locked species/island fixed margin geographic null check.

CI checks source identities, arithmetic, and published results; it cannot
recover private larval-host associations or re-run author-side evidence here.
"""
from pathlib import Path
import hashlib,json,ast,math
import pytest
ROOT=Path(__file__).resolve().parents[1]
D=ROOT/"provenance/reviewer_defenses/results/region_preserving_host_overlap_20261031"
F=D/"host_overlap_geography_null_audit.json"
CODE=D/"analyze_regional_host_overlap_null.py"

def blob(p):
    data=p.read_bytes()
    return hashlib.sha1(b"blob "+str(len(data)).encode()+b"\0"+data).hexdigest()

def test_source_pins_and_preserved_degrees():
    if "_anonymous_stage" in ROOT.parts:
        pytest.skip("Exploratory private-host source intentionally omitted from anonymous GEB review")
    x=json.loads(F.read_text(encoding="utf8"))
    ast.parse(CODE.read_text(encoding="utf8"))
    for path,sha in x["source_public_git_blobs"].items():
        assert blob(ROOT/path)==sha
    assert x["observation_cohort_n"]==113
    assert x["butterfly_species_count"]==60
    assert x["island_count"]==17
    assert x["observed_indexed_compatible_host_n"]==34
    assert math.isclose(x["observed_indexed_compatible_host_fraction"],34/113)
    assert sum(len(v) for v in x["geographic_near_network_components"])==17
    assert len(x["geographic_near_network_components"])==6
    assert len({id for comp in x["geographic_near_network_components"] for id in comp})==17
    assert all(x["preservation"].values())

def test_global_and_region_marginal_null_not_causal():
    if "_anonymous_stage" in ROOT.parts:
        pytest.skip("Original private HOSTS links not redistributed for blind GEB review")
    x=json.loads(F.read_text(encoding="utf8"))
    g=x["nulls"]["global_species_island_fixed_margins"]
    r=x["nulls"]["within_six_geographic_components_fixed_margins"]
    for n in (g,r):
        assert n["draws"]==1497
        assert len(n["chains"])==3
        assert sum(c["draws"] for c in n["chains"])==n["draws"]
        assert sum(c["draws_ge_observed"] for c in n["chains"])==n["inclusive_tail_ge_observed"]
        assert math.isclose(n["inclusive_fraction"],n["inclusive_tail_ge_observed"]/n["draws"])
        assert math.isclose(n["mean_overlap_fraction"],n["mean_overlap_n"]/113)
        assert n["interval95_overlap_n"][0]<=34<=n["interval95_overlap_n"][1]
        assert all(c["accepted_double_edge_swaps"]>10000 for c in n["chains"])
        assert all(-1<c["lag1_score_autocorrelation"]<1 for c in n["chains"])
    assert g["inclusive_tail_ge_observed"]==49
    assert r["inclusive_tail_ge_observed"]==83
    assert g["inclusive_fraction"]<r["inclusive_fraction"]
    assert 28<g["mean_overlap_n"]<29
    assert 28<r["mean_overlap_n"]<29
    assert x["inference_limits"]["not_preregistered"]
    assert x["inference_limits"]["no_host_indexed_not_true_absence"]
    assert x["inference_limits"]["mcmc_draws_not_independent_samples"]
    assert x["inference_limits"]["global_cross_region_null_is_geographically_unrealistic"]
