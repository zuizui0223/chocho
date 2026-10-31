"""Reconcile published near-island record holdout snapshot without private raw data.

Public CI verifies locked original *published* pair/botanical sources, logical
support and code syntax, not the empirical predictions from a private Parquet.
"""
from __future__ import annotations
import ast, csv, hashlib, json, math
from collections import defaultdict
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]
B=ROOT/"provenance/reviewer_defenses/results"
F=B/"maritime_host_steppingstone_20261030/near_island_host_opportunity_holdout_20261031.json"
CODE=B/"maritime_host_steppingstone_20261030/recompute_near_island_relay_holdout_20261031.py"
PAIRS=B/"island_host_network_20261030/island_host_network_pairs.csv"
PLANTS=B/"island_dated_host_match_20261030/island_port_early_plant_availability.csv"
LINKS=B/"island_dated_host_match_20261030/island_dated_host_match_links.csv"

def git_blob(p):
    b=p.read_bytes()
    return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()

def test_source_locked_botanical_opportunity_and_risk_support():
    if "_anonymous_stage" in ROOT.parts:
        pytest.skip("Anonymous GEB manuscript excludes exploratory original island sources")
    x=json.loads(F.read_text(encoding="utf8"))
    refs=x["source_public_git_blobs"]
    for p,k in [(PAIRS,"pairs"),(PLANTS,"early_plants"),(LINKS,"host_links")]:
        assert git_blob(p)==refs[k]
    assert x["source_private_parquet_sha256"]=="8362ec4538d28916b2f9aaf98d1512b343b15f8f2eefe060ee61798c9a67406f"
    assert not (ROOT/"work/barriers/global_original_species_island_year.parquet").exists()
    ast.parse(CODE.read_text(encoding="utf8"))
    pairrows=list(csv.DictReader(PAIRS.open(encoding="utf8",newline="")))
    assert len(pairrows)==1275
    near=[r for r in pairrows if float(r["port_distance_km"])<200]
    assert len(near)==x["pair_cohort"]["near_pairs"]==14
    assert sum(int(r["shared_new_species"]) for r in near)==x["pair_cohort"]["shared_new_species_pair_rows"]==41
    n=x["pair_cohort"]
    assert n["islands"]==17 and n["chain_rows"]==21
    assert (n["focal_prior_host_early"]["n"],n["focal_prior_host_early"]["chains"],n["focal_prior_host_early"]["unique_pairs"])==(3,3,2)
    assert (n["other_prior_host_unclassified"]["n"],n["other_prior_host_unclassified"]["chains"])==(38,18)
    assert (n["botanical_opportunity_min100_focal"]["n"],n["botanical_opportunity_min100_focal"]["chains"])==(3,3)
    assert (n["botanical_opportunity_min100_other_host_unclassified"]["n"],n["botanical_opportunity_min100_other_host_unclassified"]["chains"])==(9,3)
    early=defaultdict(list)
    for row in csv.DictReader(PLANTS.open(encoding="utf8",newline="")):
        early[int(row["landmass_id"])].append(int(row["early_plant_records"]))
    e={i:max(v) for i,v in early.items()}
    near_high=[r for r in near if min(e[int(r["island_a"])],e[int(r["island_b"])])>=100]
    assert len(near_high)==9
    assert sum(int(r["shared_new_species"]) for r in near_high)==12
    assert x["guardrails"]["other_unknown_links_not_negative_hosts"]
    assert x["guardrails"]["private_parquet_not_published"]

def test_time_holdout_predictions_not_overstated_as_causation():
    if "_anonymous_stage" in ROOT.parts:
        pytest.skip("Exploratory island prediction snapshot excluded from anonymous GEB manuscript")
    x=json.loads(F.read_text(encoding="utf8"))
    a=x["annual_risk"]
    assert (a["rows"],a["events"],a["butterfly_species"],a["destination_islands"])==(865,91,70,17)
    assert a["unrepeated_prior_near_source"]["events"]+a["repeated_prior_near_source"]["events"]==91
    assert a["unrepeated_prior_near_source"]["rows"]+a["repeated_prior_near_source"]["rows"]==865
    assert a["recent_5yr_unrepeated_source"]["events"]+a["recent_5yr_repeated_source"]["events"]==91
    assert a["recent_5yr_unrepeated_source"]["rows"]+a["recent_5yr_repeated_source"]["rows"]==865
    splits=x["prediction"]["holdouts"]
    assert len(splits)==4
    for s in splits:
        assert s["train_through"]<s["test"][0]<=s["test"][1]<=2022
        assert 0<s["test_events"]<s["test_rows"]
        for v in s["logloss"].values():
            assert math.isfinite(v) and v>=0
    main=splits[-1]
    assert (main["train_through"],main["test"],main["test_rows"],main["test_events"],main["test_islands"])==(2017,[2018,2022],372,28,14)
    assert main["logloss"]["ever_repeat"]>main["logloss"]["baseline"]
    assert main["logloss"]["recent_repeat"]<main["logloss"]["baseline"]
    b=x["island_bootstrap"]
    assert b["n"]==1500
    assert "no refit" in b["interval_type"]
    for key in ("ever_repeat","recent_repeat"):
        assert b[key]["ci95"][0]<0<b[key]["ci95"][1]
    assert x["guardrails"]["first_record_not_arrival"]
    assert x["guardrails"]["repeat_not_establishment"]
    assert x["guardrails"]["shipping_not_observed"]
