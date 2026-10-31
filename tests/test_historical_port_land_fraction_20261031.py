"""Protect land-area comparability of historical 38-port GBIF plant contrasts.

A 1-km radius at a coast may include much more sea at a port than at an
offport comparison. The land-cover fraction sensitivity is descriptive,
and matching sensitivity cutoffs were post-hoc, not a causal experiment.
"""
from __future__ import annotations

import csv
import json
import math
from collections import defaultdict
from pathlib import Path

D=Path(__file__).resolve().parents[1]/"provenance/reviewer_defenses/results/historical_port_nonport_botany_20261031"


def read(name):
    with (D/name).open(newline="",encoding="utf8") as f:
        return list(csv.DictReader(f))


def test_land_fraction_site_and_pair_support():
    r=json.loads((D/"land_fraction_quality_results.json").read_text(encoding="utf8"))
    sites=read("land_fraction_sites.csv")
    pairs=read("land_fraction_pairs.csv")
    assert r["status"]=="SOURCE_LOCKED_HISTORICAL_BOTANY_GEOMETRIC_AREA_QA_NOT_CAUSAL_EFFECT"
    assert (r["source_islands"],r["ports"],r["sites"])==(17,38,114)
    assert len(sites)==114 and len(pairs)==38
    assert len({s["site_id"] for s in sites})==114
    assert len({p["port_id"] for p in pairs})==38
    by_port=defaultdict(list)
    for s in sites:
        assert -1e-10<=float(s["land_area_fraction_1km"])<=1+1e-10
        by_port[s["port_id"]].append(s)
    assert all(len(v)==3 for v in by_port.values())
    for p in pairs:
        v=by_port[p["port_id"]]
        port=[x for x in v if x["site_type"]=="port"]
        nonport=[x for x in v if x["site_type"]=="offport"]
        assert len(port)==1 and len(nonport)==2
        vals=[float(x["land_area_fraction_1km"]) for x in nonport]
        portf=float(port[0]["land_area_fraction_1km"])
        assert math.isclose(portf,float(p["port_land_fraction"]),rel_tol=0,abs_tol=1e-10)
        assert math.isclose(sum(vals)/2,float(p["mean_nonport_land_fraction"]),rel_tol=0,abs_tol=1e-10)
        maxdiff=max(abs(portf-x) for x in vals)
        assert math.isclose(maxdiff,float(p["max_abs_port_control_land_fraction_diff"]),rel_tol=0,abs_tol=1e-10)
    assert r["land_fraction_mean_port"]<r["land_fraction_mean_nonport"]


def test_land_area_matched_reanalysis_fails_to_establish_port_effect():
    report=json.loads((D/"land_fraction_quality_results.json").read_text(encoding="utf8"))
    pairs=read("land_fraction_pairs.csv")
    assert len(report["tolerance_sensitivity"])==12
    for cutoff in (1,5,10,20):
        for tolerance in (.15,.25,.4):
            key=f"min{cutoff}_maxabsland{tolerance}"
            eligible=[p for p in pairs
                      if p[f"botanical_support_{cutoff}"]=="True"
                      and float(p["max_abs_port_control_land_fraction_diff"])<=tolerance]
            saved=report["tolerance_sensitivity"][key]
            assert len(eligible)==saved["matched_ports"]
            by_island=defaultdict(list)
            value_key=f"host_share_difference_{cutoff}"
            for p in eligible:
                by_island[p["island_id"]].append(float(p[value_key]))
            assert len(by_island)==saved["matched_islands"]
            island_equal=sum(sum(v)/len(v) for v in by_island.values())/len(by_island)
            assert math.isclose(island_equal,saved["island_equal_top10_port_minus_nonport_share"],
                                rel_tol=0,abs_tol=1e-11)
    strict=report["tolerance_sensitivity"]["min1_maxabsland0.15"]
    assert (strict["matched_ports"],strict["matched_islands"])==(4,4)
    assert strict["island_equal_top10_port_minus_nonport_share"]<0
    assert any("not year-specific" in s for s in report["notes"])
    assert any("causal" in s.lower() for s in report["notes"])
