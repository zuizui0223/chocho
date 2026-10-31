"""Guard 1-km source-locked historical harbor/nonharbor land-fraction confounding.

The comparison is retrospective and can never establish that ships delivered
the candidate foodplants or that butterflies bred on these islands.
"""
from __future__ import annotations

import csv
import json
import math
from collections import defaultdict
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/"provenance/reviewer_defenses/results/historical_port_nonport_botany_20261031"

def csv_rows(name):
    with (DATA/name).open(encoding="utf8",newline="") as f:
        return list(csv.DictReader(f))

def mean(v):
    return sum(v)/len(v) if v else None

def test_original_historical_land_area_sensitivity_recomputed():
    files=["land_fraction_sites.csv","land_fraction_pairs.csv","land_fraction_quality_results.json"]
    if not all((DATA/x).is_file() for x in files):
        pytest.skip("Independent GEB anonymous review bundle omits exploratory historical island and harbor sources")
    sites=csv_rows(files[0])
    pairs=csv_rows(files[1])
    audit=json.loads((DATA/files[2]).read_text(encoding="utf8"))
    assert audit["status"]=="SOURCE_LOCKED_HISTORICAL_BOTANY_GEOMETRIC_AREA_QA_NOT_CAUSAL_EFFECT"
    assert (len(sites),len(pairs),audit["source_islands"])==(114,38,17)
    assert len({r["site_id"] for r in sites})==114
    assert len({r["landmass_id"] for r in sites})==17
    for r in sites:
        val=float(r["land_area_fraction_1km"])
        assert 0<=val<=1
    for kind,expected in [("port",audit["land_fraction_mean_port"]),
                          ("offport",audit["land_fraction_mean_nonport"])]:
        subset=[float(r["land_area_fraction_1km"]) for r in sites if r["site_type"]==kind]
        assert math.isclose(mean(subset),expected,rel_tol=0,abs_tol=1e-12)
    assert audit["land_fraction_mean_port"] < audit["land_fraction_mean_nonport"]
    med=sorted(float(r["max_abs_port_control_land_fraction_diff"]) for r in pairs)[len(pairs)//2]
    assert math.isclose(med,audit["port_control_maxabs_landfrac_diff_median"],rel_tol=0,abs_tol=1e-12)
    for threshold in (1,5,10,20):
        for tol in (.15,.25,.4):
            key=f"min{threshold}_maxabsland{tol:g}"
            reference=audit["tolerance_sensitivity"][key]
            eligible=[r for r in pairs if r[f"botanical_support_{threshold}"]=="True"
               and float(r["max_abs_port_control_land_fraction_diff"])<=tol]
            by_island=defaultdict(list)
            for r in eligible:
                by_island[r["island_id"]].append(float(r[f"host_share_difference_{threshold}"]))
            value=mean([mean(v) for v in by_island.values()])
            assert reference["matched_ports"]==len(eligible)
            assert reference["matched_islands"]==len(by_island)
            if value is None:
                assert reference["island_equal_top10_port_minus_nonport_share"] is None
            else:
                assert math.isclose(value,reference["island_equal_top10_port_minus_nonport_share"],
                                    rel_tol=0,abs_tol=1e-12)
            assert reference["positive_port_count"]==sum(float(x[f"host_share_difference_{threshold}"])>0 for x in eligible)
            assert reference["negative_port_count"]==sum(float(x[f"host_share_difference_{threshold}"])<0 for x in eligible)
    assert audit["tolerance_sensitivity"]["min10_maxabsland0.25"]["matched_islands"]==3
    assert audit["tolerance_sensitivity"]["min10_maxabsland0.25"]["island_equal_top10_port_minus_nonport_share"]<0

def test_historical_plant_source_has_no_butterfly_colonization_claim():
    target=DATA/"README.md"
    if not target.exists():
        pytest.skip("Separate port-first paper results not part of the standalone GEB blinded bundle")
    t=target.read_text(encoding="utf8")
    assert "61.4%" in t
    assert "88.5%" in t
    assert "3 / 3" in t or "3港・3島" in t
    assert "cannot assume" in t or "できない" in t or "not assume" in t
