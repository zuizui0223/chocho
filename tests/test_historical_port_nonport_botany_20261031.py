"""Guard historical matched port-vs-nonport botany against false zero-species facets.

The 2001–2010 source is independent of later butterfly indexing but contains
neither verified historical island-specific alien status nor known ship routes.
"""
from __future__ import annotations

from collections import defaultdict
import csv
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/"provenance/reviewer_defenses/results/historical_port_nonport_botany_20261031"


def rows(name):
    with (DATA/name).open(encoding="utf8",newline="") as stream:
        return list(csv.DictReader(stream))


def number(row,key):
    x=row[key]
    return float(x) if x!="" else None


def test_faceted_plant_records_not_silently_zeroed():
    s=rows("historical_plant_site_facets.csv")
    r=json.loads((DATA/"historical_port_offport_results.json").read_text(encoding="utf8"))
    assert len(s)==r["historical_plant_query_complete"]+r["historical_plant_query_failed"]
    assert r["sites_geometrically_matched"]==len(s)
    assert len({x["site_id"] for x in s})==len(s)
    assert r["port_buffer_km"]==r["offport_buffer_km"]==1
    assert r["window"].startswith("2001-2010")
    assert r["source_ports_in_design"]>=17
    assert r["source_ports_in_design"]<=r["original_port_candidates"]==38
    successful=[x for x in s if x["query_status"]=="COMPLETE"]
    assert successful
    recorded=[x for x in successful if float(x["plant_occurrence_index_count"])>0]
    resolved=[x for x in recorded if float(x["plant_indexed_species_keys"])>0]
    assert len(recorded)>0 and len(resolved)>0, "Facet parsing failure: botanical records but all zero species keys"
    assert sum(x["site_type"]=="port" for x in s)>0
    assert sum(x["site_type"]=="offport" for x in s)>0
    for x in successful:
        n=float(x["plant_indexed_species_keys"])
        host=float(x["n_target_hosts"])
        assert 0<=host<=n
        assert n>=0
        if n>0:
            assert len([z for z in x["species_keys"].split(";") if z])==int(n)
        if host>0:
            assert len([z for z in x["target_host_species_keys"].split(";") if z])==int(host)
    code=(DATA/"acquire_matched_botany.py").read_text(encoding="utf8")
    assert 'facets.get("SPECIES_KEY",[])' in code
    assert "no SPECIES_KEY facet" in code


def test_port_nonport_contrasts_recomputed_from_site_level_indices():
    source=rows("historical_plant_site_facets.csv")
    pairs=rows("port_offport_historical_pairs.csv")
    report=json.loads((DATA/"historical_port_offport_results.json").read_text(encoding="utf8"))
    by_port=defaultdict(list)
    for r in source:
        by_port[r["port_id"]].append(r)
    assert len(pairs)==len(by_port)==report["source_ports_in_design"]
    for p in pairs:
        grp=by_port[p["port_id"]]
        port=[x for x in grp if x["site_type"]=="port"]
        off=[x for x in grp if x["site_type"]=="offport"]
        assert len(port)==1 and len(off)>=2
        port=port[0]
        assert p["island_id"]==port["landmass_id"]
        ready=all(x["query_status"]=="COMPLETE" for x in grp)
        assert (p["status"]=="PAIRED_INDEXED")==ready
        for th in (1,5,10,20,50):
            key=f"port_minus_control_top10_index_fraction_{th}"
            eligible=ready and all(float(x["plant_indexed_species_keys"])>=th for x in grp)
            assert (p[f"meets_min_species_{th}"]=="True")==eligible
            if not eligible:
                assert p[key]==""
                continue
            fx=lambda x:float(x["n_target_hosts"])/float(x["plant_indexed_species_keys"])
            expected=fx(port)-sum(fx(x) for x in off)/len(off)
            assert math.isclose(float(p[key]),expected,rel_tol=0,abs_tol=1e-12)
    for th in (1,5,10,20,50):
        eligible=[p for p in pairs if p[f"meets_min_species_{th}"]=="True"]
        stats=report["results_by_min_background_plant_species_keys"][str(th)]
        assert stats["eligible_port_pairs"]==len(eligible)
        assert stats["eligible_islands"]==len(set(p["island_id"] for p in eligible))
        if eligible:
            island=defaultdict(list)
            for p in eligible:
                island[p["island_id"]].append(float(p[f"port_minus_control_top10_index_fraction_{th}"]))
            mean=sum(sum(v)/len(v) for v in island.values())/len(island)
            assert math.isclose(stats["island_equal_mean_port_minus_offport_host_fraction"],mean,rel_tol=0,abs_tol=1e-12)
            assert sum(float(p[f"port_minus_control_top10_index_fraction_{th}"])>0 for p in eligible)==stats["positive_port_pairs"]
            assert sum(float(p[f"port_minus_control_top10_index_fraction_{th}"])<0 for p in eligible)==stats["negative_port_pairs"]
        else:
            assert stats["island_equal_mean_port_minus_offport_host_fraction"] is None
    assert any("not verified introduced" in item or "NOT verified introduced" in item for item in report["limits"])
