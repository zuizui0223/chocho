#!/usr/bin/env python3
"""Independently verify published full38 GBIF port/nonport botanical aggregates.

No GBIF network or GSHHG download required. Replays every ratio, source
speciesKey membership, site/port geometry provenance, and island-equal means.
It does NOT verify the original GBIF point records or actual plant introductions.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import random
from collections import Counter,defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[4]
HERE=Path(__file__).resolve().parent

SOURCES={
    "sites":HERE/"historical_plant_site_facets.csv",
    "design":HERE/"fixed_geometric_sites.csv",
    "admission":HERE/"geometry_admissions.csv",
    "pairs":HERE/"port_offport_historical_pairs.csv",
    "receipt":HERE/"historical_port_offport_results.json",
    "top10":ROOT/"provenance/reviewer_defenses/results/geb_port_island_realization_20261031/GEB_top10_pre2011_port_plant_overlap.json",
}

def csv_rows(path):
    with path.open(encoding="utf8",newline="") as f:
        return list(csv.DictReader(f))

def close(x,y,tol=1e-10):
    if x is None or y is None:
        return x is None and y is None
    return math.isclose(float(x),float(y),rel_tol=0,abs_tol=tol)

def fmt_float(cell):
    if cell is None or cell=="":return None
    return float(cell)

def mean(vals):
    return sum(vals)/len(vals) if vals else None

def quantile(v,p):
    if not v:return None
    v=sorted(v)
    i=(len(v)-1)*p
    lo=math.floor(i);hi=math.ceil(i)
    return v[lo]+(v[hi]-v[lo])*(i-lo)

def source_blob(path):
    data=path.read_bytes()
    return hashlib.sha1(b"blob "+str(len(data)).encode()+b"\0"+data).hexdigest()

def compute(draws=4000):
    src=SOURCES
    receipt=json.loads(src["receipt"].read_text(encoding="utf8"))
    known=json.loads(src["top10"].read_text(encoding="utf8"))
    sites=csv_rows(src["sites"])
    designs=csv_rows(src["design"])
    pairs=csv_rows(src["pairs"])
    decisions=csv_rows(src["admission"])
    top={str(x["GBIF_accepted_speciesKey"]) for x in known["source_contributor_plants"]}
    assert len(top)==10
    assert source_blob(src["top10"])==receipt["source_git_blobs"]["top_hosts"]["sha"]
    assert receipt["gshhg_zip_sha256"]=="8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41"
    assert receipt["window"]=="2001-2010 historical BEFORE 2011-2022 butterfly index"
    assert receipt["port_buffer_km"]==receipt["offport_buffer_km"]==1
    assert len(sites)==len(designs)==114
    assert len(pairs)==len(decisions)==38
    assert len({p["landmass_id"] for p in sites})==17
    assert len({p["port_id"] for p in sites})==38
    assert len({p["site_id"] for p in sites})==114
    assert set(p["site_id"] for p in sites)==set(p["site_id"] for p in designs)
    assert all(p["query_status"]=="COMPLETE" for p in sites)
    assert all(q["status"]=="MATCHED_GEOMETRIC_DESIGN_ONLY" and q["n_controls"]=="2"
               for q in decisions)
    by_site={x["site_id"]:x for x in sites}
    grouped=defaultdict(list)
    for row in sites:
        grouped[row["port_id"]].append(row)
        spp=set(filter(None,row["species_keys"].split(";")))
        matches=set(filter(None,row["target_host_species_keys"].split(";")))
        assert len(spp)==int(row["plant_indexed_species_keys"])
        assert matches==spp&top
        assert len(matches)==int(row["n_target_hosts"])
        assert int(row["plant_occurrence_index_count"])>=len(spp)
        assert int(row["provider_datasets_minimum"])>=0
    for d in designs:
        assert d["site_id"] in by_site
        s=by_site[d["site_id"]]
        assert s["latitude"]==d["latitude"] and s["longitude"]==d["longitude"]
        if d["site_type"]=="offport":
            assert float(d["nearest_wpi_km"])>=5-1e-5
            assert .15-1e-6<=float(d["shoreline_km"])<=1.25+1e-6
            assert float(d["port_distance_km"])>=5
    for port,group in grouped.items():
        assert len(group)==3
        assert sum(x["site_type"]=="port" for x in group)==1
        assert sum(x["site_type"]=="offport" for x in group)==2

    pp={str(x["port_id"]):x for x in pairs}
    for port,group in grouped.items():
        row=pp[port]
        p=next(x for x in group if x["site_type"]=="port")
        c=[x for x in group if x["site_type"]=="offport"]
        assert row["status"]=="PAIRED_INDEXED"
        assert int(row["n_controls"])==2
        assert int(row["port_plant_species_keys"])==int(p["plant_indexed_species_keys"])
        assert int(row["port_top10_count"])==int(p["n_target_hosts"])
        assert close(row["offport_plant_species_keys_mean"],
                      mean([int(x["plant_indexed_species_keys"]) for x in c]))
        assert close(row["offport_top10_count_mean"],
                      mean([int(x["n_target_hosts"]) for x in c]))
        for threshold in (1,5,10,20,50):
            ok=int(p["plant_indexed_species_keys"])>=threshold and all(
                int(x["plant_indexed_species_keys"])>=threshold for x in c)
            assert (row["meets_min_species_"+str(threshold)]=="True")==ok
            observed=fmt_float(row["port_minus_control_top10_index_fraction_"+str(threshold)])
            expected=(int(p["n_target_hosts"])/int(p["plant_indexed_species_keys"])
                -mean([int(x["n_target_hosts"])/int(x["plant_indexed_species_keys"])
                       for x in c])) if ok else None
            assert close(observed,expected)

    scan={}
    for threshold in (1,5,10,20,50):
        fld="port_minus_control_top10_index_fraction_"+str(threshold)
        eligible=[p for p in pairs if p["meets_min_species_"+str(threshold)]=="True"]
        island_map=defaultdict(list)
        for p in eligible:island_map[p["island_id"]].append(float(p[fld]))
        island_means=[mean(v) for v in island_map.values()]
        effect=mean(island_means)
        r=receipt["results_by_min_background_plant_species_keys"][str(threshold)]
        assert r["eligible_port_pairs"]==len(eligible)
        assert r["eligible_islands"]==len(island_map)
        assert close(r["island_equal_mean_port_minus_offport_host_fraction"],effect)
        assert r["positive_port_pairs"]==sum(float(x[fld])>0 for x in eligible)
        assert r["negative_port_pairs"]==sum(float(x[fld])<0 for x in eligible)
        assert r["equal_port_pairs"]==sum(float(x[fld])==0 for x in eligible)
        rng=random.Random(20261031)
        boot=[mean([rng.choice(island_means) for _ in island_means])
              for draw in range(draws)] if island_means else []
        scan[str(threshold)]={
            "ports":len(eligible),"islands":len(island_map),
            "island_equal_mean":effect,
            "island_bootstrap_95":[quantile(boot,.025),quantile(boot,.975)] if boot else None,
            "positive_ports":r["positive_port_pairs"],
            "negative_ports":r["negative_port_pairs"],
            "equal_ports":r["equal_port_pairs"],
        }
    coverage={}
    for typ in ("port","offport"):
        rows=[x for x in sites if x["site_type"]==typ]
        coverage[typ]={
            "sites":len(rows),
            "sites_without_2001_2010_botanical_records":sum(int(x["plant_occurrence_index_count"])==0 for x in rows),
            "sites_without_resolved_specieskey":sum(int(x["plant_indexed_species_keys"])==0 for x in rows),
            "sites_indexing_any_global_top10_host_candidate":sum(int(x["n_target_hosts"])>0 for x in rows),
            "mean_botanical_specieskey_count":mean([int(x["plant_indexed_species_keys"]) for x in rows]),
            "mean_botanical_occurrence_index":mean([int(x["plant_occurrence_index_count"]) for x in rows]),
            "mean_distinct_publishing_datasets_minimum":mean([int(x["provider_datasets_minimum"]) for x in rows])
        }
    return {
        "schema":"chocho_full38_historical_botany_independent_source_replay_v1",
        "status":"PASS_PUBLISHED_AGGREGATE_REPLAY_NOT_CAUSAL_OR_RAW_GBIF_VALIDATION",
        "exact_17_islands_38_ports_76_offport_sites":True,
        "all_historical_sites_complete":True,
        "observed_2001_2010_biology_only":True,
        "source_blob_hashes":{k:source_blob(p) for k,p in src.items()},
        "coverage":coverage,
        "threshold_sensitivity":scan,
        "bootstrap_draws":draws,
        "limitations":[
            "Bootstrap resamples at most 10 independent islands and cannot measure historical transport or reporting error.",
            "Globally important top10 host plants are not confirmed alien or larval hosts on the target island at the target time.",
            "Even exact WPI IDs are 2019 retrospectively mapped geographic proxies.",
            "Off-port coastal points exclude catalogued modern WPI within 5km but can include unlisted human infrastructure.",
            "Zero GBIF indexed botanical speciesKey means missing observed records, never biological plant absence.",
            "Unequal terrestrial area fraction of 1km harbor versus inland controls has not been matched.",
            "No butterfly observations, larval feeding, breeding, or island relay measured in this H1 source audit."
        ]
    }

def main():
    output=compute()
    target=HERE/"full38_replay_and_bias_audit.json"
    target.write_text(json.dumps(output,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print(json.dumps({"status":output["status"],"coverage":output["coverage"],
                      "thresholds":output["threshold_sensitivity"]},
        ensure_ascii=False,indent=2))

if __name__=="__main__":main()
