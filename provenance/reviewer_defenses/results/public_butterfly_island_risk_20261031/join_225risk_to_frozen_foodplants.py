#!/usr/bin/env python3
"""Join independent, outcome-unselected butterfly × island GBIF index risk to
historically indexed documented foodplants in harbor and non-harbor coasts.

CONDITIONS: only 42 *published* genuine frozen HOSTS name edges (35 butterfly
species and 6 plant taxa), never proposed BCE-addition candidates or the 2,524
unpublished author edges. The species are selected because they appeared in a
previous exploratory source, so this is a sensitivity cohort, NOT confirmatory.

GBIF current index first record != first arrival, and a candidate foodplant
unindexed at a site cannot be coded as biological absence. Contemporary
accepted keys do not prove historical island-specific introduced status, or
whether the plant was moved by ship. No in situ larval feeding or breeding.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import random
from collections import defaultdict
from datetime import datetime,timezone
from pathlib import Path
from urllib.request import Request,urlopen
from urllib.parse import urlencode
from urllib.error import HTTPError,URLError
import time

ROOT=Path(__file__).resolve().parents[4]
HERE=Path(__file__).resolve().parent
SOURCES={
 "edges":("provenance/reviewer_defenses/results/top_shared_port_hosts_20261030/top_shared_host_link_provenance.csv",
         "22099e2198c33a683c3ba4e0c19a81b0e40a18fb"),
 "butterfly_index":("provenance/reviewer_defenses/results/public_butterfly_island_risk_20261031/full17_butterfly_island_225_atrisk.csv",None),
 "botany":("provenance/reviewer_defenses/results/historical_port_nonport_botany_20261031/historical_plant_site_facets.csv",
           "eb09a3f3ec053c00b43e8bc0678d6643ce02697b"),
 "migration":("provenance/reviewer_defenses/results/butterfly_maritime_trait_gap_20261031/migration_gbif_author2025_crosswalk.csv",
            "439dd1034ca16f96ff94222695d98dc39747c2e4"),
 "island_group":("provenance/reviewer_defenses/results/historical_botanical_17island_20261031/near_pair_pre2011_plant_composition.json",
                 "54c4515947d2a3ee7545e4007985ef82a2af6868")
}
GBIF="https://api.gbif.org/v1/species/match"
UA="chocho-2026-independently-published-FOODPLANT-35species-temporal-at-risk-join"
DEFAULT_ISLAND_MIN_SPECIES_INDEX=5

def load(name):
    path,sha=SOURCES[name]
    raw=(ROOT/path).read_bytes()
    actual=hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
    if sha is not None and actual!=sha:
        raise ValueError("Pinned Git source differs: "+path+" "+actual)
    if path.endswith(".json"):
        return json.loads(raw),actual
    return list(csv.DictReader(io.StringIO(raw.decode("utf8")))),actual

def get_species(name):
    p={"name":name,"rank":"SPECIES","kingdom":"Plantae","strict":"true"}
    url=GBIF+"?"+urlencode(p)
    last=None
    for attempt in range(4):
        try:
            with urlopen(Request(url,headers={"Accept":"application/json","User-Agent":UA}),timeout=40) as resp:
                data=json.load(resp)
            if str(data.get("canonicalName") or "").lower()!=name.lower():
                raise ValueError("GBIF canonical name differs: "+name+" => "+str(data)[:400])
            if data.get("matchType")!="EXACT" or data.get("rank")!="SPECIES" or int(data.get("confidence",0))<90:
                raise ValueError("Unverified GBIF species: "+name+" "+str(data)[:400])
            accepted=data.get("acceptedUsageKey") or data.get("usageKey")
            if not accepted:
                raise ValueError("Matched but no accepted species key")
            return {"name":name,"accepted_plant_speciesKey":str(accepted),
                "source_match_type":data.get("matchType"),"source_status":data.get("status"),
                "source_rank":data.get("rank"),"source_confidence":data.get("confidence"),
                "name_applies_to_historic_island_alien_status":False}
        except (HTTPError,URLError,TimeoutError,ValueError) as ex:
            last=ex
            if isinstance(ex,HTTPError) and ex.code not in (408,425,429,500,502,503,504):
                break
            time.sleep(min(5,.5*(2**attempt)))
    raise ValueError(f"GBIF strict source name resolution failed for {name}: {type(last).__name__} {str(last)[:250]}")

def hypergeom_island_exact_margin(rows,exposure,period):
    # Margins per island: preserve how many species were eligible, how many had
    # observed new GBIF indices, and how many had a named host indexed. Repeated
    # butterfly species across islands are NOT independent, so this is only a
    # stricter island sampling control, not an exact causal null.
    eligible="eligible_"+period
    outcome="first_indexed_in_"+period
    by=defaultdict(list)
    for r in rows:
        if int(r[eligible]):by[r["landmass_id"]].append(r)
    observed=0;expected=0;variance=0
    summary=[]
    for island in sorted(by,key=int):
        group=by[island]
        N=len(group)
        K=sum(int(r[exposure]) for r in group)
        n=sum(int(r[outcome]) for r in group)
        x=sum(int(r[exposure])*int(r[outcome]) for r in group)
        mu=n*K/N
        var=n*(K/N)*(1-K/N)*(N-n)/(N-1) if N>1 else 0
        observed+=x;expected+=mu;variance+=var
        summary.append({"island":island,"risk_cells":N,"host_indexed_risk":K,
            "first_indexes":n,"host_indexed_first_indexes":x,"null_expectation":mu})
    return {"observed_coindexed_new_records":observed,
            "conditional_island_margin_expected":expected,
            "island_margin_variance":variance,
            "island_standardized_score":(observed-expected)/math.sqrt(variance) if variance>0 else None,
            "island_contributions":summary,
            "method":"hypergeometric exact island first-index and host-index margins, no species or region constraint",
            "uncertainty_is_not_corrected_for_repeated_species_or_recording_bias":True}

def analyze():
    edges,source_edges=load("edges")
    index,source_index=load("butterfly_index")
    plants,source_plants=load("botany")
    migration,source_migration=load("migration")
    islands,source_group=load("island_group")
    if len(edges)!=49 or len(index)!=3825 or len(plants)!=114 or len(migration)!=225:
        raise ValueError("Published source coverage has changed")
    genuine=[x for x in edges if x["in_frozen_HOSTS"]=="True"]
    assert len(genuine)==42
    hostnames=sorted({x["host_name"] for x in genuine})
    assert len(hostnames)==6
    taxonomy=[get_species(name) for name in hostnames]
    plant_keys={r["name"]:r["accepted_plant_speciesKey"] for r in taxonomy}
    if plant_keys["Medicago sativa"]!="9151957":
        raise ValueError("Known previously verified Medicago key differs")
    by_species=defaultdict(set)
    for edge in genuine:by_species[edge["species"]].add(edge["host_name"])
    assert len(by_species)==35

    rownames={x["species"] for x in index}
    if not set(by_species)<=rownames:
        raise ValueError("Not all published verified host-edge butterflies are in the original 225 candidates")
    by_island=defaultdict(lambda:{"port":[],"offport":[]})
    for plant in plants:
        assert plant["query_status"]=="COMPLETE"
        by_island[plant["landmass_id"]][plant["site_type"]].append(plant)
    assert len(by_island)==17
    migration_by={r["species"]:r["author_2025_GBIF_requested_taxon_key_proxy"]
                  for r in migration}
    joined=[]
    for r in index:
        if r["species"] not in by_species:continue
        island=r["landmass_id"]
        known=by_species[r["species"]]
        keys={plant_keys[x] for x in known}
        pos=by_island[island]["port"]
        off=by_island[island]["offport"]
        assert len(pos)>0 and len(off)==2*len(pos)
        def indexed(site):
            return set(filter(None,site["species_keys"].split(";")))
        known_port=set().union(*(indexed(x)&keys for x in pos))
        known_off=set().union(*(indexed(x)&keys for x in off))
        port_plant=set().union(*(indexed(x) for x in pos))
        off_plant=set().union(*(indexed(x) for x in off))
        row=dict(r)
        row.update({
            "known_confirmed_foodplant_edges_of_6_total_names":len(known),
            "known_foodplant_source_names":";".join(sorted(known)),
            "known_foodplant_GBIF_keys":";".join(sorted(keys,key=int)),
            "pre2011_port_compatible_speciesKey_recorded":int(bool(known_port)),
            "pre2011_offport_compatible_speciesKey_recorded":int(bool(known_off)),
            "pre2011_any_compatible_speciesKey_recorded":int(bool(known_port|known_off)),
            "pre2011_port_known_foodplants_indexed_n":len(known_port),
            "pre2011_offport_known_foodplants_indexed_n":len(known_off),
            "port_known_foodplant_keys_recorded":";".join(sorted(known_port,key=int)),
            "offport_known_foodplant_keys_recorded":";".join(sorted(known_off,key=int)),
            "port_site_count":len(pos),"offport_site_count":len(off),
            "pre2011_port_all_plant_species_keys_n":len(port_plant),
            "pre2011_offport_all_plant_species_keys_n":len(off_plant),
            "pre2011_port_sampled_sites_with_plant_index":sum(bool(indexed(x)) for x in pos),
            "pre2011_offport_sampled_sites_with_plant_index":sum(bool(indexed(x)) for x in off),
            "pre2011_both_area_botanical_index_five_species_each":int(len(port_plant)>=5 and len(off_plant)>=5),
            "migration_selected_positive_only":int(migration_by.get(r["species"])==
                                              "REVIEW_SELECTED_NAME_MATCHED_GBIF_REQUEST"),
            "not_selected_35_on_future_first_index":True,
            "unrecorded_named_foodplant_not_biological_absence":True,
        })
        joined.append(row)
    assert len(joined)==35*17
    assert len({(x["speciesKey"],x["landmass_id"]) for x in joined})==len(joined)

    result={
     "schema":"chocho_fixed225_unselected_risk_host35_6plants_2001_2010_island_joint_v1",
     "status":"EXPLORATORY_NAME_SPECIFIC_HISTORICAL_HOST_BUTTERFLY_INDEX_JOIN_NOT_CAUSAL",
     "checked_utc":datetime.now(timezone.utc).isoformat(),
     "source_git_blobs":{"edges":source_edges,"original225risk":source_index,
       "pre2011_botany":source_plants,"migration":source_migration,
       "original17_islands":source_group},
     "taxonomy":taxonomy,"host_species_butterflies":35,
     "original_225_predeclared_species":225,
     "source_TRUE_Host_edges":42,"foodplant_species_names":6,
     "candidate_species_island_cells":len(joined),"islands":17,
     "source_direct_HOSTS_butterflies_are_historically_posthoc_selected":True,
     "periods":{},
     "sample_quality":{}
    }
    for period in ("2011_2017","2018_2022"):
        e="eligible_"+period;o="first_indexed_in_"+period
        risk=[x for x in joined if int(x[e])]
        cases=sum(int(x[o]) for x in risk)
        result["periods"][period]={
            "risk_cells":len(risk),"new_indexed_butterfly_island_pairs":cases,
            "by_early_host_presence":{}
        }
        for exposure in ("pre2011_port_compatible_speciesKey_recorded",
                         "pre2011_offport_compatible_speciesKey_recorded",
                         "pre2011_any_compatible_speciesKey_recorded"):
            block={}
            for label in (0,1):
                subset=[x for x in risk if int(x[exposure])==label]
                k=sum(int(x[o]) for x in subset)
                block[str(label)]={"risk_cells":len(subset),"first_indexed":k,
                       "indexed_first_rate":k/len(subset) if subset else None}
            result["periods"][period]["by_early_host_presence"][exposure]=block
            result["periods"][period]["island_margin_null_"+exposure]=(
                hypergeom_island_exact_margin(joined,exposure,period))
        for quality in (0,1):
            q=[x for x in risk if int(x["pre2011_both_area_botanical_index_five_species_each"])==quality]
            result["periods"][period]["botanical_5key_quality_"+str(quality)]={
                "risk_cells":len(q),"new_indexes":sum(int(x[o]) for x in q),
                "islands":len({x["landmass_id"] for x in q}),
                "host_port_detected_cells":sum(int(x["pre2011_port_compatible_speciesKey_recorded"]) for x in q)
            }
    # These sources already leave most per-island first records unclassified by
    # full HOSTS, island-specific historical native ranges and GBIF precision.
    result["limits"]=[
        "Only 35 species already selected in exploratory shared-host work, 42 frozen genuine HOSTS relations, not complete host portfolios.",
        "All 225 butterfly candidates were frozen before post-2010 GBIF island records; 35-species host subset was prior exploratory selection, not independent validation.",
        "Missing GBIF foodplant record is not botanical absence or no larval host, especially for partly known host portfolios.",
        "2001-2010 GBIF host records do not establish historical introduced status on particular islands or travel by ship.",
        "2019 WPI harbor points and later coastline are spatial proxies rather than verified historic harbor activity.",
        "17 islands were themselves originally selected for near-neighbor retrospective island-study design, not random global sample.",
        "GBIF current 2026 butterfly index is different from the original precision-constrained private butterfly-year Parquet.",
        "All-species candidate nonindexes include many regionally implausible butterfly-island units.",
        "The island-only fixed-margin null does not preserve butterfly species or region margins; do not assign causal p-values.",
        "Island and offport GBIF plant recording density is unequal; true breeding and larval feeding not measured.",
        "Independent documented migration is positive-only; not-listed is migration status UNKNOWN, not resident.",
        "No 225-species, 17-island projection of actual shipping or island stepping-stone transition."
    ]
    return joined,result

def write_csv(file,rows):
    with file.open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader();w.writerows(rows)

def main():
    rows,obj=analyze()
    HERE.mkdir(parents=True,exist_ok=True)
    write_csv(HERE/"public_35host_butterfly_x_17island_risk.csv",rows)
    (HERE/"public_35host_named_plant_join.json").write_text(
        json.dumps(obj,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print(json.dumps({"status":obj["status"],"taxonomy":obj["taxonomy"],
      "cells":obj["candidate_species_island_cells"],"periods":obj["periods"]},
      ensure_ascii=False,indent=2))

if __name__=="__main__":main()
