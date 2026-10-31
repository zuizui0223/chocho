#!/usr/bin/env python3
"""Historical host-index vs *butterfly* first-index: source-audited public panel.

Target: frozen candidate 225 species × original 17 GSHHG islands, with 2001-2010
GBIF observations prior to 2011-2017 and 2018-2022 butterfly index windows.
Only the 42 explicitly verified in-frozen-HOSTS foodplant edges (35 species,
6 plant names) may contribute to the primary host-compatibility diagnostic.

Do NOT interpret a non-indexed plant as absent; matched rows are at-risk GBIF
index records, not documented colonizations. A 2019 WPI location does not prove
a ship path or 2001-era harbor. The source host edges were selected for earlier
port exploration, not a representative randomized sample of 225 butterflies.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import random
import time
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
from urllib.error import HTTPError,URLError
from urllib.parse import urlencode
from urllib.request import Request,urlopen

ROOT=Path(__file__).resolve().parents[4]
HERE=Path(__file__).resolve().parent
SOURCES={
 "risk":("provenance/reviewer_defenses/results/public_butterfly_island_risk_20261031/full17_butterfly_island_225_atrisk.csv",
         "60e38e31d1d1b8896fb925aaa5ee0a405fbef2df"),
 "risk_receipt":("provenance/reviewer_defenses/results/public_butterfly_island_risk_20261031/full17_query_receipt.json",
                 "6ef70869578e1bde3bca13bb49e44ef3b4fb4360"),
 "host_edges":("provenance/reviewer_defenses/results/top_shared_port_hosts_20261030/top_shared_host_link_provenance.csv",
               "22099e2198c33a683c3ba4e0c19a81b0e40a18fb"),
 "flora":("provenance/reviewer_defenses/results/historical_port_nonport_botany_20261031/historical_plant_site_facets.csv",
          "eb09a3f3ec053c00b43e8bc0678d6643ce02697b"),
 "plant_receipt":("provenance/reviewer_defenses/results/historical_port_nonport_botany_20261031/historical_port_offport_results.json",
                  "53ef6dddbe1a8efceaf518812bf2b3c77cf285d6"),
 "land_frac":("provenance/reviewer_defenses/results/historical_port_nonport_botany_20261031/land_fraction_sites.csv",
              "b2f7a7597550b61abbc182b6c5c6f2f22ae6a898"),
 "mobility":("provenance/reviewer_defenses/results/butterfly_maritime_trait_gap_20261031/migration_gbif_author2025_crosswalk.csv",
             "439dd1034ca16f96ff94222695d98dc39747c2e4"),
 "regional_null":("provenance/reviewer_defenses/results/geb_port_island_realization_20261031/GEB_port_island_realization_risk_20261031.json",
                  "25103acdd5cd0e8e28624ec92c340544f7a622ae")
}
GBIF="https://api.gbif.org/v1/species/match"
UA="chocho-butterfly-first-index-verified-HOSTS-public-join/20261031"
MIN_PLANT_SPECIES=(1,5,10,20)
SEED=20261031

def read(k):
    path,expected=SOURCES[k]
    raw=(ROOT/path).read_bytes()
    sha=hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
    if sha!=expected:raise ValueError(f"Frozen source changed {path}: {sha} expected {expected}")
    if path.endswith(".json"):return json.loads(raw)
    return list(csv.DictReader(io.StringIO(raw.decode("utf8"))))

def resolve_name(name):
    url=GBIF+"?"+urlencode({"name":name,"rank":"SPECIES","kingdom":"Plantae","strict":"true"})
    last=None
    for trial in range(4):
        try:
            with urlopen(Request(url,headers={"Accept":"application/json","User-Agent":UA}),timeout=40) as f:
                m=json.load(f)
            if not isinstance(m,dict) or m.get("error"):raise ValueError("Invalid GBIF match data")
            accepted=m.get("acceptedUsageKey") or m.get("usageKey")
            if (m.get("matchType")!="EXACT" or str(m.get("rank","")).upper()!="SPECIES"
                or str(m.get("canonicalName","")).lower()!=name.lower()
                or not accepted or int(m.get("confidence",0))<90):
                raise ValueError("Cannot accept strict name matching "+name+" "+repr(m)[:200])
            return {"name":name,"accepted_speciesKey":str(accepted),
                    "matchType":m.get("matchType"),"rank":m.get("rank"),
                    "confidence":m.get("confidence"),"index_status":m.get("status"),
                    "original_usageKey":str(m.get("usageKey")),
                    "accepted_usageKey":str(accepted),
                    "GBIF_name_match_source":"GBIF species/match Plantae SPECIES exact strict"}
        except (HTTPError,URLError,TimeoutError,ValueError,json.JSONDecodeError) as e:
            last=e
            if isinstance(e,HTTPError) and e.code not in (408,425,429,500,502,503,504):
                break
            time.sleep(min(2,0.35*2**trial))
    raise RuntimeError("Cannot resolve independent host plant name "+name+" "+repr(last))

def cache_taxonomy(host_names):
    path=HERE/"verified_HOSTS_six_plant_taxonomy.json"
    if path.exists():
        receipt=json.loads(path.read_text(encoding="utf8"))
        if (receipt.get("host_plant_names")==sorted(host_names)
            and receipt.get("status")=="ALL_SIX_GBIF_ACCEPTED_SPECIESKEY_MATCHED"):
            return receipt,False
    names=[resolve_name(x) for x in sorted(host_names)]
    if len(names)!=6 or len(set(x["accepted_speciesKey"] for x in names))!=6:
        raise ValueError("Six source plant species do not map to six independent accepted GBIF names")
    receipt={
       "schema":"chocho_verified_HOSTS_six_plant_accepted_GBIF_taxon_keys_v1",
       "status":"ALL_SIX_GBIF_ACCEPTED_SPECIESKEY_MATCHED",
       "source_host_edges_git_blob":SOURCES["host_edges"][1],
       "host_plant_names":sorted(host_names),
       "source_match_utc":datetime.now(timezone.utc).isoformat(),
       "species":names,
       "historically_introduced_status_by_island_NOT_KNOWN":True,
       "local_larval_feeding_NOT_ESTABLISHED":True
    }
    path.write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    return receipt,True

def geographic_components(data):
    z=data.get("geographic_components")
    # Prior GEB risk receipt contains original six-region adjacency,
    # but shape can differ; derive region explicitly from 2011 risk report.
    periods=data["periods"]
    rows=periods["2011_2017_first_index_among_no_pre2011_butterfly_index"]["recording_risk"]["by_geographic_component"]
    groups=[sorted([str(i) for i in entry["islands"]],key=int) for entry in rows]
    if len(groups)!=6 or sorted(map(len,groups))!=[2,2,2,2,4,5]:
        raise ValueError("Source island geographical components changed")
    flat=[i for group in groups for i in group]
    if len(flat)!=17 or len(set(flat))!=17:raise ValueError("Missing original source island IDs")
    return {island:k for k,group in enumerate(groups) for island in group},groups

def join_rows():
    panel=read("risk");report=read("risk_receipt");edges=read("host_edges")
    flora=read("flora");plants=read("plant_receipt");land=read("land_frac")
    mobility=read("mobility")
    regions,_=geographic_components(read("regional_null"))
    if len(panel)!=3825 or len(edges)!=49 or len(flora)!=114 or len(land)!=114 or len(mobility)!=225:
        raise ValueError("Source input panel rows drift")
    if report["risk_cells"]!=3825 or plants["historical_plant_query_complete"]!=114:
        raise ValueError("Original frozen source query receipts contradict row number")
    accepted=[e for e in edges if e["in_frozen_HOSTS"]=="True" and e["BCE_added_candidate"]=="False"]
    if len(accepted)!=42 or len({e["species"] for e in accepted})!=35:
        raise ValueError("Expected 42 strictly frozen HOSTS verified edges for 35 butterfly species")
    names=set(e["host_name"] for e in accepted)
    if len(names)!=6:raise ValueError("Expected exactly six independently evidenced plant taxa")
    taxonomy,created=cache_taxonomy(names)
    key_by_name={e["name"]:e["accepted_speciesKey"] for e in taxonomy["species"]}
    if len(key_by_name)!=6:raise ValueError("Plant accepted taxonomy not unique")
    by_butterfly=defaultdict(set)
    for edge in accepted:
        by_butterfly[edge["species"]].add(key_by_name[edge["host_name"]])
    taxa_panel={x["species"] for x in panel}
    if not set(by_butterfly)<=taxa_panel:raise ValueError("Original verified butterfly outside 225 frozen candidate panel")
    if len({(x["species"],x["landmass_id"]) for x in panel})!=3825:
        raise ValueError("Duplicate species×island source risk rows")
    mig={x["species"]:x["author_2025_GBIF_requested_taxon_key_proxy"]
         for x in mobility}

    by_site=defaultdict(list)
    lf={x["site_id"]:x for x in land}
    for x in flora:
        if x["query_status"]!="COMPLETE":raise ValueError("All source plant queries must be complete")
        if x["site_id"] not in lf:raise ValueError("Missing source land fraction")
        if x["site_type"] not in ("port","offport"):raise ValueError("Unknown site type")
        by_site[x["landmass_id"]].append({
            "site_id":x["site_id"],"type":x["site_type"],
            "keys":set(filter(None,x["species_keys"].split(";"))),
            "n_flora_species":int(x["plant_indexed_species_keys"]),
            "flora_records":int(x["plant_occurrence_index_count"]),
            "land_fraction":float(lf[x["site_id"]]["land_area_fraction_1km"])
        })
    if len(by_site)!=17 or set(by_site)!=set(regions):
        raise ValueError("Source island botanical coverage does not match region mapping")
    assert sum(len(q) for q in by_site.values())==114
    assert all(0<=s["land_fraction"]<=1 for v in by_site.values() for s in v)

    rows=[]
    for x in panel:
        if x["species"] not in by_butterfly:continue
        sp=x["species"];iid=x["landmass_id"]
        compatible=by_butterfly[sp]
        src=by_site[iid]
        port=[z for z in src if z["type"]=="port"]
        off=[z for z in src if z["type"]=="offport"]
        if not port or len(off)!=2*len(port):raise ValueError("Missing matched harbor and two offport sites")
        hit_p=[z for z in port if compatible&z["keys"]]
        hit_o=[z for z in off if compatible&z["keys"]]
        row={
            "species":sp,"speciesKey":x["speciesKey"],"landmass_id":iid,
            "component":regions[iid],
            "is_migrant_roster_documented":int(mig.get(sp)=="REVIEW_SELECTED_NAME_MATCHED_GBIF_REQUEST"),
            "number_original_verified_HOSTS_host_species":len(compatible),
            "accepted_host_keys":";".join(sorted(compatible,key=int)),
            "pre2001_2010_indexed_records":int(x["pre2001_2010_indexed_records"]),
            "eligible_2011_2017":int(x["eligible_2011_2017"]),
            "first_2011_2017":int(x["first_indexed_in_2011_2017"]),
            "eligible_2018_2022":int(x["eligible_2018_2022"]),
            "first_2018_2022":int(x["first_indexed_in_2018_2022"]),
            "ever_host_plant_indexed_at_port_pre2011":int(bool(hit_p)),
            "ever_host_plant_indexed_at_offport_pre2011":int(bool(hit_o)),
            "ever_host_plant_indexed_anywhere_pre2011":int(bool(hit_p or hit_o)),
            "indexed_host_keys_port":";".join(sorted(set().union(*(compatible&z["keys"] for z in hit_p)),key=int)),
            "indexed_host_keys_offport":";".join(sorted(set().union(*(compatible&z["keys"] for z in hit_o)),key=int)),
            "number_port_sites_with_matched_host":len(hit_p),
            "number_offport_sites_with_matched_host":len(hit_o),
            "number_port_sites_with_any_botanical_record":sum(z["n_flora_species"]>0 for z in port),
            "number_offport_sites_with_any_botanical_record":sum(z["n_flora_species"]>0 for z in off),
            "total_port_flora_species_index_summed_NOT_unique":sum(z["n_flora_species"] for z in port),
            "total_offport_flora_species_index_summed_NOT_unique":sum(z["n_flora_species"] for z in off),
            "butterfly_index_NOT_established_colonization":True,
            "nonindexed_botanical_host_NOT_biological_absence":True,
            "source_HOSTS_host_identity_NOT_in_situ_local_feeding":True
        }
        rows.append(row)
    assert len(rows)==595
    assert len(set(x["species"] for x in rows))==35
    assert all(len([x for x in rows if x["landmass_id"]==i])==35 for i in regions)
    return rows,taxonomy,regions,created

def run_null(rows,outcome,block_region,draws=999,attempts_per_draw=70):
    """Fixed butterfly × original region event and island event margins.

    A contingency-table double-edge switch preserves species event count within
    region, island first-record count, and eligibility. Fixed host exposure is
    not shuffled across species with different known host resource portfolios.
    """
    eligible=[r for r in rows if r["eligible_2011_2017"] if outcome=="first_2011_2017"] if False else [
       r for r in rows if (r["eligible_2011_2017"] if outcome=="first_2011_2017" else r["eligible_2018_2022"])]
    occupied={(r["species"],r["landmass_id"]) for r in eligible if r[outcome]==1}
    allowed={(r["species"],r["landmass_id"]) for r in eligible}
    exposure={(r["species"],r["landmass_id"]):r["ever_host_plant_indexed_anywhere_pre2011"]
        for r in eligible}
    region={(r["species"],r["landmass_id"]):r["component"] for r in eligible}
    nobs=sum(exposure[e] for e in occupied)
    if not occupied:return {"observed":0,"null_mean":None,"null_support":[],"accepted_swaps":0,"status":"NO_EVENTS"}
    rng=random.Random(SEED+9*(outcome=="first_2018_2022")+13*block_region)
    groups=defaultdict(list)
    for pair in sorted(occupied):
        idx=region[pair] if block_region else 0
        groups[idx].append(pair)
    nonempty=[g for g in groups if len(groups[g])>=2]
    if not nonempty:
        return {"observed":nobs,"null_mean":None,"null_support":[],"accepted_swaps":0,
                "status":"NO_SWAPPABLE_GROUPS"}
    accepted=0
    def swaps(n):
        nonlocal accepted
        for _ in range(n):
            g=rng.choice(nonempty);v=groups[g]
            i=rng.randrange(len(v));j=rng.randrange(len(v))
            if i==j:continue
            a,p=v[i];b,q=v[j]
            if a==b or p==q:continue
            xy=(a,q);yx=(b,p)
            if xy not in allowed or yx not in allowed or xy in occupied or yx in occupied:continue
            occupied.remove((a,p));occupied.remove((b,q))
            occupied.add(xy);occupied.add(yx)
            v[i]=xy;v[j]=yx
            accepted+=1
    swaps(5000)
    vals=[]
    for _ in range(draws):
        swaps(attempts_per_draw)
        vals.append(sum(exposure[e] for e in occupied))
    if accepted<100:
        verdict="UNDERIDENTIFIED_FEWER_THAN_100_ACCEPTED_SWAP_MOVES"
    else:verdict="DESCRIPTIVE_CONDITIONAL_MARGINAL_NULL"
    return {"observed":nobs,"null_mean":sum(vals)/len(vals),
        "null_min":min(vals),"null_max":max(vals),
        "tail_draws_ge_observed":sum(x>=nobs for x in vals),
        "tail_fraction_ge_observed":sum(x>=nobs for x in vals)/draws,
        "draws":draws,"accepted_swaps":accepted,"status":verdict,
        "preserved_species_event_margins":True,
        "preserved_island_event_margins":True,
        "preserved_within_region_margins":bool(block_region),
        "not_independent_MCMC_pvalue_or_causal_effect":True}

def summarize(rows):
    result={"schema":"chocho_35_footprint_HOSTS6_historical_port_plant_butterfly_first_index_v1",
        "status":"PUBLIC_EVENT_RISK_AUDIT_NO_INDEPENDENT_HARBOR_HOST_COLONIZATION_CAUSAL_INFERENCE",
        "all_unselected_butterfly_panel_species":225,"all_original_islands":17,
        "all_unselected_butterfly_risk_cells":3825,
        "selected_explicit_public_verified_HOSTS_only_species":35,
        "selected_explicit_public_verified_HOSTS_edges":42,
        "selected_explicit_public_verified_host_plant_taxa":6,
        "joined_species_island_cells":len(rows),
        "source_git_blobs":{name:{"path":v[0],"sha":v[1]} for name,v in SOURCES.items()},
        "classification_no_host_index_means_UNKNOWN_not_absent":True,
        "plant_data_2001_2010_before_butterfly_outcomes":True,
        "butterfly_selection_before_future_records_only_FOR_225_not_host_subset":True,
        "host_subset_selected_from_port_related_previous_analysis_NOT_random":True,
        "not_independent_actual_H2_establishment_or_ship_transport":True,
        "periods":{},
        "source_limitation":[
          "35 host-linked butterflies are source-selected from prior port study and not representative of all 225 source butterflies.",
          "2001-2010 plant speciesKey observations are 1km WPI port and two offport buffers, NOT complete historical island flora.",
          "No indexed host at sampled sites means unobserved, NOT biological absence.",
          "Country-level POWO introduced status exists for only named focal plants; not all 6 taxa are island-specific historically alien.",
          "2019 WPI facilities are retrospective 2001-2010 proxies, not original shipment pathways.",
          "GBIF 2026 index uses revised family taxonKey search and simplified original GSHHG geometry, so it is not the private original record chronology.",
          "A 2011-2017 butterfly index is a GBIF presence record, not first actual arrival or breeding.",
          "2018-2022 source cohort is selected by 2011-2017 nonindex and is not species holdout.",
          "Observed host-matching counts reuse 42 selected host edges, and sparse plant effort/region confound comparisons.",
          "Region-blocked double-edge-swap null tests conditional event placement, not a causal independent p-value."
        ]
    }
    for outcome,eligible in [
       ("first_2011_2017","eligible_2011_2017"),
       ("first_2018_2022","eligible_2018_2022")]:
        at=[r for r in rows if r[eligible]==1]
        n_events=sum(r[outcome] for r in at)
        estimates={}
        for exp in ("ever_host_plant_indexed_at_port_pre2011",
                    "ever_host_plant_indexed_at_offport_pre2011",
                    "ever_host_plant_indexed_anywhere_pre2011"):
            yes=[x for x in at if x[exp]]
            no=[x for x in at if not x[exp]]
            y=sum(x[outcome] for x in yes);n=sum(x[outcome] for x in no)
            estimates[exp]={
              "host_indexed_risk_cells":len(yes),"host_indexed_future_first_events":y,
              "not_indexed_UNKNOWN_risk_cells":len(no),"not_indexed_UNKNOWN_future_first_events":n,
              "raw_future_first_index_rate_at_host_index":y/len(yes) if yes else None,
              "raw_future_first_index_rate_at_unindexed_unknown":n/len(no) if no else None,
              "unobserved_host_NOT_absent":True
            }
        by_island=defaultdict(list)
        for r in at:by_island[r["landmass_id"]].append(r)
        result["periods"][outcome]={
          "risk_cells":len(at),"first_indexed_events":n_events,
          "eligible_islands":len(by_island),
          "eligible_species":len(set(r["species"] for r in at)),
          "documented_migrant_positive_only_count":sum(
             r["is_migrant_roster_documented"] for r in at),
          "compatible_host_index":estimates,
          "by_island":[{"island":island,"risk_cells":len(xs),
                       "first_indexed_events":sum(r[outcome] for r in xs),
                       "host_port_positive_risk_cells":sum(r["ever_host_plant_indexed_at_port_pre2011"] for r in xs),
                       "host_offport_positive_risk_cells":sum(r["ever_host_plant_indexed_at_offport_pre2011"] for r in xs)}
                      for island,xs in sorted(by_island.items(),key=lambda q:int(q[0]))],
          "null_global_species_and_island_margins":run_null(rows,outcome,False),
          "null_within_six_original_regions_species_and_island_margins":run_null(rows,outcome,True)
        }
    return result

def write_csv(rows,path):
    with path.open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader();w.writerows(rows)

def main():
    rows,taxa,regions,created=join_rows()
    report=summarize(rows)
    HERE.mkdir(parents=True,exist_ok=True)
    write_csv(rows,HERE/"full17_verified_42HOSTS_butterfly_atrisk.csv")
    (HERE/"full17_verified_HOSTS_butterfly_risk_result.json").write_text(
        json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps({
       "status":report["status"],
       "rows":len(rows),"taxonomic_source_resolved_newly":created,
       "plant_accepted_taxa":taxa["species"],
       "periods":{k:{"risk_cells":v["risk_cells"],"first_indexed_events":v["first_indexed_events"],
            "compatible_host_index":v["compatible_host_index"],
            "region_null":v["null_within_six_original_regions_species_and_island_margins"]}
            for k,v in report["periods"].items()}
    },ensure_ascii=False,indent=2),flush=True)

if __name__=="__main__":main()
