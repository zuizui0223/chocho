#!/usr/bin/env python3
"""Source-frozen identifiability audit for island port-host-butterfly and 2020 AIS.

Do not fit a route/host/butterfly causal or predictive model if original
same-unit temporal, route-identity and observation-effort gates fail.
The 2020 voyages POSTDATE the 2001-2010 plants and 2011-2017 butterflies.
"""
from __future__ import annotations
import argparse, csv, hashlib, json
from collections import Counter,defaultdict
from pathlib import Path

ORIGINAL_GIT_BLOBS={
 "sites":"efb4e6c8648c3acd79f4416abc91b074eba1a07c",
 "observations":"ffc353c75c05e9e5fefa81342cd5be7314d035c7",
 "coverage":"82fbbdeda2e5d350654802ead129c2743082ea27",
 "direct_voyages":"b165c00d70c8ca3c1897745513c850247e44bec3",
 "voyage_source":"ed0b74df50ed6fe3ec02735fd43508e625691c83",
 "original_site_audit":"cb349f6359603def39db18b36a7208cadac530bc",
}
ORIGINAL_SITE_PAIRS=24
ORIGINAL_ISLANDS=17
ORIGINAL_SITE_ROWS=48
ORIGINAL_WPI_PORTS=194
ORIGINAL_VERIFIED_OD_EDGES=94
GATE_MIN_PAIRS=15
GATE_MIN_ISLANDS=10

def blob_sha(p):
    b=p.read_bytes()
    return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()

def read_csv(p):
    with p.open(encoding="utf-8-sig",newline="") as f:
        return list(csv.DictReader(f))

def unique_index(rs,cols):
    d={}
    for r in rs:
        k=tuple(r[c] for c in cols)
        if k in d:raise ValueError(f"Nonunique frozen {cols}: {k}")
        d[k]=r
    return d

def num(r,col):
    return int(r[col])

def effort_ratio_match(port,nonport,field,at_least=5):
    pa=num(port,field);np=num(nonport,field)
    return pa>=at_least and np>=at_least and .25<=pa/np<=4

def main():
    p=argparse.ArgumentParser()
    for k in ("sites","observations","coverage","direct-voyages","voyage-source","original-site-audit","out-dir"):
        p.add_argument("--"+k,required=True,type=Path)
    a=p.parse_args()
    paths={k.replace("-","_"):getattr(a,k.replace("-","_")) for k in
           ("sites","observations","coverage","direct-voyages","voyage-source","original-site-audit")}
    for k,path in paths.items():
        want=ORIGINAL_GIT_BLOBS[k]
        if blob_sha(path)!=want:
            raise ValueError(f"Fail closed: original data file {k} blob changed")
    sites=read_csv(a.sites)
    obs=read_csv(a.observations)
    covered=read_csv(a.coverage)
    edges=read_csv(a.direct_voyages)
    site_audit=json.loads(a.original_site_audit.read_text())
    ship=json.loads(a.voyage_source.read_text())
    if len(sites)!=ORIGINAL_SITE_PAIRS or len(obs)!=ORIGINAL_SITE_ROWS:
        raise RuntimeError("Frozen source 24 matched port/nonport pairs or 48 site rows changed")
    if len(covered)!=ORIGINAL_WPI_PORTS or len(edges)!=ORIGINAL_VERIFIED_OD_EDGES:
        raise RuntimeError("Frozen 194 WPI port source coverage or 94 verified directed links changed")
    if len(set(r["landmass_id"] for r in sites))!=ORIGINAL_ISLANDS:
        raise RuntimeError("Original cohort 17 islands changed")
    if site_audit["cohort"]["matched_port_nonport_pairs"]!=24 or site_audit["cohort"]["matched_islands"]!=17:
        raise RuntimeError("Original 2001-10 to 2011-22 island-source audit not reproduced")
    if (ship["verified_original_WPI_ports_with_one_or_more_external_endpoints"]!=46 or
        ship["distinct_original_landmass_unordered_pairs_with_observed_2020_voyage"]!=12 or
        ship["source_2020_departure_range"][0][:4]!="2020"):
        raise RuntimeError("Independent 2020 Zenodo WPI coverage source changed")
    obs_i=unique_index(obs,("landmass_id","port_id","side"))
    coverage_i=unique_index(covered,("original_WPI_port_id",))
    by_port=defaultdict(lambda:{"cross_voyage_rows":0,"within_voyage_rows":0,"other_islands":set()})
    for e in edges:
        c=num(e,"source_observed_2020_triplet_count")
        if c<=0:raise ValueError("Verified observed voyage source count nonpositive")
        same=e["same_original_island"]=="True"
        if same !=(e["origin_island_id"]==e["dest_island_id"]):
            raise RuntimeError("WPI original landmass ID inconsistent with its own record")
        for col,other in (("origin_WPI_port_id","dest_island_id"),("dest_WPI_port_id","origin_island_id")):
            port=e[col]
            if same:by_port[port]["within_voyage_rows"]+=c
            else:
                by_port[port]["cross_voyage_rows"]+=c
                by_port[port]["other_islands"].add(e[other])

    result=[]
    for site in sites:
        island=site["landmass_id"];port=site["port_id"]
        pkey=(island,port,"port");nkey=(island,port,"nonport")
        if pkey not in obs_i or nkey not in obs_i or (port,) not in coverage_i:
            raise RuntimeError(f"Original port/nonport/2020 source key no longer aligns: {pkey}")
        wp=coverage_i[(port,)]
        if wp["island_id"]!=island:
            raise RuntimeError(f"Original coastal site and verified WPI source have different landmass IDs for {port}")
        j=by_port[port]
        verified=wp["verified_external_2020_endpoint"]=="True"
        if not verified and j["cross_voyage_rows"]>0:
            raise RuntimeError("A source-unmatched WPI port was assigned real shipping OD")
        x=obs_i[pkey];y=obs_i[nkey]
        plant_ok=(num(x,"gbif_plant_pre_total_records")>=10 and
                  num(y,"gbif_plant_pre_total_records")>=10)
        pre_ok=effort_ratio_match(x,y,"gbif_butterfly_pre_total_records",5)
        early_ok=effort_ratio_match(x,y,"gbif_butterfly_early_total_records",10)
        result.append({
           "landmass_id":island,"port_id":port,"source_WPI_port_name":wp["source_port_name"],
           "GBIF_both_plant_sources_at_least_10":plant_ok,
           "GBIF_prior_butterfly_effort_at_least_5_within_fourfold":pre_ok,
           "GBIF_2011_2017_butterfly_effort_at_least_10_within_fourfold":early_ok,
           "combined_potentially_evaluable_effort":plant_ok and pre_ok and early_ok,
           "pre_plant_records_port":num(x,"gbif_plant_pre_total_records"),
           "pre_plant_records_nonport":num(y,"gbif_plant_pre_total_records"),
           "pre_butterfly_records_port":num(x,"gbif_butterfly_pre_total_records"),
           "pre_butterfly_records_nonport":num(y,"gbif_butterfly_pre_total_records"),
           "early_butterfly_records_port":num(x,"gbif_butterfly_early_total_records"),
           "early_butterfly_records_nonport":num(y,"gbif_butterfly_early_total_records"),
           "early_225_new_species_port":num(x,"gbif_butterfly_new_2011_2017_225_species"),
           "early_225_new_species_nonport":num(y,"gbif_butterfly_new_2011_2017_225_species"),
           "WPI_2020_exact_port_name_plus_distance_verified":verified,
           "WPI_2020_verified_direct_cross_island_source_episode":j["cross_voyage_rows"]>0,
           "WPI_2020_indexed_cross_island_voyage_rows":j["cross_voyage_rows"],
           "WPI_2020_indexed_within_island_voyage_rows":j["within_voyage_rows"],
           "WPI_2020_other_original_island_ids":"|".join(sorted(j["other_islands"],key=int)),
           "AIS2020_after_2001_2010_plants_and_2011_2017_butterfly_response":True,
        })
    result.sort(key=lambda r:(int(r["landmass_id"]),int(r["port_id"])))
    if len(result)!=24:raise RuntimeError("Source index drift")
    cohorts={
      "all_original_port_nonport":lambda r:True,
      "plant_effort_both_at_least10":lambda r:r["GBIF_both_plant_sources_at_least_10"],
      "prior_butterfly_effort_both_at_least5_and_4fold":lambda r:r["GBIF_prior_butterfly_effort_at_least_5_within_fourfold"],
      "later_butterfly_effort_both_at_least10_and_4fold":lambda r:r["GBIF_2011_2017_butterfly_effort_at_least_10_within_fourfold"],
      "plant_and_prior_effort":lambda r:r["GBIF_both_plant_sources_at_least_10"] and r["GBIF_prior_butterfly_effort_at_least_5_within_fourfold"],
      "plant_and_later_effort":lambda r:r["GBIF_both_plant_sources_at_least_10"] and r["GBIF_2011_2017_butterfly_effort_at_least_10_within_fourfold"],
      "plant_and_prior_and_later_effort":lambda r:r["combined_potentially_evaluable_effort"],
      "2020_original_port_verified":lambda r:r["WPI_2020_exact_port_name_plus_distance_verified"],
      "2020_cross_island_voyage_source_confirmed":lambda r:r["WPI_2020_verified_direct_cross_island_source_episode"]
    }
    stats={}
    for name,condition in cohorts.items():
        s=[r for r in result if condition(r)]
        stats[name]={
          "matched_port_nonport_pairs":len(s),
          "independent_original_landmass_groups":len(set(r["landmass_id"] for r in s)),
          "landmass_ids":sorted(set(r["landmass_id"] for r in s),key=int),
          "gate_minimum_15pairs_and_10islands_satisfied":
                   len(s)>=GATE_MIN_PAIRS and len(set(r["landmass_id"] for r in s))>=GATE_MIN_ISLANDS
        }
    if (stats["2020_original_port_verified"]["matched_port_nonport_pairs"]!=6 or
        stats["2020_cross_island_voyage_source_confirmed"]["matched_port_nonport_pairs"]!=2 or
        stats["plant_and_prior_and_later_effort"]["matched_port_nonport_pairs"]!=3):
        raise RuntimeError("Observed source-cohort overlap changed: inspect upstream data before proceeding")
    if stats["2020_cross_island_voyage_source_confirmed"]["gate_minimum_15pairs_and_10islands_satisfied"]:
        raise RuntimeError("Guard failed: biological OD model threshold unexpectedly attained")
    if stats["plant_and_prior_and_later_effort"]["gate_minimum_15pairs_and_10islands_satisfied"]:
        raise RuntimeError("Guard failed: effort-matched biological cohort unexpectedly adequate")
    a.out_dir.mkdir(parents=True,exist_ok=True)
    with (a.out_dir/"original24_ports_WPI_AIS2020_GBIF_detection_join.csv").open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=list(result[0]));w.writeheader();w.writerows(result)
    audit={
     "schema":"chocho_original_24_site_17_island_port_host_butterfly_ship_2020_identifiability_20261117_v1",
     "status":"SAME_PORT_ID_SOURCE_FEASIBILITY_GATES_FAILED_NO_SHIPPING_OR_ESTABLISHMENT_MODEL",
     "source_git_blobs":ORIGINAL_GIT_BLOBS,
     "GBIF_port_nonport_cohort":"24 matched WPI ports, 17 original landmass groups; 2001-10 plant and butterflies, 2011-17 and 2018-22 indexed butterflies",
     "AIS_source":"Duvallet et al. 2024 Zenodo 10.5281/zenodo.10783698; source departures 2020-01-01 through 2020-10-07 only",
     "source_WPI_matching":"2020 external port exact normalized name AND <=25km original WPI position, no inferred geographic alias",
     "ship_2020_later_than_early_butterfly_response_2011_2017":True,
     "original_source_island_port_cohort":stats,
     "minimum_model_support":{"min_distinct_islands":GATE_MIN_ISLANDS,"min_matched_ports":GATE_MIN_PAIRS,
         "port_host_butterfly_and_2020_cross_island_route_model":"NOT_IDENTIFIABLE",
         "effort_adjusted_host_to_butterfly_lag_model":"NOT_IDENTIFIABLE"},
     "recorded_voyage_source_is_not_historical_absence":True,
     "no_model_was_fitted":True,
     "important_limits":[
        "No link from 2020 ship OD to a biological event indexed before 2020 can be given a causal time order.",
        "Only 6 original ports are conservatively identified in AIS2020 and only 2 have original-island-crossing WPI-matched episodes.",
        "Missing WPI-2020 link is unknown coverage, never zero real ship traffic.",
        "Pre-2011 plant and butterfly recording effort does not have a target-group recording-day denominator at each local site.",
        "Conditioning on 2011-17 butterfly record counts is a diagnostic of future observation coverage, not a valid prospective predictor or causal adjustment.",
        "A coarse host candidate list is not proof of locally available host or larval feeding.",
        "Earlier butterfly GBIF first index is not biological new colonisation or establishment.",
        "Do not treat three matched port samples on one island as three independent spatial replicates.",
        "The original globally positive host geography and 155-species GBIF prediction results remain separate validated estimands."
     ],
    }
    (a.out_dir/"source_matched_route_effort_identifiability_audit.json").write_text(
         json.dumps(audit,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print("SOURCE_PORT_HOST_VOYAGE_EFFORT_AUDIT",json.dumps({
        "status":audit["status"],"cohorts":stats,
        "no_model_was_fitted":True
    },ensure_ascii=False),flush=True)

if __name__=="__main__": main()
