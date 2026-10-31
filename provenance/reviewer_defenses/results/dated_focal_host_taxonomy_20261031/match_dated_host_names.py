#!/usr/bin/env python3
"""Cross known pre-window larval-host NAME references with new GBIF plant speciesKey×year facets.

This validates exact taxonomic source identity and indexed chronology only.
The actual original butterfly migration, reproduction, in-situ plant feeding,
historical introduced-status and ship transport are NOT measured here.
"""
from __future__ import annotations
import argparse,csv,gzip,hashlib,io,json,time
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen,Request
from urllib.error import HTTPError,URLError

ROOT=Path(__file__).resolve().parents[4]
HOST_LINKS=Path("provenance/reviewer_defenses/results/island_dated_host_match_20261030/island_dated_host_match_links.csv")
HOST_LINKS_GIT_BLOB="88db773c74dea9cd2d58bde0b67a394d93226869"
SOURCE_GZIP_SHA256="b436ba8668d67d05a247f9cafd59b846db2327031e59f4ad95c78d6870bf113e"
SOURCE_GZIP_MEMBER="local_all23ports_specieskey_year.csv.gz"
GBIF_MATCH_API="https://api.gbif.org/v1/species/match"

def fetch_json(url,params):
    full=url+"?"+urlencode(params)
    req=Request(full,headers={"User-Agent":"chocho-taxonomic-historical-port-index/20261031",
                              "Accept":"application/json"})
    with urlopen(req,timeout=50) as f:
        v=json.load(f)
    if not isinstance(v,dict):raise ValueError("Invalid GBIF taxonomic match response")
    return v

def pinned_links(root):
    raw=(root/HOST_LINKS).read_bytes()
    sha=hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
    if sha!=HOST_LINKS_GIT_BLOB:raise ValueError("Original locked butterfly-host name links changed")
    ls=list(csv.DictReader(io.StringIO(raw.decode("utf8"))))
    if len(ls)!=5:raise ValueError("Original five source host links changed")
    return ls

def taxon_match(name):
    d=fetch_json(GBIF_MATCH_API,{"name":name,"rank":"SPECIES","kingdom":"Plantae","strict":"true"})
    canonical=str(d.get("canonicalName","")).lower()
    if (d.get("matchType")!="EXACT" or d.get("rank")!="SPECIES"
        or canonical!=name.lower() or int(d.get("confidence",0))<90
        or not d.get("usageKey")):
        raise ValueError("Taxonomy match ambiguous for "+name+": "+json.dumps(d)[:500])
    accepted=int(d.get("acceptedUsageKey") or d["usageKey"])
    return {"query_name":name,"gbif_species_key":accepted,
            "exact_name_match":True,"status_at_query":d.get("status"),
            "rank":d.get("rank"),"confidence":d["confidence"],
            "usage_key":d["usageKey"],"accepted_usage_key":d.get("acceptedUsageKey")}

def run(root,original_gzip,target):
    if hashlib.sha256(original_gzip.read_bytes()).hexdigest()!=SOURCE_GZIP_SHA256:
        raise ValueError("All-23-port original historical plant species-year artifact hash changed")
    original=list(csv.DictReader(io.StringIO(gzip.decompress(original_gzip.read_bytes()).decode("utf8"))))
    if len(original)!=9788:raise ValueError("Original plant facet row count changed")
    links=pinned_links(root)
    focal=[x for x in links if sorted([x["island_a"],x["island_b"]]) in
            [["169","276"],["276","715"]]]
    if len(focal)!=3 or sorted(set(x["host_name"] for x in focal))!=["Alnus incana","Medicago sativa"]:
        raise ValueError("Expected original source-host focal 3 links")
    taxa={}
    for name in sorted(set(x["host_name"] for x in focal)):
        taxa[name]=taxon_match(name)
        time.sleep(.1)
    # Checked source endpoints are the *original* focal island host record pairs.
    results=[]
    for name,info in sorted(taxa.items()):
        key=str(info["gbif_species_key"])
        for iid in sorted({s for x in focal if x["host_name"]==name for s in (x["island_a"],x["island_b"])},key=int):
            matched=[x for x in original if x["island_id"]==iid and x["gbif_speciesKey"]==key]
            years=sorted({int(x["year"]) for x in matched})
            ports=sorted({x["port_id"] for x in matched},key=int)
            results.append(dict(host_name=name,gbif_species_key=key,
              island_id=iid,first_indexed_year=years[0] if years else None,
              indexed_years=years,number_annual_port_facets=len(matched),
              number_ports_with_matching_species_key=len(ports),
              port_ids=ports,
              sum_indexed_records_across_potentially_overlapping_port_buffers=sum(int(x["indexed_record_count"]) for x in matched),
              no_matching_row_means_NO_GBIF_INDEX_EVIDENCE_not_plant_absence=not years))
    output={
     "schema":"chocho_original_host_name_gbif_specieskey_annual_facets_v1",
     "status":"SOURCE_SPECIES_ID_AND_PRE_2011_BOTANICAL_RECORDS_NOT_CAUSAL_HOST_RESCUE",
     "source_original_host_links_git_blob":HOST_LINKS_GIT_BLOB,
     "source_all_23port_specieskey_year_gzip_sha256":SOURCE_GZIP_SHA256,
     "source_gbif_taxonomic_api":GBIF_MATCH_API,
     "checked_at_utc":datetime.now(timezone.utc).isoformat(),
     "known_historical_host_links_checked":len(focal),
     "unique_host_names_checked":len(taxa),
     "names":list(taxa.values()),
     "island_taxon_observations":results,
     "original_source_early_port_record_time_window":[2001,2010],
     "nonfocal_host_status_not_determined":True,
     "no_local_larval_feeding_inferred":True,
     "historical_native_vs_alien_status_NOT_VALIDATED":True,
     "GBIF_backbone_name_matching_not_original_WCVP_host_edge_certification":True,
     "not_documented_shipping_or_butterfly_source":True,
     "do_not_replace_missing_source_row_with_ecological_absence":True,
    }
    target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(output,indent=2,ensure_ascii=False,sort_keys=True)+"\n",encoding="utf8")
    print("FOCAL_TAXON_MATCH_AUDIT",json.dumps({"names":len(taxa),"island_records":len(results),
         "observed_prior_host_island_keys":sum(r["first_indexed_year"] is not None for r in results),
         "source_facets":sum(r["number_annual_port_facets"] for r in results)},sort_keys=True))

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root",type=Path,default=Path("."))
    p.add_argument("--source-gzip",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    run(args.root,args.source_gzip,args.output)
