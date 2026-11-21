#!/usr/bin/env python3
"""Read-only 2026 NEW GBIF provider-facet audit for seven original WPI ports.

Original 1,596 complete port cache: GBIF Plantae speciesKey facet pooled
2001-2024, historical original sources have no occurrence gbifID/datasetKey.

This pilot uses seven ports from the original 1,596 completed source and the
independent Nov 7 archived exact 5km geoDistance source of 194 port centres.
Two periods: 2001-2010 (compares to Nov7 at the same query parameters) and
2001-2024 (new provider index only; old 1,596 cache used 5km WKT geometry
rather than necessarily equivalent GBIF geoDistance).

NEVER substitute current datasetKeys for historical publishers or invent DOI.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, time, uuid
from pathlib import Path
from datetime import datetime, timezone
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError

ORIGINAL="3dbf038dab64e8bbc9ba6c0f5b0e19ce1cb2b250"
P194=Path("results/island51_pre2011_public_general_flora_20261107/source/GBIF_original194_port_2001_2010_coverage.csv")
P1596=Path("provenance/reviewer_defenses/results/port_alien_broad_inventory_20261030/port_alien_host_broad_all_ports.csv")
SOURCE_REC=Path("results/island51_pre2011_public_general_flora_20261107/source/source_receipt.json")
SOURCE_COV_SHA="ec172084494d73284fdbaeaacd59d23c030d390884ca8d409f2f23e96d8559fb"
SOURCE_PORT_BLOB="fdf73b11ddcb1644f46730123350ae52cc61cae6"
IDS=(2,28,33,145,425,428,609)
PERIODS={"early2001_2010": "2001,2010", "pooled2001_2024": "2001,2024"}
GBIF="https://api.gbif.org/v1/occurrence/search"

def sha(raw:bytes)->str:return hashlib.sha256(raw).hexdigest()
def blob(raw:bytes)->str:
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()

def read_csv(path):
    with path.open(encoding="utf-8-sig",newline="") as f:return list(csv.DictReader(f))

def query(params):
    url=GBIF+"?"+urlencode(params,doseq=True)
    error=None
    for attempt in range(5):
        try:
            request=Request(url,headers={"Accept":"application/json",
                "User-Agent":"chocho-20261120-plant-port-provider-facet-source-probe"})
            with urlopen(request,timeout=75) as response:
                if response.status!=200:raise RuntimeError("GBIF non-200")
                j=json.load(response)
            if not isinstance(j,dict) or j.get("error"):
                raise ValueError("Unusable GBIF JSON")
            return j
        except (HTTPError,URLError,TimeoutError,ValueError,RuntimeError) as e:
            error=e
            if isinstance(e,HTTPError) and e.code not in (408,425,429,500,502,503,504):
                break
            if attempt<4:time.sleep(min(12,2**attempt))
    raise RuntimeError(f"GBIF query unavailable (not a zero record): {type(error).__name__}: {str(error)[:180]}")

def facet_data(params):
    values={}
    total=None;offset=0;pages=0
    for _ in range(80):
        j=query({**params,"facetOffset":offset})
        n=int(j["count"])
        if total is not None and n!=total:
            raise RuntimeError("GBIF indexed count changed during provider paging")
        total=n
        facets={str(x["field"]).upper():x.get("counts",[]) for x in j.get("facets",[])}
        if total>0 and "DATASET_KEY" not in facets:
            raise RuntimeError("GBIF positive count but no datasetKey facet")
        part=facets.get("DATASET_KEY",[])
        for r in part:
            k=str(r.get("name",""));cnt=int(r["count"])
            if str(uuid.UUID(k))!=k.lower() or cnt<0 or k in values:
                raise RuntimeError("Invalid or duplicated datasetKey facet")
            values[k]=cnt
        pages+=1
        if len(part)<1000:break
        offset+=len(part)
    else:raise RuntimeError("Too many GBIF source contributor pages")
    if sum(values.values())!=total:
        raise RuntimeError("DatasetKey contribution counts fail exact total reconciliation")
    return total,values,pages

def audit(src:Path):
    coverbytes=(src/P194).read_bytes()
    portbytes=(src/P1596).read_bytes()
    if sha(coverbytes)!=SOURCE_COV_SHA:
        raise RuntimeError("194 original port-geometry/early-count immutable CSV changed")
    if blob(portbytes)!=SOURCE_PORT_BLOB:
        raise RuntimeError("Original 3630-port, 1596 complete-cohort Git blob changed")
    ports=read_csv(src/P1596)
    original=read_csv(src/P194)
    rr=json.loads((src/SOURCE_REC).read_text(encoding="utf8"))
    if rr["source_port_original_geo_count"]!=194 or rr["port_count_mismatches"]!=100 or rr["source_query_radius_km"]!=5:
        raise RuntimeError("Original Nov7 general-flora cohort changed")
    if rr["source_sha256"]["GBIF_original194_port_2001_2010_coverage.csv"]!=SOURCE_COV_SHA:
        raise RuntimeError("Nov7 source report does not agree with source CSV")
    if len(ports)!=3630 or sum(r["query_complete"]=="True" for r in ports)!=1596 or len(original)!=194:
        raise RuntimeError("Original 1596 / 3630 WPI or 194 baseline cohort changed")
    allport={int(r["port_id"]):r for r in ports}
    panel={int(r["port_id"]):r for r in original}
    if len(allport)!=len(ports) or len(panel)!=len(original):
        raise RuntimeError("Duplicate original port ID")
    for pid in IDS:
        if pid not in panel or pid not in allport or allport[pid]["query_complete"]!="True":
            raise RuntimeError("Original 194 cohort port not among 1596 completed GBIF ports")
    start=datetime.now(timezone.utc).isoformat()
    output=[]
    for pid in IDS:
        row=panel[pid]
        geo=row["query_5km_from_original_WPI_centroid"]
        if not geo.endswith(",5km") or len(geo.split(","))!=3:
            raise RuntimeError("Original 194 source geoDistance format changed")
        lat,lon,r=geo.split(",")
        if len(lat.split(".")[-1])!=7 or len(lon.split(".")[-1])!=7 or r!="5km":
            raise RuntimeError("Original pinned WPI port coordinate precision/radius changed")
        for period,years in PERIODS.items():
            params={"taxonKey":"6","year":years,
                    "occurrenceStatus":"PRESENT","hasCoordinate":"true",
                    "hasGeospatialIssue":"false","geoDistance":geo,
                    "facet":"datasetKey","facetLimit":"1000","limit":"0"}
            count,bykey,pages=facet_data(params)
            entry={
              "WPI_original_port_id":pid,
              "original_island_id":row["island_id"],
              "period":period,
              "source_2001_2010_GBIF_2026_11_07_index_count":int(row["live_GBIF_2001_2010_plant_records"]) if period=="early2001_2010" else None,
              "source_original_complete_1596_cache_query_count":"NOT_AVAILABLE_AS_RECORD_TOTAL",
              "source_2019_WPI_port_geoDistance_5km":geo,
              "NEW_2026_11_20_GBIF_index_records":count,
              "NEW_datasetKey_contributor_count":len(bykey),
              "NEW_datasetKey_index_counts":[{"datasetKey":k,"count":v} for k,v in sorted(bykey.items())],
              "pages":pages,
              "NEW_sum_datasetKey_records":sum(bykey.values()),
              "NEW_minus_Nov7_same_query_count":count-int(row["live_GBIF_2001_2010_plant_records"]) if period=="early2001_2010" else None,
              "original_1596_host_intro_candidate_taxa":float(allport[pid]["contemporary_not_native_in_catalogue"]),
              "not_directly_historical_host_subset":True
            }
            output.append(entry)
            print("PORT_NEW_DATASETKEY_FACETS",pid,period,"current",count,"providers",len(bykey),"earlyNov7",entry["source_2001_2010_GBIF_2026_11_07_index_count"],flush=True)
            time.sleep(0.12)
    finish=datetime.now(timezone.utc).isoformat()
    unique={q["datasetKey"] for x in output for q in x["NEW_datasetKey_index_counts"]}
    return {
      "schema":"chocho_1596_source_original_ports_194_subset_20261120_NEW_provider_facets_pilot_v1",
      "status":"SOURCE_PINNED_PORT_IDS_NEW_SNAPSHOT_PROVIDER_COUNTS_NOT_ORIGINAL_DATASET_ATTRIBUTION",
      "original_git_commit":ORIGINAL,
      "source_original_complete_WPI_ports":1596,
      "source_original_WPI_candidates":3630,
      "source_original_island_subset_ports":194,
      "source_original_nov7_early_plant_query_changed_ports":100,
      "GBIF_API":GBIF,
      "NEW_facet":"datasetKey",
      "original_1596_cached_facet":"speciesKey",
      "source_GBIF_year_periods":PERIODS,
      "NEW_query_time_start_UTC":start,"NEW_query_time_end_UTC":finish,
      "selected_WPI_port_ids":list(IDS),
      "source_immutable_194_geo_CSV_sha256":SOURCE_COV_SHA,
      "source_immutable_1596_port_cohort_git_blob":SOURCE_PORT_BLOB,
      "NEW_port_period_queries":len(output),
      "NEW_unique_datasetKeys_in_selected_ports":len(unique),
      "results":output,
      "DOI_for_2026_original_1596_ports":None,
      "NEW_2026_index_not_2026_original_index":True,
      "all_NEW_provider_counts_match_NEW_GBIF_index_total":True,
      "not_comparable_to_original_cache_geometry":"The original 1596 cached buffer used WKT geometry; the Nov7 public subset and this new probe use source-identical 194-point GBIF geoDistance 5km, not guaranteed exactly same sampled buffer points.",
      "limits":[
       "The seven selected ports are a feasibility pilot (not 1596 representative independent samples)",
       "A 2001-2024 pooled query is not an independent timestamp of plant introduction, flora origin, larval host quality or actual ship transport",
       "Ports near one another can share GBIF records; summed port counts do not represent unique occurrences",
       "Current datasetKey list is not 2026 original port inventory source attribution and cannot mint a retrospective original derived DOI",
       "Attribution for 155/224 butterfly source DOI is a distinct acquisition, never borrowed for this botanical source"
      ]
    }

def main():
    a=argparse.ArgumentParser()
    a.add_argument("--source-root",type=Path,required=True)
    a.add_argument("--out",type=Path,required=True)
    args=a.parse_args()
    j=audit(args.source_root)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(j,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print("SEVEN_PORT_TWO_PERIODS_NEW_GBIF_DATASETKEY_FACET_SOURCE_CHECK_PASS",
        json.dumps({"query_count":j["NEW_port_period_queries"],
                    "unique_new_providers":j["NEW_unique_datasetKeys_in_selected_ports"]}))
if __name__=="__main__":
    main()
