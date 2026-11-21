#!/usr/bin/env python3
"""NEW read-only all 1,596 source WPI complete-port GBIF provider facets.

Only port IDs/source cohort are reproduced from pinned original 2026 source.
Current public WPI mirror point centres match the user's private original
2019 port point data when latitude/longitude are formatted to 6 decimals
(max 0.053 m difference independently confirmed). Original 1596 cached
facets used buffered WKT geometry, not necessarily exactly geoDistance 5km.

These 2026-11-20 queries measure NEW GBIF index, NOT original 2026 historic
GBIF publisher attribution. Do not mint a historic derived DOI, use full
GBIF original 155/224 butterfly DOI, or claim biological introduction.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,time,uuid
from datetime import datetime,timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request,urlopen
from urllib.error import HTTPError,URLError

BASE=Path("provenance/reviewer_defenses/results/port_alien_broad_inventory_20261030/port_alien_host_broad_all_ports.csv")
PORT_BLOB="fdf73b11ddcb1644f46730123350ae52cc61cae6"
ARCGIS_1596_SIX_DP_SHA256="91a1dd2a2d38ceb7244b9158a38e1ae9862298c644892d26fdf9d5318de2e98b"
ARCGIS_1596_CSV_SHA256="3783827ef5a0f95205bf51b0b7f186e055075d799f3ffeb2c0a7a4686f03446a"
GBIF="https://api.gbif.org/v1/occurrence/search"

def sha(b):return hashlib.sha256(b).hexdigest()

def get_json(params):
    url=GBIF+"?"+urlencode(params,doseq=True)
    error=None
    for attempt in range(6):
        try:
            req=Request(url,headers={"Accept":"application/json",
                "User-Agent":"chocho-GEB-20261120-source-matched-1596-WPI-port-GBIF-provider-facet"})
            with urlopen(req,timeout=90) as response:
                if response.status!=200:raise RuntimeError("GBIF non-200 response")
                j=json.load(response)
            if not isinstance(j,dict) or j.get("error"):
                raise ValueError("Unexpected GBIF data response")
            return j
        except (HTTPError,URLError,TimeoutError,ValueError,RuntimeError) as e:
            error=e
            if isinstance(e,HTTPError) and e.code not in (408,425,429,500,502,503,504):
                break
            if attempt<5:time.sleep(min(20,2**attempt))
    raise RuntimeError("GBIF provider-facet request failure, not a zero-source response: "
                       +type(error).__name__+" "+str(error)[:160])

def dataset_key_facets(params):
    all_counts={};total=None;offset=0;pages=0
    for _ in range(70):
        response=get_json({**params,"facetOffset":str(offset)})
        n=int(response["count"])
        if total is not None and n!=total:
            raise RuntimeError("Source index count drift within GBIF provider paging")
        total=n
        facets={str(x.get("field","")).upper():x.get("counts",[])
                for x in response.get("facets",[])}
        if n>0 and "DATASET_KEY" not in facets:
            raise RuntimeError("Positive GBIF index count without DATASET_KEY facet")
        records=facets.get("DATASET_KEY",[])
        for item in records:
            key=str(item.get("name",""));count=int(item.get("count",0))
            if str(uuid.UUID(key))!=key.lower() or count<0 or key in all_counts:
                raise RuntimeError("Invalid or duplicate GBIF provider datasetKey")
            all_counts[key]=count
        pages+=1
        if len(records)<1000:break
        offset+=len(records)
    else:raise RuntimeError("GBIF port has too many source publishers to facet")
    if sum(all_counts.values())!=total:
        raise RuntimeError("GBIF indexed record count not fully attributable to provider datasets")
    return total,all_counts,pages

def read_source(source_root,mirror):
    raw=(source_root/BASE).read_bytes()
    gitsha=hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
    if gitsha!=PORT_BLOB:raise RuntimeError("Original 2026 fixed 1596-cohort Git blob drift")
    rows=list(csv.DictReader(raw.decode("utf-8-sig").splitlines()))
    ids=sorted(int(r["port_id"]) for r in rows if r["query_complete"]=="True")
    if len(rows)!=3630 or len(ids)!=1596 or len(set(ids))!=1596:
        raise RuntimeError("Original WPI 3630/1596 status population drift")
    raw=mirror.read_bytes()
    if sha(raw)!=ARCGIS_1596_CSV_SHA256:
        raise RuntimeError("Public WPI checked 1596 mirror CSV is not exact immutable source artifact")
    d={};readrows=list(csv.DictReader(raw.decode("utf8").splitlines()))
    for r in readrows:
        p=int(r["port_id"])
        if p in d:raise RuntimeError("WPI mirror duplicate port ID")
        la,lo=float(r["latitude_7dp"]),float(r["longitude_7dp"])
        if not (-90<=la<=90 and -180<=lo<=180):raise RuntimeError("WPI point outside Earth")
        d[p]=(la,lo)
    if set(d)!=set(ids):raise RuntimeError("WPI mirror 1596 complete source port roster mismatch")
    six="port_id,lat,lon\n"+"".join(f"{k},{d[k][0]:.6f},{d[k][1]:.6f}\n" for k in ids)
    if sha(six.encode("utf8"))!=ARCGIS_1596_SIX_DP_SHA256:
        raise RuntimeError("WPI mirror does not match privately verified original WPI2019 6-decimal source")
    return ids,d

def audit(root,mirror,part,total_shards,checkpoint_dir):
    ids,coords=read_source(root,mirror)
    if not (0<=part<total_shards) or total_shards!=4:
        raise RuntimeError("Original source deterministic 4-way shard selection changed")
    selected=ids[part::total_shards]
    if len(selected)!=399:
        raise RuntimeError("Original complete-port subgroup population is not 399")
    checkpoint_dir.mkdir(parents=True,exist_ok=True)
    ts=datetime.now(timezone.utc).isoformat()
    output=[]
    for idx,k in enumerate(selected,start=1):
        lat,lon=coords[k]
        q={"taxonKey":"6","year":"2001,2024",
           "occurrenceStatus":"PRESENT","hasCoordinate":"true",
           "hasGeospatialIssue":"false",
           "geoDistance":f"{lat:.7f},{lon:.7f},5km",
           "facet":"datasetKey","facetLimit":"1000","limit":"0"}
        count,facets,pages=dataset_key_facets(q)
        output.append({
          "original_WPI2019_OBJECTID":k,
          "new_2026_GBIF_index_count":count,
          "new_datasetKey_count":len(facets),
          "current_GBIF_record_counts_by_datasetKey":[
             {"datasetKey":id,"count":n} for id,n in sorted(facets.items())],
          "facet_pages":pages,
          "new_source_query_geoDistance":q["geoDistance"]
        })
        if idx%40==0 or idx==len(selected):
            print("NEW_1596_PORT_PROVIDER_PROGRESS",part,idx,len(selected),
                  "records",sum(x["new_2026_GBIF_index_count"] for x in output),
                  "datasets",len({z["datasetKey"] for x in output for z in x["current_GBIF_record_counts_by_datasetKey"]}),
                  flush=True)
            # Internal checkpoint: partial shard is NEVER promoted as 1596 completion.
            (checkpoint_dir/("PARTIAL_port_shard_"+str(part)+".json")).write_text(
              json.dumps({"source":"PARTIAL_ONLY_DO_NOT_USE_FOR_1596_ATTRIBUTION",
                          "source_checked_rows":len(output),"expected_rows":399,
                          "new_index_timestamp_start":ts,"data":output},
                          ensure_ascii=False)+"\n",encoding="utf8")
        time.sleep(.16)
    if len(output)!=399 or len({x["original_WPI2019_OBJECTID"] for x in output})!=399:
        raise RuntimeError("Source selected port shard not complete")
    out={
       "schema":"chocho_new_20261120_all1596_WPI_plants_GBIF_datasetkey_facets_shard_v1",
       "status":"SOURCE_COMPLETED_ONE_OF_FOUR_NEW_PUBLIC_GBIF_PROVIDER_SHARDS_NOT_ORIGINAL_2026_QUERY",
       "original_git_source_commit":"3dbf038dab64e8bbc9ba6c0f5b0e19ce1cb2b250",
       "original_all_ports":3630,"original_complete_ports":1596,
       "source_original_private_WPI2019_location_match_6dp":True,
       "new_source_query_api":GBIF,
       "query_years":"2001,2024","query_radius":"5km",
       "query_filter":"Plantae taxonKey=6; PRESENT; hasCoordinate true; hasGeospatialIssue false",
       "new_provider_data_not_original_citation":True,
       "original_1596_cache_WKT_geometry_not_exactly_equivalent_to_new_geoDistance":True,
       "source_shard_index":part,"source_shard_total":total_shards,
       "new_index_timestamp_start_utc":ts,
       "new_index_timestamp_end_utc":datetime.now(timezone.utc).isoformat(),
       "n_port_queries":len(output),
       "new_indexed_port_buffers_sum_not_unique_occurrence_records":sum(x["new_2026_GBIF_index_count"] for x in output),
       "new_unique_datasetKey_in_shard":len({z["datasetKey"] for x in output for z in x["current_GBIF_record_counts_by_datasetKey"]}),
       "new_datasetKey_facets":output,
       "historical_2026_GBIF_publisher_datasetKeys_recovered":False,
       "new_GBIF_derived_DOI_minted":False,
       "caveat":"Provider shares are based on live new GBIF index (not the source's historical 1596-port speciesKey facet), and port buffers overlap."
    }
    return out

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-root",type=Path,required=True)
    p.add_argument("--public-mirror",type=Path,required=True)
    p.add_argument("--part",type=int,required=True)
    p.add_argument("--shards",type=int,default=4)
    p.add_argument("--checkpoint-dir",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    a=audit(args.source_root,args.public_mirror,args.part,args.shards,args.checkpoint_dir)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(a,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print("ALL1596_NEW_PROVIDER_FACET_SHARD_PASS",
          a["source_shard_index"],a["n_port_queries"],a["new_unique_datasetKey_in_shard"])
if __name__=="__main__":
    main()
