#!/usr/bin/env python3
"""All 23 historically indexed ports across nine preselected island-study nodes.

Reads the same source-locked 2019 WPI IDs and GBIF query code as the 9-port pilot.
Port buffers may OVERLAP: only set-union of accepted species keys, never summed
record counts, is used as island botanical *indexed* richness. Not true flora.
"""
from __future__ import annotations
import argparse,csv,gzip,hashlib,json,time
from collections import defaultdict
from datetime import datetime,timezone
from pathlib import Path
from acquire_dated_botanical_facets import choose,get_ports,year_facets,PLANTS,PORTS,PINS,WPI,GBIF

def run(root,outfile,artifact):
    allports,chosen,countries=choose(root)
    geo=get_ports(allports,countries)
    portrows=[];yearrows=[];facets=[]
    union=defaultdict(set)
    years=defaultdict(set)
    indexed_record_sum=defaultdict(int)
    for row in sorted(allports,key=lambda x:int(x["port_id"])):
        iid=row["landmass_id"]
        pid=int(row["port_id"])
        port=geo[pid]
        n_sum=0
        for year in range(2001,2011):
            n,keys,pages=year_facets(port,year)
            n_sum+=n
            indexed_record_sum[iid]+=n
            union[iid].update(keys)
            if n:years[iid].add(year)
            yearrows.append({"island_id":iid,"port_id":pid,"year":year,
              "plant_records_indexed":n,"faceted_species_key_count":len(keys),
              "records_without_species_key":n-sum(keys.values()),"facet_pages":pages})
            for key,count in sorted(keys.items()):
                facets.append((iid,pid,year,key,count))
            time.sleep(0.18)
        portrows.append({"island_id":iid,"port_id":pid,
          "port_name":port["name"],"country":port["expected_country"],
          "original_early_record_index":int(row["early_plant_records"]),
          "new_annual_query_index_sum":n_sum,
          "difference_from_historic_source":n_sum-int(row["early_plant_records"])})
        print("ALL_PORT_BOTANICAL",iid,pid,n_sum,"earlier",row["early_plant_records"],flush=True)
    if len(yearrows)!=230 or len(portrows)!=23 or len(union)!=9:
        raise ValueError("23 × 10 year/port source census incomplete")
    artifact.parent.mkdir(parents=True,exist_ok=True)
    with gzip.open(artifact,"wt",newline="",encoding="utf8") as file:
        writer=csv.writer(file)
        writer.writerow(["island_id","port_id","year","gbif_speciesKey","indexed_record_count"])
        writer.writerows(facets)
    obs_by_island=[{
        "island_id":iid,
        "number_source_ports":sum(x["island_id"]==iid for x in portrows),
        "indexed_species_key_union_across_all_ports":len(union[iid]),
        "years_with_any_botanical_indexed_record":len(years[iid]),
        "sum_of_overlapping_port_buffer_record_queries_NOT_unique_occurrences":indexed_record_sum[iid],
        "not_original_island_wide_flora_inventory":True,
        } for iid in sorted(union,key=int)]
    manifest={
        "schema":"chocho_all23_early_port_botanical_coverage_v1",
        "status":"HISTORICAL_PORT_SPECIES_YEAR_INDEX_ACQUISITION_NOT_HOST_ABSENCE",
        "checked_at_utc":datetime.now(timezone.utc).isoformat(),
        "git_blob_pins":{str(path):sha for path,sha in PINS.items()},
        "wpi_mirror_identity_checked_ports":23,
        "source_ports_all":23,
        "source_islands_all":9,
        "annual_port_year_rows":230,
        "year_range":[2001,2010],
        "wpi_source":WPI,
        "gbif_occurrence_search":GBIF,
        "query_params":{"taxonKey":6,"hasCoordinate":True,
           "hasGeospatialIssue":False,"occurrenceStatus":"PRESENT",
           "geoDistance":"5km from 2019-derived WPI port point",
           "facet":"speciesKey","individual_years":True},
        "by_port":portrows,
        "by_island_species_key_unions":obs_by_island,
        "by_port_year_recording_support":yearrows,
        "specieskey_year_artifact":{"filename":artifact.name,
          "rows":len(facets),"sha256":hashlib.sha256(artifact.read_bytes()).hexdigest(),
          "not_committed_to_git":True},
        "overlapping_port_buffers_not_independent":True,
        "no_historical_introduced_host_status_inferred":True,
        "no_botanical_nonrecord_interpreted_as_absence":True,
        "GBIF_reindexed_and_geodistance_geometry_not_necessarily_old_cache_compatible":True,
        "no_butterfly_host_edges_or_original_occurrence_coordinates":True
    }
    outfile.parent.mkdir(parents=True,exist_ok=True)
    outfile.write_text(json.dumps(manifest,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf8")
    print("ALL_PORT_AUDIT",json.dumps({"ports":23,"islands":9,"portYears":230,"speciesKeyYearRows":len(facets),
      "speciesKeyUnions":{x["island_id"]:x["indexed_species_key_union_across_all_ports"] for x in obs_by_island}},sort_keys=True),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--root",type=Path,default=Path("."))
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--artifact",type=Path,required=True)
    a=p.parse_args()
    run(a.root,a.output,a.artifact)
