#!/usr/bin/env python3
"""Exact 2019 WPI port ID + landmass joint source coverage for island-first paper.

This is a join-support audit, not a botanical 'port vs non-port' comparison.
Existing GBIF flora was searched within a 5km port radius; 2021 ESA land-cover
sample is from 1km port/remote-coast probes. Matching WPI IDs DOES NOT
harmonize footprints, planting dates, recording effort or island introduced status.
"""
from __future__ import annotations
from collections import defaultdict
from pathlib import Path
import argparse,csv,hashlib,io,json,math

ROOT=Path("provenance/reviewer_defenses/results")
EARLY=ROOT/"island_dated_host_match_20261030/island_port_early_plant_availability.csv"
VEGETATION=ROOT/"coastal_port_vegetation_20261030/coastal_control_vegetation_pairs.csv"
ISLAND=ROOT/"island_host_network_20261030/island_host_network_islands.csv"
DATED=ROOT/"historical_botanical_17island_20261031/near_pair_pre2011_plant_composition.json"
SOURCE={
 str(EARLY):"0a01831056c2cd59a4bd41bd1dbd4c84d77cf922",
 str(VEGETATION):"10360b39e20a720000c5ac8f9f156073c67e7c9b",
 str(ISLAND):"18cecde5bf3066f276380167f451f1fa5cc1caa1",
 str(DATED):"54c4515947d2a3ee7545e4007985ef82a2af6868",
}
BINS=(0,1,10,30,100,300,1000)

def load(path,root):
 p=root/path;b=p.read_bytes()
 actual=hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()
 if actual!=SOURCE[str(path)]:raise ValueError("Pinned source changed "+str(path))
 return list(csv.DictReader(io.StringIO(b.decode("utf8")))) if p.suffix==".csv" else json.loads(b.decode("utf8"))

def summarize(root):
 flora=load(EARLY,root);cover=load(VEGETATION,root)
 sites=load(ISLAND,root);dated=load(DATED,root)
 eligible={s["landmass_id"] for s in sites if s["eligible"]=="True"}
 if len(eligible)!=51:raise ValueError("Original selected island group changed")
 datedids={str(x["island_id"]) for x in dated["islands"]}
 if len(datedids)!=17:raise ValueError("Historical island geography changed")
 original_ports=[p for p in flora if p["landmass_id"] in eligible]
 if len(original_ports)!=194:raise ValueError("Original 51 island 194 plant port inventory changed")
 cover1=[x for x in cover if x["landmass_id"] in eligible and x["near_radius_km"]=="1"]
 if len(cover1)!=273:raise ValueError("Original 51-island 1km control port census changed")
 if any(int(r["near_points"])<3 or int(r["far_points"])<3 for r in cover1):
  raise ValueError("Matched same island land-cover sites are inadequately sampled")
 cov={(r["landmass_id"],r["port_id"]):r for r in cover1}
 if len(cov)!=len(cover1):raise ValueError("Duplicate land-cover 1km port IDs")
 joint=[]
 for p in original_ports:
  key=(p["landmass_id"],p["port_id"])
  if key not in cov:continue
  c=cov[key];early=int(p["early_plant_records"])
  if early<0:raise ValueError("Negative indexed botanical recording effort")
  joint.append({
   "landmass_id":p["landmass_id"],"port_id":p["port_id"],
   "has_resolved_2001_2010_botanical_specieskey_year_data":p["landmass_id"] in datedids,
   "original_port_2001_2010_GBIF_indexed_plants_query_count":early,
   "original_port_2001_2024_GBIF_indexed_plants_query_count":int(p["all_plant_records"]),
   "landcover_2021_built_difference_1km_near_minus_own_island_nonport":float(c["built_difference"]),
   "landcover_2021_tree_difference_1km_near_minus_own_island_nonport":float(c["tree_difference"]),
   "landcover_2021_grass_difference_1km_near_minus_own_island_nonport":float(c["grass_difference"]),
   "all_controls_at_least_3_raster_samples_per_side":True
  })
 joint.sort(key=lambda x:(int(x["landmass_id"]),int(x["port_id"])))
 if len(joint)!=79:raise ValueError("Exact source WPI port join cohort changed")
 if len({x["port_id"] for x in joint})!=79:raise ValueError("Port ID collision across eligible landmass join")
 def ngrp(rows):return len({x["landmass_id"] for x in rows})
 bins=[]
 for t in BINS:
  cohort=[x for x in joint if x["original_port_2001_2010_GBIF_indexed_plants_query_count"]>=t]
  dated_cohort=[x for x in cohort if x["has_resolved_2001_2010_botanical_specieskey_year_data"]]
  bins.append({"minimum_2001_2010_indexed_botanical_query_records":t,
   "same_port_ID_plant_and_landcover_sites":len(cohort),
   "distinct_islands":ngrp(cohort),
   "historically_source_resolved_ports":len(dated_cohort),
   "historically_source_resolved_islands":ngrp(dated_cohort)})
 if [(r["same_port_ID_plant_and_landcover_sites"],r["distinct_islands"]) for r in bins][:5] != [
  (79,36),(67,32),(52,24),(44,19),(37,14)]:
  raise ValueError("Port-plant source support quality gate changed")
 return {
   "schema":"chocho_exact_WPI79_historical_port_plant_vs_same_port_landcover_source_join_v1",
   "status":"ID_LINKAGE_ONLY_NOT_PLANT_SITE_CONTROL",
   "git_blob_source_pins":SOURCE,
   "original_eligible_islands":len(eligible),
   "botanical_2001_2024_indexed_port_catalogue_ids":len(original_ports),
   "same_island_2021_landcover_matched_1km_port_ids":len(cover1),
   "identical_WPI_port_and_landmass_ID_intersection":len(joint),
   "identical_WPI_ports_in_distinct_islands":ngrp(joint),
   "historical_species_key_facets_from_2001_2010_available_port_IDs":sum(x["has_resolved_2001_2010_botanical_specieskey_year_data"] for x in joint),
   "historical_species_key_facets_from_2001_2010_available_islands":ngrp([x for x in joint if x["has_resolved_2001_2010_botanical_specieskey_year_data"]]),
   "pre2011_indexed_botanical_coverage_thresholds":bins,
   "exact_linked_source_2021_landcover_and_2001_2024_botanical_port_IDs":joint,
   "limits":{
     "original_botanical_port_radius_km":5,
     "same_port_2021_worldcover_coastal_probe_radius_km":1,
     "source_port_ID_match_does_not_mean_same_spatial_footprint":True,
     "nonport_site_plant_specieskey_records_not_acquired":True,
     "nonport_built_coast_landcover_not_equivalent_to_nonport_urbanized_coast":True,
     "original_botanical_taxonomy_species_by_year_not_available_at_all_79_port_IDs":True,
     "plant_record_count_not_alien_host_species_richness_or_butterfly_impact":True,
     "host_source_network_remains_separately_licensed":True,
     "unknown_island_local_plant_introduction_year_and_shipping_route":True,
     "no_independent_79_ports_claimed_when_multiple_same_island":True
   }}
if __name__=="__main__":
 p=argparse.ArgumentParser()
 p.add_argument("--root",type=Path,default=Path("."))
 p.add_argument("--out",type=Path,required=True)
 args=p.parse_args()
 z=summarize(args.root)
 args.out.parent.mkdir(parents=True,exist_ok=True)
 args.out.write_text(json.dumps(z,sort_keys=True,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
 print("JOINT_WPI_PORT_SOURCE_ID",json.dumps({"n":z["identical_WPI_port_and_landmass_ID_intersection"],"islands":z["identical_WPI_ports_in_distinct_islands"],"dated_port_ids":z["historical_species_key_facets_from_2001_2010_available_port_IDs"],"coverage":z["pre2011_indexed_botanical_coverage_thresholds"]}))
