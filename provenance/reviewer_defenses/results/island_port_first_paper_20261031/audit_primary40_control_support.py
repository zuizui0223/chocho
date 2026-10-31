#!/usr/bin/env python3
"""Audit the 51-island butterfly + candidate host cohort against same-island
near-port and off-port *LAND COVER* controls.

This pilot does NOT measure larval hosts at non-port sites. Port botanical
records are pooled 2001-2024 and ESA WorldCover is from 2021, so no temporal
gateway effect or historical introductions can be claimed. Paired land-cover
samples are deliberately equally weighted by matched port, then by island.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,math
from pathlib import Path

BASE=Path("provenance/reviewer_defenses/results")
ISLANDS=BASE/"island_host_network_20261030/island_host_network_islands.csv"
COAST=BASE/"coastal_port_vegetation_20261030/coastal_control_vegetation_pairs.csv"
HISTORY=BASE/"historical_botanical_17island_20261031/near_pair_pre2011_plant_composition.json"
SOURCE={
 str(ISLANDS):"18cecde5bf3066f276380167f451f1fa5cc1caa1",
 str(COAST):"10360b39e20a720000c5ac8f9f156073c67e7c9b",
 str(HISTORY):"54c4515947d2a3ee7545e4007985ef82a2af6868"
}

def load(root,p):
 data=(root/p).read_bytes()
 sha=hashlib.sha1(b"blob "+str(len(data)).encode()+b"\0"+data).hexdigest()
 if sha!=SOURCE[str(p)]:raise ValueError("Source Git blob changed: "+str(p))
 return (list(csv.DictReader(io.StringIO(data.decode("utf8"))))
         if p.suffix==".csv" else json.loads(data.decode("utf8")))

def ranks(a):
 s=sorted(range(len(a)),key=lambda i:a[i]);r=[0.]*len(a)
 i=0
 while i<len(s):
  j=i+1
  while j<len(s) and a[s[j]]==a[s[i]]:j+=1
  for k in s[i:j]:r[k]=(i+j+1)/2
  i=j
 return r

def pearson(a,b):
 n=len(a);m=sum(a)/n;k=sum(b)/n
 s=sum((x-m)*(y-k) for x,y in zip(a,b))
 den=math.sqrt(sum((x-m)**2 for x in a)*sum((y-k)**2 for y in b))
 return s/den if den else None

def spearman(rows,x,y):return pearson(ranks([float(r[x]) for r in rows]),ranks([float(r[y]) for r in rows]))

def summarize(root):
 islands=load(root,ISLANDS);coast=load(root,COAST);history=load(root,HISTORY)
 source51={r["landmass_id"]:r for r in islands if r["eligible"]=="True"}
 if len(source51)!=51:raise ValueError("Original 51-island cohort changed")
 dated17={r["island_id"] for r in history["islands"]}
 if len(dated17)!=17:raise ValueError("Source-dated island cohort changed")
 pairs=[r for r in coast if r["near_radius_km"]=="1" and r["landmass_id"] in source51]
 if len(pairs)!=273:raise ValueError("Original 1 km port/non-port paired landcover rows changed")
 if len({r["port_id"] for r in pairs})!=273:raise ValueError("One source port appears with duplicates in 1km cohort")
 # Previous vegetation pipeline filters to at least three real raster points per matched side.
 if any(int(r["near_points"])<3 or int(r["far_points"])<3 for r in pairs):
  raise ValueError("Port or remote same-island landscape control has inadequate sampled points")
 g={}
 for r in pairs:
  key=r["landmass_id"]
  g.setdefault(key,[]).append(r)
 if len(g)!=40:raise ValueError("Eligible island port/nonport overlap changed")
 rows=[]
 for lid in sorted(g,key=int):
  source=source51[lid];ports=g[lid]
  def avg(x):return sum(float(v[x]) for v in ports)/len(ports)
  rows.append({"island_id":lid,"matched_port_sites":len(ports),
   "source_eligible_ports":int(source["ports"]),
   "contemporary_alien_host_candidate_species":int(source["alien_host_species"]),
   "prior_butterfly_species_indexed":int(source["earlier_butterfly_species"]),
   "later_new_butterfly_species_indexed":int(source["new_butterfly_species"]),
   "near_port_minus_same_island_off_port_built_fraction":avg("built_difference"),
   "near_port_minus_same_island_off_port_tree_fraction":avg("tree_difference"),
   "near_port_minus_same_island_off_port_grass_fraction":avg("grass_difference"),
   "near_port_minus_same_island_off_port_crop_fraction":avg("crop_difference"),
   "also_has_2001_2010_historical_plant_species_key_facets":lid in dated17})
 if sum(x["matched_port_sites"] for x in rows)!=273:raise ValueError("Original source port counts changed")
 if sum(x["also_has_2001_2010_historical_plant_species_key_facets"] for x in rows)!=12:
  raise ValueError("Historic 17 island subset / 40 island landcover overlap changed")
 tests=[
  ("built_fraction_delta_vs_aliens","near_port_minus_same_island_off_port_built_fraction","contemporary_alien_host_candidate_species"),
  ("built_fraction_delta_vs_later_butterflies","near_port_minus_same_island_off_port_built_fraction","later_new_butterfly_species_indexed"),
  ("built_fraction_delta_vs_prior_butterflies","near_port_minus_same_island_off_port_built_fraction","prior_butterfly_species_indexed"),
  ("total_island_aliens_vs_later_butterflies","contemporary_alien_host_candidate_species","later_new_butterfly_species_indexed"),
  ("total_island_aliens_vs_prior_butterflies","contemporary_alien_host_candidate_species","prior_butterfly_species_indexed"),
  ("built_fraction_delta_vs_total_source_ports","near_port_minus_same_island_off_port_built_fraction","source_eligible_ports")
 ]
 corr={k:spearman(rows,x,y) for k,x,y in tests}
 assert len(rows)==40
 return {"schema":"chocho_island_port_first_paper_control_support_and_exploratory_ranks_v1",
 "status":"LANDCOVER_VS_ISLAND_RICHNESS_ASSOCIATIONS_NOT_HOST_PORT_CONTROL",
 "git_blob_source_pins":SOURCE,
 "preexisting_eligible_islands":51,
 "islands_with_1km_matched_port_nonport_2021_landcover":len(rows),
 "sampled_1km_port_sites_in_eligible_islands":len(pairs),
 "historically_dated_plant_islands_2001_2010":len(dated17),
 "historically_dated_plant_islands_with_1km_controls":sum(r["also_has_2001_2010_historical_plant_species_key_facets"] for r in rows),
 "all_pairs_have_at_least_three_near_and_three_far_landcover_samples":True,
 "matched_1km_port_aggregate_island_rows":rows,
 "island_weighted_exploratory_spearman_correlations":corr,
 "limits":{
  "not_botanical_host_matched_port_vs_off_port":True,
  "not_local_plant_introduced_status":True,
  "2021_worldcover_not_historical_evidence":True,
  "port_plant_indices_pool_2001_2024_and_may_postdate_butterfly_records":True,
  "not_independent_port_replicates_when_same_island":True,
  "no_climate_habitat_reporting_effort_adjustment":True,
  "does_not_test_actual_ship_supply_or_local_larval_breeding":True,
  "does_not_establish_55_percent_global_geb_resource_realization":True},
 "recommended_primary_next_data":"Same sampled 40 islands: original-coordinate and time-window matched port-versus-nonport botanical speciesKey records, accepted larval host identities, historical alien status, provider effort and local urban controls. First restrict to 12 islands with dated prior botanical source; expand only if controls feasible."}

if __name__=="__main__":
 p=argparse.ArgumentParser()
 p.add_argument("--root",type=Path,default=Path("."))
 p.add_argument("--output",type=Path,required=True)
 args=p.parse_args();audit=summarize(args.root)
 args.output.parent.mkdir(parents=True,exist_ok=True)
 args.output.write_text(json.dumps(audit,sort_keys=True,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
 print("ISLAND_PORT_PRIMARY_GATE",json.dumps({
   "islands":audit["islands_with_1km_matched_port_nonport_2021_landcover"],
   "ports":audit["sampled_1km_port_sites_in_eligible_islands"],
   "dated_islands":audit["historically_dated_plant_islands_with_1km_controls"],
   "rho":audit["island_weighted_exploratory_spearman_correlations"]},sort_keys=True))
