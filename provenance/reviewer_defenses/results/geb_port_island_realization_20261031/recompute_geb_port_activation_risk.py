#!/usr/bin/env python3
"""Source-locked 2011–17 vs 2018–22 butterfly-island first-record risk-set comparison.

The butterfly×island×year PARQUET belongs to the original user and is NOT
included in GitHub; this script is a controlled private-source replay recipe.
Original host interactions (2,524) and individual plant records remain private.

Question: Are first-indexed butterflies seen more often at islands where a
catalogue-compatible larval host plant was indexed at a 2019-WPI proxy
port buffer in 2001–2010? Compare to species×geographic-group + island
marginal fixed permutations among previously non-indexed risk cells.

NO inference of shipping, local larval use, actual introduction, arrival or
breeding. The 60 butterflies are selected for >=1 first-index event in the
2011–17 study areas; the 2018–22 window is not a fully blind species cohort.
"""
from __future__ import annotations
from pathlib import Path
from collections import defaultdict,Counter
import argparse,csv,gzip,hashlib,io,json,random,math

ROOT=Path(__file__).resolve().parents[4]
BASE=Path("provenance/reviewer_defenses/results")
ORIGINAL="8362ec4538d28916b2f9aaf98d1512b343b15f8f2eefe060ee61798c9a67406f"
HOST_EDGES="e55012d6b0eefd2b3e9e4b8f13fc4786ddc1e66c21f73e705bb20df92d0b5e22"
TAXON="d8d21e6686e3f069d4374d4be1c4e47744788f9244c46a1bbbbb9015e5cd7fb9"
PLANTS23="b436ba8668d67d05a247f9cafd59b846db2327031e59f4ad95c78d6870bf113e"
PLANTS15="44dc11ad3a0b73aef50e92ca28f49db3927f217de9496398f1b9ba9687a2bf15"
PUBLIC={
 BASE/"island_host_network_20261030/island_host_network_all225.csv":"39fe84981cab2b0e47063399f0a5c0c9731e800d",
 BASE/"historical_botanical_17island_20261031/near_pair_pre2011_plant_composition.json":"54c4515947d2a3ee7545e4007985ef82a2af6868",
 Path("data/frozen/figure_sources/host_contribution_metrics.csv"):"4cf4851499ca2233ff2adadf1c4ee7f0e863f23d",
 BASE/"island_host_continuity_20261030/island_host_continuity_events.csv":"efbc770469a720e2dea82478db2f6196ba03b353",
}
SEEDS=(20261031,20261032,20261033)
SAMPLES=499

def verify_private(path,sha):
 raw=path.read_bytes()
 if hashlib.sha256(raw).hexdigest()!=sha:
  raise ValueError("Private source SHA-256 mismatch: "+str(path))
 return raw

def public(root,path):
 raw=(root/path).read_bytes()
 actual=hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
 if actual!=PUBLIC[path]:raise ValueError("Source published Git blob mismatch "+str(path))
 return raw.decode("utf8")

def csv_rows(s):
 if isinstance(s,bytes):s=s.decode("utf8")
 return list(csv.DictReader(io.StringIO(s)))

def components(graph):
 ids={str(v["island_id"]) for v in graph["islands"]}
 if len(ids)!=17 or len(graph["source_pairs"])!=14:
  raise ValueError("Expected original 17 islands / 14 near island pairs")
 parent={i:i for i in ids}
 def find(i):
  if parent[i]!=i:parent[i]=find(parent[i])
  return parent[i]
 for p in graph["source_pairs"]:
  a,b=find(str(p["island_a"])),find(str(p["island_b"]))
  if a!=b:parent[b]=a
 comps=defaultdict(set)
 for i in ids:comps[find(i)].add(i)
 ordered=sorted([sorted(v,key=int) for v in comps.values()],key=lambda z:int(z[0]))
 if sorted(map(len,ordered))!=[2,2,2,2,4,5]:
  raise ValueError("Original 6-island regional graph changed")
 return ordered,{i:n for n,g in enumerate(ordered) for i in g}

def host_catalogue(edges,tax):
 accepted=defaultdict(set);hosts=defaultdict(set)
 for e in edges:
  accepted[e["current_accepted_name"]].add(e["original_name"])
  if e["query_eligible"]=="True" and e["query_name"]:
   k=tax[e["query_name"]]["gbif_speciesKey"]
   if k:hosts[e["original_name"]].add(str(k))
 if len(hosts)!=228:
  raise ValueError("Original species-level strict GBIF larval-host coverage changed")
 return hosts,accepted

def load_plants(p23,p15,ids):
 plants=defaultdict(set)
 for path,sha,key,nkey,count in [
   (p23,PLANTS23,"gbif_speciesKey","indexed_record_count",9788),
   (p15,PLANTS15,"speciesKey","indexed_records",327)]:
  raw=verify_private(path,sha)
  source=csv_rows(gzip.decompress(raw))
  if len(source)!=count:raise ValueError("Original plant GBIF species-by-year facets changed")
  for q in source:
   if str(q["island_id"]) not in ids or not(2001<=int(q["year"])<=2010):
    raise ValueError("Plant source is not pre-outcome or foreign island")
   if int(q[nkey])<=0:raise ValueError("Invalid historical plant indexed count")
   plants[str(q["island_id"])].add(str(q[key]))
 return plants

def read_original(path,local_pickle):
 if hashlib.sha256(path.read_bytes()).hexdigest()!=ORIGINAL:
  raise ValueError("Original uploaded butterfly-year Parquet SHA-256 changed")
 import pandas as pd
 frame=pd.read_pickle(local_pickle) if local_pickle else pd.read_parquet(path)
 if (len(frame),frame.speciesKey.nunique(),frame.landmass_id.nunique(),
     int(frame.records.sum()),int(frame.year.min()),int(frame.year.max()))!=(83573,224,3384,14425536,2001,2024):
  raise ValueError("Original precise butterfly species-island-year source changed")
 if frame.duplicated(["speciesKey","landmass_id","year"]).any() or not frame.records.eq(frame.precise1km_records).all():
  raise ValueError("Original precision or uniqueness changed")
 frame=frame.copy()
 frame["landmass_id"]=frame.landmass_id.astype(str)
 return frame

def make_risk(data,lookup,geb,hosts,alias,flora,by_area,units):
 first=data[data.landmass_id.isin(by_area)].groupby(["speciesKey","landmass_id"]).year.min().to_dict()
 candidate=sorted(lookup.keys())
 original_new=[]
 for species in candidate:
  k=lookup[species]
  for island in sorted(by_area,key=int):
   fy=first.get((k,island),9999)
   if 2011<=fy<=2017:original_new.append((species,island))
 if len(original_new)!=121 or len(set(a for a,b in original_new))!=64:
  raise ValueError("Original 2011–2017 first-record cohort changed")
 select=set(s for s,i in original_new)
 if sum(int(float(geb[s]["introduced_added_units"])>0) for s,i in original_new)!=120:
  raise ValueError("GEB 239-source opportunity gain across island events differs")
 rows=[];lost=[]
 for species in sorted(select):
  aliases={species} if species in hosts else alias.get(species,set()) & set(hosts)
  if len(aliases)!=1:
   lost.append(species);continue
  taxon=next(iter(aliases))
  for island in sorted(by_area,key=int):
   fy=int(first.get((lookup[species],island),9999))
   if fy<=2010:continue
   rows.append({"species":species,"island":island,"component":by_area[island],
     "compatible_host_indexed_2001_2010":int(bool(hosts[taxon]&flora[island])),
     "first_2011_2017":int(2011<=fy<=2017),
     "first_2018_2022":int(2018<=fy<=2022),
     "before2023_first_year":fy,
     "GEB_additional_WGSRPD3_resource_units":int(float(geb[species]["introduced_added_units"])),
     "GBIF_port_botanical_species_keys_indexed":len(flora[island])})
 import pandas as pd
 panel=pd.DataFrame(rows)
 if (len(panel),panel.species.nunique(),int(panel.first_2011_2017.sum()),int((panel.first_2011_2017*panel.compatible_host_indexed_2001_2010).sum()))!=(905,60,113,34):
  raise ValueError("Eligible original-year butterfly×island source cohort changed")
 if sum(s in lost for s,i in original_new)!=8:
  raise ValueError("Original 8 source taxonomy-unresolved first events changed")
 return panel,len(select),len(lost),original_new

def risk_rates(panel,outcome,regions):
 selected=panel.copy()
 if outcome=="first_2018_2022":
  selected=selected[selected.first_2011_2017==0].copy()
 selected["outcome"]=selected[outcome].astype(int)
 host=selected["compatible_host_indexed_2001_2010"]
 a=selected[host==1];b=selected[host==0]
 def block(d):
  return {"risk_cells":int(len(d)),"new_indexed_events":int(d.outcome.sum())}
 return selected,{
  "cells":int(len(selected)),"indexed_first_events":int(selected.outcome.sum()),
  "butterfly_species":int(selected.species.nunique()),
  "islands":int(selected.island.nunique()),
  "host_early_recorded":block(a),
  "host_early_NOT_indexed":block(b),
  "raw_rate_recorded_host":float(a.outcome.mean()),
  "raw_rate_no_host_index":float(b.outcome.mean()),
  "no_index_is_not_botanical_absence":True,
  "by_geographic_component":[
   {"islands":regions[k],"risk_cells":len(g),
    "first_index_events":int(g.outcome.sum()),
    "early_compatible_host_indexed_risk_cells":int(g.compatible_host_indexed_2001_2010.sum()),
    "first_index_events_at_preindexed_compatible_host":int((g.outcome*g.compatible_host_indexed_2001_2010).sum())}
   for k,g in selected.groupby("component")]}

def blocked_swaps(panel,seed):
 """Original at-risk cells, fixed species×component and island event margins."""
 rng=random.Random(seed)
 allowed={(x.species,x.island) for x in panel.itertuples()}
 occupied={(x.species,x.island) for x in panel.itertuples() if x.outcome}
 host={(x.species,x.island):x.compatible_host_indexed_2001_2010 for x in panel.itertuples()}
 groups=defaultdict(list)
 for x in panel.itertuples():
  if x.outcome:groups[x.component].append((x.species,x.island))
 keys=list(groups)
 original=sum(host[x] for x in occupied)
 matches=Counter(x[0] for x in occupied)
 islandmargins=Counter(x[1] for x in occupied)
 partmargins=Counter((x[0],next(z.component for z in panel.itertuples() if z.island==x[1])) for x in occupied)
 accepted=0
 def steps(N):
  nonlocal accepted
  for _ in range(N):
   group=rng.choice(keys);v=groups[group]
   i=rng.randrange(len(v));j=rng.randrange(len(v))
   if i==j:continue
   a,p=v[i];b,q=v[j]
   if a==b or p==q:continue
   if (a,q) not in allowed or (b,p) not in allowed:continue
   if (a,q) in occupied or (b,p) in occupied:continue
   occupied.remove((a,p));occupied.remove((b,q))
   occupied.add((a,q));occupied.add((b,p))
   v[i]=(a,q);v[j]=(b,p);accepted+=1
 steps(60*len(occupied))
 seq=[]
 for _ in range(SAMPLES):
  steps(10*len(occupied))
  seq.append(sum(host[k] for k in occupied))
 if Counter(x[0] for x in occupied)!=matches or Counter(x[1] for x in occupied)!=islandmargins:
  raise ValueError("At-risk null record margins violated")
 if Counter((sp,next(z.component for z in panel.itertuples() if z.island==iid)) for sp,iid in occupied)!=partmargins:
  raise ValueError("Butterfly×geography margins violated")
 return {"seed":seed,"draws":len(seq),
   "mean_matched_event_count":sum(seq)/len(seq),
   "min":min(seq),"max":max(seq),
   "interval95":[empirical_quantile(seq,.025),empirical_quantile(seq,.975)],
   "draws_at_least_observed":sum(v>=original for v in seq),
   "accepted_swaps":accepted}

def empirical_quantile(seq,q):
 s=sorted(seq)
 i=int(q*(len(s)-1))
 return float(s[i])

def null_summary(panel):
 original=int((panel.outcome*panel.compatible_host_indexed_2001_2010).sum())
 chain=[blocked_swaps(panel,seed) for seed in SEEDS]
 draws=sum(x["draws"] for x in chain)
 exceeds=sum(x["draws_at_least_observed"] for x in chain)
 return {"method":"within-geographic-component bipartite double-edge swaps, conditional on eligibility and fixed butterfly species x component / island first-record margins",
   "chains":chain,"total_draws":draws,"observed_matched_events":original,
   "null_mean_matched_events":sum(x["mean_matched_event_count"] for x in chain)/len(chain),
   "reference_tail_fraction":exceeds/draws,
   "reference_tail_draws":exceeds,
   "no_independent_mc_pvalue_or_causal_interpretation":True}

def execute(args):
 d=read_original(args.original_year_parquet,args.verified_local_decoded_pickle)
 rows=csv_rows(public(args.root,BASE/"island_host_network_20261030/island_host_network_all225.csv"))
 lookup={r["species"]:int(r["taxon_key"]) for r in rows}
 if len(lookup)!=225:raise ValueError("225 source butterflies changed")
 gebrows=csv_rows(public(args.root,Path("data/frozen/figure_sources/host_contribution_metrics.csv")))
 geb={r["species"]:r for r in gebrows}
 if len(geb)!=239:raise ValueError("GEB 239 primary source changed")
 geo=json.loads(public(args.root,BASE/"historical_botanical_17island_20261031/near_pair_pre2011_plant_composition.json"))
 components_list,by_area=components(geo)
 original_events=csv_rows(public(args.root,BASE/"island_host_continuity_20261030/island_host_continuity_events.csv"))
 if len(original_events)!=5665:raise ValueError("Historical source butterfly event census changed")
 edges=csv_rows(verify_private(args.author_original_host_edges,HOST_EDGES))
 if len(edges)!=2524:raise ValueError("Original butterfly host source changed")
 names=json.loads(gzip.decompress(verify_private(args.original_name_to_gbifkey_gzip,TAXON)))["name_resolution"]
 if len(names)!=1642:raise ValueError("Original plant name source changed")
 hosts,alias=host_catalogue(edges,names)
 flora=load_plants(args.plant23_gzip,args.plant15_gzip,set(by_area))
 for island in by_area:
  expected=next(z["n_indexed_species_keys"] for z in geo["islands"] if z["island_id"]==island)
  if len(flora[island])!=expected:
   raise ValueError("Pre-2011 botanical index species union differs from archive")
 study,nspecies,lost,original_new=make_risk(d,lookup,geb,hosts,alias,flora,by_area,components_list)
 p1,desc1=risk_rates(study,"first_2011_2017",components_list)
 p2,desc2=risk_rates(study,"first_2018_2022",components_list)
 output={
 "schema":"chocho_geb_port_island_species_at_risk_2011_2022_v1",
 "status":"SOURCE_VERIFIED_EXPLORATORY_TWO_WINDOWS_NO_CAUSAL_GATEWAY_EVIDENCE",
 "original_private_species_island_year_parquet_SHA256":ORIGINAL,
 "private_author_original_230_edges_SHA256":HOST_EDGES,
 "frozen_name_to_gbif_key_gzip_SHA256":TAXON,
 "original_23port_historical_plant_GBIF_gzip_SHA256":PLANTS23,
 "additional_15port_historical_plant_GBIF_gzip_SHA256":PLANTS15,
 "source_GEB_contribution_metrics_git_blob":PUBLIC[Path("data/frozen/figure_sources/host_contribution_metrics.csv")],
 "source_225_butterfly_name_GBIF_key_git_blob":PUBLIC[BASE/"island_host_network_20261030/island_host_network_all225.csv"],
 "geographic_components":components_list,
 "all_island_first_2011_2017_before_original230_name_screen":{
   "events":len(original_new),"species":nspecies,
   "GEB_species_with_positive_global_resource_expansion":sum(float(geb[s]["introduced_added_units"])>0 for s in set(s for s,i in original_new)),
   "GEB_events_for_species_with_positive_global_resource_expansion":sum(float(geb[s]["introduced_added_units"])>0 for s,i in original_new)},
 "selected_name_join":{"source_event_species":nspecies,"source_events":len(original_new),
   "resolved_larval_host_species":study.species.nunique(),"resolved_event_count":int(study.first_2011_2017.sum()),
   "unresolved_source_events":len(original_new)-int(study.first_2011_2017.sum()),
   "selected_species_required_at_least_one_2011_2017_original_first_index_in_study_islands":True},
 "early_host_index_year_window":"2001–2010",
 "early_host_plant_status":"exact GBIF accepted species key, indexed within one of the 38 2019 WPI 5km port buffers",
 "periods":{
   "2011_2017_first_index_among_no_pre2011_butterfly_index":{"recording_risk":desc1,
       "geographic_margin_null":null_summary(p1)},
   "2018_2022_first_index_among_still_unindexed_after_2017":{"recording_risk":desc2,
       "geographic_margin_null":null_summary(p2)}},
 "study_limits":{
   "only_60_butterfly_species_preselected_by_having_at_least_one_2011_2017_first_index":True,
   "historical_botanical_nonindex_is_NOT_plant_absence":True,
   "GBIF_first_record_not_arrival_or_establishment":True,
   "five_year_butterfly_recording_not_larval_use_or_successful_breeding":True,
   "GEB_resource_opportunity_is_modern_regional_not_historical_island_specific":True,
   "no_port_vs_nonport_flora_control":True,
   "no_actual_island_specific_plant_alien_status_or_shipping_mode":True,
   "null_MCMC_draws_correlated_and_exchangeability_unproven":True,
   "no_proof_of_positive_or_negative_butterfly_colonization_effect":True,
   "private_source_rows_and_links_not_published":True},
 "scientific_readout":"Raw preindexed compatible host presence tracks later butterfly indexing, particularly 2018-2022, but after preserving butterfly x regional component first-record margins, island first-record margins, and risk-set prior-history eligibility, the observed host-index overlap is nearly identical to the conditional expected overlap in both windows. No independent port realization of the +54.9% GEB opportunity is identified."
 }
 args.out.parent.mkdir(parents=True,exist_ok=True)
 args.out.write_text(json.dumps(output,sort_keys=True,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
 print("GEB_PORT_ISLAND_RISK",json.dumps({key:{"n":v["recording_risk"]["cells"],"events":v["recording_risk"]["indexed_first_events"],"observed":v["geographic_margin_null"]["observed_matched_events"],"expected":v["geographic_margin_null"]["null_mean_matched_events"],"ge":v["geographic_margin_null"]["reference_tail_draws"]} for key,v in output["periods"].items()}))
if __name__=="__main__":
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument("--root",type=Path,default=ROOT)
 p.add_argument("--original-year-parquet",type=Path,required=True)
 p.add_argument("--author-original-host-edges",type=Path,required=True)
 p.add_argument("--original-name-to-gbifkey-gzip",type=Path,required=True)
 p.add_argument("--plant23-gzip",type=Path,required=True)
 p.add_argument("--plant15-gzip",type=Path,required=True)
 p.add_argument("--verified-local-decoded-pickle",type=Path,default=None,
   help="Optional private local decoded Parquet verified by its actual original-source SHA")
 p.add_argument("--out",type=Path,required=True)
 execute(p.parse_args())
