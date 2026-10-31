#!/usr/bin/env python3
"""Test historical compatible host co-indexing under plant-geography-conditioned nulls.

Uses only SHA-pinned ORIGINAL author-side host edges, taxonomic name crosswalk,
17-island GBIF botanical speciesKey/year facets, and source-indexed butterflies.
All individual host edges and raw occurrence coordinates remain UNPUBLISHED.
The key 34/113 is detection co-indexing, not host-aided breeding or colonization.
"""
from __future__ import annotations
import argparse,csv,gzip,hashlib,importlib.util,io,json,math,random
from collections import Counter,defaultdict
from pathlib import Path
ROOT=Path(__file__).resolve().parents[4]
M=Path("provenance/reviewer_defenses/results/original230_dated_host_join_20261031/run_dated_host_recurrence.py")
M_BLOB="0688e362d80f23d83d8d5d25ba9cefb24a705167"
FROZEN=Path("provenance/reviewer_defenses/results/original230_dated_host_join_20261031/original230_prior_host_recurrence_audit.json")
FROZEN_BLOB="0272020aad8967bcd1a8eb8a6ca3995e33a15ccc"
SOURCE_SCHEMA="chocho_source_original230_dated_plant_host_fiveyear_butterfly_recording_v1"
TAX_SHA="d8d21e6686e3f069d4374d4be1c4e47744788f9244c46a1bbbbb9015e5cd7fb9"
N_DRAWS=499
SEED=20261031
FREQ_BIN_NAMES=("0","1","2to3","4plus")

def gitblob(path):
 b=path.read_bytes()
 return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def load_original(root):
 if gitblob(root/M)!=M_BLOB or gitblob(root/FROZEN)!=FROZEN_BLOB:
  raise ValueError("Original fixed source algorithm or verified result changed")
 spec=importlib.util.spec_from_file_location("original230_frozen",root/M)
 module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
 old=json.loads((root/FROZEN).read_text(encoding="utf8"))
 if old["schema"]!=SOURCE_SCHEMA:raise ValueError("Original source result schema changed")
 return module,old

def freq_bin(n):
 if n==0:return "0"
 if n==1:return "1"
 if n<=3:return "2to3"
 return "4plus"

def read_sources(m,args):
 if sha(args.author_roster_dir/"original230_host_edges.csv")!=m.HOSTS_SHA:
  raise ValueError("Original author edge source SHA changed")
 if sha(args.author_roster_dir/"plant_taxon_query_candidates.csv")!=m.PLANT_CANDIDATES_SHA:
  raise ValueError("Original plant taxonomic candidate identity changed")
 if sha(args.taxonomy_gzip)!=TAX_SHA:
  raise ValueError("GBIF accepted plant-name-to-speciesKey crosswalk changed")
 payload=json.loads(gzip.decompress(args.taxonomy_gzip.read_bytes()))
 if payload["source_author_plant_candidates_sha256"]!=m.PLANT_CANDIDATES_SHA:
  raise ValueError("Plant name mapping source mismatch")
 names=payload["name_resolution"]
 if len(names)!=1642:raise ValueError("Incomplete frozen exact accepted name matches")
 rows=m.readcsv((args.author_roster_dir/"original230_host_edges.csv").read_bytes())
 if len(rows)!=2524:raise ValueError("Author host roster count changed")
 edges={(z["original_name"],names[z["query_name"]]["gbif_speciesKey"],z["family"])
   for z in rows if z["query_eligible"]=="True" and z["query_name"]
      and names[z["query_name"]]["gbif_speciesKey"] is not None}
 if len(edges)!=2396:raise ValueError("Original resolved 2396 host edges changed")
 source=json.loads(m.source(args.root,m.FLORA))
 historical=(m.read_source_plant_gzip(args.original23,m.ORIGINAL23_SHA,True)
     +m.read_source_plant_gzip(args.new15,source["new_15port_gzip_sha256"],False))
 plantsets=defaultdict(set)
 for island,k,year,port,count in historical:
  plantsets[island].add(k)
 for z in source["islands"]:
  if len(plantsets[z["island_id"]])!=z["n_indexed_species_keys"]:
   raise ValueError("Island source flora union count changed")
 raw_events=m.readcsv(m.source(args.root,m.EVENTS))
 study=[z for z in raw_events if z["landmass_id"] in plantsets
          and 2011<=int(z["first_record_year"])<=2017]
 if len(study)!=121:raise ValueError("Butterfly event source cohort changed")
 aliases=defaultdict(set)
 for e in rows:aliases[e["current_accepted_name"]].add(e["original_name"])
 butterflies={b for b,k,f in edges}
 events=[]
 for e in study:
  name=e["species"]
  if name in butterflies:src=name
  elif len(aliases[name]&butterflies)==1:src=next(iter(aliases[name]&butterflies))
  else:continue
  events.append({"original_butterfly_name":src,"landmass_id":e["landmass_id"],
     "any_re_record":e["any_followup_record"]=="True",
     "record_year_five":e["record_in_year5"]=="True",
     "_island_plants":plantsets[e["landmass_id"]]})
 if len(events)!=113:raise ValueError("Source event identity mismatch")
 def match(ed):
  by=defaultdict(set)
  for s,k,f in ed:by[s].add(k)
  return [{**e,"has_prior_recorded_compatible_host":
      bool(by[e["original_butterfly_name"]]&e["_island_plants"])} for e in events]
 target=match(edges)
 if sum(e["has_prior_recorded_compatible_host"] for e in target)!=34:
  raise ValueError("Original original230 34/113 source co-indexing failed")
 return m,source,plantsets,edges,events,target,match

def graph_components(source):
 parent={}
 def find(a):
  if a not in parent:parent[a]=a
  if parent[a]!=a:parent[a]=find(parent[a])
  return parent[a]
 for pair in source["source_pairs"]:
  a,b=pair["island_a"],pair["island_b"]
  if find(a)!=find(b):parent[find(b)]=find(a)
 groups=defaultdict(set)
 for i in sorted(parent,key=int):groups[find(i)].add(i)
 sorted_groups=sorted(groups.values(),key=lambda v:min(map(int,v)))
 if (len(sorted_groups),sum(map(len,sorted_groups)))!=(6,17):
  raise ValueError("Fixed 14-near-pair graph disconnected differently")
 cluster_of={i:n for n,group in enumerate(sorted_groups) for i in group}
 return sorted_groups,cluster_of

def run_null(m,edges,events,plantsets,cluster_of,mode):
 """Fixed ORIGINAL degrees and butterfly×plant-family host counts plus geography bands."""
 all_host_keys={p for _,p,_ in edges}
 prevalence={p:len({i for i in plantsets if p in plantsets[i]}) for p in all_host_keys}
 regions={p:tuple(sorted({cluster_of[i] for i in plantsets if p in plantsets[i]})) for p in all_host_keys}
 if mode not in ("family_plus_island_prevalence_band","family_plus_region_signature_and_prevalence_band"):
  raise ValueError("Unknown fixed null")
 ed=[list(x) for x in sorted(edges)]
 orig=Counter((s,f) for s,p,f in ed)
 butterfly_degree=Counter(s for s,p,f in ed)
 plant_degree=Counter(p for s,p,f in ed)
 groups=defaultdict(list)
 for idx,(_,p,f) in enumerate(ed):
  grp=(f,freq_bin(prevalence[p]))
  if mode.endswith("_and_prevalence_band"):
   grp=(f,regions[p],freq_bin(prevalence[p]))
  groups[grp].append(idx)
 eligible=[v for v in groups.values()
   if len(v)>1 and len({ed[i][0] for i in v})>1 and len({ed[i][1] for i in v})>1]
 rng=random.Random(SEED+{"family_plus_island_prevalence_band":1,
                         "family_plus_region_signature_and_prevalence_band":2}[mode])
 used={(s,p) for s,p,f in ed}
 accepted=0
 def swap(n):
  nonlocal accepted
  for _ in range(n):
   if not eligible:return
   pool=rng.choice(eligible)
   i=rng.choice(pool);j=rng.choice(pool)
   a,x,f=ed[i];b,y,_=ed[j]
   if a==b or x==y or (a,y) in used or (b,x) in used:continue
   used.remove((a,x));used.remove((b,y))
   used.add((a,y));used.add((b,x))
   ed[i][1]=y;ed[j][1]=x
   accepted+=1
 def score():
  by=defaultdict(set)
  for s,p,f in ed:by[s].add(p)
  obs=[bool(by[e["original_butterfly_name"]]&e["_island_plants"]) for e in events]
  matched=sum(obs)
  a=sum(e["any_re_record"] for e,b in zip(events,obs) if b)
  b=sum(e["any_re_record"] for e,t in zip(events,obs) if not t)
  diff=a/matched-b/(len(events)-matched) if 0<matched<len(events) else None
  return matched/len(events),diff
 swap(25*len(ed))
 arr=[];deltas=[]
 for _ in range(N_DRAWS):
  swap(4*len(ed))
  f,d=score()
  arr.append(f)
  if d is not None:deltas.append(d)
 if (Counter(s for s,p,f in ed)!=butterfly_degree
     or Counter(p for s,p,f in ed)!=plant_degree
     or Counter((s,f) for s,p,f in ed)!=orig
     or len(used)!=len(ed)):
  raise ValueError("Null broke exact source host degrees/family")
 for s,p,f in ed:
  # Individual host plant degrees are preserved and each host always appears in same stratum.
  if not isinstance(prevalence[p],int):raise ValueError("Plant prevalence invalid")
 def quant(v,p):
  x=sorted(v)
  return x[int(p*(len(x)-1))]
 obs=34/113
 actual_delta=26/34-70/79
 return {
  "null_mode":mode,"draws":N_DRAWS,"random_seed":SEED+(1 if mode.endswith("_prevalence_band") and "region" not in mode else 2),
  "edges":len(ed),"swaps_accepted":accepted,
  "stratum_count":len(groups),"strata_with_swappable_hosts":len(eligible),
  "edges_in_swappable_strata":sum(map(len,eligible)),
  "plant_presence_pre2011_17islands_preserved_by_band":True,
  "plant_geographic_region_signature_exact_preserved":"region_signature" in mode,
  "each_plant_degree_and_each_butterfly_plant_family_count_preserved":True,
  "observed_prior_host_index_fraction":obs,
  "null_mean_fraction":sum(arr)/len(arr),
  "null_fraction_q025_q975":[quant(arr,.025),quant(arr,.975)],
  "null_draws_ge_observed_fraction":sum(v>=obs for v in arr),
  "null_lag1_fraction_autocorrelation":(
    sum((a-sum(arr)/len(arr))*(b-sum(arr)/len(arr)) for a,b in zip(arr[:-1],arr[1:]))
    /max(1e-30,sum((a-sum(arr)/len(arr))**2 for a in arr))
  ),
  "observed_reindex_delta":actual_delta,
  "null_rate_delta_mean":sum(deltas)/len(deltas) if deltas else None,
  "null_draws_ge_observed_reindex_delta":sum(v>=actual_delta for v in deltas),
  "valid_rate_delta_draws":len(deltas),
  "not_actual_plant_absence_or_larval_host_use":True,
  "not_actual_colonization_shipping_or_independent_species":True,
  "no_climate_or_observer_effort_control":True
 }

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument("--root",type=Path,default=Path("."))
 p.add_argument("--author-roster-dir",type=Path,required=True)
 p.add_argument("--original23",type=Path,required=True)
 p.add_argument("--new15",type=Path,required=True)
 p.add_argument("--taxonomy-gzip",type=Path,required=True)
 p.add_argument("--output",type=Path,required=True)
 args=p.parse_args()
 m,source,plantsets,edges,events,matched,match=read_sources(*((load_original(args.root)[0],args),))
 # IMPORTANT: redo original-source baseline to prove the 10/299 initial result matches.
 original=m.source_network_null(matched,edges)
 if original["null_draws_ge_observed_fraction"]!=10 or original["null_draws_ge_observed_delta"]!=194:
  raise ValueError("Original structural host null did not replay exactly")
 groups,cluster_of=graph_components(source)
 profiles={}
 for name in ("family_plus_island_prevalence_band","family_plus_region_signature_and_prevalence_band"):
  profiles[name]=run_null(m,edges,events,plantsets,cluster_of,name)
 output={
   "schema":"chocho_original230_geo_prevalence_host_swap_null_v1",
   "status":"EXPLORATORY_SOURCE_PINNED_HOST_INDEX_CO_OCCURRENCE_NOT_CAUSAL",
   "source_original_algo_git_blob_sha":M_BLOB,
   "source_original_frozen_audit_git_blob_sha":FROZEN_BLOB,
   "original_author_host_edges_sha256":m.HOSTS_SHA,
   "original_author_plant_candidates_sha256":m.PLANT_CANDIDATES_SHA,
   "taxon_mapping_gzip_sha256":TAX_SHA,
   "plant23_gzip_sha256":m.ORIGINAL23_SHA,
   "plant15_gzip_sha256":source["new_15port_gzip_sha256"],
   "source_species_island_events":113,
   "source_butterfly_species":len({x["original_butterfly_name"] for x in events}),
   "source_geographic_islands":17,
   "source_geographic_near_island_pairs":14,
   "host_source_original_edges":2524,
   "resolved_butterfly_host_edges":len(edges),
   "observed_compatible_preindexed_host_events":34,
   "observed_any_five_year_reindex_host_yes":26,
   "observed_any_five_year_reindex_host_no_index":70,
   "source_species_and_geographic_regions":[
      {"region_id":i,"islands":sorted(g,key=int)} for i,g in enumerate(groups)],
   "prior_family_and_degree_only_original_null_replayed":{
       "n_draws":original["n_draws"],"mean":original["null_fraction_mean"],
       "draws_ge_observed":original["null_draws_ge_observed_fraction"],
       "draws_ge_observed_recorder_delta":original["null_draws_ge_observed_delta"]
   },
   "more_constrained_structural_nulls":profiles,
   "original_host_non_detection_never_treated_as_host_absence":True,
   "plant_appearance_not_historical_intro_or_vessel_route":True,
   "not_tested_direct_mainland_source_or_butterfly_breeding":True,
   "none_of_modes_preserves_original_plot_level_botanical_recorders_or_yearly_weather":True,
   "original_host_edges_and_original_GBIF_points_unpublished":True
 }
 args.output.parent.mkdir(parents=True,exist_ok=True)
 args.output.write_text(json.dumps(output,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf8")
 print("GEO_PREVALENCE_HOST_NULL_AUDIT",json.dumps({
    "observed":34/113,
    "original":output["prior_family_and_degree_only_original_null_replayed"],
    "new":{k:{"mean":v["null_mean_fraction"],"ge":v["null_draws_ge_observed_fraction"],
       "accepted":v["swaps_accepted"]} for k,v in profiles.items()}
 },sort_keys=True))

if __name__=="__main__":main()
