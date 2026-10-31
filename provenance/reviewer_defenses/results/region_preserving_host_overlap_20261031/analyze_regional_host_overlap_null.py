#!/usr/bin/env python3
"""Swap observed butterfly×island first-indexing edges, preserving source geography.

Null contrasts:
 (A) global fixed butterfly event island degrees (crosses continents; diagnostic),
 (B) fixed species × *near-island connected component* and island event margins.
Plant speciesKey-by-island history and ORIGINAL butterfly larval-host identities
are kept FIXED. No ship routes, local larval survival or true host absences inferred.

This is post-hoc macro-ecological observation testing; MCMC draws are correlated
and do not supply a preregistered causal p-value.
"""
from __future__ import annotations
import argparse
from collections import Counter,defaultdict
import csv,gzip,hashlib,io,json,random,math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[4]
BASE=Path("provenance/reviewer_defenses/results")
BUTTERFLIES=BASE/"island_host_continuity_20261030/island_host_continuity_events.csv"
GEOGRAPHY=BASE/"historical_botanical_17island_20261031/near_pair_pre2011_plant_composition.json"
CANDIDATES=BASE/"island_host_network_20261030/island_host_network_all225.csv"
SOURCE_SHA={
  BUTTERFLIES:"efbc770469a720e2dea82478db2f6196ba03b353",
  GEOGRAPHY:"54c4515947d2a3ee7545e4007985ef82a2af6868",
  CANDIDATES:"39fe84981cab2b0e47063399f0a5c0c9731e800d",
}
HOST_SHA="e55012d6b0eefd2b3e9e4b8f13fc4786ddc1e66c21f73e705bb20df92d0b5e22"
TAXON_SHA="d8d21e6686e3f069d4374d4be1c4e47744788f9244c46a1bbbbb9015e5cd7fb9"
PLANTS23_SHA="b436ba8668d67d05a247f9cafd59b846db2327031e59f4ad95c78d6870bf113e"
DRAWS_PER_CHAIN=499
CHAINS=(20261031,20261032,20261033)
BURN_EDGE_MULTIPLIER=60
THIN_EDGE_MULTIPLIER=10

def strict_file(p,sha):
 b=p.read_bytes()
 if hashlib.sha256(b).hexdigest()!=sha:
  raise ValueError("Private source SHA-256 mismatch: "+str(p))
 return b

def git_file(root,p):
 b=(root/p).read_bytes()
 sha=hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()
 if sha!=SOURCE_SHA[p]:
  raise ValueError("Source Git blob differs: "+str(p))
 return b.decode("utf8")

def read_csv(data):
 if isinstance(data,bytes):data=data.decode("utf8")
 return list(csv.DictReader(io.StringIO(data)))

def source_data(root,hosts_path,plants23_path,plants15_path,taxonomy_path):
 geo=json.loads(git_file(root,GEOGRAPHY))
 events=read_csv(git_file(root,BUTTERFLIES))
 candidates=read_csv(git_file(root,CANDIDATES))
 if (len(events),len(candidates))!=(5665,225):
  raise ValueError("Original source candidate butterfly/event census changed")
 butterfly_names={r["species"] for r in candidates}
 names={x["island_id"] for x in geo["islands"]}
 if len(names)!=17 or len(geo["source_pairs"])!=14:
  raise ValueError("Near-island original geography changed")
 orig=read_csv(strict_file(hosts_path,HOST_SHA))
 if len(orig)!=2524:
  raise ValueError("Original author 2,524 larval host edges changed")
 taxon_blob=strict_file(taxonomy_path,TAXON_SHA)
 taxon=json.loads(gzip.decompress(taxon_blob))
 if taxon.get("schema")!="chocho_accepted_original_host_species_gbif_crosswalk_v1":
  raise ValueError("Author host accepted-name crosswalk source schema mismatch")
 mapping=taxon["name_resolution"]
 if len(mapping)!=1642:
  raise ValueError("Expected exact 1642 host plant source name resolutions")
 # Protect original source accepted plant name mapping and avoid ambiguity.
 host=defaultdict(set)
 byaccepted=defaultdict(set)
 for h in orig:
  current=h["current_accepted_name"];original=h["original_name"]
  byaccepted[current].add(original)
  if h["query_eligible"]!="True" or not h["query_name"]:
   continue
  key=mapping[h["query_name"]]["gbif_speciesKey"]
  if key:host[original].add(str(key))
 if not 220<=len(host)<=230:
  raise ValueError("Original butterfly larval host resolution changed")
 plant_rows=[]
 for p,sha,old in [(plants23_path,PLANTS23_SHA,True),
                    (plants15_path,geo["new_15port_gzip_sha256"],False)]:
  raw=strict_file(p,sha)
  rec=read_csv(gzip.decompress(raw).decode("utf8"))
  if len(rec)!=(9788 if old else 327):
   raise ValueError("Archived GBIF plant historical facet row count differs")
  for r in rec:
   iid=r["island_id"];yr=int(r["year"])
   if iid not in names or not 2001<=yr<=2010:
    raise ValueError("Botanical record outside declared islands and years")
   key=(r["gbif_speciesKey"] if old else r["speciesKey"])
   plant_rows.append((iid,str(key)))
 flora=defaultdict(set)
 for iid,key in plant_rows:flora[iid].add(key)
 if any(len(flora[s])!=n["n_indexed_species_keys"] for n in geo["islands"] for s in [n["island_id"]]):
  raise ValueError("GBIF original-island historic flora richness mismatch")
 selected=[]
 for e in events:
  iid=e["landmass_id"]
  if iid not in names or not 2011<=int(e["first_record_year"])<=2017:
   continue
  sp=e["species"]
  if sp not in butterfly_names:
   raise ValueError("Butterfly not in fixed source species roster")
  if sp in host:
   source=sp
  else:
   eligible=byaccepted.get(sp,set()) & set(host)
   source=next(iter(eligible)) if len(eligible)==1 else None
  if source is not None:
   selected.append((source,iid))
 if len(selected)!=113 or len(set(selected))!=113 or len({sp for sp,i in selected})!=60:
  raise ValueError("Unselected 113 butterfly x island first-record cohort mismatch")
 observed=sum(bool(host[sp]&flora[iid]) for sp,iid in selected)
 if observed!=34:raise ValueError("Observed prior indexed larval hosts 34/113 changed")
 return selected,host,flora,geo,observed

def components(geo):
 nodes={n["island_id"] for n in geo["islands"]}
 parent={n:n for n in nodes}
 def find(n):
  if parent[n]!=n:parent[n]=find(parent[n])
  return parent[n]
 for p in geo["source_pairs"]:
  a,b=find(p["island_a"]),find(p["island_b"])
  if a!=b:parent[b]=a
 comp={i:find(i) for i in nodes}
 members=defaultdict(set)
 for k,v in comp.items():members[v].add(k)
 if len(members)!=6 or sorted(len(g) for g in members.values())!=[2,2,2,2,4,5]:
  raise ValueError("Unexpected geographic connectivity of 17 island nodes")
 return comp,sorted([sorted(v,key=int) for v in members.values()],key=lambda a:int(a[0]))

def empirical_quantile(x,q):
 a=sorted(x)
 return a[int(q*(len(a)-1))] if a else None

def spearman_chain_lag1(xs):
 if len(xs)<4:return None
 n=len(xs)-1
 a=xs[:-1];b=xs[1:]
 ma=sum(a)/n;mb=sum(b)/n
 z=sum((x-ma)*(y-mb) for x,y in zip(a,b))
 den=math.sqrt(sum((x-ma)**2 for x in a)*sum((y-mb)**2 for y in b))
 return z/den if den>0 else None

def make_markov(edges,plant,flora,groups,seed,observed):
 rng=random.Random(seed)
 a=[(s,i) for s,i in sorted(edges,key=lambda z:(z[1],z[0]))]
 base=set(a)
 occupied=set(a)
 by_component=[list(g) for g in groups if len(g)>=2]
 before_species=Counter(s for s,i in a)
 before_island=Counter(i for s,i in a)
 before_comp_species=Counter((s,z) for z,g in enumerate(groups) for k in g for s,i in [a[k]])
 changes=0
 attempts=0
 def swap(k):
  nonlocal changes,attempts
  for _ in range(k):
   attempts+=1
   g=rng.choice(by_component)
   j1=rng.choice(g);j2=rng.choice(g)
   if j1==j2:continue
   s1,i1=a[j1];s2,i2=a[j2]
   if s1==s2 or i1==i2 or (s1,i2) in occupied or (s2,i1) in occupied:
    continue
   occupied.remove((s1,i1));occupied.remove((s2,i2))
   occupied.add((s1,i2));occupied.add((s2,i1))
   a[j1]=(s1,i2);a[j2]=(s2,i1)
   changes+=1
 swap(BURN_EDGE_MULTIPLIER*len(edges))
 score=[]
 retention=[]
 for t in range(DRAWS_PER_CHAIN):
  swap(THIN_EDGE_MULTIPLIER*len(edges))
  score.append(sum(bool(plant[s]&flora[i]) for s,i in a))
  retention.append(len(base&occupied)/len(base))
 after_species=Counter(s for s,i in a)
 after_island=Counter(i for s,i in a)
 after_comp_species=Counter((s,z) for z,g in enumerate(groups) for k in g for s,i in [a[k]])
 if before_species!=after_species or before_island!=after_island or before_comp_species!=after_comp_species:
  raise ValueError("Observed butterfly×island or regional species degrees drifted")
 if len(occupied)!=len(edges):
  raise ValueError("Null generated duplicate butterfly island source pairs")
 return {
  "seed":seed,"draws":len(score),"accepted_double_edge_swaps":changes,
  "swap_attempts":attempts,"source_pair_retention_mean":sum(retention)/len(retention),
  "scores":score,"mean":sum(score)/len(score),
  "lag1_score_autocorrelation":spearman_chain_lag1(score),
  "draws_ge_observed":sum(x>=observed for x in score),
  "range":[min(score),max(score)],
 }

def run(args):
 e,hosts,flora,geo,observed=source_data(args.root,args.original_roster,args.original23,args.new15,args.taxonomy)
 comp,members=components(geo)
 # Selected events may be absent in some geographic components, but fixed graph controls remain.
 global_groups=[list(range(len(e)))]
 regional_groups=[sorted([k for k,v in enumerate(e) if comp[v[1]]==comp[g[0]]])
  for g in members]
 regional_groups=[z for z in regional_groups if z]
 if sum(len(g) for g in regional_groups)!=len(e):
  raise ValueError("Regional groups don't partition all eligible first-record observations")
 out={"schema":"chocho_17island_species_island_fixed_margin_geography_null_v1",
  "status":"POSTHOC_SOURCE_GEOGRAPHY_BLOCKED_RECORD_MATCH_NULL_NOT_CAUSAL",
  "source_author_edges_sha256":HOST_SHA,"source_taxonomy_crosswalk_sha256":TAXON_SHA,
  "source_public_git_blobs":{str(k):v for k,v in SOURCE_SHA.items()},
  "original23_plant_gzip_sha256":PLANTS23_SHA,
  "new15_plant_gzip_sha256":geo["new_15port_gzip_sha256"],
  "observation_cohort_n":len(e),"butterfly_species_count":len({s for s,i in e}),
  "island_count":len({i for s,i in e}),
  "observed_indexed_compatible_host_n":observed,
  "observed_indexed_compatible_host_fraction":observed/len(e),
  "geographic_near_network_components":members,
  "preservation":{
    "row_species_event_counts":True,
    "column_island_event_counts":True,
    "species_x_near_component_event_counts_for_regional_null":True,
    "original_plant_species_island_presence_and_host_identity_fixed":True},
  "settings":{
    "draws_per_chain":DRAWS_PER_CHAIN,"seeds":list(CHAINS),
    "burn_attempts_per_observed_event":BURN_EDGE_MULTIPLIER,
    "attempts_between_retained_draws_per_event":THIN_EDGE_MULTIPLIER
  },
  "nulls":{},
  "inference_limits":{
    "fixed_2011_2017_indexed_first_recording_not_true_colonization":True,
    "host_plant_indexed_pre2011_not_confirmed_larval_feeding":True,
    "source_name_taxonomy_acceptance_not_historical_island_introduction_status":True,
    "no_host_indexed_not_true_absence":True,
    "within_geographic_connected_group_swaps_do_not_control_climate_or_habitat":True,
    "global_cross_region_null_is_geographically_unrealistic":True,
    "mcmc_draws_not_independent_samples":True,
    "not_preregistered":True,
    "actual_ship_routes_and_source_populations_unknown":True
  }}
 for label,groups in [("global_species_island_fixed_margins",global_groups),
                     ("within_six_geographic_components_fixed_margins",regional_groups)]:
  result=[make_markov(e,hosts,flora,groups,seed,observed) for seed in CHAINS]
  pooled=[v for item in result for v in item["scores"]]
  out["nulls"][label]={
    "draws":len(pooled),"mean_overlap_n":sum(pooled)/len(pooled),
    "mean_overlap_fraction":sum(pooled)/(len(pooled)*len(e)),
    "interval95_overlap_n":[empirical_quantile(pooled,.025),empirical_quantile(pooled,.975)],
    "inclusive_tail_ge_observed":sum(v>=observed for v in pooled),
    "inclusive_fraction":sum(v>=observed for v in pooled)/len(pooled),
    "observed_minus_null_mean_count":observed-sum(pooled)/len(pooled),
    "chains":[{k:v for k,v in result_item.items() if k!="scores"} for result_item in result],
    "samples_are_autocorrelated_markov_chain_draws":True
  }
  print("SOURCE_GEOGRAPHY_NULL",label,json.dumps({
    k:v for k,v in out["nulls"][label].items() if k not in ("chains",)
  }),flush=True)
 args.out.parent.mkdir(parents=True,exist_ok=True)
 args.out.write_text(json.dumps(out,ensure_ascii=False,sort_keys=True,indent=2)+"\n",encoding="utf8")
 print("GEO_NULL_SOURCE_COMPLETE",json.dumps({
    "events":len(e),"observed":observed,"groups":len(members),
    "regional_support":[len(group) for group in regional_groups]},sort_keys=True),flush=True)

if __name__=="__main__":
 p=argparse.ArgumentParser()
 p.add_argument("--root",type=Path,default=Path("."))
 p.add_argument("--original-roster",type=Path,required=True)
 p.add_argument("--original23",type=Path,required=True)
 p.add_argument("--new15",type=Path,required=True)
 p.add_argument("--taxonomy",type=Path,required=True)
 p.add_argument("--out",type=Path,required=True)
 run(p.parse_args())
