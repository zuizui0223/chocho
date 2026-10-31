#!/usr/bin/env python3
"""Reconcile GEB 239-species resource gain to 225-species port and 17-island recording support.

All inputs are committed and Git blob SHA locked. Do not treat WGSRPD3
resource envelope or later port host co-indexing as measured species
colonisation, historical vessel-mediated introduction, or larval feeding.
"""
from __future__ import annotations
from pathlib import Path
import csv,hashlib,io,json,math
from collections import defaultdict

ROOT=Path(__file__).resolve().parents[4]
OUT=Path("provenance/reviewer_defenses/results/geb_port_island_realization_20261031/GEB_to_port_species_bridge.json")
PATHS={
 "data/frozen/figure_sources/host_contribution_metrics.csv":"4cf4851499ca2233ff2adadf1c4ee7f0e863f23d",
 "provenance/reviewer_defenses/results/port_connectivity_hosts_20261030/port_connectivity_host_all225.csv":"239dc82992fa3374e4862a1e2202663e1aeef7da",
 "provenance/reviewer_defenses/results/island_host_continuity_20261030/island_host_continuity_events.csv":"efbc770469a720e2dea82478db2f6196ba03b353",
 "provenance/reviewer_defenses/results/historical_botanical_17island_20261031/near_pair_pre2011_plant_composition.json":"54c4515947d2a3ee7545e4007985ef82a2af6868",
 "provenance/reviewer_defenses/results/island_host_network_20261030/island_host_network_all225.csv":"39fe84981cab2b0e47063399f0a5c0c9731e800d",
 "provenance/reviewer_defenses/results/butterfly_host_contribution_concentration_v0.1.json":"201a762de94b9ecf42a766c8cb5d8d1739a2ca3c",
}
def load(root):
 d={}
 for path,expected in PATHS.items():
  raw=(root/path).read_bytes()
  digest=hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
  if digest!=expected:raise ValueError("Git source blob identity has changed: "+path)
  d[path]=(list(csv.DictReader(io.StringIO(raw.decode("utf8"))))
    if path.endswith(".csv") else json.loads(raw.decode("utf8")))
 return d
def rank(a):
 order=sorted(range(len(a)),key=lambda i:a[i])
 z=[0.]*len(a);i=0
 while i<len(a):
  j=i+1
  while j<len(a) and a[order[j]]==a[order[i]]:j+=1
  for k in order[i:j]:z[k]=(i+j+1)/2
  i=j
 return z
def corr(a,b):
 n=len(a);am=sum(a)/n;bm=sum(b)/n
 x=[z-am for z in a];y=[z-bm for z in b]
 den=math.sqrt(sum(v*v for v in x)*sum(v*v for v in y))
 return sum(u*v for u,v in zip(x,y))/den if den else None
def solve(matrix,vector):
 n=len(vector);a=[[float(v) for v in r]+[float(vector[i])] for i,r in enumerate(matrix)]
 for i in range(n):
  j=max(range(i,n),key=lambda k:abs(a[k][i]))
  a[i],a[j]=a[j],a[i]
  t=a[i][i]
  if abs(t)<1e-12:raise ValueError("Nonsingular covariate rank matrix required")
  for k in range(i,n+1):a[i][k]/=t
  for j in range(n):
   if j==i:continue
   t=a[j][i]
   for k in range(i,n+1):a[j][k]-=t*a[i][k]
 return [r[n] for r in a]
def partial_corr(x,y,covs):
 n=len(x);k=len(covs[0])
 gram=[[sum(covs[i][a]*covs[i][b] for i in range(n)) for b in range(k)] for a in range(k)]
 def residual(v):
  rhs=[sum(covs[i][j]*v[i] for i in range(n)) for j in range(k)]
  betas=solve(gram,rhs)
  return [v[i]-sum(covs[i][j]*betas[j] for j in range(k)) for i in range(n)]
 return corr(residual(x),residual(y))
def build(root=ROOT):
 d=load(root)
 geb=d["data/frozen/figure_sources/host_contribution_metrics.csv"]
 ports=d["provenance/reviewer_defenses/results/port_connectivity_hosts_20261030/port_connectivity_host_all225.csv"]
 events=d["provenance/reviewer_defenses/results/island_host_continuity_20261030/island_host_continuity_events.csv"]
 islands=d["provenance/reviewer_defenses/results/historical_botanical_17island_20261031/near_pair_pre2011_plant_composition.json"]
 previous=d["provenance/reviewer_defenses/results/island_host_network_20261030/island_host_network_all225.csv"]
 contributors=d["provenance/reviewer_defenses/results/butterfly_host_contribution_concentration_v0.1.json"]
 g={r["species"]:r for r in geb}
 i={r["species"]:r for r in previous}
 if len(g)!=239 or len(ports)!=len(i)!=225 or len(set(r["species"] for r in ports))!=225:raise ValueError("GEB S1 and island source universes changed")
 if any(r["species"] not in g or r["species"] not in i for r in ports):raise ValueError("Species-key cohort crosswalk incomplete")
 areas={str(q["island_id"]) for q in islands["islands"]}
 assert len(areas)==17
 e=[r for r in events if r["landmass_id"] in areas and 2011<=int(r["first_record_year"])<=2017]
 if len(e)!=121 or len({r["species"] for r in e})!=64:raise ValueError("Original species×island first-index cohort changed")
 if any(r["species"] not in g for r in e):raise ValueError("GEB S1 differs from butterfly original")
 added=[int(float(r["introduced_added_units"])) for r in geb]
 native=[int(float(r["native_resource_units"])) for r in geb]
 rows=[]
 for p in ports:
  m=g[p["species"]];old=i[p["species"]]
  rows.append({"gain":float(m["introduced_added_units"]),
     "n_hosts":float(m["resolved_host_species"]),
     "native":float(m["native_resource_units"]),
     "previous":float(old["earlier_records_on_eligible_islands"]),
     "island_ports":float(p["island_ports_with_matching_alien_host"])})
 x=rank([r["gain"] for r in rows])
 y=rank([r["island_ports"] for r in rows])
 raw=corr(x,y)
 fit_host=partial_corr(x,y,[[1,math.log1p(r["n_hosts"])] for r in rows])
 fit_all=partial_corr(x,y,[[1,math.log1p(r["n_hosts"]),
   math.log1p(r["native"]),math.log1p(r["previous"])] for r in rows])
 base={
 "schema":"chocho_GEB239_to_port225_to_island64_observation_bridge_v1",
 "status":"VALIDATED_PUBLIC_SOURCE_COUNTS_RETROSPECTIVE_ASSOCIATIONS_NOT_GATEWAY_CAUSALITY",
 "source_git_blobs":PATHS,
 "GEB_resource_opportunity":{
   "butterfly_species":len(geb),"expanded_species":sum(v>0 for v in added),
   "native_butterfly_x_WGSRPD3_units":sum(native),
   "introduced_added_units":sum(added),
   "expanded_percent_relative_to_native":100*sum(added)/sum(native),
   "wcvp_areas_are_not_historical_port_sources":True},
 "aligned_port_candidate_comparison":{
   "species":len(ports),"with_global_geb_added_resource_units":sum(r["gain"]>0 for r in rows),
   "species_with_any_2001_2024_island_port_alien_candidate_host_record":sum(r["island_ports"]>0 for r in rows),
   "spearman_global_added_units_vs_count_of_island_ports_with_candidate_hosts":raw,
   "partial_rank_assoc_after_adjusting_host_count_only":fit_host,
   "partial_rank_assoc_after_adjusting_log1p_host_count_native_units_and_prior_recorded_islands":fit_all,
   "temporal_leak":"2001–2024 port flora and modern regional WCVP used; not prior-year forecast",
   "not_independent_from_host_catalogue":True},
 "indexed_2011_2017_island_first_records":{
   "island_count":len(areas),"events":len(e),"species":len(set(r["species"] for r in e)),
   "events_with_positive_global_resource_gain":sum(float(g[r["species"]]["introduced_added_units"])>0 for r in e),
   "species_with_positive_global_resource_gain":sum(float(g[s]["introduced_added_units"])>0 for s in {r["species"] for r in e}),
   "any_reindexed_during_next_five_years":sum(r["any_followup_record"]=="True" for r in e),
   "cases_selected_by_observed_indexing_not_a_risk_set":True},
 "top10_plant_contributors_from_geb":{
   "contributors":len(contributors["plant_species"]["top10"]),
   "fraction_of_all_added_resource_units":sum(float(z["share"]) for z in contributors["plant_species"]["top10"]),
   "shares_derived_from_source_github":True,
   "first_year_plant_port_coverage_is_separate_audit":True},
 "primary_limits":{
   "gep_gain_is_across_butterfly_x_region_units_not_physical_range_expansion":True,
   "port_index_does_not_establish_transport_or_larval_use":True,
   "port_plant_catalogue_and_global_geb_geography_share_plant_host_network_and_native_status":True,
   "correlations_descriptive_not_causal":True,
   "species_subset_and_island_subset_posthoc":True,
   "unrecorded_hosts_not_absent":True,
   "not_a_estimate_of_fraction_of_55_percent_realized":True}
 }
 return base
if __name__=="__main__":
 import argparse
 p=argparse.ArgumentParser()
 p.add_argument("--root",type=Path,default=ROOT)
 p.add_argument("--out",type=Path,default=ROOT/OUT)
 a=p.parse_args()
 data=build(a.root)
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(data,sort_keys=True,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
 print("GEB_PORT_BRIDGE:",json.dumps({"S1":data["GEB_resource_opportunity"],"ports":data["aligned_port_candidate_comparison"],"island":data["indexed_2011_2017_island_first_records"]}))
