"""Recalculate Spearman/tie average ranks, island graph component drop, quality gates.

Uses only public locked-source 17-node summaries and original 1275 pair
Jaccards. Botanical observation not larval-host identity or port introduction.
"""
import csv,json,math,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[4]
BASE=Path("provenance/reviewer_defenses/results")
H=BASE/"historical_botanical_17island_20261031/near_pair_pre2011_plant_composition.json"
P=BASE/"island_host_network_20261030/island_host_network_pairs.csv"
GH="54c4515947d2a3ee7545e4007985ef82a2af6868"
GP="5663d2ffe42419f5c4496b92c9f3bbf76175907b"
def sha(p):
 d=p.read_bytes()
 return hashlib.sha1(b"blob "+str(len(d)).encode()+b"\0"+d).hexdigest()
def ranks(a):
 idx=sorted(range(len(a)),key=lambda i:a[i]);out=[0.0]*len(a);k=0
 while k<len(a):
  j=k+1
  while j<len(a) and a[idx[k]]==a[idx[j]]:j+=1
  for n in idx[k:j]:out[n]=(k+j+1)/2
  k=j
 return out
def pearson(a,b):
 if len(a)<3:return None
 m=sum(a)/len(a);n=sum(b)/len(b)
 v=sum((x-m)**2 for x in a);w=sum((x-n)**2 for x in b)
 if not v or not w:return None
 return sum((x-m)*(y-n) for x,y in zip(a,b))/math.sqrt(v*w)
def rho(a,b):
 return pearson(ranks(a),ranks(b))
def summarize(root=ROOT):
 if sha(root/H)!=GH or sha(root/P)!=GP:raise ValueError("Source blob changed")
 src=json.loads((root/H).read_text(encoding="utf8"))
 pairs=list(csv.DictReader((root/P).open(encoding="utf8",newline="")))
 lookup={tuple(sorted([r["island_a"],r["island_b"]],key=int)):r for r in pairs}
 rows=[]
 for r in src["source_pairs"]:
  key=tuple(sorted([r["island_a"],r["island_b"]],key=int))
  r=dict(r)
  r["modern_host_jaccard"]=float(lookup[key]["host_jaccard"])
  rows.append(r)
 if len(rows)!=14 or sum(x["shared_new_butterfly_species_pairs"] for x in rows)!=41:
  raise ValueError("Original near butterfly population changed")
 def corr(v,a,b):return rho([float(z[a]) for z in v],[float(z[b]) for z in v])
 parent={}
 def find(x):
  if x not in parent:parent[x]=x
  if parent[x]!=x:parent[x]=find(parent[x])
  return parent[x]
 for r in rows:
  x,y=find(r["island_a"]),find(r["island_b"])
  if x!=y:parent[y]=x
 groups={}
 for r in rows:groups.setdefault(find(r["island_a"]),[]).append(r)
 def quality(r,t,y):
  return (r["plant_taxon_union_a"]>=t and r["plant_taxon_union_b"]>=t
      and r["plant_recorded_years_a"]>=y and r["plant_recorded_years_b"]>=y)
 def filtered(value):
  return {"n":len(value),"rho":corr(value,"historical_plant_jaccard","new_butterfly_jaccard")}
 return {
  "all":{"n":14,
    "rho_historical_broad_plants_vs_new_butterflies":corr(rows,"historical_plant_jaccard","new_butterfly_jaccard"),
    "rho_historical_broad_plants_vs_early_butterflies":corr(rows,"historical_plant_jaccard","prior_butterfly_jaccard"),
    "rho_modern_catalogue_host_jaccard_vs_new_butterflies":corr(rows,"modern_host_jaccard","new_butterfly_jaccard")},
  "nonzero_new_pairs":filtered([r for r in rows if r["shared_new_butterfly_species_pairs"]>0]),
  "acquisition_primary_quality":filtered([r for r in rows if quality(r,20,2)]),
  "exploratory_stricter_quality":filtered([r for r in rows if quality(r,50,3)]),
  "regional_components":[
    {"islands":sorted({i for r in v for i in [r["island_a"],r["island_b"]]},key=int),"n_pairs":len(v)}
    for v in groups.values()],
  "leave_one_geography_group_out":[
    {"omitted_islands":sorted({i for r in v for i in [r["island_a"],r["island_b"]]},key=int),
    "remaining_pairs":len(rows)-len(v),
    "rho":corr([r for r in rows if r not in v],"historical_plant_jaccard","new_butterfly_jaccard")}
    for v in groups.values()],
 }
def main():
 saved=json.loads((ROOT/"provenance/reviewer_defenses/results/historical_botanical_17island_20261031/historical_butterfly_similarity_sensitivity.json").read_text())
 actual=summarize()
 for k in actual:
  def equiv(a,b):
   if isinstance(a,float):return abs(a-b)<1e-10
   if isinstance(a,dict):return set(a)==set(b) and all(equiv(v,b[t]) for t,v in a.items())
   if isinstance(a,list):return len(a)==len(b) and all(equiv(v,w) for v,w in zip(a,b))
   return a==b
  if not equiv(actual[k],saved[k]):raise ValueError("Source sensitivity snapshot mismatched "+k)
 print("SOURCE17_HISTORICAL_BOTANICAL_SENSITIVITY_PASS")
if __name__=="__main__":main()
