#!/usr/bin/env python3
"""Source-locked historical host-plant records × butterfly species-island indexing.

Original author-derived butterfly–HOSTS plant edges are read transiently from
a frozen Actions artifact; only aggregate model scores are saved. The 2001–10
plant records precede 2011–22 butterfly indexing but are NOT proven larval
feeding, biological butterfly arrivals, introduced plants or real sea currents.
"""
from __future__ import annotations
import argparse,csv,gzip,hashlib,io,json,math
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.linalg import lsqr
from scipy.stats import t as t_distribution
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

PINS={
 "hosts":"e55012d6b0eefd2b3e9e4b8f13fc4786ddc1e66c21f73e705bb20df92d0b5e22",
 "taxonomy":"d8d21e6686e3f069d4374d4be1c4e47744788f9244c46a1bbbbb9015e5cd7fb9",
 "plantkeys":"42bca05d45e97dd38f1b879fc0934828fc9b7451a2b7848a29eac11c12e23061",
 "port_coverage":"ec172084494d73284fdbaeaacd59d23c030d390884ca8d409f2f23e96d8559fb",
 "butterfly_index":"d8fbb3fbf936bfdbedc6c91f1e19f11deda882d3125ccb6600a9c8d800a57887",
 "marine_pairs":"133dc6d6a91203d09385a4ffd4db0b712dc3309130d1bf2a09bcee6bafa36bf0",
 "marine_prep":"5641711fc7ab2bdcd4626e1f483faf46c945a5d11466f0deada092db456a0363",
 "longitude":"2e4b8b1d4bda18f7193f843378da124b597ead0f6c49751a547cab10f65e2ef3"}
def get(path,name):
 b=Path(path).read_bytes()
 observed=hashlib.sha256(b).hexdigest()
 if observed!=PINS[name]:raise ValueError("Original input SHA256 changed: "+name+" "+observed)
 return b
def csvrows(blob):return list(csv.DictReader(io.StringIO(blob.decode("utf-8"))))

def construct(a):
 edges=csvrows(get(a.hosts,"hosts"))
 taxonomy=json.loads(gzip.decompress(get(a.taxonomy,"taxonomy")))["name_resolution"]
 plants=csvrows(gzip.decompress(get(a.plants,"plantkeys")))
 ports=csvrows(get(a.port_coverage,"port_coverage"))
 butterfly=pd.read_csv(io.BytesIO(get(a.butterfly,"butterfly_index")),
                       dtype={"species":str,"speciesKey":str,"landmass_id":str})
 marine=pd.read_csv(io.BytesIO(get(a.marine,"marine_pairs")),dtype={"island_a":str,"island_b":str})
 prep=json.loads(get(a.prep,"marine_prep"))
 groups=pd.read_csv(io.BytesIO(get(a.longitude,"longitude")),dtype={"island_id":str})
 if len(edges)!=2524 or len(taxonomy)!=1642 or len(plants)!=27598 or len(ports)!=194 or len(butterfly)!=11475 or len(marine)!=1275 or len(groups)!=51:
  raise ValueError("Original source dimensions do not match archived records")
 valid_keys={str(v["gbif_speciesKey"]) for v in taxonomy.values()
             if v.get("gate")=="EXACT_ACCEPTED" and v.get("gbif_speciesKey")}
 if len(valid_keys)!=1581:raise ValueError("Accepted species mapping changed")
 seen={v for v in butterfly.species.unique()}
 host=defaultdict(set);direct=set();alias=set();alias_target=defaultdict(set)
 for e in edges:
  if e["query_eligible"]!="True":continue
  status=taxonomy.get(e["query_name"],{})
  if status.get("gate")!="EXACT_ACCEPTED" or not status.get("gbif_speciesKey"):continue
  if e["current_usage_status"]=="ACCEPTED":
   name=e["current_accepted_name"];kind="direct"
  elif e["current_usage_status"]=="SYNONYM":
   name=e["original_name"];kind="alias"
  else:continue
  if name not in seen:continue
  host[name].add(str(status["gbif_speciesKey"]))
  if kind=="direct":direct.add(name)
  else:alias.add(name);alias_target[name].add(e["current_accepted_name"])
 if any(len(v)!=1 for v in alias_target.values()) or direct&alias:
  raise ValueError("Ambiguous original butterfly species alias")
 matched={s for s in seen if host.get(s)}
 if (len(matched),len(direct),len(alias))!=(218,203,15):
  raise ValueError("Conservative source butterfly crosswalk changed")
 allplant=defaultdict(set);allrec=defaultdict(int)
 for p in plants:
  if int(p["indexed_records_2001_2010"])<=0:raise ValueError("Invalid dated plant facet count")
  allplant[str(p["island_id"])].add(str(p["GBIF_speciesKey"]))
 original_ports=defaultdict(int)
 for p in ports:
  original_ports[str(p["island_id"])]+=1
  allrec[str(p["island_id"])]+=int(p["live_GBIF_2001_2010_plant_records"])
 if sum(allrec.values())!=173113 or sum(original_ports.values())!=194:raise ValueError("Original 194-port GBIF counts changed")
 islands=sorted(prep["by_island"],key=int)
 if len(islands)!=51 or set(islands)!=set(groups.island_id) or set(islands)!=set(original_ports):raise ValueError("Original island IDs mismatch")
 if any(original_ports[k]!=int(prep["by_island"][k]["original_ports"]) for k in islands):raise ValueError("Original WPI port identities mismatch")
 S=butterfly[butterfly.species.isin(matched)].copy().reset_index(drop=True)
 sp=sorted(matched);ms=len(sp);n=len(islands)
 smap={v:i for i,v in enumerate(sp)}
 imap={v:i for i,v in enumerate(islands)}
 S["sp"]=S.species.map(smap).astype(int);S["island"]=S.landmass_id.map(imap).astype(int)
 S["pre"]=S.pre2001_2010_indexed_records.astype(int)
 S["y2011"]=S.first_indexed_in_2011_2017.astype(int)
 S["y2018"]=S.first_indexed_in_2018_2022.astype(int)
 S["y"]=((S.y2011+S.y2018)>0).astype(int)
 S=S[S.pre.eq(0)].reset_index(drop=True)
 S["has_prior_host"]=np.array([int(bool(host[b]&allplant[k])) for b,k in zip(S.species,S.landmass_id)],dtype=int)
 S["n_prior_hosts"]=np.array([len(host[b]&allplant[k]) for b,k in zip(S.species,S.landmass_id)],dtype=int)
 S["n_known_hosts"]=S.species.map(lambda b:len(host[b]))
 S["n_plant_species"]=S.landmass_id.map(lambda k:len(allplant[k]))
 S["n_plant_records"]=S.landmass_id.map(lambda k:allrec[k])
 S["n_ports"]=S.landmass_id.map(lambda k:original_ports[k])
 before=butterfly.assign(pre=butterfly.pre2001_2010_indexed_records.astype(int))
 prior_by_species=before[before.pre.gt(0)].groupby("species").landmass_id.agg(list).to_dict()
 geod=np.zeros((n,n));ocean=np.zeros((n,n))
 for p in marine.itertuples(index=False):
  i=imap[p.island_a];j=imap[p.island_b]
  geod[i,j]=geod[j,i]=float(p.port_distance_km)
  ocean[i,j]=ocean[j,i]=float(p.log1p_normalized_uniform_circuit)
 if not (np.all(geod[np.triu_indices(n,1)]>0) and np.all(ocean[np.triu_indices(n,1)]>0)):
  raise ValueError("Original marine circuit or WPI source distances changed")
 psets={b:[imap[k] for k in ks] for b,ks in prior_by_species.items() if b in matched}
 S["n_prior_sources"]=S.species.map(lambda k:len(psets.get(k,[])))
 S["nearest_prior_geo_km"]=np.array([min((geod[int(i),j] for j in psets.get(b,[]) if j!=int(i)),default=0) for b,i in zip(S.species,S.island)])
 S["nearest_prior_marine"]=np.array([min((ocean[int(i),j] for j in psets.get(b,[]) if j!=int(i)),default=0) for b,i in zip(S.species,S.island)])
 S["log_source_km"]=np.log1p(S.nearest_prior_geo_km)
 S["log_uniform"]=np.log1p(S.nearest_prior_marine)
 locations=groups.set_index("island_id")
 S["lat"]=S.landmass_id.map(locations.representative_latitude)
 S["lon"]=S.landmass_id.map(locations.representative_longitude)
 S["geo_group"]=S.landmass_id.map(locations.longitude_block).astype(int)
 early_by_island=before.groupby("landmass_id").pre.agg(lambda x:int(sum(x>0))).to_dict()
 S["n_prior_butterflies"]=S.landmass_id.map(early_by_island)
 summary={
  "original_HOSTS_butterfly_edges_read_transiently":2524,
  "accepted_HOSTS_plant_name_keys":len(valid_keys),
  "historically_indexed_accepted_HOSTS_plant_keys":len(valid_keys&{str(r["GBIF_speciesKey"]) for r in plants}),
  "matched_butterflies":len(matched),"strict_name_matched_species":len(direct),
  "verified_original_source_synonym_matches":len(alias),
  "original_ports":194,"original_port_speciesKey_rows":len(plants),
  "original_early_plant_observation_sum_overlapping_buffers":sum(allrec.values()),
  "islands_with_early_general_plant_keys":sum(bool(allplant[i]) for i in islands),
  "islands_with_early_known_HOSTS_plant_record":sum(bool(valid_keys&allplant[i]) for i in islands),
  "at_risk_species_island_cells":len(S),"future_first_record_cells":int(S.y.sum()),
  "has_pre2011_compatible_HOSTS_plant":int(S.has_prior_host.sum()),
  "new_first_indices_with_prior_host":int((S.y*S.has_prior_host).sum()),
  "new_first_indices_without_prior_host":int((S.y*(1-S.has_prior_host)).sum())}
 if (summary["at_risk_species_island_cells"],summary["future_first_record_cells"],
     summary["has_pre2011_compatible_HOSTS_plant"],summary["new_first_indices_with_prior_host"])!=(10431,306,1075,53):
  raise ValueError("Original temporally prior exact species×island join changed")
 return S,ms,n,summary,direct

def two_way_lpm(S,ms,ni,controls,output=None):
 S=S.reset_index(drop=True);N=len(S)
 y=np.asarray(S.y if output is None else S[output],dtype=float)
 x=np.asarray(S.has_prior_host,dtype=float)
 a=S.sp.to_numpy(int);b=S.island.to_numpy(int)
 sa=sparse.csr_matrix((np.ones(N),(np.arange(N),a)),shape=(N,ms))
 ia=sparse.csr_matrix((np.ones(N),(np.arange(N),b)),shape=(N,ni))
 other=[sparse.csr_matrix(S[k].to_numpy(float)[:,None]) for k in controls]
 W=sparse.hstack([sa,ia]+other,format="csr")
 rx=x-W@lsqr(W,x,atol=1e-9,btol=1e-9,iter_lim=800)[0]
 ry=y-W@lsqr(W,y,atol=1e-9,btol=1e-9,iter_lim=800)[0]
 denom=float(rx@rx)
 beta=float(rx@ry/denom)
 err=ry-beta*rx;val=rx*err
 va=np.bincount(a,weights=val,minlength=ms)
 vb=np.bincount(b,weights=val,minlength=ni)
 var=(np.dot(va,va)+np.dot(vb,vb)-np.dot(val,val))/denom**2
 se=float(np.sqrt(max(0.,var)))
 dof=min(len(set(a)),len(set(b)))-1
 qt=float(t_distribution.ppf(.975,dof))
 return {"n":N,"events":int(y.sum()),"n_species":len(set(a)),"n_islands":len(set(b)),
         "estimated_risk_difference_pp":beta*100,"two_way_cluster_SE_pp":se*100,
         "descriptive_95pct_interval_pp":[(beta-qt*se)*100,(beta+qt*se)*100],
         "controls":["butterfly_species_FE","island_FE"]+controls,
         "interval_status":"descriptive two-way cluster robust LPM, not independent biological replicate"}

def heldout(S,ms):
 S=S.reset_index(drop=True);rids=S.sp.to_numpy(int);y=S.y.to_numpy(int)
 attrs=[
  np.log1p(S.nearest_prior_geo_km.to_numpy(float)),
  np.log1p(S.nearest_prior_marine.to_numpy(float)),
  S.n_prior_sources.eq(0).to_numpy(float),
  np.log1p(S.n_prior_sources.to_numpy(float)),
  np.log1p(S.n_plant_records.to_numpy(float)),
  np.log1p(S.n_plant_species.to_numpy(float)),
  S.n_plant_species.gt(0).to_numpy(float),
  np.log1p(S.n_prior_butterflies.to_numpy(float)),
  np.log1p(S.n_ports.to_numpy(float)),
  S.lat.to_numpy(float),np.sin(np.radians(S.lon.to_numpy(float))),
  np.cos(np.radians(S.lon.to_numpy(float))),
  np.log1p(S.n_known_hosts.to_numpy(float))]
 X=np.column_stack(attrs);h=S.has_prior_host.to_numpy(float)[:,None]
 onehot=sparse.csr_matrix((np.ones(len(S)),(np.arange(len(S)),rids)),shape=(len(S),ms))
 allfold={}
 for label,keys in [("LOIO",sorted(S.landmass_id.unique(),key=int)),
                    ("geography_5_groups",[0,1,2,3,4])]:
  rows=[]
  for island in keys:
   test=(S.landmass_id==island).to_numpy() if label=="LOIO" else (S.geo_group==island).to_numpy()
   train=~test
   if not test.any() or len(np.unique(y[train]))!=2:raise ValueError("Invalid held island split")
   scale=StandardScaler().fit(X[train])
   xt=scale.transform(X[train]);xv=scale.transform(X[test])
   fr={"fold":str(island),"n_train":int(train.sum()),"n_test":int(test.sum()),"n_first_indices":int(y[test].sum())}
   for model in ("B","B_plus_prior_host"):
    if model=="B":
     tr=sparse.hstack([sparse.csr_matrix(xt),onehot[train]],format="csr")
     te=sparse.hstack([sparse.csr_matrix(xv),onehot[test]],format="csr")
    else:
     tr=sparse.hstack([sparse.csr_matrix(xt),sparse.csr_matrix(h[train]),onehot[train]],format="csr")
     te=sparse.hstack([sparse.csr_matrix(xv),sparse.csr_matrix(h[test]),onehot[test]],format="csr")
    model_fit=LogisticRegression(solver="liblinear",C=1.0,random_state=20261117,max_iter=500)
    model_fit.fit(tr,y[train])
    p=model_fit.predict_proba(te)[:,1]
    fr["Brier_"+model]=float(np.mean((y[test]-p)**2))
    fr["Logloss_"+model]=float(np.mean(-(y[test]*np.log(p)+(1-y[test])*np.log1p(-p))))
   rows.append(fr)
  frame=pd.DataFrame(rows)
  weighted=label!="LOIO"
  summary={}
  for metric in ("Brier","Logloss"):
   w=frame.n_test.to_numpy(float) if weighted else np.ones(len(frame))
   a=np.average(frame[metric+"_B"],weights=w)
   b=np.average(frame[metric+"_B_plus_prior_host"],weights=w)
   summary[metric]={"base":float(a),"with_prior_known_host":float(b),
                    "relative_error_reduction_pct":float(100*(a-b)/a)}
  allfold[label]=(summary,frame)
 return allfold

def main():
 ap=argparse.ArgumentParser()
 for arg in ("hosts","taxonomy","plants","port-coverage","butterfly","marine","prep","longitude","output"):
  ap.add_argument("--"+arg,required=True,type=Path)
 a=ap.parse_args();out=a.output;out.mkdir(parents=True,exist_ok=True)
 S,ms,ni,src,direct=construct(a)
 field=lambda t,ss:S[ss].assign(y=t)[["y","sp","island","has_prior_host","log_source_km","log_uniform"]]
 primary={}
 for name,population in [
   ("all_2011_2022",S),
   ("general_plant_indexed",S[S.n_plant_species.gt(0)]),
   ("original_current_accepted_name_only",S[S.species.isin(direct)])]:
  primary[name]={
   "species_and_island_FE":two_way_lpm(population,ms,ni,[]),
   "species_island_FE_plus_geographic_and_marine_prior_conspecific_access":
     two_way_lpm(population,ms,ni,["log_source_km","log_uniform"])}
 later=S[S.y2011.eq(0)].copy()
 primary["2011_2017_only"]={"geographic_controls":two_way_lpm(S.assign(y=S.y2011),ms,ni,["log_source_km","log_uniform"])}
 primary["2018_2022_only"]={"geographic_controls":two_way_lpm(later.assign(y=later.y2018),ms,ni,["log_source_km","log_uniform"])}
 folds=heldout(S,ms)
 summary={"schema":"chocho_original230_HOSTS_194port_2001_2010_exact_known_host_GBIF_225by51_2011_2022_v1",
  "status":"SOURCE_PINNED_PRIOR_HOST_SPECIES_ISLAND_DIAGNOSTIC_NOT_COLONISATION",
  "provenance_actions":{"host_edges":40979986160,"host_taxonomy":41177682289,
    "dated_194port_plants":42361257543,"public_butterflies_225x51":42556346640,
    "uniform_marine_circuit":42193019385,"source_only_geography":42325176725},
  "input_sha256":PINS,"sample":src,"associations":primary,
  "out_of_island_prediction":{k:{"n_folds":len(v[1]),"n_test_cells":int(v[1].n_test.sum()),
     "n_indexed_events":int(v[1].n_first_indices.sum()),"scores":v[0],
     "folds_improved_Brier":int((v[1]["Brier_B_plus_prior_host"]<v[1]["Brier_B"]).sum())}
    for k,v in folds.items()},
  "method":"Historical 2001–2010 speciesKey observations at 194 original port 5km buffers intersect an exact accepted GBIF taxon-key crosswalk derived from butterfly–HOSTS, evaluated against independent later 2011–2022 GBIF first index, at-risk only where focal butterfly unindexed before 2011.",
  "limits":["The author-derived individual HOSTS butterfly edges are not included in this output; do not publish raw confidential source or joined butterfly-host species×island table",
    "Pre-2011 plant GBIF indexing does not imply plant introduction, local feeding, larval success, biological arrival, or true absence when unindexed",
    "The 51 sampled landmasses and 225 butterflies differ from 104 WCVP regions and global 239-species resource opportunity",
    "Measured undirected GSHHG electrical resistance does not represent directed currents or vessel carriage",
    "Prior conspecific sources only among the sampled 51 islands, not true nearest biological populations",
    "Predictive folds share species and histories; two-way cluster intervals are descriptive and not confirmatory causal inference"]}
 for k,(r,fs) in folds.items():
  fs.to_csv(out/("heldout_"+k+"_fold_scores.csv"),index=False)
 (out/"SOURCE_VERIFIED_PRIOR_2010_HOSTS_FIRSTINDEX_20261117.json").write_text(json.dumps(summary,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
 print("PRIOR_HOST_SOURCE_VERIFIED",json.dumps({"sample":src,"FE_adjusted":primary["all_2011_2022"]["species_island_FE_plus_geographic_and_marine_prior_conspecific_access"],"pred":summary["out_of_island_prediction"]},ensure_ascii=False),flush=True)
if __name__=="__main__":main()
