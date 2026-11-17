#!/usr/bin/env python3
"""Does independent global pre-2011 butterfly COUNTRY breadth add geographically
transferable prediction of subsequent island GBIF FIRST INDEXING?

The dated 2001–2010 global country facets are OBSERVATION-index breadth, not
true native range or capacity to disperse. First-indexed 2011–2022 butterflies
may have been present but undetected; this is NOT a colonisation risk model.

Core out-of-fold B vs B+country comparisons use only source-period GBIF counts,
2001–2010 island source record geography and 2001–2010 early plant indexed
record effort, plus static WPI port and native HOSTS known-host catalogue.

Ancillary +H adds contemporary plant-defined number of island resource gains
*only as a retrospective explanation*, NOT a prior-in-time predictor. Keep
different GBIF country record supports and WCVP island regions separate.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np,pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score
from scipy.stats import spearmanr

ROOT=Path(__file__).resolve().parents[1]
PINS={
"species_index":("data/exploratory/public51_butterfly_index_20261109/full51_225butterfly_period_index.csv","b6a3c98eb76864052d45d520a55b7d830816126c"),
"global_breadth":("results/butterfly_global_GBIF_country_breadth_20261117/225_butterflies_observed_global_GBIF_country_breadth_and_HOSTS_benefit.csv","f263d30655a698a361afe0f1c7fd4c6bc6c3ed49"),
"historical_effort":("results/island51_pre2011_effort_control_20261107/historical_botanical_coverage_by_island.csv","84f4818bd83eb8f7c3ff745ef7aa4234310397aa"),
"island_geography":("results/island51_host_circuit_atlas_20261106/island_geography.csv","f0c5b1139e54b4872d28682f0991f89f2368ab43"),
"marine_pairs":("results/island51_observed_ocean_circuit_20261106/island51_observed_circuit_pairs.csv","7832b1d9d8b1047027fd5197d7bc0d31c22cad3f"),
"longitude_groups":("results/island51_double_unseen_blocks_20261107/longitude_island_blocks.csv","ef2cfcd0b3133d308c18b517001b33a9b46e9630")
}

def gitblob(path):
 b=path.read_bytes()
 return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()

def assemble():
 for name,(rel,blob) in PINS.items():
  got=gitblob(ROOT/rel)
  if got!=blob:raise ValueError("Frozen source changed "+name+" "+got)
 index=pd.read_csv(ROOT/PINS["species_index"][0],dtype={"species":str,"speciesKey":str,"landmass_id":str})
 breadth=pd.read_csv(ROOT/PINS["global_breadth"][0],dtype={"species":str})
 early=pd.read_csv(ROOT/PINS["historical_effort"][0],dtype={"island_id":str})
 island=pd.read_csv(ROOT/PINS["island_geography"][0],dtype={"island_id":str})
 marine=pd.read_csv(ROOT/PINS["marine_pairs"][0],dtype={"island_a":str,"island_b":str})
 groups=pd.read_csv(ROOT/PINS["longitude_groups"][0],dtype={"island_id":str})
 if (len(index),len(breadth),len(early),len(island),len(marine),len(groups))!=(11475,225,51,51,1275,51):
  raise ValueError("Source panel shape not original")
 if index.set_index(["species","landmass_id"]).index.has_duplicates:raise ValueError("Original species-island cells duplicated")
 if not set(island.island_id)==set(index.landmass_id)==set(early.island_id)==set(groups.island_id):
  raise ValueError("Original 51 landmasses no longer align")
 if not set(breadth.species)==set(index.species):
  raise ValueError("Global country facets and source butterfly index do not match exact 225 species")
 if early.original_port_count.sum()!=194 or island.original_ports.sum()!=194:
  raise ValueError("Original 194 WPI port counts drifted")
 all_butterfly_early=(index.assign(early=index.pre2001_2010_indexed_records.gt(0))
                      .groupby("landmass_id").early.sum().to_dict())
 earliest=(index.assign(early=index.pre2001_2010_indexed_records.gt(0))
           .pivot(index="species",columns="landmass_id",values="early"))
 dist={}
 for r in marine.itertuples():
  a,b=r.island_a,r.island_b
  dist[(a,b)]=dist[(b,a)]=(float(r.port_distance_km),float(r.log1p_normalized_uniform_circuit))
 if len(dist)!=2550:raise ValueError("Original complete pairwise marine source not recovered")
 group=groups.set_index("island_id").longitude_block.to_dict()
 geo=island.set_index("island_id")
 sampling=early.set_index("island_id")
 species= breadth.set_index("species")
 arr=[]
 for z in index.itertuples():
  if int(z.pre2001_2010_indexed_records)>0:continue
  sp=z.species;loc=z.landmass_id
  met=species.loc[sp]
  # Exclusion requires a known sufficient *independent* country count for
  # fixed-30 rarefaction; never assign unrecorded as local endemic.
  if not np.isfinite(float(met.rarefied_30_countries)):continue
  prior=[k for k,yes in earliest.loc[sp].items() if yes and k!=loc]
  nearestGeo=min((dist[(loc,k)][0] for k in prior),default=0.)
  nearestSea=min((dist[(loc,k)][1] for k in prior),default=0.)
  station=geo.loc[loc]
  effort=sampling.loc[loc]
  y=int(int(z.first_indexed_in_2011_2017)==1 or int(z.first_indexed_in_2018_2022)==1)
  row={
   "species":sp,"landmass_id":loc,"longitude_group":int(group[loc]),
   "genus":sp.split(" ",1)[0],"future_first_GBIF_index":y,
   "rarefied_30_global_country_breadth":float(met.rarefied_30_countries),
   "log_global_country_record_count":np.log1p(float(met.n_country_facet_records)),
   "log_source_species_early_51islands":np.log1p(float(met.pre2011_51_island_indexed_count)),
   "known_native_HOSTS_geography_log":np.log1p(float(met.native_resource_units)),
   "known_NATIVE_island_HOSTS_geography_log":np.log1p(float(met.native_island_regions)),
   "known_HOSTS_species_count_log":np.log1p(float(met.resolved_host_species)),
   "known_HOSTS_families_log":np.log1p(float(met.host_family_count)),
   "source_lacking_early_conspecific_index":int(len(prior)==0),
   "log_n_prior_conspecific_islands":np.log1p(float(len(prior))),
   "log_geodesic_km_nearest_early_conspecific_index":np.log1p(nearestGeo),
   "log_undirected_marine_resistance_nearest_early_conspecific_index":np.log1p(nearestSea),
   "log_pre2011_botanical_record_effort_194_original_ports":np.log1p(float(effort.early_plant_buffer_record_count_sum)),
   "log_original_WPI_port_count":np.log1p(float(station.original_ports)),
   "log_early_butterflies_recorded_on_island":np.log1p(float(all_butterfly_early[loc])),
   "latitude":float(station.representative_latitude),
   "sin_longitude":float(np.sin(np.radians(station.representative_longitude))),
   "cos_longitude":float(np.cos(np.radians(station.representative_longitude))),
   "EXPLORATORY_CURRENT_HOST_island_resource_gain_log":np.log1p(float(met.added_island_regions))}
  arr.append(row)
 data=pd.DataFrame(arr)
 if data.future_first_GBIF_index.sum()<290 or len(data)<10500:
  raise ValueError("Wrong original 225x51 2001–2022 risk cohort on 30-record country source")
 return data

def fitpredict(X,y,train,test):
 scale=StandardScaler().fit(X[train])
 x=scale.transform(X[train]);v=scale.transform(X[test])
 m=LogisticRegression(solver="liblinear",penalty="l2",C=1.,max_iter=400,random_state=20261117)
 m.fit(x,y[train])
 p=m.predict_proba(v)[:,1]
 if not np.isfinite(p).all() or np.any(p<=0) or np.any(p>=1):raise ValueError("Invalid first-index probability")
 return p

def evalsource(data):
 B=["log_global_country_record_count",
   "log_source_species_early_51islands","known_native_HOSTS_geography_log",
   "known_NATIVE_island_HOSTS_geography_log","known_HOSTS_species_count_log",
   "known_HOSTS_families_log","source_lacking_early_conspecific_index",
   "log_n_prior_conspecific_islands","log_geodesic_km_nearest_early_conspecific_index",
   "log_undirected_marine_resistance_nearest_early_conspecific_index",
   "log_pre2011_botanical_record_effort_194_original_ports",
   "log_original_WPI_port_count","log_early_butterflies_recorded_on_island",
   "latitude","sin_longitude","cos_longitude"]
 C="rarefied_30_global_country_breadth"
 H="EXPLORATORY_CURRENT_HOST_island_resource_gain_log"
 models={"B":B,"B_plus_global_pre2011_country_breadth":B+[C],
         "B_plus_CONCURRENT_host_resource_gain_retrospective":B+[H],
         "B_plus_both_geographies_retrospective":B+[C,H]}
 y=data.future_first_GBIF_index.to_numpy(int)
 x={k:data[v].to_numpy(float) for k,v in models.items()}
 labels={
  "LOIO_51_islands":data.landmass_id.astype(str).to_numpy(),
  "source_frozen_5_geographic_groups":data.longitude_group.to_numpy(),
  "heldout_butterfly_genus_5_hashes":np.array([
      int(hashlib.sha256(("genus20261117_"+a).encode()).hexdigest()[:8],16)%5
      for a in data.genus.astype(str)])}
 result={};folds=[]
 for label,cluster in labels.items():
  klist=sorted(np.unique(cluster))
  pred={m:np.full(len(data),np.nan) for m in models}
  for k in klist:
   test=cluster==k;train=~test
   if test.sum()<5 or train.sum()<7000 or y[train].sum()<15:
    raise ValueError("Invalid heldout butterfly source fold")
   row={"design":label,"fold":str(k),"train_cells":int(train.sum()),"test_cells":int(test.sum()),
        "events":int(y[test].sum())}
   for name in models:
    p=fitpredict(x[name],y,train,test)
    pred[name][test]=p
    row["Brier_"+name]=float(np.mean((y[test]-p)**2))
    row["logloss_"+name]=float(np.mean(-(y[test]*np.log(p)+(1-y[test])*np.log1p(-p))))
   folds.append(row)
  scores={}
  for name in models:
   p=pred[name]
   if not np.isfinite(p).all():raise ValueError("No OOF prediction for "+label+"/"+name)
   scores[name]={"OOF_Brier":float(np.mean((y-p)**2)),
                 "OOF_logloss":float(np.mean(-(y*np.log(p)+(1-y)*np.log1p(-p)))),
                 "OOF_AUC":float(roc_auc_score(y,p))}
  base=scores["B"]
  inc={}
  for name in models:
   if name=="B":continue
   inc[name]={
    "Brier_error_reduction_pct":100*(base["OOF_Brier"]-scores[name]["OOF_Brier"])/base["OOF_Brier"],
    "logloss_error_reduction_pct":100*(base["OOF_logloss"]-scores[name]["OOF_logloss"])/base["OOF_logloss"],
    "AUC_delta":scores[name]["OOF_AUC"]-base["OOF_AUC"]}
  result[label]={"n_folds":len(klist),"test_total":len(data),"events":int(y.sum()),
                 "source_time":"GBIF global country and island early source 2001–2010; fixed present day HOSTS catalogue and WPI geography",
                 "scores":scores,"incremental_over_base":inc}
 return result,folds,models

def main():
 p=argparse.ArgumentParser();p.add_argument("--out-dir",type=Path,required=True);a=p.parse_args()
 data=assemble()
 scores,folds,models=evalsource(data)
 a.out_dir.mkdir(parents=True,exist_ok=True)
 pd.DataFrame(folds).to_csv(a.out_dir/"GLOBAL_PRE2011_BUTTERFLY_BREADTH_OOF_LOSSES.csv",index=False)
 score={
  "schema":"chocho_butterfly_global_GBIF_pre2011_country_breadth_vs_late51_island_first_index_record_OOF_v1",
  "status":"SOURCE_VERIFIED_RETROSPECTIVE_FIRST_GBIF_INDEX_PREDICTION_NOT_ESTABLISHMENT",
  "source_git_blobs":{name:{"path":p,"sha":sha} for name,(p,sha) in PINS.items()},
  "n_risk_cells":len(data),"n_islands":51,"n_butterflies":data.species.nunique(),
  "n_future_first_index_events":int(data.future_first_GBIF_index.sum()),
  "n_excluded_butterflies_fewer_than30_country_records":225-data.species.nunique(),
  "heldout_scores":scores,
  "models":models,
  "limits":[
   "Global 2001–2010 country facets estimate observed breadth, not native butterfly range/endemicity/dispersal ability",
   "Current HOSTS native plant geography and catalogued host counts are not source dated before the first butterfly record; models are retrospective even though independent country-index covariates are restricted to early dates",
   "The HOST gain predictor is CONTEMPORARY WCVP 2024 and cannot be interpreted as historical pre-2011 local food; only auxiliary retrospective model",
   "Independent global GBIF country facet and island first-index response share GBIF record-provider and effort biases",
   "The original 51 selected landmasses are not the 126 botanical WGSRPD3 islands",
   "First GBIF indexing is not biological arrival, successful butterfly population growth or true community homogenization",
   "Historic conspecific source distances use only the original 51 sampled landmasses, excluding external potential sources",
   "Geographic/group/generalization folds share species and islands; these are out-of-fold predictive comparisons not causal exposure experiments",
   "No tested term establishes shipborne carriage, ocean-current movement, or plant-introduction-mediated colonization"]}
 (a.out_dir/"GLOBAL_COUNTRY_BREADTH_AND_HOST_RESOURCE_FIRSTINDEX_OOF_20261117.json").write_text(json.dumps(score,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
 print("INDEPENDENT_BUTTERFLY_GBIF_BREADTH_TRANSFER="+json.dumps(
     {"cells":score["n_risk_cells"],"events":score["n_future_first_index_events"],"species":score["n_butterflies"],
       "gains":{design:z["incremental_over_base"] for design,z in scores.items()}},ensure_ascii=False,sort_keys=True),flush=True)

if __name__=="__main__":
 main()
