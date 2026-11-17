#!/usr/bin/env python3
"""Does the INDEPENDENT 225-by-51 GBIF butterfly species index become more similar?

Compare cumulative species records indexed in 2001-2010, through 2017, through
2022 on the same 51 WPI-linked sampled landmasses. Newly indexed records are
not actual colonisations; later sources may reflect effort and historical
presence. Test 2010-to-2022 Jaccard gain against constrained binary source
record rearrangements keeping each butterfly's later new island count, each
island's number of newly indexed butterfly species and the structural absence
of duplicated first-index events. Also hold butterfly-by-frozen-longitude-group
new indices fixed in the more demanding geography-null.

Read only original source-pinned public 225*51 cells, not unavailable 83k
private records. Optional HOSTS resource-footprint quartiles are supplied by the
existing independent 239-butterfly resource-reassembly analysis. Never equate
host-plant geographic breadth to the butterfly's true native range.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr,rankdata

ROOT=Path(__file__).resolve().parents[1]
PINS={
 "butterflies":("data/exploratory/public51_butterfly_index_20261109/full51_225butterfly_period_index.csv","b6a3c98eb76864052d45d520a55b7d830816126c"),
 "longitude":("results/island51_double_unseen_blocks_20261107/longitude_island_blocks.csv","ef2cfcd0b3133d308c18b517001b33a9b46e9630"),
 "beneficiary":("results/butterfly_resource_beneficiaries_20261117/beneficiary_species_native_resource_breadth.csv","c9e409441dcd2bb7f9edc5a813fc13490bf60d8f"),
}
N_DRAWS=399
SEED=20261117

def git_blob(path):
 b=path.read_bytes()
 return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\x00"+b).hexdigest()

def source():
 digest={}
 for key,(path,sha) in PINS.items():
  p=ROOT/path
  observed=git_blob(p)
  if observed!=sha:raise ValueError("Original source Git blob mismatch "+path+": "+observed)
  digest[key]={"path":path,"git_blob_sha":sha}
 return digest

def jaccard(A):
 X=A.astype(np.int32)
 sizes=np.sum(X,axis=0)
 overlap=X.T@X
 union=sizes[:,None]+sizes[None,:]-overlap
 idx=np.triu_indices(A.shape[1],1)
 valid=union[idx]>0
 out=np.full(len(valid),np.nan)
 out[valid]=overlap[idx][valid]/union[idx][valid]
 return out,valid

def first_index_null(E0,B,group,rng,regional,draws=N_DRAWS):
 E=E0.copy()
 positions=np.argwhere(E).astype(np.int32)
 n_edges=len(positions)
 n_species,n_islands=E.shape
 if not n_edges:raise ValueError("No original indexed future records")
 before_species=E.sum(axis=1)
 before_island=E.sum(axis=0)
 before_by_region=np.stack([E[:,group==g].sum(axis=1) for g in range(5)],axis=1)
 prev_edges=B.sum()
 proposals=accepted=0
 def mix(propose):
  nonlocal proposals,accepted
  picks=rng.integers(0,n_edges,size=(propose,2))
  for x,y in picks:
   proposals+=1
   if x==y:continue
   s1,i1=positions[x]
   s2,i2=positions[y]
   if s1==s2 or i1==i2:continue
   if regional and group[i1]!=group[i2]:continue
   if B[s1,i2] or B[s2,i1] or E[s1,i2] or E[s2,i1]:continue
   E[s1,i1]=False;E[s2,i2]=False
   E[s1,i2]=True;E[s2,i1]=True
   positions[x,1]=i2;positions[y,1]=i1
   accepted+=1
 def validate():
  if np.any(E&B):raise AssertionError("Original butterflies reindexed in 2001-2010")
  if not np.array_equal(E.sum(axis=1),before_species) or not np.array_equal(E.sum(axis=0),before_island):
   raise AssertionError("Fixed butterfly or island new-record margin drift")
  if regional:
   for g in range(5):
    if not np.array_equal(E[:,group==g].sum(axis=1),before_by_region[:,g]):
     raise AssertionError("Butterfly by longitude-block new-index margin drift")
 mix(65*n_edges)
 validate()
 draws=[]
 for k in range(draws):
  mix(9*n_edges)
  if k%99==0 or k==draws-1:
   validate()
  z,_=jaccard(B|E)
  draws.append(z)
  if (k+1)%100==0:print("NULL_SOURCE_COMPOSITION_PROGRESS",("within_longitude" if regional else "global"),k+1,draws,flush=True)
 D=np.asarray(draws,dtype=np.float64)
 return D,{"draws":len(D),"proposals":proposals,"accepted":accepted,"swap_acceptance_fraction":accepted/proposals,
           "burn_proposals":65*n_edges,"per_draw_proposals":9*n_edges,
           "source_margin_checks":"PASS",
           "fixed_butterfly_new_island_counts":True,
           "fixed_island_new_butterfly_species":True,
           "butterfly_by_original_longitude_group_counts_fixed":bool(regional),
           "source_period":"2001-2010 pre-indexed; 2011-2022 new record-index events"}

def summary(obs_early,obs_late,valid,subset):
 a=obs_early[valid&subset];b=obs_late[valid&subset]
 d=b-a
 return {"pairs":len(a),"mean_Jaccard_2001_2010":float(a.mean()),
         "mean_Jaccard_by_2022":float(b.mean()),
         "mean_increment":float(d.mean()),
         "median_increment":float(np.median(d)),
         "improved_pairs":int((d>0).sum()),
         "worsened_pairs":int((d<0).sum()),
         "unchanged_pairs":int((d==0).sum()),
         "share_pairs_increasing":float((d>0).mean())}

def measure_null(obs,base,nulldraws,which):
 low=np.asarray(nulldraws,dtype=float)
 inc=np.nanmean(low[:,which]-base[which][None,:],axis=1)
 return {
  "observed_mean_Jaccard_increment":float(obs),
  "null_mean":float(inc.mean()),
  "null_descriptive_interval95":[float(x) for x in np.quantile(inc,[.025,.975])],
  "extra_beyond_null":float(obs-inc.mean()),
  "descriptive_upper_tail":float((1+(inc>=obs).sum())/(len(inc)+1)),
  "descriptive_lower_tail":float((1+(inc<=obs).sum())/(len(inc)+1)),
  "draws":int(len(inc)),
  "null_CV_warning":"MCMC swaps are dependent, observation process is uneven and sample is selected; not an inferential p-value for true island butterfly homogenization"}

def analyze(args):
 pins=source()
 df=pd.read_csv(ROOT/PINS["butterflies"][0],dtype={"species":str,"speciesKey":str,"landmass_id":str})
 geo=pd.read_csv(ROOT/PINS["longitude"][0],dtype={"island_id":str})
 benefits=pd.read_csv(ROOT/PINS["beneficiary"][0])
 butterflies=sorted(df.species.unique())
 islands=sorted(df.landmass_id.unique(),key=int)
 if len(df)!=11475 or len(butterflies)!=225 or len(islands)!=51 or len(geo)!=51 or len(benefits)!=239:
  raise ValueError("Original GBIF source dimensions changed")
 if set(islands)!=set(geo.island_id) or not df.set_index(["species","landmass_id"]).index.is_unique:
  raise ValueError("Source island group or unique species island key mismatch")
 if not df.source_not_future_selected.astype(str).str.lower().eq("true").all():
  raise ValueError("Original GBIF reindex selected later records")
 if not df.index_not_true_absence_or_colonization.astype(str).str.lower().eq("true").all():
  raise ValueError("Source absence caveat removed")
 idx=pd.MultiIndex.from_product([butterflies,islands],names=["species","landmass_id"])
 df=df.set_index(["species","landmass_id"]).reindex(idx)
 if df.isna().any().any():
  raise ValueError("Original 225x51 source not dense")
 early=df.pre2001_2010_indexed_records.to_numpy().reshape(225,51)>0
 early_added=(df.future2011_2017_indexed_records.to_numpy().reshape(225,51)>0)&(~early)
 mid=early|early_added
 late_added=(df.future2018_2022_indexed_records.to_numpy().reshape(225,51)>0)&(~mid)
 future=early_added|late_added
 later=early|future
 if (int(early.sum()),int(early_added.sum()),int(late_added.sum()),int(later.sum()))!=(736,221,96,1053):
  raise ValueError("Original 2001-2010 and later indexing totals changed")
 if np.any(early&future):raise ValueError("Non-first indexing misclassified")
 labels=geo.set_index("island_id").loc[islands].longitude_block.to_numpy(int)
 if sorted(np.bincount(labels)) != [10,10,10,10,11]:raise ValueError("Original outcome-blind geography population changed")
 j0,valid=jaccard(early)
 j1,v1=jaccard(mid)
 j2,v2=jaccard(later)
 if not np.array_equal(valid,v1) or not np.array_equal(valid,v2) or int(valid.sum())!=1275:
  raise ValueError("Original early Jaccard pair support changed")
 xy=np.triu_indices(51,1)
 within=labels[xy[0]]==labels[xy[1]]
 scopes={"all_51_original_landmass_pairs":np.ones(len(valid),dtype=bool),
         "within_original_longitude_groups":within,
         "between_original_longitude_groups":~within}
 observed={scope:summary(j0,j2,valid,mask) for scope,mask in scopes.items()}
 observed["time_series"]={"indexed_2001_2010":{"species_island_cells":736,"mean_Jaccard":float(np.nanmean(j0))},
                          "through_2017":{"species_island_cells":957,"mean_Jaccard":float(np.nanmean(j1))},
                          "through_2022":{"species_island_cells":1053,"mean_Jaccard":float(np.nanmean(j2))},
                          "new_2011_2017":221,"new_2018_2022":96,
                          "original_unique_species_island_cells":11475}
 nulls={}
 for regional in [False,True]:
  name="source_fixed_butterfly_island_new_index_margins"+("_plus_each_butterfly_longitude_group" if regional else "")
  samples,details=first_index_null(future,early,labels,np.random.default_rng(SEED+int(regional)),regional)
  nulls[name]={"mixing":details,"evaluations":{
     scope:measure_null(observed[scope]["mean_increment"],j0,samples,valid&mask)
     for scope,mask in scopes.items()}}
 # Exact fixed resource quartiles; this is an independent descriptive corroboration,
 # NOT a matched 126-region resource-to-51-landmass butterfly causality test.
 benefit=benefits.set_index("species").loc[butterflies]
 index_quartiles=[]
 for q in (1,2,3,4):
  sel=(benefit.native_host_footprint_quartile.to_numpy(int)==q)
  n=int(sel.sum())
  available=int((~early[sel]).sum())
  new=int(future[sel].sum())
  index_quartiles.append({
    "native_host_plant_resource_footprint_quartile":q,
    "matched_butterfly_species":n,
    "sampled_prior_2001_2010_species_island_indices":int(early[sel].sum()),
    "at_risk_unindexed_species_island_pairs":available,
    "new_2011_2022_species_island_indices":new,
    "share_all_317_new_butterfly_island_indices":float(new/317),
    "recorded_new_index_rate_among_at_risk_cells":float(new/available),
    "island_host_resource_opportunities_added_in_DIFFERENT_126_biogeographic_regions":int(benefit.loc[sel,"added_island_regions"].sum())
   })
 if sum(x["new_2011_2022_species_island_indices"] for x in index_quartiles)!=317:
  raise ValueError("New indexed species lost in quartile link")
 index_rho=spearmanr(benefit.added_island_regions.to_numpy(),future.sum(axis=1))
 pre_index_rho=spearmanr(benefit.native_resource_units.to_numpy(),early.sum(axis=1))
 # Rank partial relative association of source-fixed added HOSTS geography with
 # actual later recording indices, controlling cataloged known-host count,
 # native host footprint and earlier observed butterfly 51-landmass breadth.
 def rankresid(v,C):
  R=rankdata(v)
  rc=np.column_stack([np.ones(len(v))]+[rankdata(z) for z in C])
  return R-rc@np.linalg.lstsq(rc,R,rcond=None)[0]
 controls=[benefit.native_resource_units.to_numpy(),benefit.resolved_host_species.to_numpy(),early.sum(axis=1)]
 t=rankresid(future.sum(axis=1),controls)
 p=rankresid(benefit.added_island_regions.to_numpy(),controls)
 adjusted_rho=float(np.corrcoef(t,p)[0,1])
 # Descriptive 225-butterfly bootstrap for quartile first-index difference.
 qarr=benefit.native_host_footprint_quartile.to_numpy(int)
 obs_rate=index_quartiles[3]["recorded_new_index_rate_among_at_risk_cells"]-index_quartiles[0]["recorded_new_index_rate_among_at_risk_cells"]
 rng=np.random.default_rng(SEED+19)
 band=[]
 for _ in range(4999):
  vals=[]
  for q in (1,4):
   candidate=np.flatnonzero(qarr==q)
   take=candidate[rng.integers(0,len(candidate),len(candidate))]
   b=future[take].sum()
   trials=(~early[take]).sum()
   vals.append(float(b/trials))
  band.append(vals[1]-vals[0])
 result={
  "schema":"chocho_cumulative_GBIF_butterfly_225x51_2001_2010_2017_2022_source_locked_margin_control_v1",
  "status":"COMPUTATION_SOURCE_VERIFIED_BUTTERFLY_INDEXING_NOT_TRUE_HOMOGENIZATION",
  "source_git_blobs":pins,
  "source_periods":{"early":"2001-2010","mid":"2011-2017","late":"2018-2022"},
  "cohort":{"butterflies":225,"fixed_landmass_groups":51,"species_island_cells":11475,
            "original_pairwise_dyads":1275,
            "earlier_indexed_cells":736,"new_first_indexed_2011_2017":221,
            "new_first_indexed_2018_2022":96,
            "new_first_indexed_2011_2022":317},
  "cumulative_recording_Jaccard":observed,
  "constrained_recording_nulls":nulls,
  "resource_beneficiary_correspondence":{
    "breadth_type":"Native distribution of known larval host plants, NOT native butterfly range or endemicity",
    "quartiles":index_quartiles,
    "earlier_butterfly_index_breadth_vs_native_HOSTS_breadth_spearman":float(pre_index_rho.statistic),
    "observed_new_index_count_vs_126_island_resource_additions_spearman":float(index_rho.statistic),
    "partial_rank_new_index_vs_added_host_resource_after_native_host_breadth_HOSTS_size_early_butterfly_index":adjusted_rho,
    "Q4_minus_Q1_at_risk_first_index_rate_pp":float(100*obs_rate),
    "Q4_minus_Q1_species_bootstrap4999_descriptive95_pp":[float(100*x) for x in np.quantile(band,[.025,.975])]},
  "critical_boundaries":[
    "Cumulative GBIF index is observed detection, not true historical butterfly occurrence, colonization or abundance; zero indexing not absence",
    "225 butterfly species and 51 nonrandom WPI-linked landmass groups differ from 239 butterflies and 126 WGSRPD3 floristic regions",
    "All 317 post-2010 first GBIF indexed butterfly-island cells are fixed in null species and island margins; uneven observer effort not fully reconstructed",
    "Longitude-5 geographic group null is a sensitivity not actual biogeographic historic butterfly source range",
    "Dependent Markov-chain null configurations do not provide independent biological replicates or confirmatory p values",
    "Resource quartiles reflect native HOSTS plant geographic range, not native or endemic butterfly distributions, and HOSTS host count confounds breadth",
    "No dated local introduction or verified local larval survival; the resource layer and observed GBIF index not a causal matched-pair panel",
    "Never promote an increase in indexed species overlap alone to actual butterfly biotic homogenization, ecological competitive replacement, ship transfer or realized establishment"]
 }
 args.out_dir.mkdir(parents=True,exist_ok=True)
 (args.out_dir/"CUMULATIVE_GBIF_BUTTERFLY_ISLAND_INDEX_SIMILARITY_20261117.json").write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
 pd.DataFrame(index_quartiles).to_csv(args.out_dir/"butterfly_early_index_and_later_index_by_native_host_resource_quartile.csv",index=False)
 pd.DataFrame([{"comparison":k,**v} for k,v in observed.items() if k!="time_series"]).to_csv(args.out_dir/"cumulative_recorded_butterfly_Jaccard_by_geographic_group.csv",index=False)
 print("ORIGINAL_BUTTERFLY_CUMULATIVE_INDEX_RESULT="+json.dumps({
    "n":result["cohort"],"jaccard":observed,"nulls":{k:{scope:z["extra_beyond_null"] for scope,z in v["evaluations"].items()} for k,v in nulls.items()},
    "benefit":result["resource_beneficiary_correspondence"]},ensure_ascii=False,sort_keys=True),flush=True)

if __name__=="__main__":
 p=argparse.ArgumentParser()
 p.add_argument("--out-dir",type=Path,required=True)
 analyze(p.parse_args())
