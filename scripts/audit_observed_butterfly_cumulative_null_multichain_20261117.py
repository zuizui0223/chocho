#!/usr/bin/env python3
"""Four independent chains to stress-test descriptive 225x51 GBIF record
indexing convergence beyond constrained addition margins.

NO native butterfly range, establishment or insect-population inference.
Only original GBIF 225×51 source and preselected five longitude groups.
The original single-chain 399-draw experiment remains untouched; this is a
diagnostic for dependence/mixing, not new independent ecological replicates.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
from audit_observed_butterfly_cumulative_index_homogenization_20261117 import source,first_index_null,jaccard,summary,SEED,PINS

CHAINS=4
DRAWS_PER_CHAIN=199
MAX_PREDEFINED_ACROSS_CHAIN_MEAN_SPREAD=.004

def main():
 ap=argparse.ArgumentParser()
 ap.add_argument("--out-dir",type=Path,required=True)
 a=ap.parse_args()
 receipts=source()
 src=pd.read_csv(PINS["butterflies"][0][0],dtype={"species":str,"landmass_id":str})
 geo=pd.read_csv(PINS["longitude"][0][0],dtype={"island_id":str})
 names=sorted(src.species.unique())
 islands=sorted(src.landmass_id.unique(),key=int)
 key=pd.MultiIndex.from_product((names,islands),names=("species","landmass_id"))
 d=src.set_index(["species","landmass_id"]).reindex(key)
 if len(d)!=225*51 or d.isna().any().any():raise ValueError("Original public225x51 input changed")
 early=d.pre2001_2010_indexed_records.to_numpy().reshape(225,51)>0
 later=((d.future2011_2017_indexed_records.to_numpy().reshape(225,51)>0) |
        (d.future2018_2022_indexed_records.to_numpy().reshape(225,51)>0))
 new=later&(~early)
 if (int(early.sum()),int(new.sum()))!=(736,317):raise ValueError("Original early/new indexing scope changed")
 groups=geo.set_index("island_id").loc[islands].longitude_block.to_numpy(dtype=int)
 j0,valid=jaccard(early)
 j2,_=jaccard(early|new)
 orig=float(np.nanmean(j2-j0))
 if abs(orig-.05680961963536349)>1e-12:raise ValueError("Original 51 pair-Jaccard increase not reproduced")
 i,j=np.triu_indices(51,1)
 mask={"all_pairs":valid,"within_longitude":valid&(groups[i]==groups[j]),"between_longitude":valid&(groups[i]!=groups[j])}
 observed={scope:float(np.mean((j2-j0)[v])) for scope,v in mask.items()}
 mode={}
 for regional in (False,True):
  name="geography_preserving" if regional else "species_and_island_margins"
  chain_rows=[]
  full={}
  for k in range(CHAINS):
   rng=np.random.default_rng(SEED+10000*(k+1)+int(regional))
   samples,receipt=first_index_null(new,early,groups,rng,regional,draws=DRAWS_PER_CHAIN)
   if receipt["source_margin_checks"]!="PASS":raise ValueError("Margin preservation lost")
   measures={scope:np.mean(samples[:,ind]-j0[ind][None,:],axis=1) for scope,ind in mask.items()}
   chain_rows.append({"chain":k+1,"seed":SEED+10000*(k+1)+int(regional),
                      "acceptance_fraction":receipt["swap_acceptance_fraction"],
                      **{scope+"_mean":float(z.mean()) for scope,z in measures.items()},
                      **{scope+"_sd":float(z.std(ddof=1)) for scope,z in measures.items()}})
   for scope,vals in measures.items():
    full.setdefault(scope,[]).extend(float(x) for x in vals)
  pooled={}
  for scope,arr in full.items():
   vals=np.asarray(arr)
   means=np.asarray([row[scope+"_mean"] for row in chain_rows])
   gap=float(means.max()-means.min())
   pooled[scope]={"observed_increment":observed[scope],
     "null_mean":float(vals.mean()),"null_descriptive_2p5_97p5":[float(x) for x in np.quantile(vals,[.025,.975])],
     "observed_excess":float(observed[scope]-vals.mean()),
     "chain_mean_spread":gap,
     "predeclared_mixing_stability_threshold":MAX_PREDEFINED_ACROSS_CHAIN_MEAN_SPREAD,
     "chain_means_stable":bool(gap<MAX_PREDEFINED_ACROSS_CHAIN_MEAN_SPREAD),
     "n_chain_draws_total":len(vals),
     "status":"EXPLORATORY_NOT_INDEPENDENT_P_VALUES"}
  mode[name]={"chains":chain_rows,"diagnostic":pooled}
 out={"schema":"chocho_225x51_early_to_late_record_jaccard_fixed_margins_four_chain_stability_v1",
      "status":"FOUR_CHAIN_DESCRIPTIVE_SOURCE_LOCKED_OBSERVATION_NOT_COLONISATION",
      "source_pins":receipts,"source_island_group":"five original 2026-11-07 outcome-blind longitude bins",
      "seed_base":SEED,"chains_per_null":CHAINS,"draws_per_chain":DRAWS_PER_CHAIN,
      "original_indexed_2001_2010":736,"post2010_new_first_indexed":317,
      "same_original_51_islands_pairwise_Jaccard_change":observed,
      "nulls":mode,
      "limitations":[
       "A separate 4-chain mixing diagnostic, not four independent ecological cohorts or source tests",
       "Scored on the same originally selected 51 WPI port-linked island groups",
       "Preserves butterfly and island recorded new-species richness, not complete observer effort or source-list bias",
       "Within-longitude-group constraint preserves taxon-year region allocation but latitude, island area, historical observation effort can still confound",
       "Cumulative GBIF first indexing is NOT true butterfly biotic homogenization, colony formation, abundance change or evidence of plant-mediated range shift",
       "Fixed source first-record indexing is vulnerable to historical re-identification; modern biodiversity data not assumed complete"] }
 a.out_dir.mkdir(parents=True,exist_ok=True)
 (a.out_dir/"MULTICHAIN_225_51_RECORD_SIMILARITY_SOURCE_NULL_20261117.json").write_text(json.dumps(out,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
 print("FOUR_CHAIN_SOURCE_NULL_VERIFICATION="+json.dumps({"observed":observed,"nulls":{n:v["diagnostic"] for n,v in mode.items()}},ensure_ascii=False,sort_keys=True),flush=True)

if __name__=="__main__":
 main()
