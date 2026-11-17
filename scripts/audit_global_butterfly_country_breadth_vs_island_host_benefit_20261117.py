#!/usr/bin/env python3
"""Does independent pre-2011 GBIF country footprint of butterflies predict who
benefits from human-redistributed HOSTS plant opportunities on islands?

Three distinct and noninterchangeable "geographic breadth" concepts:
 A. native HOSTS-plant resource regions for a butterfly (WCVP, not butterfly range);
 B. GBIF 2001-2010 COORDINATE-POSITIVE countries actually indexed for a butterfly,
    with exact finite-sample 30/100 record country rarefaction, NOT native range;
 C. 2001-2010 GBIF indexed *island* landmass groups, among 51 selected ports.

A and B are independent botanical- and butterfly-sourced geographical footprints.
C and B come from GBIF, but geographic supports differ and record-effort biases
persist. This is not butterfly endemicity, transport causation or colonization.
"""
from __future__ import annotations
import argparse,hashlib,json,math,statistics
from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd
from scipy.stats import spearmanr,rankdata

ROOT=Path(__file__).resolve().parents[1]
FACET_PATH="results/butterfly_global_GBIF_country_breadth_20261117/SOURCE_GBIF2001_2010_225COUNTRY_FACETS_20261117.json"
HOST_PATH="results/butterfly_resource_beneficiaries_20261117/beneficiary_species_native_resource_breadth.csv"
INDEX_PATH="data/exploratory/public51_butterfly_index_20261109/full51_225butterfly_period_index.csv"
FROZEN={
 HOST_PATH:"c9e409441dcd2bb7f9edc5a813fc13490bf60d8f",
 INDEX_PATH:"b6a3c98eb76864052d45d520a55b7d830816126c"
}
SEED=20261117

def gitblob(path):
 b=path.read_bytes()
 return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()

def expected_unique_countries(countries,n=30):
 """Without-replacement exact expected number of country categories in n draws.
 Countries with faceted records 0 not included; records lacking country excluded
 from N; occurrence records may be unevenly sampled among countries.
 """
 N=sum(int(a["indexed_occurrences"]) for a in countries)
 if N<n:return float("nan")
 expct=0.
 for row in countries:
  k=int(row["indexed_occurrences"])
  if k<1:continue
  if N-k<n:
   gone=0.
  else:
   gone=math.exp(math.lgamma(N-k+1)+math.lgamma(N-n+1)-math.lgamma(N-k-n+1)-math.lgamma(N+1))
  expct+=1-gone
 if expct<0 or expct>len(countries)+1.e-4:raise ValueError("Invalid country rarefaction expectation")
 return float(expct)

def residual_rank(v,controls):
 R=rankdata(np.asarray(v,dtype=float))
 C=np.column_stack([np.ones(len(R))]+[rankdata(np.asarray(x,dtype=float)) for x in controls])
 return R-C@np.linalg.lstsq(C,R,rcond=None)[0]

def partial_rho(x,y,control_cols):
 a=residual_rank(x,control_cols)
 b=residual_rank(y,control_cols)
 if a.std()<1.e-10 or b.std()<1.e-10:return float("nan")
 return float(np.corrcoef(a,b)[0,1])

def corr(x,y):
 z=spearmanr(x,y).statistic
 return float(z)

def ridge_fit_predict(X,y,train,test,alpha=10.):
 mu=X[train].mean(axis=0); sd=X[train].std(axis=0)
 sd=np.where(sd>1.e-9,sd,1)
 z=(X[train]-mu)/sd
 w=(X[test]-mu)/sd
 slope=np.linalg.solve(z.T@z+alpha*np.eye(z.shape[1]),z.T@(y[train]-y[train].mean()))
 return y[train].mean()+w@slope

def run_oof(rows,full,include,alpha=10.):
 d=rows.copy()
 # Species-level random fixed splits and genus-block holdout are both
 # predetermined by names, independent of the response values.
 names=d.species.astype(str).tolist()
 foldA=np.array([int(hashlib.sha256(("chocho20261117_"+x).encode()).hexdigest()[:8],16)%5 for x in names])
 foldB=np.array([int(hashlib.sha256(("genus20261117_"+x.split()[0]).encode()).hexdigest()[:8],16)%5 for x in names])
 target=np.log1p(d.added_island_regions.to_numpy(float))
 base=[np.log1p(d.native_resource_units.to_numpy(float)),
       np.log1p(d.native_island_regions.to_numpy(float)),
       np.log1p(d.resolved_host_species.to_numpy(float)),
       np.log1p(d.host_family_count.to_numpy(float)),
       np.log1p(d.n_country_facet_records.to_numpy(float)),
       np.log1p(d.pre2011_51_island_indexed_count.to_numpy(float))]
 bX=np.column_stack(base)
 cX=np.column_stack(base+[d.rarefied_30_countries.to_numpy(float)])
 allscores={}
 for key,fold in [("species_stable_hash_5fold",foldA),("genus_block_5fold",foldB)]:
  pred0=np.full(len(d),np.nan);pred1=np.full(len(d),np.nan)
  sizes=[]
  for k in range(5):
   test=fold==k;train=~test
   sizes.append({"fold":k,"held_species":int(test.sum()),"train_species":int(train.sum())})
   if test.sum()<2 or train.sum()<30:raise ValueError("Insufficient genus/species holdout")
   pred0[test]=ridge_fit_predict(bX,target,train,test,alpha)
   pred1[test]=ridge_fit_predict(cX,target,train,test,alpha)
  if not (np.isfinite(pred0).all() and np.isfinite(pred1).all()):raise ValueError("Invalid OOF prediction")
  mse0=float(np.mean((target-pred0)**2))
  mse1=float(np.mean((target-pred1)**2))
  allscores[key]={"folds":sizes,"baseline_MSE_log1p_added_island_resources":mse0,
                 "add_rarefied_30_countries_MSE":mse1,
                 "relative_MSE_reduction_pct":100*(mse0-mse1)/mse0,
                 "response":"log1p(number of new potential butterfly host-covered botanical islands)",
                 "baseline_controls":["log1p(native HOSTS resource regions)","log1p(native HOSTS island regions)",
                     "log1p(HOSTS known plant species)","log1p(host family count)",
                     "log1p(GBIF 2001-2010 georeferenced country facet records)",
                     "log1p(earlier indexed butterflies on the original 51 sample islands)"],
                 "added_predictor":"Expected distinct country categories in 30 independent-country-facet observation records, exact finite-sample without replacement",
                 "ridge_alpha_fixed":alpha}
 return allscores

def main():
 p=argparse.ArgumentParser();p.add_argument("--out-dir",type=Path,required=True);a=p.parse_args()
 original={}
 for name,pin in FROZEN.items():
  digest=gitblob(ROOT/name)
  if digest!=pin:raise ValueError("Frozen upstream source changed "+name)
  original[name]=pin
 facet_file=ROOT/FACET_PATH
 if not facet_file.exists():raise FileNotFoundError("Original country-facet source artifact not yet committed "+FACET_PATH)
 facets=json.loads(facet_file.read_text())
 if (facets["requested"],facets["completed"],facets["failed"])!=(225,225,0) or facets["pilot"] or facets["status"]!="ALL_REQUESTS_COMPLETE":
  raise ValueError("GBIF 225 butterfly country facet source is not COMPLETE")
 if facets["fixed_parameters"]["year"]!="2001,2010" or facets["fixed_parameters"]["hasCoordinate"]!="true":
  raise ValueError("Invalid year/georeference GBIF facet source")
 benefits=pd.read_csv(ROOT/HOST_PATH)
 records=pd.read_csv(ROOT/INDEX_PATH,dtype={"species":str,"landmass_id":str})
 if len(benefits)!=239 or len(records)!=225*51:raise ValueError("Source cohorts changed")
 b=benefits.set_index("species")
 if b.index.has_duplicates:raise ValueError("Duplicated HOSTS resource butterfly")
 record_early=records.assign(indexed=records.pre2001_2010_indexed_records.gt(0)).groupby("species").indexed.sum().to_dict()
 record_future=records.assign(new=(records.first_indexed_in_2011_2017.eq(1)|records.first_indexed_in_2018_2022.eq(1))).groupby("species").new.sum().to_dict()
 if sum(record_future.values())!=317:raise ValueError("Original public 317 future indices changed")
 results=[]
 for r in facets["countries_per_butterfly"]:
  s=r["species"];key=r["speciesKey"]
  if s not in b.index or s not in record_early:raise ValueError("GBIF butterfly not in frozen original resources or early observations")
  if r["status"]!="COMPLETE_GBIF_COUNTRY_FACETS":raise ValueError("Unresolved source failure for "+s)
  countries=r["country_facets"]
  n_country=sum(z["indexed_occurrences"] for z in countries)
  if n_country!=r["country_records_in_facets"]:raise ValueError("Source country facet sum mismatch")
  result=b.loc[s].to_dict()
  results.append({"species":s,"speciesKey":key,
    "n_country_facet_records":n_country,
    "n_GBIF_period_records_total":r["GBIF_count_indexed_2001_2010"],
    "observed_countries_any":r["n_countries_any_country_facet"],
    "observed_countries_5plus":r["n_countries_5plus_records"],
    "observed_countries_10plus":r["n_countries_10plus_records"],
    "rarefied_30_countries":expected_unique_countries(countries,30),
    "rarefied_100_countries":expected_unique_countries(countries,100),
    "native_resource_units":int(result["native_resource_units"]),
    "native_island_regions":int(result["native_island_regions"]),
    "added_island_regions":int(result["added_island_regions"]),
    "fraction_native_island_gaps_filled":float(result["fraction_native_island_gaps_filled"]),
    "resolved_host_species":int(result["resolved_host_species"]),
    "host_family_count":float(result["host_family_count"]),
    "source_native_HOSTS_resource_quartile":int(result["native_host_footprint_quartile"]),
    "pre2011_51_island_indexed_count":int(record_early[s]),
    "later_51_island_first_index_count":int(record_future[s]),
    "later_51_atrisk_cells":int(51-record_early[s])})
 d=pd.DataFrame(results).sort_values("species").reset_index(drop=True)
 if len(d)!=225 or d.species.nunique()!=225:raise ValueError("Original species join missing")
 cohort={"butterfly_species":len(d),
         "pre2011_country_facet_georeferenced_records_total":int(d.n_country_facet_records.sum()),
         "countries_with_2001_2010_any_observation_species":int((d.observed_countries_any>0).sum()),
         "eligible_30_record_rarefaction_species":int(d.rarefied_30_countries.notna().sum()),
         "eligible_100_record_rarefaction_species":int(d.rarefied_100_countries.notna().sum()),
         "GBIF_earlier_51_landmass_indexed_cells":int(d.pre2011_51_island_indexed_count.sum()),
         "later_51_landmass_first_index_cells":int(d.later_51_island_first_index_count.sum()),
         "source_239_unique_species":len(benefits)}
 if (cohort["GBIF_earlier_51_landmass_indexed_cells"],cohort["later_51_landmass_first_index_cells"])!=(736,317):
  raise ValueError("Original index source totals drifted")
 d30=d[d.rarefied_30_countries.notna()].copy()
 d100=d[d.rarefied_100_countries.notna()].copy()
 if len(d30)<80 or len(d100)<50:raise ValueError("Insufficient global GBIF records for defensible breadth sensitivity")
 primary={}
 for label,g,field in [
  ("all_225_country_count",d,"observed_countries_any"),
  ("all_225_country_5plus",d,"observed_countries_5plus"),
  ("rarefied_30_indexed_records",d30,"rarefied_30_countries"),
  ("rarefied_100_indexed_records",d100,"rarefied_100_countries")]:
  x=g[field].to_numpy(float)
  y=g.added_island_regions.to_numpy(float)
  risk=g.fraction_native_island_gaps_filled.to_numpy(float)
  conf=[g.native_resource_units.to_numpy(float),
        g.resolved_host_species.to_numpy(float),
        g.pre2011_51_island_indexed_count.to_numpy(float),
        np.log1p(g.n_country_facet_records.to_numpy(float))]
  primary[label]={"n":len(g),
    "range_proxy":field,
    "Spearman_countrybreadth_vs_island_host_resource_gains":corr(x,y),
    "Spearman_countrybreadth_vs_fraction_of_initial_native_host_gaps_filled":corr(x,risk),
    "partial_rank_with_native_host_footprint_catalogued_HOSTS_count_early_51_index_and_country_observation_effort":partial_rho(x,y,conf),
    "Spearman_countrybreadth_vs_later_51_landmass_first_GBIF_indices":corr(x,g.later_51_island_first_index_count.to_numpy(float)),
    "partial_rank_later_index_count_adjusted_host_footprint_HOSTS_count_early_51_index_country_observation_effort":
      partial_rho(x,g.later_51_island_first_index_count.to_numpy(float),conf)}
 # Independent country-footprint quartiles with no native-host region classification
 # Better use sample-size-standardised rarefied 30-record country breadth.
 d30=d30.sort_values(["rarefied_30_countries","species"]).copy()
 if d30.rarefied_30_countries.nunique()<8:raise ValueError("Insufficient rarefied country breadth resolution")
 d30["observed_butterfly_country_breadth_quartile"]=np.minimum(4,(4*np.arange(len(d30))//len(d30))+1)
 groups=[]
 for quart in (1,2,3,4):
  sub=d30[d30.observed_butterfly_country_breadth_quartile==quart]
  gains=int(sub.added_island_regions.sum())
  events=int(sub.later_51_island_first_index_count.sum())
  atrisk=int(sub.later_51_atrisk_cells.sum())
  groups.append({
   "record_country_breadth_quartile":quart,"species":len(sub),
   "rarefied_30_country_breadth_min":float(sub.rarefied_30_countries.min()),
   "rarefied_30_country_breadth_max":float(sub.rarefied_30_countries.max()),
   "mean_country_breadth_rarefied30":float(sub.rarefied_30_countries.mean()),
   "mean_known_HOSTS_resource_breadth":float(sub.native_resource_units.mean()),
   "mean_known_host_catalog_size":float(sub.resolved_host_species.mean()),
   "new_island_host_resource_opportunities":gains,
   "share_island_added_covered_by_rarefaction_eligible_species":float(gains/d30.added_island_regions.sum()),
   "new_later_51_indexed_species_island_cells":events,
   "fraction_new_51_index_events_in_rarefied_sample":float(events/d30.later_51_island_first_index_count.sum()),
   "new_index_rate_within_early_unindexed_cells":float(events/atrisk)})
 # Bootstrap main ecological sensitivity by butterfly species as descriptive.
 g=d30.copy()
 from scipy.stats import spearmanr
 rng=np.random.default_rng(SEED)
 raw_rho=[];partial=[]
 for _ in range(2999):
  ix=rng.integers(0,len(g),len(g))
  sub=g.iloc[ix]
  v=sub.rarefied_30_countries.to_numpy(float)
  w=sub.added_island_regions.to_numpy(float)
  raw_rho.append(corr(v,w))
  con=[sub.native_resource_units.to_numpy(float),sub.resolved_host_species.to_numpy(float),
       sub.pre2011_51_island_indexed_count.to_numpy(float),
       np.log1p(sub.n_country_facet_records.to_numpy(float))]
  partial.append(partial_rho(v,w,con))
 intervals={
  "rarefied30_raw_Spearman_2999_descriptive95":[float(z) for z in np.quantile(raw_rho,[.025,.975])],
  "rarefied30_effort_and_host_footprint_adjusted_partial_rank_2999_descriptive95":
    [float(z) for z in np.nanquantile(partial,[.025,.975])]}
 # Absolute amount of resource beneficiaries should not follow a rare-country
 # proxy if resource gain independently explained by host-plant distribution and
 # catalog coverage; heldout prediction quantifies rather than guesses this.
 oof=run_oof(d30,full=True,include=True)
 a.out_dir.mkdir(parents=True,exist_ok=True)
 d.to_csv(a.out_dir/"225_butterflies_observed_global_GBIF_country_breadth_and_HOSTS_benefit.csv",index=False)
 d30.to_csv(a.out_dir/"eligible_rarefied_30_butterfly_countries_by_benefit.csv",index=False)
 pd.DataFrame(groups).to_csv(a.out_dir/"independent_GBIF_country_breadth_quartile_beneficiaries.csv",index=False)
 report={
  "schema":"chocho_225_butterfly_independent_observed_country_GBIF_breadth_vs_126_island_HOSTS_resource_benefit_20261117_v1",
  "status":"SOURCE_VERIFIED_INDEPENDENT_BUTTERFLY_OBSERVED_BREADTH_NOT_TRUE_NATIVE_RANGES",
  "source_git_blob_checks":{**original,FACET_PATH:gitblob(facet_file)},
  "GBIF_country_snap_acquired_UTC":facets["obtained_at_UTC"],
  "cohort":cohort,
  "methods":{
   "country_breadth":"Global GBIF occurrence/search, speciesKey and year 2001–2010, coordinates only, country facet, not true native butterflies or island endemicity",
   "sample_effort_correction":"Exact expected unique observed countries in 30 (main) or 100 (sensitivity) records sampled without replacement from each butterfly's original 2001–2010 country occurrence-count vector",
   "rarefied_under_30":"EXCLUDED from rarefied comparison not treated as narrow butterfly",
   "country_observation_5plus":"Separate sensitivity to taxon country reporting count",
   "quartiles":"Equal-species rank groups of sample-size standardized 30-record country breadth. Ties broken lexicographically without using outcomes; small ties may cross groups",
   "partial_rank_controls":["native HOSTS plant-range region count","documented larval host catalogue size",
                            "previous 51-landmass indexed butterfly breadth","log1p(global historical GBIF country facet indexed record count)"],
   "sensitivity_prediction":"Fixed-alpha ridge, OOF 5 butterfly groups by stable SHA256 species-name hashes and independently by genus-group hashes, train-only scaling",
   "bootstrap":"2,999 butterfly-species resamples, descriptive intervals, not phylogenetically independent butterfly species"},
  "primary_correlations":primary,
  "independent_country_breadth_quartiles":groups,
  "descriptive_species_bootstraps":intervals,
  "cross_validated_increment_beyond_native_HOSTS_and_record_effort":oof,
  "interpretation_limits":[
   "GBIF observed-country breadth is INDEPENDENT butterfly recording geography but not verified true native range, endemism, dispersal ability, population abundance or actual island occupancy",
   "GBIF 2001-2010 is a new 2026 retrieval; time-filtering records is not independent source-archival temporal reconstruction, reindex drift can change snapshots",
   "Some zero/low pre-2011 GBIF species have no country coverage: do NOT call them range-restricted; exclude those with <30 records from standardized main",
   "Country sizes and within-country sampling differ; rarefaction adjusts number of records, not recorder locations/habitat/taxonomic coverage or global country border size",
   "The butterfly GBIF country facet and original 51-island butterfly panel share opportunistic GBIF observation processes, despite a distinct global geographic query",
   "The 239 HOSTS-WCVP botanical resource panel and 225 GBIF butterfly panel represent different sampling, no species-level independent plant introduction date or local larval feeding verified",
   "Correlations or predictive improvements do not establish shipborne dispersal, native butterfly community homogenization, or causal butterfly colonization by introduced plants",
   "Source datasetKey/GBIF DOI citation for new country facet queries still needs separate archive protocol before journal submission"]}
 (a.out_dir/"BUTTERFLY_GBIF_COUNTRY_BREADTH_VS_ISLAND_HOST_GAIN_20261117.json").write_text(json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
 print("ACTUAL_BUTTERFLY_COUNTRY_BREADTH_RESULTS="+json.dumps({
   "cohort":cohort,"primary":primary,"country_quartiles":groups,
   "out_of_genus_and_species_fold":oof,"descriptive_intervals":intervals},ensure_ascii=False,sort_keys=True),flush=True)

if __name__=="__main__":
 p=argparse.ArgumentParser();p.add_argument("--out-dir",type=Path,required=True);main_args=p.parse_args();main()
