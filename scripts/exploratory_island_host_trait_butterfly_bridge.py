#!/usr/bin/env python3
"""Butterfly end of proposed island port -> introduced plant -> butterfly bridge.

SOURCE LOCKS:
- Q2 v02: exactly 11,573 confirmed *potential-resource* butterfly×introduced
  plant×strict island-region edges; 5,741 unique species×island opportunity cells.
- WCVP13/GIFT3.2 exact accepted plant ID -> life-cycle and growth form crosswalk
  from source-checked run 41946522894.
- Original frozen 239 butterfly panel resource-gain data + LepTraits disturbance.
No new external requests. Never infer observed feeding, plant introduction route,
nutrient-site preference or butterfly reproduction from potential resources.

Test A: Do disturbance-associated butterflies have a higher fraction of
annual(-herb) introduced host plants among the *distinct plants supporting
new island opportunity*?
Test B: Do these host traits explain species-level potential resource gain
after host inventory size, native opportunity and host breadth are controlled?
This is an exploratory association and selection-by-opportunity is material.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd

EXPECTED_EDGES=11573
EXPECTED_CELLS=5741
EXPECTED_INTRO_PLANTS=564
SCORED={"Disturbance-associated (strong)":1, "Disturbance-avoidant (strong)":0}
SEED=20261105


def file_sha(path:Path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_and_check(edges_path:Path,crosswalk_path:Path,panel_path:Path,desc_path:Path):
    edges=pd.read_csv(edges_path,dtype={"plant_id":str,"species":str,"island_code":str})
    fields={"species","island_code","plant_id"}
    if not fields.issubset(edges.columns):raise ValueError("Unexpected Q2 introduced-host source fields")
    if len(edges)!=EXPECTED_EDGES or edges[list(fields)].duplicated().any():
        raise ValueError(f"Q2 introduced-host source edges changed: {len(edges)}")
    if (edges[["species","island_code"]].drop_duplicates().shape[0]!=EXPECTED_CELLS
        or edges.plant_id.nunique()!=EXPECTED_INTRO_PLANTS):
        raise ValueError("Q2 239-butterfly island resource or distinct introduced host counts changed")
    cross=pd.read_csv(crosswalk_path,dtype={"plant_name_id":str},keep_default_na=False)
    if not {"plant_name_id","taxon_key","lifecycle","growth_form"}.issubset(cross.columns):
        raise ValueError("WCVP×GIFT source crosswalk missing required columns")
    # Source Q2 edges were exported by R as integer IDs, whereas the
    # pyreadr/pandas GIFT crosswalk may export numeric R IDs with '.0'.
    # Normalize only lossless integer-like WCVP keys; never fuzzy-match taxa.
    def id_key(v):
        t=str(v).strip()
        if re.fullmatch(r"[0-9]+(?:\.0+)?",t):
            return t.split(".")[0]
        return t
    # Python source code must use a single slash in the regex; verify first.
    cross["plant_name_id"]=cross.plant_name_id.map(id_key)
    edges["plant_id"]=edges.plant_id.map(id_key)
    if cross.plant_name_id.duplicated().any():
        raise ValueError("Trait crosswalk has conflicting accepted WCVP IDs")
    ids=set(edges.plant_id)
    present=ids.intersection(set(cross.plant_name_id))
    print("Q2_TO_GIFT_WCVP_ACCEPTED_ID_JOIN",
          json.dumps({"q2_distinct_ids":len(ids),"matched_trait_crosswalk":len(present),
                      "q2_sample":sorted(ids)[:12],
                      "crosswalk_sample":cross.plant_name_id.head(12).to_list()},sort_keys=True),flush=True)
    if len(present)<250:
        raise ValueError("Too few exact accepted WCVP IDs crossmapped to Q2 introduced hosts; do not infer no known traits")
    cross["known_lifecycle"]=cross.lifecycle.isin(("annual","biennial","perennial"))
    cross["known_both"]=cross["known_lifecycle"]&cross.growth_form.isin(("herb","shrub","tree"))
    cross["annual"]=cross.lifecycle.eq("annual")
    cross["annual_herb"]=cross["annual"]&cross.growth_form.eq("herb")
    candidates=edges[["plant_id"]].drop_duplicates().merge(cross,left_on="plant_id",
        right_on="plant_name_id",how="left",validate="one_to_one")
    for flag in ("known_lifecycle","known_both","annual","annual_herb"):
        # NA means unknown taxon–trait match; never recode as a negative fact
        # in denominators. Explicit flags for summary numerator only.
        candidates[flag]=candidates[flag].fillna(False).astype(bool)
    plant_mapping=candidates.set_index("plant_id")[
        ["known_lifecycle","known_both","annual","annual_herb"]]
    expanded=edges[["species","island_code","plant_id"]].merge(
        plant_mapping,left_on="plant_id",right_index=True,
        how="left",validate="many_to_one")
    if expanded[["known_lifecycle","known_both"]].isna().any().any():
        raise ValueError("Unaccounted missing trait classification mapping")
    panel=pd.read_csv(panel_path)
    desc=pd.read_csv(desc_path)
    if len(panel)!=239 or panel.species.nunique()!=239:
        raise ValueError("Frozen 239 butterfly scientific panel changed")
    if desc.species.duplicated().any():
        raise ValueError("Species descriptor duplication")
    if edges.species.isin(panel.species).sum()!=len(edges):
        raise ValueError("Q2 edges include butterflies absent from frozen panel")
    return expanded,candidates,panel,desc


def plant_level_species_summary(expanded,panel,desc):
    def summarize(g,unit):
        total=len(g)
        life=int(g.known_lifecycle.sum())
        both=int(g.known_both.sum())
        return {"total":int(total),
            "known_lifecycle":life,
            "known_both":both,
            "annual_n":int((g.annual & g.known_lifecycle).sum()),
            "annual_herb_n":int((g.annual_herb & g.known_both).sum()),
            "annual_fraction":float((g.annual & g.known_lifecycle).sum()/life) if life else np.nan,
            "annual_herb_fraction":float((g.annual_herb & g.known_both).sum()/both) if both else np.nan,
            "lifecycle_coverage":life/total if total else np.nan,
            "both_trait_coverage":both/total if total else np.nan,
            "unit":unit}
    rows=[]
    for sp,g in expanded.groupby("species",sort=True):
        uniq=g.drop_duplicates("plant_id")
        a=summarize(uniq,"distinct_introduced_host_taxa")
        b=summarize(g,"butterfly_introduced_host_island_edges")
        rows.append({"species":sp,"genus":sp.split()[0],
                     "island_added_cells":int(g.island_code.nunique()) if False else
                       int(g[["species","island_code"]].drop_duplicates().shape[0]),
                     **{"distinct_"+k:v for k,v in a.items()},
                     **{"edge_"+k:v for k,v in b.items()}})
    joined=pd.DataFrame(rows).merge(panel,on="species",how="left",
                validate="one_to_one").merge(desc[["species","disturbance_affinity"]],
                on="species",how="left",validate="one_to_one")
    joined["disturbance_strong_binary"]=joined.disturbance_affinity.map(SCORED)
    joined["log_native"]=np.log1p(joined.native_resource_units)
    joined["log_host_families"]=np.log1p(joined.host_family_count)
    joined["log_introduced_host_count"]=np.log1p(joined.distinct_total)
    joined["log_gain"]=joined.log_resource_expansion
    if len(joined)>239:
        raise ValueError("Source Q2 expanded butterfly count exceeds S1 panel")
    return joined


def adjusted_beta(data,response,predictor,controls):
    fields=[predictor,*controls]
    d=data[["genus",response,*fields]].dropna().copy()
    if len(d)<15:return None
    if d[predictor].nunique()<2:return None
    X=np.column_stack([np.ones(len(d)),*[d[x].to_numpy(dtype=float) for x in fields]])
    y=d[response].to_numpy(dtype=float)
    if np.linalg.matrix_rank(X)<X.shape[1]:
        return None
    coef=np.linalg.lstsq(X,y,rcond=None)[0]
    resid=y-X@coef
    return {"n":int(len(d)),"genera":int(d.genus.nunique()),"beta":float(coef[1]),
            "R2":float(1-np.sum(resid*resid)/np.sum((y-y.mean())**2)),
            "source_estimator":"species_level_OLS_with_genus_cluster_bootstrap",
            "predictor":predictor,"response":response,"controls":controls}


def genus_bootstrap(data,response,predictor,controls,*,B=1999,seed=SEED):
    original=adjusted_beta(data,response,predictor,controls)
    if original is None:return {"status":"UNIDENTIFIABLE_FIXED_SOURCE"}
    d=data[["genus",response,predictor,*controls]].dropna()
    groups=[g.copy() for _,g in d.groupby("genus",sort=True)]
    rng=np.random.default_rng(seed)
    vals=[]
    for _ in range(B):
        ix=rng.integers(0,len(groups),len(groups))
        take=pd.concat([groups[i] for i in ix],ignore_index=True)
        val=adjusted_beta(take,response,predictor,controls)
        if val is not None:vals.append(val["beta"])
    if len(vals)<B*.95:return {"status":"TOO_MANY_BOOTSTRAP_FAILURES",
                               "success":len(vals),"attempts":B,"original":original}
    v=np.asarray(vals)
    out={**original,
        "status":"EXPLORATORY_NONCAUSAL_GENUS_BOOTSTRAP",
        "genus_bootstrap_draws":B,"successes":len(vals),
        "beta_boot95":list(map(float,np.quantile(v,[.025,.975]))),
        "boot_positive_fraction":float((v>0).mean())}
    return out


def analyze(expanded,candidates,panel,desc,*,B=1999):
    species=plant_level_species_summary(expanded,panel,desc)
    plants_total=len(candidates)
    known_both=int(candidates.known_both.sum())
    known_life=int(candidates.known_lifecycle.sum())
    if known_life<100 or known_both<100:
        raise ValueError(f"Insufficient trait-matched introduced host taxa: {known_life=} {known_both=}; source join/trait gate requires investigation")
    scored=species[species.disturbance_strong_binary.notna()].copy()
    if len(scored)<60:raise ValueError("Insufficient disturbance-trait butterflies")
    controls=["log_introduced_host_count","log_native","log_host_families",
              "distinct_both_trait_coverage"]
    results={"design":{
       "source_intro_plant_island_edges":len(expanded),
       "source_species_island_added_cells":expanded[["species","island_code"]].drop_duplicates().shape[0],
       "distinct_introduced_hosts":plants_total,
       "plant_lifecycle_known":known_life,
       "plant_lifecycle_plus_growth_known":known_both,
       "distinct_butterflies_with_added_island_resource":len(species),
       "disturbance_status_scored_butterflies_with_added_island_resource":len(scored),
       "disturbance_associated_strong":int(scored.disturbance_strong_binary.eq(1).sum()),
       "disturbance_avoidant_strong":int(scored.disturbance_strong_binary.eq(0).sum()),
       "unexpanded_species_excluded_not_coded_zero_host_trait":int(len(panel)-len(species))},
       "predeclared_primary":{
           "trait":"GIFT annual+herb fraction among DISTINCT introduced host plants used in new island opportunity",
           "direction":"strong disturbance-associated > strong disturbance-avoidant",
           "minimum_known_host_species":2,"minimum_coverage":.30,
           "not_perfect_diagnostic":"These host plants *already produced binary resource opportunity* in source data; not independent evidence of real colonization",
       },
       "cohorts":{},
       "limitations":[
          "The identified introduced hosts are selected *because they create potential added island resource cells*; test is conditional on observed reconstructed gain, not population-level interaction mechanism.",
          "Life-cycle annual or herbaceous habit does not imply high nitrogen requirement, ruderal preference, host palatability or butterfly feeding.",
          "Source GIFT coverage varies by plant genus and region; unscored host taxa are omitted from trait fraction denominators, not assigned negative.",
          "Disturbance affinity labels derive from LepTraits source groups, not quantitative butterfly fecundity or actual urban adaptation.",
          "Genus cluster bootstrap does not handle host-plant reuse among butterflies or plant phylogeny.",
          "Controlling total host taxon richness and native resource breadth is exploratory adjustment; no causal mediation is estimated.",
          "Butterfly host taxa may be used by multiple butterfly species and many islands; host-plant trait prevalence can be nonindependent.",
          "No positive evidence for port enrichment of *introduced* annual/herb flora was found in previous primary all-flora test.",
       ]}
    # Primary: compare butterflies with >=2 known plants and >=30% coverage.
    for label,minknown,cov,adequate in [
        ("primary_2plants_30pct",2,.30,False),
        ("min_3plants_50pct",3,.50,False),
        ("min_2plants_50pct",2,.50,False),
        ("source_host_taxonomy_adequate_only",2,.30,True),
        ("all_covered_any_scored",1,0,False),
    ]:
        sub=scored[(scored.distinct_known_both>=minknown)&
                   (scored.distinct_both_trait_coverage>=cov)].copy()
        if adequate:
            sub=sub[sub.host_taxonomy_lower_bound_adequate.astype(str).str.casefold().eq("true")]
        group=sub.groupby("disturbance_strong_binary").distinct_annual_herb_fraction.agg(["count","mean"])
        res={"n":len(sub),"genera":sub.genus.nunique(),
             "disturbance_group_n":sub.disturbance_strong_binary.value_counts().to_dict(),
             "annual_herb_means_by_group":group.to_dict(orient="index"),
             "source_trait_coverage_median":float(sub.distinct_both_trait_coverage.median()) if len(sub) else None}
        if len(sub)>=30 and sub.disturbance_strong_binary.nunique()==2:
            res["disturbance_to_annual_herb_fraction"]=genus_bootstrap(
              sub,"distinct_annual_herb_fraction","disturbance_strong_binary",controls,B=B,seed=SEED)
            res["trait_to_potential_gain_adjusting_disturbance"]=genus_bootstrap(
              sub,"log_gain","distinct_annual_herb_fraction",
              ["disturbance_strong_binary",*controls],B=B,seed=SEED+1)
            res["edge_weighted_trait_sensitivity"]=genus_bootstrap(
              sub,"edge_annual_herb_fraction","disturbance_strong_binary",
              ["log_introduced_host_count","log_native","log_host_families",
               "edge_both_trait_coverage"],B=B,seed=SEED+2)
        else:res["status"]="INSUFFICIENT_GROUPS"
        results["cohorts"][label]=res
        print("ISLAND_HOST_TRAIT_TO_BUTTERFLY_DISTURBANCE",label,json.dumps({
          "n":res["n"],"group_n":res["disturbance_group_n"],
          "means":res["annual_herb_means_by_group"],
          "primary":res.get("disturbance_to_annual_herb_fraction"),
          "gain":res.get("trait_to_potential_gain_adjusting_disturbance"),
          "edge_weighted":res.get("edge_weighted_trait_sensitivity")
          },sort_keys=True)[:14000],flush=True)
    return species,results


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for v in ("q2-edges","wcvp-gift-crosswalk","panel","descriptors","out","species-out"):
        p.add_argument("--"+v,type=Path,required=True)
    p.add_argument("--bootstrap",type=int,default=1999)
    a=p.parse_args()
    q,plants,panel,desc=read_and_check(a.q2_edges,a.wcvp_gift_crosswalk,a.panel,a.descriptors)
    s,results=analyze(q,plants,panel,desc,B=a.bootstrap)
    summary={
      "schema":"chocho_20261105_island_introduced_host_GIFT_annual_herb_butterfly_disturbance_bridge_v1",
      "status":"EXPLORATORY_NOT_CAUSAL_BUTTERFLY_COLONIZATION_OR_NITROGEN_SELECTION",
      "source":{
         "Q2_github_action":41398149093,
         "WCVP_GIFT_crosswalk_action":41946522894,
         "239_butterfly_panel_SHA256":file_sha(a.panel),
         "disturbance_source_SHA256":file_sha(a.descriptors),
         "introduced_island_host_edges_SHA256":file_sha(a.q2_edges),
         "GIFT_exact_trait_crosswalk_SHA256":file_sha(a.wcvp_gift_crosswalk)
      },**results}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(summary,indent=2,ensure_ascii=False,default=str)+"\n",encoding="utf-8")
    a.species_out.parent.mkdir(parents=True,exist_ok=True)
    s.to_csv(a.species_out,index=False)
    print("WCVP_GIFT_HOST_BUTTERFLY_TRAIT_BRIDGE_QA",json.dumps(summary["design"],default=str),flush=True)

if __name__=="__main__":main()
