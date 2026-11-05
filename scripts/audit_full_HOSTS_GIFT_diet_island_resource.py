#!/usr/bin/env python3
"""Independent ALL-documented-HOSTS butterfly diet × GIFT 3.2 plant life-cycle.

Key improvement: butterfly host-plant annual/herb trait fraction is estimated
BEFORE selecting plants that generated newly possible island resource cells.
A source-locked full 239-butterfly HOSTS×WCVP host list defines trait
fractions; the island Q2 source defines a separate 0-inclusive outcome.
No claim about real butterfly colonisation, annuality = nitrophily, or port
transport mechanism. This is an association and a source-quality audit.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from exploratory_island_host_trait_butterfly_bridge import genus_bootstrap

GOOD={"Disturbance-associated (strong)":1.0,
      "Disturbance-avoidant (strong)":0.0}

def id_key(v):
    s=str(v).strip()
    if s.endswith(".0") and s[:-2].isdigit():return s[:-2]
    return s

def load(host_path,cross_path,panel_path,desc_path,edges_path):
    h=pd.read_csv(host_path,dtype={"accepted_plant_name_id":str},
                   keep_default_na=False)
    if not {"insect_species","accepted_plant_name_id","accepted_name","family"}.issubset(h.columns):
        raise ValueError("Not original pinned WCVP–HOSTS accepted source sidecar")
    h["plant_id"]=h.accepted_plant_name_id.map(id_key)
    h=h[["insect_species","plant_id","family"]].drop_duplicates()
    if h.duplicated(["insect_species","plant_id"]).any():
        raise ValueError("HOSTS same insect × plant ID duplicated")

    c=pd.read_csv(cross_path,dtype={"plant_name_id":str},
                  keep_default_na=False)
    c["plant_id"]=c.plant_name_id.map(id_key)
    if c.plant_id.duplicated().any():raise ValueError("Source GIFT accepted WCVP trait key not unique")
    c["known_lifecycle"]=c.lifecycle.isin(["annual","biennial","perennial"])
    c["known_both"]=c.known_lifecycle & c.growth_form.isin(["herb","tree","shrub"])
    c["annual"]=c.lifecycle.eq("annual")
    c["annual_herb"]=c.annual & c.growth_form.eq("herb")

    p=pd.read_csv(panel_path)
    d=pd.read_csv(desc_path)
    e=pd.read_csv(edges_path,dtype={"plant_id":str})
    if len(p)!=239 or p.species.nunique()!=239:
        raise ValueError("Original scientific panel not 239 butterflies")
    if len(e)!=11573 or e[["species","island_code"]].drop_duplicates().shape[0]!=5741:
        raise ValueError("Original corrected Q2 island add source not restored")
    h=h[h.insect_species.isin(p.species)]
    if h.insect_species.nunique()<210:raise ValueError("Full WCVP-host reconstruction lost panel members")
    known_count=len(set(h.plant_id)&set(c.plant_id))
    if known_count<250:
        raise ValueError(f"Full source WCVP-host taxon crosswalk too few accepted IDs: {known_count}")
    # Missing trait IDs remain unavailable; never interpret an unmatched
    # accepted host taxon as a perennial or non-herb.
    resolved=h.merge(c[["plant_id","known_lifecycle","known_both","annual","annual_herb"]],
                     on="plant_id",how="left",validate="many_to_one")
    for col in ("known_lifecycle","known_both","annual","annual_herb"):
        resolved[col]=resolved[col].fillna(False).astype(bool)
    g=resolved.groupby("insect_species",sort=False)
    rows=g.agg(total_hosts=("plant_id","nunique"),
               known_lifecycle=("known_lifecycle","sum"),
               known_both=("known_both","sum"),
               annual=("annual","sum"),
               annual_herb=("annual_herb","sum")).reset_index().rename(columns={"insect_species":"species"})
    rows["annual_fraction"]=rows.annual/rows.known_lifecycle.replace(0,np.nan)
    rows["annual_herb_fraction"]=rows.annual_herb/rows.known_both.replace(0,np.nan)
    rows["coverage"]=rows.known_both/rows.total_hosts
    island=e[["species","island_code"]].drop_duplicates().groupby("species").size().rename("island_added")
    if island.sum()!=5741:raise ValueError("Invalid corrected Q2 island opportunities")
    rows=p.merge(rows,on="species",how="left",validate="one_to_one")
    rows=rows.merge(d[["species","disturbance_affinity"]],on="species",
                    how="left",validate="one_to_one")
    rows["island_added"]=rows.species.map(island).fillna(0).astype(int)
    rows["log_island_added"]=np.log1p(rows.island_added)
    rows["log_global_gain"]=rows.log_resource_expansion
    rows["log_native"]=np.log1p(rows.native_resource_units)
    rows["log_host_families"]=np.log1p(rows.host_family_count)
    rows["log_diet_hosts"]=np.log1p(rows.total_hosts)
    rows["disturbance_binary"]=rows.disturbance_affinity.map(GOOD)
    rows["genus"]=rows.species.str.split().str[0]
    rows["taxonomic_adequacy"]=rows.host_taxonomy_lower_bound_adequate.astype(str).str.casefold().eq("true")
    rowdiag={
      "fixed_butterfly_panel_species":len(rows),
      "butterfly_panel_species_with_resolved_hosts":int(rows.total_hosts.notna().sum()),
      "butterfly_panel_species_with_2plus_scored_plants":int((rows.known_both>=2).sum()),
      "distinct_accepted_host_plants_in_panel":int(h.plant_id.nunique()),
      "distinct_accepted_host_plants_with_any_GIFT_trait":known_count,
      "distinct_accepted_host_plant_island_Q2_source":564,
      "Q2_island_new_opportunity_cells":int(rows.island_added.sum()),
      "butterflies_with_island_added_cells":int((rows.island_added>0).sum()),
      "butterflies_with_zero_island_added_cells":int((rows.island_added==0).sum()),
      "all_host_traits_constructed_without_selecting_positive_island_resource_gain":True}
    return rows, rowdiag, resolved

def model_summaries(rows,B=799):
    out={}
    for name,known,cover,adequate in (
        ("primary_2hosts_30pct",2,.3,False),
        ("tight_3hosts_50pct",3,.5,False),
        ("taxonomy_adequate_2hosts_30pct",2,.3,True),
        ("any_1host_0pct",1,0,False)
    ):
        q=rows[(rows.known_both>=known)&(rows.coverage>=cover)].copy()
        if adequate:q=q[q.taxonomic_adequacy]
        scored=q[q.disturbance_binary.notna()].copy()
        res={"n_all_with_trait":len(q),"n_scored_disturbance":len(scored),
             "genera":scored.genus.nunique(),
             "scored_disturbance_groups":scored.disturbance_binary.value_counts().to_dict(),
             "includes_zero_island_gain":int((scored.island_added==0).sum()),
             "median_GIFT_trait_coverage":float(scored.coverage.median()) if len(scored) else None}
        if len(scored)<25 or scored.disturbance_binary.nunique()<2:
            res["status"]="INSUFFICIENT_SOURCE";out[name]=res;continue
        ctrl=["log_diet_hosts","log_native","log_host_families","coverage"]
        res["disturbance_to_diet_annual_herb"]=genus_bootstrap(
            scored,"annual_herb_fraction","disturbance_binary",ctrl,
            B=B,seed=20261105)
        for response in ("log_island_added","log_global_gain"):
            res["diet_annual_herb_to_"+response]=genus_bootstrap(
                scored,response,"annual_herb_fraction",
                ["disturbance_binary",*ctrl],B=B,seed=20261115)
        # The all species comparator does not require LepTraits disturbance status.
        res["all_butterflies_diet_to_island_cells"]=genus_bootstrap(
             q,"log_island_added","annual_herb_fraction",ctrl,
             B=B,seed=20261125)
        out[name]=res
        print("FULL_HOST_TRAIT_BUTTERFLY_GAIN",name,json.dumps({
            "n":res["n_scored_disturbance"],
            "n_zero":res["includes_zero_island_gain"],
            "disturbance_vs_diet":res["disturbance_to_diet_annual_herb"],
            "diet_to_island":res["diet_annual_herb_to_log_island_added"],
            "diet_to_global":res["diet_annual_herb_to_log_global_gain"],
            "all_panel_diet_island":res["all_butterflies_diet_to_island_cells"]
        },sort_keys=True)[:16000],flush=True)
    return out

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("hosts","crosswalk","panel","descriptors","edges","out","species-out"):
        p.add_argument("--"+key,type=Path,required=True)
    p.add_argument("--draws",type=int,default=799)
    a=p.parse_args()
    d,diagnosis,_=load(a.hosts,a.crosswalk,a.panel,a.descriptors,a.edges)
    res={
      "schema":"chocho_WCVP13_HOSTS_full_diet_GIFT32_annualherb_butterfly_island_resource_v1",
      "status":"FULL_LARVAL_DIET_NO_GAIN_SELECTION_EXPLORATORY_NOT_REAL_COLONIZATION",
      "source_gift":41946522894,"source_q2":41398149093,
      "source_host_database_commit":"808e0b869f9ec1adf8efff87cf6a395adda103e0",
      "source_wcvp_commit":"65bed76bae9d644ccb6ad200c05f9f5071d89e05",
      "original_sources_sha256":{k:hashlib.sha256(getattr(a,k).read_bytes()).hexdigest()
         for k in ("hosts","crosswalk","panel","descriptors","edges")},
      "source_population":diagnosis,
      "scenarios":model_summaries(d,B=a.draws),
      "limitations":[
        "Diet traits include all pinned, exactly resolved host plant taxa, independent of whether they introduced opportunity to a specific island.",
        "Host list from HOSTS is not proof of plant utilization in a specific local island or butterfly reproduction.",
        "GIFT recorded life cycle and herb growth not direct ruderal or nutrient-site preference.",
        "Genus bootstrap does not account for shared host-plant phylogenetic and geographical dependence across butterflies.",
        "Global potential resource expansion differs from island-only added resource count. Count lacks native-island denominator; adjusted for global native breadth.",
        "Host taxa lacking life-history data remain unknown and enter coverage diagnostics, not negative class.",
        "Country, port infrastructure and historical introductions are not explicitly fitted in this test.",
        "Source butterfly panel trait and host completeness nonrandom, with taxonomic sensitivity subsets."
      ]}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(res,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    a.species_out.parent.mkdir(parents=True,exist_ok=True)
    d.to_csv(a.species_out,index=False)
    print("FULL_WCVP_HOSTS_239_BUTTERFLY_SOURCE_DIAGNOSTIC",json.dumps(diagnosis),flush=True)
if __name__=="__main__":main()
