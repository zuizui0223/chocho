#!/usr/bin/env python3
"""Conservative structural null for full (not island-selected) butterfly host diets.

Source: WCVP13 accepted HOSTS butterfly×plant IDs across fixed 239 butterflies;
GIFT3.2 life-cycle/growth form, WCVP13 global plant native/alien range records,
Q2 v02 5,741 distinct *potential* new butterfly×island resource opportunities.

Fixed-source plant-trait label shuffle stratified by host family and shared
butterfly degree, optionally global native range. Assesses trait association
BEYOND those network/taxonomic factors, not direct causal trait influence.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import pyreadr

from audit_full_HOSTS_GIFT_diet_island_resource import load,id_key

SEED=20261105


def native_range_counts(path:Path,plant_ids:set[str],strict_islands:set[str]):
    d=pyreadr.read_r(str(path))["wcvp_distributions"]
    req={"plant_name_id","area_code_l3","introduced","extinct","location_doubtful"}
    if not req.issubset(d.columns):
        raise ValueError("Pinned WCVP regional distribution schema missing fields")
    d=d[list(req)].copy()
    d["plant_id"]=d.plant_name_id.map(id_key)
    d=d[d.plant_id.isin(plant_ids)].copy()
    for c in ("introduced","extinct","location_doubtful"):
        d[c]=pd.to_numeric(d[c],errors="coerce")
    d=d[(d.extinct==0)&(d.location_doubtful==0)&d.introduced.isin([0,1])]
    d=d[["plant_id","area_code_l3","introduced"]].drop_duplicates()
    if d.duplicated(["plant_id","area_code_l3"]).any():
        raise ValueError("Conflicting WCVP native introduced status at plant×region")
    d["island"]=d.area_code_l3.isin(strict_islands)
    out=d.groupby("plant_id").agg(
        global_native_regions=("introduced",lambda x:int((x==0).sum())),
        global_introduced_regions=("introduced",lambda x:int((x==1).sum())),
        island_introduced_regions=("introduced",lambda x:int(((x==1)&d.loc[x.index,"island"]).sum())),
        island_native_regions=("introduced",lambda x:int(((x==0)&d.loc[x.index,"island"]).sum()))
    ).reset_index()
    return out


def prepared(full_host,trait_cross,range_data,summary,*,scored_only=True,minimum_known=2,coverage=.30):
    host=pd.read_csv(full_host,dtype={"accepted_plant_name_id":str},keep_default_na=False)
    host["plant_id"]=host.accepted_plant_name_id.map(id_key)
    if host[["insect_species","plant_id"]].duplicated().any():
        raise ValueError("Repeated full host edges")
    c=pd.read_csv(trait_cross,dtype={"plant_name_id":str},keep_default_na=False)
    c["plant_id"]=c.plant_name_id.map(id_key)
    if c.plant_id.duplicated().any():raise ValueError("Duplicated GIFT accepted IDs")
    c["known"]=c.lifecycle.isin(["annual","biennial","perennial"])&c.growth_form.isin(["herb","shrub","tree"])
    c=c[c.known].copy()
    c["annual_herb"]=(c.lifecycle.eq("annual")&c.growth_form.eq("herb")).astype(float)
    h=host.merge(c[["plant_id","annual_herb"]],on="plant_id",how="inner",validate="many_to_one")
    plant=h[["plant_id","annual_herb","family"]].drop_duplicates("plant_id")
    r=pd.read_csv(range_data,dtype={"plant_id":str})
    plant=plant.merge(r,on="plant_id",how="left",validate="one_to_one")
    for col in ("global_native_regions","global_introduced_regions","island_introduced_regions"):
        plant[col]=plant[col].fillna(0).astype(float)
    plant["degree"]=plant.plant_id.map(h.groupby("plant_id").insect_species.nunique()).astype(int)
    plant["family"]=plant.family.fillna("UNKNOWN").astype(str)
    plant["degree_bin"]=np.where(plant.degree==1,"degree1",
                        np.where(plant.degree<=3,"degree2to3","degree4plus"))
    plant["native_bin"]=np.where(plant.global_native_regions<=3,"native0to3",
                        np.where(plant.global_native_regions<=20,"native4to20","native21plus"))

    d=summary.copy()
    d=d[(d.known_both>=minimum_known)&(d.coverage>=coverage)]
    if scored_only:d=d[d.disturbance_binary.notna()]
    d=d.copy()
    idx={p:i for i,p in enumerate(plant.plant_id)}
    sidx={s:i for i,s in enumerate(d.species)}
    M=np.zeros((len(d),len(plant)),dtype=np.uint8)
    for row in h[["insect_species","plant_id"]].itertuples(index=False):
        if row.insect_species in sidx:M[sidx[row.insect_species],idx[row.plant_id]]=1
    denom=M.sum(axis=1).astype(float)
    if (denom<minimum_known).any():
        raise ValueError("Full host diet identity does not match primary butterfly trait counts")
    ytrait=M.astype(float)@plant.annual_herb.to_numpy()/denom
    if not np.allclose(ytrait,d.annual_herb_fraction):
        raise ValueError("All original 239 butterfly host trait fractions disagree")
    return d,plant,M,denom


def controls(d,with_disturbance):
    terms=["log_diet_hosts","log_native","log_host_families","coverage"]
    if with_disturbance:terms+=["disturbance_binary"]
    Z=np.column_stack([np.ones(len(d)),*(
      d[x].to_numpy(dtype=float) for x in terms)])
    if np.linalg.matrix_rank(Z)!=Z.shape[1]:
        raise ValueError("Covariate rank failure in full HOSTS diet network null")
    q,_=np.linalg.qr(Z,mode="reduced")
    return q,terms


def pnull(d,p,M,denom,*,outcome,stratify,B=799,seed=SEED):
    if stratify not in ("family_degree","family_degree_native"):
        raise ValueError("Unexpected plant network shuffle")
    groupcols=["family","degree_bin"]+(
        ["native_bin"] if stratify=="family_degree_native" else [])
    groups=p.groupby(groupcols,observed=True).indices
    vals=p.annual_herb.to_numpy(dtype=float)
    movable=[np.asarray(g,dtype=int) for g in groups.values()
             if len(g)>=2 and 0<np.sum(vals[g])<len(g)]
    q,term=controls(d,"disturbance_binary" in d.columns and d.disturbance_binary.notna().all())
    y=d[outcome].to_numpy(dtype=float)
    ry=y-q@(q.T@y)
    Mf=M.astype(float)
    def slope(lab):
        x=Mf@lab/denom
        x=x-q@(q.T@x)
        den=float(x@x)
        return float(x@ry/den) if den>1e-10 else None
    observed=slope(vals)
    if observed is None:raise ValueError("No residual plant-trait variation")
    rng=np.random.default_rng(seed)
    numbers=[]
    for _ in range(B):
        label=vals.copy()
        for group in movable:
            label[group]=vals[rng.permutation(group)]
        b=slope(label)
        if b is not None:numbers.append(b)
    if len(numbers)!=B:
        raise ValueError("Network null is not identifiable")
    mu=float(np.mean(numbers))
    return {
      "response":outcome,"shuffle":stratify,
      "butterfly_species":len(d),"butterfly_genera":int(d.genus.nunique()),
      "original_beta":observed,"null_mean":mu,
      "null_q025_q975":list(map(float,np.quantile(numbers,[.025,.975]))),
      "p_two_centered_descriptive":(1+sum(abs(b-mu)>=abs(observed-mu)-1e-12 for b in numbers))/(B+1),
      "shuffle_species_count":int(len(set(i for g in movable for i in g))),
      "known_plant_species_count":len(p),"shuffle_strata":len(movable),
      "plant_network_preserved":True,
      "source_flora_traits_missing_preserved":True,
      "controls":term,
      "no_biological_introduction_or_breeding_inference":True}


def run(host,cross,wdist,panel,desc,edges,q1,*,B=799):
    d,diag,_=load(host,cross,panel,desc,edges)
    q1source=json.loads(Path(q1).read_text(encoding="utf-8"))
    island=set(q1source["classification"]["strict_island_codes"])
    if len(island)!=126:raise ValueError("Original strict island source drift")
    hostp=pd.read_csv(host,dtype={"accepted_plant_name_id":str})
    names=set(hostp.accepted_plant_name_id.map(id_key))
    ranges=native_range_counts(wdist,names,island)
    result={
      "schema":"chocho_full_HOSTS_GIFT_plant_family_degree_native_range_shuffle_v1",
      "status":"OBSERVATIONAL_ALL_DIET_NETWORK_NULL_NOT_CAUSAL_PLANT_TRAIT_MEDIATION",
      "source_precedent_R_full_rebuild":41974576917,
      "source_python_stream_rebuild":41975502719,
      "source_Q2":41398149093,
      "source_WCVP_GIFT":41946522894,
      "source_range_sha256":hashlib.sha256(Path(wdist).read_bytes()).hexdigest(),
      "n_distinct_wcvp_hosts":len(names),
      "n_plant_hosts_with_observed_range":len(ranges),
      "cohorts":{}
    }
    for label,scored,n,cover in (("all_191_source_positive_and_zero_gain",False,2,.30),
                                 ("disturbance_146_species",True,2,.30),
                                 ("disturbance_tight_116_species",True,3,.50)):
        subset,plants,M,denom=prepared(host,cross,ranges,d,
                    scored_only=scored,minimum_known=n,coverage=cover)
        rec={"n_species":len(subset),"n_genera":int(subset.genus.nunique()),
             "source_plant_N":len(plants),"source_known_host_diet_island_zero_count":int((subset.island_added==0).sum()),
             "nulls":{}}
        for response in ("log_island_added","log_global_gain"):
            for strat in ("family_degree","family_degree_native"):
                a=pnull(subset,plants,M,denom,outcome=response,stratify=strat,
                         B=B,seed=SEED+len(label)+len(response)+len(strat))
                rec["nulls"][response+"__"+strat]=a
        result["cohorts"][label]=rec
        print("FULL_DIET_GIFT_NETWORK_NULL",label,json.dumps({
           "n":len(subset),
           "nulls":{k:{"beta":v["original_beta"],"expected":v["null_mean"],
                       "p":v["p_two_centered_descriptive"],
                       "movable":v["shuffle_species_count"]}
                    for k,v in rec["nulls"].items()}},sort_keys=True),flush=True)
    # descriptive plant trait vs WCVP geographic introduced island range
    plants=p
    result["limitations"]=[
      "Family × butterfly degree permutations preserve broad phylogenetic proxy but not exact plant phylogeny or trait syndromes.",
      "Further matching on native regional breadth is a conservative constraint that can remove a real trait-to-range association; not a definitive no-effect test.",
      "All larval host plant taxa are from a nonrandom fixed 239 butterfly panel, not entire terrestrial plant diversity.",
      "Annual-herb is not actual nitrophily, field ruderal tolerance, direct larval performance or sewage nutrient delivery.",
      "Butterfly island resource count is geographically possible plant host occupancy, not observed butterfly establishment.",
      "No historical ship route, plant introduction year or local colonization sequence is identified.",
      "Country, botanical observational error and shared plant-host interactions can remain confounded even after this shuffle."
    ]
    return result


def main():
    pa=argparse.ArgumentParser()
    for key in ("host","cross","wcvp-dist","panel","desc","edges","q1","out"):
        pa.add_argument("--"+key,type=Path,required=True)
    pa.add_argument("--draws",type=int,default=799)
    a=pa.parse_args()
    r=run(a.host,a.cross,a.wcvp_dist,a.panel,a.desc,a.edges,a.q1,B=a.draws)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(r,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print("SOURCE_PINNED_FULL_DIET_GIFT_NETWORK_NULL_COMPLETE",r["status"],flush=True)
if __name__=="__main__":main()
