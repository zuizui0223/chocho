#!/usr/bin/env python3
"""Plant-network-stratified negative control for introduced-host trait × butterfly gain.

Uses frozen WCVP-v13 × GIFT-v3.2 accepted species and Q2 v02
introduced-host × butterfly × island edges, no network calls.
The TWO outcome grains are distinct:
  (a) actual count of newly POSSIBLE island butterfly×region RESOURCE CELLS,
  (b) frozen GLOBAL potential resource log gain.
Neither is observed colonization or native-island-normalized gain.

Conditioning on a trait computed solely among hosts that already generate
new island opportunities can yield selection bias; shuffled plant labels
preserving each plant's network and regional frequency supply a structural
null. This is not a confirmatory randomized ecological experiment.
"""
from __future__ import annotations
import argparse, hashlib, json, math
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
from exploratory_island_host_trait_butterfly_bridge import (
    read_and_check,plant_level_species_summary,genus_bootstrap)

SOURCE_RUN=41973320821
SEED=20261105
def sh(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def build(inputs, cohort="primary"):
    exp,candidates,panel,desc=inputs
    species=plant_level_species_summary(exp,panel,desc)
    species["log1p_added_island_cells"]=np.log1p(species.island_added_cells)
    species["log_global_gain"]=species.log_resource_expansion
    scored=species[species.disturbance_strong_binary.notna()].copy()
    if cohort=="primary":
        scored=scored[(scored.distinct_known_both>=2)&(scored.distinct_both_trait_coverage>=.3)].copy()
    elif cohort=="adequate":
        scored=scored[(scored.distinct_known_both>=2)&(scored.distinct_both_trait_coverage>=.3)&
          (scored.host_taxonomy_lower_bound_adequate.astype(str).str.casefold()=="true")].copy()
    elif cohort=="tight":
        scored=scored[(scored.distinct_known_both>=3)&(scored.distinct_both_trait_coverage>=.5)].copy()
    else:raise ValueError("Unknown cohort")
    valid=candidates.loc[candidates.known_both,["plant_id","annual_herb"]].copy()
    valid["plant_id"]=valid.plant_id.astype(str)
    if valid.plant_id.duplicated().any():raise ValueError("Repeated plant traits after key normalization")
    degree=exp.groupby("plant_id").agg(
        n_butterflies=("species","nunique"),
        n_island_regions=("island_code","nunique"),
        n_edges=("species","size")).reset_index()
    valid=valid.merge(degree,on="plant_id",how="left",validate="one_to_one")
    if valid.n_edges.isna().any():raise ValueError("Trait plant absent from network")
    # Degree and island range define *fixed, outcome-independent* bins,
    # formed without looking at annual labels or the butterfly outcome.
    valid["degree_class"]=pd.cut(valid.n_butterflies,bins=[0,1,3,10000],
        labels=["one_butterfly","2-3_butterflies","4plus_butterflies"],include_lowest=True)
    valid["island_class"]=pd.qcut(valid.n_island_regions.rank(method="first"),
        3,labels=["low_island_range","mid_island_range","high_island_range"])
    valid["stratum"]=valid.degree_class.astype(str)+"|"+valid.island_class.astype(str)
    plants=valid.plant_id.astype(str).tolist()
    idx={s:i for i,s in enumerate(plants)}
    links=exp[["species","plant_id"]].drop_duplicates()
    matrix=np.zeros((len(scored),len(plants)),dtype=np.uint8)
    sidx={s:i for i,s in enumerate(scored.species)}
    for item in links.itertuples(index=False):
        if item.species in sidx and str(item.plant_id) in idx:
            matrix[sidx[item.species],idx[str(item.plant_id)]]=1
    denom=matrix.sum(axis=1).astype(float)
    if (denom<2).any():
        raise ValueError("Primary cohort has butterfly with fewer than 2 distinct trait-known introduced island hosts")
    actual=valid.annual_herb.astype(float).to_numpy()
    fractions=matrix.astype(float)@actual/denom
    if not np.allclose(fractions,scored.distinct_annual_herb_fraction.to_numpy(),atol=1e-9):
        raise ValueError("Plant-edge trait metric and original species trait fraction diverge")
    groups=[np.asarray(v,dtype=int) for v in valid.groupby("stratum",observed=True).indices.values()]
    movable=[g for g in groups if len(g)>=2 and 0<sum(actual[g])<len(g)]
    return scored,valid,matrix,denom,actual,movable

def design(scored, response, with_disturbance=True):
    control=["log_introduced_host_count","log_native","log_host_families",
             "distinct_both_trait_coverage"]
    if with_disturbance:control+=["disturbance_strong_binary"]
    Z=np.column_stack([np.ones(len(scored)),*(
         scored[c].astype(float).to_numpy() for c in control)])
    if np.linalg.matrix_rank(Z)<Z.shape[1]:
        raise ValueError("Rank-deficient island host trait structural-null controls")
    q,_=np.linalg.qr(Z,mode="reduced")
    y=scored[response].astype(float).to_numpy()
    yr=y-q@(q.T@y)
    def slope(x):
        xr=x-q@(q.T@x)
        den=float(xr@xr)
        return float(xr@yr/den) if den>1e-12 else None
    return slope,control

def structural_null(scored,valid,matrix,denom,actual,movable,*,response,B=1999,seed=SEED):
    f,controls=design(scored,response)
    observed_fraction=(matrix.astype(float)@actual)/denom
    obs=f(observed_fraction)
    if obs is None:raise ValueError("No residual trait variance after controls")
    rng=np.random.default_rng(seed)
    d=np.asarray(actual,dtype=float)
    perm=[]
    for _ in range(B):
        shuffled=d.copy()
        for g in movable:
            shuffled[g]=d[rng.permutation(g)]
        trait=(matrix.astype(float)@shuffled)/denom
        b=f(trait)
        if b is not None:perm.append(b)
    if len(perm)<.99*B:raise ValueError("Null permutations not estimable")
    mu=float(np.mean(perm))
    return {
       "response":response,
       "island_butterfly_species":len(scored),
       "independent_butterfly_genera":int(scored.genus.nunique()),
       "network_plant_taxa_known_both":len(valid),
       "strata_total":len(valid.groupby("stratum")),
       "network_plant_taxa_exchangeable":int(len(set(i for g in movable for i in g))),
       "original_beta_trait_fraction":obs,
       "permutation_mean":mu,
       "permutation_q025_q975":list(map(float,np.quantile(perm,[.025,.975]))),
       "observed_minus_null_mean":obs-mu,
       "two_sided_centered_p_descriptive":(1+sum(abs(x-mu)>=abs(obs-mu)-1e-12 for x in perm))/(len(perm)+1),
       "positive_tail_p_descriptive":(1+sum(x>=obs-1e-12 for x in perm))/(len(perm)+1),
       "permutations":len(perm),
       "held_fixed":"each introduced plant identity, butterfly links, island occurrence count, GIFT trait missingness, butterfly gains, and network controls",
       "shuffle":"annual-herb labels among known plants in degree × island-range strata",
       "controls":controls,
       "causal_inference":False,
       "null_is_not_prospective_and_does_not_destroy_other_phylogenetic_associations":True
    }

def run(inputs,B=1999):
    out={"schema":"chocho_island_Q2_introduced_host_annual_trait_network_structural_null_v1",
         "status":"SOURCE_PINNED_EXPLORATORY_NOT_PORT_SELECTION_NITROGEN_OR_BUTTERFLY_COLONIZATION",
         "source_Q2_GIFT_bridge_run":SOURCE_RUN,
         "cohorts":{}}
    for label in ("primary","adequate","tight"):
        scored,valid,matrix,denom,actual,movable=build(inputs,label)
        if len(scored)<40:raise ValueError("Too few host-scored butterfly species "+label)
        cohort={"n_species":len(scored),"n_genera":int(scored.genus.nunique()),
                "group_count":scored.disturbance_strong_binary.value_counts().to_dict(),
                "structural_null":{}}
        for response in ("log1p_added_island_cells","log_global_gain"):
            # Genus-bootstrap sensitivity and plant structural null are distinct.
            boot=genus_bootstrap(scored,response,"distinct_annual_herb_fraction",
                ["disturbance_strong_binary","log_introduced_host_count","log_native",
                 "log_host_families","distinct_both_trait_coverage"],B=B,seed=SEED+len(label))
            null=structural_null(scored,valid,matrix,denom,actual,movable,
                 response=response,B=B,seed=SEED+len(label)+len(response))
            cohort["structural_null"][response]={"genus_bootstrap":boot,"plant_trait_label_null":null}
        out["cohorts"][label]=cohort
        print("HOST_TRAIT_ISLAND_GAIN_STRUCTURAL_NULL",label,
            json.dumps({k:{"coef":v["genus_bootstrap"].get("beta"),
                           "CI":v["genus_bootstrap"].get("beta_boot95"),
                           "plant_null":v["plant_trait_label_null"]["two_sided_centered_p_descriptive"],
                           "plant_null_beta":v["plant_trait_label_null"]["permutation_mean"]}
                        for k,v in cohort["structural_null"].items()},sort_keys=True),flush=True)
    out["limits"]=[
      "The raw island outcome is count of new butterfly×island *resource* cells, not realized abundance or normalized percentage island gain.",
      "Introduced hosts are selected for having generated added island potential: collider/selection and source geography can create effects.",
      "Distinct plant annual-herb traits not direct nitrogen preference nor butterfly host performance.",
      "Network shuffle preserves plant degree and island-range terciles but not exact phylogenetic relatedness or within-family trait distribution.",
      "Main plant annual-herb port-selection stage was not positively supported; ecological bridge is not a multi-step causal mechanism.",
      "Different butterflies share host plants and are not independent; genus resampling cannot solve shared plant pseudo-replication."
    ]
    return out

def main():
    p=argparse.ArgumentParser()
    for key in ("edges","crosswalk","panel","descriptors","out"):
        p.add_argument("--"+key,type=Path,required=True)
    p.add_argument("--draws",type=int,default=799)
    a=p.parse_args()
    data=read_and_check(a.edges,a.crosswalk,a.panel,a.descriptors)
    out=run(data,a.draws)
    out["sources_SHA256"]={
       key:sh(getattr(a,key)) for key in ("edges","crosswalk","panel","descriptors")
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print("ISLAND_HOST_TRAIT_NETWORK_STRUCTURAL_AUDIT_COMPLETE",out["status"],flush=True)

if __name__=="__main__":main()
