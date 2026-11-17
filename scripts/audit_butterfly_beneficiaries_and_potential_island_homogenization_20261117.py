#!/usr/bin/env python3
"""Which butterflies gain potential resources on islands, and does potential
consumer composition become more similar across islands?

Original 239 butterflies x 369 WGSRPD3 source resource cells, 126 island
regions, retained native and contemporary WCVP-HOSTS opportunity maps.
The "breadth" of native HOSTS-plant availability is NOT a butterfly's range.
Historical GBIF breadth on 51 selected landmasses is a distinct secondary
indexed-observation proxy, NOT butterfly island endemism.

Exploratory, source-locked, noncausal; source was previously inspected.
Fixed-margin swaps preserve every island added-consumer richness, every
butterfly number of gained islands, and exclude native-covered cells.
A second null additionally preserves each butterfly's gains by Level-1 region.
"""
import argparse, hashlib, json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT=Path(__file__).resolve().parents[1]
SOURCES={
 "cells":("results/island_resource_gap_filling/resource_cells.csv.gz","1ba2c5b5377248237235074b15d6950cc8655447"),
 "species":("results/island_resource_gap_filling/by_species.csv","a230f6d214c96f874a3ab184ee9823cc9023c1bb"),
 "metrics":("data/frozen/figure_sources/anthropogenic_species_metrics.csv","9f83d2d9ac8e46f1cfd4ded572c7a023273503dd"),
 "regions":("results/island-mainland-resource-q1/region_metrics.csv","4aba89565354edc6c44d566523864cec66a30cbe"),
 "gbif51":("data/exploratory/public51_butterfly_index_20261109/full51_225butterfly_period_index.csv","b6a3c98eb76864052d45d520a55b7d830816126c")
}
SEED=20261117
DRAWS=199
BURN_MULT=50
PER_DRAW_MULT=8

def gitblob(path):
    b=path.read_bytes()
    return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()

def pin_sources():
    details={}
    for name,(rel,pin) in SOURCES.items():
        path=ROOT/rel
        sha=gitblob(path)
        if sha != pin:
            raise ValueError(f"Source {name} modified: {sha}, expected {pin}")
        details[name]={"path":rel,"git_blob_sha":pin,"size_bytes":path.stat().st_size}
    return details

def rank_quartiles(species):
    sort=species.sort_values(["native_resource_units","species"]).copy()
    sort["native_host_footprint_quartile"]=np.minimum(4, np.floor(4*np.arange(len(sort))/len(sort)).astype(int)+1)
    return sort

def consume_jaccards(mask,rows,columns,valid=None):
    mask=mask.astype(np.int32)
    richness=mask.sum(axis=0).astype(np.float64)
    inter=(mask.T@mask).astype(np.float64)
    u=richness[:,None]+richness[None,:]-inter
    tri=np.triu_indices(mask.shape[1],1)
    union=u[tri]
    good=union>0
    if valid is not None:
        good &= valid
    out=np.full(len(union),np.nan,dtype=np.float64)
    out[good]=inter[tri][good]/union[good]
    return out,good,inter[tri],union

def swap_preserving_margins(original,native,level1,regional,rng,draws=DRAWS):
    A=original.copy().astype(np.bool_)
    E=np.argwhere(A).astype(np.int32)
    L=A.sum()
    if not np.array_equal(E.shape,[L,2]):raise ValueError("Wrong added edge list")
    nsp,nis=A.shape
    if len(level1)!=nis:raise ValueError("Wrong origin stratum roster")
    proposals=0;accepted=0
    def do_swaps(count):
        nonlocal proposals,accepted
        pairs=rng.integers(0,L,size=(count,2))
        for x,y in pairs:
            proposals+=1
            if x==y:continue
            sp1,isl1=E[x];sp2,isl2=E[y]
            if sp1==sp2 or isl1==isl2:continue
            if regional and level1[isl1]!=level1[isl2]:continue
            if native[sp1,isl2] or native[sp2,isl1]:continue
            if A[sp1,isl2] or A[sp2,isl1]:continue
            A[sp1,isl1]=False;A[sp2,isl2]=False
            A[sp1,isl2]=True;A[sp2,isl1]=True
            E[x,1]=isl2;E[y,1]=isl1
            accepted+=1
    # Burn in, then continue one fixed-margin Markov chain. Null draws are
    # dependent; intervals and tail counts are descriptive Monte Carlo only.
    do_swaps(BURN_MULT*L)
    null=[]
    for i in range(draws):
        do_swaps(PER_DRAW_MULT*L)
        if i in (0,draws-1):
            if not np.array_equal(A.sum(axis=1),original.sum(axis=1)) or not np.array_equal(A.sum(axis=0),original.sum(axis=0)):
                raise AssertionError("Row/column degree-preserving null failed")
            if np.any(A&native):raise AssertionError("Structural native-cell exclusion failed")
            if regional:
                for c in np.unique(level1):
                    sel=level1==c
                    if not np.array_equal(A[:,sel].sum(axis=1),original[:,sel].sum(axis=1)):
                        raise AssertionError("Species-by-bioregion margin changed")
        null.append(A.copy())
    return null,{"seed":SEED+(1 if regional else 0),
                 "regional_preserves_species_by_level1":bool(regional),
                 "draws":draws,"proposals":proposals,"accepted":accepted,
                 "acceptance_fraction":accepted/proposals,
                 "source_margin_preservation_checks":"PASS"}

def descriptive_null(obs,null):
    x=np.asarray(null,dtype=float)
    return {"observed":float(obs),"null_mean":float(x.mean()),
            "null_sd":float(x.std(ddof=1)),
            "null_2p5_97p5":[float(v) for v in np.quantile(x,[.025,.975])],
            "observed_minus_null_mean":float(obs-x.mean()),
            "upper_tail_descriptive_mc":float((1+(x>=obs).sum())/(len(x)+1)),
            "lower_tail_descriptive_mc":float((1+(x<=obs).sum())/(len(x)+1))}

def main():
    pa=argparse.ArgumentParser()
    pa.add_argument("--out-dir",type=Path,required=True)
    pa.add_argument("--n-draws",type=int,default=DRAWS)
    args=pa.parse_args()
    if args.n_draws!=DRAWS:raise ValueError("Frozen null 199 draws must not be changed after outcome inspection")
    pinned=pin_sources()
    species=pd.read_csv(ROOT/SOURCES["species"][0])
    metrics=pd.read_csv(ROOT/SOURCES["metrics"][0])
    regions=pd.read_csv(ROOT/SOURCES["regions"][0])
    raw=pd.read_csv(ROOT/SOURCES["cells"][0])
    gbif=pd.read_csv(ROOT/SOURCES["gbif51"][0])
    assert len(metrics)==239 and len(species)==478 and len(raw)==239*369 and len(gbif)==225*51
    assert set(species.group)=={"island","mainland"}
    assert (regions["class"].value_counts().to_dict()=={"mainland":237,"island":126,"mixed":6})
    assert regions.area_code_l3.is_unique and not raw.duplicated(["species","area_code_l3"]).any()
    native=raw.native_hosts.gt(0).astype(np.uint8)
    added=raw.contemporary.astype(int)-native
    assert added.isin((0,1)).all()
    assert int(native.sum())==26530 and int(added.sum())==14553

    island=regions[regions["class"]=="island"].sort_values("area_code_l3")
    island_codes=island.area_code_l3.tolist()
    sps=sorted(metrics.species)
    j=raw[raw.area_code_l3.isin(island_codes)].copy()
    assert len(j)==239*126
    j["native"]=j.native_hosts.gt(0).astype(np.uint8)
    j["added"]=j.contemporary.astype(int)-j["native"]
    N=(j.pivot(index="species",columns="area_code_l3",values="native")
         .reindex(index=sps,columns=island_codes)).to_numpy(dtype=bool)
    A=(j.pivot(index="species",columns="area_code_l3",values="added")
         .reindex(index=sps,columns=island_codes)).to_numpy(dtype=bool)
    C=N|A
    assert N.shape==(239,126) and int(N.sum())==4829 and int(A.sum())==5741 and int(C.sum())==10570
    assert not (A & N).any()
    bs=species[species.group=="island"][["species","native","added","fill","gain"]].rename(
            columns={"native":"native_island_regions",
                     "added":"added_island_regions",
                     "fill":"fraction_native_island_gaps_filled"})
    data=metrics[["species","native_resource_units","resolved_host_species",
                  "host_family_count","host_taxonomy_lower_bound_adequate"]].merge(
                  bs,on="species",validate="one_to_one")
    assert len(data)==239 and data.added_island_regions.sum()==5741
    data=rank_quartiles(data)
    qsum=[]
    for q in range(1,5):
        z=data[data.native_host_footprint_quartile==q]
        total=z.added_island_regions.sum()
        qsum.append({
           "quartile":q, "butterflies":len(z),
           "min_native_HOSTS_resource_regions":int(z.native_resource_units.min()),
           "max_native_HOSTS_resource_regions":int(z.native_resource_units.max()),
           "new_island_resource_cells":int(total),
           "fraction_all_island_resource_additions":float(total/5741),
           "mean_added_island_regions_per_butterfly":float(z.added_island_regions.mean()),
           "proportion_species_with_added_island_resources":float((z.added_island_regions>0).mean()),
           "fraction_all_original_native_host_gaps_filled":
               float(total/(126*len(z)-z.native_island_regions.sum())),
           "mean_known_HOSTS_plant_species":float(z.resolved_host_species.mean())})
    rho=spearmanr(data.native_resource_units,data.added_island_regions)
    rho_fill=spearmanr(data.native_resource_units,data.fraction_native_island_gaps_filled)
    gbif_early=(gbif.assign(has_early=gbif.pre2001_2010_indexed_records.gt(0))
                    .groupby("species").has_early.sum()).rename("early_sampled_island_occupancy_proxy")
    obs=data.join(gbif_early,on="species",how="inner")
    assert len(obs)==225
    rho_gbif=spearmanr(obs.early_sampled_island_occupancy_proxy,obs.added_island_regions)
    rho_gbif_baseline=spearmanr(obs.early_sampled_island_occupancy_proxy,obs.native_resource_units)

    geo=np.array([int(x) for x in island.level1])
    tri=np.triu_indices(len(island),1)
    is_same=geo[tri[0]]==geo[tri[1]]
    nj,valid,_,_=consume_jaccards(N,N,A)
    cj,validC,_,_=consume_jaccards(C,N,A,valid=valid)
    assert validC.sum()==valid.sum()
    pair_delta=cj[valid]-nj[valid]
    obs_pair={
      "island_pairs_all":int(len(valid)),
      "pairs_comparable_both_native_and_contemporary_nonempty":int(valid.sum()),
      "pairs_with_native_union_zero_and_excluded":int((~valid).sum()),
      "mean_potential_butterfly_Jaccard_native":float(np.nanmean(nj)),
      "mean_potential_butterfly_Jaccard_contemporary":float(np.nanmean(cj)),
      "absolute_mean_Jaccard_change":float(pair_delta.mean()),
      "mean_Jaccard_change_percentage_points":float(100*pair_delta.mean()),
      "fraction_pairs_with_increased_potential_Jaccard":float((pair_delta>0).mean()),
      "fraction_pairs_with_decreased_potential_Jaccard":float((pair_delta<0).mean()),
      "fraction_pairs_unchanged":float((pair_delta==0).mean())}
    strata={}
    for name,mask in [("within_Level1",is_same),("between_Level1",~is_same)]:
        selected=valid&mask
        strata[name]={"pairs":int(selected.sum()),
                      "native_mean_Jaccard":float(np.nanmean(nj[selected])),
                      "contemporary_mean_Jaccard":float(np.nanmean(cj[selected])),
                      "mean_Jaccard_change":float(np.nanmean(cj[selected]-nj[selected]))}
    observed_delta=obs_pair["absolute_mean_Jaccard_change"]
    nulls={}
    for regional in (False,True):
        label="species_island_addition_margins_and_within_Level1_species_margins" if regional else "species_island_addition_margins_global"
        rng=np.random.default_rng(SEED+(1 if regional else 0))
        configurations,swap_summary=swap_preserving_margins(A,N,geo,regional,rng,draws=DRAWS)
        metrics_null=[]
        for sample in configurations:
            new,_good,_,_=consume_jaccards(N|sample,N,sample,valid=valid)
            metrics_null.append(float(np.nanmean(new-nj)))
        nulls[label]={"comparison":"raw added-resource Jaccard gain against constrained assignment preserving column and row gain totals",
                      "null":descriptive_null(observed_delta,metrics_null),"swaps":swap_summary}
    # Species-level descriptive boostrapping of quartile allocation; the
    # sampled catalog is not a random sample of all butterflies globally.
    rng=np.random.default_rng(SEED+100)
    low=data[data.native_host_footprint_quartile==1].added_island_regions.to_numpy()
    high=data[data.native_host_footprint_quartile==4].added_island_regions.to_numpy()
    diff=np.empty(4999)
    for k in range(len(diff)):
        lo=low[rng.integers(0,len(low),len(low))].mean()
        hi=high[rng.integers(0,len(high),len(high))].mean()
        diff[k]=hi-lo

    args.out_dir.mkdir(parents=True,exist_ok=True)
    keep=["species","native_resource_units","native_island_regions","added_island_regions",
          "fraction_native_island_gaps_filled","resolved_host_species",
          "host_family_count","native_host_footprint_quartile"]
    data[keep].sort_values(["native_host_footprint_quartile","species"]).to_csv(
           args.out_dir/"beneficiary_species_native_resource_breadth.csv",index=False)
    pd.DataFrame(qsum).to_csv(args.out_dir/"beneficiary_four_quartiles.csv",index=False)
    # Region and pair summaries contain no private butterfly-HOSTS edges.
    pd.DataFrame([{"scope":name,**v} for name,v in strata.items()]).to_csv(
          args.out_dir/"potential_butterfly_island_Jaccard_by_geographic_scope.csv",index=False)
    report={
      "schema":"chocho_butterfly_resource_beneficiary_native_footprint_and_potential_homogenization_v1",
      "status":"ORIGINAL_239x369_SOURCE_VERIFIED_EXPLORATORY_NOT_BUTTERFLY_ENDISM_OR_ESTABLISHMENT",
      "cohort":{"butterflies":239,"original_botanical_regions":369,"strict_islands":126,
                "original_native_resource_cells":26530,"original_global_added_cells":14553,
                "original_island_added_cells":5741,"gbif_early_51island_subset_species":225},
      "source_git_blobs":pinned,
      "who_benefits":{
         "breadth_type":"number of WGSRPD3 regions with at least one native documented HOSTS host; NOT butterfly geographic range",
         "native_host_footprint_quartiles":qsum,
         "spearman_native_HOSTS_breadth_vs_island_added":float(rho.statistic),
         "spearman_native_HOSTS_breadth_vs_gap_fill_fraction":float(rho_fill.statistic),
         "upper_vs_lower_quartile_mean_added_islands_difference":float(high.mean()-low.mean()),
         "upper_vs_lower_95pct_species_resampling_descriptive":[float(x) for x in np.quantile(diff,[.025,.975])],
         "early_GBIF_51island_breadth_sensitivity":{
            "matched_species":len(obs),
            "earlier_indexed_in_at_least_one_of_51_islands":int((obs.early_sampled_island_occupancy_proxy>0).sum()),
            "spearman_indexed_2001_2010_islands_vs_added_resource_islands":float(rho_gbif.statistic),
            "spearman_indexed_2001_2010_islands_vs_native_HOSTS_breadth":float(rho_gbif_baseline.statistic),
            "warning":"GBIF observed presence in 51 selected port-linked landmasses is recording effort affected, NOT validated native geographic range"}}
       ,
      "potential_island_consumer_composition":{
          "presence_definition":"At least one accepted known larval host in each island region, regardless of butterfly occurrence",
          "observed":obs_pair,
          "within_vs_between_Level1":strata,
          "degree_constrained_nulls":nulls,
          "warning":"Potential community Jaccard may rise arithmetically as resource opportunities increase; two structural nulls preserve per-island added richness and per-butterfly added-island count; fixed-width CIs and null tails descriptive."},
      "limits":[
       "Butterfly species outside the frozen 239 HOSTS panel are not assessed; no verified island-endemic/narrow-range butterfly classification",
       "Native host-resource breadth (plant distribution) is NOT observed butterfly geographic distribution",
       "Early butterfly GBIF footprint restricted to 51 WPI-selected landmass groups and uneven opportunistic detection",
       "Potential consumer assemblages are not observed butterfly communities or population abundances",
       "Contemporary WCVP introduced host ranges do not establish plant arrival years, local host quality, reproduction, or biological invasions",
       "A constrained reshuffling null tests allocation structure, not randomized causal plant movement; repeated swap configurations are Markov dependent",
       "Different sample and geometry from 104-region actual botanical and 51-landmass indexed butterfly marine-circuit similarity"
      ]}
    (args.out_dir/"BENEFICIARY_AND_POTENTIAL_HOMOGENIZATION_20261117.json").write_text(
          json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("BENEFICIARY_POTENTIAL_HOMOGENIZATION_RESULT="+json.dumps({
       "sample":report["cohort"],
       "quartiles":qsum,
       "rho_native_resource_footprint":report["who_benefits"]["spearman_native_HOSTS_breadth_vs_island_added"],
       "rho_recorded_early_butterfly":report["who_benefits"]["early_GBIF_51island_breadth_sensitivity"]["spearman_indexed_2001_2010_islands_vs_added_resource_islands"],
       "potential_Jaccard":obs_pair,
       "potential_nulls":{k:v["null"] for k,v in nulls.items()}},sort_keys=True),flush=True)

if __name__=="__main__":
    main()
