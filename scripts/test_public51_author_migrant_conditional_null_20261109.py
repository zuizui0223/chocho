#!/usr/bin/env python3
"""Independent author-listed migratory butterfly subset: conditioned island-index null.

Tests source-record proximity (not true dispersal), preserving positive
first-index counts per butterfly species and island, then species x region.
The publication-linked roster is positive-only; unlisted status is UNKNOWN,
not a biological sedentary control.
"""
from __future__ import annotations
import argparse,importlib.util,json,random
from collections import Counter
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
PARENT=ROOT/"scripts/run_public51_butterfly_trait_access_20261109.py"
MIGRANT=ROOT/"scripts/run_public51_author_selected_migrant_access_20261109.py"
NAMES=("nearest_prior_conspecific_geo_access","residual_uniform_ocean_access",
       "centered_wing_x_geo_access","centered_wing_x_ocean_access")
SEED=20261109

def module(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    if spec is None or spec.loader is None:raise RuntimeError("Expected pinned original-source code missing")
    x=importlib.util.module_from_spec(spec);spec.loader.exec_module(x)
    return x

def swap(edges,current,elig,rng,regional=None):
    a,b=rng.randrange(len(edges)),rng.randrange(len(edges))
    if a==b:return False
    u,v=edges[a],edges[b]
    if u[0]==v[0] or u[1]==v[1]:return False
    if regional is not None and regional[u[1]]!=regional[v[1]]:return False
    x,y=(u[0],v[1]),(v[0],u[1])
    if x not in elig or y not in elig or x in current or y in current:
        return False
    edges[a],edges[b]=x,y
    current.difference_update((u,v));current.update((x,y))
    return True

def synthetic():
    edges=[(1,1),(2,2),(3,3)]
    es=set(edges);rng=random.Random(50)
    allowed={(i,j) for i in range(1,4) for j in range(1,4)}
    for _ in range(10000):swap(edges,es,allowed,rng)
    assert es==set(edges)
    assert Counter(i for i,j in es)=={1:1,2:1,3:1}
    assert Counter(j for i,j in es)=={1:1,2:1,3:1}

def sample_stat(edges,lookup,wing):
    vals=np.asarray([lookup[e] for e in edges],float)
    sz=np.asarray([wing[e[0]] for e in edges],float)
    centered=sz-sz.mean()
    return np.array([np.mean(vals[:,0]),np.mean(vals[:,1]),
                     np.mean(centered*vals[:,0]),np.mean(centered*vals[:,1])])

def switch_support(events,elig,region=None):
    paircount=0;involved=set()
    ss=set(events)
    for ii,x in enumerate(events):
        for y in events[ii+1:]:
            if x[0]==y[0] or x[1]==y[1]:continue
            if region is not None and region[x[1]]!=region[y[1]]:continue
            if (x[0],y[1]) in elig and (y[0],x[1]) in elig and (x[0],y[1]) not in ss and (y[0],x[1]) not in ss:
                paircount+=1;involved.add(x);involved.add(y)
    return {"eligible_2x2_switch_pairs_initial":paircount,
            "positive_cells_initially_in_swaps":len(involved),
            "positive_cells":len(events),
            "fraction_initially_swappable":len(involved)/len(events),
            "support_does_not_prove_MCMC_mix":True}

def test_null(events,elig,lookup,wing,*,regional,seed,draws,chains=4):
    support=switch_support(events,elig,regional)
    if support["eligible_2x2_switch_pairs_initial"]<12 or support["fraction_initially_swappable"]<.35:
        return {"status":"REGION_CONDITIONED_NULL_UNDERPOWERED","support":support,
                "no_biological_negative_inference":True}
    observed=sample_stat(events,lookup,wing)
    out=[];diag=[]
    species_counts=Counter(k for k,i in events)
    target_counts=Counter(i for k,i in events)
    reg_counts=Counter((k,regional[i]) for k,i in events) if regional is not None else None
    for ch in range(chains):
        rng=random.Random(seed+911*ch)
        arr=list(events);current=set(arr)
        accepted=0;proposed=0
        burn=4500
        while accepted<burn:
            accepted+=int(swap(arr,current,elig,rng,regional))
            proposed+=1
            if proposed>1500000:raise ValueError("Conditional design stalled; do not infer p-values")
        chain=[]
        for _ in range(draws):
            moves=0
            while moves<60:
                moves+=int(swap(arr,current,elig,rng,regional));proposed+=1
                if proposed>12_000_000:raise ValueError("Conditional simulation failed to mix")
            chain.append(sample_stat(arr,lookup,wing))
        if Counter(k for k,i in arr)!=species_counts or Counter(i for k,i in arr)!=target_counts:
            raise ValueError("Fixed species/island positive count invariant lost")
        if regional is not None and Counter((k,regional[i]) for k,i in arr)!=reg_counts:
            raise ValueError("Fixed species/region positive count invariant lost")
        aa=np.asarray(chain,float);out.extend(chain)
        diag.append({"chain":ch,"proposed":proposed,"accepted":burn+60*draws,
                     "acceptance_fraction":(burn+60*draws)/proposed,
                     "statistics_mean":dict(zip(NAMES,aa.mean(axis=0).tolist()))})
    data=np.asarray(out,float);tests={}
    for k,name in enumerate(NAMES):
        v=data[:,k];center=float(v.mean());obs=float(observed[k])
        tests[name]={"observed":obs,"null_mean":center,
                     "observed_minus_null":obs-center,
                     "null_95pct":list(map(float,np.quantile(v,[.025,.975]))),
                     "p_two_sided_monte_carlo":float((1+np.sum(abs(v-center)>=abs(obs-center)))/(len(v)+1)),
                     "p_one_sided_high":float((1+np.sum(v>=obs))/(len(v)+1)),
                     "num_draws":len(v)}
    return {"status":"EXPLORATORY_AUTHOR_SELECTED_MIGRANT_CONDITIONAL_GBIF_INDEX",
            "support":support,"chains":diag,"tests":tests,
            "not_arrival_colonisation_or_actual_migration":True}

def run(outdir):
    synthetic()
    base=module(PARENT,"existing_public51_source")
    migr=module(MIGRANT,"original_author_migrant_proxy")
    panel,receipt=base.load_source(
        ROOT/"data/exploratory/public51_butterfly_index_20261109")
    roster=migr.migration_crosswalk(panel)
    src,trait,data,pins,bridge,bmeta=base.data_bridge()
    islands={z["landmass_id"] for z in panel}
    regional=base.geographic_blocks(data,islands)
    report={"schema":"chocho_public51_author_selected_migrant_two_margin_20261109_v1",
            "source_sha256":receipt["merged_table_sha256"],
            "roster_git_sha1":migr.CROSSWALK_BLOB,
            "status":"EXPLORATORY_WITHIN_POSITIVELY_AUTHOR_REVIEW_SELECTED_SPECIES_NOT_MIGRANT_VS_NONMIGRANT",
            "migrant_proxy_source":"Author-released GBIF predicate derived from published migration review; original Table S1 inaccessible",
            "nonlisted_migration":"UNKNOWN_NOT_NONMIGRANT",
            "null_method":"2x2 species/island case count fixed, + species/region case counts for constrained null",
            "scientific_target":"GBIF first INDEX of butterflies with previous indexed conspecific source within the 51 landmass sample; not real dispersal",
            "periods":{}}
    for idx,period in enumerate(("2011_2017","2018_2022")):
        rr,cov=base.make_rows(panel,period,data,bridge,regional)
        rr=[z for z in rr if roster[z["speciesKey"]]==migr.KNOWN]
        if len(rr)<200:raise ValueError("Author-selected candidates source undercovered")
        keys=[(int(z["speciesKey"]),int(z["landmass_id"])) for z in rr]
        elig=set(keys)
        if len(elig)!=len(keys):raise ValueError("Duplicate selected positive source")
        events=[key for key,z in zip(keys,rr) if z["record_index_case"]]
        if len(events)<15:raise ValueError("Insufficient positive records in author-selected migration roster")
        geo=np.asarray([z["features"]["min_log_geodesic_km"] for z in rr],float)
        sea=np.asarray([z["features"]["min_log_uniform_effective_resistance"] for z in rr],float)
        X=np.column_stack((np.ones(len(rr)),geo,geo*geo))
        beta=np.linalg.lstsq(X,sea,rcond=None)[0]
        access=np.column_stack((-geo,-(sea-X@beta)))
        sd=access.std(axis=0)
        if min(sd)<1e-9:raise ValueError("No distance-independent original ocean residual variation")
        access=(access-access.mean(axis=0))/sd
        lookup={key:access[i] for i,key in enumerate(keys)}
        wing={int(z["speciesKey"]):float(z["wing_size_z"]) for z in rr}
        within={int(k):int(v) for k,v in regional.items()}
        first=test_null(events,elig,lookup,wing,regional=None,
                        seed=SEED+idx*100,draws=1000,chains=4)
        spatial=test_null(events,elig,lookup,wing,regional=within,
                          seed=SEED+idx*100+41,draws=750,chains=4)
        report["periods"][period]={
            "risk_cells":len(rr),"indexed_first_cases":len(events),
            "candidate_species":len({z["speciesKey"] for z in rr}),
            "source_geographic_islands":len({z["landmass_id"] for z in rr}),
            "ocean_from_geo_source_only_quadratic_beta":beta.tolist(),
            "two_margin":first,"within_original_geographic_region":spatial}
    report["limitations"]=[
      "Author-review linked GBIF selected taxa are indirect documented migrant candidates; not a direct life-history dispersal measurement.",
      "Species absent from selected list are UNKNOWN, not confirmed non-migrants.",
      "Non-index does not prove biological absence; first index does not prove immigration, breeding or establishment.",
      "Region, species and island positive case totals alone do not control geographic visibility and seasonal observation effort.",
      "Sea resistance is uniform Circuitscape, not verified ship travel or self-propelled butterfly ocean crossings.",
      "No historic species-specific larval-host presence or plant introduction dates are used.",
      "Same recorded species can occur at several sampled islands and does not represent independent biological crossing events.",
      "Conditional MCMC p-values are exploratory under dependent swap draws and posthoc comparisons."
    ]
    outdir.mkdir(parents=True,exist_ok=True)
    path=outdir/"author_selected_migrant_two_margin_index_null.json"
    path.write_text(json.dumps(report,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf8")
    print(json.dumps({"output":str(path),"observations":{
      p:{"risk":q["risk_cells"],"events":q["indexed_first_cases"],
         "geo_global_p":q["two_margin"].get("tests",{}).get("nearest_prior_conspecific_geo_access",{}).get("p_two_sided_monte_carlo"),
         "sea_global_p":q["two_margin"].get("tests",{}).get("residual_uniform_ocean_access",{}).get("p_two_sided_monte_carlo"),
         "spatial_status":q["within_original_geographic_region"]["status"],
         "spatial_tests":{k:v["p_two_sided_monte_carlo"] for k,v in q["within_original_geographic_region"].get("tests",{}).items()}
        } for p,q in report["periods"].items()}},indent=2,ensure_ascii=False))
    return report

if __name__=="__main__":
    a=argparse.ArgumentParser()
    a.add_argument("--out-dir",type=Path,required=True)
    z=a.parse_args()
    run(z.out_dir)
