#!/usr/bin/env python3
"""Two-margin and within-region null for butterfly traits and indexed island access.

Outcome = subsequent public GBIF FIRST INDEX on 51 preselected source landmasses.
Never infer biological first arrival, successful flight, shipment or breeding.
A species by island first-index case is one positive cell (not an independent
replication of all other species and islands). Source and 51-island cohort are
validated byte-for-byte by companion public51 source importer.
"""
from __future__ import annotations
import argparse, importlib.util, json, random
from collections import Counter
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
SOURCE_SCRIPT=ROOT/"scripts/run_public51_butterfly_trait_access_20261109.py"
SEED=20261109
PERIODS=("2011_2017","2018_2022")
STAT_NAMES=(
 "geo_access_mean", "sea_residual_access_mean",
 "wing_proxy_x_geo", "wing_proxy_x_sea_residual",
 "small_minus_large_geo", "small_minus_large_sea_residual",
 "disturbance_x_geo", "disturbance_x_sea_residual"
)

def module(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    if spec is None or spec.loader is None:raise RuntimeError("Source module missing")
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    return m

def switch(edges,event_set,eligible,rng,regional=None):
    # A reversible 2x2 row-column switch preserves BOTH butterfly species
    # and target island event counts; regional also holds species×region.
    a,b=rng.randrange(len(edges)),rng.randrange(len(edges))
    if a==b:return False
    u,v=edges[a],edges[b]
    if u[0]==v[0] or u[1]==v[1]:return False
    if regional and regional[u[1]]!=regional[v[1]]:return False
    x,y=(u[0],v[1]),(v[0],u[1])
    if x not in eligible or y not in eligible or x in event_set or y in event_set:
        return False
    edges[a],edges[b]=x,y
    event_set.difference_update((u,v))
    event_set.update((x,y))
    return True

def synthetic_test():
    es=[(1,1),(2,2),(3,3),(4,4)]
    allowed={(s,i) for s in range(1,5) for i in range(1,5)}
    initial=set(es);rng=random.Random(23)
    for _ in range(8000):switch(es,initial,allowed,rng)
    assert set(es)==initial
    assert Counter(x for x,_ in es)=={1:1,2:1,3:1,4:1}
    assert Counter(i for _,i in es)=={1:1,2:1,3:1,4:1}

def get_scores(events,lookup,traits):
    v=np.asarray([lookup[e] for e in events],float)
    wing=np.asarray([traits[e[0]][0] for e in events])
    grp=np.asarray([traits[e[0]][1] for e in events])
    dis=np.asarray([traits[e[0]][2] for e in events])
    if not (np.any(grp==0) and np.any(grp==1)):
        raise ValueError("Wing groups missing among source positive index records")
    if not (np.any(dis==1) and np.any(dis==-1)):
        raise ValueError("Disturbance groups missing among positive index records")
    return np.array([
      np.mean(v[:,0]),np.mean(v[:,1]),
      # Under fixed species event margins, the event-weighted trait means
      # are constant. Centering them here removes the shared access main effect
      # from each trait-interaction statistic. Without this correction, an
      # overwhelmingly disturbance-associated positive cohort can fake a
      # significant disturbance × sea effect whenever sea has a main effect.
      np.mean((wing-wing.mean())*v[:,0]),
      np.mean((wing-wing.mean())*v[:,1]),
      np.mean(v[grp==0,0])-np.mean(v[grp==1,0]),
      np.mean(v[grp==0,1])-np.mean(v[grp==1,1]),
      np.mean((dis-dis.mean())*v[:,0]),
      np.mean((dis-dis.mean())*v[:,1]),
    ],dtype=float)

def valid_switch_support(events,eligible,regional=None):
    sw=0;involved=set()
    for idx,a in enumerate(events):
        for b in events[idx+1:]:
            if a[0]==b[0] or a[1]==b[1]:continue
            if regional and regional[a[1]]!=regional[b[1]]:continue
            if (a[0],b[1]) in eligible and (b[0],a[1]) in eligible and (
                a[0],b[1]) not in events and (b[0],a[1]) not in events:
                sw+=1;involved.add(a);involved.add(b)
    return {"initial_eligible_2x2_switch_pairs":sw,
            "initial_positive_events_in_at_least_one_switch":len(involved),
            "positive_events":len(events),
            "fraction_initially_swappable":len(involved)/len(events),
            "initial_support_not_markov_chain_ergodicity_proof":True}

def conditional_null(events,eligible,lookup,traits,*,seed,regional=None,draws=1000,chains=4):
    support=valid_switch_support(events,eligible,regional)
    if support["initial_positive_events_in_at_least_one_switch"]<0.35*len(events) or (
        support["initial_eligible_2x2_switch_pairs"]<10):
        return {"status":"NOT_IDENTIFIABLE_WITHIN_FIXED_STRATUM",
                "support":support,"not_a_negative_biological_result":True}
    truth=get_scores(events,lookup,traits)
    results=[]
    chain_diag=[]
    ss=Counter(k for k,v in events);ii=Counter(v for k,v in events)
    groups0=Counter((k,regional[v]) for k,v in events) if regional else None
    for ch in range(chains):
        rng=random.Random(seed+501*ch)
        a=events.copy();state=set(a)
        attempt=accept=0
        burn=6500
        while accept<burn:
            accept+=int(switch(a,state,eligible,rng,regional))
            attempt+=1
            if attempt>1500000:raise ValueError("Insufficient valid 2x2 switches under this design")
        vals=[]
        for it in range(draws):
            got=0
            while got<60:
                got+=int(switch(a,state,eligible,rng,regional));attempt+=1
                if attempt>10000000:raise ValueError("Conditional null chain stalled")
            vals.append(get_scores(a,lookup,traits))
        if Counter(k for k,v in a)!=ss or Counter(v for k,v in a)!=ii:
            raise ValueError("Species or island marginal broke")
        if regional and Counter((k,regional[v]) for k,v in a)!=groups0:
            raise ValueError("Species by geographic cluster margins broke")
        arr=np.stack(vals)
        results.extend(vals)
        chain_diag.append({
            "chain":ch,"accepted_switches":burn+60*draws,
            "proposed_switches":attempt,
            "acceptance_fraction":(burn+60*draws)/attempt,
            "sampled_mean":dict(zip(STAT_NAMES,[float(v) for v in arr.mean(axis=0)]))})
    arr=np.stack(results)
    out={}
    for j,name in enumerate(STAT_NAMES):
        nv=arr[:,j];obs=truth[j];center=float(nv.mean())
        out[name]={
          "observed":float(obs),"null_mean":center,
          "observed_minus_null":float(obs-center),
          "null_sd":float(nv.std(ddof=1)),
          "null_95":[float(z) for z in np.quantile(nv,[.025,.975])],
          "p_two_sided_monte_carlo":float((1+(abs(nv-center)>=abs(obs-center)).sum())/(len(nv)+1)),
          "p_direction_higher":float((1+(nv>=obs).sum())/(len(nv)+1)),
          "draws":int(len(nv)),
        }
    return {"status":"EXPLORATORY_CONDITIONAL_RECORD_INDEX_NULL",
            "support":support,"chains":chain_diag,"tests":out,
            "MCMC_NOT_TRUE_ECOLOGICAL_RANDOM_ASSIGNMENT":True}

def prepare(source_data,period,base,bridge,groups):
    rows,cov=base.make_rows(source_data,period,base.data_bridge()[2],bridge,groups)
    # "Wing size" is the frozen S1 mean proxy. It does NOT measure
    # actual dispersal ability, self-propelled crossing or ship transport.
    n=len(rows);keys=[(int(r["speciesKey"]),int(r["landmass_id"])) for r in rows]
    if len(set(keys))!=n:raise ValueError("Duplicate GBIF 51 risk row")
    idx={key:k for k,key in enumerate(keys)}
    lat=np.array([r["features"]["min_log_geodesic_km"] for r in rows])
    sea=np.array([r["features"]["min_log_uniform_effective_resistance"] for r in rows])
    if not np.isfinite(lat).all() or not np.isfinite(sea).all():
        raise ValueError("Source geography nonfinite")
    z=np.column_stack((np.ones(n),lat,lat**2))
    coeff=np.linalg.lstsq(z,sea,rcond=None)[0]
    access=np.column_stack((-lat,-(sea-z@coeff)))
    sd=access.std(axis=0)
    if np.any(sd<1e-9):raise ValueError("No independent source geography/resistance variation")
    access=(access-access.mean(axis=0))/sd
    lookup={key:access[i] for i,key in enumerate(keys)}
    sizes=np.array([r["wing_size_z"] for r in rows])
    sizes=(sizes-sizes.mean())/sizes.std()
    median=base.data_bridge()[5]["wing_size_proxy_median"]
    traits={}
    for index,r in enumerate(rows):
        k=int(r["speciesKey"])
        wing_original=bridge[r["species"]]["wing_size_proxy"]
        group=0 if wing_original<=median else 1
        tup=(float(sizes[index]),group,int(r["disturbance"]))
        if k in traits and traits[k]!=tup:
            # Centered z is constant within species, so should not drift
            raise ValueError("Trait values change between islands")
        traits[k]=tup
    events=[key for key,r in zip(keys,rows) if r["record_index_case"]]
    eligible=set(keys)
    if len(events)!=cov["prior_source_and_trait_events"]:
        raise ValueError("Reconstructed case count mismatch")
    regional={int(i):int(g) for i,g in groups.items()}
    info={"rows":len(rows),"positive_cases":len(events),
          "source_species":len({r["species"] for r in rows}),
          "islands":len({r["landmass_id"] for r in rows}),
          "lat_min":float(lat.min()),"lat_max":float(lat.max()),
          "sea_residual_quadratic_fit_coefficients":list(map(float,coeff)),
          "source_prior_only":True}
    return events,eligible,lookup,traits,regional,info

def run(source_dir,outdir):
    synthetic_test()
    base=module(SOURCE_SCRIPT,"fixed_source51_trait")
    panel,receipt=base.load_source(source_dir)
    src,trait,data,pins,bridge,bmeta=base.data_bridge()
    groups=base.geographic_blocks(data,{r["landmass_id"] for r in panel})
    result={
      "schema":"chocho_public51_butterfly_trait_two_margin_geo_conditional_20261109_v2",
      "status":"EXPLORATORY_FIRST_GBIF_INDEX_NOT_BIOLOGICAL_COLONIZATION",
      "source_sha256":receipt["merged_table_sha256"],
      "traits_git_sha1":trait.TRAITS[1],
      "data_selection":"199 LepTraits matched taxon names among 225 original butterflies; prior indexed same-species sources on original 51 only",
      "design":"2x2 species × island first-index switches; spatial sensitivity also constrains species × KMeans-7 geographic cluster counts",
      "interaction_calculation":"center each trait by its event-weighted mean (fixed by species event totals) before multiplying by access, removing the main access effect",
      "no_truth_assumptions":["no historical larval host at endpoint","no actual butterfly migration","no documented boat transport","non-index not ecological absence"],
      "periods":{},
    }
    outdir.mkdir(parents=True,exist_ok=True)
    for i,period in enumerate(PERIODS):
        # Use same frozen panel / geography / trait data for both epochs.
        rows,meta=base.make_rows(panel,period,data,bridge,groups)
        keys=[(int(r["speciesKey"]),int(r["landmass_id"])) for r in rows]
        wing_by_species={int(r["speciesKey"]):bridge[r["species"]] for r in rows}
        geo=np.asarray([r["features"]["min_log_geodesic_km"] for r in rows])
        sea=np.asarray([r["features"]["min_log_uniform_effective_resistance"] for r in rows])
        X=np.column_stack((np.ones(len(rows)),geo,geo**2))
        beta=np.linalg.lstsq(X,sea,rcond=None)[0]
        access=np.column_stack((-geo,-(sea-X@beta)))
        access=(access-access.mean(axis=0))/access.std(axis=0)
        lookup={k:access[j] for j,k in enumerate(keys)}
        w=np.asarray([r["wing_size_z"] for r in rows])
        wmu,wsd=float(w.mean()),float(w.std())
        median=float(bmeta["wing_size_proxy_median"])
        traits={k:(float((t["wing_size_z"]-wmu)/wsd),
                   int(t["wing_size_proxy"]>median),
                   int(t["disturbance_category"]))
               for k,t in wing_by_species.items()}
        eligible=set(keys)
        events=[key for key,r in zip(keys,rows) if r["record_index_case"]]
        if len(eligible)!=len(rows) or len(events)!=meta["prior_source_and_trait_events"]:
            raise ValueError("Wrong original source event risk cohort")
        region={int(k):int(v) for k,v in groups.items()}
        primary=conditional_null(events,eligible,lookup,traits,seed=SEED+i*100)
        regional=conditional_null(events,eligible,lookup,traits,seed=SEED+i*100+77,
                                  regional=region,draws=750,chains=4)
        result["periods"][period]={
          "source_risk_cells":len(rows),
          "source_first_index_event_cells":len(events),
          "source_species":len(wing_by_species),
          "source_landmasses":len({r["landmass_id"] for r in rows}),
          "primary_two_margin":primary,
          "within_geographic_cluster_two_margin":regional,
          "sea_residual_on_quadratic_geodesic":list(map(float,beta)),
        }
    path=outdir/"public51_butterfly_conditional_trait_marine.json"
    path.write_text(json.dumps(result,indent=2,ensure_ascii=False,sort_keys=True)+"\n",encoding="utf8")
    print(json.dumps({
      "source_sha":result["source_sha256"],
      "periods":{p:{
         "rows":z["source_risk_cells"],"events":z["source_first_index_event_cells"],
         "global":{name:{"effect":v["observed_minus_null"],
                          "p_two":v["p_two_sided_monte_carlo"]}
                   for name,v in z["primary_two_margin"].get("tests",{}).items()},
         "regional_status":z["within_geographic_cluster_two_margin"]["status"],
         "regional_support":z["within_geographic_cluster_two_margin"]["support"]
      } for p,z in result["periods"].items()}},indent=2))
    return result

if __name__=="__main__":
    a=argparse.ArgumentParser()
    a.add_argument("--source-dir",type=Path,required=True)
    a.add_argument("--out-dir",type=Path,required=True)
    x=a.parse_args()
    run(x.source_dir,x.out_dir)
