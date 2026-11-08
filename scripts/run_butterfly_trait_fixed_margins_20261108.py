#!/usr/bin/env python3
"""Species x island fixed-margin null: butterfly wing traits and prior ocean access.

Observations are indexed GBIF first records, not colonisation, breeding, nor
directly measured dispersal. Inputs are the original *public* 225 x 17 island
source and frozen LepTraits S1, each protected by exact original git-blob SHA.
Two-by-two switches preserve each species and each island first-index count.
Retrospective posthoc exploratory analysis; null draws are autocorrelated.
"""
from __future__ import annotations
import argparse
from collections import Counter
import importlib.util
import json
import random
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
SOURCE_SCRIPT=ROOT/"scripts/run_ncomms_public17_ocean_access_holdout.py"
TRAIT_SCRIPT=ROOT/"scripts/run_butterfly_trait_marine_access_20261108.py"
PERIODS=("2011_2017","2018_2022")
NAMES=("overall_geo","overall_sea_resid","small_geo","small_sea_resid",
       "large_geo","large_sea_resid","difference_small_minus_large_geo",
       "difference_small_minus_large_sea_resid",
       "continuous_wing_x_geodesic_access","continuous_wing_x_sea_residual_access",
       "disturbance_assoc_minus_avoid_geo","disturbance_assoc_minus_avoid_sea_resid")
SEED=20261108

def import_script(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    if spec is None or spec.loader is None:raise RuntimeError("Missing source module")
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module

def swap(edges,event_set,eligible,rng):
    i,j=rng.randrange(len(edges)),rng.randrange(len(edges))
    if i==j:return False
    a,b=edges[i],edges[j]
    if a[0]==b[0] or a[1]==b[1]:return False
    c,d=(a[0],b[1]),(b[0],a[1])
    if c not in eligible or d not in eligible or c in event_set or d in event_set:
        return False
    edges[i],edges[j]=c,d
    event_set.difference_update((a,b))
    event_set.update((c,d))
    return True

def stats(events,lookup,small,wing,disturbance):
    a=np.asarray([lookup[e] for e in events],dtype=float)
    group=np.asarray([e[0] in small for e in events])
    w=np.asarray([wing[e[0]] for e in events])
    dis=np.asarray([disturbance[e[0]]==1 for e in events])
    if not (group.any() and (~group).any() and dis.any() and (~dis).any()):
        raise ValueError("Missing positive event in a group")
    return np.asarray([
        a[:,0].mean(),a[:,1].mean(),
        a[group,0].mean(),a[group,1].mean(),
        a[~group,0].mean(),a[~group,1].mean(),
        a[group,0].mean()-a[~group,0].mean(),
        a[group,1].mean()-a[~group,1].mean(),
        np.mean(a[:,0]*w),np.mean(a[:,1]*w),
        a[dis,0].mean()-a[~dis,0].mean(),
        a[dis,1].mean()-a[~dis,1].mean()])

def test_synthetic_swap():
    edges=[(1,1),(2,2),(3,3)]
    eligible={(s,i) for s in range(1,4) for i in range(1,4)}
    s=set(edges);rng=random.Random(12)
    for _ in range(15000):swap(edges,s,eligible,rng)
    assert s==set(edges)
    assert Counter(x[0] for x in edges)=={1:1,2:1,3:1}
    assert Counter(x[1] for x in edges)=={1:1,2:1,3:1}

def summarize(obs,draws):
    out={}
    for k,name in enumerate(NAMES):
        v=draws[:,k];mu=float(v.mean());o=float(obs[k])
        out[name]={
          "observed":o, "null_mean":mu,
          "observed_minus_null":o-mu,
          "null_sd":float(v.std(ddof=1)),
          "null_95":[float(a) for a in np.quantile(v,[0.025,0.975])],
          "p_one_sided_high":float((1+(v>=o).sum())/(len(v)+1)),
          "p_one_sided_low":float((1+(v<=o).sum())/(len(v)+1)),
          "p_two_sided_centered":float((1+(np.abs(v-mu)>=abs(o-mu)).sum())/(len(v)+1)),
        }
    return out

def run_period(m,data,bridge,period):
    original,_=m.prepare(data,period)
    rows=[r for r in original if int(r["prior_source_islands"])>0
          and bridge[r["species"]].get("trait_available")]
    expected=(719,52) if period=="2011_2017" else (868,14)
    events=[(int(r["speciesKey"]),int(r["landmass_id"])) for r in rows
            if int(r["record_index_case"])==1]
    if (len(rows),len(events))!=expected:
        raise ValueError("Frozen original index-year sample changed")
    eligible={(int(r["speciesKey"]),int(r["landmass_id"])) for r in rows}
    assert len(eligible)==len(rows) and len(set(events))==len(events)
    near=np.asarray([r["features"]["min_log_geodesic_km"] for r in rows],float)
    sea=np.asarray([r["features"]["min_log_uniform_effective_resistance"] for r in rows],float)
    assert np.isfinite(near).all() and np.isfinite(sea).all()
    fit=np.linalg.lstsq(np.column_stack((np.ones(len(rows)),near,near*near)),sea,rcond=None)[0]
    residual=sea-(fit[0]+fit[1]*near+fit[2]*near*near)
    val=np.column_stack((-near,-residual))
    val=(val-val.mean(axis=0))/val.std(axis=0)
    lookup={(int(r["speciesKey"]),int(r["landmass_id"])):val[i] for i,r in enumerate(rows)}
    traits={int(r["speciesKey"]):bridge[r["species"]] for r in rows}
    small={k for k,z in traits.items() if z["wing_median_group"]=="small_50pct"}
    large={k for k,z in traits.items() if z["wing_median_group"]=="large_50pct"}
    assert small.isdisjoint(large)
    wing={k:float(z["wing_size_z"]) for k,z in traits.items()}
    disturbance={k:int(z["disturbance_category"]) for k,z in traits.items()}
    obs=stats(events,lookup,small,wing,disturbance)
    # Additional geographic-component blocking is NOT a valid full null when
    # very few event edges are swappable. Diagnose support before testing.
    component={int(r["landmass_id"]):int(r["component"]) for r in rows}
    can_swap=[]
    for ii,a in enumerate(events):
        for b in events[ii+1:]:
            if a[0]==b[0] or a[1]==b[1] or component[a[1]]!=component[b[1]]:
                continue
            cross1,cross2=(a[0],b[1]),(b[0],a[1])
            if cross1 in eligible and cross2 in eligible and cross1 not in events and cross2 not in events:
                can_swap.append((a,b))
    involved=set(x for pair in can_swap for x in pair)
    regional_support={
       "distinct_original_components_with_events":len({component[e[1]] for e in events}),
       "initially_valid_within_component_2x2_event_pairs":len(can_swap),
       "events_initially_involved_in_eligible_component_swaps":len(involved),
       "events_all":len(events),
       "fraction_initially_swappable":len(involved)/len(events),
       "regional_matched_null_not_promoted_without_support":True,
       "note":"Initial state structural support only, not a uniform randomization inference on regional groups"
    }
    if (period=="2011_2017" and (len(can_swap),len(involved))!=(6,8)) or (
        period=="2018_2022" and (len(can_swap),len(involved))!=(5,6)):
        raise ValueError("Original within-component swap support changed")
    samples=[]
    diagnostics=[]
    for ch in range(4):
        rng=random.Random(SEED+ch*1000+(0 if period=="2011_2017" else 600))
        arr=events.copy();es=set(arr)
        attempted=accepted=0
        while accepted<6000:
            accepted+=int(swap(arr,es,eligible,rng))
            attempted+=1
            if attempted>800000:raise RuntimeError("No MCMC mixing during burn-in")
        chain=[]
        for _ in range(1500):
            successes=0
            while successes<60:
                successes+=int(swap(arr,es,eligible,rng));attempted+=1
                if attempted>10_000_000:raise RuntimeError("2x2 null chain stalled")
            chain.append(stats(arr,lookup,small,wing,disturbance))
        assert Counter(x[0] for x in arr)==Counter(x[0] for x in events)
        assert Counter(x[1] for x in arr)==Counter(x[1] for x in events)
        samples.extend(chain)
        diagnostics.append({"chain":ch,"accepted_swaps":96000,
                            "proposed_swaps":attempted,
                            "acceptance_rate":96000/attempted,
                            "means":{key:float(np.mean([r[j] for r in chain]))
                                     for j,key in enumerate(NAMES)}})
    draws=np.asarray(samples)
    # Assert the six decisive contrasts are structurally present;
    # NOT asserting positive ecological findings in source code.
    tests=summarize(obs,draws)
    return {
      "eligible_species_x_island_cells":len(rows),
      "positive_record_first_indices":len(events),
      "small_positive_count":sum(e[0] in small for e in events),
      "large_positive_count":sum(e[0] in large for e in events),
      "disturbance_assoc_positive_count":sum(disturbance[e[0]]==1 for e in events),
      "disturbance_avoid_positive_count":sum(disturbance[e[0]]==-1 for e in events),
      "sea_log_resistance_on_geo_log_dist_quadratic":list(map(float,fit)),
      "regional_component_null_support_audit":regional_support,
      "four_chain_diagnostics":diagnostics,
      "tests":tests,
    }

def run(output:Path):
    test_synthetic_swap()
    source=import_script(SOURCE_SCRIPT,"original_public17_module")
    trait=import_script(TRAIT_SCRIPT,"original_trait_module")
    data,pins=source.load(ROOT)
    bridge,coverage=trait.build_bridge(data,trait.read_pinned_csv(*trait.TRAITS),
                                       trait.read_pinned_csv(*trait.RESOURCES))
    result={
       "schema":"chocho_trait_x_marine_fixed_species_and_island_first_index_margins_v1",
       "status":"EXPLORATORY_POST_HOC_RECORD_INDEX_NULL_NOT_DISPERSAL",
       "source_pins":pins,
       "trait_coverage":coverage,
       "methods":{
           "fixed_margins":"2x2 event/non-event cell swaps preserving BOTH species and target-island case counts",
           "eligible_at_risk":"strict historically not-yet-indexed cells with a prior conspecific other-island index; no biological absence",
           "access_metrics":"z(-log1p nearest source distance); z(-[log1p sea resistance - quadratic prediction from log1p distance])",
           "wing_test":"continuous wing x access, and 50/50 wing group differences",
           "chain":"4 independently seeded 1500-draw chains, 6000 accepted burn-in, 60 accepted 2x2 swaps between draws",
           "p_values":"Monte Carlo conditional two-sided deviation from sampled mean; exploratory multiple tests",
           "caution":"No verified movement, reproduction, host resources or sea shipping; 17 selected island groups; no site-level historical effort",
       },
       "periods":{p:run_period(source,data,bridge,p) for p in PERIODS}
    }
    output.mkdir(parents=True,exist_ok=True)
    p=output/"BUTTERFLY_TRAIT_SOURCE_ISLAND_FIXED_NULL.json"
    p.write_text(json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps({
      "output":str(p),
      "result":{per:{"n":r["eligible_species_x_island_cells"],
                     "events":r["positive_record_first_indices"],
                     "sea_resid_p":r["tests"]["overall_sea_resid"]["p_two_sided_centered"],
                     "small_large_sea_diff_p":r["tests"]["difference_small_minus_large_sea_resid"]["p_two_sided_centered"],
                     "continuous_wing_geo_p":r["tests"]["continuous_wing_x_geodesic_access"]["p_two_sided_centered"],
                     "continuous_wing_sea_p":r["tests"]["continuous_wing_x_sea_residual_access"]["p_two_sided_centered"]}
                for per,r in result["periods"].items()}},indent=2,ensure_ascii=False))
    return result

if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir",type=Path)
    ap.add_argument("--self-test",action="store_true")
    args=ap.parse_args()
    if args.self_test:
        test_synthetic_swap();print("PASS fixed species and island margins for 15,000 synthetic attempts")
    else:
        if args.out_dir is None:ap.error("--out-dir is required for actual run")
        run(args.out_dir)
