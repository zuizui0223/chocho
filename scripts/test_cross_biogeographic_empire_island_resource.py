#!/usr/bin/env python3
"""Q2-F: Cross-biogeographic-region shared-colonial-history resource network test.

The scientific estimand is descriptive: among island pairs from DIFFERENT
WGSRPD Level1 regions, is a shared >=50 year published empire history related
to higher Jaccard similarity of newly enabled butterfly resource opportunities?
Unlike global-only shuffling, exposure is permuted WITHIN each Level1 area,
and origin/destination Level1 PAIR fixed effects are included in the outcome
regression. This preserves the number of nodes affiliated with each empire
within each Level1 in every null draw and does not estimate a causal effect.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import random
from collections import defaultdict
from pathlib import Path

import numpy as np
from pyproj import Geod
from shapely.geometry import shape

EMPIRES = ("Great_Britain", "Spain", "Portugal", "Netherlands")
SEED = 20261102
B = 4999

def rows(path):
    with path.open(newline="", encoding="utf-8") as h:
        yield from csv.DictReader(h)

def jaccard(a, b):
    union = a | b
    return len(a & b) / len(union) if union else 0.

def geo_coordinates(path):
    data = json.loads(path.read_text(encoding="utf-8"))
    out = {}
    for f in data["features"]:
        props = f.get("properties") or {}
        code = str(props.get("LEVEL3_COD") or "").strip()
        if code and f.get("geometry"):
            p = shape(f["geometry"]).representative_point()
            out[code] = (float(p.x), float(p.y))
    if len(out) != 369:
        raise RuntimeError(f"Wrong official Level3 geometry codes: {len(out)}")
    return out

def make_model(nodes, hosts, butterflies, coords, *, level_relation="different"):
    n = len(nodes)
    if n < 12:
        return {"status": "TOO_FEW_ISLANDS", "nodes": n}
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)
             if (nodes[i]["level1"] != nodes[j]["level1"]) == (level_relation == "different")]
    if len(pairs) < 20:
        return {"status":"TOO_FEW_ISLAND_PAIRS","nodes":n,"pairs":len(pairs)}
    u=np.array([i for i,j in pairs],dtype=int)
    v=np.array([j for i,j in pairs],dtype=int)
    empires=np.array([[int(node["years"][e]>=50) for e in EMPIRES] for node in nodes],dtype=np.uint8)
    levs=[x["level1"] for x in nodes]
    keys=["|".join(sorted((levs[i],levs[j]))) for i,j in pairs]
    levels=sorted(set(levs))
    keygroups=defaultdict(list)
    for k,key in enumerate(keys):
        keygroups[key].append(k)
    geo=Geod(ellps="WGS84")
    distances=[]
    for i,j in pairs:
        a=coords[nodes[i]["code"]]
        b=coords[nodes[j]["code"]]
        _,_,m=geo.inv(*a,*b)
        distances.append(m/1000.)
    distances=np.array(distances)
    areas=np.log(np.array([x["area_km2"] for x in nodes]))
    hostsets=[hosts[x["code"]] for x in nodes]
    butterflysets=[butterflies[x["code"]] for x in nodes]
    hostcnt=np.array([len(x) for x in hostsets])
    buttercnt=np.array([len(x) for x in butterflysets])
    y_host=np.array([jaccard(hostsets[i],hostsets[j]) for i,j in pairs])
    y_butter=np.array([jaccard(butterflysets[i],butterflysets[j]) for i,j in pairs])
    X=np.column_stack((
        np.log1p(distances),
        np.abs(areas[u]-areas[v]),
        np.abs(np.log1p(hostcnt[u])-np.log1p(hostcnt[v])),
        np.abs(np.log1p(buttercnt[u])-np.log1p(buttercnt[v])),
    ))
    # First remove Level1-pair fixed effects from X and each outcome.
    def demean_group(z):
        result=z.astype(float).copy()
        for indices in keygroups.values():
            result[indices] -= np.mean(result[indices], axis=0)
        return result
    Z=demean_group(X)
    rank=np.linalg.matrix_rank(Z)
    if rank != X.shape[1]:
        raise RuntimeError(f"Geographic nuisance controls not full rank: {rank}")
    Q,_=np.linalg.qr(Z,mode="reduced")
    def residualized(v):
        w=demean_group(v)
        return w-Q@(Q.T@w)
    yh=residualized(y_host)
    yb=residualized(y_butter)
    def shared(matrix):
        return np.any((matrix[u]&matrix[v])>0,axis=1).astype(float)
    def statistic(x, y):
        xr=residualized(x)
        denominator=float(xr@xr)
        if denominator<=1e-10:return None
        return float(xr@y/denominator)
    obs_x=shared(empires)
    exposed=obs_x==1
    if not (exposed.any() and (~exposed).any()):
        return {"status":"NO_SHARED_VERSUS_NON_SHARED_HISTORY_VARIATION","nodes":n,"pairs":len(pairs)}
    rkeys={k: {"pairs":len(indices),
               "same_empire_pairs":int(np.sum(exposed[indices])),
               "different_empire_pairs":len(indices)-int(np.sum(exposed[indices]))}
           for k,indices in sorted(keygroups.items())}
    comparison=[k for k,x in rkeys.items() if x["same_empire_pairs"]>0 and x["different_empire_pairs"]>0]
    def one_test(kind, obs, fixed_y, labelblocks):
        buckets=defaultdict(list)
        for i,b in enumerate(labelblocks):
            buckets[b].append(i)
        movable=[ids for ids in buckets.values() if len(ids)>=2 and
                 len({tuple(empires[i]) for i in ids})>1]
        eligible=set(i for ids in movable for i in ids)
        rng=random.Random(SEED)
        vals=[]
        identical=0
        invalid=0
        for _ in range(B):
            perm=empires.copy()
            for ids in movable:
                order=ids.copy()
                rng.shuffle(order)
                perm[ids]=empires[order]
            x=shared(perm)
            if np.array_equal(x,obs_x):identical+=1
            s=statistic(x,fixed_y)
            if s is None:invalid+=1
            else: vals.append(s)
        if obs is None or not vals:
            return {"status":"NO_IDENTIFIABLE_HISTORY_CONTRAST",
                    "effective_islands":len(eligible),
                    "exchangeable_blocks":len(movable),"invalid_permutations":invalid}
        null=np.array(vals)
        m=float(np.mean(null))
        return {
            "status":"DESCRIPTIVE_NODE_PERMUTATION",
            "adjusted_shared_history_Jaccard_difference":obs,
            "permutation_mean":m,
            "observed_minus_null_mean":obs-m,
            "null_q025":float(np.quantile(null,.025)),
            "null_q975":float(np.quantile(null,.975)),
            "p_high":(1+int(np.sum(null>=obs-1e-12)))/(len(null)+1),
            "p_two_sided_centered":(1+int(np.sum(np.abs(null-m)>=abs(obs-m)-1e-12)))/(len(null)+1),
            "effective_islands":len(eligible),
            "exchangeable_blocks":len(movable),
            "identical_exposure_pair_draws":identical,
            "invalid_permutations":invalid,
            "valid_permutations":len(null)
        }
    groupings={
      "global":["all"]*n,
      "within_level1":[x["level1"] for x in nodes],
      "within_level1_area_half":[x["level1"]+"|"+("small" if x["area_km2"]<=
           float(np.median([q["area_km2"] for q in nodes])) else "large") for x in nodes],
    }
    observed_host=statistic(obs_x,yh)
    observed_butter=statistic(obs_x,yb)
    # Fixed Y/noise controls and history cohorts: labels are swapped by NODES, not pairs.
    tests={}
    for name, blocks in groupings.items():
        tests[name]={
          "host_plant_network":one_test(name,observed_host,yh,blocks),
          "butterfly_resource_opportunity_network":one_test(name,observed_butter,yb,blocks),
        }
    # Pair outcomes share islands, so report leave-one-ISLAND-out sign changes
    # instead of treating 798 dyads as 798 independent replicates.
    def leave_one_island_out_y(outcome):
        values=[]
        for removed in range(n):
            keep=(u != removed)&(v != removed)
            if int(np.sum(keep))<20:
                continue
            z=X[keep].astype(float).copy()
            x=obs_x[keep].astype(float).copy()
            y=outcome[keep].astype(float).copy()
            groups=defaultdict(list)
            for ix,k in enumerate(np.array(keys)[keep]):
                groups[k].append(ix)
            for ind in groups.values():
                z[ind]-=np.mean(z[ind],axis=0)
                x[ind]-=np.mean(x[ind])
                y[ind]-=np.mean(y[ind])
            if np.linalg.matrix_rank(z)!=z.shape[1]:
                continue
            q,_=np.linalg.qr(z,mode="reduced")
            xx=x-q@(q.T@x)
            yy=y-q@(q.T@y)
            denom=float(xx@xx)
            if denom<=1e-10:
                continue
            values.append({"excluded_code":nodes[removed]["code"],
                           "adjusted_beta":float(xx@yy/denom)})
        return {
            "valid_leave_one_island_out":len(values),
            "n_positive":sum(v["adjusted_beta"]>0 for v in values),
            "minimum_coefficient":min((v["adjusted_beta"] for v in values),default=None),
            "maximum_coefficient":max((v["adjusted_beta"] for v in values),default=None),
            "range_is_not_confidence_interval":True,
            "extreme_influence_cases":(
                sorted(values,key=lambda x:x["adjusted_beta"])[:3]+
                sorted(values,key=lambda x:x["adjusted_beta"])[-3:]
            ) if values else [],
        }
    loocv={
        "host_plant_network": leave_one_island_out_y(y_host),
        "butterfly_resource_opportunity_network": leave_one_island_out_y(y_butter),
    }
    exposure_by_region={
      reg:{
        "islands":sum(x["level1"]==reg for x in nodes),
        "affiliated_islands_by_empire":{
          e:sum(x["level1"]==reg and x["years"][e]>=50 for x in nodes) for e in EMPIRES
        },
      } for reg in levels
    }
    raw_host=float(np.mean(y_host[exposed])-np.mean(y_host[~exposed]))
    raw_butter=float(np.mean(y_butter[exposed])-np.mean(y_butter[~exposed]))
    return {
      "status":"DESCRIPTIVE_CROSS_REGION_PAIR_TEST",
      "n_islands":n,"pairs_nonindependent":len(pairs),
      "within_level1_pairs_excluded":n*(n-1)//2-len(pairs) if level_relation=="different" else len(pairs),
      "pairs_same_empire":int(np.sum(exposed)),
      "pairs_no_shared_empire":int(np.sum(~exposed)),
      "level1_pairs":len(keygroups),
      "level1_pairs_with_both_exposure_classes":len(comparison),
      "support_region_pairs_with_both_classes":comparison,
      "region_pair_support":rkeys,
      "regions_by_empire_presence":exposure_by_region,
      "mean_host_plant_jaccard":float(np.mean(y_host)),
      "mean_new_butterfly_resource_jaccard":float(np.mean(y_butter)),
      "raw_same_minus_no_shared_host_jaccard":raw_host,
      "raw_same_minus_no_shared_butterfly_jaccard":raw_butter,
      "controlled_shared_history_host_coefficient":observed_host,
      "controlled_shared_history_butterfly_coefficient":observed_butter,
      "island_leave_one_out_influence":loocv,
      "controls":"Unordered Level1 pair fixed effects, log great-circle distance, abs log WGSRPD area difference, abs log host richness difference, abs log butterfly resource opportunity richness difference.",
      "label_permutations":tests,
    }

def main():
    p=argparse.ArgumentParser()
    for k in ("q2-edges","q2-audit","empire-islands","region-area","wgsrpd3-geojson","out"):
        p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();a.out.parent.mkdir(parents=True,exist_ok=True)
    audit=json.loads(a.q2_audit.read_text(encoding="utf-8"))
    edges=list(rows(a.q2_edges))
    if audit["scope"]["introduced_butterfly_plant_island_edges"]!=len(edges):
        raise RuntimeError("Q2 full resource-host edge source changed")
    plants=defaultdict(set)
    butterflies=defaultdict(set)
    for e in edges:
        plants[e["island_code"]].add(e["plant_id"])
        butterflies[e["island_code"]].add(e["species"])
    if len(plants)!=122 or sum(map(len,butterflies.values()))!=5674:
        raise RuntimeError("Q1/Q2 island resource census changed")
    hist={r["area_code_l3"]:r for r in rows(a.empire_islands)}
    area={r["code"]:r for r in rows(a.region_area)}
    if len(hist)!=125:raise RuntimeError("Lenzner 125-island historical map changed")
    coords=geo_coordinates(a.wgsrpd3_geojson)
    nodes=[]
    for code in sorted(butterflies):
        h=hist[code]
        if code not in area or code not in coords:raise RuntimeError("Unmatched island "+code)
        nodes.append({"code":code,"level1":h["level1"],"area_km2":float(area[code]["area_km2"]),
                      "coverage":float(h["lenzner_polygon_coverage"]),
                      "strict":int(h["strict_historical_empire"])==1,
                      "years":{e:float(h[e+"_weighted_exposure_years"]) for e in EMPIRES}})
    def sample(coverage=.8,strict=False,no_europe=False):
        return [x for x in nodes if x["coverage"]>=coverage and
                any(y>=50 for y in x["years"].values()) and
                (not strict or x["strict"]) and
                (not no_europe or x["level1"]!="1")]
    populations={
        "primary_80pct_cover_50yrs":sample(),
        "sensitivity_50pct_cover_50yrs":sample(.5),
        "sensitivity_strict_single_empire":sample(.8,True),
        "sensitivity_non_european":sample(.8,False,True),
    }
    result={
        "schema":"chocho_cross_biogeographic_empire_butterfly_resource_v0.1",
        "status":"CROSS_REGION_HISTORICAL_RESOURCE_ASSOCIATION_NOT_CAUSATION",
        "scope":{"original_butterflies":239,"original_island_added_opportunities":5674,
                 "native_host_taxa_source":audit["scope"]["distinct_intro_host_taxa"],
                 "history_empires":list(EMPIRES)},
        "question":"Are shared historical empires associated with introduced-host butterfly-resource similarity even across different Level1 biogeographic regions?",
        "sample_definition":"Same source cohort as Q2-E, preserve all eligible historical empire affiliations >=50 years; only pairs between DIFFERENT WGSRPD Level1",
        "cohorts":{key:make_model(cohort,plants,butterflies,coords) for key,cohort in populations.items()},
        "limits":[
            "WGSRPD3 Level1 is a broad botanical region, not a randomized or naturally exogenous historical assignment.",
            "Permuting island history within Level1 preserves group-level empire prevalence, but does not eliminate within-region colonial-selection confounding.",
            "Controlled association and permutation tails from non-independent island dyads are descriptive, not proof of colonial causal influence.",
            "Only four European colonial empires and high-coverage islands are included; France and other histories cannot be inferred from this source.",
            "Same-empire resources can come from broadly native plant hosts present across several world regions; this alone does not demonstrate ship-mediated introductions.",
            "Jaccard butterfly opportunity is inferred from WCVP contemporary hosts, not actual island butterfly occupancy.",
            "No observed-larval-feeding/port-specific plant/dated shipping route evidence is supplied.",
        ],
    }
    a.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({"status":result["status"],"cohorts":{
      k:{"n_islands":v.get("n_islands"),
         "pairs":v.get("pairs_nonindependent"),"shared":v.get("pairs_same_empire"),
         "matched_region_pair_support":v.get("level1_pairs_with_both_exposure_classes"),
         "butterfly":v.get("label_permutations",{}).get("within_level1",{}).get("butterfly_resource_opportunity_network"),
         "host":v.get("label_permutations",{}).get("within_level1",{}).get("host_plant_network")}
      for k,v in result["cohorts"].items()}},indent=2,ensure_ascii=False))

if __name__ == "__main__":
    main()
