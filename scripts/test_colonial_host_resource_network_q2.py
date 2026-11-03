#!/usr/bin/env python3
"""Q2-E: Shared former-empire history and *potential* butterfly resource networks.

Population: 239 fixed butterfly species, WCVP-introduced larval-host opportunities
on 122 strict island WGSRPD3 regions. No GBIF butterfly observation is needed.

Historical comparison is restricted to independently geolocated islands with
>=50 weighted colonial years and published map coverage; full empire vectors
are permuted *as units* among island nodes within region/area strata. All pair
statistics are descriptive, not independent samples or evidence of ship routes.
"""
from __future__ import annotations
import argparse,csv,json,math,random
from collections import defaultdict,Counter
from pathlib import Path
import numpy as np
from pyproj import Geod
from shapely.geometry import shape


EMPIRES=("Great_Britain","Spain","Portugal","Netherlands")
SEED=20261102


def readcsv(path):
    with path.open(newline="",encoding="utf-8") as h:yield from csv.DictReader(h)


def rankquantile(seq,q):
    xs=sorted(seq)
    return xs[int((len(xs)-1)*q)]


def jaccard(a,b):
    size=len(a|b)
    return 0.0 if size==0 else len(a&b)/size


def region_representative_lonlat(geojson_path):
    features=json.loads(geojson_path.read_text(encoding="utf-8"))["features"]
    coords={}
    for f in features:
        code=str((f.get("properties") or {}).get("LEVEL3_COD") or "")
        if not code or not f.get("geometry"):continue
        poly=shape(f["geometry"])
        if poly.is_empty:continue
        rp=poly.representative_point()
        coords[code]=(float(rp.x),float(rp.y))
    return coords


def pair_model(cohort,plant_sets,butterfly_sets,coords,geod):
    n=len(cohort)
    if n<12:return {"status":"TOO_FEW_NODES","nodes":n}
    level=[r["level1"] for r in cohort]
    areas=np.array([math.log(r["area_km2"]) for r in cohort])
    yrs=np.array([[int(r["history"][e]>=50) for e in EMPIRES] for r in cohort],dtype=np.uint8)
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    u=np.array([i for i,j in pairs]);v=np.array([j for i,j in pairs])
    if len(pairs)<20:raise RuntimeError("Insufficient distinct island pairs")
    distances=[]
    for i,j in pairs:
        ax,ay=coords[cohort[i]["code"]]
        bx,by=coords[cohort[j]["code"]]
        _,_,distance_m=geod.inv(ax,ay,bx,by)
        distances.append(distance_m/1000)
    distance=np.array(distances)
    levelmatch=np.array([int(level[i]==level[j]) for i,j in pairs])
    ageom=np.abs(areas[u]-areas[v])
    host_sets=[plant_sets[r["code"]] for r in cohort]
    butterflies=[butterfly_sets[r["code"]] for r in cohort]
    host_counts=np.array([len(x) for x in host_sets])
    butter_counts=np.array([len(x) for x in butterflies])
    yplants=np.array([jaccard(host_sets[i],host_sets[j]) for i,j in pairs])
    ybutterflies=np.array([jaccard(butterflies[i],butterflies[j]) for i,j in pairs])
    # All nuisance variables determined independently of empire assignment.
    X=np.column_stack([
       np.ones(len(pairs)),
       np.log1p(distance),
       ageom,
       levelmatch,
       np.abs(np.log1p(host_counts[u])-np.log1p(host_counts[v])),
       np.abs(np.log1p(butter_counts[u])-np.log1p(butter_counts[v]))
    ])
    if np.linalg.matrix_rank(X)!=X.shape[1]:
        raise RuntimeError("Non-full-rank geographic/resource-richness control matrix")
    Q,_=np.linalg.qr(X,mode="reduced")
    resid_plant=yplants-Q@(Q.T@yplants)
    resid_butter=ybutterflies-Q@(Q.T@ybutterflies)
    def getx(yearmatrix):
        return np.any((yearmatrix[u]&yearmatrix[v])>0,axis=1).astype(float)
    def slope(yres,x):
        res=x-Q@(Q.T@x)
        denom=float(res@res)
        return None if denom<1e-12 else float((res@yres)/denom)

    original_x=getx(yrs)
    groups={
       "global":["all"]*n,
       "level1":level,
       "level1_area_half":[level[i]+"|"+("small" if areas[i]<=float(np.median(areas)) else "large") for i in range(n)]
    }
    base={
       "n_islands":n,"pairs_nonindependent":len(pairs),
       "same_empire_pairs":int(np.sum(original_x)),
       "other_empire_pairs":int(len(pairs)-np.sum(original_x)),
       "within_level1_pairs":int(sum(levelmatch)),
       "source_host_plant_jaccard_mean":float(np.mean(yplants)),
       "new_butterfly_opportunity_jaccard_mean":float(np.mean(ybutterflies)),
       "mean_distance_km":float(np.mean(distance)),
       "mean_difference_shared_minus_not_host_jaccard":
             float(np.mean(yplants[original_x==1])-np.mean(yplants[original_x==0])),
       "mean_difference_shared_minus_not_butterfly_jaccard":
             float(np.mean(ybutterflies[original_x==1])-np.mean(ybutterflies[original_x==0])),
    }
    tests={}
    for gname,gv in groups.items():
        ids=defaultdict(list)
        for i,b in enumerate(gv):ids[b].append(i)
        movable=[inds for inds in ids.values()
                 if len({tuple(yrs[i]) for i in inds})>1]
        rng=random.Random(SEED)
        obs_host=slope(resid_plant,original_x)
        obs_butter=slope(resid_butter,original_x)
        nullhost=[];nullbutter=[]
        for iteration in range(4999):
            perm=yrs.copy()
            for inds in movable:
                idx=inds.copy();rng.shuffle(idx)
                perm[inds]=yrs[idx]
            xperm=getx(perm)
            h=slope(resid_plant,xperm)
            b=slope(resid_butter,xperm)
            if h is not None and b is not None:
                nullhost.append(h);nullbutter.append(b)
        def result(obs,draws):
            if obs is None or not draws:return {"status":"NO_IDENTIFIABLE_SHARED_HISTORY_EFFECT"}
            draws.sort()
            mean=float(np.mean(draws))
            return {
                "adjusted_shared_history_coefficient":obs,
                "null_mean":mean,
                "difference_vs_permutation_mean":obs-mean,
                "null_q025":draws[int(.025*(len(draws)-1))],
                "null_q975":draws[int(.975*(len(draws)-1))],
                "p_two_sided_descriptive":(1+sum(abs(x-mean)>=abs(obs-mean)-1e-12 for x in draws))/(len(draws)+1),
                "p_high_descriptive":(1+sum(x>=obs-1e-12 for x in draws))/(len(draws)+1),
                "draws":len(draws)
            }
        tests[gname]={
          "effective_shuffled_nodes":len({i for arr in movable for i in arr}),
          "exchangeable_blocks":len(movable),
          "plant_host_similarity":result(obs_host,nullhost),
          "butterfly_new_opportunity_similarity":result(obs_butter,nullbutter),
        }
    return {**base,"history_comparisons":tests}


def main():
    p=argparse.ArgumentParser()
    for k in ("q2-edges","q2-audit","empire-islands","region-area","wgsrpd3-geojson","output"):
        p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();a.output.parent.mkdir(parents=True,exist_ok=True)
    edges=list(readcsv(a.q2_edges))
    audit=json.loads(a.q2_audit.read_text(encoding="utf-8"))
    if len(edges)!=audit["scope"]["introduced_butterfly_plant_island_edges"]:
        raise RuntimeError("Detailed original WCVP-HOSTS edges differ from frozen Q2 source audit")
    plants=defaultdict(set);butterflies=defaultdict(set)
    for row in edges:
        island=row["island_code"]
        plants[island].add(row["plant_id"])
        butterflies[island].add(row["species"])
    if len(plants)!=122 or sum(len(bs) for bs in butterflies.values())!=5674:
        raise RuntimeError("Island x butterfly panel drift")
    hist={r["area_code_l3"]:r for r in readcsv(a.empire_islands)}
    areas={r["code"]:r for r in readcsv(a.region_area)}
    coords=region_representative_lonlat(a.wgsrpd3_geojson)
    if len(hist)!=125:raise RuntimeError("Colonial historical island cohort drift")
    rows=[]
    for code in plants:
        entry=hist[code]
        if code not in coords or code not in areas:raise RuntimeError("Unmapped geographic island "+code)
        years={e:float(entry[e+"_weighted_exposure_years"]) for e in EMPIRES}
        rows.append({"code":code,"level1":entry["level1"],
                     "coverage":float(entry["lenzner_polygon_coverage"]),
                     "area_km2":float(areas[code]["area_km2"]),
                     "history":years,"strict":entry["strict_historical_empire"]=="1"})
    def get(mincoverage,strict=False):
        return sorted([r for r in rows if r["coverage"]>=mincoverage and
                any(yrs>=50 for yrs in r["history"].values()) and
                (not strict or r["strict"])],key=lambda r:r["code"])
    geod=Geod(ellps="WGS84")
    groups={
      "primary_geographic_cover80_50yr":get(.8),
      "sensitivity_cover50_50yr":get(.5),
      "sensitivity_strict_single_empire_cover80_50yr":get(.8,True)
    }
    result={
      "schema":"chocho_q2_historic_empire_shared_resource_similarity_v0.1",
      "status":"EXPLORATORY_PAIRWISE_HISTORICAL_RESOURCE_GEOGRAPHY",
      "sources":{"WCVP":"65bed76bae9d644ccb6ad200c05f9f5071d89e05",
                 "HOSTS":"808e0b869f9ec1adf8efff87cf6a395adda103e0",
                 "Lenzner":"7aae72aadd39ae3827ece6935c6740dc0aace030",
                 "WGSRPD":"52da7828aba9d461dd133c27b3bd7a4407161f54"},
      "shared_history":"Both islands have >=50 area-weighted former colonizer years for at least one of 4 European empires",
      "butterfly_outcome":"Jaccard overlap of fixed 239-butterfly species for which introduced larval hosts create new WGSRPD3 potential opportunity in respective islands",
      "plant_outcome":"Jaccard overlap of host species whose introduced distribution supports a new butterfly opportunity",
      "controls":"OLS descriptive adjustment: great-circle distance, absolute log-area difference, same Level1, absolute log host richness difference, absolute log butterfly richness difference",
      "null":"Shuffle full four-empire history vector among island nodes, with global, Level1 or Level1×area-half restriction; all pairs are NON-INDEPENDENT",
      "cohorts":{name:pair_model(cohort,plants,butterflies,coords,geod) for name,cohort in groups.items()},
      "limitations":[
         "Exploratory avenue developed after absence of support for simple metropolitan native-source prediction.",
         "Shared empire assignments do not demonstrate historic plant movement between islands; support region-to-region association only.",
         "Historical records cover only four European empires; French and other colonial histories are absent.",
         "Fixed plant/host interaction network and region-level native/introduced WCVP data imply resource opportunity, not butterfly presence.",
         "WCVP source sites and introduced host records are unequal in size, sampling, and taxonomic resolution.",
         "Descriptive permutations of nodes are not an independent pair test or causal identification; residual spatial structure is untested."
      ]
    }
    a.output.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(result,ensure_ascii=False,indent=2))
if __name__=="__main__":
    main()
