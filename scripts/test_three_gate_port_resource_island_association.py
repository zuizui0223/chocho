#!/usr/bin/env python3
"""Descriptive 'port as plant-resource supply faucet' test on corrected 126 islands.

Tests modern indexed World Port Index count against the *plant-mediated
potential host-resource expansion* of 239 butterflies. No butterfly GBIF data
are used in this test. Modern ports are not an historical shipping/plant import
measure, and an association cannot identify ship-borne introduction.
"""
from __future__ import annotations
import argparse,csv,json,math
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr

SEED=20261102
N_PERM=4999

def readcsv(p):
    with p.open(encoding="utf-8",newline="") as f:
        return list(csv.DictReader(f))

def export(path,rows):
    if not rows:raise RuntimeError("No rows "+str(path))
    with path.open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader();w.writerows(rows)

def group_moments(a, outcome):
    if not len(a):return {"n":0}
    vals=np.asarray([x[outcome] for x in a],dtype=float)
    return {"n":len(a),"mean":float(np.mean(vals)),"median":float(np.median(vals)),
            "min":float(np.min(vals)),"max":float(np.max(vals))}

def prepare_projection(data, *, native_control=False):
    n=len(data)
    levs=sorted({d["level1"] for d in data})
    nuisance=[np.ones(n),np.array([d["log_area"] for d in data])]
    names=["intercept","ln_area_km2"]
    if native_control:
        nuisance.append(np.array([math.log1p(d["native"]) for d in data]))
        names.append("ln1p_native_host_resource_units")
    for l in levs[1:]:
        nuisance.append(np.array([int(d["level1"]==l) for d in data],float))
        names.append("level1=="+l)
    Z=np.column_stack(nuisance)
    rank=np.linalg.matrix_rank(Z)
    if rank!=Z.shape[1]:
        raise RuntimeError("Nuisance covariate matrix rank deficient "+str((rank,Z.shape)))
    Q,_=np.linalg.qr(Z,mode="reduced")
    def residualize(vector):
        x=np.asarray(vector,dtype=float)
        return x-Q@(Q.T@x)
    return residualize,names

def regression(data, outcome, *, native_control=False, exposure="log2ports",
               permutations=True, rngseed=SEED):
    n=len(data)
    y=np.array([d[outcome] for d in data],float)
    x=np.array([d[exposure] for d in data],float)
    residualize,names=prepare_projection(data,native_control=native_control)
    yr=residualize(y)
    def slope(x):
        xr=residualize(x)
        den=float(xr@xr)
        return None if den<=1e-9 else float((xr@yr)/den)
    obs=slope(x)
    if obs is None:
        return {"status":"NON_IDENTIFIABLE_EXPOSURE","n":n}
    out={"n":n,"exposure":exposure,"outcome":outcome,
         "controls":names,"coefficient":obs,
         "outcome_mean":float(np.mean(y)),"exposure_mean":float(np.mean(x))}
    if not permutations:return out
    areas=np.array([d["log_area"] for d in data])
    med_by_level={l:float(np.median([d["log_area"] for d in data if d["level1"]==l]))
                  for l in sorted({d["level1"] for d in data})}
    blocked={
       "global":["all"]*n,
       "within_level1":[d["level1"] for d in data],
       "within_level1_and_within_level_area_half":[
          d["level1"]+"|"+str(int(d["log_area"]>med_by_level[d["level1"]])) for d in data],
    }
    tests={}
    for blockname,labels in blocked.items():
        groups=defaultdict(list)
        for i,label in enumerate(labels):groups[label].append(i)
        moved=[indices for indices in groups.values()
               if len(indices)>1 and len({float(x[i]) for i in indices})>1]
        touched={i for indices in moved for i in indices}
        rng=np.random.default_rng(rngseed+len(blockname)+n)
        null=[]
        failures=0
        for i in range(N_PERM):
            xp=x.copy()
            for ind in moved:
                ind=np.array(ind,dtype=int)
                xp[ind]=x[rng.permutation(ind)]
            beta=slope(xp)
            if beta is None:failures+=1
            else:null.append(beta)
        if not null:
            tests[blockname]={"status":"NO_EXCHANGEABLE_SUPPORT",
                              "exchangeable_islands":len(touched),
                              "exchangeable_blocks":len(moved)}
            continue
        arr=np.array(null,float)
        center=float(np.mean(arr))
        tests[blockname]={
           "status":"DESCRIPTIVE_CONDITIONAL_PORT_PERMUTATION",
           "permutations":len(arr),"excluded_degenerate_permutations":failures,
           "exchangeable_islands":len(touched),"exchangeable_blocks":len(moved),
           "exposure_coefficient_observed":obs,
           "null_mean":center,
           "observed_minus_null_mean":obs-center,
           "null_q025":float(np.quantile(arr,.025)),
           "null_q975":float(np.quantile(arr,.975)),
           "one_sided_positive_p":(1+int(np.sum(arr>=obs-1e-12)))/(len(arr)+1),
           "two_sided_centered_p":(1+int(np.sum(np.abs(arr-center)>=abs(obs-center)-1e-12)))/(len(arr)+1),
           "notes":"No causal or spatial independence guarantee; port counts are shuffled among measured regions."
        }
    out["permutations"]=tests
    return out

def main():
    p=argparse.ArgumentParser()
    for x in ("q1-regions","q1-audit","wpi-island-regions","island-area",
              "wpi-audit","out-dir"):
        p.add_argument("--"+x,required=True,type=Path)
    a=p.parse_args();a.out_dir.mkdir(parents=True,exist_ok=True)
    q1=json.loads(a.q1_audit.read_text(encoding="utf-8"))
    wpi=json.loads(a.wpi_audit.read_text(encoding="utf-8"))
    by_q1={r["area_code_l3"]:r for r in readcsv(a.q1_regions)}
    by_wpi={r["area_code_l3"]:r for r in readcsv(a.wpi_island_regions)}
    by_area={r["code"]:r for r in readcsv(a.island_area)}
    if (len(by_q1),len(by_wpi))!=(369,369):raise RuntimeError("WGS/WPI source unit drift")
    if len(q1["classification"]["strict_island_codes"])!=126 or q1["panel_total"]["n_panel_mismatches"]!=0:
        raise RuntimeError("Expected source-audited corrected 126-island Q1")
    if (q1["groups"]["island"]["native_butterfly_region_units"],
        q1["groups"]["island"]["contemporary_butterfly_region_units"],
        q1["groups"]["island"]["added_butterfly_region_units"])!=(4829,10570,5741):
        raise RuntimeError("Q1 island resource counts drift")
    if wpi["corrected_indexed_WPI_island_ports"]!=982:
        raise RuntimeError("WPI corrected modern port count drift")
    all_islands={k:v for k,v in by_q1.items() if v["class"]=="island"}
    if len(all_islands)!=126:raise RuntimeError("Q1 island class incorrect")
    if sum(by_wpi[k]["geographic_group"]=="island" for k in all_islands)!=126:
        raise RuntimeError("WPI island class inconsistent")
    active=[]
    unavailable=[]
    for code, r in sorted(all_islands.items()):
        cont=int(r["contemporary"]); n=int(r["native"]); added=int(r["added"])
        port=int(by_wpi[code]["WPI_port_count"])
        if cont<=0:
            unavailable.append(code)
            continue
        if code not in by_area:
            raise RuntimeError("No Q1 geodesic area for active island "+code)
        area=float(by_area[code]["area_km2"])
        if area<=0:raise RuntimeError("Nonpositive source area "+code)
        if cont!=n+added:
            raise RuntimeError("Resource counts not additive "+code)
        if abs(added/cont-float(by_area[code]["share"]))>1e-8:
            raise RuntimeError("GIS area source resource share mismatch "+code)
        active.append({
           "code":code,"name":r["region_name"],"level1":r["level1"],
           "area_km2":area,"log_area":math.log(area),
           "native":n,"contemporary":cont,"added":added,
           "added_share":added/cont,"log_added":math.log1p(added),
           "ports":port,"has_port":int(port>0),
           "log2ports":math.log2(1+port),
        })
    if len(active)!=123 or len(unavailable)!=3:
        raise RuntimeError(f"Active source cohort changed {len(active)}; zeros {unavailable}")
    if sum(x["ports"] for x in active)!=982:
        raise RuntimeError("Ports present in non-active islands; update audit instead of assuming")
    cohorts={
       "all_active_123":active,
       "native_host_baseline_positive":[x for x in active if x["native"]>0],
       "contemporary_at_least_75":[x for x in active if x["contemporary"]>=75],
       "exclude_pacific_antarctic":[x for x in active if x["level1"] not in ("6","9")],
       "exclude_single_highest_port_region":[x for x in active
          if x["ports"]<max(z["ports"] for z in active)],
    }
    summary={
       "schema":"chocho_q4_port_index_and_host_resource_faucet_v0.1",
       "status":"SOURCE_AUDITED_DESCRIPTIVE_2023_WPI_VS_CURRENT_HOST_RESOURCES",
       "design":{
           "population":"126 strict WGSRPD3 island regions, 123 with measurable contemporary resources",
           "response_primary":"newly introduced host-resource opportunity / contemporary opportunity, species×region",
           "exposure":"log2(1+World Port Index count, 2023 published snapshot)",
           "unit":"WGSRPD3 regions, NOT individual islands or port catchments",
           "covariates":"log region polygon area and fixed Level1 botanical biogeographic regions",
           "source_source":"Pinned 239 butterfly GEB WCVP×HOSTS, official WGSRPD, public 3,630 WPI port catalogue",
           "permutations_per_test":N_PERM,"seed":SEED,
       },
       "baseline":{
          "all_islands":len(all_islands),
          "active_islands":len(active),
          "zero_contemporary_resource_islands":unavailable,
          "islands_with_indexed_port":sum(x["has_port"] for x in active),
          "islands_without_indexed_port":sum(not x["has_port"] for x in active),
          "all_island_indexed_ports":sum(x["ports"] for x in active),
          "port_count_quantiles":[float(v) for v in np.quantile(
             [x["ports"] for x in active],[0,.25,.5,.75,.9,1])],
          "port_zero_added_share":group_moments([x for x in active if x["ports"]==0],"added_share"),
          "port_positive_added_share":group_moments([x for x in active if x["ports"]>0],"added_share"),
          "spearman_port_count_versus_added_share":{
            "rho":float(spearmanr([x["ports"] for x in active],[x["added_share"] for x in active]).statistic),
            "p_naive_nonspatial":float(spearmanr([x["ports"] for x in active],[x["added_share"] for x in active]).pvalue)
          },
       },
       "cohorts":{},
       "limitations":[
         "No causal statement on ships or plant import is warranted from port counts and current WCVP host opportunities.",
         "Present-day 2023 WPI port inventory cannot assign exposure to eighteenth-/nineteenth-century introductions.",
         "WPI catalog coverage is incomplete and biased; zero indexed ports are not proof of no maritime access.",
         "Land area, economy, plant introduction effort, urbanization, tourism, climate and sampling are incompletely controlled.",
         "Area- and Level1-stratified permutations do not remove phylogenetic, spatial or colonial-history dependence.",
         "Each potential butterfly region uses a binary accepted HOSTS link, not an empirical larval establishment probability.",
       ],
    }
    for i,(key,frame) in enumerate(cohorts.items()):
        mainfit=regression(frame,"added_share",rngseed=SEED+i*100,
                            native_control=False,exposure="log2ports")
        withnative=regression(frame,"added_share",rngseed=SEED+i*100+11,
                              native_control=True,exposure="log2ports",permutations=False)
        absolute=regression(frame,"log_added",rngseed=SEED+i*100+23,
                            native_control=True,exposure="log2ports",permutations=False)
        binary=regression(frame,"added_share",rngseed=SEED+i*100+37,
                          native_control=False,exposure="has_port",permutations=False)
        summary["cohorts"][key]={
          "n":len(frame),"positive_WPI_regions":sum(x["ports"]>0 for x in frame),
          "primary_added_share":mainfit,
          "secondary_native_resource_adjusted":withnative,
          "secondary_log_absolute_added":absolute,
          "secondary_binary_any_indexed_port":binary,
        }
    primary=summary["cohorts"]["all_active_123"]["primary_added_share"]
    # Leave-one-region-out robustness keeps the same formula; excludes only
    # one island at a time and does not create independent cases.
    loo=[]
    for x in active:
        drop=[z for z in active if z["code"]!=x["code"]]
        res=regression(drop,"added_share",native_control=False,
                       exposure="log2ports",permutations=False)
        loo.append({"excluded_code":x["code"],"coefficient":res.get("coefficient")})
    v=[x["coefficient"] for x in loo if x["coefficient"] is not None]
    summary["leave_one_island_out"]={
      "n_valid":len(v),"n_positive":sum(z>0 for z in v),
      "min":min(v),"max":max(v),
      "largest_positive_coefficient_exclusions":sorted(loo,key=lambda r:r["coefficient"],reverse=True)[:5],
      "lowest_coefficient_exclusions":sorted(loo,key=lambda r:r["coefficient"])[:5],
      "range_not_CI":True
    }
    export(a.out_dir/"region_port_resource_join.csv",active)
    (a.out_dir/"island_port_host_resource_association_audit.json").write_text(
        json.dumps(summary,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({"baseline":summary["baseline"],
                      "main_added_share":primary,
                      "leave_one_out":summary["leave_one_island_out"]},
                     ensure_ascii=False,indent=2))
if __name__=="__main__":main()
