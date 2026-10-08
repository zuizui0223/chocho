#!/usr/bin/env python3
"""Post-hoc route-month stage composition audit of Erickson et al. 2025.

Egg and larval instar reports are cross-sectional. No cohort survival is measured.
"""
from __future__ import annotations
import argparse,csv,json,random,math
from collections import defaultdict,Counter
from pathlib import Path
STAGES=["Eggs","Instar_1","Instar_2","Instar_3","Instar_4","Instar_5"]
SEASONS={"winter":{1,2,3},"spring":{4,5,6},"summer":{7,8,9},"fall":{10,11,12}}

def valid(x):
    try:
        v=float(x)
        return v if math.isfinite(v) and v>=0 else None
    except (ValueError,TypeError):return None

def read(path):
    from datetime import datetime
    cleaned=[];source=0;removed=Counter()
    with path.open(newline="",encoding="utf-8-sig") as f:
        for r in csv.DictReader(f):
            source+=1
            try:d=datetime.strptime(r["Date"],"%m/%d/%y")
            except ValueError:raise RuntimeError("original survey date mismatch")
            nums={k:valid(r.get(k)) for k in STAGES+["Num_plants"]}
            if any(v is None for v in nums.values()) or nums["Num_plants"]<=0:
                removed["unusable_stage_or_plant_count"]+=1;continue
            if r["Native_status"] not in ("Exotic","Native"):raise RuntimeError("unexpected native status")
            e={**nums,"year":d.year,"month":d.month,"route":int(r["Route"]),
               "native_status":r["Native_status"],"milkweed":r["Milkweed_sp"],
               "Too_far_accessible":(str(r.get("Too_far") or "").strip().lower() in ("","no"))}
            e["late"]=e["Instar_4"]+e["Instar_5"]
            e["early"]=e["Instar_1"]+e["Instar_2"]
            e["all_larvae"]=sum(e[f"Instar_{j}"] for j in range(1,6))
            cleaned.append(e)
    if source!=4518:raise RuntimeError(f"frozen 4518 rows drift: {source}")
    return cleaned,{"source_count":source,"clean_count":len(cleaned),"excluded":dict(removed)}

def tables(rows,year,endpoint):
    grouped=defaultdict(lambda:Counter())
    for r in rows:
        if r["year"]==year:
            key=(r["route"],r["month"],r["native_status"])
            for metric in ["Eggs",endpoint,"Num_plants"]:grouped[key][metric]+=r[metric]
    pairs={}
    for (route,month,status),numbers in grouped.items():pairs.setdefault((route,month),{})[status]=numbers
    active=[];all_pairs=0
    for (route,month),group in sorted(pairs.items()):
        if "Exotic" not in group or "Native" not in group:continue
        all_pairs+=1
        a=group["Exotic"][endpoint];b=group["Exotic"]["Eggs"]
        c=group["Native"][endpoint];d=group["Native"]["Eggs"]
        if a+b>0 and c+d>0:active.append((route,month,a,b,c,d))
    return active,all_pairs

def or_mh(active):
    if not active:return None
    numerator=denominator=0.
    for _,_,a,b,c,d in active:
        n=a+b+c+d
        numerator+=a*d/n
        denominator+=b*c/n
    if denominator==0:return None
    return numerator/denominator

def quantile(values,p):
    values=sorted(values)
    z=p*(len(values)-1);i=int(z);j=math.ceil(z)
    return values[i] if i==j else values[i]*(j-z)+values[j]*(z-i)

def season(active,months):
    return [s for s in active if s[1] in months]

def block_bootstrap(active,routes,draws,seed):
    rng=random.Random(seed);routeids=sorted(set(routes))
    groups={r:[x for x in active if x[0]==r] for r in routeids}
    collected={k:[] for k in ("overall","spring","fall","spring_over_fall")}
    for _ in range(draws):
        chosen=[rng.choice(routeids) for _ in routeids]
        draw=[item for r in chosen for item in groups[r]]
        all_or=or_mh(draw)
        sp=or_mh(season(draw,SEASONS["spring"]))
        fa=or_mh(season(draw,SEASONS["fall"]))
        for key,val in (("overall",all_or),("spring",sp),("fall",fa),("spring_over_fall",sp/fa if sp is not None and fa not in (None,0) else None)):
            if val is not None and math.isfinite(val):collected[key].append(val)
    out={}
    for key,x in collected.items():
        out[key]={"finite_draws":len(x),"total_draws":draws,
                  "ci95":[quantile(x,.025),quantile(x,.975)] if len(x)>=draws*.8 else None,
                  "median":quantile(x,.5) if x else None,
                  "fraction_below_1":sum(v<1 for v in x)/len(x) if x else None}
    return out

def analyze(rows,year,endpoint="late",bootstrap=False):
    active,n_pairs=tables(rows,year,endpoint)
    totals=defaultdict(lambda:Counter())
    for r in rows:
        if r["year"]==year:
            for k in ["Eggs",endpoint,"Num_plants"]:totals[r["native_status"]][k]+=r[k]
    a,b=totals["Exotic"][endpoint],totals["Exotic"]["Eggs"]
    c,d=totals["Native"][endpoint],totals["Native"]["Eggs"]
    seasons={name:{"active_paired_route_months":len(season(active,months)),
        "exotic_vs_native_late_to_eggs_MH_OR":or_mh(season(active,months)),
        "stage_counts":{k:int(sum(x[i] for x in season(active,months))) for k,i in (("exotic_late",2),("exotic_eggs",3),("native_late",4),("native_eggs",5))}}
        for name,months in SEASONS.items()}
    sp=seasons["spring"]["exotic_vs_native_late_to_eggs_MH_OR"]
    fa=seasons["fall"]["exotic_vs_native_late_to_eggs_MH_OR"]
    out={"year":year,"endpoint":endpoint,"patch_observations":sum(r["year"]==year for r in rows),
      "route_months_both_plant_origins":n_pairs,"active_route_month_strata":len(active),
      "unstratified_stage_OR":a*d/(b*c) if b*c else None,
      "unstratified_counts":{"exotic_eggs":int(b),"exotic_late":int(a),"native_eggs":int(d),"native_late":int(c),
           "exotic_plant_reports":int(totals["Exotic"]["Num_plants"]),"native_plant_reports":int(totals["Native"]["Num_plants"])},
      "stratified_MH_OR":or_mh(active),"season":seasons,
      "spring_to_fall_OR_ratio":sp/fa if sp is not None and fa not in (None,0) else None}
    if bootstrap:
        routes=sorted({r["route"] for r in rows if r["year"]==year})
        out["route_bootstrap"]=block_bootstrap(active,routes,4999,20261008)
        out["leave_one_route_out"]={str(r):or_mh([x for x in active if x[0]!=r]) for r in routes}
    return out,active

def main():
    p=argparse.ArgumentParser()
    for name in ("source_csv","protocol_json","output_json","output_matched_csv"):
        p.add_argument("--"+name.replace("_","-"),required=True,type=Path)
    args=p.parse_args()
    protocol=json.loads(args.protocol_json.read_text(encoding="utf-8"))
    if protocol["schema"]!="chocho_urban_monarch_stage_host_origin_v0.1":
        raise RuntimeError("unexpected protocol")
    rows,audit=read(args.source_csv)
    primary,matched=analyze(rows,2022,bootstrap=True)
    if primary["patch_observations"]!=3001 or primary["route_months_both_plant_origins"]!=131 or primary["active_route_month_strata"]!=88:
        raise RuntimeError("frozen 2022 matched panel drift")
    sensitive={end:analyze(rows,2022,endpoint=end)[0] for end in ("early","all_larvae","Instar_5")}
    access,_=analyze([r for r in rows if r["Too_far_accessible"]],2022)
    partial={str(y):analyze(rows,y)[0] for y in (2023,2024)}
    species=defaultdict(set);counts=Counter()
    for r in rows:
        if r["year"]==2022:
            species[r["milkweed"]].add(r["native_status"])
            counts[(r["milkweed"],r["native_status"])]+=r["Num_plants"]
    if any(len(x)>1 for x in species.values()):
        raise RuntimeError("plant origin not perfectly taxon-specific: interpretation must change")
    result={
      "schema":"chocho_monarch_stage_origin_route_month_v0.1",
      "status":"EXPLORATORY_REANALYSIS_OF_PUBLISHED_COUNTS_NOT_TRACKED_COHORT_SURVIVAL",
      "source":"Erickson Schultz Crone 2025 DOI 10.1002/ecs2.70259; Figshare 10.6084/m9.figshare.25648644.v1",
      "protocol":"docs/exploratory/URBAN_MILKWEED_ORIGIN_STAGE_MATCHED_PROTOCOL_V01.json",
      "source_audit":audit,"primary_2022":primary,
      "sensitivity_accessible_patches":{k:access[k] for k in ("patch_observations","route_months_both_plant_origins","active_route_month_strata","stratified_MH_OR")},
      "alternate_larval_stages":{k:{"OR":v["stratified_MH_OR"],"active_strata":v["active_route_month_strata"]} for k,v in sensitive.items()},
      "incomplete_other_years":{k:{"OR":v["stratified_MH_OR"],"active_strata":v["active_route_month_strata"],"paired_route_months":v["route_months_both_plant_origins"]} for k,v in partial.items()},
      "plant_species_identity":{"plant_species_in_both_native_statuses":sum(len(x)>1 for x in species.values()),
          "top_species_by_plant_reports":[{"species":s,"status":n,"reported_plants":int(v)} for (s,n),v in counts.most_common(8)]},
      "inference_limit":["Egg and instar observations are different cross-sectional individuals; the ratio is NOT egg-to-larva survival.",
      "Native versus exotic categories fully confounded with plant species identity.",
      "15 route cluster bootstrap; garden and individual plant detection and independence may remain unaddressed.",
      "One full season (2022); 2023 incomplete and 2024 partial, NOT independent replicates.",
      "The published paper already described seasonal monarch egg-to-larva patterns; this is posthoc origin-stratified reanalysis.",
      "Do not change main GEB manuscript."]
    }
    args.output_json.parent.mkdir(parents=True,exist_ok=True)
    args.output_json.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    args.output_matched_csv.parent.mkdir(parents=True,exist_ok=True)
    with args.output_matched_csv.open("w",newline="",encoding="utf-8") as f:
        w=csv.writer(f);w.writerow(["Route","month","exotic_late","exotic_eggs","native_late","native_eggs"])
        w.writerows(matched)
    print(json.dumps({"source_clean":audit,"primary":primary["stratified_MH_OR"],
          "season_or":{k:v["exotic_vs_native_late_to_eggs_MH_OR"] for k,v in primary["season"].items()},
          "boot":primary["route_bootstrap"],"accessibility":result["sensitivity_accessible_patches"],
          "stages":result["alternate_larval_stages"]},indent=2),flush=True)
if __name__=="__main__":main()
