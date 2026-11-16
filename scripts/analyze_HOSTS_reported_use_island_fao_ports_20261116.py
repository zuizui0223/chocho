#!/usr/bin/env python3
"""Location-indexed NHM HOSTS records × WCVP regional origin × FAO crops × WPI ports.

Main inference unit: resolved original HOSTS Location, never a raw row as an
independent biological island. Source locations are a literature catalogue,
not verified field larval feeding or insect occupancy.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import mannwhitneyu

SOURCE_HOSTS_GIT_BLOB = "06749c04bd5b7a754d6f7717a31c01c51357cae8"
SOURCE_HOSTS_COMMIT = "808e0b869f9ec1adf8efff87cf6a395adda103e0"
SOURCE_WCVP_COMMIT = "65bed76bae9d644ccb6ad200c05f9f5071d89e05"
SEED = 20261116


def read(path: Path):
    with path.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def hash_file(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def simple(x) -> str:
    return re.sub(r"\s+", " ", str(x or "").strip())


def name_binomial(x) -> str:
    p=simple(x).split()
    if len(p)<2 or not re.fullmatch("[A-Z][a-zA-Z-]+",p[0]):
        return ""
    if not re.fullmatch("[a-z][a-z-]+",p[1]):
        return ""
    return " ".join(p[:2])


def location_origin(plantid: str, codes: list[str], wcvp: dict) -> tuple[str,str]:
    values=[]
    for c in codes:
        flags=wcvp.get((plantid,c),set())
        if not flags:
            return "unresolved_WCVP_location",f"missing_{c}"
        if len(flags)>1:
            return "unresolved_WCVP_location",f"mixed_origin_same_region_{c}"
        values.append(next(iter(flags)))
    if len(set(values))>1:
        return "unresolved_WCVP_location","discordant_multi_region_origin"
    return ("introduced" if values[0]=="1" else "native"),"all_source_regions_concordant"


def build_rows(args):
    geo=read(args.locations)
    ports=read(args.ports)
    areas=read(args.areas)
    full239=read(args.butterflies239)
    names=read(args.host_resolution)
    wrows=read(args.wcvp)
    crop=json.loads(args.fao.read_text(encoding="utf8"))
    if crop.get("schema")!="chocho_fao_crop_taxa_v0.1":raise RuntimeError("Wrong FAO source taxonomy")
    port={r["area_code_l3"]:r for r in ports}
    area={r["area_code_l3"]:float(r["area_km2"]) for r in areas}
    if len(ports)!=369 or len(port)!=len(ports):raise RuntimeError("WPI region universe differs from source")
    species239={x["species"] for x in full239}
    if len(species239)!=239:raise RuntimeError("Source butterfly panel not fixed at 239")
    loc={}
    used_regions=set()
    for x in geo:
        k=x["HOSTS_location"]
        if k in loc:raise RuntimeError("Duplicate original HOSTS locality mapping: "+k)
        codes=x["WGSRPD3_codes"].split("|")
        if not codes:raise RuntimeError("Missing locality WGSRPD3 codes")
        if x["mapping_tier"] not in ("exact","composite","approximate","provisional"):
            raise RuntimeError("Unknown mapping tier "+x["mapping_tier"])
        for code in codes:
            if code not in port or code not in area:
                raise RuntimeError(f"Missing WPI/area source region for {k}: {code}")
            if port[code]["geographic_group"]!=x["comparison_group"]:
                raise RuntimeError("Mismatched WCVP geographic group and original locality "+k)
            used_regions.add(code)
        if len(codes)>1 and x["mapping_tier"] not in ("composite","provisional"):
            raise RuntimeError("Multi-region HOSTS source cannot be labeled exact: "+k)
        if x["comparison_group"] not in ("island","mainland"):
            raise RuntimeError("Unknown island status")
        if x["comparison_group"]=="mainland" and x["island_origin_class"]!="mainland":
            raise RuntimeError("Mainland geological class differs")
        if x["comparison_group"]=="island" and x["island_origin_class"] not in (
                "oceanic","continental","mixed_unresolved"):
            raise RuntimeError("Unreviewed island geology")
        values={int(port[c]["WPI_port_count"]) for c in codes}
        loc[k]={**x,"codes":codes,"level1":next(iter({port[c]["level1"] for c in codes}))
                if len({port[c]["level1"] for c in codes})==1 else "mixed",
                "WPI_ports":sum(int(port[c]["WPI_port_count"]) for c in codes),
                "region_area_km2":sum(area[c] for c in codes)}
    if len(loc)<50:raise RuntimeError("Too few curated original HOSTS localities")
    resolution={}
    for x in names:
        key=x["input_host_name"]
        if key in resolution: raise RuntimeError("Duplicate WCVP taxonomic resolution")
        resolution[key]=x
    wcvp=defaultdict(set)
    for x in wrows:
        key=(x["accepted_plant_name_id"],x["area_code_l3"])
        if x["introduced"] not in ("0","1"):
            raise RuntimeError("Unknown WCVP origin flag")
        wcvp[key].add(x["introduced"])
    exact=set(crop["exact_binomials"])
    genus_wc=set(crop["genus_wildcards"])
    dedup_keys=set()
    output=[]
    stat=Counter()
    with args.hosts.open(newline="",encoding="utf-8-sig") as f:
        for x in csv.DictReader(f):
            stat["source_rows"]+=1
            raw_location=simple(x.get("Location"))
            if raw_location not in loc:
                stat["not_strictly_mapped"]+=1
                continue
            L=loc[raw_location]
            insect=simple(x.get("Insect Genus"))+" "+simple(x.get("Insect Species"))
            input_name=simple(x.get("Hostplant Genus"))+" "+simple(x.get("Hostplant Species"))
            r=resolution.get(input_name)
            if not r or r["match_status"]!="matched":
                stat["host_taxon_unmatched_or_ambiguous"]+=1
                continue
            if len(insect.split())!=2:
                stat["insect_name_incomplete"]+=1
                continue
            plantid=r["accepted_plant_name_id"]
            plant=name_binomial(r.get("accepted_name") or "")
            rawbin=name_binomial(input_name)
            origin,reason=location_origin(plantid,L["codes"],wcvp)
            if origin=="unresolved_WCVP_location":
                stat["unresolved_origin_"+reason]+=1
                continue
            lab=simple(x.get("Lab Rearing"))
            crop_exact=(plant in exact or rawbin in exact)
            crop_genus=any(z.split(" ")[0] in genus_wc for z in (plant,rawbin) if z)
            row={
                "HOSTS_row_id":x.get("_id",""),
                "location":raw_location,
                "insect_species":insect,
                "accepted_plant_id":plantid,
                "plant":plant,
                "origin":origin,
                "rearing":lab,
                "focal_butterfly239":insect in species239,
                "crop_exact":crop_exact,
                "crop_genus_wildcard":crop_genus,
                "mapping_tier":L["mapping_tier"],
                "comparison_group":L["comparison_group"],
                "origin_class":L["island_origin_class"],
                "codes":"|".join(L["codes"]),
                "WPI_ports":L["WPI_ports"],
                "region_area_km2":L["region_area_km2"],
                "level1":L["level1"],
            }
            dup=(row["HOSTS_row_id"],input_name,insect,raw_location)
            if dup in dedup_keys:
                stat["duplicate_same_source_id_plant_insect_location"]+=1
                continue
            dedup_keys.add(dup)
            output.append(row)
            stat["linked_unambiguous_origin_rows"]+=1
    return output,loc,dict(stat),{
        "HOSTS_sha256":hash_file(args.hosts),
        "wcvp_origin_sha256":hash_file(args.wcvp),
        "host_resolution_sha256":hash_file(args.host_resolution),
        "FAO_sha256":hash_file(args.fao),
        "locality_crosswalk_sha256":hash_file(args.locations),
        "WPI_ports_sha256":hash_file(args.ports),
        "WGSRPD3_area_sha256":hash_file(args.areas),
        "butterfly239_sha256":hash_file(args.butterflies239),
        "HOSTS_Git_blob":SOURCE_HOSTS_GIT_BLOB,
        "HOSTS_Git_commit":SOURCE_HOSTS_COMMIT,
        "WCVP_Git_commit":SOURCE_WCVP_COMMIT,
        "faocrop_rows":len(exact),
        "faocrop_genus_wildcards":len(genus_wc),
        "mapped_raw_locality_labels":len(loc),
        "mapped_WGSRPD3_codes":len(used_regions),
        "fixed_focal_butterflies":len(species239),
    }


def bootstrap_median(a,b):
    rng=np.random.default_rng(SEED)
    p=np.array(a,dtype=float)
    q=np.array(b,dtype=float)
    if min(len(p),len(q))<2:return None
    vals=np.array([np.median(rng.choice(p,len(p),replace=True))-
                   np.median(rng.choice(q,len(q),replace=True))
                   for _ in range(2999)])
    return [float(x) for x in np.quantile(vals,[.025,.975])]


def unweighted_lm(rows, response="empirical_logit", include_geology=False, bootstrap=1499):
    if len(rows)<10:return {"status":"TOO_FEW_LOCALITIES","n":len(rows)}
    def design(arr):
        base=[[1.,float(r["is_island"]),np.log2(1.+r["WPI_ports"]),
               np.log(r["region_area_km2"]),np.log1p(r["total_record_count"])]
              for r in arr]
        if include_geology:
            for v,r in zip(base,arr):
                v.append(float(r["origin_class"]=="oceanic"))
        return np.asarray(base,float)
    X=design(rows)
    y=np.array([r[response] for r in rows])
    rank=int(np.linalg.matrix_rank(X))
    if rank<X.shape[1]:
        return {"status":"INSUFFICIENT_COMMON_SUPPORT","n":len(rows),"rank":rank,"k":X.shape[1]}
    b=np.linalg.lstsq(X,y,rcond=None)[0]
    labels=["intercept","is_island","log2_1plus_WPI_ports","log_region_area_km2",
            "log1plus_HOSTS_record_count"]
    if include_geology:labels+=["is_oceanic_island_vs_continental"]
    rng=np.random.default_rng(SEED+int(include_geology))
    bs=[]
    for _ in range(bootstrap):
        ix=rng.integers(0,len(rows),size=len(rows))
        xx=X[ix];yy=y[ix]
        if np.linalg.matrix_rank(xx)<X.shape[1]:continue
        bs.append(np.linalg.lstsq(xx,yy,rcond=None)[0])
    interval={}
    if len(bs)>=min(300,bootstrap//2):
        vals=np.asarray(bs)
        interval={k:[float(v) for v in np.quantile(vals[:,i],[.025,.975])]
                  for i,k in enumerate(labels)}
    return {
      "status":"EXPLORATORY_EQUAL_LOCALITY_EMPIRICAL_LOGIT",
      "n":len(rows),"rank":rank,
      "coef_log_odds":{k:float(v) for k,v in zip(labels,b)},
      "source_locality_bootstrap_95pct_descriptive":interval,
      "valid_bootstrap_draws":len(bs),
      "model_unit":"one original HOSTS Location label",
      "model_description":"unweighted empirical logit of (introduced+0.5)/(resolved+1); origin as outcome, \
island and marine-geology indicators and WPI port count plus region area and recording intensity",
      "not_causal":True,
    }


def compute(runrows,localities,threshold):
    grouped=defaultdict(list)
    for x in runrows:
        grouped[x["location"]].append(x)
    records=[]
    for location,vals in sorted(grouped.items()):
        n=len(vals)
        if n<threshold:continue
        sample=vals[0]
        I=sum(x["origin"]=="introduced" for x in vals)
        N=n-I
        if n!=I+N:raise AssertionError("Unknown introduced origin leaked")
        pairs={(x["insect_species"],x["accepted_plant_id"]):x["origin"]
               for x in vals}
        unique_I=sum(x=="introduced" for x in pairs.values())
        records.append({
            "location":location,
            "n":n, "introduced":I, "native":N,
            "introduced_share":I/n,
            "pairs_n":len(pairs),"pairs_I":unique_I,
            "pairs_intro_share":unique_I/len(pairs),
            "WPI_ports":sample["WPI_ports"],
            "region_area_km2":sample["region_area_km2"],
            "is_island":int(sample["comparison_group"]=="island"),
            "origin_class":sample["origin_class"],
            "mapping_tier":sample["mapping_tier"],
            "level1":sample["level1"],
            "total_record_count":n,
            "empirical_logit":math.log((I+.5)/(N+.5)),
        })
    island=[r for r in records if r["is_island"]]
    mainland=[r for r in records if not r["is_island"]]
    oceanic=[r for r in island if r["origin_class"]=="oceanic"]
    continental=[r for r in island if r["origin_class"]=="continental"]
    mixed=[r for r in island if r["origin_class"]=="mixed_unresolved"]
    def med(a,key="introduced_share"):
        return float(np.median([r[key] for r in a])) if a else None
    def total(a,key="introduced"):
        return sum(r[key] for r in a)
    def weighted(a):
        return total(a,"introduced")/total(a,"n") if total(a,"n") else None
    def co(a):
        return {"localities":len(a),"median_share":med(a),
                "record_weighted_share":weighted(a),
                "accepted_HOSTS_rows":total(a,"n"),
                "unique_pair_median_share":med(a,"pairs_intro_share")}
    tests={}
    for key,a,b in [("island_vs_mainland",island,mainland),
                    ("oceanic_vs_continental",oceanic,continental)]:
        if min(len(a),len(b))>=3:
            x=[r["introduced_share"] for r in a]
            y=[r["introduced_share"] for r in b]
            u=mannwhitneyu(x,y,alternative="two-sided",method="auto")
            tests[key]={
                "median_difference_fraction":med(a)-med(b),
                "descriptive_Mann_Whitney_p":float(u.pvalue),
                "locality_resample_median_difference_95pct":bootstrap_median(x,y),
                "not_inferentially_independent_species_records":True,
            }
        else:tests[key]={"status":"FEWER_THAN_THREE_LOCALITIES_IN_GROUP"}
    result={
        "status":"SOURCE_LOCALITY_DESCRIPTIVE",
        "threshold_known_origin_records_per_location":threshold,
        "n_localities":len(records),
        "island":co(island),"mainland":co(mainland),
        "oceanic":co(oceanic),"continental":co(continental),
        "mixed_geology":co(mixed),
        "tests":tests,
        "adjusted_model_island_vs_mainland":unweighted_lm(records,include_geology=False),
        "adjusted_model_with_oceanic_continental":unweighted_lm(
            [r for r in records if r["origin_class"]!="mixed_unresolved"],
            include_geology=True),
        "details":records,
    }
    return result


def main():
    p=argparse.ArgumentParser()
    for name in ("hosts","host-resolution","wcvp","locations","fao","ports","areas",
                 "butterflies239","out-dir"):
        p.add_argument("--"+name,required=True,type=Path)
    a=p.parse_args()
    rows,localities,coverage,source=build_rows(a)
    a.out_dir.mkdir(exist_ok=True,parents=True)
    runs={}
    details=[]
    for cohort in ("fixed239_butterflies","all_HOSTS_lepidoptera"):
        for geography in ("strict_exact_or_composite","expanded_approximate"):
            for cropmode in ("all_hosts","exclude_FAO_exact","exclude_FAO_exact_plus_genus_wildcards"):
                sel=[
                    x for x in rows
                    if (cohort=="all_HOSTS_lepidoptera" or x["focal_butterfly239"])
                    and (geography=="expanded_approximate" or x["mapping_tier"] in ("exact","composite"))
                    and x["mapping_tier"]!="provisional"
                    and x["rearing"]!="Y"
                    and (cropmode=="all_hosts" or not x["crop_exact"])
                    and (cropmode!="exclude_FAO_exact_plus_genus_wildcards" or not x["crop_genus_wildcard"])
                ]
                key="__".join((cohort,geography,cropmode))
                main=compute(sel,localities,10)
                sensitivity=compute(sel,localities,20)
                runs[key]={
                    "eligible_original_HOSTS_rows":len(sel),
                    "source_localities_present":len(set(x["location"] for x in sel)),
                    "min10":{k:v for k,v in main.items() if k!="details"},
                    "min20":{k:v for k,v in sensitivity.items() if k!="details"},
                }
                details.extend({**r,"analysis_key":key,"threshold":10}
                               for r in main["details"])
    primary="fixed239_butterflies__strict_exact_or_composite__all_hosts"
    focal=runs[primary]
    extra=[r for r in rows if r["rearing"]=="Y"]
    audit={
        "schema":"chocho_HOSTS_regionally_introduced_record_use_20261116_v1",
        "status":"SOURCE_ANALYSIS_NOT_FIELD_CONFIRMED_LARVAL_USE",
        "source":source,
        "source_coverage":coverage,
        "total_exact_origin_linked_rows":len(rows),
        "rows_definitely_laboratory_reared_flag_Y":len(extra),
        "primary_key":primary,
        "results":runs,
        "interpretation_limitations":[
            "The HOSTS Location field represents literature-reported geographic attribution, not a mapped larval feeding event, and historical publication/rearing effort is highly biased.",
            "The primary butterfly panel is the 239 original chocho species; all Lepidoptera are an explicit separate contextual sensitivity, not butterfly-only data.",
            "Lab Rearing Y is omitted, but blank Lab Rearing means UNKNOWN, not field-verified feeding. Do not claim observed wild larval consumption.",
            "WCVP v13 introduced/native is current regional status and may not reflect historic introduction at the recording site or at the date of the HOSTS observation.",
            "Exact locality to one or more WGSRPD3 codes is curated; multi-region origin disagreement/missingness excluded rather than imputed.",
            "The oceanic/continental classification is operational, manually curated and geographically heterogeneous, with unresolved arcs and compound islands classified as mixed rather than force-split.",
            "FAO crop species exclusion and explicit FAO spp genus wildcard are different sensitivity scopes, not proof of plant cultivation at the locality.",
            "Raw literature rows share taxa, authors and sources. The independent unit of the descriptive comparison and bootstrap is the location, NOT the HOSTS row.",
            "Current WPI ports do not date introductions or demonstrate transport. Port-adjusted models are locality-level descriptive regressions with unmeasured geography and observation effort.",
            "Country-wide HOSTS locations are not fine enough to prove local host and butterfly co-occurrence; broad/ambiguous names are excluded from main scope.",
            "Restricted 239-species sample may not reproduce a pilot using all moths and butterflies; absence of support is not biological non-use.",
        ],
    }
    (a.out_dir/"host_reported_use_island_oceanic_fao_ports.json").write_text(
        json.dumps(audit,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    if details:
        with (a.out_dir/"per_locality_all_scopes.csv").open("w",newline="",encoding="utf8") as f:
            w=csv.DictWriter(f,fieldnames=list(details[0]))
            w.writeheader()
            for x in details:w.writerow(x)
    print("HOSTS_LOCATION_USE_MODEL",json.dumps({
        "primary":focal,
        "expanded":runs["fixed239_butterflies__expanded_approximate__all_hosts"],
        "all_lepidoptera":runs["all_HOSTS_lepidoptera__strict_exact_or_composite__all_hosts"],
        "crop_primary":runs["fixed239_butterflies__strict_exact_or_composite__exclude_FAO_exact"],
        "model_scope_count":len(runs),
        "source_coverage":coverage,
    },ensure_ascii=False),flush=True)


if __name__=="__main__":main()
