#!/usr/bin/env python3
"""Global island introduced/native life-cycle and growth-form trait association.

Pinned WCVP v13: accepted botanical names + native/introduced WGSRPD3 ranges.
GIFT v3.2: independent, category-level plant lifecycle / herb growth traits.

TEST: Are modern port-indexed WGSRPD3 island regions richer in ANNUAL plant
traits within introduced flora after accounting for total introduced flora and
botanical effort? Native annual share is a negative geographic comparator.

Neither "annual", "herb" nor their conjunction means "nitrophilous" or
"ruderal". Traits are missing not absent. Modern port index is not historical
shipping or causal plant introduction. Data versions are locked below.
"""
from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
import math
import random
import re
import statistics
import time
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
import pyreadr

GIFT_BASE="https://gift.uni-goettingen.de/api/extended/index3.2.php?"
WCVP_COMMIT="65bed76bae9d644ccb6ad200c05f9f5071d89e05"
SOURCE={"lifecycle":"2.1.1","growth_form":"1.2.1"}
B=1999
SEED=20261105

def get_page(query:dict,retries=4):
    u=GIFT_BASE+urllib.parse.urlencode(query)
    err=None
    for k in range(retries):
        try:
            req=urllib.request.Request(u,headers={"User-Agent":"chocho-island-trait-source-research/1.0","Accept":"application/json"})
            with urllib.request.urlopen(req,timeout=55) as response:
                b=response.read(25000000)
            if len(b)>=24000000:raise ValueError("GIFT payload limit, no truncated data allowed")
            j=json.loads(b.decode("utf-8-sig"))
            if not isinstance(j,list):raise ValueError("Expected GIFT source list, got "+str(type(j)))
            return j,hashlib.sha256(b).hexdigest()
        except Exception as e:
            err=str(e)[:300]
            time.sleep(2.5*(k+1))
    raise RuntimeError("GIFT official API request failed: "+u+" "+str(err))

def all_pages(query,step,max_pages,source):
    result=[]
    receipts=[]
    for n in range(max_pages):
        q={**query,"startat":n*step}
        rows,sha=get_page(q)
        for row in rows:
            if not isinstance(row,dict):raise ValueError("Unexpected GIFT trait row type")
        result.extend(rows)
        receipts.append({"source":source,"startat":n*step,"rows":len(rows),"sha256":sha})
        print("GIFT_PINNED_TRAIT_PAGE",source,n*step,len(rows),flush=True)
        if len(rows)<step:
            break
        time.sleep(.18)
    else:
        raise ValueError("GIFT source pagination max_pages reached; incomplete source forbidden")
    if not result:raise ValueError("GIFT source empty "+source)
    return result,receipts

def as_bins(name):
    if not isinstance(name,str):return None
    # Exact canonical binomials only: no fuzzy matching, authors, hybrids or
    # automatic genus-level imputation.
    x=" ".join(name.split())
    return x.casefold() if re.fullmatch(r"[A-Z][a-zA-Z-]+ [a-z][a-z-]+",x) else None

def traits_from_gift(species,rows):
    names_by_work={}
    ambiguous_work=set()
    for r in species:
        try:k=int(r["work_ID"])
        except (TypeError,ValueError):continue
        name=as_bins(r.get("work_species"))
        if name is None:continue
        if k in names_by_work and names_by_work[k]!=name:
            ambiguous_work.add(k)
        names_by_work[k]=name
    for k in ambiguous_work:names_by_work.pop(k,None)
    byname={}
    ambiguous_name=set()
    counts=collections.Counter()
    for trait,records in rows.items():
        for r in records:
            score=r.get("agreement")
            try:score=float(score) if score is not None else 1.0
            except (TypeError,ValueError):continue
            if not math.isfinite(score) or score<.66:
                counts[f"{trait}_low_agreement"]+=1
                continue
            try:work=int(r["work_ID"])
            except (TypeError,ValueError):continue
            name=names_by_work.get(work)
            if name is None:continue
            val=str(r.get("trait_value") or "").strip().casefold()
            allowed={"lifecycle":{"annual","biennial","perennial"},
                     "growth_form":{"herb","tree","shrub"}}
            if val not in allowed[trait]:
                counts[f"{trait}_ambiguous"]+=1
                continue
            x=byname.setdefault(name,{})
            if trait in x and x[trait]!=val:
                ambiguous_name.add((name,trait))
            else:x[trait]=val
            counts[f"{trait}_{val}"]+=1
    for name,trait in ambiguous_name:
        byname[name].pop(trait,None)
    for name in list(byname):
        if not byname[name]:del byname[name]
    return byname,dict(counts),{"ambiguous_GIFT_species_names":len(ambiguous_name),
          "unmatched_work_IDs":len(ambiguous_work),"known_name_keys":len(names_by_work)}

def read_wcvp(npath,dpath):
    names=pyreadr.read_r(str(npath))["wcvp_names"]
    distr=pyreadr.read_r(str(dpath))["wcvp_distributions"]
    cols=list(names.columns)
    required={"plant_name_id","taxon_name","taxon_rank","taxon_status"}
    if not required.issubset(cols):raise ValueError("WCVP names columns changed "+str(cols))
    dcols={"plant_name_id","area_code_l3","introduced","extinct","location_doubtful"}
    if not dcols.issubset(distr.columns):
        raise ValueError("WCVP distributions columns changed "+str(list(distr.columns)))
    names=names.loc[(names.taxon_rank.astype(str).str.casefold()=="species")&
                    (names.taxon_status.astype(str).str.casefold()=="accepted"),
                    ["plant_name_id","taxon_name"]].copy()
    names["plant_name_id"]=names.plant_name_id.astype(str)
    names["taxon_key"]=names.taxon_name.map(as_bins)
    names=names[names.taxon_key.notna()]
    if len(names)<200000:raise ValueError("Unexpected WCVP accepted canonical species size")
    if names.plant_name_id.duplicated().any():
        raise ValueError("WCVP accepted IDs unexpectedly ambiguous")
    return names,distr

def source_table(wcvp_names,distributions,trait,strict_codes):
    names=wcvp_names[["plant_name_id","taxon_key"]].copy()
    data=distributions[["plant_name_id","area_code_l3","introduced","extinct","location_doubtful"]].copy()
    data["plant_name_id"]=data.plant_name_id.astype(str)
    data["area_code_l3"]=data.area_code_l3.astype(str)
    data=data[data.area_code_l3.isin(strict_codes)]
    for col in ("introduced","extinct","location_doubtful"):
        data[col]=pd.to_numeric(data[col],errors="coerce")
    data=data[(data.extinct==0)&(data.location_doubtful==0)&data.introduced.isin([0,1])]
    data=data.merge(names,on="plant_name_id",how="inner",validate="many_to_one")
    data=data[["plant_name_id","area_code_l3","introduced","taxon_key"]].drop_duplicates()
    if data.duplicated(["plant_name_id","area_code_l3"]).any():
        raise ValueError("Conflicting WCVP region native versus introduced status")
    if len(data)<20000:raise ValueError("Not enough WCVP plant-island records")
    d=pd.DataFrame([{"taxon_key":k,**v} for k,v in trait.items()])
    data=data.merge(d,on="taxon_key",how="left",validate="many_to_one")
    return data

def summarize(data,geography,effort,*,coverage_threshold=.2,min_scored=20):
    regions=[]
    ef={r["code"]:r for r in effort}
    for r in geography:
        code=r["code"]
        d=data[data.area_code_l3==code]
        stats={"code":code,"level1":r["level1"],"region_name":r["region_name"],
               "log_area":float(r["log_area"]),"log2ports":float(r["log2ports"]),
               "ports":int(r["ports"]),
               "source_introduced_species":int(r["introduced"]),
               "source_native_species":int(r["native"]),
               "source_added_resource_units":int(r["contemporary_butterfly_resource_units"])}
        for origin,label in ((1,"introduced"),(0,"native")):
            z=d[d.introduced==origin]
            n=len(z)
            known=z.lifecycle.isin(["annual","biennial","perennial"])
            known_herb=known & z.growth_form.isin(["herb","tree","shrub"])
            annual=z.lifecycle=="annual"
            stats.update({
                f"{label}_N":n,
                f"{label}_lifecycle_known":int(known.sum()),
                f"{label}_lifeform_both_known":int(known_herb.sum()),
                f"{label}_annual_count":int(annual.sum()),
                f"{label}_annual_fraction_known":float(annual.sum()/known.sum()) if known.sum() else None,
                f"{label}_annual_herb_count":int((annual&z.growth_form.eq("herb")&known_herb).sum()),
                f"{label}_annual_herb_fraction_both_known":(
                    float((annual&z.growth_form.eq("herb")&known_herb).sum()/known_herb.sum())
                    if known_herb.sum() else None),
                f"{label}_lifecycle_coverage":float(known.sum()/n) if n else None,
                f"{label}_both_trait_coverage":float(known_herb.sum()/n) if n else None
            })
        if stats["introduced_N"]!=int(r["introduced"]) or stats["native_N"]!=int(r["native"]):
            raise ValueError("Complete pinned WCVP source island counts drift for "+code)
        if code in ef:
            stats["GBIF_allplant_indexed"]=float(ef[code]["gbif_all_plant_indexed_records"])
        else:stats["GBIF_allplant_indexed"]=None
        if stats["GBIF_allplant_indexed"] is not None:
            for metric in ("lifecycle","both"):
                thresh=f"{metric}_coverage" if metric=="lifecycle" else "both_trait_coverage"
                key1= f"introduced_{thresh}" if metric!="lifecycle" else "introduced_lifecycle_coverage"
                key2= f"native_{thresh}" if metric!="lifecycle" else "native_lifecycle_coverage"
                n1=f"introduced_{metric}_known" if metric=="lifecycle" else "introduced_lifeform_both_known"
                n2=f"native_{metric}_known" if metric=="lifecycle" else "native_lifeform_both_known"
                stats[f"eligible_{metric}"]=(
                    stats[n1]>=min_scored and stats[n2]>=min_scored
                    and (stats[key1] or 0)>=coverage_threshold
                    and (stats[key2] or 0)>=coverage_threshold)
        regions.append(stats)
    return regions

def fwl_coef(df,response,controls):
    d=df.copy()
    d["log_gbif"]=np.log1p(d.GBIF_allplant_indexed)
    d["log_intro"]=np.log1p(d.introduced_N)
    d["log_native"]=np.log1p(d.native_N)
    X=[np.ones(len(d)),d.log_area.to_numpy(),d.log_gbif.to_numpy(),d.log_intro.to_numpy(),d.log_native.to_numpy()]
    for key in controls:
        X.append(d[key].to_numpy())
    for lev in sorted(d.level1.unique())[1:]:
        X.append((d.level1==lev).to_numpy().astype(float))
    Z=np.column_stack(X).astype(float)
    if np.linalg.matrix_rank(Z)!=Z.shape[1]:
        raise ValueError("Rank-deficient trait covariate design")
    x=d.log2ports.to_numpy().astype(float)
    y=d[response].to_numpy().astype(float)
    # FWL residualisation uses stable least squares.
    def residual(v):
        return v-Z@np.linalg.lstsq(Z,v,rcond=None)[0]
    yr=residual(y)
    def beta(xarr):
        xr=residual(xarr)
        den=float(xr@xr)
        return None if den<1e-12 else float((xr@yr)/den)
    return beta,x,Z

def port_perm(frame,response,cover_cols,B=B,seed=SEED):
    d=frame.copy()
    n=len(d)
    if n<45 or d.level1.nunique()<5:
        return {"status":"INSUFFICIENT_EXCHANGEABLE_GLOBAL_TRAIT_SUPPORT",
                "n":n,"broad_regions":int(d.level1.nunique())}
    beta,x,_=fwl_coef(d,response,cover_cols)
    observed=beta(x)
    grouped=collections.defaultdict(list)
    med=d.groupby("level1").log_area.median().to_dict()
    for i,row in enumerate(d.itertuples()):
        grouped[str(row.level1)+"|"+str(int(row.log_area>med[row.level1]))].append(i)
    move=[v for v in grouped.values() if len(v)>1 and len(set(float(x[i]) for i in v))>1]
    if not move:return {"status":"NO_EXCHANGEABLE_PORTS","n":n}
    rng=np.random.default_rng(seed)
    null=[]
    for _ in range(B):
        xp=x.copy()
        for ids in move:
            ind=np.asarray(ids,dtype=int)
            xp[ind]=x[rng.permutation(ind)]
        z=beta(xp)
        if z is not None:null.append(z)
    if len(null)<B*.95:return {"status":"INVALID_PERMUTATIONS","n":n}
    mu=float(np.mean(null))
    p=(1+sum(abs(v-mu)>=abs(observed-mu)-1e-12 for v in null))/(len(null)+1)
    return {"status":"EXPLORATORY_REGION_LEVEL_TRAIT_ASSOCIATION",
            "n_island_botanical_regions":n,
            "broad_biogeographic_regions":int(d.level1.nunique()),
            "response":response,"observed_beta_per_log2_indexed_port":observed,
            "null_mean":mu,
            "null_q025":float(np.quantile(null,.025)),
            "null_q975":float(np.quantile(null,.975)),
            "two_sided_descriptive_permutation_p":p,
            "positive_direction":observed>0,
            "exchangeable_regions":len(set(i for inds in move for i in inds)),
            "permutations":len(null),
            "predictor_controls":[
                "ln_WGSRPD3_region_area","ln1p_GBIF_all_plant_indexing",
                "ln1p_total_introduced_flora","ln1p_native_flora",
                "trait_coverage_for_both_origin_groups",
                "WGSRPD_Level1_fixed_effects"],
            "same_species_and_genus_trait_missingness_not_random":True}

def analyse(regions):
    d=pd.DataFrame(regions)
    summary={"n_regions":len(regions),"source_regions_with_botanical_effort":int(d.GBIF_allplant_indexed.notna().sum())}
    out={}
    for typ,suffix in (("lifecycle","annual_fraction_known"),
                       ("both","annual_herb_fraction_both_known")):
        eligible=d[d[f"eligible_{typ}"].fillna(False)].copy()
        primary="introduced_"+suffix
        native="native_"+suffix
        cover=["introduced_lifecycle_coverage","native_lifecycle_coverage"] if typ=="lifecycle" else [
            "introduced_both_trait_coverage","native_both_trait_coverage"]
        eligible["intro_minus_native"]=eligible[primary]-eligible[native]
        out[typ]={
            "cohort_n":len(eligible),
            "coverage_intro_median":float(eligible[cover[0]].median()) if len(eligible) else None,
            "coverage_native_median":float(eligible[cover[1]].median()) if len(eligible) else None,
            "fraction_intro_annual_mean":float(eligible[primary].mean()) if len(eligible) else None,
            "fraction_native_annual_mean":float(eligible[native].mean()) if len(eligible) else None,
            "introduced_trait_selection":port_perm(eligible,primary,cover,seed=SEED),
            "native_geographic_comparator":port_perm(eligible,native,cover,seed=SEED+101),
            "introduced_minus_native_selector":port_perm(eligible,"intro_minus_native",cover,seed=SEED+202)
        }
    summary["models"]=out
    return summary

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ("wcvp-names","wcvp-distributions","q1-audit","island-flora","botanical-effort","out"):
        ap.add_argument("--"+name,required=True,type=Path)
    ap.add_argument("--map-out",type=Path)
    a=ap.parse_args()
    meta=json.loads(a.q1_audit.read_text())
    strict=set(meta["classification"]["strict_island_codes"])
    if len(strict)!=126:raise ValueError("Pinned 126 source island names drift")
    with a.island_flora.open(encoding="utf-8",newline="") as f:flora=list(csv.DictReader(f))
    with a.botanical_effort.open(encoding="utf-8",newline="") as f:effort=list(csv.DictReader(f))
    if (len(flora),len(effort))!=(126,98):raise ValueError("Frozen WCVP/GBIF island panel changed")
    print("LOADING_PINNED_WCVP_RDA",flush=True)
    names,distributions=read_wcvp(a.wcvp_names,a.wcvp_distributions)
    print("PINNED_WCVP_ACCEPTED_BINOMIALS",len(names),
          "DISTRIBUTION_ROWS",len(distributions),flush=True)
    species,receipt=all_pages({"query":"species"},100000,6,"species")
    traits={}
    for label,traitid in SOURCE.items():
        rows,pages=all_pages({"query":"traits","traitid":traitid,"biasref":0,"biasderiv":0},
                           10000,35,label)
        traits[label]=rows
        receipt.extend(pages)
    print("GIFT_3p2_PAGES_FETCHED",len(receipt),flush=True)
    byname,counts,identity=traits_from_gift(species,traits)
    if len(byname)<80000:raise ValueError("Unexpectedly small GIFT matched life-cycle catalogue")
    data=source_table(names,distributions,byname,strict)
    regions=summarize(data,flora,effort,coverage_threshold=.2,min_scored=20)
    result={"schema":"chocho_WCVP_v13_GIFT3p2_lifecycle_annual_herb_island_port_trait_selection_v1",
        "status":"MEASURED_OBSERVATIONAL_IF_COMPLETE_COVERAGE_GATES",
        "trait_source_version":"GIFT 3.2 2026 official API; types Lifecycle_1 and Growth_form_1",
        "WCVP_commit":WCVP_COMMIT,
        "WCVP_ID_grain":"accepted species plant_name_id × WGSRPD3 region, introduced 1/native 0",
        "GIFT_trait_coverage_is_not_complete":True,
        "trait_definition":{"annual":"GIFT 2.1.1 annual among unambiguous annual/perennial/biennial",
                            "annual_herb":"GIFT annual AND GIFT 1.2.1 herb, only where both classified",
                            "ambiguity":"variable, mixed, low agreement <0.66 and unresolvable species excluded",
                            "not_equivalent_to":"ruderal, nitrophilous, fast-growing physiology or plant performance"},
        "source":{"GIFT_pages":receipt,"WCVP_names_sha256":hashlib.sha256(a.wcvp_names.read_bytes()).hexdigest(),
                  "WCVP_distribution_sha256":hashlib.sha256(a.wcvp_distributions.read_bytes()).hexdigest(),
                  "GIFT_value_counts":counts,"name_key_audit":identity},
        "source_island_species_rows":len(data),
        "all_introduced_WCVP_species_island_rows":int((data.introduced==1).sum()),
        "all_native_WCVP_species_island_rows":int((data.introduced==0).sum()),
        "WCVP_island_trait_status":regions,
        "model_results":analyse(regions),
        "limitations":[
            "Not all WCVP accepted taxa have GIFT life-cycle or growth-form values; non-random trait coverage may drive the comparison.",
            "Annual herb is a limited life-history and morphological proxy, NOT direct nutrient affinity, weediness or ruderal occupancy.",
            "GIFT trait data may reflect flora-based observation bias; native/introduced trait coverage differs across island and family.",
            "Contemporary WCVP introduced status and 2023 WPI ports lack local timing, movement or plant establishment information.",
            "Counterfactual port-exposure shuffles control Level1 x area only, not all socioeconomics, historical shipping or island climate.",
            "Total introduced flora is a potential downstream pathway; conditioning on it tests composition, not causal mediation.",
            "Trait selection at WGSRPD Level3 region need not imply local enrichment near a physical port.",
            "The 239-butterfly larval host panel used downstream remains a separate biological mapping; no butterfly response analyzed here."
        ]}
    if a.map_out:
        tax=names[["plant_name_id","taxon_key"]].drop_duplicates()
        tax=tax[tax.taxon_key.isin(byname)]
        t=pd.DataFrame([{"taxon_key":k,**v} for k,v in byname.items()])
        joined=tax.merge(t,on="taxon_key",how="inner",validate="many_to_one")
        a.map_out.parent.mkdir(parents=True,exist_ok=True)
        joined.to_csv(a.map_out,index=False)
        result["crosswalk_accepted_WCVP_species_with_any_trait"]=len(joined)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    compact=result["model_results"]
    print("GIFT_WCVP_PORT_TRAIT_GLOBAL_RESULTS",
          json.dumps({"source":result["source_island_species_rows"],
                      "name_keys":result["source"]["name_key_audit"],
                      "coverage":{k:{"n":v["cohort_n"],
                           "intro_model":v["introduced_trait_selection"],
                           "native_comparator":v["native_geographic_comparator"],
                           "intro_minus_native":v["introduced_minus_native_selector"]}
                          for k,v in compact["models"].items()}},
                     ensure_ascii=False)[:24000],flush=True)

if __name__=="__main__":main()
