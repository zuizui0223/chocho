#!/usr/bin/env python3
"""Independent Greenstein 2022 127-monarch-host quality x pinned WCVP alien geography.

Protocol: docs/exploratory/MONARCH_HOST_QUALITY_GLOBALIZATION_PROTOCOL_V01.json
Scientific output is a conditional, posthoc plant-distribution association, never fitness.
"""
from __future__ import annotations
import argparse,csv,json,math,random,re
from collections import Counter,defaultdict
from pathlib import Path
from statistics import median

QUALITY_COUNTS={"H1":8,"H2":15,"H3":11,"L1":5,"L2":8,"L3":29,"N":33,"U":18}
def rows(path):
    with path.open(newline="",encoding="utf-8") as f:return list(csv.DictReader(f))
def save_rows(path,records,fields):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(records)
def extract(inp,out):
    tables=json.loads(inp.read_text(encoding="utf-8"))
    assert isinstance(tables,list)
    matching=[x for x in tables if x.get("table_number")==2]
    if len(matching)!=1:raise RuntimeError("Original S1 Appendix Table 2 missing")
    source=matching[0]["rows"]
    if source[0][:2]!=["Plant Species","Class"]:raise RuntimeError("Unrecognized original classification header")
    extracted=[]
    for entry in source[1:]:
        nm=str(entry[0]).strip();detail=str(entry[1]).strip()
        m=re.match(r"^(H[123]|L[123]|N|U)(?=$|\s|\[)",detail)
        if not m:raise RuntimeError(f"Unknown Greenstein source class {nm} {detail}")
        code=m.group(1)
        if len(nm.split())<2:raise RuntimeError(f"Not a species binomial: {nm}")
        extracted.append({"original_species_name":nm,"quality_subclass":code,
            "quality_bin":"high" if code.startswith("H") else ("low" if code.startswith("L") else code),
            "original_source_class_text":detail})
    c=Counter(x["quality_subclass"] for x in extracted)
    if len(extracted)!=127 or dict(c)!=QUALITY_COUNTS:raise RuntimeError(f"127-class drift: {len(extracted)} {dict(c)}")
    if len({x["original_species_name"] for x in extracted})!=127:raise RuntimeError("duplicate original source species")
    save_rows(out,extracted,list(extracted[0]))
    print(json.dumps({"source_species":len(extracted),"subclasses":dict(sorted(c.items()))}),flush=True)

def quantile(xs,p):
    vals=sorted(xs);t=p*(len(vals)-1);i=int(t);j=math.ceil(t)
    return vals[i] if i==j else vals[i]*(j-t)+vals[j]*(t-i)

def analyze_group(rows,seed,perms):
    species=sorted(rows,key=lambda x:(int(x["native_regions"]),x["accepted_id"]))
    labs=[x["quality_bin"] for x in species]
    count=Counter(labs)
    if min(count["high"],count["low"])<5:
        return {"status":"INSUFFICIENT_PERFORMANCE_COVERAGE","n":len(species),"counts":dict(count)}
    X=[math.log1p(int(v["native_regions"])) for v in species]
    Y=[math.log1p(int(v["introduced_regions"])) for v in species]
    mx=sum(X)/len(X);my=sum(Y)/len(Y)
    denom=sum((v-mx)**2 for v in X)
    slope=sum((x-mx)*(y-my) for x,y in zip(X,Y))/denom if denom>0 else 0.0
    residual=[y-(my+slope*(x-mx)) for x,y in zip(X,Y)]
    def contrast(labels):
        hi=[residual[i] for i,x in enumerate(labels) if x=="high"]
        lo=[residual[i] for i,x in enumerate(labels) if x=="low"]
        return sum(hi)/len(hi)-sum(lo)/len(lo)
    observed=contrast(labs)
    n=len(species)
    groups=[list(range(k*n//4,(k+1)*n//4)) for k in range(4)]
    eligible=sum(1 for g in groups if len({labs[i] for i in g})==2)
    if eligible<2:
        return {"status":"INSUFFICIENT_MIXED_NATIVE_RANGE_STRATA","n":n,"counts":dict(count),
                "mixed_strata":eligible}
    rng=random.Random(seed)
    draws=[]
    for _ in range(perms):
        pseudo=labs[:]
        for inds in groups:
            shuffled=[labs[i] for i in inds]
            rng.shuffle(shuffled)
            for j,idx in enumerate(inds):pseudo[idx]=shuffled[j]
        draws.append(contrast(pseudo))
    mid=median(draws)
    p=(1+sum(abs(v-mid)>=abs(observed-mid)-1e-12 for v in draws))/(perms+1)
    introduced=[int(r["introduced_regions"]) for r in species]
    native=[int(r["native_regions"]) for r in species]
    return {
       "status":"COMPUTED_EXPLORATORY_CONDITIONAL_PERMUTATION","n":n,"classes":dict(count),
       "introduced_n_high_median":median(v for v,label in zip(introduced,labs) if label=="high"),
       "introduced_n_low_median":median(v for v,label in zip(introduced,labs) if label=="low"),
       "native_n_high_median":median(v for v,label in zip(native,labs) if label=="high"),
       "native_n_low_median":median(v for v,label in zip(native,labs) if label=="low"),
       "any_introduced_high":sum(n>0 for n,label in zip(introduced,labs) if label=="high"),
       "any_introduced_low":sum(n>0 for n,label in zip(introduced,labs) if label=="low"),
       "effect_residualized_mean_log1p_intro_high_minus_low":observed,
       "null_median":mid,"null_95_interval":[quantile(draws,.025),quantile(draws,.975)],
       "excess_relative_to_conditional_null":observed-mid,
       "two_sided_monte_carlo_p":p,
       "mixed_native_range_quartiles":eligible,"permutations":perms,"seed":seed,
       "interpretation":"Species-level WCVP introduced-range association, not consumer demography, not proof of preferential introduction caused by host quality."
    }

def test(inp,out,perms,seed):
    rec=rows(inp)
    if len(rec)!=127:raise RuntimeError(f"mapped-source rows drift n={len(rec)}")
    by_id=defaultdict(list);match=Counter()
    for row in rec:
        match[row["match_status"]]+=1
        if row["match_status"]=="matched" and row["accepted_id"]:by_id[row["accepted_id"]].append(row)
    accepted=[]
    duplicate_conflicts=[]
    duplicate_consistent=[]
    for ident,group in by_id.items():
        quality={x["quality_bin"] for x in group}
        if len(quality)>1:
            duplicate_conflicts.append([ident,[x["original_species_name"] for x in group]])
        else:
            accepted.append(group[0])
            if len(group)>1:duplicate_consistent.append([ident,[x["original_species_name"] for x in group]])
    def eligible(x):
        return (x["quality_bin"] in ("high","low")
            and int(x["native_regions"])>=1 and int(x["current_regions"])>=1)
    hi_lo=[x for x in accepted if eligible(x)]
    asclepias=[x for x in hi_lo if x["original_species_name"].startswith("Asclepias ")]
    strong=[x for x in asclepias if x["quality_subclass"] in ("H3","L3")]
    sensitivity=[x for x in asclepias if x["original_species_name"]!="Asclepias curassavica"]
    apocynaceae=[x for x in hi_lo if x["family"]=="Apocynaceae"]
    result={
       "schema":"chocho_monarch_host_quality_globalization_v0.1",
       "status":"EXPLORATORY_INDEPENDENT_PERFORMANCE_VS_PLANT_REDISTRIBUTION",
       "protocol":"docs/exploratory/MONARCH_HOST_QUALITY_GLOBALIZATION_PROTOCOL_V01.json",
       "source":"Greenstein et al 2022 S1 Appendix, DOI 10.1371/journal.pone.0269701, independent of WCVP introduced distribution",
       "quality_panel":{"original_species":127,"source_counts":dict(Counter(x["quality_bin"] for x in rec)),
           "matching":dict(match),"unique_matched_accepted_taxa":len(accepted),
           "consistent_synonym_duplicates":duplicate_consistent,"conflicting_duplicates_excluded":duplicate_conflicts,
           "matched_high_low_with_native_range":len(hi_lo),
           "asclepias_high_low_with_native_range":len(asclepias)},
       "primary_Asclepias_H_vs_L":analyze_group(asclepias,seed,perms),
       "secondary_strong_H3_vs_L3_Asclepias":analyze_group(strong,seed,perms),
       "sensitivity_without_Asclepias_curassavica":analyze_group(sensitivity,seed,perms),
       "descriptive_Apocynaceae_H_vs_L":analyze_group(apocynaceae,seed,perms),
       "selected_Asclepias_species":[{"name":x["original_species_name"],"quality":x["quality_subclass"],
         "accepted_name":x["accepted_name"],"native_regions":int(x["native_regions"]),
         "introduced_regions":int(x["introduced_regions"]),"match_type":x["matched_as"]}
          for x in sorted(asclepias,key=lambda r:(-int(r["introduced_regions"]),r["original_species_name"]))],
       "interpretation_boundary":"Evidence classes reflect literature coverage and are not uniform larval experimental survival probabilities. Macro-geographic plant ranges cannot diagnose realized butterfly adaptation or occurrence. Comparisons are exploratory and post-hoc; the original GEB manuscript is unchanged."
    }
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in result.items() if k!="selected_Asclepias_species"},indent=2),flush=True)

def main():
    parser=argparse.ArgumentParser()
    sub=parser.add_subparsers(dest="cmd",required=True)
    p=sub.add_parser("extract");p.add_argument("--original-tables-json",type=Path,required=True)
    p.add_argument("--output-csv",type=Path,required=True)
    p=sub.add_parser("test");p.add_argument("--mapped-csv",type=Path,required=True)
    p.add_argument("--output-json",type=Path,required=True);p.add_argument("--permutations",type=int,default=9999)
    p.add_argument("--seed",type=int,default=20261008)
    a=parser.parse_args()
    if a.cmd=="extract":extract(a.original_tables_json,a.output_csv)
    if a.cmd=="test":test(a.mapped_csv,a.output_json,a.permutations,a.seed)
if __name__=="__main__":main()
