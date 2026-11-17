#!/usr/bin/env python3
"""Validate uniform 33-species Canary larval-host source and quantify ascertainment.

Manual 24-species host genus annotations are accepted ONLY when each genus
appears verbatim within that source's larval CICLO VITAL section and the page
was independently found for that butterfly. Exact fixed-margin source-coverage
tests are descriptive of document availability, never of ecological host lack.
"""
import argparse
import csv
import json
import math
import re
from collections import defaultdict
from pathlib import Path

ANNOTATIONS = Path("provenance/reviewer_defenses/results/island_endemic_synonym_origin_20261118/uniform_source_larval_genera_24.json")

def read_rows(path):
    with path.open(encoding="utf-8-sig",newline="") as f:
        return list(csv.DictReader(f))

def hypergeom_tail_low(endo,other,total_present,observed_endo_present):
    low=max(0,total_present-other)
    high=min(endo,total_present)
    denom=math.comb(endo+other,total_present)
    dist={k:math.comb(endo,k)*math.comb(other,total_present-k)/denom
          for k in range(low,high+1)}
    obs_p=dist[observed_endo_present]
    return {
        "one_sided_endemic_less_or_equal":sum(p for k,p in dist.items() if k<=observed_endo_present),
        "two_sided_fisher_probability_ordered":sum(p for p in dist.values() if p<=obs_p+1e-15),
        "null_expected_endemic_present":endo*total_present/(endo+other),
    }

def family_conditional_one_sided(by_family):
    """Exact convolution retaining each family's endemic/native counts and page total."""
    distribution={0:1.0}
    details=[]
    for family,parts in sorted(by_family.items()):
        endo=parts["endemic"]
        native=parts["native"]
        e,n=len(endo),len(native)
        if not (e+n):
            continue
        total=sum(r["host_genus_documented_uniform"] for r in endo+native)
        denom=math.comb(e+n,total)
        dist={k:math.comb(e,k)*math.comb(n,total-k)/denom
              for k in range(max(0,total-n),min(e,total)+1)}
        nxt=defaultdict(float)
        for current,p in distribution.items():
            for k,prob in dist.items():
                nxt[current+k]+=p*prob
        distribution=dict(nxt)
        details.append({
            "family":family,"endemic_n":e,"native_nonendemic_n":n,
            "source_pages_endemic":sum(r["host_genus_documented_uniform"] for r in endo),
            "source_pages_native":sum(r["host_genus_documented_uniform"] for r in native),
            "source_pages_total":total,
        })
    obs=sum(x["source_pages_endemic"] for x in details)
    return {
        "observed_endemic_pages":obs,
        "one_sided_conditional_fixed_family_origin_page_totals":sum(
            p for k,p in distribution.items() if k<=obs
        ),
        "family_tables":details,
        "interpretation":"Source account availability association conditional on family, NOT a biological resource-absence effect; species share phylogeny."
    }

def run(source_csv, annotation_json, output_dir):
    rows=read_rows(source_csv)
    raw=json.loads(annotation_json.read_text(encoding="utf-8"))
    manually=raw["annotations"]
    if len(rows)!=33 or len({r["independent_checklist_species"] for r in rows})!=33:
        raise ValueError("Original complete INBO 33-species denominator changed")
    if len(manually)!=24 or sum(len(v) for v in manually.values())!=86:
        raise ValueError("Curated source annotation counts changed")
    names={r["independent_checklist_species"] for r in rows}
    if not set(manually)<=names:
        raise ValueError("Curated source not a subset of independent INBO butterflies")
    for row in rows:
        name=row["independent_checklist_species"]
        documented=int(row["source_ciclo_vital_larva_text_present"])
        approved=manually.get(name,[])
        if documented!=int(bool(approved)):
            raise ValueError("Reported life section and audited host genera disagree: " + name)
        if documented:
            if row["fetch_status"] not in ("http_200","cached") or int(row["source_page_identity_ok"])!=1:
                raise ValueError("No validated page source for " + name)
            section=row["source_life_section_verbatim"]
            for genus in approved:
                if not re.search(r"\b"+re.escape(genus)+r"\b",section,re.I):
                    raise ValueError("Manually asserted larval host genus absent from life section: "+name+": "+genus)
            if len(set(approved))!=len(approved):
                raise ValueError("Repeated manually annotated larval genus")
        row["host_genus_documented_uniform"]=int(bool(approved))
        row["larval_host_genera_uniform"]=";".join(approved)
        row["larval_host_genus_count_uniform"]=len(approved)
        row["union_original_exact_HOSTS_and_uniform_source"]=int(
            int(row["NHM_HOSTS_exact_genus"])==1 or bool(approved)
        )
    if sum(r["host_genus_documented_uniform"] for r in rows)!=24:
        raise ValueError("Expected 24 independently sourced host-section reports")
    if sum(r["union_original_exact_HOSTS_and_uniform_source"] for r in rows)!=26:
        raise ValueError("Expected 26 document-recovered taxa from original exact HOSTS or independent source")
    by_origin=defaultdict(list)
    by_family=defaultdict(lambda:{"endemic":[],"native":[]})
    for row in rows:
        by_origin[row["canary_origin"]].append(row)
        if row["canary_origin"]=="endemic_confirmed":
            by_family[row["family"]]["endemic"].append(row)
        elif row["canary_origin"].startswith("native_nonendemic"):
            by_family[row["family"]]["native"].append(row)
    def report(sub):
        return {
            "n":len(sub),
            "uniform_source_larval_genus":sum(r["host_genus_documented_uniform"] for r in sub),
            "exact_original_HOSTS":sum(int(r["NHM_HOSTS_exact_genus"]) for r in sub),
            "union_original_exact_HOSTS_plus_uniform_larval_genus":sum(
                r["union_original_exact_HOSTS_and_uniform_source"] for r in sub
            ),
            "page_missing_not_biological_host_absence":sum(
                r["fetch_status"] not in ("http_200","cached") for r in sub
            ),
            "larval_host_genus_mentions":sum(r["larval_host_genus_count_uniform"] for r in sub),
        }
    strata={stratum:report(sub) for stratum,sub in sorted(by_origin.items())}
    endemic=by_origin["endemic_confirmed"]
    native=by_origin["native_nonendemic_secure"]+by_origin["native_nonendemic_probable"]
    if (len(endemic),len(native))!=(13,19):
        raise ValueError("Origin comparison denominator drift")
    fixed=hypergeom_tail_low(
        len(endemic),len(native),
        sum(r["host_genus_documented_uniform"] for r in endemic+native),
        sum(r["host_genus_documented_uniform"] for r in endemic)
    )
    fam=family_conditional_one_sided(by_family)
    # A related Canary Hipparchia endemic radiation accounts for five of the
    # seven missing endemic account pages. Do not claim five independent
    # evolutionary replications or an origin effect after removing the clade.
    endemic_without_hipparchia=[
        r for r in endemic if not r["independent_checklist_species"].startswith("Hipparchia ")
    ]
    if len(endemic_without_hipparchia)!=8:
        raise ValueError("Canary endemic Hipparchia species-complex sensitivity count changed")
    page_endo_without=sum(x["host_genus_documented_uniform"] for x in endemic_without_hipparchia)
    if page_endo_without!=6:
        raise ValueError("Unexpected non-Hipparchia endemic page coverage")
    excluding_clade=hypergeom_tail_low(
        len(endemic_without_hipparchia),len(native),
        sum(x["host_genus_documented_uniform"] for x in endemic_without_hipparchia+native),
        page_endo_without
    )
    # Separate potentially disputed regional-native classifications from a
    # source-asserted origin association. This uses the same response status.
    discordant={"Cacyreus marshalli","Pieris rapae"}
    native_conservative=[x for x in native if x["independent_checklist_species"] not in discordant]
    if len(native_conservative)!=17:
        raise ValueError("Regional-origin sensitivity cohort changed")
    excluding_discordant=hypergeom_tail_low(
        len(endemic),len(native_conservative),
        sum(x["host_genus_documented_uniform"] for x in endemic+native_conservative),
        sum(x["host_genus_documented_uniform"] for x in endemic)
    )
    strict_native=by_origin["native_nonendemic_secure"]
    if len(strict_native)!=5:
        raise ValueError("Native-secure control category changed")
    strict_comparison=hypergeom_tail_low(
        len(endemic),len(strict_native),
        sum(x["host_genus_documented_uniform"] for x in endemic+strict_native),
        sum(x["host_genus_documented_uniform"] for x in endemic)
    )
    family_two_sided={}
    for family,part in by_family.items():
        if not part["endemic"] or not part["native"]:
            continue
        x=part["endemic"];y=part["native"]
        family_two_sided[family]=hypergeom_tail_low(
            len(x),len(y),
            sum(z["host_genus_documented_uniform"] for z in x+y),
            sum(z["host_genus_documented_uniform"] for z in x)
        )
    result={
        "schema":"chocho_source_balanced_canary33_larval_host_genus_audit_v1",
        "selected_denominator":"All 33 original INBO-present Canarian butterflies, regardless of HOSTS or origin",
        "uniform_source":"https://polinizadoresdecanarias.org/repositorio/",
        "source_section":"CICLO VITAL, NOT adult PLANTAS VISITADAS",
        "n":len(rows),
        "source_species_pages_with_named_larval_host_genus":24,
        "named_larval_host_species_genus_pairs":86,
        "distinct_larval_host_genera":len({
            genus for genera in manually.values() for genus in genera
        }),
        "combined_original_exact_HOSTS_plus_uniform":26,
        "origin_strata":strata,
        "native_nonendemic_combined":report(native),
        "endemic_combined":report(endemic),
        "source_availability_fisher_fixed_margin":fixed,
        "source_availability_family_conditional":fam,
        "source_availability_sensitivities":{
            "exclude_five_related_endemic_Hipparchia":{
                "endemic_n":len(endemic_without_hipparchia),
                "endemic_larval_documented":page_endo_without,
                "native_nonendemic_n":len(native),
                "native_documented":sum(x["host_genus_documented_uniform"] for x in native),
                "conditional_fisher":excluding_clade,
            },
            "exclude_native_origin_discordant_two":{
                "excluded_species":sorted(discordant),
                "endemic_n":len(endemic),
                "native_nonendemic_n":len(native_conservative),
                "native_documented":sum(x["host_genus_documented_uniform"] for x in native_conservative),
                "conditional_fisher":excluding_discordant,
            },
            "native_secure_only":{
                "native_nonendemic_n":len(strict_native),
                "native_documented":sum(x["host_genus_documented_uniform"] for x in strict_native),
                "conditional_fisher":strict_comparison
            },
            "within_family_two_sided":family_two_sided
        },
        "strict_interpretation":[
            "A source-visible butterfly entry is not a direct observation of larval feeding in the wild.",
            "All 33 were tested for presence of a source species account with a larval foodplant genus.",
            "0/7 nonretrieved endemic pages cannot be interpreted as 0 native host plants.",
            "The literature source itself covers endemic taxa less often, even after conditioning on butterfly family; phylogenetic clustering makes raw p-values optimistic as a general evolutionary replication.",
            "Genus-level host labels cannot be translated to introduced/native WCVP species without finer plant identity and local origin.",
            "In the single uniform source, the 6 source-covered endemics and 18 source-covered native butterflies each have at least one named larval genus. Conditional on record availability this is 6/6 and 18/18, not evidence of comparable true larval host availability.",
            "Do not inflate original GEB host opportunity +54.9%/+118.9% by imputation."
        ]
    }
    output_dir.mkdir(parents=True,exist_ok=True)
    with (output_dir/"uniform_canary33_host_genus_verified.csv").open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)
    (output_dir/"uniform_canary33_host_genus_summary.json").write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print("UNIFORM_CANARY33_NAMED_LARVAL_GENUS_AUDIT_PASSED")
    print("UNIFORM_CANARY33_NAMED_LARVAL_GENUS_SUMMARY="+json.dumps({
        "n":result["n"],"genus_pairs":result["named_larval_host_species_genus_pairs"],
        "distinct_genera":result["distinct_larval_host_genera"],
        "uniform_source":result["source_species_pages_with_named_larval_host_genus"],
        "union_exactHOSTS_plus_source":result["combined_original_exact_HOSTS_plus_uniform"],
        "cohorts":strata,
        "fisher":fixed,
        "family_conditional":fam["one_sided_conditional_fixed_family_origin_page_totals"],
        "without_Hipparchia_fisher_two_sided":excluding_clade["two_sided_fisher_probability_ordered"],
        "without_native_discordant_fisher_two_sided":excluding_discordant["two_sided_fisher_probability_ordered"]
    },ensure_ascii=False))
    return result

if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--uniform-csv",required=True,type=Path)
    ap.add_argument("--annotations",type=Path,default=ANNOTATIONS)
    ap.add_argument("--out-dir",required=True,type=Path)
    a=ap.parse_args()
    run(a.uniform_csv,a.annotations,a.out_dir)
