#!/usr/bin/env python3
"""Merge independently sourced Canary pages and earlier targeted endemic rescue without doubles.

This reports the number of externally independently enumerated butterfly species
with *some named larval host genus in reviewed sources*. NOT an unbiased
estimate of the complete biological larval-host knowledge or resource ecology.
"""
import argparse
import csv
import json
from pathlib import Path

ROOT=Path("provenance/reviewer_defenses/results")
BASE=ROOT/"independent_island_butterfly_host_coverage_20261117"
NEXT=ROOT/"island_endemic_synonym_origin_20261118"
COLLISIONS={
    "Hipparchia bacchus","Hipparchia gomera","Pararge xiphia",
    "Pararge xiphioides","Thymelicus christi",
}

def rows(path):
    with path.open(encoding="utf-8-sig",newline="") as f:
        return list(csv.DictReader(f))

def audit(uniform_csv, output_dir):
    frozen=rows(BASE/"independent_42_species_vs_frozen239.csv")
    rescue=rows(NEXT/"independent_larval_host_rescue_10.csv")
    verify=json.loads((NEXT/"VERIFIED_SOURCE_OUTCOME_20261118.json").read_text())
    alias=verify["all_unique_42"]["newly_confirmed_alias_matches"]
    source=rows(uniform_csv)
    biota=rows(NEXT/"canary_BIOTA_source_status_20.csv")+rows(NEXT/"canary_BIOTA_remaining_13.csv")
    if (len(frozen),len(rescue),len(source),len(biota))!=(42,10,33,33):
        raise ValueError("Pinned source populations or uniform Canary 33 source changed")
    if len(alias)!=2 or set(alias)!={"Argynnis pandora","Pontia daplidice"}:
        raise ValueError("Source-verified two-name alias list drift")
    original={r["scientific_name"]:r for r in frozen}
    origin={r["independent_checklist_name"]:r["canary_origin_stratum"] for r in biota}
    by_source={r["independent_checklist_species"]:r for r in source}
    from_literature={r["independent_checklist_species"] for r in rescue}
    if len(set(origin))!=33 or len(by_source)!=33 or set(origin)!=set(by_source):
        raise ValueError("Canary uniform source and origin denominator differ")
    if not set(alias)<=set(original) or not from_literature<=set(original):
        raise ValueError("A rescue name is absent from independent 42")
    if set(alias)&from_literature:
        raise ValueError("Taxonomic additional HOSTS aliases overlap source rescue")
    combined=[]
    for species,original_row in sorted(original.items()):
        exact=int(original_row["hosts_named_genus"])
        alias_verified=int(species in alias)
        targeted=int(species in from_literature)
        source_row=by_source.get(species)
        uniform=int(source_row["host_genus_documented_uniform"]) if source_row else 0
        if targeted and (exact or alias_verified):
            raise ValueError("Independently sourced host species was already HOSTS reported: "+species)
        if alias_verified and exact:
            raise ValueError("Alias list includes exact original HOSTS species: "+species)
        all_sources=int(any((exact,alias_verified,targeted,uniform)))
        combined.append({
            "independent_checklist_species":species,
            "family":original_row["checklist_family"],
            "archipelagos":original_row["archipelagos"],
            "canary_independent_origin":origin.get(species,"not_in_Canary"),
            "madeira_origin":original_row["madeira_origin"],
            "exact_HOSTS_named_host_genus":exact,
            "guarded_source_confirmed_alias_only":alias_verified,
            "targeted_endemic_external_larval_host_genus":targeted,
            "uniform_Canary_source_larval_host_genus":uniform,
            "first_documented_by_uniform_beyond_prior_sources":int(
                uniform and not (exact or alias_verified or targeted)
            ),
            "multi_source_named_host_genus_confirmed":all_sources,
            "original_GEB239_exact_selected":int(original_row["in_frozen239"]),
        })
    if sum(r["exact_HOSTS_named_host_genus"] for r in combined)!=19:
        raise ValueError("Exact source baseline changed")
    if sum(r["guarded_source_confirmed_alias_only"] for r in combined)!=2:
        raise ValueError("Synonym confirmation baseline changed")
    if sum(r["targeted_endemic_external_larval_host_genus"] for r in combined)!=10:
        raise ValueError("Targeted independent literature rescue baseline changed")
    if sum(r["uniform_Canary_source_larval_host_genus"] for r in combined)!=24:
        raise ValueError("Source-balance study baseline changed")
    additions=[r["independent_checklist_species"] for r in combined
               if r["first_documented_by_uniform_beyond_prior_sources"]]
    if sorted(additions)!=sorted([
        "Aricia cramera","Euchloe charlonia",
        "Polyommatus celina","Vanessa vulcania",
    ]):
        raise ValueError("Unexpected genuinely new named larval-host reports: "+str(additions))
    count=sum(r["multi_source_named_host_genus_confirmed"] for r in combined)
    if count!=35:
        raise ValueError("Source union must recover 35/42 not "+str(count))
    def summary(sub):
        return {
            "n":len(sub),
            "source_documented_genus_n":sum(r["multi_source_named_host_genus_confirmed"] for r in sub),
            "source_undocumented_n":sum(not r["multi_source_named_host_genus_confirmed"] for r in sub),
            "original_exact_HOSTS_n":sum(r["exact_HOSTS_named_host_genus"] for r in sub),
            "prior_targeted_literature_n":sum(r["targeted_endemic_external_larval_host_genus"] for r in sub),
            "uniform_source_genus_n":sum(r["uniform_Canary_source_larval_host_genus"] for r in sub),
            "original_GEB239_exact_selection_n":sum(r["original_GEB239_exact_selected"] for r in sub),
        }
    canary=[r for r in combined if r["canary_independent_origin"]!="not_in_Canary"]
    endemic=[r for r in canary if r["canary_independent_origin"]=="endemic_confirmed"]
    native=[r for r in canary if r["canary_independent_origin"].startswith("native_nonendemic")]
    introduced=[r for r in canary if r["canary_independent_origin"]=="introduced_probable"]
    madeira_endemic=[r for r in combined if r["madeira_origin"]=="island_endemic_species"]
    extinct=[r for r in combined if r["madeira_origin"]=="historic_extinct_endemic_taxon"]
    azores_endemic=[r for r in combined if r["independent_checklist_species"] in {
        "Hipparchia azorina","Hipparchia miguelensis"
    }]
    unresolved=[r["independent_checklist_species"] for r in combined
                if not r["multi_source_named_host_genus_confirmed"]]
    if (len(endemic),len(native),len(introduced),len(extinct))!=(13,19,1,1):
        raise ValueError("Independent origin strata invalid")
    if sum(r["multi_source_named_host_genus_confirmed"] for r in native)!=19:
        raise ValueError("Uniform source plus HOSTS should resolve all 19 native nonendemics")
    if sum(r["multi_source_named_host_genus_confirmed"] for r in endemic)!=7:
        raise ValueError("Independent nonuniform rescue should resolve 7/13 endemics")
    if sum(r["multi_source_named_host_genus_confirmed"] for r in madeira_endemic)!=3:
        raise ValueError("Madeira known larval host genus external sources absent")
    if sum(r["multi_source_named_host_genus_confirmed"] for r in azores_endemic)!=2:
        raise ValueError("Azores known HOSTS endemic host genus absent")
    result={
        "schema":"chocho_42_source_rescue_union_consistent_denominator_v1",
        "independent_species_n":42,
        "original_exact_NHM_HOSTS_n":19,
        "guarded_taxonomy_NHM_HOSTS_min":21,
        "targeted_endemic_external_larval_hosts_additional_n":10,
        "uniform_source_new_larval_hosts_additional_n":4,
        "uniform_source_new_names":additions,
        "combined_documented_host_genus_lower_bound_n":count,
        "combined_documented_host_genus_lower_bound_pct":round(count/42*100,4),
        "exclude_historic_extinct_n":41,
        "exclude_historic_extinct_documented_n":count,
        "undocumented_named_host_genus_species":unresolved,
        "cohorts":{
            "all42":summary(combined),
            "Canary_all33":summary(canary),
            "Canary_endemic13":summary(endemic),
            "Canary_native_nonendemic19":summary(native),
            "Canary_introduced_probable1":summary(introduced),
            "Madeira_extant_endemic3":summary(madeira_endemic),
            "Azores_Hipparchia_endemic2":summary(azores_endemic),
        },
        "warnings":[
            "35/42 is a lower bound on host-genus documentation found in specific sources, NOT an unbiased biological true host-knowledge rate.",
            "The target-rich external literature supplied ten endemics nonuniformly, while the later single Polinizadores source queried all 33 Canarian butterflies equally.",
            "The six remaining nonextinct unreported taxa are endemic (all Canary) including a related Hipparchia species complex: taxonomic/phylogenetic sources can still be selective.",
            "The natural history larval-host evidence is genus-level and not species-level; native-versus-introduced local plant availability and performance cannot be inferred.",
            "The original GEB239 resource expansion +54.9%/+118.9% remains conditional to species-level mapped host taxa and has not been recomputed from the independent genus-only records.",
            "Original GEB239 selection count for this independent cohort requires separate taxonomically corrected crosswalk; exact-binomial count is not the same as matched species identity."
        ]
    }
    output_dir.mkdir(parents=True,exist_ok=True)
    with (output_dir/"independent42_multisource_larval_genus_union.csv").open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(combined[0]));w.writeheader();w.writerows(combined)
    (output_dir/"independent42_multisource_larval_genus_union.json").write_text(
        json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8"
    )
    print("INDEPENDENT42_MULTISOURCE_HOST_GENUS_UNION_PASSED")
    print("INDEPENDENT42_MULTISOURCE_HOST_GENUS_UNION_SUMMARY="+json.dumps({
        "n":42,"combined":count,"new":additions,"unreported":unresolved,
        "cohorts":result["cohorts"]
    },ensure_ascii=False))
    return result

if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--uniform-validated-csv",required=True,type=Path)
    ap.add_argument("--out-dir",required=True,type=Path)
    args=ap.parse_args()
    audit(args.uniform_validated_csv,args.out_dir)
