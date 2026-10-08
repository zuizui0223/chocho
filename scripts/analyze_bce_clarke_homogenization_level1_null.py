#!/usr/bin/env python3
"""BCE/Clarke source-coverage sensitivity under WGSRPD Level1-constrained null.

Replicates the original resource-homogenization Level1 edge-swap sampler.
Compares two source-conditioned networks on the SAME original 355 WGSRPD3
regions; potential resource geography is not actual consumer population biology.
"""
from __future__ import annotations

import argparse
import json
import random
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median

import analyze_butterfly_resource_homogenization as core
import analyze_butterfly_resource_homogenization_level1_null as original
from analyze_bce_clarke_homogenization_fixed_domain import validate_and_construct


def within_level1_null(native_sets, added_sets, native_masks, level1_labels, permutations, seed):
    """Reproduce original within-Level1 swap mechanism, with explicit invariants."""
    if not (len(native_sets) == len(added_sets) and len(native_masks) == len(level1_labels)):
        raise ValueError("Invalid source-matrix dimensions")
    region_count = len(level1_labels)
    if not region_count or len(native_sets) == 0:
        raise ValueError("Empty source matrix")
    native = [set(x) for x in native_sets]
    added = [set(x) for x in added_sets]
    if any(n & a for n,a in zip(native,added)):
        raise ValueError("Native and introduced-only resource cells must be disjoint")
    if any(not (0<=c<region_count) for row in native+added for c in row):
        raise ValueError("Regional cell index outside frozen domain")

    original_region_cols=Counter(c for row in added for c in row)
    original_row_totals=[len(row) for row in added]

    def level1_row_margins(rows):
        out=[]
        for row in rows:
            counts=Counter(level1_labels[c] for c in row)
            out.append(dict(counts))
        return out

    by_level1_original = level1_row_margins(added)
    # Stable traversal is REQUIRED: Python set insertion/iteration order can
    # differ when host strings are read through a hash-dependent source union.
    # A seed alone is insufficient to reproduce the exact null draw sequence.
    edges=[(i,c) for i,row in enumerate(added) for c in sorted(row)]
    groups=defaultdict(list)
    for k,(_,c) in enumerate(edges):
        groups[level1_labels[c]].append(k)
    eligible=sorted(key for key,v in groups.items() if len(v)>=2)
    if not eligible or not edges:
        raise RuntimeError("No original added edges in eligible Level1 groups")

    introduced_masks=[0]*region_count
    for i, cols in enumerate(added):
        for col in cols:
            introduced_masks[col] |= 1<<i

    rng=random.Random(seed)

    def swaps(attempts):
        accepted=0
        for _ in range(attempts):
            group=rng.choice(eligible)
            members=groups[group]
            k1=members[rng.randrange(len(members))]
            k2=members[rng.randrange(len(members))]
            if k1==k2: continue
            i,a=edges[k1]
            j,b=edges[k2]
            if i==j or a==b:continue
            if b in added[i] or a in added[j]:continue
            if b in native[i] or a in native[j]:continue
            # By group construction: a and b must share their continent.
            if level1_labels[a]!=level1_labels[b]:
                raise RuntimeError("Attempted invalid inter-continent swap")
            added[i].remove(a);added[j].remove(b)
            added[i].add(b);added[j].add(a)
            edges[k1]=(i,b);edges[k2]=(j,a)
            bi=1<<i
            bj=1<<j
            introduced_masks[a]^=bi
            introduced_masks[a]^=bj
            introduced_masks[b]^=bj
            introduced_masks[b]^=bi
            accepted+=1
        return accepted

    burn=max(10000,15*len(edges))
    step=max(2000,3*len(edges))
    accepted=swaps(burn)
    regional_samples=[]
    butterfly_samples=[]
    for _ in range(permutations):
        accepted+=swaps(step)
        masks=[native_masks[k]|introduced_masks[k] for k in range(region_count)]
        region_mean,_,_=core.mean_pairwise_jaccard(masks)
        species_mean,_,_=core.mean_pairwise_set_jaccard(
            [native[i]|added[i] for i in range(len(added))]
        )
        regional_samples.append(region_mean)
        butterfly_samples.append(species_mean)
    if level1_row_margins(added)!=by_level1_original:
        raise RuntimeError("Species × Level1 margins not invariant")
    if Counter(c for row in added for c in row)!=original_region_cols:
        raise RuntimeError("Region columns do not preserve added incidence margins")
    if [len(row) for row in added]!=original_row_totals:
        raise RuntimeError("Butterfly rows do not preserve added incidence margins")
    if any(native[i]&added[i] for i in range(len(native))):
        raise RuntimeError("Null places resource additions inside native support")

    return {
        "regional_mean_null_median":median(regional_samples),
        "regional_mean_null_ci95":[core.quantile(regional_samples,0.025),core.quantile(regional_samples,0.975)],
        "butterfly_mean_null_median":median(butterfly_samples),
        "butterfly_mean_null_ci95":[core.quantile(butterfly_samples,0.025),core.quantile(butterfly_samples,0.975)],
        "regional_null_samples":regional_samples,
        "butterfly_null_samples":butterfly_samples,
        "accepted_swaps_total":accepted,
        "added_edges_in_original_region_domain":len(edges),
        "eligible_Level1_groups":len(eligible),
        "burn_attempts":burn,
        "spacing_attempts":step,
        "sampled_margins_verified":True,
    }


def analyze_one(label, species, region_order, native_by_species, current_by_species, level1_by_region, permutations, seed):
    idx={unit:i for i,unit in enumerate(region_order)}
    labels=[level1_by_region[unit] for unit in region_order]
    ns=[{idx[unit] for unit in native_by_species[sp] if unit in idx} for sp in species]
    cs=[{idx[unit] for unit in current_by_species[sp] if unit in idx} for sp in species]
    added=[c-n for n,c in zip(ns,cs)]
    if any(not n<=c for n,c in zip(ns,cs)):
        raise RuntimeError("Native incidence not contained contemporary")
    native_masks=core.masks_from_species_units(species,region_order,native_by_species)
    current_masks=core.masks_from_species_units(species,region_order,current_by_species)
    nreg,_,pn=core.mean_pairwise_jaccard(native_masks)
    creg,_,pc=core.mean_pairwise_jaccard(current_masks)
    nsp,_,sn=core.mean_pairwise_set_jaccard(ns)
    csp,_,sc=core.mean_pairwise_set_jaccard(cs)
    if (pn,pc,sn,sc)!=(62835,62835,28441,28441):
        raise RuntimeError("Original Jaccard denominator mismatch")
    print("START_LEVEL1",label,"native_jaccard",nreg,"current_jaccard",creg,
          "added_edges",sum(map(len,added)),"reps",permutations,flush=True)
    null=within_level1_null(ns,added,native_masks,labels,permutations,seed)
    regional_samples=null.pop("regional_null_samples")
    butterfly_samples=null.pop("butterfly_null_samples")
    null.update({
        "regional_observed_minus_null_median":creg-null["regional_mean_null_median"],
        "regional_one_sided_p":(1+sum(v>=creg for v in regional_samples))/(permutations+1),
        "butterfly_observed_minus_null_median":csp-null["butterfly_mean_null_median"],
        "butterfly_one_sided_p":(1+sum(v>=csp for v in butterfly_samples))/(permutations+1),
        "permutations":permutations,"seed":seed,
        "preserves":["species total added resources","regional total added resources",
                     "species x WGSRPD Level1 added resources","native structural exclusions"],
    })
    print("LEVEL1_RESULT",label,json.dumps({
        "region_native":nreg,"region_current":creg,
        "region_null":null["regional_mean_null_median"],
        "regional_excess":null["regional_observed_minus_null_median"],
        "regional_p":null["regional_one_sided_p"],
        "butterfly_excess":null["butterfly_observed_minus_null_median"],
        "butterfly_p":null["butterfly_one_sided_p"],
        "accepted_swaps_total":null["accepted_swaps_total"]}),flush=True)
    return {
        "scenario":label,"fixed_domain_regions":len(region_order),
        "native_mean_regional_jaccard":nreg,
        "contemporary_mean_regional_jaccard":creg,
        "raw_regional_gain":creg-nreg,
        "native_mean_butterfly_resource_jaccard":nsp,
        "contemporary_mean_butterfly_resource_jaccard":csp,
        "raw_butterfly_resource_gain":csp-nsp,
        "native_incidence_in_fixed_domain":sum(map(len,ns)),
        "contemporary_incidence_in_fixed_domain":sum(map(len,cs)),
        "added_incidence_in_fixed_domain":sum(map(len,added)),
        "level1_fixed_margin_null":null
    }


def main():
    p=argparse.ArgumentParser()
    for name in (
        "protocol_json","source_crosswalk_json","insect_host_csv",
        "native_csv","contemporary_csv","added_native_csv","added_contemporary_csv",
        "frozen_metrics_csv","level3_geojson","output_json"):
        p.add_argument("--"+name.replace("_","-"),type=Path,required=True)
    args=p.parse_args()
    proto=json.loads(args.protocol_json.read_text(encoding="utf-8"))
    if proto.get("schema")!="chocho_bce_clarke_level1_geography_preserving_null_v0.1":
        raise RuntimeError("Unknown frozen protocol")
    # validate_and_construct was independently source-audited before outcome
    class AuditArgs:
        pass
    a=AuditArgs()
    a.crosswalk_json=args.source_crosswalk_json
    a.insect_host_csv=args.insect_host_csv
    a.native_csv=args.native_csv
    a.contemporary_csv=args.contemporary_csv
    a.added_native_csv=args.added_native_csv
    a.added_contemporary_csv=args.added_contemporary_csv
    a.frozen_metrics_csv=args.frozen_metrics_csv
    from pathlib import Path as _Path
    old={
        "fixed_input":{
            "primary_species_count":proto["input"]["focal_panel"],
            "primary_native_butterfly_region_units":26530,
            "primary_contemporary_butterfly_region_units":41083,
            "primary_added_units":proto["input"]["original_added_units_all_regions"],
            "original_active_native_regions":proto["input"]["orig_native_active_regions"],
        },
        "source":{
            "candidate_accepted_link_count":proto["input"]["additional_butterfly_host_links"],
            "candidate_accepted_plant_species":679
        }
    }
    sp, regions, ((oldn,oldc),(newn,newc)), summary=validate_and_construct(a,old)
    level1=original.load_level1(args.level3_geojson)
    missing=sorted(set(regions)-set(level1))
    if missing:
        raise RuntimeError("Level1 map missing original regional codes "+str(missing[:15]))
    seed=proto["design"]["random_seed"]
    reps=proto["design"]["permutations"]
    frozen=analyze_one("frozen_HOSTS",sp,regions,oldn,oldc,level1,reps,seed)
    target=proto["input"]["exact_reference"]
    for actual,refer in (
        (frozen["native_mean_regional_jaccard"],target["original_native_region_jaccard"]),
        (frozen["contemporary_mean_regional_jaccard"],target["original_contemporary_region_jaccard"]),
        (frozen["native_mean_butterfly_resource_jaccard"],target["original_native_butterfly_geography_jaccard"]),
        (frozen["contemporary_mean_butterfly_resource_jaccard"],target["original_contemporary_butterfly_geography_jaccard"]),
    ):
        if abs(actual-refer)>1e-9:
            raise RuntimeError("Original source Jaccard value not reproduced")
    aug=analyze_one("BCE_Clarke_accepted_ID_added",sp,regions,newn,newc,level1,reps,seed)
    result={
        "schema":"chocho_bce_clarke_level1_null_source_sensitivity_result_v0.1",
        "status":"POSTHOC_BCE_HOST_SOURCE_CONTINENT_PRESERVING_NULL_NOT_BIOLOGICAL_CAUSAL_TEST",
        "frozen_source_and_domain_verified":True,
        "source_receipts":{
            "accepted_crosswalk":proto["input"]["BCE_accepted_host_crosswalk_run"],
            "previous_geography":37778961442,  # Fixed provenance receipt, not a statistical input.
            "original_HOSTS_commit":proto["input"]["original_source_Hosts_commit"],
            "original_WCVP_commit":proto["input"]["botany_WCVP_commit"],
            "WGSRPD_commit":proto["input"]["geography_WGSRPD_commit"]
        },
        "shared_region_count":len(regions),
        "base_native_active":summary["baseline_active_native_regions"],
        "augmented_native_active":summary["augmented_active_native_regions"],
        "original":frozen,
        "BCE_augmented":aug,
        "comparison":{
            "regional_excess_augmented_minus_original":(
                aug["level1_fixed_margin_null"]["regional_observed_minus_null_median"]-
                frozen["level1_fixed_margin_null"]["regional_observed_minus_null_median"]),
            "butterfly_excess_augmented_minus_original":(
                aug["level1_fixed_margin_null"]["butterfly_observed_minus_null_median"]-
                frozen["level1_fixed_margin_null"]["butterfly_observed_minus_null_median"])
        },
        "boundaries":[
            "Comparison is on SAME original 355 native-active WGSRPD3 region domain; augmented domain has additional region excluded from contrast.",
            "European BCE/Clarke references are not independent original feeding observations and cannot represent worldwide local host utilization.",
            "Null is conditional on each scenario's margins and species-by-Level1 added opportunity counts.",
            "Structural butterfly resource opportunity, not observed butterfly colonization, competition, or fitness.",
            "Original GEB primary results and PR38 are unchanged."
        ]
    }
    args.output_json.parent.mkdir(parents=True,exist_ok=True)
    args.output_json.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print("LEVEL1_HOST_SOURCE_DECISION",json.dumps({
       "original_excess":frozen["level1_fixed_margin_null"]["regional_observed_minus_null_median"],
       "augmented_excess":aug["level1_fixed_margin_null"]["regional_observed_minus_null_median"],
       "original_p":frozen["level1_fixed_margin_null"]["regional_one_sided_p"],
       "augmented_p":aug["level1_fixed_margin_null"]["regional_one_sided_p"],
       "original_butterfly_excess":frozen["level1_fixed_margin_null"]["butterfly_observed_minus_null_median"],
       "augmented_butterfly_excess":aug["level1_fixed_margin_null"]["butterfly_observed_minus_null_median"]
    }),flush=True)

if __name__=="__main__":
    main()
