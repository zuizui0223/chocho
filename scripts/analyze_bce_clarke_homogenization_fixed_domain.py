#!/usr/bin/env python3
"""Post-hoc source-completeness stress test of the fixed-margin homogenization result.

Reruns the ORIGINAL fixed-margin null separately for frozen HOSTS and the
one-direction BCE/Clarke accepted-host augmentation, on EXACTLY the same
355 originally native-active WGSRPD3 regions. Values are structural potential
resource opportunity, NOT observed butterfly community change or fitness.
"""
from __future__ import annotations
import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

import analyze_butterfly_resource_homogenization as core
from analyze_bce_clarke_added_hosts_geography import load_source, additions, native_current, footprint


def readrows(path: Path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        yield from csv.DictReader(f)


def validate_and_construct(args, protocol):
    src = load_source(args.crosswalk_json)
    if src["verified_butterflies"] != 83:
        raise RuntimeError("Expected 83 evidence-bearing European pages")
    source_additions, nlinks = additions(src)
    if nlinks != protocol["source"]["candidate_accepted_link_count"]:
        raise RuntimeError("Source accepted additional host-link count mismatch")
    distinct = len(set().union(*source_additions.values()))
    if distinct != protocol["source"]["candidate_accepted_plant_species"]:
        raise RuntimeError("Source accepted plant species count mismatch")
    frozen = {r["species"].strip(): r for r in readrows(args.frozen_metrics_csv)}
    if len(frozen) != protocol["fixed_input"]["primary_species_count"]:
        raise RuntimeError("Frozen butterfly sample drift")
    hostids = defaultdict(set)
    for r in readrows(args.insect_host_csv):
        sp = r["insect_species"].strip()
        pid = r["accepted_plant_name_id"].strip()
        if sp and pid: hostids[sp].add(pid)
    orig_n, orig_c = native_current(args.native_csv, args.contemporary_csv)
    added_n, added_c = native_current(args.added_native_csv, args.added_contemporary_csv)
    for pid in set().union(*source_additions.values()):
        if pid in orig_n and orig_n[pid] != added_n[pid]:
            raise RuntimeError("WCVP native range drift for source host "+pid)
        if pid in orig_c and orig_c[pid] != added_c[pid]:
            raise RuntimeError("WCVP current range drift for source host "+pid)
        orig_n[pid] = added_n[pid]
        orig_c[pid] = added_c[pid]
    if not set(source_additions).issubset(frozen):
        raise RuntimeError("BCE host link outside frozen butterfly panel")
    species = sorted(frozen)
    original_native, original_current, augmented_native, augmented_current = {},{},{},{}
    counts = [0, 0, 0]
    for sp in species:
        ids = hostids[sp]
        if not ids:
            raise RuntimeError("No frozen HOSTS edge for butterfly "+sp)
        on = footprint(ids, orig_n)
        oc = footprint(ids, orig_c)
        if not on <= oc:
            raise RuntimeError("Old native geography not contemporary "+sp)
        vals = (len(on), len(oc), len(oc-on))
        want = tuple(int(frozen[sp][k]) for k in
                 ("native_resource_units","contemporary_resource_units","introduced_added_units"))
        if vals != want:
            raise RuntimeError("Full original baseline mismatch "+sp+" "+str(vals)+" != "+str(want))
        nn = footprint(ids | source_additions[sp], orig_n)
        nc = footprint(ids | source_additions[sp], orig_c)
        if not on <= nn or not oc <= nc or not nn <= nc:
            raise RuntimeError("Source-added resource support non-monotonic "+sp)
        original_native[sp], original_current[sp] = on, oc
        augmented_native[sp], augmented_current[sp] = nn, nc
        for i,x in enumerate(vals): counts[i] += x
    expected = protocol["fixed_input"]
    if counts != [expected["primary_native_butterfly_region_units"],
                  expected["primary_contemporary_butterfly_region_units"],
                  expected["primary_added_units"]]:
        raise RuntimeError("Frozen panel aggregate resource count drift "+repr(counts))
    old_active = sorted(set().union(*original_native.values()))
    if len(old_active) != expected["original_active_native_regions"]:
        raise RuntimeError("Original native active WGSRPD3 domain drift")
    new_active = sorted(set().union(*augmented_native.values()))
    summary = {
       "frozen_239_native_units":counts[0],
       "frozen_239_contemporary_units":counts[1],
       "frozen_239_added_units":counts[2],
       "baseline_active_native_regions":len(old_active),
       "augmented_active_native_regions":len(new_active),
       "new_native_active_region_codes":sorted(set(new_active)-set(old_active)),
       "source_additional_accepted_host_links":nlinks,
       "butterflies_with_added_links":len(source_additions),
       "source_unique_accepted_host_species":distinct,
    }
    return species, old_active, (
        (original_native, original_current),
        (augmented_native, augmented_current),
    ), summary


def analyze(scenario, species, region_order, native, contemporary, perm, seed):
    # Always use SAME *original* 355-region analysis domain for both scenarios.
    index = {u: i for i,u in enumerate(region_order)}
    ns = [set(index[r] for r in native[s] if r in index) for s in species]
    cs = [set(index[r] for r in contemporary[s] if r in index) for s in species]
    additions = [c-n for n,c in zip(ns,cs)]
    if any(not n <= c for n,c in zip(ns,cs)):
        raise RuntimeError("scenario native incidence missing from contemporary")
    native_masks = core.masks_from_species_units(species, region_order, native)
    contemp_masks = core.masks_from_species_units(species, region_order, contemporary)
    native_mean,_,npairs = core.mean_pairwise_jaccard(native_masks)
    observed,_,cpairs = core.mean_pairwise_jaccard(contemp_masks)
    if npairs != cpairs or npairs != 62835: raise RuntimeError("regional comparator drift")
    native_sp,_,sp_pairs = core.mean_pairwise_set_jaccard(ns)
    contemporary_sp,_,sp_pairs2 = core.mean_pairwise_set_jaccard(cs)
    if sp_pairs != sp_pairs2 or sp_pairs != 28441:
        raise RuntimeError("butterfly pair comparator drift")
    print("START_NULL",scenario,"fixed_regions",len(region_order),
          "introduced_added_cells_in_domain",sum(len(x) for x in additions),
          "reps",perm,flush=True)
    null = core.swap_null(ns, additions, native_masks, perm, seed)
    reg = null.pop("_null_means")
    spe = null.pop("_null_species_means")
    residual = observed-null["mean_similarity_null_median"]
    residual_sp = contemporary_sp-null["mean_butterfly_niche_overlap_null_median"]
    result = {
      "scenario":scenario,
      "fixed_domain_wgsrpd3_regions":len(region_order),
      "native_resource_incidence_in_fixed_domain":sum(len(x) for x in ns),
      "contemporary_resource_incidence_in_fixed_domain":sum(len(x) for x in cs),
      "introduced_added_in_fixed_domain":sum(map(len,additions)),
      "native_mean_regional_jaccard":native_mean,
      "contemporary_mean_regional_jaccard":observed,
      "regional_raw_gain":observed-native_mean,
      "regional_relative_gain":observed/native_mean-1 if native_mean else None,
      "native_mean_butterfly_resource_overlap":native_sp,
      "contemporary_mean_butterfly_resource_overlap":contemporary_sp,
      "butterfly_resource_overlap_raw_gain":contemporary_sp-native_sp,
      "fixed_margin_null":{
          "permutations":perm,"seed":seed,
          "regional_null_median":null["mean_similarity_null_median"],
          "regional_null_ci95":null["mean_similarity_null_ci95"],
          "regional_observed_minus_null_median":residual,
          "regional_one_sided_p":(1+sum(v>=observed for v in reg))/(len(reg)+1),
          "butterfly_null_median":null["mean_butterfly_niche_overlap_null_median"],
          "butterfly_null_ci95":null["mean_butterfly_niche_overlap_null_ci95"],
          "butterfly_observed_minus_null_median":residual_sp,
          "butterfly_one_sided_p":(1+sum(v>=contemporary_sp for v in spe))/(len(spe)+1),
          "swap_count":null["accepted_swaps_total"],
          "preserves":"Added incidence margins for each butterfly and each region, native forbidden placements",
      }
    }
    print("NULL_RESULT",scenario,json.dumps({
          "regional_n":native_mean,"regional_c":observed,
          "null_med":null["mean_similarity_null_median"],
          "excess":residual,"p":result["fixed_margin_null"]["regional_one_sided_p"],
          "butterfly_excess":residual_sp,
          "butterfly_p":result["fixed_margin_null"]["butterfly_one_sided_p"]
    }),flush=True)
    return result


def main():
    p = argparse.ArgumentParser()
    for arg in ("protocol_json","crosswalk_json","insect_host_csv","native_csv",
                "contemporary_csv","added_native_csv","added_contemporary_csv",
                "frozen_metrics_csv","reference_json","output_json"):
        p.add_argument("--"+arg.replace("_","-"),type=Path,required=True)
    a=p.parse_args()
    contract=json.loads(a.protocol_json.read_text(encoding="utf-8"))
    if contract["schema"]!="chocho_bce_clarke_resource_homogenization_fixed_domain_v0.1":
        raise RuntimeError("Unexpected protocol")
    species,domain,sets,summary=validate_and_construct(a,contract)
    rep=contract["analysis"]["permutations"]
    seed=contract["analysis"]["seed"]
    baseline=analyze("frozen_HOSTS",species,domain,*sets[0],rep,seed)
    expected=contract["fixed_input"]
    for key,name in (
         ("native_mean_regional_jaccard","original_regional_native_jaccard"),
         ("contemporary_mean_regional_jaccard","original_regional_contemporary_jaccard"),
         ("native_mean_butterfly_resource_overlap","original_butterfly_native_jaccard"),
         ("contemporary_mean_butterfly_resource_overlap","original_butterfly_contemporary_jaccard")
    ):
        if abs(baseline[key]-expected[name])>1e-9:
            raise RuntimeError("Original published spatial metric drift "+key+" "+repr(baseline[key]))
    reference=json.loads(a.reference_json.read_text(encoding="utf-8"))
    if reference.get("schema")!="chocho_butterfly_resource_homogenization_v0.1":
        raise RuntimeError("Original fixed margin receipt schema drift")
    if reference["panel"]["active_wgsrpd3_regions"]!=len(domain):
        raise RuntimeError("Original reference domain drift")
    aug=analyze("BCE_Clarke_accepted_ID_source_added",species,domain,*sets[1],rep,seed)
    result={
      "schema":"chocho_bce_clarke_fixed_domain_resource_homogenization_sensitivity_v0.1",
      "status":"POSTHOC_COMPARABLE_MARGINS_NULL_SOURCE_COMPLETENESS_NOT_CORRECTED_BIOLOGY",
      "frozen_baseline_verified":True,
      "panel":summary,
      "fixed_comparable_WGSRPD3_domain_codes":domain,
      "baseline":baseline,
      "BCE_augmented":aug,
      "difference_augmented_minus_original":{
          "native_jaccard":aug["native_mean_regional_jaccard"]-baseline["native_mean_regional_jaccard"],
          "contemporary_jaccard":aug["contemporary_mean_regional_jaccard"]-baseline["contemporary_mean_regional_jaccard"],
          "regional_raw_gain":aug["regional_raw_gain"]-baseline["regional_raw_gain"],
          "regional_excess_over_scenario_matched_fixed_margin_null":(
              aug["fixed_margin_null"]["regional_observed_minus_null_median"]-
              baseline["fixed_margin_null"]["regional_observed_minus_null_median"]),
          "butterfly_resource_overlap_excess_over_scenario_matched_fixed_margin_null":(
              aug["fixed_margin_null"]["butterfly_observed_minus_null_median"]-
              baseline["fixed_margin_null"]["butterfly_observed_minus_null_median"])
      },
      "limits":[
        "BCE/Clarke source compiled in Europe and same host may not be utilized throughout its entire WCVP range.",
        "No host biomass, butterfly range shift, feeding, fitness, competition or demographic effects measured.",
        "Original 355-native-active region domain is kept fixed for numerical comparability; augmented added native regions outside domain are excluded from paired null, but listed separately.",
        "Each null conditions on its OWN resource-incidence margins; difference of p-values is not a test of difference of causal effects.",
        "Posthoc one-direction host-list completion may exacerbate bibliography ascertainment and is not a global corrected outcome.",
        "Original GEB manuscript figures and CI remain frozen."
      ]
    }
    a.output_json.parent.mkdir(parents=True,exist_ok=True)
    a.output_json.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print("PAIRED_HOMOGENIZATION_SENSITIVITY",json.dumps({
      "native_active_original":summary["baseline_active_native_regions"],
      "native_active_augmented":summary["augmented_active_native_regions"],
      "original_regional_excess":baseline["fixed_margin_null"]["regional_observed_minus_null_median"],
      "augmented_regional_excess":aug["fixed_margin_null"]["regional_observed_minus_null_median"],
      "original_p":baseline["fixed_margin_null"]["regional_one_sided_p"],
      "augmented_p":aug["fixed_margin_null"]["regional_one_sided_p"],
    }),flush=True)


if __name__=="__main__":
    main()
