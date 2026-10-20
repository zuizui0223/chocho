#!/usr/bin/env python3
"""Freeze pre-outcome alternative real plants for glacial patch graph comparison.

The alternative plant taxon identities/groups were matched using 0 BP climate
alone in TTF v0.5, BEFORE their 21ka climate analogue was read. Here only
WCVP-native plant model availability is used. No genetic response or new
LGM patch-connectivity outcome is accessed or used for source qualification.
"""
from __future__ import annotations
import argparse,csv,hashlib,json
from collections import Counter
from pathlib import Path

SOURCE_SHA={
    "original_0BP_alternative_plan":"c8099e516a42850c8ab946c4cdc23e71afe1e8eb5360a90fcd848365b3cb8075",
    "all518_native_WCVP_models":"153189ad188a8647c49bb5401b20786098386979ca1f00f1cbc67bf55a0ad733"
}

def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def eligible_actual_and_other_plants(original,modelled,min_groups=4):
    if min_groups!=4:raise ValueError("minimum already frozen as four modeled alternative groups")
    result={}
    statuses=[]
    for species,entry in sorted(original.items()):
        actual=list(entry["known_host_taxa"])
        valid=[g for g in entry["alternative_groups"]
               if all(key in modelled for key in g["donor_taxa"])]
        host_ok=all(key in modelled for key in actual)
        passed=host_ok and len(valid)>=min_groups
        if not host_ok:status="SOURCE_ACTUAL_KNOWN_HOST_MODEL_INCOMPLETE"
        elif len(valid)<min_groups:status="SOURCE_TOO_FEW_MODELED_ALTERNATIVE_PLANT_GROUPS"
        else:status="SOURCE_QUALIFIED_BEFORE_NEW_LGM_CONTACT_RESULT"
        statuses.append({"insect_species":species,"status":status,
            "actual_host_n":len(actual),"actual_host_models_complete":host_ok,
            "original_alternative_plant_groups":len(entry["alternative_groups"]),
            "modeled_alternative_plant_groups":len(valid)})
        if passed:
            alt_groups=[{"group_index":int(g["group_index"]),
                         "donor_taxa":list(g["donor_taxa"]),
                         "modern_climate_moment_gap":float(g["mean_plant_signature_gap"])}
                        for g in valid]
            if len(set(g["group_index"] for g in alt_groups))!=len(alt_groups):
                raise ValueError("duplicate preselected alternative-plant group index")
            if any(set(g["donor_taxa"])&set(actual) for g in alt_groups):
                raise ValueError("original negative control contains a recorded host taxon")
            result[species]={"known_host_taxa":actual,"alternative_groups":alt_groups,
                             "previous_modern_matching_scenario":"ORIGINAL_TTF_v05_0BP_ONLY"}
    return result,statuses

def main():
    p=argparse.ArgumentParser()
    for name in ("source-alt-plan","all-518-plant-models","out"):
        p.add_argument("--"+name,type=Path,required=True)
    args=p.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    for key,path in (("original_0BP_alternative_plan",args.source_alt_plan),
                      ("all518_native_WCVP_models",args.all_518_plant_models)):
        if digest(path)!=SOURCE_SHA[key]:
            raise ValueError("source exact 0BP-only plant selection/model result digest drift")
    plan=json.loads(args.source_alt_plan.read_text(encoding="utf8"))
    if (plan["stage"]!="CURRENT_0BP_SELECTION_ONLY_NO_LGM_CLIMATE_OPENED"
        or not plan["no_lgm_raster_or_genetic_data_used"]
        or plan["eligible_alternative_control_n"]!=335):
        raise ValueError("original TTF modern-only alternative plant source stage changed")
    with args.all_518_plant_models.open(newline="",encoding="utf8") as f:
        rows=list(csv.DictReader(f))
    if len(rows)!=518:raise ValueError("original WCVP plant model panel changed")
    modelled={r["wcvp_id"] for r in rows
              if r["model_status"]=="COMPLETE_POTENTIAL_CLIMATE_HINDCAST"}
    if len(modelled)!=416:raise ValueError("416 source plants changed")
    retained,status=eligible_actual_and_other_plants(plan["matched_species"],modelled,4)
    if len(status)!=335 or len(retained)!=292:raise ValueError("original pre-LGM 335-to292 eligibility changed")
    if sum(x["status"]=="SOURCE_ACTUAL_KNOWN_HOST_MODEL_INCOMPLETE" for x in status)!=41:
        raise ValueError("actual modern native host model completeness changed")
    counts=[len(v["alternative_groups"]) for v in retained.values()]
    donors={x for v in retained.values() for group in v["alternative_groups"] for x in group["donor_taxa"]}
    if sum(counts)!=2108 or len(donors)!=362 or sum(x==8 for x in counts)!=146:
        raise ValueError("fully-modelled donor group counts changed")
    source={"schema":"chocho_v08_native_model_prequalified_0BP_matched_real_plant_sets",
        "status":"ELIGIBILITY_ONLY_NO_NEW_LGM_PATCH_OUTCOME_READ",
        "original_0BP_climate_matched_insects":335,
        "modelled_actual_host_insects":294,
        "complete_known_host_and_at_least_4_group_insects":292,
        "valid_frozen_alternative_group_count":2108,
        "distinct_other_plant_model_taxa":len(donors),
        "all_eight_alternatives_modelled_insects":146,
        "original_TTF_and_chocho_source_sha256":dict(SOURCE_SHA),
        "no_new_21ka_climate_or_genetic_patch_response_accessed":True,
        "matched_species":retained,
        "limits":["Source alternative plant identities selected on 0 BP observational climate moments in prior TTF v0.5",
           "No outcome-based species/group substitutions; missing all-recorded-host or <4 group plants are not backfilled",
           "Not-listed larval host plant is not a biologically proven non-host",
           "Current WCVP native models are modern-trained potential niches, never 21ka occupied ranges"]}
    (args.out/"source_fixed_292_otherplant_groups.json").write_text(
          json.dumps(source,indent=2,sort_keys=True)+"\n")
    with (args.out/"all_335_preoutcome_native_plant_model_eligibility.csv").open("w",newline="",encoding="utf8") as f:
        w=csv.DictWriter(f,fieldnames=list(status[0]),lineterminator="\n");w.writeheader();w.writerows(status)
    summary={k:v for k,v in source.items() if k!="matched_species"}
    (args.out/"summary.json").write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
    print(json.dumps(summary,sort_keys=True))

if __name__=="__main__":main()
