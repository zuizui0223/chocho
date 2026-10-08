#!/usr/bin/env python3
"""Non-confirmatory E. editha HOSTS completeness sensitivity.

Adds EXACTLY the two Haan et al. local hosts absent from the frozen worldwide
species inventory. No new plant identities, no ecological effect size or model
refit, and no changes to the primary published global reconstruction.
"""
import argparse,csv,json
from pathlib import Path
from collections import defaultdict

def rows(path):
    with open(path,encoding="utf-8-sig",newline="") as f:
        yield from csv.DictReader(f)

def dist(path):
    out=defaultdict(set)
    for r in rows(path):
        id=(r.get("accepted_plant_name_id") or "").strip()
        region=(r.get("area_code_l3") or "").strip()
        if id and region:
            out[id].add(region)
    return out

def footprint(ids,native,contemporary):
    return {
        "native":set().union(*(native[id] for id in ids)),
        "contemporary":set().union(*(contemporary[id] for id in ids))
    }

def describe(x):
    n=x["native"]
    c=x["contemporary"]
    if not n.issubset(c):raise RuntimeError("Native-to-contemporary monotonicity broken")
    return {"native":len(n),"contemporary":len(c),
            "introduced_added":len(c-n)}

def plantago_loss(ids,plantago_id,native,contemporary):
    fp=footprint(ids,native,contemporary)
    after=footprint(ids-{plantago_id},native,contemporary)["contemporary"] | native[plantago_id]
    return fp["contemporary"]-after

def main():
    parser=argparse.ArgumentParser()
    for name in ["protocol_json","original_interaction_csv","original_native_csv",
                 "original_contemporary_csv","augmented_taxa_csv","augmented_native_csv",
                 "augmented_contemporary_csv","output_json"]:
        parser.add_argument("--"+name.replace("_","-"),type=Path,required=True)
    a=parser.parse_args()
    p=json.loads(a.protocol_json.read_text())
    if p["schema"]!="chocho_euphydryas_two_omitted_hosts_augmentation_v0.1":
        raise RuntimeError("Incorrect prespecified sensitivity schema")
    sp=p["frozen_inputs"]["butterfly"]
    local=p["literature_host_additions"]
    if [x["species"] for x in local]!=["Castilleja hispida","Castilleja levisecta"]:
        raise RuntimeError("Post-hoc species set must not change")
    interactions=[r for r in rows(a.original_interaction_csv)
                  if (r.get("insect_species") or "").strip()==sp]
    baseline_ids={r["accepted_plant_name_id"].strip() for r in interactions}
    if len(baseline_ids)!=p["frozen_inputs"]["baseline_taxa"]:
        raise RuntimeError("Original accepted host inventory differs")
    plantago={r["accepted_plant_name_id"].strip() for r in interactions
              if r.get("input_host_name")=="Plantago lanceolata" or
                  r.get("accepted_name")=="Plantago lanceolata"}
    if len(plantago)!=1:raise RuntimeError("Expected one exact original Plantago")
    pl=next(iter(plantago))
    n=dist(a.original_native_csv)
    c=dist(a.original_contemporary_csv)
    base=footprint(baseline_ids,n,c)
    metric=describe(base)
    frozen=p["frozen_inputs"]
    if metric!={"native":frozen["native_units"],
                "contemporary":frozen["contemporary_units"],
                "introduced_added":frozen["added_units"]}:
        raise RuntimeError("Frozen baseline mismatch: "+repr(metric))
    base_loss=plantago_loss(baseline_ids,pl,n,c)
    if len(base_loss)!=frozen["Plantago_lanceolata_sole_added_units"]:
        raise RuntimeError("Frozen original Plantago dependency mismatch")
    focal=[]
    additions=set()
    matching=list(rows(a.augmented_taxa_csv))
    if len(matching)!=2 or {x["host_name"] for x in matching}!={z["species"] for z in local}:
        raise RuntimeError("Botanical taxon match input was not exactly the two prespecified hosts")
    an=dist(a.augmented_native_csv)
    ac=dist(a.augmented_contemporary_csv)
    for r in matching:
        species=r["host_name"]
        status=r["match_status"]
        id=(r["accepted_plant_name_id"] or "").strip()
        if status=="UNIQUE_EXACT_ACCEPTED":
            if not id or id in baseline_ids or id in additions:
                raise RuntimeError("Focal new plant accepted ID not distinct")
            if an[id]-ac[id]:
                raise RuntimeError("New host native region missing in contemporary data")
            additions.add(id)
        else:
            if id:raise RuntimeError("Nonmatching taxon has resolved accepted ID")
        focal.append({
            "species":species,
            "match_state":status,
            "accepted_plant_name_id":id or None,
            "native_regions":len(an[id]) if id else None,
            "contemporary_regions":len(ac[id]) if id else None,
            "introduced_only_regions":len(ac[id]-an[id]) if id else None,
        })
    # The original all-Lepidoptera sidecar may include these plants through
    # ANOTHER butterfly/insect even when there is no E. editha link. Check its
    # frozen botanical range is identical rather than assuming it was absent.
    for id in additions:
        if (id in n and n[id] != an[id]) or (id in c and c[id] != ac[id]):
            raise RuntimeError("Target-host botanical range disagrees with frozen global sidecar")
        n[id]=an[id]
        c[id]=ac[id]
    ids=baseline_ids|additions
    aug=footprint(ids,n,c)
    metric2=describe(aug)
    newloss=plantago_loss(ids,pl,n,c)
    if newloss-base_loss:
        raise RuntimeError("Adding known hosts cannot increase absolute Plantago-sole resource cells")
    base_added=base["contemporary"]-base["native"]
    aug_added=aug["contemporary"]-aug["native"]
    result={
      "schema":"chocho_euphydryas_two_host_completeness_sensitivity_result_v0.1",
      "status":"POSTHOC_BOTANICAL_CROSSWALK_AND_STRUCTURAL_SENSITIVITY_NOT_DEMOGRAPHIC",
      "focal_butterfly":sp,
      "baseline":{"host_species":len(baseline_ids),**metric,
         "sole_Plantago_added_regions":len(base_loss),
         "sole_Plantago_share_of_added":len(base_loss)/len(base_added)},
      "augmented":{"host_species":len(ids),**metric2,
         "sole_Plantago_added_regions":len(newloss),
         "sole_Plantago_share_of_added":len(newloss)/len(aug_added) if aug_added else None},
      "original_missing_hosts":focal,
      "comparison":{
         "added_native_regions":len(aug["native"]-base["native"]),
         "added_contemporary_regions":len(aug["contemporary"]-base["contemporary"]),
         "previously_added_regions_reclassified_as_native":len(base_added-aug_added),
         "newly_introduced_only_regions":len(aug_added-base_added),
         "baseline_Plantago_sole_added_regions_now_rescued_by_Castilleja":len(base_loss-newloss),
         "baseline_Plantago_sole_added_regions_still_uncovered":len(base_loss&newloss),
         "total_delta_introduced_added":len(aug_added)-len(base_added)
      },
      "interpretation_boundaries":[
         "This is a 1-butterfly post-hoc missing-host completion, not revised 239-species findings.",
         "Castilleja association is from a specific published population and does not prove use across its entire mapped plant range.",
         "Native and introduced classifications are WCVP botanical range records, not local field resource quality.",
         "Classification as introduced-added can decrease if adding a missing host expands the native baseline.",
         "Structural host removal neither models actual plant eradication nor population survival or extinction."
      ]
    }
    a.output_json.parent.mkdir(parents=True,exist_ok=True)
    a.output_json.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False),flush=True)

if __name__=="__main__":
    main()
