#!/usr/bin/env python3
"""One-direction botanical resource sensitivity for the BCE/Clarke host gap.

An evidence-ranked European published host relation is extrapolated across that
plant's WCVP regions ONLY for this deliberately upper-bound sensitivity.
No causal ecological, demographic or global 'corrected' estimate is asserted.
"""
from __future__ import annotations
import argparse,csv,json
from collections import Counter,defaultdict
from pathlib import Path

SCHEMA="chocho_bce_clarke92_wcvp_accepted_id_gap_v0.1"
TYPE="ACCEPTED_ID_ABSENT_FROM_FROZEN_HOSTS_FOR_BUTTERFLY"

def csvrows(path):
    with path.open(encoding="utf-8-sig",newline="") as f:
        yield from csv.DictReader(f)

def load_source(path):
    r=json.loads(path.read_text(encoding="utf-8"))
    if r.get("schema")!=SCHEMA or r.get("exact_focal_butterflies")!=92:
        raise RuntimeError("WCVP reconciliation isn't the frozen full 92")
    return r

def additions(src):
    out=defaultdict(set)
    raw_count=0
    for butterfly in src["per_butterfly"]:
        b=butterfly["species"]
        for row in butterfly["details"]:
            if row["match_status"]==TYPE:
                id=row.get("accepted_plant_name_id")
                if not id:raise RuntimeError("Missing accepted ID")
                out[b].add(str(id))
                raw_count+=1
    return out,raw_count

def native_current(csv_native,csv_contemporary):
    def read(path):
        d=defaultdict(set)
        for row in csvrows(path):
            key=row["accepted_plant_name_id"].strip()
            region=row["area_code_l3"].strip()
            if key and region:d[key].add(region)
        return d
    n=read(csv_native);c=read(csv_contemporary)
    for key,v in n.items():
        if not v.issubset(c[key]):
            raise RuntimeError("Native not contained in contemporary for "+key)
    return n,c

def footprint(ids,mapping):
    result=set()
    for id in ids:result.update(mapping[id])
    return result

def source_rows(path):
    return {row["species"].strip():row for row in csvrows(path)}

def prepare(src,args):
    a,links=additions(src)
    ids=sorted(set().union(*a.values())) if a else []
    if not ids:
        raise RuntimeError("No uniquely resolved accepted additional hosts")
    args.output_csv.parent.mkdir(parents=True,exist_ok=True)
    with args.output_csv.open("w",encoding="utf-8",newline="") as f:
        w=csv.writer(f);w.writerow(["accepted_plant_name_id"])
        w.writerows([[id] for id in ids])
    print(json.dumps({"accepted_extra_link_records":links,"extra_accepted_pair_count":sum(map(len,a.values())),
                      "additional_accepted_plant_species":len(ids)}),flush=True)

def audit(src,args):
    added,original_link_count=additions(src)
    baseline=defaultdict(set)
    for row in csvrows(args.insect_host_csv):
        sp=row["insect_species"].strip()
        id=row["accepted_plant_name_id"].strip()
        if sp and id:baseline[sp].add(id)
    original_native,original_cont=native_current(args.native_csv,args.contemporary_csv)
    added_native,added_cont=native_current(args.added_native_csv,args.added_contemporary_csv)
    for id in set().union(*added.values()):
        if id in original_native and original_native[id]!=added_native[id]:
            raise RuntimeError("Botanical native range discordance for "+id)
        if id in original_cont and original_cont[id]!=added_cont[id]:
            raise RuntimeError("Botanical contemporary range discordance for "+id)
        original_native[id]=added_native[id]
        original_cont[id]=added_cont[id]
    frozen=source_rows(args.frozen_metrics_csv)
    if len(frozen)!=239:raise RuntimeError("Original 239 butterfly panel changed")
    if not set(added).issubset(frozen):
        raise RuntimeError("BCE missing edge butterfly outside frozen 239")
    species_results=[]
    mismatch=[]
    all_native_before=all_current_before=all_added_before=0
    all_native_after=all_current_after=all_added_after=0
    for sp,row in sorted(frozen.items()):
        ids=baseline[sp]
        if not ids:raise RuntimeError("Frozen species without accepted host "+sp)
        old_n=footprint(ids,original_native)
        old_c=footprint(ids,original_cont)
        old_a=old_c-old_n
        old=(len(old_n),len(old_c),len(old_a))
        expected=(int(row["native_resource_units"]),int(row["contemporary_resource_units"]),int(row["introduced_added_units"]))
        if old!=expected:
            mismatch.append({"species":sp,"observed":old,"frozen":expected})
        newids=ids | added[sp]
        new_n=footprint(newids,original_native)
        new_c=footprint(newids,original_cont)
        new_a=new_c-new_n
        if not old_n.issubset(new_n) or not old_c.issubset(new_c):
            raise RuntimeError("Augmentation lost old resource geography")
        if not new_n.issubset(new_c):
            raise RuntimeError("Augmentation native union not subset contemporary")
        all_native_before+=len(old_n);all_current_before+=len(old_c);all_added_before+=len(old_a)
        all_native_after+=len(new_n);all_current_after+=len(new_c);all_added_after+=len(new_a)
        species_results.append({
            "species":sp,"BCE_unique_accepted_additional_hosts":len(added[sp]),
            "baseline_host_species":len(ids),"augmented_host_species":len(newids),
            "baseline_native":len(old_n),"augmented_native":len(new_n),
            "baseline_contemporary":len(old_c),"augmented_contemporary":len(new_c),
            "baseline_introduced_added":len(old_a),"augmented_introduced_added":len(new_a),
            "new_native_regions":len(new_n-old_n),
            "new_contemporary_regions":len(new_c-old_c),
            "old_introduced_added_regions_reclassified_native":len(old_a-new_a),
            "newly_introduced_added_regions":len(new_a-old_a)
        })
    if mismatch:
        raise RuntimeError("Frozen 239 baseline differs from reconstructed input: "+json.dumps(mismatch[:12]))
    if all_added_before!=14553 or sum(r["baseline_introduced_added"]>0 for r in species_results)!=206:
        raise RuntimeError("Frozen 14,553/206 global checks disagree")
    observed=set(b["species"] for b in src["per_butterfly"])
    subset=[r for r in species_results if r["species"] in observed]
    changed=[r for r in species_results if r["augmented_host_species"]>r["baseline_host_species"]]
    output={
        "schema":"chocho_bce_clarke_one_direction_european_augmentation_v0.1",
        "status":"POSTHOC_STRUCTURAL_PARTIAL_EUROPE_BUTTERFLY_HOST_ADDITION_SENSITIVITY_NOT_REVISED_GEB",
        "frozen_239_baseline_verified":True,
        "frozen_species":239,"exact_European_butterflies":92,
        "verified_BCE_host_species_pages":src["verified_butterflies"],
        "additional_unique_accepted_host_links":sum(len(v) for v in added.values()),
        "additional_raw_reconciled_source_records":original_link_count,
        "butterflies_with_additional_host_links":len(changed),
        "all_239_descriptive_partial_augmented":{
            "baseline_native":all_native_before,"augmented_native":all_native_after,
            "baseline_contemporary":all_current_before,"augmented_contemporary":all_current_after,
            "baseline_introduced_added":all_added_before,"augmented_introduced_added":all_added_after,
            "baseline_resource_increase_fraction":all_added_before/all_native_before,
            "partial_augmented_increase_fraction":all_added_after/all_native_after,
            "baseline_species_with_expansion":206,
            "partial_augmented_species_with_expansion":sum(r["augmented_introduced_added"]>0 for r in species_results)
        },
        "european_overlap_subset":{
            "species":len(subset),
            "species_with_host_additions":sum(r["augmented_host_species"]>r["baseline_host_species"] for r in subset),
            "baseline_native":sum(r["baseline_native"] for r in subset),
            "augmented_native":sum(r["augmented_native"] for r in subset),
            "baseline_introduced_added":sum(r["baseline_introduced_added"] for r in subset),
            "augmented_introduced_added":sum(r["augmented_introduced_added"] for r in subset),
        },
        "species_with_changed_geo_counts":[x for x in species_results if
            (x["new_native_regions"] or x["new_contemporary_regions"])],
        "all_239_species_results":species_results,
        "limits":[
            "Europe-specific source is BCE-derived Clarke 2024 and not globally complete.",
            "Additional accepted plant IDs do not confirm local usage, developmental success or colonization in every WCVP region.",
            "The partial 239-panel recomputation is a sensitivity illustration, not a corrected representative global result.",
            "Adding native host links may make formerly introduced-added cells native, even when total coverage increases.",
            "No null significance or causal fitness/competition/consumer response inference.",
            "Original GEB metrics and figures remain unchanged."
        ]
    }
    args.output_json.parent.mkdir(parents=True,exist_ok=True)
    args.output_json.write_text(json.dumps(output,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print("GEOGRAPHIC_SENSITIVITY_SUMMARY",json.dumps({
        k:output[k] for k in ["frozen_239_baseline_verified","additional_unique_accepted_host_links","butterflies_with_additional_host_links","all_239_descriptive_partial_augmented","european_overlap_subset"]}),flush=True)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("mode",choices=("prepare","audit"))
    for name in ("crosswalk_json","output_csv","output_json","insect_host_csv","native_csv",
                 "contemporary_csv","added_native_csv","added_contemporary_csv","frozen_metrics_csv"):
        parser.add_argument("--"+name.replace("_","-"),type=Path)
    a=parser.parse_args()
    if a.crosswalk_json is None:parser.error("Need --crosswalk-json")
    src=load_source(a.crosswalk_json)
    if a.mode=="prepare":
        if a.output_csv is None:parser.error("Need --output-csv")
        prepare(src,a)
    else:
        needed=["output_json","insect_host_csv","native_csv","contemporary_csv",
                "added_native_csv","added_contemporary_csv","frozen_metrics_csv"]
        if any(getattr(a,k) is None for k in needed):parser.error("Missing audit inputs")
        audit(src,a)

if __name__=="__main__":main()
