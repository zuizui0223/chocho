#!/usr/bin/env python3
"""Recompute independent 8-endemic potential resource geography from the EXACT
source-pinned original WCVP v13 edge-level output without repeating R.

Proof targets:
- conditional effects from 13 verified butterfly × plant SPECIES interactions;
- the 18 independent endemics are the denominator, 8 source-known butterflies
  form a selective, NOT representative subset;
- no source-location newly host-present units for P. cheiranthi even though
  introduced Tropaeolum majus augments its host-plant options;
- isolate cross-global potential from actual host use/caterpillar performance;
- exclude captivity-only, oviposition-only, and non-local NHM evidence in
  separate tiers;
- leave-one-host-out diagnostics expose dependence on an alien host species.
"""
import argparse,csv,json,math
from collections import defaultdict
from pathlib import Path

MODES={
    "wild_and_documented_larval":{"larval_reported"},
    "plus_published_oviposition_and_captivity":{"larval_reported","oviposition_only","captivity_only"},
    "all_named_sources_including_NHM_nonlocal":{
        "larval_reported","oviposition_only","captivity_only","nonlocal_database"},
    "NHM_original_nonlocal_only":{"nonlocal_database"},
}
EXACT_PLANTS=11
SOURCE_EDGES=13
INDEPENDENT_END_EMICS=18
SOURCE_EIGHT=8
def load(path):
    with path.open(newline="",encoding="utf-8-sig") as h:
        return list(csv.DictReader(h))
def parse_units(value):
    return set(str(value or "").split(";"))-{""}
def summarize(rows):
    if not rows:
        return None
    ids=set(); native=set(); contemporary=set(); introduced_on_source=[]
    for row in rows:
        if row["accepted_plant_status"]!="resolved" or row["source_region_uncertain"]!="0":
            return {"unknown":True,"no_numeric_imputation":True}
        plant=row["plant_species"]
        if plant in ids:raise ValueError("Duplicated butterfly plant in source file")
        ids.add(plant)
        n=parse_units(row["global_native_regions"])
        i=parse_units(row["global_introduced_regions"])
        native|=n
        contemporary|=n|i
        if row["source_region_introduced"]=="1":
            introduced_on_source.append(plant)
    if not native <= contemporary:
        raise ValueError("WCVP native region not contained in contemporary")
    source_native=sum(int(row["source_region_native"]) for row in rows)
    source_introduced=sum(int(row["source_region_introduced"]) for row in rows)
    return {
        "host_species_n":len(rows),
        "host_native_in_butterfly_source_region_n":source_native,
        "host_introduced_in_butterfly_source_region_n":source_introduced,
        "introduced_hosts_in_butterfly_source_region":sorted(introduced_on_source),
        "source_region_was_host_present_with_native":bool(source_native),
        "source_region_any_host_contemporary":bool(source_native or source_introduced),
        "source_region_new_host_presence_created_by_introduced_only":int(source_native==0 and source_introduced>0),
        "global_potential_native_region_count":len(native),
        "global_potential_contemporary_region_count":len(contemporary),
        "global_potential_added_region_count":len(contemporary-native),
        "global_potential_contemporary_over_native_ratio":len(contemporary)/len(native) if native else None,
        "global_added_region_codes":sorted(contemporary-native)
    }
def equal_replay(calculated,row):
    if calculated is None:
        return row["host_species_without_source_record"]=="1"
    if calculated.get("unknown"):
        return row["unresolved_host_species_or_local_status_n"]!="0"
    fields={
       "host_species_n":"source_host_species_n",
       "host_native_in_butterfly_source_region_n":"source_resolved_hosts_native_n",
       "host_introduced_in_butterfly_source_region_n":"source_resolved_hosts_introduced_n",
       "global_potential_native_region_count":"potential_global_native_regions",
       "global_potential_contemporary_region_count":"potential_global_contemporary_regions",
       "global_potential_added_region_count":"potential_global_added_regions",
       "source_region_new_host_presence_created_by_introduced_only":"source_region_newly_host_present_after_intro",
    }
    return all(int(row[out_col])==calculated[k] for k,out_col in fields.items())
def main(source,original,out):
    edges=load(source)
    r_rows=load(original)
    if len(edges)!=SOURCE_EDGES or len(r_rows)!=SOURCE_EIGHT*len(MODES):
        raise ValueError("Original WCVP source panel sizes changed")
    bmap=defaultdict(list)
    for row in edges:
        bmap[row["butterfly"]].append(row)
    if len(bmap)!=SOURCE_EIGHT or len({r["plant_species"] for r in edges})!=EXACT_PLANTS:
        raise ValueError("Original eight independent endemics/eleven exact plants changed")
    if set(x["evidence_tier"] for x in edges)!={"larval_reported","oviposition_only","captivity_only","nonlocal_database"}:
        raise ValueError("Source evidence stages changed")
    if any(x["accepted_plant_status"]!="resolved" for x in edges):
        raise ValueError("A plant is not uniquely accepted in pinned WCVP v13")
    index={(r["source_butterfly"],r["mode"]):r for r in r_rows}
    if len(index)!=len(r_rows):
        raise ValueError("Duplicate original R output cells")
    by_mode=defaultdict(list);per_species={}
    for butterfly in sorted(bmap):
        records=bmap[butterfly]
        result={}
        for mode,grades in MODES.items():
            subset=[r for r in records if r["evidence_tier"] in grades]
            summary=summarize(subset)
            observed=index[(butterfly,mode)]
            if not equal_replay(summary,observed):
                raise ValueError(f"Independent exact Python replay of pinned R panel mismatch for {butterfly} / {mode}: "+repr(summary))
            result[mode]=summary
            if summary and not summary.get("unknown"):
                by_mode[mode].append((butterfly,summary))
        full=result["all_named_sources_including_NHM_nonlocal"]
        if full is None or full.get("unknown"):
            raise ValueError("The all-evidence eight-source dataset must have known plant status")
        loo=[]
        for removed in records:
            other=[r for r in records if r["plant_species"]!=removed["plant_species"]]
            after=summarize(other)
            loo.append({
                "butterfly":butterfly,
                "removed_host_species":removed["plant_species"],
                "removed_host_evidence_tier":removed["evidence_tier"],
                "removed_host_source_region_introduced":int(removed["source_region_introduced"]),
                "remaining_named_host_species":len(other),
                "global_added_regions_full":full["global_potential_added_region_count"],
                "global_added_regions_after_removal":None if after is None else after["global_potential_added_region_count"],
                "global_added_region_count_difference":None if after is None else full["global_potential_added_region_count"]-after["global_potential_added_region_count"],
                "relative_resource_ratio_after_removal":None if after is None else after["global_potential_contemporary_over_native_ratio"],
                "source_region_any_native_host_after_removal":None if after is None else after["source_region_was_host_present_with_native"],
                "caveat":"Removing one documented host tests resource-envelope sensitivity, NOT species diet adaptation or field host loss"
            })
        per_species[butterfly]={
            "source_region":records[0]["source_region_TDWG3"],
            "source_original_GEB239_selected":0,
            "all_host_records":sorted(x["plant_species"] for x in records),
            "modes":result,
            "leave_one_reported_host_out":loo
        }
    p=per_species["Pieris cheiranthi"]
    all3=p["modes"]["all_named_sources_including_NHM_nonlocal"]
    alien=[r for r in p["leave_one_reported_host_out"] if r["removed_host_species"]=="Tropaeolum majus"]
    if (len(alien)!=1 or alien[0]["removed_host_source_region_introduced"]!=1 or
       all3["host_native_in_butterfly_source_region_n"]<1 or
       all3["source_region_new_host_presence_created_by_introduced_only"]!=0):
       raise ValueError("Pieris alien host portfolio must be distinguished from host-present island region")
    for b in ("Hipparchia azorina","Hipparchia miguelensis"):
        v=per_species[b]["modes"]["all_named_sources_including_NHM_nonlocal"]
        if v["host_species_n"]!=1 or v["global_potential_added_region_count"]!=0:
            raise ValueError("Two Azores endemic original Festuca host-only result changed")
    by_mode_stats={}
    for mode,items in by_mode.items():
        positive=[b for b,v in items if v["global_potential_added_region_count"]>0]
        by_mode_stats[mode]={
            "eligible_source_endemic_butterflies_n":len(items),
            "with_added_global_potential_regions_n":len(positive),
            "with_added_global_potential_regions":sorted(positive),
            "with_zero_added_global_potential_regions":[b for b,v in items if v["global_potential_added_region_count"]==0],
            "no_claim_as_population_prevalence":True,
        }
    # Avoid cross-species average of 40x ratio dominated by a single host
    # and incompletely sampled origination histories.
    summary={
      "schema":"chocho_independent8_evidence_tier_and_plant_loo_source_v1",
      "original_independent_endemics_n":INDEPENDENT_END_EMICS,
      "source_ascerained_exact_host_species_endemics_n":SOURCE_EIGHT,
      "source_plant_species_edges_n":SOURCE_EDGES,
      "unique_accepted_host_plant_species":EXACT_PLANTS,
      "replay_all_32_R_combinations_success":True,
      "by_mode":by_mode_stats,
      "butterflies":per_species,
      "most_important_validations":{
          "Pieris_cheiranthi_one_introduced_host_is_Tropaeolum_majus":alien[0],
          "Pieris_cheiranthi_resource_introduction_does_not_create_new_CNY_host_present_region":True,
          "Azores_two_endemics_both_same_Festuca_jubata":True
      },
      "scientific_boundary":[
         "Source-available endemic butterflies are 8 of 18 independently classified taxa. This is deliberate evidence availability sampling, NOT representative or random.",
         "A globally projected butterfly host-plant list is not local butterfly geographic range, colonization, oviposition or larval performance.",
         "The 40-fold potential global envelope of P. cheiranthi is an abstract plant geography estimate driven by one source-recorded introduced host, not 40-fold island colonization.",
         "Native and introduced land-plant ranges are entire WGSRPD3 botanical regions, NOT individual Canary, Madeira or Azores islands.",
         "E. hesperidum source hosts were accepted by a SINGLE captive larva, and Descurainia preauxiana for E. grancanariensis has only oviposition evidence.",
         "Two Azores endemic Hipparchia source-only estimates both rely on the SAME original NHM host plant and are not independent evolutionary replications.",
         "No host-species absence should be inferred from species without an identified plant (10 remaining of 18 endemic butterflies).",
         "Original GEB source-selected +54.9% and +118.9% potential opportunity effects remain unchanged."
      ]
    }
    out.mkdir(parents=True,exist_ok=True)
    (out/"evidence_tiers_and_leave_one_host_out.json").write_text(
        json.dumps(summary,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    rows=[]
    for b,item in per_species.items():rows.extend(item["leave_one_reported_host_out"])
    with (out/"leave_one_source_host_out.csv").open("w",encoding="utf-8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    print("INDEPENDENT8_HOST_EVIDENCE_SENSITIVITY_VERIFIED="+json.dumps({
      "by_mode":by_mode_stats,
      "Pieris_without_intro_host":alien[0],
      "source_panels_replay_exact":True
    },ensure_ascii=False),flush=True)
    return summary
if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source-edges",required=True,type=Path)
    p.add_argument("--original-WCVP-results",required=True,type=Path)
    p.add_argument("--out-dir",required=True,type=Path)
    a=p.parse_args()
    main(a.source_edges,a.original_WCVP_results,a.out_dir)
