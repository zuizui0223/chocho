#!/usr/bin/env python3
"""Source-verified LTMN 1976–1995 vs 2005–2021 historical detection pilot.

The unit is a VERIFIED, all-taxa-compatible survey *visit*, not a GBIF record.
Non-detection means NOT RECORDED during a focal-conditioned visit with >=3
OTHER identified butterflies. It NEVER proves local biological absence.
Source transect IDs and original visit dates are in-memory only.
"""
from __future__ import annotations
import argparse,csv,json,collections,hashlib
from pathlib import Path
from audit_uk_ltmn_historical_event_join_20261028 import (
    FIELDS_VISITS,FIELDS_COUNTS,get_source,excel_rows,parse_event,integral)
from audit_uk_ltmn_site_geography_20261028 import read_crosswalk

SCHEMA="chocho_uk_LTMN_19_source_transect_conditional_butterfly_detection_20261028_v01"
# Fixed one-to-one UKBMS common-name concordance, independent of source rows.
CROSSWALK={
"Aglais urticae":"Small Tortoiseshell",
"Anthocharis cardamines":"Orange-tip",
"Aphantopus hyperantus":"Ringlet",
"Argynnis aglaja":"Dark Green Fritillary",
"Argynnis paphia":"Silver-washed Fritillary",
"Callophrys rubi":"Green Hairstreak",
"Celastrina argiolus":"Holly Blue",
"Coenonympha pamphilus":"Small Heath",
"Colias croceus":"Clouded Yellow",
"Colias hyale":"Pale Clouded Yellow",
"Cupido minimus":"Small Blue",
"Erynnis tages":"Dingy Skipper",
"Gonepteryx rhamni":"Brimstone",
"Hesperia comma":"Silver-spotted Skipper",
"Hipparchia semele":"Grayling",
"Leptidea sinapis":"Wood White",
"Lycaena dispar":"Large Copper",
"Lycaena phlaeas":"Small Copper",
"Maniola jurtina":"Meadow Brown",
"Papilio machaon":"Swallowtail",
"Pararge aegeria":"Speckled Wood",
"Pieris brassicae":"Large White",
"Pieris napi":"Green-veined White",
"Pieris rapae":"Small White",
"Polyommatus icarus":"Common Blue",
"Vanessa atalanta":"Red Admiral",
"Vanessa cardui":"Painted Lady",
}
MOTH_OR_INSUFFICIENT={
"Burnet Companion","Chimney Sweeper","Clouded Buff","Common Heath",
"Five-spot Burnet","Humming-bird Hawk-moth","Latticed Heath","Mother Shipton",
"Narrow-bordered 5-spot Burnet","Narrow-bordered Bee Hawk-moth",
"Scarlet Tiger","Silver Y","Six-spot Burnet","Speckled Yellow",
"The Cinnabar","The Forester","The Vapourer","Wood Tiger","Yellow Belle",
"Unknown","Unknown Blue","Unknown Fritillary","Unknown White",
"No butterflies seen","Essex/Small Skipper","Small/Essex Skipper",
}
PERIODS={"early":(1976,1995),"late":(2005,2021)}
OTHER_SPECIES_MIN=3
MIN_QUALIFYING_SITES=10
MIN_QUALIFIED_VISITS_PER_ERA=10
MIN_OBS_DETECTED_ANY_ERA=20

def original_taxonomy_ok(path):
    src=json.loads(path.read_text(encoding="utf-8"))
    if src.get("schema")!="chocho_GBIF_butterfly_239_WGSRPD3_bounded_2010_2025_v01" or len(src.get("species_receipts",[]))!=239:
        raise ValueError("Chocho original 239-name source altered")
    lookup={x["species"]:x for x in src["species_receipts"]}
    if not set(CROSSWALK)<=set(lookup) or any(
         not isinstance(lookup[n].get("taxon",{}).get("key"),int)
         for n in CROSSWALK):
        raise ValueError("Target UK name crosswalk not in original accepted chocho taxon keys")
    if len(set(CROSSWALK.values()))!=len(CROSSWALK):
        raise ValueError("Duplicated source common-name concordance")
    if lookup["Lasiommata megera"]["status"]!="TAXON_UNRESOLVED":
        raise ValueError("Original source unresolved butterfly exclusion changed")
    return {"source_original_239_accepted_key_manifest_sha256":hashlib.sha256(path.read_bytes()).hexdigest(),
            "mapped_target_species":len(CROSSWALK),
            "unresolved_wall_taxon_is_excluded":True}

def site_effort(visits,geocoded):
    """Fixed, outcome-independent *original visit* site panel."""
    original=set();site_years=collections.defaultdict(lambda:{p:set() for p in PERIODS})
    period_events=collections.defaultdict(lambda:{p:set() for p in PERIODS})
    for r in visits:
        ev=parse_event(r)
        if ev is None or ev in original:continue
        original.add(ev)
        site,year,month,day,clock=ev
        if site not in geocoded or not 4<=month<=9:continue
        for era,(a,b) in PERIODS.items():
            if a<=year<=b:
                site_years[site][era].add(year)
                period_events[site][era].add(ev)
    panel={s for s,v in period_events.items() if
         all(len(v[p])>=10 and len(site_years[s][p])>=3 for p in PERIODS)}
    return original,panel,period_events

def count_at_visits(source_counts,verified_visits,panel):
    """Positives only; missing count is not a zero until visit protocol gate."""
    observed=collections.defaultdict(set)
    seen={}
    qc=collections.Counter()
    known_butterflies=set(CROSSWALK.values())
    for r in source_counts:
        ev=parse_event(r)
        if ev is None or ev not in verified_visits or ev[0] not in panel:continue
        name=str(r["COMMON_NAME"] or "").strip()
        try:count=float(r["VISIT_TOTAL"])
        except (ValueError,TypeError):qc["non_numeric_counts"]+=1;continue
        if not count>0:continue
        if name in MOTH_OR_INSUFFICIENT:continue
        # 2023 original UK source named all common butterflies. The remaining
        # named taxa are non-ambiguous UK butterflies; unknown source strings are
        # retained as other-target taxonomic evidence with source-id QC.
        key=(ev,name)
        if key in seen:
            if count!=seen[key]:qc["conflicting_source_count"]+=1
            else:qc["duplicate_count_row"]+=1
            continue
        seen[key]=count
        observed[ev].add(name)
        known_butterflies.add(name)
    qc["n_positive_original_named_visit_species_pairs"]=len(seen)
    qc["n_source_visit_events_with_at_least_one_named_butterfly"]=len(observed)
    qc["n_positive_source_common_names"]=len(known_butterflies)
    return observed,dict(qc)

def tally(visit_events,panel,observed,threshold=OTHER_SPECIES_MIN):
    answer=[]
    for sci,common in sorted(CROSSWALK.items()):
        bysite=collections.defaultdict(lambda:{p:[0,0] for p in PERIODS})
        seen_focal_era=collections.Counter()
        nsource_eligible=collections.Counter()
        for s in sorted(panel):
            for era,(start,end) in PERIODS.items():
                for ev in visit_events[s][era]:
                    names=observed.get(ev,set())
                    if sum(x!=common for x in names)<threshold:continue
                    nsource_eligible[era]+=1
                    bysite[s][era][0]+=1
                    if common in names:
                        bysite[s][era][1]+=1
                        seen_focal_era[era]+=1
        eligible=[s for s in panel if all(bysite[s][era][0]>=MIN_QUALIFIED_VISITS_PER_ERA for era in PERIODS)]
        stats={}
        for era in PERIODS:
            visits=sum(bysite[s][era][0] for s in eligible)
            nseen=sum(bysite[s][era][1] for s in eligible)
            fractions=[bysite[s][era][1]/bysite[s][era][0] for s in eligible]
            stats[era]={
                "source_qualified_focal_compatible_visits":visits,
                "focal_detected_visits":nseen,
                "site_equal_detected_fraction":sum(fractions)/len(fractions) if fractions else None,
                "n_original_site_panel_with_focal_detection":sum(bysite[s][era][1]>0 for s in eligible),
                "zero_is_visit_specific_not_occupancy_absence":True,
            }
        gate=(len(eligible)>=MIN_QUALIFYING_SITES and
              max(stats["early"]["focal_detected_visits"],stats["late"]["focal_detected_visits"])>=MIN_OBS_DETECTED_ANY_ERA)
        a,b=(stats[k]["site_equal_detected_fraction"] for k in PERIODS)
        answer.append({
            "species":sci,"original_UKBMS_common_name":common,
            "n_preselected_sites":len(panel),
            "n_sites_with_minimum_era_survey_coverage":len(eligible),
            "target_specific_eligible_source_visits_all_panel_by_era":dict(nsource_eligible),
            "periods":stats,
            "site_equal_late_minus_early_visit_detection_pp":100*(b-a) if a is not None and b is not None else None,
            "information_gate":"QUALIFIED_SOURCE_REPEAT_SITE_VISIT_DETECTION" if gate else "INSUFFICIENT_FOCAL_DETECTION_SUPPORT",
            "source_species_ambiguity_flag":"UKBMS Wood White potentially species complex" if sci=="Leptidea sinapis" else None,
            "migratory_species_not_evidence_of_local_recruitment":sci in {"Vanessa atalanta","Vanessa cardui","Colias croceus"},
        })
    return answer

def run(args):
    taxon=original_taxonomy_ok(args.original_manifest)
    geocoded,gq=read_crosswalk()
    original,panel,visit_events=site_effort(excel_rows(get_source("visits"),FIELDS_VISITS),geocoded)
    if len(panel)!=19:raise ValueError(f"Fixed source 19 long-term crosswalk panel changed: {len(panel)}")
    observed,qc=count_at_visits(excel_rows(get_source("counts"),FIELDS_COUNTS),original,panel)
    primary=tally(visit_events,panel,observed,3)
    strict=tally(visit_events,panel,observed,5)
    broad=tally(visit_events,panel,observed,1)
    result={
       "schema":SCHEMA,"status":"ORIGINAL_UK_VERIFIED_SURVEY_BASED_DETECTION_CHANGES_NOT_CAUSAL",
       "original_taxon_source":taxon,
       "source_original_LTMN_xls_sha256":{
          "visits":"c23945feb1568b4b60ac139136c485c234699710eae395a104e0be92a958a84d",
          "counts":"6afcc7d03eba0b657ae0e3a0d9da217b88df3ff39da66d7359df7fbb062ead3f",
          "site_crosswalk":"e38853b775ecadf39332a8412e2a4a2fa3c83512cf7218e4f1eda172dd83c9cd",
       },
       "source":{"n_original_exact_visits_all_years":len(original),
          "n_original_fixed_longterm_geocoded_transects":len(panel),
          "source_positive_count_qc":qc,"geography_source_qc":gq},
       "original_site_comparable_focal_butterflies":primary,
       "source_only_min_other_butterfly_labels_sensitivity":{
           "primary_min3":primary,"strict_min5":strict,"broad_min1":broad},
       "inference":{
          "real_ukbms_source_visit_denominator_not_GBIF_record_total":True,
          "event_with_missing_focal_name_is_nondetection_only_under_all_species_survey_assumption":True,
          "UKBMS_targeted_transects_cannot_be_independently_labelled_for_each_original_site":True,
          "not_biological_absence_or_colonization":True,
          "historical_local_host_plant_resource_not_measured":True,
          "human_maintained_host_plant_status_not_measured":True,
          "historical_annual_site_temperature_not_yet_joined":True,
          "no_host_mediated_warming_tracking_or_population_growth_inferred":True,
          "no_original_survey_ids_transect_coordinates_date_count_rows_exported":True,
       },
    }
    from collections import Counter
    result["source"]["n_primary_qualified_focal_species"]=sum(x["information_gate"].startswith("QUALIFIED") for x in primary)
    result["source"]["n_primary_source_taxa_with_at_least_one_detection_early"]=sum(x["periods"]["early"]["focal_detected_visits"]>0 for x in primary)
    result["source"]["n_primary_source_taxa_with_at_least_one_detection_late"]=sum(x["periods"]["late"]["focal_detected_visits"]>0 for x in primary)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({"status":result["status"],
                      "sites":len(panel),"species":len(CROSSWALK),
                      "qualified_species":result["source"]["n_primary_qualified_focal_species"],
                      "qualifying_source_deltas":[{"species":r["species"],
                        "delta_pp":r["site_equal_late_minus_early_visit_detection_pp"],
                        "sites":r["n_sites_with_minimum_era_survey_coverage"]}
                         for r in primary if r["information_gate"].startswith("QUALIFIED")]},
                      ensure_ascii=False),flush=True)
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--original-manifest",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    run(p.parse_args())
