#!/usr/bin/env python3
"""Recover public standardized Caterpillars Count! survey ID/site/year filters.

Field outcome is all Lepidoptera larvae, not butterfly species specific.
Survey sites, full coordinates and individual observer names never exported.
This extends the previously completed Table S2 cross-HOSTS join.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math
from collections import defaultdict,Counter
from pathlib import Path
import numpy as np
from analyze_field_caterpillar_host_realization_filter_20261027 import (
    normalize_binomial, accepted_global_links, solve_logit, expit, RIDGE,
)
from analyze_hosts_location_wcvp_origin_20261027 import (
    accepted_pairs, read_rows, sha256,
)

SCHEMA="chocho_external_field_plant_site_year_pairing_20261027_v01"
SOURCE_COMMIT="2c5690ac7c6dba7eca506da886c49781c37617bd"
FIELDS={"ID","Name","PlantFK","sciName","Species","julianday",
        "Longitude","Latitude","WetLeaves","Year","Group","Quantity",
        "ObservationMethod"}
COWEETA={"Coweeta - BB","Coweeta - BS","Coweeta - RK"}


def read_inferred(path):
    fields={}
    with Path(path).open(encoding="utf-8-sig",newline="") as f:
        it=csv.DictReader(f)
        need={"PlantFK","InferredSciName","NameConfidence"}
        if not need<=set(it.fieldnames or []):
            raise ValueError("Original plant inference source header drift")
        for r in it:
            pid=r["PlantFK"].strip()
            if not pid:continue
            try:confidence=float(r["NameConfidence"] or "nan")
            except ValueError:confidence=math.nan
            name=" ".join((r["InferredSciName"] or "").strip().split())
            if pid in fields and fields[pid]!=(name,confidence):
                raise ValueError("Original inferred plant FK has multiple contradictory identities")
            fields[pid]=(name,confidence)
    return fields


def read_origin(path):
    statuses={}
    with Path(path).open(encoding="utf-8-sig",newline="") as f:
        it=csv.DictReader(f)
        need={"scientificName","nativeStatus","plantOrigin"}
        if not need<=set(it.fieldnames or []):
            raise ValueError("Original USDA plant source schema drift")
        for r in it:
            name=" ".join((r["scientificName"] or "").strip().split())
            status=str(r["plantOrigin"] or "").strip().lower()
            if status not in ("native","alien"):
                continue
            if name in statuses and statuses[name]!=status:
                raise ValueError("Original field USDA plant origin conflicted")
            statuses[name]=status
    return statuses


def raw_filter(path, inferred, plant_origins, accepted_239):
    qc=Counter()
    surveys={}
    header=None
    all_years=Counter()
    with Path(path).open(encoding="utf-8-sig",newline="") as f:
        it=csv.DictReader(f)
        header=it.fieldnames
        if not header or not FIELDS<=set(header):
            raise ValueError("Original field fullDataset lacks expected raw survey source columns")
        for r in it:
            qc["original_raw_long_format_sighting_rows"]+=1
            site=(r["Name"] or "").strip()
            if site in COWEETA:qc["Coweeta_source_rows_excluded"]+=1;continue
            try:
                jd=int(float(r["julianday"]))
                lon=float(r["Longitude"])
                year=int(float(r["Year"]))
            except (ValueError,TypeError):
                qc["invalid_calendar_longitude"]+=1;continue
            if not (152<=jd<=212 and lon>=-100 and 2010<=year<=2025):
                qc["out_of_JuneJuly_eastNA_or_year"]+=1;continue
            wet=(r["WetLeaves"] or "").strip().lower()
            if wet not in ("false","f","0"):
                qc["wet_or_missing_leaf_condition"]+=1;continue
            surveyid=str(r["ID"]).strip()
            plantfk=str(r["PlantFK"] or "").strip()
            if not surveyid or not site or not plantfk:
                qc["missing_source_survey_id_site_or_branch"]+=1;continue
            scientific=" ".join((r.get("sciName") or "").strip().split())
            inference=inferred.get(plantfk)
            if (str(r.get("Species") or "").strip()=="N/A"
                and inference and inference[1]>=2.):
                scientific=inference[0]
            if not scientific or scientific.lower() in ("na","n/a","unknown"):
                qc["no_source_scientific_plant_name"]+=1;continue
            botanic=normalize_binomial(scientific)
            taxon_status=("SPECIES_BINOMIAL" if botanic else "SOURCE_NONBINOMIAL_TAXON")
            plant_name=botanic or scientific
            origin=plant_origins.get(plant_name)
            if origin not in ("native","alien"):
                qc["no_original_USDA_origin_for_source_plant"]+=1;continue
            if not (-180<=lon<=180):
                qc["invalid_longitude"]+=1;continue
            check=(site,year,jd,plantfk,plant_name,origin,taxon_status)
            if surveyid not in surveys:
                surveys[surveyid]={"key":check,"caterpillar":False,
                                   "record_rows":0}
            elif surveys[surveyid]["key"]!=check:
                qc["survey_ID_conflicting_source_keys"]+=1
                continue
            p=surveys[surveyid]
            p["record_rows"]+=1
            if str(r.get("Group") or "").strip().lower()=="caterpillar":
                try:
                    quantity=float(str(r.get("Quantity") or "0"))
                except ValueError:
                    quantity=0.
                if quantity>0:
                    p["caterpillar"]=True
            all_years[year]+=1
    source_counts=Counter()
    surveyrows=[]
    for surveyid,info in surveys.items():
        site,year,jd,plantfk,plant_name,origin,taxon_status=info["key"]
        source_counts["survey_IDs_final_quality_gate"]+=1
        source_counts[f"USDA_{origin}_survey_IDs"]+=1
        source_counts["source_caterpillar_positive_survey_IDs"]+=int(info["caterpillar"])
        match=accepted_239.get(plant_name)
        if taxon_status=="SPECIES_BINOMIAL" and isinstance(match,int):
            source_counts["accepted_239_linked_survey_IDs"]+=1
        else:
            match=None
        surveyrows.append({
            "site":site,"year":year,"week":jd//7,
            "origin":origin,"plant":plant_name,
            "taxon_status":taxon_status,
            "caterpillar":int(info["caterpillar"]),
            "host_239":match,
        })
    source_counts["unique_source_sites"]=len({r["site"] for r in surveyrows})
    source_counts["unique_source_plant_names"]=len({r["plant"] for r in surveyrows})
    source_counts["source_distinct_calendar_years"]=len({r["year"] for r in surveyrows})
    return surveyrows,dict(qc),dict(source_counts),header


def site_year_matched(surveys, *, min_survey=5, week=False, source_links=False):
    groups=defaultdict(list)
    for row in surveys:
        if source_links and row["host_239"] is None:
            continue
        key=(row["site"],row["year"],row["week"]) if week else (row["site"],row["year"])
        groups[key].append(row)
    screened=[]
    for key,group in groups.items():
        native=[r for r in group if r["origin"]=="native"]
        alien=[r for r in group if r["origin"]=="alien"]
        if len(native)<min_survey or len(alien)<min_survey:
            continue
        if not {r["plant"] for r in native} or not {r["plant"] for r in alien}:
            continue
        screened.append({
            "site":key[0],"site_year_key":key,
            "native_surveys":len(native),
            "native_positive":sum(r["caterpillar"] for r in native),
            "alien_surveys":len(alien),
            "alien_positive":sum(r["caterpillar"] for r in alien),
            "native_distinct_plants":len({r["plant"] for r in native}),
            "alien_distinct_plants":len({r["plant"] for r in alien}),
        })
    return screened


def aggregate_matched(rows):
    n=len(rows)
    sites=len({x["site"] for x in rows})
    native_n=sum(x["native_surveys"] for x in rows)
    alien_n=sum(x["alien_surveys"] for x in rows)
    native_pos=sum(x["native_positive"] for x in rows)
    alien_pos=sum(x["alien_positive"] for x in rows)
    info={
        "matched_site_year_or_week_strata":n,
        "independent_field_sites_with_matching_both_plant_origins":sites,
        "source_native_survey_events_in_matched_contexts":native_n,
        "source_alien_survey_events_in_matched_contexts":alien_n,
        "native_larval_positive_events":native_pos,
        "alien_larval_positive_events":alien_pos,
        "descriptive_native_detection_frac":native_pos/native_n if native_n else None,
        "descriptive_alien_detection_frac":alien_pos/alien_n if alien_n else None,
        "source_site_year_paired_evidence_information_gate":(
            "SUFFICIENT_CONTEXT_STRATA_FOR_DESCRIPTIVE_PAIRING"
            if n>=20 and sites>=10 else "INSUFFICIENT_MATCHED_SITES_OR_YEARS"),
    }
    # Mantel-Haenszel source-ascertained within-site/year odds ratio (descriptive).
    if n>=20 and sites>=10:
        upper=lower=0.
        for r in rows:
            a=r["native_positive"]; b=r["native_surveys"]-a
            c=r["alien_positive"]; d=r["alien_surveys"]-c
            t=a+b+c+d
            upper+=a*d/t; lower+=b*c/t
        info["descriptive_Mantel_Haenszel_native_vs_alien_larval_presence_odds_ratio"]=(upper/lower if lower else None)
    return info


def run(args):
    author=json.loads((args.author_source/"SOURCE_MANIFEST.json").read_text())
    if author.get("inputs_count")!=11:raise ValueError("Not original WCVP/HOSTS manifest")
    for name,digest in author["sha256"].items():
        if sha256(args.author_source/name)!=digest:
            raise ValueError("Original author-only input member SHA mismatch")
    original=read_rows(args.author_source/"inputs/anthropogenic_species_metrics.csv")
    focal={r["species"] for r in original}
    if len(original)!=239 or len(focal)!=239:
        raise ValueError("Original 239 frozen original species changed")
    by_name,by_id,ambiguous,meta=accepted_global_links(
        args.author_source/"inputs/original/insect_host_accepted.csv",focal)
    selected={
        name:len(by_id[next(iter(ids))])
        for name,ids in by_name.items()
        if len(ids)==1 and name not in ambiguous
    }
    original_origin=read_origin(args.origin)
    original_inferred=read_inferred(args.inferred)
    surveys,qc,source,cols=raw_filter(args.original_raw,original_inferred,
                                     original_origin,selected)
    siteyear=site_year_matched(surveys,min_survey=5)
    siteweek=site_year_matched(surveys,min_survey=2,week=True)
    siteyear_239=site_year_matched(surveys,min_survey=5,source_links=True)
    matched_plant_names={r["plant"] for r in surveys if r["host_239"] is not None}
    native_plants={r["plant"] for r in surveys if r["host_239"] is not None
                   and r["origin"]=="native"}
    alien_plants={r["plant"] for r in surveys if r["host_239"] is not None
                  and r["origin"]=="alien"}
    new={
        "schema":SCHEMA,
        "status":"RAW_FIELD_SOURCE_SITE_YEAR_PAIRED_EVIDENCE_QC_NOT_BUTTERFLY_RECRUITMENT",
        "source":{
            "external_commit":SOURCE_COMMIT,
            "original_raw_field_csv_sha256":sha256(args.original_raw),
            "original_USDA_plant_status_sha256":sha256(args.origin),
            "original_inferred_botanical_names_sha256":sha256(args.inferred),
            "original_field_raw_columns":cols,
            "original_239_butterfly_host_graph_sha256":sha256(
                args.author_source/"inputs/original/insect_host_accepted.csv"),
            "original_239_accepted_host_links":meta["original_239_accepted_species_plant_pairs"],
        },
        "raw_long_format_quality":qc,
        "source_unique_surveys":source,
        "matched_site_year_both_origins":aggregate_matched(siteyear),
        "matched_site_7day_week_both_origins":aggregate_matched(siteweek),
        "original_239_source_exact_plant_linked_site_year_both_origins":aggregate_matched(siteyear_239),
        "source_original239_accepted_plant_match_coverage":{
            "species_name_linked_plants_at_any_site":len(matched_plant_names),
            "native_plant_names_with_original_239_link":len(native_plants),
            "alien_plant_names_with_original_239_link":len(alien_plants),
            "source_original_239_linked_species_model_gate":(
                "POTENTIALLY_TESTABLE_FROM_SITE_YEAR" if len(matched_plant_names)>=15
                and len(native_plants)>=3 and len(alien_plants)>=3
                and aggregate_matched(siteyear_239)["matched_site_year_or_week_strata"]>=10
                and aggregate_matched(siteyear_239)[
                    "independent_field_sites_with_matching_both_plant_origins"]>=10
                else "INSUFFICIENT_SOURCE_239_MATCHED_PLANT_SITE_YEAR_SUPPORT"),
        },
        "red_lines":[
            "The raw source field branch surveys include moths and butterflies; not a measure of invasive butterfly species host use.",
            "Native/alien survey effects were published in Whitener and Hurlbert 2026; this analysis tests feasibility of source HOSTS link transfer, not novel estimation of published effects.",
            "Site/year pairing does not randomize botanical origin, planted host choice or observer effort and does not establish ecological causality.",
            "Source negative occurrence is one survey's nondetection, not field larval survival failure or host unsuitability.",
            "A source genus-only plant cannot be linked to a unique original accepted WCVP host species.",
            "No exact site names, observer personal details, coordinates, individual arthropod observations, or raw per-species site-year panels are exported.",
        ],
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(new,sort_keys=True,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps({
        "status":new["status"],
        "qc":qc,
        "source":source,
        "site_year":new["matched_site_year_both_origins"],
        "week":new["matched_site_7day_week_both_origins"],
        "source_239_sites":new["original_239_source_exact_plant_linked_site_year_both_origins"],
        "source_239_coverage":new["source_original239_accepted_plant_match_coverage"],
    },sort_keys=True))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--author-source",required=True,type=Path)
    p.add_argument("--original-raw",required=True,type=Path)
    p.add_argument("--origin",required=True,type=Path)
    p.add_argument("--inferred",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    run(p.parse_args())
