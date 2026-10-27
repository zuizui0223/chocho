#!/usr/bin/env python3
"""Positive-only match: published 2024 established alien insects vs chocho 239.

Zenodo 2024 v0.2.1 occurrence_table_20230317.csv has 13047
establishment species×region positive rows; NOT a nonestablishment census.
Conservative exact binomial match to original 239 butterfly HOSTS/WCVP
species summaries. Missing Zenodo records are NOT arrival/establishment
failures. Both source SHA and original 239 cohort metrics verified.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import re
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path

import numpy as np

from audit_zenodo_2024_nonnative_insect_establishment_source_20261027 import get

SCHEMA="chocho_239_butterfly_original_Zenodo_positive_alien_establishment_20261027_v01"
ZENODO_RECORD=11519262
ZENODO_META_URL=f"https://zenodo.org/api/records/{ZENODO_RECORD}"
DATA={
    "occurrence_table_20230317.csv":"05483a2a9a7c926fb5105c44058ee52a4598039dc82c5bb10d1db284b3b66fe7",
    "attribute_table_20230317.csv":"d184315959746f867fef7f377d5532c383cc7d7fc1b3e135f6bf967c586275c5",
    "source_table_20230317.csv":"ad2258907d07e6216b3a708d8bdee8780c6ba66a0b8046952b9342575e273ec4",
}
RAW_ROWS={"occurrence_table_20230317.csv":13047,
          "attribute_table_20230317.csv":8514,
          "source_table_20230317.csv":13210}
EXPECTED_ORIG_SPECIES=239
ORIG_RESOURCE_GAIN=14553
SOURCE_UNLISTED="NOT_DOCUMENTED_AS_ALIEN_ESTABLISHED_IN_2024_SOURCE_REGIONS"


def sha(content):
    return hashlib.sha256(content).hexdigest()


def binomial(text):
    r=" ".join(str(text or "").strip().split())
    if re.fullmatch(r"[A-Z][a-zA-Z-]+ [a-z][a-z-]+",r):
        return r
    return None


def read_zenodo_member(name,item):
    link=(item.get("links") or {}).get("self") or (
        item.get("links") or {}).get("content")
    if not link:
        raise ValueError("Zenodo original species table missing authorized download")
    payload,url,ctype=get(link,max_bytes=9_000_000)
    if sha(payload)!=DATA[name]:
        raise ValueError("Zenodo 2024 published original member SHA mismatch "+name)
    if item.get("size") and len(payload)!=int(item["size"]):
        raise ValueError("Zenodo original member size mismatch "+name)
    declared=item.get("checksum")
    if declared and declared.startswith("md5:") and \
       hashlib.md5(payload).hexdigest()!=declared.partition(":")[2]:
        raise ValueError("Zenodo original member checksum mismatch "+name)
    try:
        text=payload.decode("utf-8-sig")
    except UnicodeDecodeError:
        text=payload.decode("latin-1")
    rows=list(csv.DictReader(io.StringIO(text)))
    if len(rows)!=RAW_ROWS[name]:
        raise ValueError("Original 2024 established insect positive row number changed "+name)
    if not rows or any(None in r for r in rows[:20]):
        raise ValueError("Bad original source table shape")
    return rows,{
        "source_member":name,
        "sha256":sha(payload),
        "size_bytes":len(payload),
        "record_url":url,
        "source_columns":list(rows[0]),
        "source_rows":len(rows),
        "source_declared_md5":declared,
    }


def load_establishment_source():
    payload,final,ctype=get(ZENODO_META_URL)
    record=json.loads(payload)
    if record.get("id")!=ZENODO_RECORD:
        raise ValueError("Different original Zenodo 2024 dataset record")
    meta=record.get("metadata") or {}
    if str(meta.get("version"))!="0.2.1":
        raise ValueError("Published 2024 v0.2.1 establishment source version drifted")
    files=record.get("files") or []
    if isinstance(files,dict):
        files=list(files.values())
    byname={f.get("key") or f.get("filename"):f for f in files}
    if not set(DATA) <= set(byname):
        raise ValueError("Zenodo 2024 original 3 CSV file names changed")
    sources={}
    receipts={}
    for name,item in DATA.items():
        rows,receipt=read_zenodo_member(name,byname[name])
        sources[name]=rows
        receipts[name]=receipt
    return sources,{
        "record_id":record["id"],
        "doi":meta.get("doi"),
        "publication_date":meta.get("publication_date"),
        "original_version":meta.get("version"),
        "record_metadata_SHA256":sha(payload),
        "sources":receipts,
    }


def original239(path):
    raw=Path(path).read_bytes()
    with Path(path).open(encoding="utf-8-sig",newline="") as f:
        rows=list(csv.DictReader(f))
    if len(rows)!=239 or len({r["species"] for r in rows})!=239:
        raise ValueError("Not exactly 239 original accepted source butterflies")
    n=sum(int(r["introduced_added_units"]) for r in rows)
    expanded=sum(int(r["introduced_added_units"])>0 for r in rows)
    if n!=ORIG_RESOURCE_GAIN or expanded!=206:
        raise ValueError("Original accepted HOSTS/WCVP source 14553/206 margins changed")
    for r in rows:
        if not binomial(r["species"]):
            raise ValueError("Invalid original 239 butterfly taxon name")
    return {r["species"]:r for r in rows},{
        "original_239_metrics_sha256":sha(raw),
        "original_239_butterflies":len(rows),
        "original_added_host_resource_units":n,
        "original_expanded_butterflies":expanded,
    }


def analyze(zen,original):
    occurrence=zen["occurrence_table_20230317.csv"]
    attrs=zen["attribute_table_20230317.csv"]
    sources=zen["source_table_20230317.csv"]
    req={"genus_species","family","order","rank","year","region",
         "intentional_release","established_indoors_or_outdoors","eradicated"}
    if not req<=set(occurrence[0]) or \
       not {"genus_species","order","plant_feeding","total"}<=set(attrs[0]):
        raise ValueError("Original 2024 established alien insect source field grain changed")
    known=set(original)
    occ_by_taxon=defaultdict(list)
    positive_identified=0
    regime=Counter()
    all_order=Counter()
    all_rank=Counter()
    all_regional=Counter()
    years=Counter()
    warning=Counter()
    pair_seen=set()
    for r in occurrence:
        order=(r.get("order") or "").strip()
        rank=(r.get("rank") or "").strip()
        sp=binomial(r["genus_species"])
        area=r["region"].strip()
        all_order[order]+=1
        all_rank[rank]+=1
        all_regional[area]+=1
        yr=(r.get("year") or "").strip()
        if yr and yr!="NA":
            years["available"]+=1
        else:
            years["unknown"]+=1
        if order=="Lepidoptera" and rank=="species" and sp:
            positive_identified+=1
            if (sp,area) in pair_seen:
                warning["duplicate_lepidoptera_species_area"]+=1
            pair_seen.add((sp,area))
            occ_by_taxon[sp].append({
                "recipient_zone":area,
                "year_source_first_report":yr if yr not in ("NA","") else None,
                "source_intentional_release_flag":(r.get("intentional_release") or "").strip(),
                "source_established_indoor_outdoor_flag":(
                    r.get("established_indoors_or_outdoors") or "").strip(),
                "source_eradication_flag":(r.get("eradicated") or "").strip(),
                "source_family":(r.get("family") or "").strip(),
                "source_2024_row_means_establishment_positive":True,
            })
            regime[(r.get("intentional_release") or "").strip()]+=1
    unique_lepid_species=set(occ_by_taxon)
    attrs_by_name=defaultdict(list)
    for r in attrs:
        attrs_by_name[(r["genus_species"] or "").strip()].append(r)
    source_name_taxonomic_overlaps=[]
    for sp in sorted(known & unique_lepid_species):
        rows=occ_by_taxon[sp]
        originalr=original[sp]
        attribute=attrs_by_name.get(sp,[])
        zones=sorted({z["recipient_zone"] for z in rows})
        yrs=[int(z["year_source_first_report"]) for z in rows
             if z["year_source_first_report"] and
             z["year_source_first_report"].isdigit() and
             len(z["year_source_first_report"])==4]
        source_name_taxonomic_overlaps.append({
            "species":sp,
            "original_239_family_breadth":float(originalr["host_family_count"]),
            "original_239_host_species_breadth":int(originalr["resolved_host_species"]),
            "original_239_source_native_resource_units":int(originalr["native_resource_units"]),
            "original_239_source_added_WCVP_resource_units":int(originalr["introduced_added_units"]),
            "original_239_source_introduced_share_contemporary":
                float(originalr["introduced_share_of_contemporary"]),
            "source_documented_non_native_recipient_area_count":len(zones),
            "source_documented_non_native_recipient_areas":zones,
            "original_source_region_and_year_records":rows,
            "source_earliest_documented_establishment_year":min(yrs) if yrs else None,
            "source_attribute_taxon_match_count":len(attribute),
            "source_attribute_plant_feeding_label":(
                attribute[0].get("plant_feeding") if len(attribute)==1 else None
            ),
            "source_species_is_established_positive_NOT_general_arrival_success_rate":True,
        })
    by_recipient=defaultdict(set)
    for z in source_name_taxonomic_overlaps:
        for area in z["source_documented_non_native_recipient_areas"]:
            by_recipient[area].add(z["species"])
    species_with_time=sum(z["source_earliest_documented_establishment_year"] is not None
                          for z in source_name_taxonomic_overlaps)
    # Nonnative species list is positive-only. Any original 239 species
    # not in the occurrence sheet is not a legitimate "failed establishment".
    no_record=sorted(known-set(z["species"] for z in source_name_taxonomic_overlaps))
    source_status={
        "original_239_butterflies_matched_as_established_non_native_in_at_least_one_source_region":
            len(source_name_taxonomic_overlaps),
        "original_239_butterflies_not_in_source_positive_list":len(no_record),
        "positive_original239_butterflies_with_documented_year":species_with_time,
        "positive_original239_butterflies_without_establishment_year":
            len(source_name_taxonomic_overlaps)-species_with_time,
        "original_239_nonlisted_species_cannot_be_called_failed_or_untransported":True,
        "original_239_matched_positive_by_recipient":{
            k:len(v) for k,v in sorted(by_recipient.items())
        },
        "original_239_positive_species_non_native_area_count_distribution":dict(
            Counter(z["source_documented_non_native_recipient_area_count"]
                    for z in source_name_taxonomic_overlaps)
        ),
        "all_source_Lepidoptera_positive_species_unique":len(unique_lepid_species),
        "all_source_Lepidoptera_positive_species_zone_rows":positive_identified,
        "all_source_all_taxa_order_positive_rows":dict(all_order),
        "all_source_all_taxa_rank_positive_rows":dict(all_rank),
        "all_source_positive_region_counts_all_taxa":dict(all_regional),
        "all_source_occurrence_records_by_year_status":dict(years),
        "all_source_Lepidoptera_intentional_release_flag_values":dict(regime),
        "source_duplicate_status":dict(warning),
        "source_positive_original239_species":source_name_taxonomic_overlaps,
    }
    # Positive-only associations are descriptive. No use of absent list as zero.
    if len(source_name_taxonomic_overlaps)>=10:
        series=[
            (z["source_documented_non_native_recipient_area_count"],
             z["original_239_source_added_WCVP_resource_units"],
             z["original_239_family_breadth"])
            for z in source_name_taxonomic_overlaps
        ]
        def spearman(a,b):
            aa=np.array(a,dtype=float);bb=np.array(b,dtype=float)
            def rank(x):
                out=np.empty(len(x),float)
                for v in np.unique(x):
                    mask=x==v
                    ix=np.where(x<v)[0].size
                    out[mask]=(2*ix+mask.sum()+1)/2
                return out
            if len(set(a))<2 or len(set(b))<2:
                return None
            return float(np.corrcoef(rank(aa),rank(bb))[0,1])
        source_status["positive_established_species_only_descriptive"]={
            "rho_non_native_area_count_vs_global_host_resource_added_units":
                spearman([x[0] for x in series],[x[1] for x in series]),
            "rho_non_native_area_count_vs_global_host_family_breadth":
                spearman([x[0] for x in series],[x[2] for x in series]),
            "correlation_is_NOT_colonization_probability_effect":True,
            "original_sample_size_species":len(series),
        }
    else:
        source_status["positive_established_species_only_descriptive"]={
            "status":"INSUFFICIENT_MATCHED_ESTABLISHMENT_TAXA",
            "n":len(source_name_taxonomic_overlaps),
        }
    source_status["interpretive_limits"]=[
        "The original published positive-only nonnative establishment rows are not taxon×region arrival attempts or systematically sampled failed establishments.",
        "The absence of a focal butterfly from this 2024 species list is NOT a biologically established absence, immigration failure, host-resource obstacle or negative survey.",
        "The 2024 source includes indoor/intentional/recently eradicated statuses that require filtering and historical interpretation, not proof of current reproducing wild populations.",
        "First documented establishment years are not individual border interceptions and are not linked to a founding shipment of known host plants.",
        "A taxon-wide 14553 resource-gain index may be high far from a particular recipient region: global resource opportunity cannot substitute local habitat availability.",
        "The original 239 butterfly names and Zenodo taxonomic species names may contain synonyms or unresolved subspecies not rescued by an unsafe fuzzy join.",
    ]
    return source_status


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--original239",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    data,receipt=load_establishment_source()
    original,source=original239(args.original239)
    result=analyze(data,original)
    status={
        "schema":SCHEMA,
        "status":"REAL_2024_ALIEN_ESTABLISHMENT_POSITIVE_LIST_MATCHED_TO_ORIGINAL_239_NO_FAILURE_ABSENCES",
        "source_zenodo":receipt,
        "source_original239":source,
        "source_positive_matches":result,
        "computed_utc":datetime.now(timezone.utc).isoformat(),
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(status,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "status":status["status"],
        "source_total_lepidoptera_unique":
            result["all_source_Lepidoptera_positive_species_unique"],
        "established_non_native_in_at_least_one_area":
            result["original_239_butterflies_matched_as_established_non_native_in_at_least_one_source_region"],
        "original_239_matched_by_receiving_area":
            result["original_239_matched_positive_by_recipient"],
        "positive_species_details":result["source_positive_original239_species"],
        "positive_only_resource_contrast":
            result["positive_established_species_only_descriptive"],
    },ensure_ascii=False,sort_keys=True))


if __name__=="__main__":
    main()
