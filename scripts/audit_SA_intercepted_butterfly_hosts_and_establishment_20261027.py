#!/usr/bin/env python3
"""Actual butterfly taxa among Lepidoptera in South African plant import samples.

All original Figshare CSV and author-only WCVP/HOSTS files are SHA-audited.
A butterfly identified on a crop shipment is an intercepted insect-commodity
ASSOCIATION, not independently verified larval feeding or future establishment.
"""
from __future__ import annotations
import argparse,json,re
from collections import Counter,defaultdict
from pathlib import Path

from analyze_sa_import_lepidoptera_detection_by_commodity_20261027 import (
    load_primary_file,commodity_class,lepidoptera_result
)
from analyze_hosts_location_wcvp_origin_20261027 import (
    accepted_pairs, read_rows, sha256,
)

SCHEMA="chocho_SA_original_intercepted_butterflies_host_commodity_pair_20261027_v01"
EXACT_SOURCE_BUTTERFLIES={"Lampides boeticus","Pieris brassicae"}
SOURCE_2021_META_LATER_ESTABLISHED_MOTH="Spodoptera frugiperda"


def parsed_species(label):
    s=str(label or "").strip()
    if "cf." in s or "aff." in s or "?" in s:
        return None
    m=re.fullmatch(r"([A-Z][a-z]+ [a-z][a-z-]+)(?: \([A-Z][A-Za-z\- ]+\))?",s)
    if not m:return None
    if m.group(1).split()[0] in ("Lepidoptera","Pyralidae","Tortricidae"):
        return None
    return m.group(1)


def explicit_stage_label(label):
    t=str(label).lower()
    if "larva" in t or "caterpillar" in t:
        return "EXPLICIT_LARVAL_STAGE"
    if "egg" in t or "ovipos" in t:
        return "EXPLICIT_EGG_OR_OVIPOSITION"
    if "pupa" in t or "chrysalis" in t:
        return "EXPLICIT_PUPAL_STAGE"
    if "adult" in t or "imago" in t:
        return "EXPLICIT_ADULT_STAGE"
    return "LIFE_STAGE_UNSPECIFIED"


def source_host_names(accepted_file,original239_file):
    original=read_rows(original239_file)
    species={r["species"] for r in original}
    if len(original)!=239 or len(species)!=239 or \
       not EXACT_SOURCE_BUTTERFLIES <= species:
        raise ValueError("Original 239-butterfly source taxonomy changed")
    pairs,quality,_=accepted_pairs(accepted_file,species)
    if quality["source_ambiguous_accepted_id_names"]:
        raise ValueError("Original accepted plant pairs ambiguous")
    known=defaultdict(set)
    for (insect,plant),accepted_id in pairs.items():
        if insect in EXACT_SOURCE_BUTTERFLIES:
            known[insect].add(plant)
    return known,quality


def interpret(rows,metadata,known):
    m=[r for r in metadata if (r.get("scientificName") or "").strip()==
       SOURCE_2021_META_LATER_ESTABLISHED_MOTH and
       (r.get("order") or "").strip()=="Lepidoptera"]
    if len(m)!=1:
        raise ValueError("Original Spodoptera status metadata source changed")
    moth=m[0]
    if (moth["SA.occurenceStatus.historical"].strip().lower()!="absent" or
        not moth["SA.occurenceStatus.current"].strip().lower().startswith("present")):
        raise ValueError("Source Spodoptera historical/modern status changed")
    records=defaultdict(Counter)
    stage=defaultdict(Counter)
    stage_by_commodity=defaultdict(Counter)
    stage_by_material=defaultdict(Counter)
    originals=defaultdict(lambda:{"commodity":Counter(),"year":Counter(),
                                  "crop":Counter(),"regime":Counter(),
                                  "stage":Counter(),"exact_crop_host_yes":0,
                                  "exact_crop_binomial_unknown":0})
    for sample in rows:
        value=(sample.get("Insecta: Lepidoptera") or "").strip()
        if lepidoptera_result(value)!="DETECTED_LEPIDOPTERA":
            continue
        commodity=(sample.get("crop.commodity") or "").strip() or "<EMPTY>"
        group=commodity_class(commodity)
        year=(sample.get("year") or "").strip() or "<MISSING>"
        regime=(sample.get("sample.source") or "").strip() or "<MISSING>"
        cp=(sample.get("crop.scientificName") or "").strip()
        tax=parsed_species(value)
        group_type="ORIGINAL_239_EXACT_BUTTERFLY" if tax in EXACT_SOURCE_BUTTERFLIES else (
            "SOURCE_MOTH_SPODOPTERA" if tax==SOURCE_2021_META_LATER_ESTABLISHED_MOTH else
            "OTHER_SOURCE_TAXON_OR_UNRESOLVED"
        )
        records[group_type][commodity]+=1
        stage[group_type][explicit_stage_label(value)]+=1
        stage_by_commodity[commodity][explicit_stage_label(value)]+=1
        stage_by_material[group][explicit_stage_label(value)]+=1
        if tax in EXACT_SOURCE_BUTTERFLIES or tax==SOURCE_2021_META_LATER_ESTABLISHED_MOTH:
            d=originals[tax]
            d["commodity"][commodity]+=1
            d["year"][year]+=1
            d["crop"][cp or "<MISSING_CROP_NAME>"]+=1
            d["regime"][regime]+=1
            d["stage"][explicit_stage_label(value)]+=1
            if tax in EXACT_SOURCE_BUTTERFLIES:
                if cp in known[tax]:
                    d["exact_crop_host_yes"]+=1
                elif not cp or len(cp.split())<2:
                    d["exact_crop_binomial_unknown"]+=1
    if set(originals)&EXACT_SOURCE_BUTTERFLIES!=EXACT_SOURCE_BUTTERFLIES:
        raise ValueError("External intercepted exact-source GEB butterfly set changed")
    if sum(sum(v.values()) for v in records.values())!=178:
        raise ValueError("Original positive Lepidoptera source cell count changed")
    named=[]
    for tax in sorted(originals):
        d=originals[tax]
        is_butterfly=tax in EXACT_SOURCE_BUTTERFLIES
        named.append({
            "identified_lepidoptera_name":tax,
            "in_original_GEB_239_butterflies":is_butterfly,
            "positive_original_plant_sample_rows":sum(d["commodity"].values()),
            "original_imported_commodity_type_counts":dict(d["commodity"]),
            "original_imported_commodity_class_counts":dict(Counter({
                kind:sum(v for commodity,v in d["commodity"].items()
                         if commodity_class(commodity)==kind)
                for kind in ("LIVE_PROPAGATION","SEED_GRAIN",
                             "FRESH_OR_DRY_PRODUCE","OTHER_OR_UNRESOLVED")
            })),
            "original_source_year_counts":dict(d["year"]),
            "original_sampling_regime_counts":dict(d["regime"]),
            "original_crop_scientific_name_counts":dict(d["crop"]),
            "stage_metadata_source_annotation_counts":dict(d["stage"]),
            "original_crop_binomial_exact_to_original_known_HOSTS_plant_pair":(
                d["exact_crop_host_yes"] if is_butterfly else None),
            "crop_taxonomy_unknown_or_genus_only":(
                d["exact_crop_binomial_unknown"] if is_butterfly else None),
            "recording_scope":"intercepted_order_insect_positive_specimen_not_founder_or_wild_larval_feeding",
        })
    sp=next((z for z in named if z["identified_lepidoptera_name"]==
            SOURCE_2021_META_LATER_ESTABLISHED_MOTH),None)
    if sp is not None:
        source_year=[int(y) for y in sp["original_source_year_counts"] if y.isdigit()]
        before2017=sum(n for yr,n in sp["original_source_year_counts"].items()
                       if yr.isdigit() and int(yr)<2017)
        sp["source_positive_samples_pre_2017"]=before2017
        sp["source_earliest_positive_year"]=min(source_year) if source_year else None
        sp["source_2017_year_in_current_status_is_metadata_not_specimen_date"]=True
    return {
        "all_positive_Lepidoptera_commodity_types_by_source_taxon_category":{
            typ:dict(items) for typ,items in records.items()},
        "all_positive_Lepidoptera_explicit_stages_by_source_taxon_category":{
            typ:dict(items) for typ,items in stage.items()},
        "all_positive_Lepidoptera_stages_by_original_commodity":{
            kind:dict(stages) for kind,stages in stage_by_commodity.items()},
        "all_positive_Lepidoptera_stages_by_material_class":{
            kind:dict(stages) for kind,stages in stage_by_material.items()},
        "intercepted_exact_original_239_butterfly_names":
            sorted(EXACT_SOURCE_BUTTERFLIES),
        "intercepted_species_and_metadata_recent_establishment":named,
        "historical_absent_to_present_moth":moth["scientificName"],
        "moth_SA_historical_status":moth["SA.occurenceStatus.historical"],
        "moth_SA_current_status":moth["SA.occurenceStatus.current"],
        "moth_has_unambiguous_original_sample_positive":sp is not None,
        "original_SA_taxon_data_without_link_to_specific_established_descendants":True,
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--original239",type=Path,required=True)
    p.add_argument("--original-host-source",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    original=json.loads((args.original_host_source/"SOURCE_MANIFEST.json").read_text())
    if original.get("inputs_count")!=11:
        raise ValueError("Wrong original author-only 239 WCVP source")
    for name,digest in original["sha256"].items():
        if sha256(args.original_host_source/name)!=digest:
            raise ValueError("Original accepted HOSTS source member digest changed")
    known,source_qc=source_host_names(
        args.original_host_source/"inputs/original/insect_host_accepted.csv",
        args.original239)
    samples,receipt=load_primary_file("sample")
    taxonomy,tax_receipt=load_primary_file("metadata")
    r=interpret(samples,taxonomy,known)
    report={
        "schema":SCHEMA,
        "status":"EXACT_TWO_ORIGINAL_BUTTERFLIES_INTERCEPTED_SOURCE_POSITIVE_BY_COMMODITY_NOT_ESTABLISHMENT",
        "source":{
            "original_public_figshare_DOI":"10.6084/m9.figshare.c.5180681",
            "original_SA_samples_sha256":receipt["member_SHA256"],
            "original_SA_taxon_metadata_sha256":tax_receipt["member_SHA256"],
            "original_author_only_butterfly_HOSTS_source_sha256":
                sha256(args.original_host_source/"inputs/original/insect_host_accepted.csv"),
            "original_239_source":sha256(args.original239),
            "original_WCVP_accepted_join_qc":source_qc,
        },
        "data":r,
        "inferential_gates":[
            "Two precisely identified original 239 butterflies are too few for a host-vector-to-establishment coefficient.",
            "A detected butterfly taxon on the cargo is not evidence that the cargo plant is its larval host.",
            "A matched original HOSTS binomial crop/plant link does not show observed larval feeding during this interception.",
            "The establishment of Spodoptera frugiperda in South Africa is a species-wide historical metadata status, not traced ancestry from any intercepted specimen.",
            "A 178-positive Lepidoptera sample may include moths or higher larval labels, not necessarily butterflies.",
            "No original full importer/exporter identities, exact sample catalogue IDs or private localities appear in this output."
        ],
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "status":report["status"],
        "source_original_239_exact_butterfly_positive_samples":[
            r for r in report["data"]["intercepted_species_and_metadata_recent_establishment"]
            if r["in_original_GEB_239_butterflies"]
        ],
        "source_moth_historical_establishment_case":[
            r for r in report["data"]["intercepted_species_and_metadata_recent_establishment"]
            if r["identified_lepidoptera_name"]==SOURCE_2021_META_LATER_ESTABLISHED_MOTH
        ],
        "positive_taxonomic_stage_distribution":
            r["all_positive_Lepidoptera_explicit_stages_by_source_taxon_category"],
        "positive_stage_by_commodity":r["all_positive_Lepidoptera_stages_by_original_commodity"],
        "positive_stage_by_material":r["all_positive_Lepidoptera_stages_by_material_class"]
    },sort_keys=True))


if __name__=="__main__":
    main()
