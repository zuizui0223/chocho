#!/usr/bin/env python3
"""Audit two previously discarded original HOSTS fields: Lab Rearing, Damage.

The contents of these literal source columns are inspected BEFORE attaching
any ecological semantics. This is data availability and quality profiling,
NOT wild larval breeding validation or independent source confirmation.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from collections import Counter,defaultdict
from pathlib import Path

SCHEMA="chocho_pinned_hosts_lab_rearing_damage_metadata_audit_20261027_v1"
ORIGINAL_239="data/exploratory/gbif_country_20261020/original_239_taxon_manifest.json"
FOCAL=set(["Pieris rapae","Thymelicus lineola","Pieris brassicae","Danaus plexippus",
           "Papilio demoleus","Papilio xuthus","Cacyreus marshalli","Chilades pandava"])
REQUIRED={"Lab Rearing","Damage","Location","Insect Genus","Insect Species",
          "Hostplant Genus","Hostplant Species"}
COUNTERS=("original_all_records","original_complete_species_rows","original_239_taxon_rows")


def value(raw):
    return str(raw or "").strip() or "<EMPTY>"


def load_focal(path):
    d=json.loads(Path(path).read_text())
    if "species_receipts" not in d or len(d["species_receipts"])!=239:
        raise ValueError("Not original accepted butterfly taxon manifest")
    names=[x["species"] for x in d["species_receipts"]]
    if len(set(names))!=239:
        raise ValueError("Duplicate original 239 source butterflies")
    if not (FOCAL-{"Cacyreus marshalli","Chilades pandava","Papilio xuthus"}) <=set(names):
        raise ValueError("Original butterfly original-239 focus changed")
    return set(names)


def inventory(path,original):
    totals=Counter()
    field_stats=defaultdict(lambda:{"Lab Rearing":Counter(),"Damage":Counter()})
    sampled=defaultdict(lambda:{"Lab Rearing":Counter(),"Damage":Counter()})
    meta={}
    with Path(path).open(encoding="utf-8-sig",newline="") as f:
        reader=csv.DictReader(f)
        if not REQUIRED<=set(reader.fieldnames or []):
            raise ValueError("Pinned original HOSTS complete field header changed")
        meta["columns"]=reader.fieldnames
        for row in reader:
            totals["original_all_records"]+=1
            complete=all((row.get(k) or "").strip() for k in
                        ("Insect Genus","Insect Species","Hostplant Genus","Hostplant Species"))
            if not complete:
                totals["source_rows_missing_full_species_binomials"]+=1
                continue
            totals["original_complete_species_rows"]+=1
            insect=row["Insect Genus"].strip()+" "+row["Insect Species"].strip()
            field_stats["all_lepidopteran_source_with_full_binomials"]["Lab Rearing"][
                value(row["Lab Rearing"])]+=1
            field_stats["all_lepidopteran_source_with_full_binomials"]["Damage"][
                value(row["Damage"])]+=1
            if insect in original:
                totals["original_239_taxon_rows"]+=1
                group="original_239_source_butterflies"
                for col in ("Lab Rearing","Damage"):
                    field_stats[group][col][value(row[col])]+=1
            if insect in FOCAL:
                totals["candidate_eight_raw_taxon_rows"]+=1
                for col in ("Lab Rearing","Damage"):
                    sampled[insect][col][value(row[col])]+=1
    if (totals["original_all_records"],totals["original_complete_species_rows"],
        totals["original_239_taxon_rows"])!=(140485,86826,4695):
        raise ValueError("Original HOSTS global or fixed original 239 row counts drifted")
    result={
        "schema":SCHEMA,
        "status":"SOURCE_FIELD_COMPLETENESS_AND_VALUE_INVENTORY_NO_LARVAL_HOST_USE_ASSERTION",
        "source":{
            "HOSTS_original_commit":"globalbioticinteractions/HOSTS@808e0b869f9ec1adf8efff87cf6a395adda103e0",
            "original_resource_SHA256":hashlib.sha256(Path(path).read_bytes()).hexdigest(),
            "original_239_taxon_manifest_SHA256":hashlib.sha256(
                Path(ORIGINAL_239).read_bytes()).hexdigest() if Path(ORIGINAL_239).exists() else None,
            "original_source_columns":meta["columns"],
        },
        "row_flow":dict(totals),
        "by_population":{
            group:{
                col:{
                    "all_value_frequencies":dict(sorted(colcount.items(),
                                                        key=lambda x:(-x[1],x[0]))),
                    "distinct_original_literal_values":len(colcount),
                    "original_rows_with_nonempty_field":sum(
                        n for val,n in colcount.items() if val!="<EMPTY>"),
                    "original_rows_empty":colcount["<EMPTY>"]
                }
                for col,colcount in cols.items()
            }
            for group,cols in field_stats.items()
        },
        "eight_invasion_candidates_source_stage_field_counts":{
            species:{
                col:{"original_raw_taxon_rows":sum(counts.values()),
                     "distinct_literal_values":len(counts),
                     "nonempty_rows":sum(n for v,n in counts.items() if v!="<EMPTY>"),
                     "top_12_literal_values":[{"source_value":k,"original_rows":v}
                         for k,v in counts.most_common(12)]}
                for col,counts in cols.items()
            }
            for species,cols in sorted(sampled.items())
        },
        "limits":[
            "The HOSTS Lab Rearing and Damage columns are source text metadata with potentially undocumented encoding; values are not accepted larval survival or development until original HOSTS schema defines them.",
            "Even non-empty Lab Rearing does not imply a free-ranging larva or a novel geographical host adoption.",
            "Original HOSTS has no per-row original DOI/study citation columns; separately curated papers may overlap source records.",
            "Off-panel butterflies are presented separately, not counted as original 239.",
            "All empty original values are missing evidence, not the absence of ecological interaction."
        ]
    }
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--hosts-resource",required=True,type=Path)
    p.add_argument("--239-manifest",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    args=p.parse_args()
    original=load_focal(args.__dict__["239_manifest"])
    result=inventory(args.hosts_resource,original)
    result["source"]["original_239_taxon_manifest_SHA256"]=hashlib.sha256(
        args.__dict__["239_manifest"].read_bytes()).hexdigest()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print(json.dumps({
        "source":result["status"],
        "source_counts":result["row_flow"],
        "239_field_counts":{
            field:{k:v for k,v in data.items() if k in (
                "original_rows_with_nonempty_field","original_rows_empty",
                "distinct_original_literal_values")}
            for field,data in result["by_population"]["original_239_source_butterflies"].items()
        },
        "239_literal_values":{
            field:data["all_value_frequencies"]
            for field,data in result["by_population"]["original_239_source_butterflies"].items()
        },
        "eight_candidates":{
            sp:{col:record["nonempty_rows"] for col,record in val.items()}
            for sp,val in result["eight_invasion_candidates_source_stage_field_counts"].items()
        }
    },sort_keys=True))


if __name__=="__main__":
    main()
