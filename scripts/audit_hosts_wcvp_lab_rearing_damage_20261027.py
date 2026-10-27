#!/usr/bin/env python3
"""Read original HOSTS Lab Rearing / Damage text at WCVP report-area grain.

Missing metadata does NOT mean butterfly larvae failed to develop.
Source Y is not direct field host use; Damage contains adult/egg labels too.
No original butterfly x plant x locality link is publicly exported.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

from analyze_hosts_location_wcvp_origin_20261027 import (
    SOURCE_ZIP_SHA, accepted_pairs, classify_location, load_country_atoms,
    read_rows, sha256, wcvp_range, STATUS as SOURCE_ORIGINS
)

SCHEMA="chocho_hosts_location_wcvp_239_lab_rearing_damage_source_audit_20261027_v1"
ORIGINAL_SOURCE_ROW_COUNTERS=(140485,86826,4695,3898,3586)


def collect(hosts, accepted, focal):
    totals=Counter()
    triples=defaultdict(lambda:{
        "lab":set(),"damage":set(),"rows":0
    })
    for r in read_rows(hosts):
        totals["all_host_source_rows"]+=1
        fields=("Insect Genus","Insect Species","Hostplant Genus","Hostplant Species")
        if not all(str(r.get(f) or "").strip() for f in fields):
            continue
        totals["full_species_source_rows"]+=1
        butterfly=r["Insect Genus"].strip()+" "+r["Insect Species"].strip()
        if butterfly not in focal:
            continue
        totals["focal_239_source_rows"]+=1
        host=r["Hostplant Genus"].strip()+" "+r["Hostplant Species"].strip()
        accepted_id=accepted.get((butterfly,host))
        if accepted_id is None:
            totals["missing_frozen_accepted_pair"]+=1
            continue
        totals["accepted_focal_source_rows"]+=1
        label=(r.get("Location") or "").strip()
        lab=(r.get("Lab Rearing") or "").strip() or "<EMPTY>"
        damage=(r.get("Damage") or "").strip() or "<EMPTY>"
        if lab not in {"<EMPTY>","Y","?"}:
            raise ValueError("Unknown original HOSTS laboratory rearing encoding")
        record=triples[butterfly,accepted_id,label]
        record["lab"].add(lab)
        record["damage"].add(damage)
        record["rows"]+=1
    if (totals["all_host_source_rows"], totals["full_species_source_rows"],
        totals["focal_239_source_rows"], totals["accepted_focal_source_rows"],
        len(triples))!=ORIGINAL_SOURCE_ROW_COUNTERS:
        raise ValueError("Original 140485/86826/4695/3898/3586 source identity drift")
    return triples,dict(totals)


def summarize(triples, statuses):
    counts=defaultdict(Counter)
    damage_text=defaultdict(Counter)
    all_lab_values=Counter()
    for key, info in triples.items():
        status=statuses.get(key)
        if status not in SOURCE_ORIGINS:
            raise ValueError("Missing expected source botanical origin category")
        q=counts[status]
        q["source_distinct_butterfly_host_Location_triples"]+=1
        q["underlying_original_HOSTS_rows"]+=info["rows"]
        if "Y" in info["lab"]:
            q["triple_contains_any_Lab_Rearing_Y"]+=1
        if "?" in info["lab"]:
            q["triple_contains_any_Lab_Rearing_question_mark"]+=1
        if info["lab"]=={"<EMPTY>"}:
            q["triple_all_Lab_Rearing_blank"]+=1
        meaningful={s for s in info["damage"] if s!="<EMPTY>"}
        if meaningful:
            q["triple_contains_any_nonempty_Damage"]+=1
            for label in meaningful:
                damage_text[status][label]+=1
        if len(info["lab"])>1:
            q["triple_contains_multiple_rearing_literal_values"]+=1
        if len(meaningful)>1:
            q["triple_contains_multiple_damage_literal_values"]+=1
        all_lab_values.update(info["lab"])
    return {
        "counts_by_source_origin":{
            status:dict(c) for status,c in sorted(counts.items())
        },
        "introduced_area_top_damage_literal_values":[
            {"original_literal_damage":label,"distinct_triples":n}
            for label,n in damage_text["INTRODUCED_CONFIRMED"].most_common(20)
        ],
        "lab_literal_values_across_source_triples":dict(all_lab_values),
        "caution":"Y/?/Damage are original HOSTS metadata, not independently validated larval-rearing outcomes, field surveys or ecological absences."
    }


def run(args):
    source=args.author_source
    manifest=json.loads((source/"SOURCE_MANIFEST.json").read_text())
    if manifest.get("inputs_count")!=11:
        raise ValueError("Not original 11-member WCVP/HOSTS offline source")
    for name, digest in manifest["sha256"].items():
        if sha256(source/name)!=digest:
            raise ValueError("Original author-only WCVP/HOSTS member SHA mismatch")
    root=source/"inputs"
    focal_rows=read_rows(root/"anthropogenic_species_metrics.csv")
    focal={x["species"] for x in focal_rows}
    if len(focal)!=239 or len(focal_rows)!=239:
        raise ValueError("Original butterfly taxon universe has changed")
    accepted,meta,_=accepted_pairs(root/"original/insect_host_accepted.csv",focal)
    if meta["source_ambiguous_accepted_id_names"]:
        raise ValueError("Original botanical accepted host pair ambiguity")
    native,current=wcvp_range(
        root/"original/native_extant_nondoubtful_wgsrpd3.csv",
        root/"original/extant_nondoubtful_wgsrpd3.csv")
    atoms, atom_metadata, regions, geo = load_country_atoms(
        args.country_atoms, args.country_atoms_manifest, args.country_summary)
    triples,flow=collect(args.hosts_resource,accepted,focal)
    result={}
    for cutoff,label in ((.95,"0.95"),(.80,"0.8")):
        status_cache={}
        by_status={}
        for (sp,pid,loc) in triples:
            k=pid,loc
            if k not in status_cache:
                status_cache[k]=classify_location(
                    pid,loc,threshold=cutoff,
                    atoms=atoms,atom_metadata=atom_metadata,
                    original_regions=regions,native=native,current=current)
            by_status[sp,pid,loc]=status_cache[k]
        report=summarize(triples,by_status)
        original={
            "0.95":{"NATIVE_CONFIRMED":399,
                    "INTRODUCED_CONFIRMED":179,"WCVP_UNLISTED":303},
            "0.8":{"NATIVE_CONFIRMED":424,
                   "INTRODUCED_CONFIRMED":197,"WCVP_UNLISTED":475}
        }[label]
        for x, n in original.items():
            if report["counts_by_source_origin"].get(x,{}).get(
                "source_distinct_butterfly_host_Location_triples",0)!=n:
                raise ValueError("Underlying original documented origin triple count changed")
        if sum(v["source_distinct_butterfly_host_Location_triples"] for v
               in report["counts_by_source_origin"].values())!=3586:
            raise ValueError("Source accepted triple count conservation failed")
        result[label]=report
    out={
        "schema":SCHEMA,
        "status":"ORIGINAL_HOSTS_LAB_REARING_DAMAGE_LITERAL_METADATA_ONLY",
        "source":{
            "original_source_zip_sha256":SOURCE_ZIP_SHA,
            "original_HOSTS_sha256":sha256(args.hosts_resource),
            "original_239_butterfly_species":239,
            "original_accepted_host_join":meta,
            "geographic_source":geo,
            "original_accepted_pairs_sha256":sha256(root/"original/insect_host_accepted.csv"),
        },
        "source_row_flow":flow,
        "source95_and_sensitivity80":result,
        "inferential_stop":[
            "Lab Rearing 'Y' is a source-encoded method flag with no documented proof of local wild adult emergence.",
            "Lab Rearing '?' and blanks are unknown; do not call them larval failure.",
            "Damage values such as on leaves, adult association or oviposition only do not demonstrate completed development.",
            "An original HOSTS Location native/introduced source category classifies PLANT region, not butterfly nativity.",
            "These flags are not independent external data and cannot validate de novo host adoption or colonization.",
            "No butterfly-plant-source Location pair names or original copyrighted rows are exported."
        ]
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    headline=out["source95_and_sensitivity80"]["0.95"]["counts_by_source_origin"]
    print(json.dumps({
        "status":out["status"],
        "source_total":out["source_row_flow"],
        "area95_intro_rearing_damage":headline.get("INTRODUCED_CONFIRMED"),
        "area95_native_rearing_damage":headline.get("NATIVE_CONFIRMED"),
        "area95_unlisted_rearing_damage":headline.get("WCVP_UNLISTED"),
        "area80_intro_rearing_damage":
            result["0.8"]["counts_by_source_origin"].get("INTRODUCED_CONFIRMED"),
    },sort_keys=True))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--hosts-resource",type=Path,required=True)
    p.add_argument("--author-source",type=Path,required=True)
    p.add_argument("--country-atoms",type=Path,required=True)
    p.add_argument("--country-atoms-manifest",type=Path,required=True)
    p.add_argument("--country-summary",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    run(p.parse_args())
