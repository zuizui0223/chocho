#!/usr/bin/env python3
"""Quality-stratify mapped monarch resources using independent 127-host review.

SOURCE: Greenstein, Steele & Taylor 2022 PLOS ONE DOI:
10.1371/journal.pone.0269701, supplementary S1 CC BY.
STATUS: post-hoc independent host-quality evidence pilot for Danaus plexippus,
NOT cross-taxon fitness evidence or measured recruitment. The PLOS table must
reproduce 34 H, 42 L, 33 N, 18 U before any map is calculated. Exact taxon
matches only; unrecognized and synonymous names are NOT silently declared poor.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

from audit_external_monarch_host_quality_appendix import (
    tables_in_docx, code_in_row,
)

PLOS_SHA="bc19e33823a3e860443b0086af129e39f4218663233fdceeccb6b08a15b3af79"
PLOS_DOI="10.1371/journal.pone.0269701.s001"
SPECIES="Danaus plexippus"
OBSERVED=(40,181,263,82)
LATIN=re.compile(r"\b([A-Z][a-z]{2,}\s+[a-z][a-z.-]+)\b")
EXPECTED={"high":34,"low":42,"non_host":33,"unsubstantiated":18}


def identity(s:str)->str:
    s=s.strip().replace("×","x")
    mt=LATIN.search(s)
    return mt.group(1).casefold() if mt else ""


def read_csv(path:Path)->list[dict]:
    with path.open(newline="",encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def source_quality(docx:Path)->tuple[dict[str,str],dict]:
    sha=hashlib.sha256(docx.read_bytes()).hexdigest()
    if sha!=PLOS_SHA:
        raise ValueError("PLOS S1 publication supplementary bytes changed")
    choices=[]
    for idx,rows in enumerate(tables_in_docx(docx)):
        records=[]
        for row in rows:
            cat=code_in_row(row)
            if not cat:
                continue
            # Published plant-name column is expected first, but allow fallback
            # from the remaining cells only when the row has ONE binomial.
            first=identity(row[0] if row else "")
            all_candidates={identity(v) for v in row}
            all_candidates.discard("")
            if first:
                plant=first
            elif len(all_candidates)==1:
                plant=next(iter(all_candidates))
            else:
                raise ValueError(f"plant name ambiguous in source appendix table {idx}")
            records.append((plant,cat))
        if len(records)==127:
            choices.append((idx,records))
    if len(choices)!=1:
        raise ValueError("not exactly one 127-record host-quality table in supplement")
    idx,rows=choices[0]
    labels={}
    for name,cat in rows:
        if not name or name in labels:
            raise ValueError("duplicate/blank taxon in 127-host source; review synonym manually")
        labels[name]=cat
    grouped=Counter(
        "high" if v.startswith("H")
        else "low" if v.startswith("L")
        else "non_host" if v=="N" else "unsubstantiated"
        for v in labels.values()
    )
    if dict(grouped)!=EXPECTED:
        raise ValueError(f"published 34/42/33/18 classes drift: {dict(grouped)}")
    return labels,{
        "doi":PLOS_DOI,
        "source_sha256":sha,
        "host_classification_rows":len(rows),
        "table_index_zero_based":idx,
        "reported_counts":dict(grouped),
        "source_type":"Independent published performance-class compilation; global transfer is NOT validated",
        "license":"CC BY; cite Greenstein et al. (2022)",
    }


def geography(path:Path)->dict[str,set[str]]:
    d=defaultdict(set)
    for row in read_csv(path):
        key=row["accepted_plant_name_id"].strip()
        unit=row["area_code_l3"].strip()
        if key and unit:
            d[key].add(unit)
    return dict(d)


def sum_envelope(hosts:set[str],geo:dict[str,set[str]])->set[str]:
    union=set()
    for host in hosts:
        union.update(geo.get(host,set()))
    return union


def reconstruct_quality_pilot(
    quality:dict[str,str],
    host_links:list[dict],
    native:dict[str,set[str]],
    contemp:dict[str,set[str]],
) -> dict:
    # Fixed study focal species, not all monarchs or all butterflies.
    this=[x for x in host_links if x["insect_species"].strip()==SPECIES]
    all_ids={x["accepted_plant_name_id"].strip() for x in this}
    if len(all_ids)!=OBSERVED[0]:
        raise ValueError(f"frozen Danaus plexippus accepted-host count drift {len(all_ids)}")
    ids={}
    names={}
    matched_by_status=defaultdict(set)
    unmatched=[]
    for row in this:
        hid=row["accepted_plant_name_id"].strip()
        if not hid:
            raise ValueError("empty accepted plant ID")
        preferred=identity(row.get("accepted_name",""))
        original=identity(row.get("input_host_name",""))
        vals={quality[k] for k in (preferred,original) if k in quality}
        if len(vals)>1:
            raise ValueError("conflicting quality classes across accepted/original botanical names")
        class_status=next(iter(vals)) if vals else None
        if hid in ids and ids[hid]!=class_status:
            raise ValueError("ambiguous status for same accepted host ID")
        is_first_row=hid not in ids
        ids[hid]=class_status
        names[hid]=row.get("accepted_name","") or row.get("input_host_name","")
        if class_status is None and is_first_row:
            unmatched.append({"accepted_host_id":hid,"plant_name":names[hid]})
        elif class_status is not None:
            matched_by_status[class_status].add(hid)
    main_native=sum_envelope(all_ids,native)
    main_cont=sum_envelope(all_ids,contemp)
    main_added=main_cont-main_native
    if (len(all_ids),len(main_native),len(main_cont),len(main_added))!=OBSERVED:
        raise ValueError("original Danaus resource map did NOT replicate 40/181/263/82")

    def stats(class_set:set[str])->dict:
        n=sum_envelope(class_set,native)
        c=sum_envelope(class_set,contemp)
        add=c-n
        introduced_host_map=set()
        for hid in class_set:
            introduced_host_map.update(contemp.get(hid,set())-native.get(hid,set()))
        support=main_added & introduced_host_map
        return {
            "accepted_host_species":len(class_set),
            "native_regions_under_quality_filter":len(n),
            "contemporary_regions_under_quality_filter":len(c),
            "new_regions_under_quality_filter":len(add),
            "fraction_of_original_82_still_classified_new":len(main_added & add)/82,
            "original_82_cells_supported_by_any_source_quality_class_introduced_host":len(support),
            "source_quality_supported_fraction_of_original_82":len(support)/82,
            "note":"Presence and source quality are not region-specific actual larval survival",
        }
    high={h for h,status in ids.items() if status and status.startswith("H")}
    strict_high={h for h,status in ids.items() if status=="H3"}
    low={h for h,status in ids.items() if status and status.startswith("L")}
    nonhost={h for h,status in ids.items() if status=="N"}
    unstated={h for h,status in ids.items() if status=="U"}
    unmatched_ids=set(all_ids)-set().union(high,low,nonhost,unstated)
    if len(unmatched_ids)!=len(unmatched):
        raise ValueError("unmatched accepted host count mismatch")

    sets={
        "all_original_HOSTS":all_ids,
        "literature_high_H1_H2_H3":high,
        "literature_confirmed_high_H3_only":strict_high,
        "literature_low_L1_L2_L3":low,
        "literature_nonhosts_N":nonhost,
        "literature_unsubstantiated_U":unstated,
        "not_exactly_matched_in_published_127_host_review":unmatched_ids,
        "quality_high_or_low":high|low,
    }
    output={
        "schema":"chocho_danaus_plexippus_larval_host_quality_geography_v0.1",
        "status":"EXPLORATORY_SOURCE_QUALITY_POTENTIAL_GEOGRAPHY_ONLY",
        "focal_butterfly":SPECIES,
        "original_frozen_geography":{
            "accepted_HOSTS_species":40,
            "native_regions":181,
            "contemporary_regions":263,
            "distinct_added_regions":82,
        },
        "classes":{
            key:stats(values) for key,values in sets.items()
        },
        "host_class_counts":dict(Counter(
            status or "NO_EXACT_NAME_MATCH" for status in ids.values()
        )),
        "taxa_requiring_manual_synonym_review":unmatched,
        "support_rule":{
            "within_original_82_cells":
                "At least one introduced HOSTS host with an independently documented source-quality class; not demonstration of functional larval use in that region",
            "quality_filtered_new_regions":
                "Recompute BOTH native and contemporary unions after restricting host list; denominator can change and new-region identities may be reclassified",
        },
        "prohibited_claims":[
            "PLOS inferred high-performance host in New Orleans or other study setting != globally viable larval host in every WCVP plant range cell",
            "Excluding unmatched or unsubstantiated hosts from a conditional high-quality map does not establish they are unsuitable",
            "Source quality classes reflect 127 candidate monarch host taxa only and cannot be exported to the other 238 focal butterfly species",
            "No oviposition diversion, breeding, reproductive failure, source-sink fitness or new butterfly colonization was directly observed",
        ],
    }
    return output


def main():
    p=argparse.ArgumentParser()
    for key in ("plos_s1_docx","insect_host_csv","native_distribution_csv",
                "contemporary_distribution_csv","output_json"):
        p.add_argument("--"+key.replace("_","-"),type=Path,required=True)
    a=p.parse_args()
    quality,source=source_quality(a.plos_s1_docx)
    report=reconstruct_quality_pilot(
        quality,read_csv(a.insect_host_csv),
        geography(a.native_distribution_csv),
        geography(a.contemporary_distribution_csv),
    )
    report["external_performance_source"]=source
    a.output_json.parent.mkdir(parents=True,exist_ok=True)
    a.output_json.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in report.items()
        if k not in ("taxa_requiring_manual_synonym_review",)},indent=2))


if __name__=="__main__":
    main()
