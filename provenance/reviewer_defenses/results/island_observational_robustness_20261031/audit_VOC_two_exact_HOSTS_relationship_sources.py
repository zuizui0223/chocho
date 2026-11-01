#!/usr/bin/env python3
"""Check public original HOSTS evidence for two exact VOС/GloNAF source triads.

No species is reclassified based on the label "HOSTS". We report exact
original interaction rows and whether the source provides life-stage and
bibliographic source annotations. This check does NOT certify a butterfly
actually ate the host locally.
"""
from __future__ import annotations

import csv
from collections import Counter,defaultdict
import hashlib
import json
from pathlib import Path

BASE=Path(__file__).resolve().parent
SOURCE=Path("_external/HOSTS/resource.csv")
OUT=BASE/"voc_cargo_two_exact_HOSTS_relationship_provenance_v1.json"
TARGET={
    ("zizina","otis","indigofera","zollingeriana"),
    ("pieris","rapae","nicotiana","tabacum"),
}


def classify_rows(rows):
    columns=set()
    selected=defaultdict(list)
    n=0
    for row in rows:
        columns.update(row)
        n+=1
        taxon=(str(row.get("Insect Genus") or "").strip().casefold(),
               str(row.get("Insect Species") or "").strip().casefold(),
               str(row.get("Hostplant Genus") or "").strip().casefold(),
               str(row.get("Hostplant Species") or "").strip().casefold())
        if taxon in TARGET:
            # Only two publicly documented source records are exported.
            # Keep annotations for data-quality review, with bounded length.
            meta={k:str(v or "")[:220] for k,v in row.items()
                  if k not in ("Insect Genus","Insect Species",
                               "Hostplant Genus","Hostplant Species") and v}
            selected[taxon].append(meta)
    return {
        "source_rows":n,
        "all_source_field_names":sorted(columns),
        "source_provenance_field_names":[k for k in sorted(columns)
            if any(t in k.casefold() for t in
                   ("reference","source","stage","life","interaction","note","paper","cite","record"))],
        "matched":{
            " ".join(k):{"original_exact_match_count":len(selected.get(k,[])),
                        "source_annotations":selected.get(k,[])}
            for k in sorted(TARGET)
        },
    }


def run():
    with SOURCE.open(encoding="utf-8-sig",newline="") as f:
        x=classify_rows(csv.DictReader(f))
    if x["source_rows"]!=140485:
        raise RuntimeError("Frozen HOSTS row count drifted; stop provenance inference")
    if any(x["matched"][" ".join(k)]["original_exact_match_count"]<1 for k in TARGET):
        raise RuntimeError("One of the source HOSTS exact triplets was not found")
    h=hashlib.sha256()
    with SOURCE.open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    receipt={
      "schema":"chocho_17c_VOC_two_exact_modern_HOSTS_record_data_quality_v1",
      "status":"ORIGINAL_HOSTS_ROW_PROVENANCE_NOT_LOCAL_LARVAL_USE_CONFIRMATION",
      "pinned_HOSTS_original_resource_sha256":h.hexdigest(),
      **x,
      "limits":[
        "The original HOSTS interaction may include copied/incorrect or unusual ecological associations; independent stage and biological validity still require external primary literature.",
        "A 17th-century commodity genus and this host source row cannot establish the exact plant species was live cargo on an original island route.",
        "Indexed first butterfly records and GloNAF current naturalized plant species do not establish local butterfly larval feeding or establishment.",
      ],
    }
    OUT.write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print(json.dumps(receipt,ensure_ascii=False,indent=2))


if __name__=="__main__":
    run()
