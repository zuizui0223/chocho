#!/usr/bin/env python3
"""Source-locked PUBLIC label/field coverage only; NO phenology effect inspection."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from collections import Counter
from pathlib import Path
from tempfile import TemporaryDirectory
from zipfile import ZipFile

from audit_monarch_annotated_host_source_20261028 import source_download, SHARE_URL

ORIGINAL_ZIP_SHA256="4e4def63669b34738bcbfb8ba32ccf6c61d02cc6992b8f1ebfa8bc350a3987db"
CORE="v1/v1_core_larval_host_records.csv"
THINNED="v1/v1_thinned_larval_host_records.csv"
FULL="v1/v1_full_larval_host_records.csv"
FIELDS=("host","confidence","setting","cultivated","rearing","native_status","EW")
MISSING={"", "na", "n/a", "nan", "none", "null", "unknown", "?", "[empty]"}


def count_file(archive, member):
    distributions={x:Counter() for x in FIELDS}
    distinct={x:set() for x in ("host","hex_id","year")}
    missing=Counter()
    format_audit=Counter()
    host_names=Counter()
    n=0
    with archive.open(member) as stream:
        csvr=csv.DictReader(io.TextIOWrapper(stream,encoding="utf-8-sig",newline=""))
        columns=csvr.fieldnames or []
        for row in csvr:
            n+=1
            for field in FIELDS:
                if field not in columns:
                    continue
                value=str(row.get(field) or "").strip()
                normalized=value.lower() or "[empty]"
                if normalized in MISSING:
                    missing[field]+=1
                distributions[field][normalized]+=1
            # Inspect original format only; deliberately do NOT use phenology outcomes.
            import re
            host=str(row.get("host") or "").strip()
            state=str(row.get("State") or "").strip()
            hexagon=str(row.get("hex_id") or "").strip()
            format_audit["valid_host_label"]+=int(host.lower() not in MISSING)
            format_audit["exact_binomial"]+=int(bool(re.fullmatch(r"[A-Z][a-z]+ [a-z][a-z-]*",host)))
            format_audit["state_nonempty"]+=int(state.lower() not in MISSING)
            format_audit["hex_nonempty"]+=int(hexagon.lower() not in MISSING)
            format_audit["host_state_hex_present"]+=int(host.lower() not in MISSING and state.lower() not in MISSING and hexagon.lower() not in MISSING)
            if host.lower() not in MISSING:
                host_names[host[:100]]+=1
            for field in distinct:
                if field not in columns:
                    continue
                value=str(row.get(field) or "").strip()
                if value.lower() not in MISSING:
                    distinct[field].add(value)
    return {
        "rows":n,
        "field_label_categories":{field:dict(sorted(v.items())) for field,v in distributions.items()
                                  if field in columns and field!="host"},
        "host_unique_taxon_labels":len(distinct["host"]),
        "geographic_hexagon_labels":len(distinct["hex_id"]),
        "source_calendar_years":len(distinct["year"]),
        "missing_label_counts":dict(sorted(missing.items())),
        "format_diagnostics":dict(sorted(format_audit.items())),
        "public_host_taxon_names_and_counts":dict(sorted(host_names.items())),
        "source_columns_present":{f:(f in columns) for f in (FIELDS+("hex_id","week","year","month","observed_on","latitude","longitude","observer"))},
    }


def run(outfile):
    result={"schema":"chocho_monarch_public_source_label_coverage_20261028_v01",
            "source_zip_sha256":ORIGINAL_ZIP_SHA256,
            "raw_observations_or_people_exported":False}
    with TemporaryDirectory(prefix="chocho-host-coverage-") as d:
        path=Path(d)/"origin.zip"
        result["transport"]=source_download(SHARE_URL,path)
        sha=hashlib.sha256(path.read_bytes()).hexdigest()
        if sha!=ORIGINAL_ZIP_SHA256:
            raise ValueError("Original public annotated source archive changed")
        with ZipFile(path) as z:
            if not all(x in z.namelist() for x in (CORE,THINNED,FULL)):
                raise ValueError("Frozen public original row sources are incomplete")
            result["core"]=count_file(z,CORE)
            result["thinned"]=count_file(z,THINNED)
            result["metadata_filenames"]=[x for x in z.namelist() if x.lower().endswith(".rtf")]
    if result["core"]["rows"]!=30320 or result["thinned"]["rows"]!=12942:
        raise ValueError("Published core or thinned cohort no longer matches")
    result["status"]="COMPLETE_SOURCE_CATEGORICAL_COVERAGE_ONLY"
    outfile.parent.mkdir(parents=True,exist_ok=True)
    outfile.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"status":result["status"],"core_rows":result["core"]["rows"],
                      "thinned_rows":result["thinned"]["rows"],
                      "thinned_categories":result["thinned"]["field_label_categories"]},sort_keys=True))
    return result


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",required=True,type=Path)
    a=p.parse_args()
    run(a.output)
