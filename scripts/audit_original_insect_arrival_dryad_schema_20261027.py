#!/usr/bin/env python3
"""Source-only public Dryad schema audit for distinct stages of insect invasion.

Fenn-Moltu 2023 DOI dryad.dv41ns243: source-country x destination x
commodity class x decade intercepted INSECT SPECIES RICHNESS, versus
native-bioregion x destination x decade established species richness.
Mally et al. 2022 DOI dryad.31zcrjdnw: LEPIDOPTERA species-level
nonnative establishment LIST plus family-level commodity records.

These have different taxonomic and geographical grains. Source preview
does not construct any invalid species×shipment×establishment join.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import time
import zipfile
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen

DOIS={
  "insect_flows_2023":"doi:10.5061/dryad.dv41ns243",
  "lepidoptera_establishment_2022":"doi:10.5061/dryad.31zcrjdnw",
}
EXPECTED={
 "insect_flows_2023":{"transport_flows_updated.csv",
                     "establishment_flows.csv"},
 "lepidoptera_establishment_2022":{
    "S5_list_of_non-native_species.txt","S8_commodities.txt",
    "S9_plant-commodities.txt",
 },
}
SCHEMA="chocho_original_public_dryad_insect_two_stage_grain_audit_20261027_v01"


def get_bytes(url):
    req=Request(url,headers={
        "User-Agent":"chocho-source-invasion-macro/0.1 (scientific data reuse; GitHub zuizui0223/chocho)",
        "Accept":"application/zip, application/json, text/plain, */*",
    })
    errors=[]
    for i in range(4):
        try:
            with urlopen(req,timeout=90) as r:
                if r.status!=200:
                    raise RuntimeError(f"HTTP {r.status}")
                raw=r.read(50_000_000)
                if len(raw)==50_000_000:
                    raise ValueError("Unexpected >50MB public source archive")
                return raw,r.url,r.headers.get("Content-Type")
        except (HTTPError,URLError,TimeoutError,RuntimeError) as exc:
            errors.append(str(exc)[:170])
            if i<3:
                time.sleep(min(2**(i+1),10))
    raise RuntimeError("Public Dryad source unavailable: "+str(errors))


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def schema_rows(name,raw):
    # Scientific metadata, never assume a text file is comma separated.
    text=raw.decode("utf-8-sig",errors="replace")
    sep="\t" if name.endswith(".txt") else ","
    first_lines=[line for line in text.splitlines()
                 if line.strip() and not line.startswith("#")]
    first=first_lines[0] if first_lines else ""
    if sep not in first and "," in first:
        sep=","
    if sep not in first and ";" in first:
        sep=";"
    reader=csv.reader(io.StringIO(text),delimiter=sep)
    parsed=list(reader)
    head=parsed[0] if parsed else []
    return {
        "filename":name,
        "byte_size":len(raw),
        "source_sha256":sha(raw),
        "source_nonempty_lines":len(first_lines),
        "source_detected_delimiter":repr(sep),
        "first_line_columns":head[:40],
        "first_5_raw_data_line_previews":[line[:350] for line in first_lines[1:6]],
        "source_duplicate_header_names":[h for h in set(head) if head.count(h)>1],
        "source_header_column_count":len(head),
        "source_ragged_row_lengths":sorted({
           len(x) for x in parsed[:min(1000,len(parsed))]})[:30],
    }


def acquire(doi):
    target="https://datadryad.org/api/v2/datasets/"+quote(doi,safe="")+"/download"
    content,url,mime=get_bytes(target)
    if not zipfile.is_zipfile(io.BytesIO(content)):
        raise ValueError(
            f"DOI {doi}: Dryad response is not a ZIP {mime}, first bytes {content[:130]!r}"
        )
    with zipfile.ZipFile(io.BytesIO(content)) as z:
        names=z.namelist()
        members={}
        for filename in names:
            path=Path(filename)
            if path.is_absolute() or ".." in path.parts:
                raise ValueError("ZIP path traversal in public original archive")
            if filename.endswith("/") or path.name.startswith("."):
                continue
            if path.name in members:
                raise ValueError("Ambiguous same-name public source files")
            members[path.name]=z.read(filename)
    return {
        "dataset_doi":doi,
        "original_public_source_download_api":target,
        "effective_source_redirect":url,
        "source_content_type":mime,
        "archive_size_bytes":len(content),
        "public_source_archive_SHA256":sha(content),
        "file_names":sorted(members),
        "files":[schema_rows(name,raw)
                 if name.lower().endswith((".txt",".csv",".md"))
                 else {"filename":name,"byte_size":len(raw),
                       "source_sha256":sha(raw),"nontext_source":True}
                 for name,raw in sorted(members.items())],
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output",required=True,type=Path)
    args=p.parse_args()
    out={}
    for alias,doi in DOIS.items():
        result=acquire(doi)
        actual=set(result["file_names"])
        missing=EXPECTED[alias]-actual
        if missing:
            raise ValueError(f"External source file missing: {alias}: {sorted(missing)}; actually: {sorted(actual)}")
        out[alias]=result
    report={
        "schema":SCHEMA,
        "status":"COMPLETE_ORIGINAL_TWO_PUBLIC_DATASET_SCHEMA_GRAINS_NO_SPECIES_COMMODITY_JOIN",
        "source_datasets":out,
        "frozen_source_ontology":{
            "insect_flows_2023":"Transport origin country x destination x commodity class x decade species RICHNESS; establishment native bioregion x destination x decade species RICHNESS. Not same species.",
            "lepidoptera_establishment_2022":"Species names by established recipient in S5; border interception commodity may be aggregated by Lepidoptera FAMILY or commodity. Do not infer individual butterfly arrival modes.",
        },
        "point_of_audit":"Read original ZIP member names and original schema before constructing any statistical model or joining datasets. Pure diagnostic, no invasion-rate estimates.",
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2,ensure_ascii=False,sort_keys=True)+"\n")
    print(json.dumps({
        "status":report["status"],
        "sources":{k:{
            "source_DOI":v["dataset_doi"],
            "source_sha256":v["public_source_archive_SHA256"],
            "files":v["file_names"],
            "source_schemas":v["files"],
        } for k,v in out.items()}
    },ensure_ascii=False,sort_keys=True))


if __name__=="__main__":
    main()
