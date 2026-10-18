#!/usr/bin/env python3
"""Pre-outcome source audit of published monarch host-performance classifications.

PLOS ONE 2022 S1 appendix (CC BY), DOI 10.1371/journal.pone.0269701.s001.
This script reports table shapes and whether performance-coded rows are present
but never automatically asserts a 40-HOSTS/WCVP accepted-ID quality match or
identifies larval breeding at any global locality. All classification uses
require manual source-table review and taxonomic concordance.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import re
import zipfile
from collections import Counter
from pathlib import Path
from xml.etree import ElementTree as ET

NS="{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
CODES={"H1","H2","H3","L1","L2","L3","N","U"}
NAMES=re.compile(r"\b[A-Z][a-z]+\s+[a-z][a-z\-]+\b")


def tables_in_docx(path:Path)->list[list[list[str]]]:
    with zipfile.ZipFile(path) as zipf:
        if "word/document.xml" not in zipf.namelist():
            raise ValueError("not a readable OOXML supplement")
        root=ET.fromstring(zipf.read("word/document.xml"))
    tables=[]
    for table in root.iter(NS+"tbl"):
        cells=[]
        for row in table.findall(NS+"tr"):
            rec=[]
            for td in row.findall(NS+"tc"):
                text=" ".join(
                    (x.text or "").strip()
                    for x in td.iter(NS+"t")
                    if x.text and x.text.strip()
                )
                rec.append(" ".join(text.split()))
            cells.append(rec)
        tables.append(cells)
    return tables


def code_in_row(row:list[str])->str | None:
    found=[]
    for cell in row:
        for token in re.split(r"[;,/:\s()]+",cell):
            token=token.strip().upper()
            if token in CODES:
                found.append(token)
    choices=set(found)
    if len(choices)==1:
        return next(iter(choices))
    return None


def audit_tables(tables:list[list[list[str]]])->dict:
    items=[]
    for i,tab in enumerate(tables):
        codes=Counter(filter(None,(code_in_row(r) for r in tab)))
        species_like=sum(
            1 for row in tab if any(
                NAMES.search(c or "") for c in row
            )
        )
        items.append({
            "table_index":i,
            "rows":len(tab),
            "max_columns":max((len(x) for x in tab),default=0),
            "binomial_like_rows":species_like,
            "single_status_code_rows":sum(codes.values()),
            "detected_singleton_status_code_counts":dict(codes),
            # No raw taxonomy/labels extracted at this data-feasibility stage.
        })
    return {
        "schema":"chocho_monarch_2022_host_quality_supplement_schema_v0.1",
        "status":"PUBLISHED_DOCUMENT_TABLE_AUDIT_NOT_GEOGRAPHIC_FITNESS_RESULT",
        "reference":"Greenstein et al. 2022, PLOS ONE, 10.1371/journal.pone.0269701",
        "source_appendix_doi":"10.1371/journal.pone.0269701.s001",
        "peer_reviewed_reported_host_classes":{
            "high_performance":34,
            "low_performance":42,
            "nonhost":33,
            "unsubstantiated":18,
            "total":127,
        },
        "tables":items,
        "eligibility":{
            "at_least_one_table_with_100_plus_binomial_rows":any(
                x["binomial_like_rows"]>=100 for x in items
            ),
            "exact_127_status_row_table_present":any(
                x["single_status_code_rows"]>=120 and x["binomial_like_rows"]>=100
                for x in items
            ),
            "independent_40_HOSTS_accepted_taxa_crosswalk":"NOT_PERFORMED",
            "butterfly_new_region_quality_filter":"NOT_PERFORMED",
            "observed_larval_recruitment":"NOT_MEASURED",
        },
        "limitations":[
            "Row counts alone cannot validate botanical names, synonymy, evidence class or the reported 34/42/33/18 split.",
            "Host-quality categories cannot be projected as fixed survival rates at every WGSRPD3 region.",
            "The dataset concerns Danaus plexippus, not global performance measured across 239 butterflies.",
        ]
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--supplement-docx",type=Path,required=True)
    p.add_argument("--output-json",type=Path,required=True)
    a=p.parse_args()
    result=audit_tables(tables_in_docx(a.supplement_docx))
    result["source_content_sha256"]=hashlib.sha256(
        a.supplement_docx.read_bytes()
    ).hexdigest()
    a.output_json.parent.mkdir(parents=True,exist_ok=True)
    a.output_json.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__=="__main__":
    main()
