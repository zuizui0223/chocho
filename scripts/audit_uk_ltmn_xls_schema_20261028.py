#!/usr/bin/env python3
"""Read original UKBMS LTMN OLE XLS in-memory, inspect ONLY safe schemas.

The original site, transect, recorder, location and per-visit butterfly count
rows never leave the GitHub Action. No absence is derived from missing rows.
"""
from __future__ import annotations
import argparse,collections,hashlib,io,json,re
from pathlib import Path
from urllib.request import Request,urlopen

SOURCES={
"visits":("https://publications.naturalengland.org.uk/file/6020843017338880","c23945feb1568b4b60ac139136c485c234699710eae395a104e0be92a958a84d",9_000_000),
"counts":("https://publications.naturalengland.org.uk/file/6589256320155648","6afcc7d03eba0b657ae0e3a0d9da217b88df3ff39da66d7359df7fbb062ead3f",36_000_000)}
HEADER_RX=re.compile(r"(site|transect|visit|species|butterfly|count|week|date|year|start|route|section|number|location|record|name|id$)",re.I)

def get_source(key):
    url,sha,limit=SOURCES[key]
    req=Request(url,headers={"User-Agent":"chocho-original-UKBMS-LTMN-reanalysis/1.0"})
    with urlopen(req,timeout=65) as r:data=r.read(limit+1)
    if len(data)>limit or hashlib.sha256(data).hexdigest()!=sha:
        raise ValueError("Exact original LTMN source digest/size mismatch for "+key)
    if data[:8]!=bytes.fromhex("D0CF11E0A1B11AE1"):
        raise ValueError("Source no longer original OLE XLS "+key)
    return data

def headers_only(sheet):
    """Pick rows that resemble English field headers, never publish observations."""
    candidates=[]
    for i in range(min(20,sheet.nrows)):
        row=sheet.row_values(i,0,min(60,sheet.ncols))
        cells=[(col,str(value).strip()[:75]) for col,value in enumerate(row)
               if isinstance(value,str) and value.strip()]
        score=sum(bool(HEADER_RX.search(value)) for col,value in cells)
        nonempty=len(cells)
        if score>=2 and score>=max(2,nonempty//4):
            candidates.append((score,i,cells))
    candidates.sort(key=lambda x:(-x[0],x[1]))
    if not candidates:return {"source_header_status":"UNRESOLVED","header_candidates":[]}
    best=candidates[0]
    # Do not export any actual row values that did not match a known field name.
    matched=[{"column_number":col+1,"field_name":val} for col,val in best[2] if HEADER_RX.search(val)]
    return {"source_header_status":"LIKELY_FIELD_HEADER",
            "header_row_1based":best[1]+1,
            "recognized_fields":matched[:60],
            "recognized_field_count":len(matched),
            "header_row_total_string_cells":len(best[2])}

def run(path):
    try:import xlrd
    except ImportError as exc:raise ValueError("Read-only source OLE XLS decoder missing") from exc
    result={"schema":"chocho_LTMN_ukbms_original_visit_count_schema_20261028_v01",
            "inputs":{},"no_original_site_ids_observers_locations_counts_exported":True,
            "no_nondetections_inferred":True}
    for key in SOURCES:
        data=get_source(key)
        wb=xlrd.open_workbook(file_contents=data,on_demand=True)
        sheets=[]
        for sh in wb.sheets():
            m={"sheet_number":len(sheets)+1,"original_rows":sh.nrows,
               "original_columns":sh.ncols}
            m.update(headers_only(sh))
            # Only names of original field headers, never original source rows.
            if m.get("source_header_status")=="LIKELY_FIELD_HEADER":
                head=m["header_row_1based"]-1
                m["field_names_all_columns"]=[
                    str(v).strip()[:80] for v in sh.row_values(head)
                ]
            sheets.append(m)
        result["inputs"][key]={"sha256":SOURCES[key][1],"bytes":len(data),
                              "format":"OLE_XLS","sheets":sheets}
        wb.release_resources()
    result["status"]="SOURCE_VERIFIED_SCHEMA_NO_SURVEY_ABSENCES_INFERRED"
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps(result,sort_keys=True))
    return result

if __name__=="__main__":
    a=argparse.ArgumentParser();a.add_argument("--output",required=True,type=Path)
    run(a.parse_args().output)
