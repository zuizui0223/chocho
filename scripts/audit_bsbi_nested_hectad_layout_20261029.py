#!/usr/bin/env python3
"""Audit original *nested* BSBI 2020 historical hectad layout, source only.

Original 2024 Zenodo ZIP has a smaller nested hectadDistributions.zip. This
inspection publishes only source file names, CSV header fields and row counts;
never taxa×hectad occurrence bodies or original UKBMS transect locations.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,zipfile
from pathlib import Path
from collections import Counter
from audit_bsbi_historical_hectad_source_20261029 import (
    obtain,EXPECTED_MD5,inspect_archive,classify_header,MAX_UNPACKED
)

NESTED="PA2020_10by10kmSquareData/hectadDistributions.zip"
SCHEMA="chocho_bsbi_PA2020_original_nested_hectad_date_classes_schema_20261029_v01"

def inspect_nested(source:bytes):
    with zipfile.ZipFile(io.BytesIO(source)) as z:
        if NESTED not in z.namelist():raise ValueError("Expected original 2024 nested hectad file absent")
        blob=z.read(NESTED)
    with zipfile.ZipFile(io.BytesIO(blob)) as nested:
        if nested.testzip() is not None:raise ValueError("Nested original source CRC invalid")
        members=[x for x in nested.infolist() if not x.is_dir()]
        if not members or sum(x.file_size for x in members)>650_000_000:
            raise ValueError("Nested original file decompression size unsupported")
        extensions=Counter(Path(x.filename).suffix.lower() for x in members)
        schemas=[]
        for i,info in enumerate(members):
            path=Path(info.filename)
            if path.is_absolute() or ".." in path.parts:
                raise ValueError("Original nested data path is unsafe")
            if Path(info.filename).suffix.lower() not in {".csv",".tsv",".txt"}:
                continue
            if len(schemas)>=18:
                continue
            delimiter="\t" if path.suffix==".tsv" else ","
            with nested.open(info) as stream:
                rdr=csv.reader(io.TextIOWrapper(stream,encoding="utf-8-sig",
                        errors="replace",newline=""),delimiter=delimiter)
                header=next(rdr,[])
                # Count rows by streaming so never serialize individual records.
                n=sum(1 for _ in rdr)
            schemas.append({
                "original_member_filename":info.filename[:180],
                "original_member_uncompressed_bytes":info.file_size,
                "row_count_no_original_records_exported":n,
                **classify_header(header),
            })
        return {
            "original_nested_filename":NESTED,
            "original_nested_zip_sha256":hashlib.sha256(blob).hexdigest(),
            "original_nested_member_count":len(members),
            "original_nested_uncompressed_bytes_total":sum(x.file_size for x in members),
            "original_nested_extension_counts":dict(sorted(extensions.items())),
            "nested_original_member_schema_sample":schemas,
            "no_hectad_species_period_original_rows_exported":True,
        }

def run(output:Path):
    result={
        "schema":SCHEMA,"original_zenodo_source":"https://zenodo.org/records/10998003",
        "original_Zenodo_expected_MD5":EXPECTED_MD5,
        "individual_plant_or_butterfly_observations_exported":False
    }
    try:
        source,url,errors=obtain()
        result["original_source_zip_sha256"]=hashlib.sha256(source).hexdigest()
        result["original_source_url"]=url
        result.update(inspect_nested(source))
        result["status"]="COMPLETE_NESTED_ORIGINAL_BSBI_DATED_GRID_SCHEMA"
    except (ValueError,OSError,zipfile.BadZipFile) as exc:
        result["status"]="ORIGINAL_SOURCE_UNAVAILABLE_OR_UNREADABLE"
        result["blocker_class"]=type(exc).__name__
        result["blocker"]=str(exc)[:320]
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,sort_keys=True,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({k:result.get(k) for k in ("status","original_nested_member_count",
          "original_nested_extension_counts","original_nested_uncompressed_bytes_total",
          "nested_original_member_schema_sample","blocker")},ensure_ascii=False))
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",required=True,type=Path)
    x=run(p.parse_args().output)
    if x["status"]!="COMPLETE_NESTED_ORIGINAL_BSBI_DATED_GRID_SCHEMA":raise SystemExit(2)
