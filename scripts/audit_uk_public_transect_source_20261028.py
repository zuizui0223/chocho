#!/usr/bin/env python3
"""Identify ORIGINAL UK transect visit+species count download availability.

No source rows, recorder IDs, site IDs, grid coordinates, or historical
butterfly outcomes are exported. At this stage, Excel files are only format-
checked; content is NOT interpreted as confirmed zero butterfly surveys.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,re,zipfile
from pathlib import Path
from tempfile import TemporaryDirectory
from urllib.request import Request,urlopen
from urllib.error import URLError,HTTPError

SOURCES={
"natural_england_ltmn_visits":{
"url":"https://publications.naturalengland.org.uk/file/6020843017338880",
"expected":"xls","max_bytes":9_000_000},
"natural_england_ltmn_weekly_counts":{
"url":"https://publications.naturalengland.org.uk/file/6589256320155648",
"expected":"xls","max_bytes":35_000_000},
"natural_england_ltmn_transect_site_crosswalk":{
"url":"https://publications.naturalengland.org.uk/file/6033934228652032",
"expected":"xlsx","max_bytes":2_000_000},
"eidc_ecn_weekly_butterflies_1993_2015":{
"url":"https://data-package.ceh.ac.uk/data/5aeda581-b4f2-4e51-b1a6-890b6b3403a3",
"expected":"zipcsv","max_bytes":60_000_000},
"eidc_ukbms_site_indices_1976_2024":{
"url":"https://data-package.ceh.ac.uk/data/53a1c705-4437-48bb-9bcd-b68977f8edd9",
"expected":"zipcsv","max_bytes":140_000_000},
}
MAX_MEMBER_UNPACKED=120_000_000

def describe(data:bytes,kind:str):
    """Strict magic-byte and public CSV schema inspection, never source rows."""
    ans={"format_detected":"UNRECOGNIZED"}
    if data[:8]==bytes.fromhex("D0CF11E0A1B11AE1"):
        ans["format_detected"]="LEGACY_OLE_XLS"
        return ans
    if data[:4]==b"PK\x03\x04":
        ans["format_detected"]="ZIP"
        with zipfile.ZipFile(io.BytesIO(data)) as z:
            if z.testzip() is not None:raise ValueError("Corrupt original ZIP CRC")
            files=[x for x in z.infolist() if not x.is_dir()]
            if sum(x.file_size for x in files)>400_000_000:
                raise ValueError("Unbounded compressed original source")
            csv_meta=[]
            suffixes={}
            for f in files:
                suffix=Path(f.filename).suffix.lower()
                suffixes[suffix]=suffixes.get(suffix,0)+1
                if suffix==".csv" and f.file_size<=MAX_MEMBER_UNPACKED:
                    with z.open(f) as c:
                        rd=csv.reader(io.TextIOWrapper(c,encoding="utf-8-sig",errors="replace",newline=""))
                        header=next(rd,[])
                        rows=sum(1 for _ in rd)
                    csv_meta.append({"original_file_index":len(csv_meta)+1,
                                      "rows":rows,
                                      "columns":[str(h).strip()[:80] for h in header][:80]})
            ans.update({"member_suffix_counts":suffixes,"csv_schema_metadata":csv_meta[:30],
                        "total_csv_source_rows":sum(x["rows"] for x in csv_meta)})
        return ans
    if data[:1] in (b"<",b"{") or b"<html" in data[:800].lower():
        ans["format_detected"]="HTML_OR_JSON_NOT_TABLE"
        return ans
    if b"," in data[:1000] and b"\n" in data[:1000]:
        ans["format_detected"]="TEXT_CSV"
        r=csv.reader(io.StringIO(data[:12000].decode("utf-8-sig",errors="replace")))
        ans["preview_header_only"]=next(r,[])[:80]
        return ans
    return ans

def fetch_source(key,source):
    limit=source["max_bytes"]
    req=Request(source["url"],headers={"User-Agent":"chocho-UK-transect-public-source/1.0",
                                       "Accept":"application/zip,application/octet-stream,*/*"})
    result={"source_key":key,"source_url":source["url"],"declared_format":source["expected"],
            "raw_site_count_or_observer_data_exported":False}
    try:
        with urlopen(req,timeout=65) as r:
            result["status_http"]=r.status
            result["final_domain"]=r.url.split("/")[2].lower()
            result["source_content_type"]=r.headers.get("Content-Type","")
            data=r.read(limit+1)
        if len(data)>limit:raise ValueError("Source exceeds bounded size")
        result.update({"source_bytes":len(data),"source_sha256":hashlib.sha256(data).hexdigest()})
        result.update(describe(data,source["expected"]))
        accepted=(
           source["expected"]=="xls" and result["format_detected"]=="LEGACY_OLE_XLS"
           or source["expected"]=="xlsx" and result["format_detected"]=="ZIP"
           or source["expected"]=="zipcsv" and
              (result["format_detected"]=="TEXT_CSV" or
               (result["format_detected"]=="ZIP" and result.get("total_csv_source_rows",0)>0))
        )
        result["source_stage"]= "ORIGINAL_SOURCE_FORMAT_AND_SCHEMA_VERIFIED" if accepted else "INCOMPATIBLE_OR_ACCESS_PAGE"
    except (OSError,URLError,HTTPError,ValueError,zipfile.BadZipFile) as e:
        result["source_stage"]="ORIGINAL_SOURCE_UNAVAILABLE_OR_BLOCKED"
        result["error_type"]=type(e).__name__
        result["error_details"]=str(e)[:250]
    return result

def run(path):
    outcomes=[fetch_source(key,obj) for key,obj in SOURCES.items()]
    complete=[x["source_key"] for x in outcomes if x["source_stage"]=="ORIGINAL_SOURCE_FORMAT_AND_SCHEMA_VERIFIED"]
    result={"schema":"chocho_uk_public_transect_visits_counts_source_feasibility_20261028_v01",
            "status":"PUBLIC_SOURCE_AT_LEAST_ONE_AVAILABLE" if complete else "NO_PUBLIC_ORIGINAL_SOURCE_VERIFIED",
            "verified_sources":complete,"source_count":len(outcomes),"sources":outcomes,
            "any_source_rows_or_precise_transect_sites_exported":False,
            "actual_survey_absences_verified":False,
            "source_linkage_or_true_nondetections_inferred":False}
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({"status":result["status"],
       "sources":[{"key":x["source_key"],"status":x["source_stage"],
                   "format":x.get("format_detected"),
                   "bytes":x.get("source_bytes"),
                   "csv_rows":x.get("total_csv_source_rows"),
                   "error":x.get("error_details")} for x in outcomes]},ensure_ascii=False))
    return result

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",required=True,type=Path)
    args=parser.parse_args()
    result=run(args.output)
    if result["status"]!="PUBLIC_SOURCE_AT_LEAST_ONE_AVAILABLE":
        raise SystemExit(2)
