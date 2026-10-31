#!/usr/bin/env python3
"""Fetch already-published butterfly occurrence source (14,637,152 GBIF rows).

This uses the EXISTING successful GBIF DOI https://doi.org/10.15468/dl.umqmsf
and downloads only its existing public ZIP. No account credentials, new GBIF
download request, or per-island GBIF search is required.

Archive is streamed to file with SHA256 check against the live GBIF archive
size, not extracted into runner disk. This stage introspects the CSV schema
and original predicate; it does NOT claim 14m usable precise island records.
"""
from __future__ import annotations
import argparse,hashlib,json,os,urllib.request,urllib.error,zipfile,io,time
from datetime import datetime,timezone
from pathlib import Path

HERE=Path(__file__).resolve().parent
API="https://api.gbif.org/v1/occurrence/download/10.15468/dl.umqmsf"
DOI="10.15468/dl.umqmsf"
EXPECTED_DOWNLOAD_KEY="0013264-261023074918325"
EXPECTED_TOTAL_RECORDS=14_637_152
UA="chocho-global-public-GH-workflow-original-GBIF-DOI-reuse/20261031"

def fetch_metadata():
    req=urllib.request.Request(API,headers={"Accept":"application/json","User-Agent":UA})
    with urllib.request.urlopen(req,timeout=60) as r:
        j=json.load(r)
    if (j.get("doi")!=DOI or j.get("key")!=EXPECTED_DOWNLOAD_KEY
        or j.get("status")!="SUCCEEDED"
        or int(j.get("totalRecords",0))!=EXPECTED_TOTAL_RECORDS):
        raise ValueError("Published original GBIF download metadata changed: "+
                         str({k:j.get(k) for k in ("doi","key","status","totalRecords")}))
    link=j.get("downloadLink")
    if not link or not link.startswith("https://api.gbif.org/"):
        raise ValueError("Unexpected original download source URL")
    size=int(j.get("size",0))
    if not 10**9<=size<=2*10**9:
        raise ValueError("Expected known original ~1.42 GB ZIP, unexpected size "+str(size))
    return j,size,link

def ensure_archive(path,size,link):
    if path.is_file() and path.stat().st_size==size and zipfile.is_zipfile(path):
        return {"cached_original_existing":True,"size":size}
    path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix(".partial")
    if temp.exists():temp.unlink()
    h=hashlib.sha256()
    written=0
    req=urllib.request.Request(link,headers={"User-Agent":UA,"Accept":"application/zip,*/*"})
    with urllib.request.urlopen(req,timeout=150) as resp, temp.open("wb") as out:
        if resp.status!=200:raise ValueError("Download source status="+str(resp.status))
        while True:
            chunk=resp.read(8*1024*1024)
            if not chunk:break
            written+=len(chunk)
            if written>2*10**9:raise ValueError("More ZIP bytes than size cap")
            h.update(chunk)
            out.write(chunk)
            if written%(128*1024*1024)<len(chunk):
                print("ORIGINAL_SOURCE_GBIF_ZIP_STREAMED",written,"of",size,flush=True)
    if written!=size or not zipfile.is_zipfile(temp):
        raise ValueError("Original GBIF ZIP download incomplete or wrong bytes "+
                         str(written)+" expected "+str(size))
    temp.replace(path)
    return {"cached_original_existing":False,"size":written,
            "sha256_newly_downloaded":h.hexdigest()}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--zip",type=Path,required=True)
    args=p.parse_args()
    metadata,size,link=fetch_metadata()
    receipt=ensure_archive(args.zip,size,link)
    digest=hashlib.sha256(args.zip.read_bytes()).hexdigest()
    with zipfile.ZipFile(args.zip) as z:
        info=[{"filename":i.filename,"compressed_size":i.compress_size,
               "uncompressed_size":i.file_size,
               "compression_method":i.compress_type,
               "has_large_uncompressed_size":i.file_size>2*10**9}
              for i in z.infolist()]
        possible=[x for x in z.infolist() if x.filename.lower().endswith(
            ("occurrence.txt","occurrence.csv","simple.csv","occurrence.tsv"))]
        if not possible:
            possible=sorted([x for x in z.infolist() if x.filename.lower().endswith((".csv",".tsv",".txt"))],
                            key=lambda x:x.file_size,reverse=True)[:1]
        if len(possible)!=1 or possible[0].file_size<1*10**9:
            raise ValueError("Unexpected GBIF ZIP no unambiguous occurrence table "+
                             repr(info)[:1400])
        mainfile=possible[0]
        with z.open(mainfile) as r:
            raw=[r.readline(1_000_000) for _ in range(3)]
        first=[x.decode("utf8",errors="replace").strip()[:1500] for x in raw]
    req=metadata.get("request") or {}
    predicate=req.get("predicate") if isinstance(req,dict) else {}
    if isinstance(predicate,str):
        try:predicate=json.loads(predicate)
        except Exception:predicate={"predicate_format":"unknown_text"}
    source={
      "schema":"chocho_original_public_14m_butterfly_GBIF_archive_source_audit_v1",
      "status":"EXISTING_ORIGINAL_GBIF_ZIP_RETRIEVED_AND_INSPECTED_NOT_YET_SPATIAL_JOINED",
      "checked_at_utc":datetime.now(timezone.utc).isoformat(),
      "doi":"https://doi.org/"+DOI,
      "download_key":metadata.get("key"),"total_original_download_records":metadata.get("totalRecords"),
      "zip_bytes":size,"zip_sha256":digest,"cached_file_path_not_committed":str(args.zip),
      "download_source":link,
      "source_request_format":req.get("format"),
      "source_request_predicate":predicate,
      "ZIP_entry_infos":info,
      "main_occurrence_zip_member":mainfile.filename,
      "main_occurrence_uncompressed_bytes":mainfile.file_size,
      "first_three_csv_rows_preview":first,
      "source_license":metadata.get("license"),
      "records_are_NOT_all_confirmed_geo_accurate_or_on_target_islands":True,
      "original_private_14425536_record_panel_not_reconstructed_yet":True,
      "download_metadata_frozen_at_utc":metadata.get("created"),
      "download_archive_retained_until":metadata.get("eraseAfter"),
      "checked_zip_bytes":receipt["size"],
      "from_existing_local_cached_zip":receipt["cached_original_existing"],
    }
    HERE.mkdir(exist_ok=True,parents=True)
    (HERE/"original_GBIF_14m_source_zip_audit.json").write_text(
      json.dumps(source,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps({k:source[k] for k in
      ("status","total_original_download_records","zip_bytes","zip_sha256",
       "main_occurrence_zip_member","main_occurrence_uncompressed_bytes",
       "source_request_format","first_three_csv_rows_preview")},
      indent=2,ensure_ascii=False),flush=True)

if __name__=="__main__":main()
