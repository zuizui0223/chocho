#!/usr/bin/env python3
"""Recover publisher-blob-verified GLSN2018 Excel PORT and EDGE original research sources.

No openpyxl dependency: extract XLSX XML with Python stdlib, preserve all original rows.
Research source: Network-Maritime-Complexity/GLSN-efficiency-vulnerability (MIT).
Neither edge means direct cargo journey nor butterfly dispersal.
"""
from __future__ import annotations
import argparse,base64,csv,hashlib,io,json,posixpath,re,time,urllib.request,urllib.error,zipfile
from datetime import datetime,timezone
from pathlib import Path
import xml.etree.ElementTree as ET

ROOT=Path(__file__).resolve().parents[1]
PUB="Network-Maritime-Complexity/GLSN-efficiency-vulnerability"
GIT_API="https://api.github.com/repos/"+PUB+"/git/blobs/"
FILES={
    "original_GLNS_ports_2018.xlsx":{"path":"data/port 2018.xlsx","sha":"5486d8c204a7c49b375a14762fe5cb36c7877d96","size":34675,"role":"PORTS"},
    "original_GLNS_edges_2018.xlsx":{"path":"data/Edge 2018_with great circle distance (GC) and weight.xlsx","sha":"18098cab432ddbcd7e0a7091fc1c4cedbda40fb7","size":1018504,"role":"EDGES"},
}
NS="http://schemas.openxmlformats.org/spreadsheetml/2006/main"
REL="http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PKG="http://schemas.openxmlformats.org/package/2006/relationships"

def github_blob_sha(data):
    return hashlib.sha1(b"blob "+str(len(data)).encode()+b"\0"+data).hexdigest()

def blob_download(sha):
    url=GIT_API+sha
    err=None
    for retry in range(4):
        try:
            req=urllib.request.Request(url,headers={
                "User-Agent":"chocho-GLSN2018-source-recovery/1.0",
                "Accept":"application/vnd.github+json",
                "X-GitHub-Api-Version":"2022-11-28"})
            with urllib.request.urlopen(req,timeout=60) as res:
                data=res.read(6_000_000)
            doc=json.loads(data)
            if doc.get("encoding")!="base64":
                raise ValueError("GitHub blob API source encoding not base64")
            binary=base64.b64decode(doc["content"])
            if github_blob_sha(binary)!=sha or doc.get("sha")!=sha:
                raise ValueError("Original publisher Git blob identity mismatch")
            return binary,{"url":url,"http_status":200,"retrieved_utc":datetime.now(timezone.utc).isoformat()}
        except (urllib.error.HTTPError,urllib.error.URLError,TimeoutError,ConnectionError) as exc:
            err=str(exc)
            if retry<3:time.sleep(3*(retry+1))
        except Exception:
            raise
    raise RuntimeError("Original publisher Git blob could not be retrieved: "+str(err))

def column_index(reference):
    name=re.match(r"^[A-Z]+",reference or "")
    if name is None:raise ValueError("Excel cell missing alphabetic reference")
    result=0
    for c in name.group(0):
        result=result*26+ord(c)-64
    return result-1

def sheet_tables(data):
    z=zipfile.ZipFile(io.BytesIO(data))
    wb=ET.fromstring(z.read("xl/workbook.xml"))
    relationships=ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
    by_rel={item.get("Id"):item.get("Target") for item in relationships}
    sst=[]
    if "xl/sharedStrings.xml" in z.namelist():
        ss=ET.fromstring(z.read("xl/sharedStrings.xml"))
        for si in ss.findall(f"{{{NS}}}si"):
            sst.append("".join(e.text or "" for e in si.iter(f"{{{NS}}}t")))
    result=[]
    for sh in wb.find(f"{{{NS}}}sheets"):
        rel=sh.get(f"{{{REL}}}id")
        path=by_rel.get(rel)
        if not path:raise ValueError("Missing XLSX worksheet relationship")
        path=(path.lstrip("/") if path.startswith("/") else posixpath.normpath(posixpath.join("xl",path)))
        if path not in z.namelist():raise ValueError("XLSX sheet target missing: "+path)
        sheet=ET.fromstring(z.read(path))
        rows=[]
        body=sheet.find(f"{{{NS}}}sheetData")
        if body is None:raise ValueError("XLSX sheetData not found")
        for line in body.findall(f"{{{NS}}}row"):
            fields={}
            for cell in line.findall(f"{{{NS}}}c"):
                x=column_index(cell.get("r"))
                kind=cell.get("t")
                value=cell.find(f"{{{NS}}}v")
                if kind=="inlineStr":
                    inline=cell.find(f"{{{NS}}}is")
                    text="" if inline is None else "".join(x.text or "" for x in inline.iter(f"{{{NS}}}t"))
                elif value is None:
                    text=""
                elif kind=="s":
                    idx=int(value.text)
                    if idx<0 or idx>=len(sst):raise ValueError("Invalid XLSX shared string index")
                    text=sst[idx]
                else:
                    text=value.text or ""
                fields[x]=text
            if fields:
                n=max(fields)+1
                if n>100:raise ValueError("Source XLSX extremely wide unexpected columns")
                rows.append([fields.get(i,"") for i in range(n)])
            if len(rows)>200000:
                raise ValueError("Original XLSX exceeds source-only source size limit")
        if not rows:continue
        width=max(map(len,rows))
        padded=[row+[""]*(width-len(row)) for row in rows]
        result.append({"name":sh.get("name"),"id":sh.get("sheetId"),"rows":padded,"width":width})
    if not result:raise ValueError("No rows in publisher XLSX")
    return result

def write_csv(path,rows):
    with Path(path).open("w",encoding="utf-8",newline="") as stream:
        writer=csv.writer(stream)
        writer.writerows(rows)

def run(outdir):
    dest=Path(outdir).resolve()
    if any(dest.is_relative_to((ROOT/s).resolve()) for s in ("data","results","provenance")):
        raise ValueError("GLSN source recovery cannot overwrite original source")
    dest.mkdir(parents=True,exist_ok=True)
    receipts=[]
    for local,info in FILES.items():
        raw_path=dest/local
        if raw_path.exists():
            b=raw_path.read_bytes();origin={"retrieved_utc":"verified_local_cache","http_status":None}
        else:
            b,origin=blob_download(info["sha"])
            raw_path.write_bytes(b)
        if len(b)!=info["size"] or github_blob_sha(b)!=info["sha"]:
            raise ValueError("GLSN 2018 XLSX original publisher size or Git SHA-1 changed")
        sheets=sheet_tables(b)
        sheet_info=[]
        for k,sheet in enumerate(sheets,1):
            path=dest/f"GLSN2018_{info['role'].lower()}_original_sheet_{k}.csv"
            write_csv(path,sheet["rows"])
            sheet_info.append({"sheet_index":k,"name":sheet["name"],
                "all_nonempty_rows":len(sheet["rows"]),"columns":sheet["width"],
                "column_header_candidates":sheet["rows"][0][:60],
                "example_rows":sheet["rows"][1:4],
                "extracted_csv_file":path.name,
                "extracted_csv_sha256":hashlib.sha256(path.read_bytes()).hexdigest()})
        receipts.append({"role":info["role"],"original_publisher_path":info["path"],
            "original_publisher_git_blob":info["sha"],
            "original_recovered_bytes":len(b),"original_sha256":hashlib.sha256(b).hexdigest(),
            "publisher_license":"MIT","source_url":"https://github.com/"+PUB+"/blob/main/"+info["path"].replace(" ","%20"),
            "download":origin,"sheets":sheet_info})
        print("GLSN2018_SOURCE_FILE="+json.dumps({"role":info["role"],"git_blob_sha":info["sha"],
                "sheets":[{"name":x["name"],"rows":x["all_nonempty_rows"],"headers":x["column_header_candidates"][:16],
                           "samples":x["example_rows"][:2]} for x in sheet_info]},sort_keys=True),flush=True)
    receipt={"schema":"chocho_GLSN2018_recovered_original_github_excel_source_v1",
        "status":"ORIGINAL_PUBLISHER_GIT_BLOB_XLSX_AND_FULL_EXTRACTED_TABLES_VERIFIED",
        "publisher_repo":PUB,"license":"MIT",
        "observed_real_port_to_port_voyages":False,
        "liner_network_edges_not_necessarily_consecutive_port_calls":True,
        "source_retrieval_utc":datetime.now(timezone.utc).isoformat(),
        "sources":receipts}
    (dest/"GLSN2018_recovered_source_receipt.json").write_text(json.dumps(receipt,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    return receipt

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--out-dir",type=Path,required=True)
    run(parser.parse_args().out_dir)
