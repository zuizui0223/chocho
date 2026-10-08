#!/usr/bin/env python3
"""Audit Figshare 25648644 public files before any stage-vs-host-composition model.

Do not infer missing survey counts, run downloaded code, or treat route-level counts
as individual egg-to-larva survival.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,re,time,urllib.error,urllib.parse,urllib.request,zipfile
from pathlib import Path

ARTICLE_ID=25648644
MAX_FILE_SIZE=50_000_000
MAX_TOTAL=120_000_000
RECOGNIZED={".csv",".tsv",".txt",".r",".rmd",".xlsx",".rds",".rdata",".rda",".zip",".md",".json",".pdf"}
def get(url,max_bytes=MAX_FILE_SIZE):
    req=urllib.request.Request(url,headers={
        "User-Agent":"chocho-open-research-source-gate/0.1 (github.com/zuizui0223/chocho)",
        "Accept":"application/json,application/octet-stream;q=0.8,*/*;q=0.5"
    })
    with urllib.request.urlopen(req,timeout=70) as response:
        data=response.read(max_bytes+1)
        if len(data)>max_bytes:raise ValueError("file exceeds size gate")
        return data
def payload(url):
    data=get(url,2_000_000)
    return json.loads(data)
def endpoint_probe():
    urls=[
      "https://api.figshare.com/v2/articles/25648644",
      "https://api.figshare.com/v2/articles/25648644/versions/1"
    ]
    errors=[]
    for url in urls:
        for attempt in range(3):
            try:
                obj=payload(url)
                if not isinstance(obj,dict) or not isinstance(obj.get("files"),list):
                    errors.append({"url":url,"error":"metadata missing files array"})
                    break
                if str(obj.get("id"))!="25648644":
                    raise RuntimeError("Figshare returned unexpected article id")
                return obj,errors
            except Exception as e:
                errors.append({"url":url,"attempt":attempt+1,"error":str(e)[:400]})
                if attempt<2:time.sleep(min(15,3*(attempt+1)))
    return None,errors
def schema_extract(path,data):
    ext=path.suffix.lower()
    out={"name":path.name,"bytes":len(data),"sha256":hashlib.sha256(data).hexdigest(),"extension":ext}
    if ext in (".csv",".tsv",".txt"):
        text=data[:300000].decode("utf-8-sig",errors="replace")
        delimiter="\t" if ext==".tsv" else ","
        if ext==".txt":
            try:delimiter=csv.Sniffer().sniff(text[:10000],delimiters=",\t;").delimiter
            except csv.Error:return out
        rows=csv.reader(io.StringIO(text),delimiter=delimiter)
        head=next(rows,[])
        out["header"]=head[:100]
        out["first_rows"]=[next(rows,[])[:15] for _ in range(2)]
        out["approx_total_newlines"]=data.count(b"\n")
    elif ext==".zip":
        try:
            with zipfile.ZipFile(io.BytesIO(data)) as z:
                out["archive_members"]=[{"name":v.filename,"size":v.file_size} for v in z.infolist()[:200]]
        except (zipfile.BadZipFile,ValueError) as e:
            out["bad_archive"]=str(e)
    elif ext in (".r",".rmd",".md"):
        out["first_text"]=data[:3000].decode("utf-8",errors="replace")
    return out
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--protocol-json",type=Path,required=True)
    ap.add_argument("--output-dir",type=Path,required=True)
    a=ap.parse_args()
    p=json.loads(a.protocol_json.read_text(encoding="utf8"))
    if p["figshare_article_id"]!=ARTICLE_ID:raise RuntimeError("Figshare study id drift")
    a.output_dir.mkdir(parents=True,exist_ok=True)
    metadata,errors=endpoint_probe()
    result={
      "schema":"chocho_figshare_urban_milkweed_raw_data_source_gate_v0.1",
      "study":"Erickson et al. 2025, Ecosphere 10.1002/ecs2.70259",
      "expected_figshare_article_id":ARTICLE_ID,
      "source_metadata_accessible":metadata is not None,
      "source_errors":errors,
      "files":[],
      "files_downloaded":[],
      "data_completeness_unknown":True,
      "does_not_fit_stage_model":True,
      "decision":"SOURCE_UNAVAILABLE" if metadata is None else "SOURCE_METADATA_ACCESSIBLE_NOT_YET_MODELLED"
    }
    if metadata:
        result["figshare_metadata"]={k:metadata.get(k) for k in ("id","title","doi","published_date","modified_date","version","url_public_html","license","description")}
        filelist=metadata["files"]
        result["files"]=[{k:f.get(k) for k in ("id","name","size","is_link_only","download_url","computed_md5")} for f in filelist]
        total=0
        for f in filelist:
            name=Path(str(f.get("name") or "unnamed")).name
            size=int(f.get("size") or 0)
            suffix=Path(name).suffix.lower()
            if suffix not in RECOGNIZED or size<=0 or size>MAX_FILE_SIZE or total+size>MAX_TOTAL:
                result["files_downloaded"].append({"name":name,"status":"not_downloaded_size_or_extension"})
                continue
            try:
                url=f.get("download_url") or f"https://ndownloader.figshare.com/files/{int(f['id'])}"
                data=get(url)
                if len(data)!=size:
                    raise RuntimeError(f"download mismatch: expected {size} got {len(data)}")
                path=a.output_dir/("figshare_"+str(f["id"])+"_"+re.sub(r"[^a-zA-Z0-9._-]","_",name))
                path.write_bytes(data)
                total+=len(data)
                result["files_downloaded"].append({"name":name,"status":"downloaded","source_file_id":f["id"],"local_basename":path.name})
                result.setdefault("file_schemas",[]).append(schema_extract(path,data))
            except Exception as e:
                result["files_downloaded"].append({"name":name,"status":"download_failed","error":str(e)[:400]})
            time.sleep(.5)
        result["decision"]="SOURCE_READY_FOR_SCHEMA_REVIEW" if any(x.get("status")=="downloaded" for x in result["files_downloaded"]) else "SOURCE_METADATA_ONLY"
    (a.output_dir/"figshare_source_inventory_v01.json").write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps({"decision":result["decision"],
        "metadata":result["source_metadata_accessible"],
        "n_files":len(result["files"]),
        "downloaded":[x["name"] for x in result["files_downloaded"] if x.get("status")=="downloaded"],
        "errors":result["source_errors"][-4:]},indent=2),flush=True)
if __name__=="__main__":main()
