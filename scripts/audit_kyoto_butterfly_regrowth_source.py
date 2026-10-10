#!/usr/bin/env python3
"""Audit original Kyoto native/exotic butterfly competition dataset, WITHOUT new effects."""
from __future__ import annotations
import argparse,csv,hashlib,io,json
from pathlib import Path
from urllib.request import urlopen,Request
from urllib.error import HTTPError,URLError

ARTICLE_ID=23170898
ARTICLE_URL=f"https://api.figshare.com/v2/articles/{ARTICLE_ID}"
EXPECTED=["Hashimoto_and_Ohgushi_2023_larvaldata.csv","Hashimoto_and_Ohgushi_2023_plantdata.csv"]

def get(url,limit):
    with urlopen(Request(url,headers={"User-Agent":"chocho-original-butterfly-scientific-audit/1.0"}),timeout=40) as res:
        body=res.read(limit+1)
    if not body or len(body)>limit:
        raise ValueError("missing or oversized source")
    return body

def csv_metadata(data):
    text=data.decode("utf-8-sig")
    start=text.splitlines()[:15]
    try:
        dialect=csv.Sniffer().sniff("\n".join(start),delimiters=",;\t")
        delim=dialect.delimiter
    except csv.Error:
        delim=","
    dr=csv.DictReader(io.StringIO(text),delimiter=delim)
    header=dr.fieldnames
    if not header or len(header)<3 or len(header)>150:
        raise ValueError("invalid original field structure")
    count=0
    missing={}
    categories={}
    for r in dr:
        if None in r:
            raise ValueError("inconsistent CSV row widths")
        count+=1
        for name in header:
            if not str(r.get(name) or "").strip():
                missing[name]=missing.get(name,0)+1
    return {"records":count,"columns":header,"missing_by_field":missing,
            "delimiter":"TAB" if delim=="\t" else delim}

def audit():
    metadata=json.loads(get(ARTICLE_URL,1_000_000))
    if metadata.get("id")!=ARTICLE_ID:
        raise ValueError("not the original source article")
    files=metadata.get("files",[])
    receipt={
       "schema":"chocho_original_kyoto_butterfly_competition_source_v01",
       "article_id":ARTICLE_ID,"title":metadata.get("title"),
       "doi":metadata.get("doi"),"licence":metadata.get("license"),
       "files":[],"all_expected_files_verified":False,
       "new_regrowth_timing_effect_estimated":False,
       "temporal_resource_pulse_randomization_verified":False,
       "GEB_PR38_untouched":True
    }
    successes=set()
    for f in files:
        name=f.get("name") or ""
        item={"name":name,"size":f.get("size"),"file_id":f.get("id"),"download_url":f.get("download_url")}
        if not name.lower().endswith((".csv",".tsv")):
            item["status"]="UNPARSED_NONCSV"
            receipt["files"].append(item);continue
        try:
            raw=get(f["download_url"],4_000_000)
            actual_md5=hashlib.md5(raw).hexdigest()
            expected=f.get("computed_md5") or f.get("supplied_md5")
            if expected and expected!=actual_md5:
                raise ValueError("original md5 mismatch")
            item.update({"sha256":hashlib.sha256(raw).hexdigest(),"md5":actual_md5,
                         "schema":csv_metadata(raw),"status":"ORIGINAL_CSV_VERIFIED"})
            if name in EXPECTED:successes.add(name)
        except (ValueError,UnicodeError,HTTPError,URLError,OSError,TimeoutError,KeyError,csv.Error) as e:
            item.update({"status":"SOURCE_INVALID_OR_BLOCKED","reason":f"{type(e).__name__}: {str(e)[:180]}"})
        receipt["files"].append(item)
    receipt["all_expected_files_verified"]=set(EXPECTED).issubset(successes)
    receipt["status"]="SOURCE_FEASIBILITY_ONLY" if receipt["all_expected_files_verified"] else "SOURCE_INCOMPLETE"
    return receipt

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--receipt",type=Path,required=True)
    args=ap.parse_args()
    try:
        result=audit()
    except (ValueError,HTTPError,URLError,OSError,TimeoutError,json.JSONDecodeError) as e:
        result={"schema":"chocho_original_kyoto_butterfly_competition_source_v01",
                "status":"SOURCE_ACCESS_BLOCKED","reason":f"{type(e).__name__}: {str(e)[:180]}",
                "new_regrowth_timing_effect_estimated":False,"GEB_PR38_untouched":True}
    args.receipt.parent.mkdir(parents=True,exist_ok=True)
    args.receipt.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2,ensure_ascii=False))

if __name__=="__main__":
    main()
