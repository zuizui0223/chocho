#!/usr/bin/env python3
"""Recover 2022 peer-reviewed 127-plant monarch host-performance source tables.

No guessed taxonomy. Produces raw source table inventory and candidate classification
rows; a later reviewed crosswalk must match independently to pinned HOSTS/WCVP.
"""
from __future__ import annotations
import argparse,hashlib,json,re,time,urllib.error,urllib.parse,urllib.request
from pathlib import Path

DOI="10.1371/journal.pone.0269701"
PLOS="https://journals.plos.org/plosone/article/file?id="+DOI+".s001&type=supplementary"
DRYAD="https://datadryad.org/api/v2/datasets/doi%3A10.5061%2Fdryad.7wm37pvvp"
MAX_BYTES=30_000_000

def retrieve(url,max_bytes=MAX_BYTES):
    request=urllib.request.Request(url,headers={
      "Accept":"application/octet-stream,application/vnd.openxmlformats-officedocument.wordprocessingml.document,application/json,*/*",
      "User-Agent":"chocho-reproducible-butterfly-quality-audit/0.1 (github.com/zuizui0223/chocho)"
    })
    with urllib.request.urlopen(request,timeout=65) as res:
        data=res.read(max_bytes+1)
        if len(data)>max_bytes:raise RuntimeError("download too large")
        return data,{"status":res.status,"url":res.geturl(),"content_type":res.headers.get("Content-Type")}
def safe_fetch(url,tries=3):
    attempts=[]
    for i in range(tries):
        try:
            data,meta=retrieve(url)
            return data,meta,attempts
        except Exception as e:
            attempts.append({"attempt":i+1,"error":str(e)[:450]})
            if i+1<tries:time.sleep(min(15,3+4*i))
    return None,None,attempts
def table_inventory(path):
    from docx import Document
    document=Document(path)
    inventories=[]
    for i,t in enumerate(document.tables):
        rows=[[re.sub(r"\s+"," ",cell.text).strip() for cell in row.cells] for row in t.rows]
        inventories.append({
          "table_number":i+1,"row_count":len(rows),
          "max_columns":max((len(x) for x in rows),default=0),
          "first_rows":rows[:9],
          "rows":rows
        })
    return inventories
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--protocol-json",type=Path,required=True)
    p.add_argument("--output-dir",type=Path,required=True)
    a=p.parse_args()
    config=json.loads(a.protocol_json.read_text())
    if config["independent_study"]["doi"]!=DOI:raise RuntimeError("DOI protocol drift")
    outdir=a.output_dir;outdir.mkdir(parents=True,exist_ok=True)
    result={"schema":"chocho_monarch_independent_quality_source_v0.1",
       "status":"SOURCE_PROBE", "article_doi":DOI,
       "original_quality_study":"Greenstein Steele Taylor (2022)",
       "source_access":{},"tables":[],
       "quality_classification_matched":False,
       "no_imputed_taxon_names":True,"no_model_fitted":True}
    for name,url in (("plos_s1_docx",PLOS),("dryad_dataset_metadata",DRYAD)):
        data,meta,errors=safe_fetch(url)
        result["source_access"][name]={"url":url,"success":data is not None,"metadata":meta,"attempt_errors":errors}
        if data is None:continue
        ext=".docx" if name=="plos_s1_docx" else ".json"
        source=outdir/("greenstein_2022_"+name+ext)
        source.write_bytes(data)
        result["source_access"][name].update({"size":len(data),"sha256":hashlib.sha256(data).hexdigest()})
        if name=="plos_s1_docx":
            try:
                tables=table_inventory(source)
                result["table_summary"]=[{k:v for k,v in t.items() if k!="rows"} for t in tables]
                inventory=outdir/"greenstein_2022_s1_tables_full_v01.json"
                inventory.write_text(json.dumps(tables,indent=2,ensure_ascii=False)+"\n")
                result["table_count"]=len(tables)
                result["status"]="PLOS_SUPPLEMENT_EXTRACTED_REQUIRES_CLASSIFICATION_REVIEW"
            except Exception as e:
                result["source_access"][name]["parse_error"]=str(e)[:800]
        else:
            try:
                z=json.loads(data)
                result["dryad_version_links"]=z.get("_links",{})
                result["dryad_metadata_title"]=z.get("title")
            except Exception as e:result["source_access"][name]["parse_error"]=str(e)[:500]
    if not result.get("table_summary"):
        result["status"]="NO_INDEPENDENT_CLASSIFICATION_SOURCE_EXTRACTED"
    (outdir/"monarch_quality_source_gate_result_v01.json").write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps({"status":result["status"],
        "plos_source":result["source_access"]["plos_s1_docx"]["success"],
        "dryad_metadata":result["source_access"]["dryad_dataset_metadata"]["success"],
        "table_count":result.get("table_count",0),
        "tables":result.get("table_summary",[])[:10],
        "error_summary":{k:v["attempt_errors"][-1:] for k,v in result["source_access"].items()}},ensure_ascii=False,indent=2),flush=True)
if __name__=="__main__":main()
