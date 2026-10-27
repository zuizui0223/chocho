#!/usr/bin/env python3
"""Public DOI-identified South African plant import interception source schema audit.

Figshare COLLECTION DOI 10.6084/m9.figshare.c.5180681; source original paper
Scientific Data 2021 Saccaggi et al DOI 10.1038/s41597-021-00869-z.
Pre-analysis source schema only. Confidential importer/exporter record IDs,
raw sample rows and exact identifiers are not published in the output.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,time
from collections import Counter
from pathlib import Path
from urllib.error import URLError,HTTPError
from urllib.request import Request,urlopen

SCHEMA="chocho_southafrica_plant_import_Lepidoptera_source_schema_20261027_v01"
PUBLIC_COLLECTION_ID=5180681
PUBLIC_API=f"https://api.figshare.com/v2/collections/{PUBLIC_COLLECTION_ID}/articles"
EXPECTED_FILES=("List of contaminants","Metadata of contaminants")


def source_get(url,maximum=60_000_000):
    errors=[]
    for i in range(3):
        req=Request(url,headers={
            "User-Agent":"chocho-invasive-insect-pathways-research/0.1",
            "Accept":"application/json,text/csv,*/*",
        })
        try:
            with urlopen(req,timeout=100) as r:
                if r.status!=200:
                    raise RuntimeError("HTTP "+str(r.status))
                data=r.read(maximum+1)
                if len(data)>maximum:
                    raise ValueError("Source exceeds bounded original public ZIP/CSV size")
                return data,r.url,r.headers.get("Content-Type")
        except (HTTPError,URLError,TimeoutError,RuntimeError) as e:
            errors.append(str(e)[:150])
            if i<2:time.sleep(2**i)
    raise RuntimeError("Figshare public source fetch error: "+str(errors))


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def field_values(path,raw):
    # Original documentation explicitly says semicolon-delimited sheets.
    sample=raw[:4096].decode("utf-8-sig",errors="replace")
    delimiter=";" if sample.splitlines()[0].count(";") > sample.splitlines()[0].count(",") else ","
    text=raw.decode("utf-8-sig",errors="replace")
    lines=csv.DictReader(io.StringIO(text),delimiter=delimiter)
    head=lines.fieldnames or []
    if len(set(head))!=len(head):
        raise ValueError("Ambiguous duplicated original sample/contaminant field names")
    stats=Counter()
    probes={k:Counter() for k in head
            if any(q in k.lower() for q in (
                "commodity","sample.source","lepidopter","status",
                "taxon","sample.type","order","class","kind"
            ))}
    for r in lines:
        stats["source_rows"]+=1
        for k in probes:
            val=str(r.get(k) or "").strip()
            if val:
                probes[k][val]+=1
        if any(r.get(k) is None for k in head):
            stats["ragged_rows"]+=1
    report={
        "source_original_filename":path,
        "original_source_member_sha256":sha(raw),
        "original_source_file_size_bytes":len(raw),
        "original_delimiter":delimiter,
        "column_count":len(head),
        "columns":head,
        "source_row_count":stats["source_rows"],
        "source_ragged_rows":stats["ragged_rows"],
        "selected_original_raw_column_vocabulary":{
            key:{
                "distinct_nonempty_values":len(cnt),
                "missing_or_empty_count":stats["source_rows"]-sum(cnt.values()),
                "top8_values":[{"value":v[:110],"n":n} for v,n in cnt.most_common(8)]
            }
            for key,cnt in probes.items()
            if len(cnt)<20000
        },
    }
    return report


def run():
    raw,url,_=source_get(PUBLIC_API)
    obj=json.loads(raw)
    articles=obj if isinstance(obj,list) else obj.get("articles",obj.get("items",[]))
    if not isinstance(articles,list) or not articles:
        raise ValueError("Original Figshare collection c.5180681 has no articles")
    member_articles=[]
    files=[]
    for article in articles:
        article_id=article.get("id")
        if not isinstance(article_id,int):
            raise ValueError("Invalid original Figshare collection member article id")
        raw_details,_,_=source_get(f"https://api.figshare.com/v2/articles/{article_id}")
        details=json.loads(raw_details)
        member_articles.append({
            "figshare_article_ID":article_id,
            "article_title":details.get("title"),
            "article_original_file_names":[v.get("name") for v in details.get("files",[])],
        })
        for item in details.get("files",[]):
            if str(item.get("name","")).lower().endswith(".csv"):
                files.append({**item,"origin_article_ID":article_id,
                              "origin_article_title":details.get("title")})
    if len(files)<2:
        raise ValueError("Original source expects two CSV files")
    outputs=[]
    for prefix in EXPECTED_FILES:
        chosen=[f for f in files if f["name"].lower().startswith(prefix.lower())]
        if len(chosen)!=1:
            raise ValueError("Original published pair CSV missing or ambiguous: "+prefix+"; actual original Figshare CSV names: "+repr([f["name"] for f in files]))
        item=chosen[0]
        payload,actual,mime=source_get(item["download_url"])
        if item.get("size") and len(payload)!=item["size"]:
            raise ValueError("Original Figshare listed file size mismatch "+prefix)
        if item.get("computed_md5") and hashlib.md5(payload).hexdigest() != item["computed_md5"]:
            raise ValueError("Original Figshare source digest mismatch "+prefix)
        info=field_values(item["name"],payload)
        info.update({
            "original_figshare_file_id":item["id"],
            "original_source_article_ID":item["origin_article_ID"],
            "original_source_article_title":item["origin_article_title"],
            "original_figshare_file_download_url":item["download_url"],
            "effective_download_url":actual,
            "original_figshare_file_md5":item.get("computed_md5") or item.get("supplied_md5"),
            "content_type":mime,
        })
        outputs.append(info)
    return {
        "schema":SCHEMA,
        "status":"PUBLIC_FIGSHARE_CONSIGNMENT_AND_CONTAMINANT_SOURCE_SCHEMAS_COMPLETE",
        "original_public_paper_DOI":"10.1038/s41597-021-00869-z",
        "original_public_dataset_DOI":"10.6084/m9.figshare.c.5180681",
        "figshare_collection_members_SHA256":sha(raw),
        "figshare_collection_ID":PUBLIC_COLLECTION_ID,
        "figshare_collection_article_inventory":member_articles,
        "metadata_download_API":PUBLIC_API,
        "original_files":outputs,
        "grain_warnings":[
            "Original first sheet rows are inspected plant samples and contain 79 taxon results; NA means not tested, 0 means tested negative.",
            "Original second sheet rows are taxonomic descriptions/statuses; these are not plant import shipments.",
            "Whether a contaminant has later established is not a causal fraction of interceptions and may be absent at species or Lepidoptera level.",
            "Intervention, audit and post-entry quarantine differ in source sampling design; do not pool them into unbiased invasion arrivals.",
            "No raw sample importer/exporter identifiers, full source names in each consignment, or original record matrix exported."
        ],
    }


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()
    report=run()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2,ensure_ascii=False,sort_keys=True)+"\n")
    print(json.dumps({
        "status":report["status"],
        "source_articles":report["figshare_collection_article_inventory"],
        "source_sheets":report["original_files"]
    },sort_keys=True,ensure_ascii=False))
