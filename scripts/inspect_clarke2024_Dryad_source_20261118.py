#!/usr/bin/env python3
"""Inspect and source-lock Clarke 2024 Dryad butterfly foodplant CSV sources.

This step retrieves published source data directly from the Dryad REST API.
It does not assume the correct foreign keys before inspecting the original CSV
headers; no plant/butterfly association or rank is invented if the API changes.
Data files are kept in ephemeral CI artifacts with per-file SHA256 receipts.
"""
import argparse,csv,hashlib,json,re,time,urllib.parse,urllib.request,urllib.error
from pathlib import Path

API="https://datadryad.org/api/v2"
DOI="doi:10.5061/dryad.1vhhmgr12"
DATASET=API+"/datasets/"+urllib.parse.quote(DOI,safe="")
DESIRED={
 "butterfly_synonym.csv","larval_host_plant.csv","larval_host_plant_status.csv",
 "larval_plants.csv","wcvp_extendedb.csv","reference.csv","region.csv","README.md"
}

def request_json(url):
    last=None
    for retry in range(4):
        try:
            req=urllib.request.Request(url,headers={
              "Accept":"application/json",
              "User-Agent":"chocho-independent-butterfly-Dryad-2024-source-audit"
            })
            with urllib.request.urlopen(req,timeout=50) as f:
                content=f.read(5_000_000+1)
                if len(content)>5_000_000: raise ValueError("Dryad JSON too large")
                data=json.loads(content)
            return data
        except (OSError,ValueError,urllib.error.HTTPError,urllib.error.URLError) as e:
            last=str(e);time.sleep(2*(retry+1))
    raise RuntimeError("Source Dryad API unavailable: "+str(last))

def resource_url(href):
    if not href: return None
    if href.startswith("https://datadryad.org"): return href
    if href.startswith("http"): raise ValueError("Unexpected API target "+href)
    return "https://datadryad.org/"+href.lstrip("/")

def discover_version():
    dataset=request_json(DATASET)
    links=dataset.get("_links",{})
    latest=links.get("stash:latestVersion",{}) or links.get("stash:version",{})
    if isinstance(latest,dict):
        href=latest.get("href")
        if href:
            obj=request_json(resource_url(href))
            number=obj.get("versionNumber")
            identifier=obj.get("id")
            if not identifier:
                identifier=href.rstrip("/").split("/")[-1]
            if str(identifier).isdigit(): return dataset,obj,int(identifier)
    # Explicit versions endpoint and source-ranked max published version.
    versions=request_json(DATASET+"/versions")
    items=versions.get("_embedded",{}).get("stash:versions",[])
    if not items and isinstance(versions,dict):
        items=versions.get("versions",[])
    if not items: raise ValueError("Dryad version list missing")
    for item in items:
        if not item.get("id"):
            href=(item.get("_links",{}).get("self") or {}).get("href","")
            suffix=href.rstrip("/").split("/")[-1]
            if suffix.isdigit():item["id"]=int(suffix)
    acceptable=[x for x in items if x.get("id")]
    if not acceptable:
        print("DRYAD_DIAGNOSTIC_DATASET_LINKS="+json.dumps(
          {"identifier":dataset.get("identifier"),"links":dataset.get("_links"),"keys":list(dataset)},ensure_ascii=False),flush=True)
        print("DRYAD_DIAGNOSTIC_VERSIONS_STRUCTURE="+json.dumps({
          "keys":list(versions) if isinstance(versions,dict) else "not_object",
          "embedded_keys":list(versions.get("_embedded",{})) if isinstance(versions,dict) else [],
          "first_versions":items[:2],
        },ensure_ascii=False)[:5000],flush=True)
        raise ValueError("No Dryad version IDs; see API structure diagnostic")
    obj=max(acceptable,key=lambda x:(int(x.get("versionNumber") or 0),int(x["id"])))
    return dataset,obj,int(obj["id"])

def fetch_files(version_id):
    url=API+"/versions/"+str(version_id)+"/files?per_page=100"
    objects=[]
    visited=set()
    while url:
        if url in visited or len(visited)>5: raise ValueError("Cyclic Dryad pagination")
        visited.add(url)
        block=request_json(url)
        items=block.get("_embedded",{}).get("stash:files",[])
        if not items and isinstance(block.get("files"),list): items=block["files"]
        objects.extend(items)
        href=(block.get("_links",{}).get("next") or {}).get("href")
        url=resource_url(href) if href else None
    if not objects: raise ValueError("Dryad source version has zero listed files")
    return objects

def download(uri,output):
    req=urllib.request.Request(uri,headers={"User-Agent":"chocho-source-pinned-Dryad-2024"})
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req,timeout=120) as response:
                payload=response.read(15_000_000)
            if not payload: raise ValueError("Empty Dryad source file")
            output.write_bytes(payload)
            return payload
        except (OSError,urllib.error.HTTPError,urllib.error.URLError,ValueError):
            if attempt==2:raise
            time.sleep(2*(attempt+1))

def main(dest, metadata_only=False):
    dest.mkdir(parents=True,exist_ok=True)
    dataset,version,version_id=discover_version()
    files=fetch_files(version_id)
    manifest={
      "doi":DOI,"dataset_api":DATASET,
      "dryad_dataset_id":dataset.get("id"),
      "dryad_version_id":version_id,
      "dryad_version_number":version.get("versionNumber"),
      "file_inventory":[{"filename":i.get("path") or i.get("filename"),"id":i.get("id"),"size":i.get("size")} for i in files],
      "retrieved":{}
    }
    available={Path(str(x.get("path") or x.get("filename") or "")).name:x for x in files}
    print("DRYAD_SOURCE_FILES_AVAILABLE="+json.dumps(sorted(available),ensure_ascii=False))
    required={"butterfly_synonym.csv","larval_host_plant.csv","larval_host_plant_status.csv","larval_plants.csv"}
    if not required<=set(available):
        raise ValueError("Missing necessary original Dryad source files: "+str(required-set(available)))
    if metadata_only:
        manifest["status"]="PUBLIC_METADATA_ONLY_NO_SOURCE_ROWS_RETRIEVED"
        manifest["download_attempt_status"]="Prior source-identity matched run 44204959410 returned HTTP 401 Unauthorized on file bytes; none may be analyzed without authorized access"
        manifest["caution"]="No butterfly-host relationships or evidence ranks were extracted or analyzed from this DOI"
        (dest/"SOURCE_MANIFEST.json").write_text(
            json.dumps(manifest,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
        print("DRYAD_METADATA_ONLY_NO_FILE_BYTES_RETRIEVED",len(files),flush=True)
        return manifest
    for filename in sorted(DESIRED & set(available)):
        record=available[filename]
        href=(record.get("_links",{}).get("stash:download") or {}).get("href")
        if href:
            uri=resource_url(href)
        else:
            ident=record.get("id")
            if not ident:raise ValueError("Missing file id/URL "+filename)
            uri=API+"/files/"+str(ident)+"/download"
        data=download(uri,dest/filename)
        sha=hashlib.sha256(data).hexdigest()
        info={"file_id":record.get("id"),"download_url":uri,"size":len(data),"sha256":sha}
        if filename.endswith(".csv"):
            with (dest/filename).open(encoding="utf-8-sig",newline="") as f:
                reader=csv.DictReader(f)
                fields=list(reader.fieldnames or [])
                head=[dict(x) for _,x in zip(range(3),reader)]
            info["headers"]=fields
            info["sample_records"]=head
        manifest["retrieved"][filename]=info
        print("DRYAD_SOURCE_HEADER="+json.dumps({"file":filename,"sha256":sha,
                "size":len(data),"headers":info.get("headers"),
                "sample":info.get("sample_records",[])},ensure_ascii=False),flush=True)
    (dest/"SOURCE_MANIFEST.json").write_text(
       json.dumps(manifest,ensure_ascii=False,indent=2)+"\n",encoding="utf-8"
    )
    print("DRYAD_2024_SOURCE_INSPECTION_SUCCESS")
    return manifest

if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir",type=Path,required=True)
    ap.add_argument("--metadata-only",action="store_true")
    args=ap.parse_args()
    main(args.out_dir,args.metadata_only)
