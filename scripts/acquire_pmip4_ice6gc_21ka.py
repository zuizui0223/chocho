#!/usr/bin/env python3
"""Fetch the authoritative PMIP4 ICE-6G_C 21ka/0ka land+ice source if accessible.

Fail closed: never silently substitute modern bathymetry, generic present-day ice
or static sea-level lowering for reconstructed 21ka land and ice.
"""
from __future__ import annotations
import argparse,hashlib,json,os,tempfile
from pathlib import Path
from urllib.parse import quote
import requests

MAX_SOURCE_BYTES=250_000_000

def candidates(token,filename):
    fname=quote(filename,safe="")
    base="https://sharebox.lsce.ipsl.fr"
    return [
        ("owncloud_public_selection",
         f"{base}/index.php/s/{token}/download?path=%2F&files={fname}",None),
        ("owncloud_public_webdav",
         f"{base}/public.php/webdav/{fname}",(token,"")),
        ("owncloud_dav_public_files",
         f"{base}/remote.php/dav/public-files/{token}/{fname}",None),
    ]

def fetch_exact(session,token,filename,md5,out_dir):
    out_dir=Path(out_dir);out_dir.mkdir(parents=True,exist_ok=True)
    errors=[]
    for name,url,auth in candidates(token,filename):
        temporary=None
        try:
            with session.get(url,auth=auth,timeout=(15,60),stream=True,allow_redirects=True) as response:
                response.raise_for_status()
                ctype=response.headers.get("Content-Type","").lower()
                if "text/html" in ctype or "application/json" in ctype:
                    raise RuntimeError(f"received {ctype} rather than a NetCDF file")
                header_size=response.headers.get("Content-Length")
                if header_size and int(header_size)>MAX_SOURCE_BYTES:
                    raise RuntimeError("PMIP source exceeds declared source cap")
                fd,temporary=tempfile.mkstemp(dir=out_dir,prefix=".partial_pmip_",suffix=".nc")
                m=hashlib.md5();sha=hashlib.sha256();count=0
                with os.fdopen(fd,"wb") as file:
                    for chunk in response.iter_content(chunk_size=1024*1024):
                        if not chunk:continue
                        count+=len(chunk)
                        if count>MAX_SOURCE_BYTES:raise RuntimeError("PMIP source cap exceeded")
                        m.update(chunk);sha.update(chunk);file.write(chunk)
                if m.hexdigest()!=md5:
                    raise RuntimeError(f"MD5 mismatch actual={m.hexdigest()} expected={md5} bytes={count}")
                if count<1_000_000:raise RuntimeError("implausibly small PMIP global NetCDF")
                target=out_dir/filename;os.replace(temporary,target)
                return {"filename":filename,"retrieval_method":name,"canonical_source_url":url,
                        "md5":m.hexdigest(),"sha256":sha.hexdigest(),"byte_count":count,
                        "source_url_after_redirect":response.url}
        except (requests.RequestException,ValueError,RuntimeError,OSError) as exc:
            errors.append({"method":name,"error":type(exc).__name__+": "+str(exc)[:240]})
        finally:
            if temporary and os.path.exists(temporary):
                os.unlink(temporary)
    raise RuntimeError(json.dumps({"filename":filename,"errors":errors},sort_keys=True))

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--contract",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    spec=json.loads(args.contract.read_text())
    if spec["status"]!="SOURCE_ONLY_PRE_BARRIER_RESULT" or len(spec["time_files"])!=2:
        raise ValueError("PMIP4 source contract status drift")
    args.out.mkdir(parents=True,exist_ok=True)
    session=requests.Session()
    session.headers.update({"User-Agent":"chocho-21ka-native-host-connectivity-physics/0.1 research reproducibility"})
    source_receipts=[];issues=[]
    order=sorted(spec["time_files"],key=lambda x:-x["age_ka"])
    for file in order:
        try:
            result=fetch_exact(session,"1RSwVf7afj38cmG",file["filename"],file["expected_md5"],args.out)
            source_receipts.append(result)
            print(json.dumps({"status":"SOURCE_MD5_VERIFIED","filename":file["filename"],
                 "sha256":result["sha256"],"size":result["byte_count"]}),flush=True)
        except RuntimeError as exc:
            issues.append({"filename":file["filename"],"failure":str(exc)})
            print(json.dumps({"status":"SOURCE_NOT_AVAILABLE","filename":file["filename"],
                 "error":str(exc)}),flush=True)
            break
    receipt={"schema":"chocho_ice6gc_source_retrieval_v01","verified_sources":source_receipts,
             "unresolved":issues,"complete":len(source_receipts)==2 and not issues,
             "PMIP4_source":"https://pmip4.lsce.ipsl.fr/doku.php/data:ice_ice6g_c",
             "genetic_or_host_LGM_model_opened":False,
             "unverified_substitute_fabricated":False}
    (args.out/"source_receipt.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    if not receipt["complete"]:
        raise SystemExit("Required authoritative 21ka source unavailable; stop before connectivity scoring")
if __name__=="__main__":main()
