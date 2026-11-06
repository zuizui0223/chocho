#!/usr/bin/env python3
"""Read ORIGINAL World Bank / IMF 2015–2021 commercial AIS density ZIP
central-directory metadata using HTTP byte ranges, without downloading
480 MB or mistaking AIS positions for ship voyages or plant transport.

This preparatory audit never fabricates a raster resistance field. Inspect
member compression, uncompressed byte size, raster name and licence before
planning reproducible 51-island georeferenced calculations.
"""
from __future__ import annotations
import argparse,hashlib,io,json,time,urllib.request,urllib.error,zipfile
from pathlib import Path

URL="https://datacatalogfiles.worldbank.org/ddh-published/0037580/5/DR0045405/shipdensity_commercial_.zip"
EXPECTED_LENGTH=480353369

def head():
    req=urllib.request.Request(URL,method="HEAD",headers={"User-Agent":"chocho-2026-island51-AIS-archive-source-audit/1.0"})
    with urllib.request.urlopen(req,timeout=35) as r:
        return {k:v for k,v in r.headers.items()}

class RemoteRangeZip(io.RawIOBase):
    """Minimal exact remote read/seek; strict HTTP 206 verification."""
    def __init__(self,url,size):
        self.url=url;self.size=size;self.cursor=0
        self.requests=[]
    def readable(self):return True
    def seekable(self):return True
    def seek(self,offset,whence=0):
        if whence==0:pos=offset
        elif whence==1:pos=self.cursor+offset
        elif whence==2:pos=self.size+offset
        else:raise ValueError("Invalid seek mode")
        if pos<0 or pos>self.size:raise ValueError("Requested ZIP offset outside official archive")
        self.cursor=pos
        return pos
    def tell(self):return self.cursor
    def read(self,n=-1):
        if n<0:n=self.size-self.cursor
        n=min(n,self.size-self.cursor)
        if n<=0:return b""
        start=self.cursor;end=start+n-1
        req=urllib.request.Request(self.url,headers={
            "Range":f"bytes={start}-{end}",
            "User-Agent":"chocho-2026-island51-AIS-archive-source-audit/1.0"})
        error=None
        for attempt in range(4):
            try:
                with urllib.request.urlopen(req,timeout=55) as r:
                    if r.status!=206:
                        raise ValueError(f"Expected 206 byte-range, got HTTP {r.status}")
                    received=r.read(n+1)
                    if len(received)!=n:
                        raise ValueError(f"Truncated byte range {start}-{end}: {len(received)}")
                    cr=r.headers.get("Content-Range","")
                    if not cr.startswith(f"bytes {start}-{end}/"):
                        raise ValueError("Official response has wrong Content-Range "+cr)
                self.cursor+=n
                self.requests.append({"start":start,"end":end,"sha256":hashlib.sha256(received).hexdigest()})
                return received
            except Exception as e:
                error=repr(e)[:300];time.sleep(attempt+1)
        raise RuntimeError("HTTP range ZIP central-directory read failed: "+str(error))

def run():
    h=head()
    size=int(h.get("Content-Length","0"))
    if size!=EXPECTED_LENGTH or "bytes" not in h.get("Accept-Ranges","").lower():
        raise ValueError(f"World Bank AIS source changed: content length {size}, accept ranges {h.get('Accept-Ranges')}")
    stream=RemoteRangeZip(URL,size)
    with zipfile.ZipFile(stream,allowZip64=True) as z:
        member=[]
        for v in z.infolist():
            member.append({
              "filename":v.filename,
              "uncompressed_bytes":v.file_size,
              "archive_compressed_bytes":v.compress_size,
              "compress_type":int(v.compress_type),
              "crc32_hex":hex(v.CRC),
              "is_dir":v.is_dir()
            })
        bad=z.testzip() if False else None  # Never fetch/decompress entire 480MB archive for audit
    return {
       "schema":"chocho_WorldBank_IMF_commercial_AIS_ZIP_HTTP_range_manifest_20261106_v1",
       "status":"PUBLIC_OFFICIAL_AIS_ZIP_CENTRAL_DIRECTORY_VERIFIED_NO_RASTER_VALUES_EXTRACTED",
       "catalogue":"https://datacatalog.worldbank.org/search/dataset/0037580/global-shipping-traffic-density",
       "archive_url":URL,
       "source_origin":"IMF hourly AIS broadcasts 2015-01 through 2021-02",
       "source_units":"AIS positions reported per 0.005 by 0.005 degree cell, NOT ship voyages or inter-island cargo quantity",
       "source_licence":"World Bank published CC BY 4.0",
       "source_zip_bytes":size,
       "HTTP_headers":{k:h.get(k) for k in ("Content-Type","Content-Length","Accept-Ranges","Last-Modified")},
       "central_directory_members":member,
       "HTTP_range_requests":stream.requests,
       "HTTP_byte_range_download_total_bytes":int(sum(x["end"]-x["start"]+1 for x in stream.requests)),
       "actual_raster_pixels_decoded":0,
       "archive_was_not_fully_downloaded":True,
       "go_forward_gate":"Need verified original 51 WPI port coordinates, raster extraction metadata and ocean mask before true AIS resistance matrices.",
       "not_evidence":"Modern AIS density cannot establish 19th-century colonial plant exchanges, cargo content or butterfly establishment."
    }

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    d=run()
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(d,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print("WORLDBANK_COMMERCIAL_AIS_2015_2021_ORIGINAL_ZIP_MANIFEST",json.dumps({
      "status":d["status"],
      "total_http_bytes":d["HTTP_byte_range_download_total_bytes"],
      "members":d["central_directory_members"]
    },ensure_ascii=False,sort_keys=True),flush=True)
if __name__=="__main__":main()
