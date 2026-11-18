#!/usr/bin/env python3
"""Check publisher/PMC Appendix S1 of Clarke 2024 as an independent evidence source.

Never treat the previous 401 Dryad dataset inventory as data access.
This script downloads actual licensed journal supplementary XLSX bytes or
reports source failure; original SHA256 and raw sheet/row text extracts
are retained. Reads Open XML using Python stdlib zipfile/XML only; no
external XLSX writers, pandas, or openpyxl. This audit does NOT change
the journal manuscript, HOSTS data, or original 239-butterfly cohort.
"""
import csv
import hashlib
import json
import re
import time
import urllib.error
import urllib.request
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT=Path("provenance/reviewer_defenses/results/independent_island_butterfly_host_coverage_20261117")
INBO=ROOT/"independent_42_species_vs_frozen239.csv"
BIBLIO={"article":"Clarke 2024, Ecology and Evolution e10834",
        "article_doi":"10.1002/ece3.10834",
        "supplement":"ece310834-sup-0001-AppendixS1.xlsx"}
CANDIDATES=[
 "https://pmc.ncbi.nlm.nih.gov/articles/instance/10771928/bin/ece310834-sup-0001-AppendixS1.xlsx",
 "https://pmc.ncbi.nlm.nih.gov/articles/PMC10771928/bin/ece310834-sup-0001-AppendixS1.xlsx",
 "https://europepmc.org/articles/PMC10771928/bin/ece310834-sup-0001-AppendixS1.xlsx",
 "https://onlinelibrary.wiley.com/doi/suppl/10.1002/ece3.10834/supinfo/ece310834-sup-0001-AppendixS1.xlsx"
]
NS={"x":"http://schemas.openxmlformats.org/spreadsheetml/2006/main",
    "r":"http://schemas.openxmlformats.org/officeDocument/2006/relationships",
    "p":"http://schemas.openxmlformats.org/package/2006/relationships"}

def fetch(url):
    req=urllib.request.Request(url,headers={
       "Accept":"application/vnd.openxmlformats-officedocument.spreadsheetml.sheet,application/octet-stream,*/*",
       "User-Agent":"Mozilla/5.0 (research provenance fetch; referenced publicly accessible journal supplement)"
    })
    with urllib.request.urlopen(req,timeout=40) as f:
        data=f.read(6_000_000)
        status=f.status
        ctype=f.headers.get("content-type")
        final_url=f.url
    return data,{"status":status,"content_type":ctype,"final_url":final_url,"length":len(data)}

def text_from_si(node):
    return "".join((t.text or "") for t in node.findall(".//x:t",NS))

def read_shared(z):
    if "xl/sharedStrings.xml" not in z.namelist():
        return []
    root=ET.fromstring(z.read("xl/sharedStrings.xml"))
    return [text_from_si(node) for node in root.findall("x:si",NS)]

def row_cells(row,strings):
    res={}
    for c in row.findall("x:c",NS):
        ref=c.attrib.get("r","")
        col=re.match(r"[A-Z]+",ref)
        if not col:continue
        typ=c.attrib.get("t","")
        if typ=="inlineStr":
            value=text_from_si(c.find("x:is",NS))
        else:
            v=c.find("x:v",NS)
            value="" if v is None else (v.text or "")
            if typ=="s" and value.isdigit():
                i=int(value)
                value=strings[i] if i<len(strings) else ""
        res[col.group()]=str(value)
    return res

def scan_xlsx(data,candidates):
    out=[]
    with zipfile.ZipFile(__import__("io").BytesIO(data)) as z:
        names=set(z.namelist())
        if "[Content_Types].xml" not in names or "xl/workbook.xml" not in names:
            raise ValueError("Not a real XLSX spreadsheet")
        shared=read_shared(z)
        sheet_xml=sorted((x for x in names if re.fullmatch(r"xl/worksheets/sheet\d+\.xml",x)),
            key=lambda x:int(re.search(r"(\d+)\.xml",x).group(1)))
        if not sheet_xml:raise ValueError("No XLSX worksheets found")
        for fn in sheet_xml:
            root=ET.fromstring(z.read(fn))
            rows=root.findall(".//x:sheetData/x:row",NS)
            examples=[]
            matched=[]
            for i,row in enumerate(rows):
                cells=row_cells(row,shared)
                if i<6:
                    examples.append(cells)
                joined=" ".join(x for x in cells.values() if x).lower()
                found=[name for name in candidates if name.lower() in joined]
                if found:
                    matched.append({"row":int(row.attrib.get("r",i+1)),
                                   "butterfly_names":found,"cells":cells})
            out.append({"sheet_xml":fn,"rows":len(rows),"top_rows":examples,
                        "matched_rows_n":len(matched),
                        "matched_rows":matched[:800]})
    return out

def main(output):
    output.mkdir(parents=True,exist_ok=True)
    with INBO.open(encoding="utf-8-sig",newline="") as f:
        source=list(csv.DictReader(f))
    if len(source)!=42:raise ValueError("Original independent 42 changed")
    names=sorted({x["scientific_name"] for x in source})
    audit={"schema":"chocho_clarke2024_publisher_supplement_source_probe_v1",
           "requested_source":BIBLIO,"source_urls_considered":CANDIDATES,
           "independent_taxa_n":len(names),"download_attempts":[],"file_status":"NOT_RETRIEVED",
           "source_rows_used_for_scientific_conclusions":0}
    data=None
    for url in CANDIDATES:
        try:
            raw,metadata=fetch(url)
            if not raw.startswith(b"PK"):
                raise ValueError("Downloaded body does not start with XLSX ZIP magic")
            if len(raw)<1000 or len(raw)>5_000_000:
                raise ValueError("Unreasonable original appendix size")
            data=raw
            metadata["sha256"]=hashlib.sha256(raw).hexdigest()
            audit["download_attempts"].append({"url":url,**metadata})
            audit["file_status"]="SOURCE_XLSX_BYTES_VERIFIED"
            audit["verified_original_url"]=metadata["final_url"]
            audit["original_xlsx_sha256"]=metadata["sha256"]
            break
        except Exception as exc:
            audit["download_attempts"].append({"url":url,"error":str(exc)[:250]})
            print("CLARKE_PUBLISHER_ATTEMPT_FAILED="+json.dumps(audit["download_attempts"][-1]),flush=True)
            time.sleep(0.3)
    if data is not None:
        xlsx_path=output/"source_clarke2024_appendixS1.xlsx"
        xlsx_path.write_bytes(data)
        sheets=scan_xlsx(data,names)
        audit["sheets"]=sheets
        audit["source_rows_used_for_scientific_conclusions"]=sum(x["matched_rows_n"] for x in sheets)
        print("CLARKE_2024_PUBLISHER_SOURCE_ROWS="+json.dumps({
          "sha256":audit["original_xlsx_sha256"],
          "sheets":[{"file":x["sheet_xml"],"rows":x["rows"],
                     "matching_original_independent_taxa_rows":x["matched_rows_n"],
                     "first_records":x["top_rows"][:3],
                     "first_matches":x["matched_rows"][:4]} for x in sheets]
        },ensure_ascii=False)[:8000],flush=True)
    else:
        audit["limitation"]="Publisher/PMC source downloads unavailable; no record-level host/evidence inference attempted"
        print("CLARKE_2024_PUBLISHER_SOURCE_UNAVAILABLE_NO_DATA",flush=True)
    (output/"appendixS1_source_audit.json").write_text(
       json.dumps(audit,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print("CLARKE_2024_PUBLISHER_SOURCE_PROBE_COMPLETED_"+audit["file_status"])
    return audit

if __name__=="__main__":
    import argparse
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out-dir",required=True,type=Path)
    args=p.parse_args()
    main(args.out_dir)
