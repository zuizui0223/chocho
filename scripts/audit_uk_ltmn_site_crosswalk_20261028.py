#!/usr/bin/env python3
"""Read OGL UK LTMN transect xlsx crosswalk metadata without exporting locations."""
from __future__ import annotations
import argparse,hashlib,io,json,re,zipfile
from pathlib import Path
from urllib.request import Request,urlopen
from xml.etree import ElementTree as ET

URL="https://publications.naturalengland.org.uk/file/6033934228652032"
SHA="e38853b775ecadf39332a8412e2a4a2fa3c83512cf7218e4f1eda172dd83c9cd"
N={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
NS=N['m']

def original_rows(body):
    if hashlib.sha256(body).hexdigest()!=SHA:raise ValueError("Original LTMN XLSX changed")
    with zipfile.ZipFile(io.BytesIO(body)) as z:
        if z.testzip() is not None:raise ValueError("Original OGL spreadsheet corrupted")
        strings=[]
        if "xl/sharedStrings.xml" in z.namelist():
            root=ET.fromstring(z.read("xl/sharedStrings.xml"))
            for si in root.findall('m:si',N):
                strings.append(''.join(v.text or '' for v in si.iter(f'{{{NS}}}t')))
        sheets=sorted(x for x in z.namelist() if re.fullmatch(r"xl/worksheets/sheet\d+\.xml",x))
        if not sheets:raise ValueError("No genuine worksheet XML")
        for path in sheets:
            root=ET.fromstring(z.read(path))
            out=[]
            for r in root.findall('.//m:sheetData/m:row',N):
                data={}
                for cell in r.findall("m:c",N):
                    ref=cell.get("r","")
                    col=re.match(r"[A-Z]+",ref)
                    if not col:continue
                    n=0
                    for ch in col.group(0):n=n*26+ord(ch)-64
                    typ=cell.get("t","")
                    obj=cell.find("m:v",N)
                    if obj is None and typ=="inlineStr":
                        inline=cell.find("m:is",N)
                        val=''.join(x.text or '' for x in inline.iter(f'{{{NS}}}t')) if inline is not None else ''
                    else:
                        val=obj.text if obj is not None else ''
                        if typ=="s" and str(val).isdigit():
                            val=strings[int(val)]
                    data[n]=val
                out.append(data)
            yield path,out

def inspect(rows):
    if not rows:return {"n_excel_rows":0}
    header=rows[0]
    fields={idx:str(name).strip()[:90] for idx,name in header.items()}
    # Report structural field names only, not any source row.
    return {"n_excel_rows_including_header":len(rows),
            "n_data_rows":len(rows)-1,
            "header_names":fields,
            "plausible_spatial_columns":[v for v in fields.values()
                  if re.search(r"easting|northing|coord|latitude|longitude|grid.?ref|post.?code|location",v,re.I)],
            "source_row_values_or_site_labels_exported":False}

def run(path):
    req=Request(URL,headers={"User-Agent":"chocho-LTMN-public-transect-crosswalk/1.0"})
    with urlopen(req,timeout=35) as r:body=r.read(350000)
    if len(body)>350000:raise ValueError("Original crosswalk unexpectedly large")
    out={"schema":"chocho_LTMN_original_site_coordinate_crosswalk_schema_20261028_v01",
          "source_sha256":hashlib.sha256(body).hexdigest(),
          "source_url":URL,"sheets":[{"sheet":p,**inspect(rows)} for p,rows in original_rows(body)],
          "original_site_ids_or_coordinates_exported":False}
    out["status"]="ORIGINAL_PUBLIC_CROSSWALK_FIELDS_VERIFIED"
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps(out,sort_keys=True))
    return out

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--output",required=True,type=Path)
    run(p.parse_args().output)
