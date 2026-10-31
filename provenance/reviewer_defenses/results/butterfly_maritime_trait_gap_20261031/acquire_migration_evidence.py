#!/usr/bin/env python3
"""Acquire/check a published independently documented butterfly migration roster.

The review's Table S1 enumerates documented migrant species; absence from the
table DOES NOT establish non-migration. No original supplement is redistributed.
The analysis must not equate migration with movement to an indexed island.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import re
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET
import zipfile
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parents[4]
HERE=Path(__file__).resolve().parent

PAPER="https://doi.org/10.1111/brv.12714"
FILENAME="brv12714-sup-0001-TableS1.xlsx"
# Each candidate must be checked by actual valid XLSX parse, not URL shape alone.
SOURCE_URLS=(
 "https://onlinelibrary.wiley.com/doi/supinfo/10.1111/brv.12714/"+FILENAME,
 "https://onlinelibrary.wiley.com/doi/suppl/10.1111/brv.12714/supinfo/"+FILENAME,
 "https://onlinelibrary.wiley.com/store/10.1111/brv.12714/asset/supinfo/"+FILENAME,
 "https://onlinelibrary.wiley.com/doi/supinfo/10.1111/brv.12714/"+FILENAME+"?download=1"
)
INPUT="provenance/reviewer_defenses/results/island_host_network_20261030/island_host_network_all225.csv"
INPUT_SHA="39fe84981cab2b0e47063399f0a5c0c9731e800d"
GAPS="provenance/reviewer_defenses/results/butterfly_maritime_trait_gap_20261031/matched_species_screen.csv"
UA="chocho-species-migration-evidence-crosswalk/1.0 (+research citation)"
NS={"m":"http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
SHEET_NS="{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
BINOMIAL=re.compile(r"^\s*([A-Z][a-zA-Z-]+)\s+([a-z][a-zA-Z-]+)(?:\b|$)")
COLUMN=re.compile(r"^([A-Z]+)")
LABELS=("species","scientific name","binomial","taxon","butterfly name")

def source_csv():
    raw=(ROOT/INPUT).read_bytes()
    sha=hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
    if sha!=INPUT_SHA:
        raise ValueError("Fixed 225 butterfly scope changed: "+sha)
    rows=list(csv.DictReader(io.StringIO(raw.decode("utf-8"))))
    if len(rows)!=225 or len({r["species"] for r in rows})!=225:
        raise ValueError("Expected unique 225 candidate butterfly species")
    return rows

def latin(cell):
    if cell is None:
        return None
    cell=re.sub(r"\s+"," ",str(cell).strip())
    m=BINOMIAL.match(cell)
    if not m or m.group(2).lower() in {"sp","spp","cf","aff"}:
        return None
    return m.group(1)+" "+m.group(2)

def read_xlsx(data):
    if not zipfile.is_zipfile(io.BytesIO(data)):
        raise ValueError("Supplied bytes are not ZIP based XLSX")
    all_columns=[]
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        names=z.namelist()
        if "[Content_Types].xml" not in names or "xl/workbook.xml" not in names:
            raise ValueError("Wrong archive, not an XLSX workbook")
        shared=[]
        if "xl/sharedStrings.xml" in names:
            root=ET.fromstring(z.read("xl/sharedStrings.xml"))
            for si in root.findall("m:si",NS):
                shared.append("".join(t.text or "" for t in si.iter(SHEET_NS+"t")))
        for filename in sorted(p for p in names
                               if re.fullmatch(r"xl/worksheets/sheet\d+\.xml",p)):
            rows=[]
            root=ET.fromstring(z.read(filename))
            for row in root.iter(SHEET_NS+"row"):
                vals={}
                for cell in row.findall("m:c",NS):
                    a=cell.attrib.get("r","")
                    m=COLUMN.match(a)
                    if not m:continue
                    col=m.group(1)
                    typ=cell.attrib.get("t")
                    value=cell.find("m:v",NS)
                    if typ=="inlineStr":
                        inline=cell.find("m:is",NS)
                        s="".join(t.text or "" for t in inline.iter(SHEET_NS+"t")) if inline is not None else ""
                    elif value is None:
                        s=""
                    elif typ=="s":
                        idx=int(value.text)
                        s=shared[idx] if 0<=idx<len(shared) else ""
                    else:s=value.text or ""
                    vals[col]=s.strip()
                if vals:rows.append(vals)
            if len(rows)<100:continue
            for column in sorted(set().union(*(set(r) for r in rows))):
                values=[r.get(column,"") for r in rows]
                headers=values[:12]
                header_match=next((x for x in headers
                    if any(token in x.lower().replace("_"," ").strip() for token in LABELS)),None)
                taxa=[latin(x) for x in values]
                nameset={s for s in taxa if s is not None}
                all_columns.append({
                    "sheet_xml":filename,"column":column,"header":header_match,
                    "documented_names":sorted(nameset),"observed_taxa":len(nameset),
                    "row_count":len(rows),
                })
    eligible=[v for v in all_columns if v["header"] and v["observed_taxa"]>=100]
    if not eligible:
        raise ValueError("No >=100-species column with recognized header; manual source QA required")
    eligible.sort(key=lambda v:(-v["observed_taxa"],v["sheet_xml"],v["column"]))
    selected=eligible[0]
    return selected, [{k:v for k,v in c.items() if k!="documented_names"} for c in all_columns]

def acquire():
    attempts=[]
    for url in SOURCE_URLS:
        try:
            req=urllib.request.Request(url,headers={"User-Agent":UA,"Accept":"application/vnd.openxmlformats-officedocument.spreadsheetml.sheet,*/*"})
            with urllib.request.urlopen(req,timeout=25) as f:
                final=f.url
                b=f.read(500_000)
            if not(16_000<=len(b)<=400_000):
                raise ValueError("Unexpected XLSX byte size")
            selected,all_columns=read_xlsx(b)
            if len(selected["documented_names"])>1000:
                raise ValueError("Implausibly many species for ~568-review roster")
            return {"status":"VALIDATED_TABLE_ROSTER","source_url":final,
                    "source_sha256":hashlib.sha256(b).hexdigest(),
                    "source_bytes":len(b),"selection":selected,"columns_audit":all_columns,
                    "fetch_attempts":attempts}
        except (urllib.error.URLError,TimeoutError,ValueError,zipfile.BadZipFile,ET.ParseError) as exc:
            attempts.append({"url":url,"error_type":type(exc).__name__,
                             "message":str(exc)[:200]})
    return {"status":"SUPPLEMENT_UNAVAILABLE_OR_UNVERIFIED","fetch_attempts":attempts}

def summarize(source):
    all225=source_csv()
    report={"schema":"chocho_independent_published_butterfly_migration_roster_crosswalk_v1",
      "source_paper":PAPER,"supplement":FILENAME,"source_225_git_blob":INPUT_SHA,
      "status":source["status"],"source_url":source.get("source_url"),
      "source_sha256":source.get("source_sha256"),
      "source_bytes":source.get("source_bytes"),
      "fetch_attempts":source["fetch_attempts"],
      "original_supplement_not_redistributed":True,
      "not_listed_does_not_prove_nonmigratory":True,
      "wing_size_is_not_direct_migration_evidence":True,
      "no_causal_maritime_host_establishment_claim":True}
    if source["status"]!="VALIDATED_TABLE_ROSTER":
        return report,[]
    selected=source["selection"]
    recorded=set(selected["documented_names"])
    candidates={x["species"] for x in all225}
    rows=[{"species":x["species"],"speciesKey":x["taxon_key"],
        "migration_evidence":"DOCUMENTED_MIGRANT_REVIEW_TABLE_S1" if x["species"] in recorded
                             else "NOT_LISTED_STATUS_UNKNOWN"}
          for x in all225]
    matched=sorted(recorded & candidates)
    report.update({"selected_column":{k:v for k,v in selected.items() if k!="documented_names"},
      "audit_columns":source["columns_audit"],
      "review_documented_unique_species":len(recorded),
      "chocho_candidate_species":225,
      "exact_documented_species_in_225":len(matched),
      "candidate_not_listed_status_unknown":225-len(matched),
      "exact_match_method":"unmodified accepted binomial, no unresolved synonym substitutions"})
    gaps=list(csv.DictReader((ROOT/GAPS).open(encoding="utf-8",newline="")))
    if len(gaps)!=193 or len({x["species"] for x in gaps})!=193:
        raise ValueError("Expected pinned 193-species public gap diagnostic")
    d=[x for x in gaps if x["species"] in recorded]
    other=[x for x in gaps if x["species"] not in recorded]
    def stats(part):
        vals=[int(x["far500"])/int(x["n"]) for x in part]
        return {"species":len(part),"record_events":sum(int(x["n"]) for x in part),
         "record_gaps_500km":sum(int(x["far500"]) for x in part),
         "species_equal_mean_gap_fraction":sum(vals)/len(vals) if vals else None}
    report["record_gap_193_species_subsets"]={
      "documented_migrant":stats(d),
      "not_listed_unknown":stats(other)
    }
    return report,rows

def main():
    source=acquire()
    report,rows=summarize(source)
    out=HERE/"migration_roster_acquisition.json"
    out.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    if rows:
        with (HERE/"migration_roster_225_crosswalk.csv").open("w",newline="",encoding="utf-8") as f:
            w=csv.DictWriter(f,fieldnames=["species","speciesKey","migration_evidence"])
            w.writeheader();w.writerows(rows)
    print(json.dumps({k:report.get(k) for k in
          ("status","source_url","review_documented_unique_species",
           "exact_documented_species_in_225","record_gap_193_species_subsets")},indent=2))

if __name__=="__main__":
    main()
