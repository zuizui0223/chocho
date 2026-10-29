#!/usr/bin/env python3
"""Check all original BSBI 3497 per-taxon CSV date flag vocabularies.

SOURCE-ONLY source-level categorical counts; no plant taxa, hectads or individual
geographies exported. Original MD5 and ZIP path/source value gates enforced.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,zipfile
from collections import Counter
from pathlib import Path

from audit_bsbi_historical_hectad_source_20261029 import MD5,MAX_ARCHIVE,nested_source_metadata

CLASSES=("to 1929","1930 - 1969","1970 - 1986","1987 - 1999","2000 - 2009","2010 - 2019")
NAMELOOKUP="PA2020_10by10kmSquareData/nameLookup.csv"
NESTED="PA2020_10by10kmSquareData/hectadDistributions.zip"

def run(archive):
    raw=archive.read_bytes()
    if len(raw)>MAX_ARCHIVE or hashlib.md5(raw).hexdigest()!=MD5:
        raise ValueError("BSBI original public ZIP changed")
    result={"schema":"chocho_bsbi_original_all_taxa_date_flag_categories_20261029_v01",
            "source_zip_md5":MD5,"date_fields":list(CLASSES),
            "raw_hectad_plant_records_exported":False}
    with zipfile.ZipFile(io.BytesIO(raw)) as outer:
        if not {NESTED,NAMELOOKUP}<=set(outer.namelist()):
            raise ValueError("Original archive missing name lookup or nested 3497 taxa")
        name_rows=list(csv.DictReader(io.TextIOWrapper(outer.open(NAMELOOKUP),encoding="utf-8-sig",newline="")))
        if len(name_rows)!=3497:raise ValueError(f"Source taxon name lookup count differs: {len(name_rows)}")
        orig_name_cols=list(name_rows[0])
        nested=outer.read(NESTED)
    count=Counter()
    flags={field:Counter() for field in CLASSES}
    headers=Counter()
    source_hectad_patterns=Counter()
    with zipfile.ZipFile(io.BytesIO(nested)) as z:
        files=[f for f in z.infolist() if f.filename.endswith(".csv")]
        if len(files)!=3497:raise ValueError("Original 3497 species source ZIP changed")
        for idx,f in enumerate(files,1):
            with z.open(f) as handle:
                rows=csv.DictReader(io.TextIOWrapper(handle,encoding="utf-8-sig",errors="strict",newline=""))
                head=tuple(rows.fieldnames or [])
                headers[head]+=1
                if not set(CLASSES)<=set(head) or "hectad" not in head:
                    raise ValueError(f"Required original date class missing, file index {idx}")
                for row in rows:
                    count["total_source_hectad_taxon_rows"]+=1
                    # Only aggregate flag vocabularies; NEVER hectad identifier.
                    val=str(row["hectad"] or "").strip()
                    if len(val)==4 and val[:2].isalpha() and val[2:].isdigit():
                        source_hectad_patterns["2letters2digits"]+=1
                    elif val:
                        source_hectad_patterns["other_nonempty"]+=1
                    else:source_hectad_patterns["missing"]+=1
                    for field in CLASSES:
                        v=str(row.get(field) or "").strip()
                        flags[field][v[:45] if v else "[empty]"]+=1
        count["nested_original_source_csv_taxa"]=len(files)
    result["source_total_hectad_taxon_rows"]=count["total_source_hectad_taxon_rows"]
    result["original_species_lookup_rows"]=len(name_rows)
    result["original_lookup_fields"]=orig_name_cols
    result["original_csv_header_profiles"]=[
        {"fields":list(k),"n_taxon_files":v} for k,v in headers.items()
    ]
    result["source_hectad_format_counts"]=dict(sorted(source_hectad_patterns.items()))
    result["original_date_field_value_vocabularies"]={k:dict(sorted(v.items())) for k,v in flags.items()}
    result["status"]="ALL_ORIGINAL_3497_TAXA_PERIOD_FLAG_SOURCE_VOCABULARY_VERIFIED"
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--original-zip",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    a=p.parse_args()
    d=run(a.original_zip)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(d,indent=2,ensure_ascii=False,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({"status":d["status"],"original_lookup_fields":d["original_lookup_fields"],
                      "source_hectad_format_counts":d["source_hectad_format_counts"],
                      "date_values":d["original_date_field_value_vocabularies"],
                      "species_files":d["original_species_lookup_rows"],
                      "source_taxon_hectad_records":d["source_total_hectad_taxon_rows"]},ensure_ascii=False))
