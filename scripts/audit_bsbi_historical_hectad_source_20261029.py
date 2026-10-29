#!/usr/bin/env python3
"""Audit original public BSBI Plant Atlas 2020 ZIP before any butterfly endpoint join.

Source first: verify publisher file MD5, ZIP members and only source field names.
No individual plant × hectad rows, original monitoring site coordinates,
recorders or other raw data are copied into artifacts.
"""
from __future__ import annotations

import argparse,csv,hashlib,io,json,stat,zipfile
from pathlib import Path

MD5="66a8bb4305df484839587b5fe8acf1d2"
DOI="10.5281/zenodo.10998003"
FILENAME="PA2020_10by10kmSquareData.zip"
SOURCE_PAGE="https://zenodo.org/records/10998003"
MAX_ARCHIVE=18_000_000
MAX_UNPACKED=3_000_000_000
MAX_MEMBER=2_600_000_000

def audit_zip(path:Path):
    if not path.is_file() or path.stat().st_size>MAX_ARCHIVE:
        raise ValueError("Expected bounded original source archive missing or too large")
    digest=hashlib.md5(path.read_bytes()).hexdigest()
    if digest!=MD5:
        raise ValueError("Publisher file MD5 differs; no biological inference allowed")

def nested_source_metadata(nested_zip:bytes):
    """Safe source nested ZIP: inspect member types and source headers only."""
    from collections import Counter
    with zipfile.ZipFile(io.BytesIO(nested_zip)) as inner:
        members=[x for x in inner.infolist() if not x.is_dir()]
        if not members or sum(x.file_size for x in members)>MAX_UNPACKED:
            raise ValueError("Original nested source exceeds uncompressed size limit")
        suffix_count=Counter()
        header_examples=[]
        source_names=[]
        for item in members:
            p=Path(item.filename)
            if p.is_absolute() or ".." in p.parts or item.file_size>MAX_MEMBER:
                raise ValueError("Unsafe original inner ZIP member")
            suffix=p.suffix.lower()
            suffix_count[suffix]+=1
            if len(source_names)<20:source_names.append(item.filename[-180:])
            if suffix in (".csv",".tsv",".txt") and len(header_examples)<8:
                with inner.open(item) as handle:
                    txt=io.TextIOWrapper(handle,encoding="utf-8-sig",errors="replace",newline="")
                    line=txt.readline()
                    delim="\\t" if "\\t" in line and line.count("\\t")>line.count(",") else ","
                    header=next(csv.reader([line],delimiter=delim),[])
                    header_examples.append({"original_file_index":len(header_examples)+1,
                        "source_member_name":item.filename[-150:],
                        "header_fields":[str(v).strip()[:100] for v in header][:75],
                        "source_text_delimiter":delim,
                        "original_file_uncompressed_bytes":item.file_size})
        return {"nested_member_count":len(members),
            "nested_uncompressed_bytes":sum(x.file_size for x in members),
            "nested_suffix_counts":dict(sorted(suffix_count.items())),
            "nested_first_member_names":source_names,
            "nested_text_header_examples":header_examples}

    output=[]
    with zipfile.ZipFile(path) as z:
        if z.testzip() is not None:
            raise ValueError("Original BSBI ZIP CRC invalid")
        files=[f for f in z.infolist() if not f.is_dir()]
        if not files or sum(f.file_size for f in files)>MAX_UNPACKED:
            raise ValueError("Original archive empty or decompression size bound exceeded")
        for item in files:
            p=Path(item.filename)
            if p.is_absolute() or '..' in p.parts or item.file_size>MAX_MEMBER:
                raise ValueError("Unsafe original member path or size")
            mode=item.external_attr>>16
            if stat.S_ISLNK(mode):
                raise ValueError("No symbolic links in source")
            suffix=p.suffix.lower()
            m={"original_member_name":item.filename[-180:],
               "suffix":suffix,"compressed_bytes":item.compress_size,
               "original_uncompressed_bytes":item.file_size}
            if suffix==".zip" and item.file_size<=MAX_ARCHIVE:
                m["nested_original_source_metadata"]=nested_source_metadata(z.read(item))
            if suffix in (".csv",".tsv",".txt"):
                with z.open(item) as raw:
                    text=io.TextIOWrapper(raw,encoding="utf-8-sig",errors="replace",newline="")
                    sample=[text.readline() for _ in range(3)]
                    sample=[x for x in sample if x]
                    if sample:
                        delim='\t' if suffix==".tsv" else (',' if ',' in sample[0] else '\t')
                        header=next(csv.reader(sample[:1],delimiter=delim),[])
                        m["first_row_field_names"]=[s.strip()[:120] for s in header][:120]
                        m["detected_delimiter"]=delim
                        m["first_three_line_lengths"]=[len(x) for x in sample]
            output.append(m)
    return {"schema":"chocho_bsbi_historical_hectad_original_source_schema_20261029_v01",
            "status":"PUBLISHER_MD5_ORIGINAL_ZIP_SCHEMA_VERIFIED",
            "source_page":SOURCE_PAGE,"doi":DOI,"archive_name":FILENAME,
            "publisher_original_archive_md5":digest,
            "original_archive_bytes":path.stat().st_size,
            "source_members":output,
            "source_member_count":len(output),
            "original_total_uncompressed_bytes":sum(m["original_uncompressed_bytes"] for m in output),
            "raw_plant_hectad_records_exported":False,
            "UKBMS_site_ids_or_coordinates_exported":False,
            "requires_source_date_class_and_license_interpretation":True}

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--original-zip",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    result=audit_zip(args.original_zip)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,ensure_ascii=False,sort_keys=True)+"\n",encoding="utf8")
    print(json.dumps({"status":result["status"],"compressed_bytes":result["original_archive_bytes"],
                      "uncompressed_bytes":result["original_total_uncompressed_bytes"],
                      "files":result["source_members"]},ensure_ascii=False))
