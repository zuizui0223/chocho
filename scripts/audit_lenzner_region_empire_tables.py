#!/usr/bin/env python3
"""Inspect exact Lenzner et al. 2022 GloNAF/colonial data schemas and join keys.

Source: https://doi.org/10.5281/zenodo.6640264, original MD5 frozen.
Does NOT republish any source row, site coordinates, colonial occupation
values or the 46 MB shapefile. Profile only keys, counts and metadata.
"""
from __future__ import annotations

from collections import Counter
import csv
import hashlib
import io
import json
from pathlib import Path
import zipfile

ZIP=Path("_temp_historical_sources/lenzner_6640264.zip")
OUT=Path("provenance/reviewer_defenses/results/global_colonial_resource_beta_20261031/"
         "lenzner_region_empire_schema.json")
EXPECTED_MD5="16f5585d2500676951a4b664bbf0ffde"
CSV_NAMES=("Empires.csv","Regions.csv","Predictors.csv","GloNAF.csv")


def csv_profile(member, data):
    raw=io.TextIOWrapper(data,encoding="utf-8-sig",errors="replace",newline="")
    sample=raw.read(7000)
    raw.seek(0)
    try:
        dialect=csv.Sniffer().sniff(sample,delimiters=",;\t|")
    except csv.Error:
        dialect=csv.excel
    reader=csv.DictReader(raw,dialect=dialect)
    columns=[str(x).strip() for x in (reader.fieldnames or [])]
    if not columns:
        raise RuntimeError("No CSV header: "+member)
    keylike=[name for name in columns if any(k in name.lower()
             for k in ("id","code","region","empire","colon","occup","duration",
                       "country","start","end","year","time","grb","glonaf"))]
    encountered={k:set() for k in keylike}
    uniques_saturated={}
    empty=Counter()
    n=0
    for row in reader:
        n+=1
        for k in keylike:
            v=(row.get(k) or "").strip()
            if not v:
                empty[k]+=1
            elif len(encountered[k])<100000:
                encountered[k].add(v)
            else:
                uniques_saturated[k]=True
    return {
        "table":member,"rows":n,"columns":columns,
        "keylike_columns":keylike,
        "candidate_key_cardinality":{k:len(v) for k,v in encountered.items()},
        "candidate_key_empty":{k:empty[k] for k in keylike},
        "uniques_truncated_at_100000":uniques_saturated,
        "file_bytes_uncompressed":None,
        "schema_only_no_source_rows":True,
    }


def profile(path=ZIP):
    if not path.exists():
        raise RuntimeError("Original Zenodo zip not downloaded")
    with path.open("rb") as f:
        md5=hashlib.md5()
        for x in iter(lambda:f.read(1<<20),b""):
            md5.update(x)
    if md5.hexdigest()!=EXPECTED_MD5:
        raise RuntimeError("Wrong Lenzner source checksum")
    out={}
    with zipfile.ZipFile(path) as z:
        archive=z.namelist()
        for name in CSV_NAMES:
            found=[k for k in archive if k.endswith("/Datasets/"+name)]
            if len(found)!=1:
                raise RuntimeError("Expected unique "+name+" in archive")
            with z.open(found[0]) as handle:
                # zipfile returns non-seekable stream; wrap in BytesIO
                p=csv_profile(name,io.BytesIO(handle.read()))
            p["file_bytes_uncompressed"]=z.getinfo(found[0]).file_size
            out[name]=p
        shp=[f for f in archive if f.lower().endswith(".shp")]
        dbf=[f for f in archive if f.lower().endswith(".dbf")]
        if len(shp)!=1 or len(dbf)!=1:
            raise RuntimeError("Expected single colonial/GloNAF region shapefile")
        # File header fields via dbf, without data rows.
        import shapefile
        reader=shapefile.Reader(dbf=io.BytesIO(z.read(dbf[0])))
        geom={"shapefile_uncompressed_bytes":z.getinfo(shp[0]).file_size,
              "geometry_record_count_from_dbf":len(reader),
              "attributes":[{"name":f[0],"dbf_type":f[1]}
                            for f in reader.fields[1:]],
              "geometry_not_crosswalked":True}
    return {
        "status":"COLONIAL_GLO_NAF_DATA_SCHEMA_VERIFIED_NOT_WGSRPD3_MAPPING",
        "source_doi":"10.5281/zenodo.6640264",
        "zip_md5":EXPECTED_MD5,
        "csv_profiles":out,
        "region_geometry":geom,
        "limitations":[
            "Schema and record counts do not establish a TDWG level3 spatial crosswalk.",
            "A region belonging to multiple colonial powers needs time-specific empire overlap, not a single empire label.",
            "GloNAF naturalized-alien records and all WCVP introduced ranges use different definitions and geographic units.",
            "The published study's colonial-flora association is not an independently novel insect result.",
        ],
    }


def main():
    p=profile()
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(p,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps({"status":p["status"],"tables":{
      k:{"rows":v["rows"],"columns":v["columns"],
         "keys":v["candidate_key_cardinality"]}
      for k,v in p["csv_profiles"].items()},
      "shape":p["region_geometry"]},indent=2))


if __name__=="__main__":
    main()
