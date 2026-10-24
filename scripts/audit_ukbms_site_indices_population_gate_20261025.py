#!/usr/bin/env python3
"""Audit source schema/feasibility of original UKBMS site indices 2024.

NO ecological effect estimated. No raw site-level rows, identifiers,
coordinates, count values, or recorder metadata are retained in output.
"""
from __future__ import annotations
import argparse
import collections
import csv
import hashlib
import io
import json
import math
import os
from pathlib import Path
import re
import tempfile
import urllib.error
import urllib.request
import zipfile

ROOT = Path(__file__).resolve().parents[1]
URL = "https://data-package.ceh.ac.uk/sd/53a1c705-4437-48bb-9bcd-b68977f8edd9.zip"
BLOB = "9f83d2d9ac8e46f1cfd4ded572c7a023273503dd"
LIMIT_BYTES = 250 * 1024 * 1024
DATASET_DOI = "10.5285/53a1c705-4437-48bb-9bcd-b68977f8edd9"
ATtribution = ("Contains UK Butterfly Monitoring Scheme (UKBMS) data © copyright "
               "and database right Butterfly Conservation, the UK Centre for Ecology & "
               "Hydrology, British Trust for Ornithology, and the Joint Nature "
               "Conservation Committee.")
SITE_ALIASES = ("site", "site_id", "siteid", "transectid", "transect_id", "transect", "sitename", "sitecode")
SPECIES_ALIASES = ("species", "speciesname", "species_name", "taxon", "taxonname", "scientificname", "scientific_name")
YEAR_ALIASES = ("year", "surveyyear", "survey_year", "recordingyear", "monitoringyear")
INDEX_ALIASES = ("siteindex", "site_index", "index", "abundance", "annualindex", "annual_index", "gai", "site_index_value")
REPEATS = 100


def clean_col(x):
    return re.sub(r"[^a-z0-9]+", "", str(x).lower())


def choose_column(headers, synonyms):
    mapping = collections.defaultdict(list)
    for x in headers:
        mapping[clean_col(x)].append(x)
    options = []
    for name in synonyms:
        options.extend(mapping.get(clean_col(name), []))
    return options[0] if len(set(options)) == 1 or options else None


def original_focal():
    file = ROOT / "data/frozen/figure_sources/anthropogenic_species_metrics.csv"
    raw = file.read_bytes()
    sha = hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\x00" + raw).hexdigest()
    if sha != BLOB:
        raise ValueError("Original 239 butterfly input differs from pinned Git blob")
    with file.open(encoding="utf-8", newline="") as handle:
        names = [row["species"].strip() for row in csv.DictReader(handle)]
    if len(names) != 239 or len(set(names)) != 239:
        raise ValueError("Expected 239 exactly distinct species")
    return set(names)


def fetch_archive(tmp):
    info = {"url":URL, "status":"NOT_REQUESTED"}
    req = urllib.request.Request(URL, headers={
        "User-Agent":"chocho-open-science-ukbms-source-audit/0.1 (contact: GitHub repository zuizui0223/chocho)",
        "Accept":"application/zip, application/octet-stream, */*",
    })
    digest = hashlib.sha256()
    try:
        with urllib.request.urlopen(req, timeout=55) as response:
            info["http_status"] = response.status
            info["content_type"] = response.headers.get("Content-Type","")
            info["content_length_header"] = response.headers.get("Content-Length")
            info["resolved_url"] = response.url.split("?")[0]
            if response.status != 200:
                info["status"] = "HTTP_NOT_OK"
                return info
            total=0
            with open(tmp,"wb") as fh:
                while True:
                    block=response.read(1024*1024)
                    if not block:
                        break
                    total+=len(block)
                    if total > LIMIT_BYTES:
                        info["status"]="TOO_LARGE_SOURCE"
                        info["bytes_read_before_cap"]=total
                        return info
                    digest.update(block)
                    fh.write(block)
        info.update(status="RETRIEVED",sha256=digest.hexdigest(),bytes=total)
    except urllib.error.HTTPError as e:
        info.update(status="HTTP_ERROR",http_status=e.code,error_type=type(e).__name__)
    except (urllib.error.URLError,TimeoutError,OSError) as e:
        info.update(status="SOURCE_TRANSPORT_ERROR",error_type=type(e).__name__,
                    error_brief=str(e)[:160])
    return info


def read_schema(zf, member, focal):
    result={"file":member.filename,"uncompressed_bytes":member.file_size}
    known_sites=set()
    known_names=set()
    years=collections.Counter()
    eligible=collections.Counter()
    repeated=collections.defaultdict(set)
    indices_nonmissing=0
    missing_indices=0
    n_rows=0
    valid_year_rows=0
    nonnumeric_year_rows=0
    bad_year_rows=0
    sample_extensions=[]
    try:
        with zf.open(member,"r") as bytes_stream:
            stream=io.TextIOWrapper(bytes_stream,encoding="utf-8-sig",errors="replace",newline="")
            head=stream.readline()
            for delimiter in [",","\t",";"]:
                try:
                    nfields=len(next(csv.reader([head],delimiter=delimiter)))
                except Exception:
                    nfields=0
                sample_extensions.append((nfields,delimiter))
            sep=max(sample_extensions)[1]
            stream.seek(0)
            reader=csv.DictReader(stream,delimiter=sep)
            headers=reader.fieldnames or []
            result["header"] = headers[:120]
            result["delimiter"] = repr(sep)
            cols={key:choose_column(headers,aliases)
                  for key,aliases in [("site",SITE_ALIASES),("species",SPECIES_ALIASES),
                                      ("year",YEAR_ALIASES),("index",INDEX_ALIASES)]}
            result["inferred_columns"]=cols
            if any(v is None for v in cols.values()):
                result["status"]="UNVERIFIED_REQUIRED_SCHEMA"
                return result
            # Scan full valid CSV and tally only aggregate unique counts.
            for row in reader:
                n_rows += 1
                site=str(row.get(cols["site"]) or "").strip()
                species=str(row.get(cols["species"]) or "").strip()
                year=str(row.get(cols["year"]) or "").strip()
                ix=str(row.get(cols["index"]) or "").strip()
                if species:known_names.add(species)
                if site:known_sites.add(site)
                try:
                    yy=int(year)
                except ValueError:
                    nonnumeric_year_rows+=1
                    continue
                if not 1970 <= yy <= 2026:
                    bad_year_rows+=1
                    continue
                years[yy]+=1
                valid_year_rows+=1
                if not ix or ix.upper() in {"NA","NULL","NAN","-"}:
                    missing_indices+=1
                    continue
                try:
                    val=float(ix)
                except ValueError:
                    missing_indices+=1
                    continue
                if not math.isfinite(val) or val < 0:
                    missing_indices+=1
                    continue
                indices_nonmissing += 1
                if site and species in focal:
                    eligible[species]+=1
                    repeated[(site,species)].add(yy)
    except (UnicodeError, csv.Error, OSError, zipfile.BadZipFile) as e:
        result.update(status="SOURCE_SCHEMA_ERROR",error_type=type(e).__name__)
        return result
    across_2=sum(len(v)>=2 for v in repeated.values())
    across_3=sum(any(y+1 in v and y+2 in v for y in v) for v in repeated.values())
    result.update(status="PARSED_LONG_SITE_YEAR",rows=n_rows,valid_year_rows=valid_year_rows,
                  year_range=[min(years),max(years)] if years else None,
                  years_distinct=len(years),
                  years_coverage_aggregate=dict(sorted(years.items())) if len(years)<=60 else {},
                  species_name_count=len(known_names),exact_239_species_matches=len(known_names & focal),
                  example_exact_focal_names=sorted(known_names & focal)[:38],
                  full_source_sites_distinct=len(known_sites),
                  valid_nonnegative_index_rows=indices_nonmissing,
                  missing_or_invalid_index_rows=missing_indices,
                  invalid_year_rows=nonnumeric_year_rows+bad_year_rows,
                  eligible_focal_species_n=len(eligible),
                  original_focal_species_with_30plus_rows=sum(v>=30 for v in eligible.values()),
                  focal_species_site_pairs_2plus_years=across_2,
                  focal_species_site_pairs_3consecutive_years=across_3,
                  source_species_name_tallies_omitted=True,
                  source_site_ids_or_individual_indices_persisted=False)
    result["time_series_gate_pass"]=bool(
        len(eligible)>=15 and len(known_sites)>=100 and across_3>=REPEATS
    )
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    original=original_focal()
    result={
        "schema":"chocho_ukbms_site_abundance_source_gate_20261025_v02",
        "dataset_doi":DATASET_DOI,
        "source_terms":"Open Government Licence UKBMS with UKBMS attribution",
        "required_attribution":ATtribution,
        "focal_species_count":len(original),
        "source_type":"observed UKBMS annual site index (relative monitored population size), not GBIF observation counts",
        "no_butterfly_abundance_effect_estimated":True,
        "no_competition_coefficient_estimated":True,
        "limit_bytes":LIMIT_BYTES,
    }
    with tempfile.TemporaryDirectory(prefix="ukbms_source_") as temp:
        tmp=Path(temp)/"original_package"
        download=fetch_archive(tmp)
        result["transport"]=download
        if download["status"]!="RETRIEVED":
            result["status"]="SOURCE_ACCESS_UNAVAILABLE_NO_BIOLOGICAL_INFERENCE"
        elif not zipfile.is_zipfile(tmp):
            result["status"]="NOT_ZIP_OR_UNVERIFIED_SOURCE"
        else:
            with zipfile.ZipFile(tmp) as z:
                error=z.testzip()
                if error:
                    result["status"]="CORRUPT_ZIP_MEMBER"
                    result["bad_member"]=error
                else:
                    files=[m for m in z.infolist() if not m.is_dir()]
                    result["archive_members"]=[
                        {"basename":Path(m.filename).name,"bytes":m.file_size}
                        for m in files if not m.filename.startswith("__MACOSX")]
                    candidate=[m for m in files if m.filename.lower().endswith((".csv",".tsv",".txt"))]
                    result["tabular_results"]=[]
                    for member in sorted(candidate,key=lambda m:-m.file_size)[:4]:
                        result["tabular_results"].append(read_schema(z,member,original))
                    valid=[a for a in result["tabular_results"] if a["status"]=="PARSED_LONG_SITE_YEAR"]
                    pass_sites=[a for a in valid if a["time_series_gate_pass"]]
                    result["status"]=(
                        "FEASIBLE_SITE_YEAR_ABUNDANCE_NEEDS_TAXONOMY_AND_MISSINGNESS_REVIEW"
                        if pass_sites else
                        "SOURCE_SCHEMA_OR_TIME_SERIES_GATE_NOT_ESTABLISHED"
                    )
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({
        "status":result["status"],
        "download_status":result["transport"]["status"],
        "http_status":result["transport"].get("http_status"),
        "archive_bytes":result["transport"].get("bytes"),
        "tabular": [
            {k:v for k,v in tab.items() if k in (
                "file","status","inferred_columns","rows","year_range",
                "exact_239_species_matches","full_source_sites_distinct",
                "focal_species_site_pairs_3consecutive_years","time_series_gate_pass")}
            for tab in result.get("tabular_results",[])
        ],
    },ensure_ascii=False,indent=2))


if __name__=="__main__":
    main()
