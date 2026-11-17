#!/usr/bin/env python3
"""Uniform independent larval-host source census for all Canary INBO butterflies.

Every one of the original 33 present Canary butterflies gets the same
deterministic single-source query, irrespective of prior HOSTS coverage and
endemism. Adult flower visits are NEVER promoted to caterpillar feeding.
HTML sections are harvested as VERBATIM evidence for manual genus adjudication;
this first pass only measures source/page availability and larval-section
presence, not the true absence of hosts in nature.
"""
import argparse
import csv
import datetime as dt
import hashlib
import json
import re
import time
import unicodedata
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen
from bs4 import BeautifulSoup

BASE = Path("provenance/reviewer_defenses/results/independent_island_butterfly_host_coverage_20261117")
INBO42 = BASE / "independent_42_species_vs_frozen239.csv"
HERE = Path("provenance/reviewer_defenses/results/island_endemic_synonym_origin_20261118")
BIOTA20 = HERE / "canary_BIOTA_source_status_20.csv"
BIOTA13 = HERE / "canary_BIOTA_remaining_13.csv"
SOURCE = "https://polinizadoresdecanarias.org/portfolio/"
# Website's independently published alternative binomials; NOT a HOSTS-driven matching.
SOURCE_NAME_ALIASES = {
    "Cyclyrius webbianus": "Leptotes webbianus",
    "Colias crocea": "Colias croceus",
}

def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))

def slug(name):
    value = unicodedata.normalize("NFKD", name).encode("ascii", "ignore").decode()
    return re.sub(r"[^a-z0-9]+", "-", value.lower()).strip("-")

def fetch_html(url, destination, offline=False):
    if destination.exists():
        data = destination.read_bytes()
        if data:
            return data, "cached"
    if offline:
        return b"", "offline_missing"
    request = Request(url, headers={
        "User-Agent":"chocho-independent-canary-host-source-audit/2026-11-18 (academic study)",
        "Accept":"text/html,application/xhtml+xml"
    })
    try:
        with urlopen(request, timeout=22) as f:
            if f.status != 200:
                return b"", "http_"+str(f.status)
            data = f.read(3_000_000+1)
            if len(data)>3_000_000:
                return b"", "oversized_page"
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(data)
        return data, "http_200"
    except HTTPError as exc:
        return b"", "http_"+str(exc.code)
    except (URLError, TimeoutError, OSError) as exc:
        return b"", "network_error_"+type(exc).__name__

def life_section(data, expected_source_binomial):
    soup = BeautifulSoup(data, "html.parser")
    for item in soup(["script","style","noscript","nav","footer"]):
        item.decompose()
    title = soup.title.get_text(" ", strip=True) if soup.title else ""
    if expected_source_binomial.lower() not in title.lower():
        return title,"","name_not_in_page_title"
    text = soup.get_text(" ",strip=True)
    canonical = re.sub(r"\s+"," ",text).strip()
    pos = re.search(r"\bCICLO\s+VITAL\b",canonical,re.I)
    if not pos:
        return title,"","no_ciclo_vital_heading"
    section = canonical[pos.end():]
    # Source page repeats heading "Biota" after the larval-life-cycle section.
    end = re.search(r"\b(?:Biota|Ficha\s+de\s+la\s+especie|Descubre\s+m[aá]s)\b",section,re.I)
    if end:
        section=section[:end.start()]
    section=section.strip()
    if not section:
        return title,"","empty_life_section"
    # Keep raw section; never elevate adult plant visits outside CICLO VITAL.
    if re.search(r"\b(?:larvas?|orugas?|oruga)\b",section,re.I):
        return title,section[:3000],"larval_section_text"
    return title,section[:3000],"life_section_without_explicit_larvae"

def run(output, offline=False):
    original=read_csv(INBO42)
    canary=[x for x in original if "Canary Islands" in x["archipelagos"].split("; ")]
    if len(canary)!=33:
        raise ValueError("Independent 33 Canary universe changed")
    statuses=read_csv(BIOTA20)+read_csv(BIOTA13)
    status={x["independent_checklist_name"]:x for x in statuses}
    if len(status)!=33 or len(set(status))!=33:
        raise ValueError("Canary source origin statuses are not exhaustive")
    if set(status)!={x["scientific_name"] for x in canary}:
        raise ValueError("The same 33 checklist butterflies must enter the uniform query")
    records=[]
    out=output / "page_receipts"
    output.mkdir(parents=True,exist_ok=True)
    for item in sorted(canary,key=lambda x:x["scientific_name"]):
        name=item["scientific_name"]
        source_name=SOURCE_NAME_ALIASES.get(name,name)
        url=SOURCE+slug(source_name)+"/"
        data,fetch_status=fetch_html(url,out/(slug(name)+".html"),offline)
        title,section,section_status=life_section(data,source_name) if data else ("","","not_fetched")
        row={
            "independent_checklist_species":name,
            "family":item["checklist_family"],
            "canary_origin":status[name]["canary_origin_stratum"],
            "origin_certainty":status[name]["origin_certainty"],
            "NHM_HOSTS_exact_genus":int(item["hosts_named_genus"]),
            "in_original239":int(item["in_frozen239"]),
            "uniform_source_name":source_name,
            "uniform_source_url":url,
            "fetch_status":fetch_status,
            "source_page_identity_ok":int(bool(data) and "name_not_in_page_title" != section_status),
            "source_page_sha256":hashlib.sha256(data).hexdigest() if data else "",
            "life_section_status":section_status,
            "source_ciclo_vital_larva_text_present":int(section_status=="larval_section_text"),
            "source_life_section_verbatim":section,
        }
        records.append(row)
        print("SOURCE_PAGE="+json.dumps({
            "name":name,"origin":row["canary_origin"],
            "source_name":source_name,"status":fetch_status,"life":section_status,
            "excerpt":section[:420]
        },ensure_ascii=False),flush=True)
        if not offline:
            time.sleep(0.25)
    if len(records)!=33 or len({x["independent_checklist_species"] for x in records})!=33:
        raise ValueError("Uniform source query did not cover every independent Canary species")
    with (output/"uniform_polinizadores_ciclo_vital_canary33.csv").open("w",encoding="utf-8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    cohorts={}
    for origin in sorted({x["canary_origin"] for x in records}):
        sub=[x for x in records if x["canary_origin"]==origin]
        cohorts[origin]={
            "original_independent_n":len(sub),
            "page_http_200_n":sum(x["fetch_status"] in ("http_200","cached") for x in sub),
            "source_identity_verified_n":sum(x["source_page_identity_ok"] for x in sub),
            "life_section_explicit_larva_n":sum(x["source_ciclo_vital_larva_text_present"] for x in sub),
            "hosts_original_exact":sum(x["NHM_HOSTS_exact_genus"] for x in sub)
        }
    summary={
        "schema":"chocho_uniform_single_source_canary_butterfly_larval_host_audit_v1",
        "date_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
        "uniform_source":"Polinizadores de Canarias species accounts, CICLO VITAL section only",
        "source_inventory":"https://polinizadoresdecanarias.org/repositorio/",
        "selected_from":"All 33 INBO Canary present butterflies, not selected from NHM HOSTS",
        "universe_n":33,
        "queried_n":len(records),
        "cohorts":cohorts,
        "all_statuses":{r["independent_checklist_species"]:{
            "origin":r["canary_origin"],"fetch":r["fetch_status"],
            "life_section":r["life_section_status"],"source_url":r["uniform_source_url"]
        } for r in records},
        "strict_limitations":[
            "Not every butterfly has a species account at this source. Missing pages are NOT no larval host.",
            "Source curated accounts describe larvae variably; the study cannot treat documentation as direct field proof.",
            "Text presence is not yet plant-genus resolution; requires inspection of EVERY section before coding reported host genera.",
            "Do NOT read PLANTAS VISITADAS/adult nectar visits as larval foodplants.",
            "Page/source archive SHA256 allows later verification; no original GEB butterfly host or flora reconstruction changed."
        ]
    }
    (output/"uniform_source_summary.json").write_text(json.dumps(summary,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print("UNIFORM_CANARY33_QUERY_SUMMARY="+json.dumps({
        "n":len(records),"cohorts":cohorts,
        "http200":sum(r["fetch_status"]=="http_200" for r in records),
        "larval_text":sum(r["source_ciclo_vital_larva_text_present"] for r in records)
    },ensure_ascii=False))
    return summary

if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir",required=True,type=Path)
    ap.add_argument("--offline",action="store_true")
    a=ap.parse_args()
    run(a.out_dir,a.offline)
