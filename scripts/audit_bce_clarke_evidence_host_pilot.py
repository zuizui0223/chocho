#!/usr/bin/env python3
"""Clarke (2024) evidence 1-3 *BCE-derived* source access and exact-link pilot.

BCE website reflects publication source, not independent new observations.
Every missing exact lexical match needs independent taxonomy/source adjudication.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

from bs4 import BeautifulSoup

TAXONOMY_URL = "https://www.bc-europe.eu/taxonomy.php"
BIOLOGICAL_LINK_SOURCE = "Clarke (2024) doi:10.1002/ece3.10834 derived BCE filtered evidence 1-3"


def download(url):
    req=urllib.request.Request(url,headers={
        "User-Agent":"chocho-scientific-evidence-audit/0.1 (published public information)",
        "Accept":"text/html,application/xhtml+xml;q=0.9,*/*;q=0.8",
    })
    with urllib.request.urlopen(req,timeout=30) as r:
        if r.status!=200:raise RuntimeError(f"HTTP {r.status} from {url}")
        data=r.read()
    if len(data)<1500:raise RuntimeError(f"Unexpectedly short page {url}")
    return data,{"url":url,"bytes":len(data),"sha256":hashlib.sha256(data).hexdigest()}


def clean(value):
    return re.sub(r"\s+"," ",value.strip())


def table_rows(html,kind):
    soup=BeautifulSoup(html,"html.parser")
    for table in soup.find_all("table"):
        parsed=[]
        for tr in table.find_all("tr"):
            cells=tr.find_all("td",recursive=False)
            if len(cells)<7:continue
            pieces=[clean(c.get_text(" ",strip=True)) for c in cells]
            if not pieces[0].isdigit():continue
            if kind=="butterflies":
                # taxonomic listing: number, family, subfamily, genus, species,...
                if not pieces[1].endswith("idae") or not pieces[3]:continue
                parsed.append((pieces[3],pieces[4]))
            else:
                # foodplant table: number, order, family, genus, species,...
                if not pieces[2].endswith("aceae") and pieces[2] not in {"Poaceae","Fabaceae"}:
                    continue
                parsed.append((pieces[3],pieces[4]))
        if len(parsed)>=(100 if kind=="butterflies" else 1):
            return parsed
    return []


def valid_species_binomial(genus,epithet):
    return bool(re.fullmatch(r"[A-Z][A-Za-z-]+",genus) and re.fullmatch(r"[a-z][a-z-]+",epithet))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--protocol",required=True,type=Path)
    ap.add_argument("--focal-csv",required=True,type=Path)
    ap.add_argument("--hosts-csv",required=True,type=Path)
    ap.add_argument("--output",required=True,type=Path)
    a=ap.parse_args()
    p=json.loads(a.protocol.read_text(encoding="utf-8"))
    if p["schema"]!="chocho_bce_clarke2024_evidence_filtered_host_pilot_v0.1":
        raise RuntimeError("Unexpected protocol")
    with a.focal_csv.open(encoding="utf-8",newline="") as f:
        focal={r["species"].strip() for r in csv.DictReader(f)}
    if len(focal)!=239:raise RuntimeError(f"Frozen butterfly panel drift: {len(focal)}")
    index,indexsource=download(TAXONOMY_URL)
    site_taxa=table_rows(index,"butterflies")
    if len(site_taxa)<400:
        raise RuntimeError(f"BCE Europe 501-species listing HTML parse mismatch; observed {len(site_taxa)}")
    site={f"{g} {s}" for g,s in site_taxa if valid_species_binomial(g,s)}
    intersection=sorted(focal&site)
    if len(intersection)<12:
        raise RuntimeError(f"Insufficient exact 239 x BCE site species matches: {len(intersection)}")
    selected=intersection[:12]
    original=defaultdict(set)
    with a.hosts_csv.open(encoding="utf-8-sig",newline="") as file:
        for r in csv.DictReader(file):
            b=clean((r.get("Insect Genus") or "")+" "+(r.get("Insect Species") or ""))
            if b not in selected:continue
            g=(r.get("Hostplant Genus") or "").strip()
            s=(r.get("Hostplant Species") or "").strip()
            if valid_species_binomial(g,s):
                original[b].add(g+" "+s)
    audited=[]
    for butterfly in selected:
        genus,epithet=butterfly.split(" ")
        url="https://www.bc-europe.eu/butterfly.php?"+urllib.parse.urlencode({"genus":genus,"species":epithet})
        status="SOURCE_UNAVAILABLE"
        source=None
        extra=None
        foodplants=set()
        excluded=[]
        for attempt in range(2):
            try:
                page,source=download(url)
                soup=BeautifulSoup(page,"html.parser")
                title=soup.find("h1")
                if title is None or butterfly not in clean(title.get_text(" ",strip=True)):
                    raise RuntimeError("BCE original named species page does not match queried taxon")
                all_pairs=table_rows(page,"plants")
                if not all_pairs:
                    status="NO_VERIFIED_FOODPLANT_TABLE"
                    break
                for g,s in all_pairs:
                    if valid_species_binomial(g,s):
                        foodplants.add(g+" "+s)
                    else:
                        excluded.append(f"{g} {s}".strip())
                status="SOURCE_VERIFIED" if foodplants else "ONLY_UNRESOLVED_PLANT_NAMES"
                break
            except Exception as err:
                extra=f"{type(err).__name__}:{str(err)[:180]}"
                if attempt==0:time.sleep(1)
        ref=original[butterfly]
        audited.append({
            "butterfly":butterfly,"source_status":status,"source_provenance":source,
            "source_error":extra if status=="SOURCE_UNAVAILABLE" else None,
            "BCE_evidence_1_to_3_exact_species_count":len(foodplants),
            "HOSTS_exact_original_species_count":len(ref),
            "exact_lexical_overlap":len(foodplants&ref),
            "BCE_only_exact_binomials":sorted(foodplants-ref),
            "HOSTS_only_exact_binomials":sorted(ref-foodplants),
            "BCE_excluded_non_species_or_infraspecific":sorted(set(excluded)),
        })
        print("BCE_SOURCE",butterfly,status,"BCE_exact",len(foodplants),
              "HOSTS_exact",len(ref),"overlap",len(foodplants&ref),flush=True)
        time.sleep(1.1)
    good=[x for x in audited if x["source_status"]=="SOURCE_VERIFIED"]
    output={
       "schema":"chocho_bce_clarke_high_evidence_species_host_pilot_v0.1",
       "status":"POSTHOC_SOURCE_PILOT_ONLY_NO_GLOBAL_INFERENCE",
       "source_label":BIOLOGICAL_LINK_SOURCE,
       "taxonomy_source":indexsource,
       "BCE_501_species_page_rows":len(site_taxa),
       "BCE_unambiguous_species":len(site),
       "frozen_chocho_butterflies":len(focal),
       "exact_species_intersection":len(intersection),
       "pilot_selection_rule":"alphabetically first 12 exact species in intersection; no host outcomes used for selection",
       "pilot_attempts":len(audited),"pilot_verified":len(good),
       "pilot_results":audited,
       "pooled_verified_only":{
           "species":len(good),
           "BCE_exact_host_links":sum(x["BCE_evidence_1_to_3_exact_species_count"] for x in good),
           "HOSTS_exact_original_host_links":sum(x["HOSTS_exact_original_species_count"] for x in good),
           "overlap_exact_spelling":sum(x["exact_lexical_overlap"] for x in good),
           "BCE_only_exact_spelling":sum(len(x["BCE_only_exact_binomials"]) for x in good),
       },
       "stop_boundaries":[
           "This is an independently curated DERIVATIVE listing, not a new field observation or original Clarke Dryad table.",
           "BCE site includes only evidence ranks 1-3 but suppresses individual tier and source citations.",
           "Exact spelling mismatch is not biological absence: WCVP accepted-name and butterfly synonyms are not resolved.",
           "No extrapolation of 12 alphabetically chosen European butterflies to all 239.",
           "No WCVP geographic counterfactual, butterfly fitness or colonization claim.",
       ]
    }
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(output,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print("SUMMARY",json.dumps({k:output[k] for k in ("BCE_501_species_page_rows",
        "exact_species_intersection","pilot_attempts","pilot_verified",
        "pooled_verified_only")}),flush=True)
    if len(good)<6:
        raise RuntimeError(f"Insufficient BCE host-page coverage: verified {len(good)}/12, no inference")


if __name__=="__main__":
    main()
