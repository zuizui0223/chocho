#!/usr/bin/env python3
"""Full 92-species *lexical* audit of BCE's Clarke-2024 derived rank-1–3 foodplants.

Not taxonomically reconciled. Does not measure independent biological observation,
actual butterfly larval performance, global native/introduced area or colonization.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,re,time,urllib.parse
from collections import defaultdict,Counter
from pathlib import Path
from bs4 import BeautifulSoup
from audit_bce_clarke_evidence_host_pilot import (
    TAXONOMY_URL, download, table_rows, valid_species_binomial, clean
)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--protocol",required=True,type=Path)
    ap.add_argument("--focal-csv",required=True,type=Path)
    ap.add_argument("--hosts-csv",required=True,type=Path)
    ap.add_argument("--output",required=True,type=Path)
    args=ap.parse_args()
    p=json.loads(args.protocol.read_text(encoding="utf-8"))
    if p["schema"]!="chocho_bce_clarke2024_all_exact_europe_species_link_audit_v0.1":
        raise RuntimeError("Wrong full sample protocol")
    with args.focal_csv.open(encoding="utf-8",newline="") as f:
        focus={r["species"].strip() for r in csv.DictReader(f)}
    if len(focus)!=239:raise RuntimeError("Frozen chocho panel no longer exactly 239")
    html,source=download(TAXONOMY_URL)
    listed=table_rows(html,"butterflies")
    if len(listed)<400:raise RuntimeError("BCE species table not verified")
    bce={f"{g} {s}" for g,s in listed if valid_species_binomial(g,s)}
    selected=sorted(focus&bce)
    if len(selected)!=92:
        raise RuntimeError("Frozen exact 92-intersection drift: "+str(len(selected)))
    host_map=defaultdict(set)
    with args.hosts_csv.open(encoding="utf-8-sig",newline="") as f:
        for r in csv.DictReader(f):
            insect=clean((r.get("Insect Genus") or "")+" "+(r.get("Insect Species") or ""))
            if insect not in focus:continue
            g=(r.get("Hostplant Genus") or "").strip()
            sp=(r.get("Hostplant Species") or "").strip()
            if valid_species_binomial(g,sp):
                host_map[insect].add(g+" "+sp)
    rows=[]
    for index,sp in enumerate(selected):
        genus,epithet=sp.split(" ")
        url="https://www.bc-europe.eu/butterfly.php?"+urllib.parse.urlencode({"genus":genus,"species":epithet})
        sourcefile=None
        status="SOURCE_UNAVAILABLE"
        plants=set()
        discarded=[]
        error=None
        for attempt in range(2):
            try:
                contents,sourcefile=download(url)
                tree=BeautifulSoup(contents,"html.parser")
                h1=tree.find("h1")
                if h1 is None or sp not in clean(h1.get_text(" ",strip=True)):
                    raise RuntimeError("BCE page scientific binomial mismatch")
                botanical=table_rows(contents,"plants")
                if not botanical:
                    status="NO_SPECIES_RANK_TABLE"
                    break
                for g,epi in botanical:
                    if valid_species_binomial(g,epi):plants.add(g+" "+epi)
                    else:discarded.append((g+" "+epi).strip())
                status="SOURCE_VERIFIED" if plants else "ONLY_UNRESOLVED_RANK"
                break
            except Exception as exc:
                error=f"{type(exc).__name__}: {str(exc)[:180]}"
                if attempt==0:time.sleep(1.5)
        raw=host_map[sp]
        rows.append({
            "species":sp,"site_status":status,"site_provenance":sourcefile,
            "source_error":error if status=="SOURCE_UNAVAILABLE" else None,
            "BCE_rank1_to3_raw_species_binomials":len(plants),
            "HOSTS_original_raw_species_binomials":len(raw),
            "intersection_exact_spelling":len(plants&raw),
            "BCE_only_exact_spelling":sorted(plants-raw),
            "HOSTS_only_exact_spelling":sorted(raw-plants),
            "BCE_excluded_non_species_taxa":sorted(set(discarded)),
        })
        print(f"SOURCE {index+1}/{len(selected)} {sp} {status} "+
              f"BCE={len(plants)} HOSTS={len(raw)} both={len(plants&raw)}",flush=True)
        time.sleep(1.0)
    verified=[r for r in rows if r["site_status"]=="SOURCE_VERIFIED"]
    failures=Counter(r["site_status"] for r in rows)
    pooled={
        "verified_butterflies":len(verified),
        "BCE_rank1to3_host_links":sum(x["BCE_rank1_to3_raw_species_binomials"] for x in verified),
        "HOSTS_raw_exact_host_links":sum(x["HOSTS_original_raw_species_binomials"] for x in verified),
        "exact_spelling_overlap":sum(x["intersection_exact_spelling"] for x in verified),
        "BCE_only_exact_names":sum(len(x["BCE_only_exact_spelling"]) for x in verified),
        "HOSTS_only_exact_names":sum(len(x["HOSTS_only_exact_spelling"]) for x in verified),
        "BCE_non_species_excluded_names":sum(len(x["BCE_excluded_non_species_taxa"]) for x in verified)
    }
    results={
      "schema":"chocho_bce_clarke_all_exact_92_source_names_v0.1",
      "status":"DESCRIPTIVE_DERIVED_CHECKLIST_EXACT_LEXICAL_AUDIT_NOT_CONFIRMED_MISSING_HOSTS",
      "source_basis":"BCE website derived from Clarke 2024 and includes only evidence levels 1-3; not original new field observations",
      "source_taxonomy":source,
      "BCE_taxonomy_rows":len(listed),
      "frozen_chocho_species":len(focus),
      "exact_bce_chocho_intersection":len(selected),
      "species_page_status_counts":dict(failures),
      "descriptive_pooled_on_verified_species":pooled,
      "per_butterfly":rows,
      "claim_limits":[
         "Nomenclatural synonyms, accepted WCVP plant IDs and butterfly synonyms are not resolved",
         "BCE site is a Clarke-derived evidence 1-3 filter without per-link evidence status or source citations",
         "This is only the exact European taxonomic intersection, not a random sample of 239 butterflies",
         "BCE-only exact spellings are NOT verified HOSTS ecological missing links",
         "No inference about local native-host alternatives, realized competition, extinction or fitness",
         "No change to frozen global resource opportunity result or GEB manuscript",
      ]
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(results,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print("FULL_SUMMARY",json.dumps({
        "exact_species":len(selected),"status_counts":dict(failures),
        "pooled":pooled
    }),flush=True)
    if len(verified)<int(0.75*len(selected)):
        raise RuntimeError("Over a quarter BCE species pages lack interpretable species-rank evidence")


if __name__=="__main__":
    main()
